"""Project orchestration (pipeline) used by both GUI and headless drivers."""

from __future__ import annotations

import os
from typing import Callable, Dict, List, Optional

from . import (batch, box as boxmod, cofactor, docking, envinfo, interactions,
               inventory, ligand, models, qc, ranking, receptor, reports,
               selfdock, structure as st, utils, visualize)
from .models import BoxDef, RunConfig, project_subdirs

LogCb = Optional[Callable[[str], None]]
ProgCb = Optional[Callable[[dict], None]]

PHASES = ["prepare", "screen", "refine", "selfdock", "analyze", "visualize", "report"]


class Pipeline:
    def __init__(self, cfg: RunConfig, on_log: LogCb = None, on_progress: ProgCb = None):
        self.cfg = cfg
        self.log = on_log
        self.progress = on_progress
        self.out = cfg.out_dir
        self.sub = project_subdirs(self.out)
        models.ensure_project_layout(self.out)
        self.state_path = os.path.join(self.out, "pipeline_state.json")
        self.state = self._load_state()
        # computed maps (persisted in state)
        self.receptors = self.state.get("receptors", {})       # name -> dict(pdbqt, clean)
        self.ligands: Dict[str, dict] = self.state.get("ligands", {})   # lig_id -> entry dict
        self.cocrystal = self.state.get("cocrystal", {})        # rec -> dict
        self.boxes = self.cfg.boxes or self.state.get("boxes", {})
        self.cfg.boxes = self.boxes   # keep a single source of truth

    # -- logging helpers ------------------------------------------------
    def _log(self, m: str):
        if self.log:
            try:
                self.log(m)
            except Exception:
                pass

    def _prog(self, d: dict):
        if self.progress:
            try:
                self.progress(d)
            except Exception:
                pass

    # -- state ------------------------------------------------------------
    def _load_state(self) -> dict:
        if os.path.exists(self.state_path):
            try:
                return utils.read_json(self.state_path)
            except Exception:
                return {}
        return {}

    def _save_state(self):
        utils.write_json({
            "receptors": self.receptors, "ligands": self.ligands,
            "cocrystal": self.cocrystal, "boxes": self.boxes,
        }, self.state_path)

    # ------------------------------------------------------------------ #
    # PHASE: prepare
    # ------------------------------------------------------------------ #
    def phase_prepare(self):
        self._log("== [prepare] input inventory + structure preparation ==")
        if self.receptors and self.ligands and self.boxes:
            self._log("[prepare] already prepared (resume) - skipping preparation")
            return
        # 1. persist config
        cfg_json = os.path.join(self.out, "config.json")
        self.cfg.boxes = self.boxes
        self.cfg.save(cfg_json)
        envinfo.write_environment_json(os.path.join(self.out, "environment.json"))

        # 2. inventory
        inv = inventory.run_inventory(self.cfg, self.out)

        # 3. boxes
        self.boxes = self._decide_boxes()
        self.cfg.boxes = self.boxes
        utils.write_json(self.boxes, os.path.join(self.out, "boxes.json"))

        # 4. receptors
        for rec in self.cfg.receptors:
            nm = rec["name"]
            boxd = self.boxes.get(nm, {})
            self._log(f"[prepare] receptor {nm}")
            self._prog({"stage": "prepare", "current": f"receptor {nm}"})
            # reference points = cocrystal ligand heavy atoms when used for box
            ref_pts = self._cocrystal_reference_points(nm, rec, boxd)
            try:
                pr = receptor.prepare_receptor(
                    rec["path"], self.sub["prepared_receptors"], nm,
                    chains=rec.get("chains") or None,
                    box_center=boxd.get("center"), box_size=boxd.get("size"),
                    reference_points=ref_pts,
                    delete_bad_res_far=self.cfg.delete_bad_res,
                )
                self.receptors[nm] = {
                    "pdbqt": pr["pdbqt"], "clean": pr["cleaned_pdb"],
                    "source": rec["path"], "atom_counts": pr["atom_counts"],
                    "n_atoms_pdbqt": pr["n_atoms_pdbqt"],
                    "deleted_bad": pr["deleted_bad_residues"],
                }
            except Exception as e:
                self._log(f"[prepare] receptor {nm} FAILED: {e}")
                raise

        # 5. ligands
        for lig in self.cfg.ligands:
            base = lig["name"] or os.path.splitext(os.path.basename(lig["path"]))[0]
            lib_dir = os.path.join(self.sub["prepared_ligands"], utils.safe_name(base))
            os.makedirs(lib_dir, exist_ok=True)
            self._log(f"[prepare] ligand library {lig['path']}")
            entries = ligand.prepare_ligand_library(lig["path"], lib_dir, ph=self.cfg.ph,
                                                    library_name=base)
            for e in entries:
                if e.pdbqt_ok:
                    self.ligands[e.name] = {"pdbqt": e.pdbqt, "prep_sdf": e.prep_sdf,
                                            "source": lig["path"]}
            # ligand_list rows already written by inventory; update pdbqt status
        # 6. cocrystal references (selfdock) and ligand pdbqt for native ligand
        self._prepare_cocrystal()

        # write updated ligand_list with pdbqt
        self._save_state()

    def _decide_boxes(self) -> Dict[str, dict]:
        if self.boxes:
            return self.boxes
        boxes = {}
        for rec in self.cfg.receptors:
            nm = rec["name"]
            b = self._auto_box_for(nm, rec)
            boxes[nm] = b.to_dict()
            self._log(f"[box] {nm}: {b.method} center={b.center} size={b.size}")
        return boxes

    def _auto_box_for(self, rec_name: str, rec: dict) -> BoxDef:
        chains = rec.get("chains") or []
        pdb = rec["path"]
        rows = inventory.inventory_receptor(pdb, rec_name, chains)
        cand = None
        for r in rows:
            if r["chain"] in chains and r.get("ligand"):
                if cand is None or (r["ligand_n_heavy"] or 0) > (cand["ligand_n_heavy"] or 0):
                    cand = r
        if cand:
            recs = st.coords_of_ligand(pdb, cand["ligand"], cand["chain"], cand["ligand_resseq"])
            return boxmod.box_from_ligand(recs, min_axis=self.cfg.min_box_axis,
                                          pad=self.cfg.ligand_pad_ang,
                                          max_axis=self.cfg.max_box_axis,
                                          ligand_tag=cand["ligand"])
        # no cocrystal: user must decide (per protocol, never guess silently)
        raise RuntimeError(
            f"receptor {rec_name}: no search box and no co-crystal ligand found in docking "
            f"chains. Please define the binding site (explicit center/size or residue list) "
            f"or explicitly choose blind docking.")

    def _cocrystal_reference_points(self, rec_name: str, rec: dict, boxd: dict) -> Optional[List[list]]:
        lig_ref = boxd.get("ligand_ref")
        if not lig_ref:
            return None
        chains = rec.get("chains") or []
        rows = inventory.inventory_receptor(rec["path"], rec_name, chains)
        for r in rows:
            if r.get("ligand") == lig_ref and (not chains or r["chain"] in chains):
                recs = st.coords_of_ligand(rec["path"], lig_ref, r["chain"], r.get("ligand_resseq"))
                pts = [[a.x, a.y, a.z] for a in recs if a.element != "H"]
                if pts:
                    return pts
        return None

    def _prepare_cocrystal(self):
        """Build crystal-reference SDF + PDBQT for self-docking of native ligands."""
        for rec in self.cfg.receptors:
            nm = rec["name"]
            chains = rec.get("chains") or []
            pdb = rec["path"]
            # identify cocrystal ligand used for the box
            box = self.boxes.get(nm, {})
            lig_ref = box.get("ligand_ref")
            if not lig_ref:
                continue
            # find chain/resseq of that ligand in chosen chains
            rows = inventory.inventory_receptor(pdb, nm, chains)
            found = None
            for r in rows:
                if r.get("ligand") == lig_ref and r["chain"] in chains:
                    found = r
                    break
            if not found:
                # maybe ligand in other chain; still allow
                for r in rows:
                    if r.get("ligand") == lig_ref:
                        found = r
                        break
            if not found:
                self._log(f"[cocrystal] {nm}: ligand {lig_ref} not located; self-dock skipped")
                continue
            ccd_file = os.path.join(os.path.dirname(pdb), f"{lig_ref}.cif")
            if not os.path.exists(ccd_file):
                # try local data dir / CCD cache next to receptor
                ccd_file = os.path.join(self.out, "ccd", f"{lig_ref}.cif")
                os.makedirs(os.path.dirname(ccd_file), exist_ok=True)
            out_sdf = os.path.join(self.sub["cocrystal"], f"{nm}_{lig_ref}.sdf")
            if not os.path.exists(ccd_file):
                # leave marker for user-provided CCD; attempt download to project dir
                self._log(f"[cocrystal] {nm}: CCD file {lig_ref}.cif not found next to PDB; "
                          f"download attempt into {ccd_file}")
                ccd_file = self._download_ccd(lig_ref, ccd_file)
            try:
                stats = cofactor.write_cocrystal_ligand(
                    pdb, lig_ref, ccd_file, out_sdf,
                    chain=found["chain"], resseq=found.get("ligand_resseq"),
                    name=f"{nm}_{lig_ref}")
                self._log(f"[cocrystal] {nm}: {lig_ref} SDF ok "
                          f"(ccd centroid dev {stats.get('ccd_centroid_dev_A')} A)")
                # prepare pdbqt from this SDF
                lig_out_dir = os.path.join(self.sub["prepared_ligands"], "cocrystal")
                os.makedirs(lig_out_dir, exist_ok=True)
                crys_id = f"{nm}_{lig_ref}"
                pdbqt_path = os.path.join(lig_out_dir, f"{crys_id}.pdbqt")
                try:
                    ligand.prepare_cocrystal_pdbqt(out_sdf, pdbqt_path, ph=self.cfg.ph)
                except Exception as e:
                    self._log(f"[cocrystal] {nm}: PDBQT prep failed: {e}")
                    pdbqt_path = ""
                self.cocrystal[nm] = {"resname": lig_ref, "lig_id": crys_id,
                                      "crystal_sdf": out_sdf, "pdbqt": pdbqt_path,
                                      "chain": found["chain"], "resseq": found.get("ligand_resseq"),
                                      "ccd_stats": stats}
            except Exception as e:
                self._log(f"[cocrystal] {nm}: skipped ({e})")
                self.cocrystal[nm] = {"error": str(e), "resname": lig_ref, "lig_id": f"{nm}_{lig_ref}"}

    def _download_ccd(self, resname: str, out_path: str) -> str:
        import urllib.request
        url = f"https://files.rcsb.org/ligands/view/{resname}.cif"
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        urllib.request.urlretrieve(url, out_path)
        return out_path

    # ------------------------------------------------------------------ #
    # PHASE: screen
    # ------------------------------------------------------------------ #
    def phase_screen(self, time_slice_s: Optional[float] = None) -> dict:
        self._log("== [screen] initial docking (exhaustiveness "
                  f"{self.cfg.exhaustiveness}) ==")
        try:
            self._log(f"[screen] docking backend: {docking.resolve_backend()}")
        except Exception as e:
            self._log(f"[screen] docking backend ERROR: {e}")
        if not self.ligands:
            raise RuntimeError("no prepared ligands available for screening")
        pairs = [(r, l) for r in self.receptors for l in self.ligands]
        return batch.run_docking_stage(
            self.cfg, "screening", self.receptors,
            {k: v["pdbqt"] for k, v in self.ligands.items()},
            pairs=pairs, time_slice_s=time_slice_s,
            progress=self._prog, log=self._log, stage_subdir=self.sub["screening"])

    # ------------------------------------------------------------------ #
    # PHASE: refine
    # ------------------------------------------------------------------ #
    def phase_refine(self, time_slice_s: Optional[float] = None) -> dict:
        self._log("== [refine] shortlist re-docking (exhaustiveness "
                  f"{self.cfg.refine_exhaustiveness}) ==")
        results = batch.collect_results(self.sub["screening"])
        if not results:
            raise RuntimeError("no screening results; refine requires screening")
        recs = sorted({r["receptor"] for r in results})
        ligs = sorted({r["ligand"] for r in results})
        mat = ranking.affinity_matrix(results, recs, ligs)
        ranking.write_affinity_matrix(mat, os.path.join(self.sub["results"], "screening_matrix.csv"), ligs)
        pairs = ranking.choose_refine_pairs(self.cfg, mat)
        if not pairs:
            raise RuntimeError("refine shortlist empty")
        return batch.run_docking_stage(
            self.cfg, "refine", self.receptors,
            {k: v["pdbqt"] for k, v in self.ligands.items()},
            pairs=pairs, time_slice_s=time_slice_s,
            progress=self._prog, log=self._log, stage_subdir=self.sub["refine"])

    # ------------------------------------------------------------------ #
    # PHASE: selfdock
    # ------------------------------------------------------------------ #
    def phase_selfdock(self, time_slice_s: Optional[float] = None) -> dict:
        self._log("== [selfdock] crystal-ligand redock validation ==")
        usable = {r: m for r, m in self.cocrystal.items() if m.get("pdbqt")}
        if not usable:
            self._log("[selfdock] no cocrystal references available; skipped")
            return {"done": 0, "skipped": 0, "errors": 0, "total": 0,
                    "remaining": [], "note": "no cocrystal"}
        return batch.run_docking_stage(
            self.cfg, "selfdock", self.receptors,
            ligand_pdbqt_map={},
            cocrystal_pdbqt_map={r: {"pdbqt": m["pdbqt"], "lig_id": m["lig_id"]}
                                 for r, m in usable.items()},
            time_slice_s=time_slice_s, progress=self._prog, log=self._log,
            stage_subdir=os.path.join(self.out, "selfdock"))

    # ------------------------------------------------------------------ #
    # PHASE: analyze (top-K, selfdock eval, PLIP)
    # ------------------------------------------------------------------ #
    def phase_analyze(self) -> dict:
        self._log("== [analyze] final Top-K + selfdock RMSD + PLIP ==")
        results = batch.collect_results(self.sub["screening"])
        recs = sorted({r["receptor"] for r in results})
        ligs = sorted({r["ligand"] for r in results})
        mat = ranking.affinity_matrix(results, recs, ligs)
        refine_res = batch.collect_results(self.sub["refine"])
        top = ranking.write_final_topk(self.cfg, refine_res, mat, self.sub["results"], ligs)

        # selfdock evaluation
        selfdock.run_selfdock_validation(self.cfg, self.cocrystal, self.out,
                                         progress=self._prog, log=self._log)

        # PLIP for final top-K (best pose from refine if available, else screening)
        best_pose_map = {}
        for item in top["rows"]:
            rec = item["receptor"]; lig = item["ligand"]
            bp = self._best_pose_path(rec, lig, prefer="refine")
            if bp:
                best_pose_map[(rec, lig)] = bp
        top_pairs = [{"receptor": r["receptor"], "ligand": r["ligand"],
                      "best_pdbqt": best_pose_map.get((r["receptor"], r["ligand"]))}
                     for r in top["rows"]]
        clean_map = {rec: self.receptors[rec]["clean"] for rec in self.receptors}
        if self.cfg.run_plip:
            try:
                interactions.run_interactions_for_topk(self.cfg, clean_map, top_pairs,
                                                       self.out, log=self._log,
                                                       progress=self._prog)
            except Exception as e:
                self._log(f"[analyze] PLIP failed: {e}")
        # stash top rows + poses for visualization/report
        utils.write_json({"top_rows": top["rows"],
                          "top_pairs": top_pairs}, os.path.join(self.sub["results"], "top_ctx.json"))
        return {"top": top, "n_top": len(top["rows"])}

    def _best_pose_path(self, rec: str, lig: str, prefer: str = "refine") -> Optional[str]:
        cand = []
        if prefer == "refine":
            cand.append(os.path.join(self.sub["refine"], f"{rec}__{lig}", "best.pdbqt"))
        cand.append(os.path.join(self.sub["screening"], f"{rec}__{lig}", "best.pdbqt"))
        for c in cand:
            if os.path.exists(c):
                return c
        return None

    # ------------------------------------------------------------------ #
    # PHASE: visualize
    # ------------------------------------------------------------------ #
    def phase_visualize(self) -> dict:
        self._log("== [visualize] PyMOL rendering ==")
        ctx_path = os.path.join(self.sub["results"], "top_ctx.json")
        if not os.path.exists(ctx_path):
            self._log("[visualize] no top context; skipped")
            return {"note": "no top context"}
        ctx = utils.read_json(ctx_path)
        top_pairs = ctx.get("top_pairs", [])
        rendered = []
        for tp in top_pairs:
            rec = tp["receptor"]; lig = tp["ligand"]
            name = f"{rec}__{lig}"
            best = tp.get("best_pdbqt")
            if not best or not os.path.exists(best):
                continue
            clean = self.receptors.get(rec, {}).get("clean")
            if not clean:
                continue
            workdir = os.path.join(self.sub["results"], "pymol_data", name)
            os.makedirs(workdir, exist_ok=True)
            complex_pdb = os.path.join(workdir, "complex.pdb")
            try:
                interactions._write_complex_pdb(clean, best, complex_pdb)
            except Exception as e:
                self._log(f"[visualize] complex build failed for {name}: {e}")
                continue
            png3 = os.path.join(self.sub["plots3d"], f"{name}.png")
            try:
                st = visualize.render_complex(workdir, name, png3,
                                              include_pse=self.cfg.write_pse)
                rendered.append({**tp, "png3d": png3 if st["ok"] else None,
                                 "log": st["log_tail"]})
                if st.get("ok"):
                    self._log(f"[visualize] {name} ok")
                else:
                    self._log(f"[visualize] {name} FAILED: {st['log_tail'][-600:]}")
                if self.cfg.write_pse and st.get("pse"):
                    pse_dst = os.path.join(self.sub["pse"], f"{name}.pse")
                    os.replace(st["pse"], pse_dst)
            except Exception as e:
                self._log(f"[visualize] {name} error: {e}")
            # 2D & composite
            prep = self.ligands.get(lig, {}).get("prep_sdf")
            if prep and os.path.exists(prep):
                png2 = os.path.join(self.sub["lig2d"], f"{name}.png")
                if visualize.ligand_2d_png(prep, png2):
                    if st.get("ok"):
                        visualize.composite_png(png3, png2, os.path.join(self.sub["composite"], f"{name}.png"))
        # tile per receptor
        for rec in self.receptors:
            imgs = [os.path.join(self.sub["plots3d"], f"{r['receptor']}__{r['ligand']}.png")
                    for r in rendered if r["receptor"] == rec and r.get("png3d")]
            if len(imgs) >= 2:
                visualize.tile_pngs(imgs, os.path.join(self.sub["topk_panels"], f"{rec}_topK.png"))
        utils.write_json({"rendered": rendered}, os.path.join(self.sub["results"], "viz_ctx.json"))
        return {"rendered": len(rendered)}

    # ------------------------------------------------------------------ #
    # PHASE: report
    # ------------------------------------------------------------------ #
    def phase_report(self) -> dict:
        self._log("== [report] methods + quality assessment ==")
        env = envinfo.gather_environment()
        env_path = os.path.join(self.out, "environment.json")
        if os.path.exists(env_path):
            try:
                env = utils.read_json(env_path)
            except Exception:
                pass
        checks = qc.run_qc(self.cfg, self.out)
        self._prog({"stage": "report"})
        selfdock_rows = None
        sd_path = os.path.join(self.sub["results"], "selfdock_summary.json")
        if os.path.exists(sd_path):
            try:
                selfdock_rows = utils.read_json(sd_path).get("selfdock", [])
            except Exception:
                pass
        reports.write_receptor_protocol(self.cfg, self.out, selfdock_rows)
        reports.write_methods_report(self.cfg, self.out, env, {}, self.boxes, self.receptors)
        reports.write_quality_report(self.cfg, self.out, checks)
        reports.write_readme(self.out, self.cfg, env)
        # directory tree listing
        tree = reports.directory_tree(self.out, maxdepth=4)
        with open(os.path.join(self.sub["reports"], "directory_tree.txt"), "w", encoding="utf-8") as fh:
            fh.write(tree)
        return {"checks": checks}


def run_project(cfg: RunConfig, on_log: LogCb = None, on_progress: ProgCb = None,
                phases: Optional[List[str]] = None,
                time_slice_s: Optional[float] = None,
                on_cancel: Optional[Callable[[], bool]] = None) -> dict:
    """Execute a project. Returns a summary dict.

    Docking stages honour ``time_slice_s``; when it expires the pipeline stops
    and returns ``remaining`` so the caller can resume later. ``on_cancel`` is
    polled between phases (safe-point cancellation for the GUI).
    """
    phases = phases or PHASES
    # Use an ASCII (Windows 8.3 short) working path internally so RDKit/meeko/
    # Vina never see non-ASCII folder names (Chinese paths break their narrow
    # char file APIs). The short path points to the SAME user-visible folder.
    if cfg.out_dir:
        cfg.out_dir = utils.short_path(cfg.out_dir)
        if any(ord(ch) > 127 for ch in cfg.out_dir):
            raise RuntimeError(
                "输出目录包含非英文字符,且当前系统无法生成 8.3 短路径。\n"
                "请改用纯英文/数字的输出目录,例如 C:\\Users\\fang\\Desktop\\dock_test 。")
    pipe = Pipeline(cfg, on_log, on_progress)
    summary = {"phases": [], "interrupted": False}
    for ph in phases:
        if on_cancel and on_cancel():
            summary["interrupted"] = True
            pipe._log(f"[pipeline] cancelled before phase {ph}")
            break
        if ph == "prepare":
            pipe.phase_prepare()
        elif ph == "screen":
            r = pipe.phase_screen(time_slice_s=time_slice_s)
            summary["screen"] = r
            for e in r.get("error_list", []):
                pipe._log(f"[screen] FAILED {e['pair']}: {e['error']}")
            if r.get("remaining"):
                summary["interrupted"] = True
                pipe._save_state()
                break
        elif ph == "refine":
            r = pipe.phase_refine(time_slice_s=time_slice_s)
            summary["refine"] = r
            for e in r.get("error_list", []):
                pipe._log(f"[refine] FAILED {e['pair']}: {e['error']}")
            if r.get("remaining"):
                summary["interrupted"] = True
                pipe._save_state()
                break
        elif ph == "selfdock":
            r = pipe.phase_selfdock(time_slice_s=time_slice_s)
            summary["selfdock"] = r
            for e in r.get("error_list", []):
                pipe._log(f"[selfdock] FAILED {e['pair']}: {e['error']}")
            if r.get("remaining"):
                summary["interrupted"] = True
                pipe._save_state()
                break
        elif ph == "analyze":
            pipe.phase_analyze()
        elif ph == "visualize":
            pipe.phase_visualize()
        elif ph == "report":
            pipe.phase_report()
        summary["phases"].append(ph)
    pipe._save_state()
    return summary
