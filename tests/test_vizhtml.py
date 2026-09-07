"""Unit tests for the interactive HTML 3D viewer / report (v2.0 feature #17)."""

from dockstudio.core import vizhtml


PDB = ("ATOM      1  N   MET A   1       1.000   2.000   3.000  1.00 20.00           N  \n"
       "HETATM 999  C1  LIG A   2       4.000   5.000   6.000  1.00 20.00           C  \n")


def test_local_snippet_is_a_script_loader():
    """_local_3dmol_snippet must return a <script> loader (CDN fallback ok)."""
    snippet = vizhtml._local_3dmol_snippet()
    assert "<script" in snippet
    assert snippet.strip()


def test_viewer_embeds_3dmol_loader_when_provided():
    """A generated viewer must load 3Dmol.js (regression: loader was dropped)."""
    loader = '<script src="https://cdn.example/3Dmol-min.js"></script>'
    html = vizhtml._viewer_html("R__L", PDB, [], [], used_cdn=True,
                                mol_snippet=loader)
    assert "https://cdn.example/3Dmol-min.js" in html
    # the viewer initialisation references $3Dmol, so the loader must be present
    assert html.index("$3Dmol") > html.index("cdn.example")


def test_viewer_reports_cdn_note_when_using_cdn():
    loader = '<script src="https://cdn.example/3Dmol-min.js"></script>'
    html = vizhtml._viewer_html("R__L", PDB, [], [], used_cdn=True,
                                mol_snippet=loader)
    assert "联网" in html
