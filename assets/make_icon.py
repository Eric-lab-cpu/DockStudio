#!/usr/bin/env python3
"""Generate DockStudio / Eric Studio brand icons (assets/icon.png, icon.ico).

Pure Pillow - no external assets needed.
"""

import os

from PIL import Image, ImageDraw, ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))


def _font(size: int, bold: bool = True):
    candidates = []
    if os.name == "nt":
        candidates += [r"C:\Windows\Fonts\seguisb.ttf", r"C:\Windows\Fonts\arialbd.ttf"]
    candidates += ["/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
                   "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"]
    for c in candidates:
        if os.path.exists(c):
            try:
                return ImageFont.truetype(c, size)
            except Exception:
                pass
    return ImageFont.load_default()


def draw_icon(size: int = 256) -> Image.Image:
    img = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    # rounded background
    r = int(size * 0.22)
    d.rounded_rectangle([0, 0, size - 1, size - 1], radius=r, fill=(21, 80, 158, 255))
    # subtle lighter bottom band
    d.rounded_rectangle([0, int(size * 0.72), size - 1, size - 1], radius=r,
                        fill=(24, 96, 190, 255))
    # white letter D
    font = _font(int(size * 0.52))
    txt = "D"
    bbox = d.textbbox((0, 0), txt, font=font)
    tw = bbox[2] - bbox[0]
    th = bbox[3] - bbox[1]
    d.text(((size - tw) / 2 - bbox[0], size * 0.13 - bbox[1]), txt,
           font=font, fill=(255, 255, 255, 255))
    # small hex "molecule" decoration at bottom
    small = max(3, int(size * 0.03))
    cx, cy = size * 0.5, size * 0.84
    r6 = size * 0.10
    pts = []
    for k in range(6):
        ang = 3.14159 / 3 * k - 3.14159 / 2
        pts.append((cx + r6 * __import__("math").cos(ang),
                    cy + r6 * __import__("math").sin(ang)))
    d.polygon(pts, outline=(190, 220, 255, 255), width=max(2, size // 48))
    # nodes on hex vertices
    for px, py in pts:
        d.ellipse([px - small, py - small, px + small, py + small],
                  fill=(190, 220, 255, 255))
    return img


def main():
    os.makedirs(HERE, exist_ok=True)
    png = os.path.join(HERE, "icon.png")
    ico = os.path.join(HERE, "icon.ico")
    img = draw_icon(256)
    img.save(png)
    # multi-resolution ico
    sizes = [16, 24, 32, 48, 64, 128, 256]
    img.save(ico, sizes=[(s, s) for s in sizes])
    print("wrote", png, ico)


if __name__ == "__main__":
    main()
