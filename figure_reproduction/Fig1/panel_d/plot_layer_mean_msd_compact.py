#!/usr/bin/env python3
"""Reproduce Fig. 1(d): layer-mean MSD at 1500 K and nominal 0 GPa."""

from __future__ import annotations

import argparse
import csv
import html
import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence, Tuple

from PIL import Image, ImageDraw, ImageFont
from matplotlib import font_manager


FIGURE_WIDTH_IN = 3.6
FIGURE_HEIGHT_IN = 1.2
DEFAULT_DPI = 600
SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT_DIR = SCRIPT_DIR / "data"
DEFAULT_OUT_STEM = SCRIPT_DIR / "output" / "layer_mean_msd_3x1p5in"

ARIAL = Path(r"C:\Windows\Fonts\arial.ttf")
if not ARIAL.is_file():
    ARIAL = Path(font_manager.findfont(font_manager.FontProperties(family=["Arial", "DejaVu Sans"])))
FONT_FAMILY = font_manager.FontProperties(fname=str(ARIAL)).get_name()
BLACK = "#202020"
WHITE = "#FFFFFF"

INTERFACE_B = "#08519C"  # deep blue
BULK_B = "#6BAED6"  # light blue
INTERFACE_FE = "#8B1E2D"  # deep red
BULK_FE = "#D9828B"  # light red


@dataclass(frozen=True)
class SeriesSpec:
    filename: str
    label: str
    color: str


@dataclass
class SeriesData:
    spec: SeriesSpec
    time_ps: List[float]
    msd: List[float]


PANELS: Tuple[Tuple[SeriesSpec, ...], ...] = (
    (
        SeriesSpec("FeB_bulk_B_layer_msd.csv", "Bulk B", BULK_B),
        SeriesSpec("FeB_interface_Fe_layer_msd.csv", "Interface Fe", INTERFACE_FE),
        SeriesSpec("FeB_bulk_Fe_layer_msd.csv", "Bulk Fe", BULK_FE),
    ),
    (
        SeriesSpec("B_interface_B_layer_msd.csv", "Interface B", INTERFACE_B),
        SeriesSpec("B_bulk_B_layer_msd.csv", "Bulk B", BULK_B),
    ),
)

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--out-stem", type=Path, default=DEFAULT_OUT_STEM)
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI)
    return parser.parse_args()


def read_series(input_dir: Path, spec: SeriesSpec) -> SeriesData:
    path = input_dir / spec.filename
    time_ps: List[float] = []
    msd: List[float] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            # One LAMMPS timestep is one fs, and 1000 fs = 1 ps.
            time_ps.append(float(row["timestep"]) / 1000.0)
            msd.append(float(row["msd"]))
    if not time_ps:
        raise ValueError(f"no data rows in {path}")
    return SeriesData(spec, time_ps, msd)


def load_panels(input_dir: Path) -> List[List[SeriesData]]:
    return [[read_series(input_dir, spec) for spec in specs] for specs in PANELS]


def nice_y_axis(maximum: float) -> Tuple[float, List[float]]:
    rough_step = maximum / 4.0
    magnitude = 10.0 ** math.floor(math.log10(rough_step))
    fraction = rough_step / magnitude
    if fraction <= 1.0:
        nice_fraction = 1.0
    elif fraction <= 2.0:
        nice_fraction = 2.0
    elif fraction <= 5.0:
        nice_fraction = 5.0
    else:
        nice_fraction = 10.0
    step = nice_fraction * magnitude
    upper = math.ceil(maximum / step) * step
    return upper, [i * step for i in range(int(round(upper / step)) + 1)]


def x_axis(maximum: float) -> Tuple[float, List[float]]:
    upper = math.ceil(maximum / 50.0) * 50.0
    return upper, [value for value in range(0, int(upper) + 1, 200)]


def panel_boxes() -> Tuple[Tuple[float, float, float, float], ...]:
    width_pt = FIGURE_WIDTH_IN * 72.0
    height_pt = FIGURE_HEIGHT_IN * 72.0
    left = 27.0
    right = 3.0
    gap = 8.0
    bottom = 22.0
    top = height_pt - 4.0
    panel_width = (width_pt - left - right - gap) / 2.0
    return (
        (left, bottom, left + panel_width, top),
        (left + panel_width + gap, bottom, width_pt - right, top),
    )


def map_point(
    x: float,
    y: float,
    box: Tuple[float, float, float, float],
    x_max: float,
    y_max: float,
) -> Tuple[float, float]:
    left, bottom, right, top = box
    return (
        left + x / x_max * (right - left),
        bottom + y / y_max * (top - bottom),
    )


def svg_text(
    x: float,
    y: float,
    text: str,
    size: float,
    anchor: str = "middle",
    transform: str = "",
) -> str:
    transform_attr = f' transform="{transform}"' if transform else ""
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}" '
        f'font-family="{FONT_FAMILY}, sans-serif" font-size="{size:.2f}pt" '
        f'fill="{BLACK}"{transform_attr}>{html.escape(text)}</text>'
    )


def render_svg(panels: Sequence[Sequence[SeriesData]], out_path: Path) -> None:
    width_pt = FIGURE_WIDTH_IN * 72.0
    height_pt = FIGURE_HEIGHT_IN * 72.0
    boxes = panel_boxes()
    x_max, x_ticks = x_axis(
        max(max(series.time_ps) for panel in panels for series in panel)
    )
    y_max, y_ticks = nice_y_axis(max(max(series.msd) for panel in panels for series in panel))

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{FIGURE_WIDTH_IN}in" '
        f'height="{FIGURE_HEIGHT_IN}in" viewBox="0 0 {width_pt:.2f} {height_pt:.2f}">',
        f'<rect width="{width_pt:.2f}" height="{height_pt:.2f}" fill="{WHITE}"/>',
        "<defs>",
    ]
    for index, box in enumerate(boxes):
        left, bottom, right, top = box
        parts.append(
            f'<clipPath id="clip-{index}"><rect x="{left:.2f}" y="{height_pt-top:.2f}" '
            f'width="{right-left:.2f}" height="{top-bottom:.2f}"/></clipPath>'
        )
    parts.append("</defs>")

    for panel_index, (panel, box) in enumerate(zip(panels, boxes)):
        left, bottom, right, top = box
        plot_top = height_pt - top
        plot_bottom = height_pt - bottom

        parts.append(
            f'<rect x="{left:.2f}" y="{plot_top:.2f}" width="{right-left:.2f}" '
            f'height="{top-bottom:.2f}" fill="none" stroke="{BLACK}" stroke-width="0.65"/>'
        )

        if panel_index == 0:
            for y_tick in y_ticks:
                _, py = map_point(0.0, y_tick, box, x_max, y_max)
                sy = height_pt - py
                parts.append(
                    f'<line x1="{left-2.5:.2f}" y1="{sy:.2f}" x2="{left:.2f}" y2="{sy:.2f}" '
                    f'stroke="{BLACK}" stroke-width="0.55"/>'
                )
                parts.append(svg_text(left - 3.5, sy + 2.4, f"{y_tick:g}", 7.0, "end"))

        for x_tick in x_ticks:
            px, _ = map_point(x_tick, 0.0, box, x_max, y_max)
            parts.append(
                f'<line x1="{px:.2f}" y1="{plot_bottom:.2f}" x2="{px:.2f}" '
                f'y2="{plot_bottom+2.5:.2f}" stroke="{BLACK}" stroke-width="0.55"/>'
            )
            if math.isclose(x_tick, x_max):
                parts.append(svg_text(px - 1.0, plot_bottom + 9.0, f"{x_tick:g}", 7.0, "end"))
            else:
                parts.append(svg_text(px, plot_bottom + 9.0, f"{x_tick:g}", 7.0))

        draw_order = sorted(
            panel,
            key=lambda series: {
                "Bulk Fe": 0,
                "Bulk B": 1,
                "Interface Fe": 2,
            }.get(series.spec.label, 3),
        )
        for series in draw_order:
            path = " ".join(
                ("M" if index == 0 else "L") + f"{px:.3f},{height_pt-py:.3f}"
                for index, (px, py) in enumerate(
                    map_point(x, y, box, x_max, y_max)
                    for x, y in zip(series.time_ps, series.msd)
                )
            )
            parts.append(
                f'<path d="{path}" clip-path="url(#clip-{panel_index})" fill="none" '
                f'stroke="{series.spec.color}" stroke-width="0.65"/>'
            )

        legend_series = list(panel)
        legend_x = left + 4.0
        legend_y = plot_top + 7.0
        legend_height = len(legend_series) * 8.5 + 2.0
        parts.append(
            f'<rect x="{legend_x-2.0:.2f}" y="{legend_y-6.0:.2f}" width="58" '
            f'height="{legend_height:.2f}" fill="{WHITE}" fill-opacity="0.9"/>'
        )
        for row, series in enumerate(legend_series):
            y = legend_y + row * 8.5
            parts.append(
                f'<line x1="{legend_x:.2f}" y1="{y:.2f}" x2="{legend_x+9.0:.2f}" '
                f'y2="{y:.2f}" stroke="{series.spec.color}" stroke-width="1.15"/>'
            )
            parts.append(svg_text(legend_x + 12.0, y + 2.4, series.spec.label, 7.0, "start"))

    y_label_x = 7.2
    y_label_y = height_pt - sum((boxes[0][1], boxes[0][3])) / 2.0
    parts.append(
        svg_text(
            y_label_x,
            y_label_y,
            "MSD (Å²)",
            7.5,
            transform=f"rotate(-90 {y_label_x:.2f} {y_label_y:.2f})",
        )
    )
    parts.append(svg_text(width_pt / 2.0, height_pt - 2.0, "Time (ps)", 7.5))

    parts.append("</svg>")
    out_path.write_text("\n".join(parts), encoding="utf-8")


def load_font(points: float, scale: float) -> ImageFont.FreeTypeFont:
    return ImageFont.truetype(str(ARIAL), max(1, round(points * scale)))


def draw_rotated_label(
    image: Image.Image,
    center: Tuple[int, int],
    value: str,
    text_font: ImageFont.FreeTypeFont,
) -> None:
    bbox = text_font.getbbox(value)
    width = bbox[2] - bbox[0] + 12
    height = bbox[3] - bbox[1] + 12
    label = Image.new("RGBA", (width, height), (255, 255, 255, 0))
    label_draw = ImageDraw.Draw(label)
    label_draw.text((6 - bbox[0], 6 - bbox[1]), value, font=text_font, fill=BLACK)
    rotated = label.rotate(90, expand=True, resample=Image.Resampling.BICUBIC)
    image.alpha_composite(rotated, (center[0] - rotated.width // 2, center[1] - rotated.height // 2))


def render_png(
    panels: Sequence[Sequence[SeriesData]],
    out_path: Path,
    dpi: int,
) -> None:
    width_px = round(FIGURE_WIDTH_IN * dpi)
    height_px = round(FIGURE_HEIGHT_IN * dpi)
    scale = dpi / 72.0
    boxes = panel_boxes()
    x_max, x_ticks = x_axis(
        max(max(series.time_ps) for panel in panels for series in panel)
    )
    y_max, y_ticks = nice_y_axis(max(max(series.msd) for panel in panels for series in panel))

    image = Image.new("RGBA", (width_px, height_px), WHITE)
    draw = ImageDraw.Draw(image)
    tick_font = load_font(7.0, scale)
    axis_font = load_font(7.5, scale)
    legend_font = load_font(7.0, scale)

    def pixel(point: Tuple[float, float]) -> Tuple[float, float]:
        return point[0] * scale, height_px - point[1] * scale

    for panel_index, (panel, box) in enumerate(zip(panels, boxes)):
        left, bottom, right, top = box
        left_px = left * scale
        right_px = right * scale
        top_px = height_px - top * scale
        bottom_px = height_px - bottom * scale
        draw.rectangle(
            (left_px, top_px, right_px, bottom_px),
            outline=BLACK,
            width=max(1, round(0.65 * scale)),
        )

        if panel_index == 0:
            for y_tick in y_ticks:
                _, py = map_point(0.0, y_tick, box, x_max, y_max)
                y_px = height_px - py * scale
                draw.line(
                    ((left - 2.5) * scale, y_px, left_px, y_px),
                    fill=BLACK,
                    width=max(1, round(0.55 * scale)),
                )
                draw.text(
                    ((left - 3.5) * scale, y_px),
                    f"{y_tick:g}",
                    font=tick_font,
                    fill=BLACK,
                    anchor="rm",
                )

        for x_tick in x_ticks:
            x_pt, _ = map_point(x_tick, 0.0, box, x_max, y_max)
            x_px = x_pt * scale
            draw.line(
                (x_px, bottom_px, x_px, bottom_px + 2.5 * scale),
                fill=BLACK,
                width=max(1, round(0.55 * scale)),
            )
            if math.isclose(x_tick, x_max):
                draw.text(
                    (x_px - 1.0 * scale, bottom_px + 6.5 * scale),
                    f"{x_tick:g}",
                    font=tick_font,
                    fill=BLACK,
                    anchor="rm",
                )
            else:
                draw.text(
                    (x_px, bottom_px + 6.5 * scale),
                    f"{x_tick:g}",
                    font=tick_font,
                    fill=BLACK,
                    anchor="mm",
                )

        draw_order = sorted(
            panel,
            key=lambda series: {
                "Bulk Fe": 0,
                "Bulk B": 1,
                "Interface Fe": 2,
            }.get(series.spec.label, 3),
        )
        for series in draw_order:
            points = [
                pixel(map_point(x, y, box, x_max, y_max))
                for x, y in zip(series.time_ps, series.msd)
            ]
            draw.line(
                points,
                fill=series.spec.color,
                width=max(1, round(0.65 * scale)),
                joint="curve",
            )

        legend_series = list(panel)
        legend_x = (left + 4.0) * scale
        legend_y = top_px + 7.0 * scale
        legend_height = (len(legend_series) * 8.5 + 2.0) * scale
        draw.rectangle(
            (
                legend_x - 2.0 * scale,
                legend_y - 6.0 * scale,
                legend_x + 56.0 * scale,
                legend_y - 6.0 * scale + legend_height,
            ),
            fill=WHITE,
        )
        for row, series in enumerate(legend_series):
            y = legend_y + row * 8.5 * scale
            draw.line(
                (legend_x, y, legend_x + 9.0 * scale, y),
                fill=series.spec.color,
                width=max(1, round(1.15 * scale)),
            )
            draw.text(
                (legend_x + 12.0 * scale, y),
                series.spec.label,
                font=legend_font,
                fill=BLACK,
                anchor="lm",
            )

    draw_rotated_label(
        image,
        (
            round(7.2 * scale),
            round(height_px - ((boxes[0][1] + boxes[0][3]) / 2.0) * scale),
        ),
        "MSD (Å²)",
        axis_font,
    )
    draw.text(
        (width_px / 2.0, height_px - 3.8 * scale),
        "Time (ps)",
        font=axis_font,
        fill=BLACK,
        anchor="mm",
    )

    image.convert("RGB").save(out_path, dpi=(dpi, dpi), optimize=True)


def render_pdf(
    panels: Sequence[Sequence[SeriesData]],
    out_path: Path,
    png_path: Path,
    dpi: int,
) -> None:
    try:
        from reportlab.pdfbase import pdfmetrics
        from reportlab.pdfbase.ttfonts import TTFont
        from reportlab.pdfgen import canvas
    except ModuleNotFoundError:
        with Image.open(png_path) as image:
            image.convert("RGB").save(out_path, "PDF", resolution=dpi)
        return

    width_pt = FIGURE_WIDTH_IN * 72.0
    height_pt = FIGURE_HEIGHT_IN * 72.0
    boxes = panel_boxes()
    x_max, x_ticks = x_axis(
        max(max(series.time_ps) for panel in panels for series in panel)
    )
    y_max, y_ticks = nice_y_axis(max(max(series.msd) for panel in panels for series in panel))

    font_name = "ArialMSD"
    if font_name not in pdfmetrics.getRegisteredFontNames():
        pdfmetrics.registerFont(TTFont(font_name, str(ARIAL)))

    pdf = canvas.Canvas(str(out_path), pagesize=(width_pt, height_pt), pageCompression=1)
    pdf.setFillColor(WHITE)
    pdf.rect(0, 0, width_pt, height_pt, stroke=0, fill=1)

    for panel_index, (panel, box) in enumerate(zip(panels, boxes)):
        left, bottom, right, top = box
        pdf.setStrokeColor(BLACK)
        pdf.setLineWidth(0.65)
        pdf.rect(left, bottom, right - left, top - bottom, stroke=1, fill=0)

        pdf.setFillColor(BLACK)
        pdf.setFont(font_name, 7.0)
        if panel_index == 0:
            for y_tick in y_ticks:
                _, py = map_point(0.0, y_tick, box, x_max, y_max)
                pdf.setLineWidth(0.55)
                pdf.line(left - 2.5, py, left, py)
                pdf.drawRightString(left - 3.5, py - 2.4, f"{y_tick:g}")

        for x_tick in x_ticks:
            px, _ = map_point(x_tick, 0.0, box, x_max, y_max)
            pdf.setLineWidth(0.55)
            pdf.line(px, bottom, px, bottom - 2.5)
            if math.isclose(x_tick, x_max):
                pdf.drawRightString(px - 1.0, bottom - 9.0, f"{x_tick:g}")
            else:
                label = f"{x_tick:g}"
                label_width = pdfmetrics.stringWidth(label, font_name, 7.0)
                pdf.drawString(px - label_width / 2.0, bottom - 9.0, label)

        pdf.saveState()
        clip = pdf.beginPath()
        clip.rect(left, bottom, right - left, top - bottom)
        pdf.clipPath(clip, stroke=0, fill=0)
        draw_order = sorted(
            panel,
            key=lambda series: {
                "Bulk Fe": 0,
                "Bulk B": 1,
                "Interface Fe": 2,
            }.get(series.spec.label, 3),
        )
        for series in draw_order:
            path = pdf.beginPath()
            for index, (x, y) in enumerate(zip(series.time_ps, series.msd)):
                px, py = map_point(x, y, box, x_max, y_max)
                if index == 0:
                    path.moveTo(px, py)
                else:
                    path.lineTo(px, py)
            pdf.setStrokeColor(series.spec.color)
            pdf.setLineWidth(0.65)
            pdf.drawPath(path, stroke=1, fill=0)
        pdf.restoreState()

        legend_series = list(panel)
        legend_x = left + 4.0
        first_y = top - 7.0
        legend_height = len(legend_series) * 8.5 + 2.0
        pdf.setFillColor(WHITE)
        pdf.rect(
            legend_x - 2.0,
            first_y - (len(legend_series) - 1) * 8.5 - 5.0,
            58.0,
            legend_height,
            stroke=0,
            fill=1,
        )
        pdf.setFont(font_name, 7.0)
        for row, series in enumerate(legend_series):
            y = first_y - row * 8.5
            pdf.setStrokeColor(series.spec.color)
            pdf.setLineWidth(1.15)
            pdf.line(legend_x, y, legend_x + 9.0, y)
            pdf.setFillColor(BLACK)
            pdf.drawString(legend_x + 12.0, y - 2.4, series.spec.label)

    pdf.setFillColor(BLACK)
    pdf.setFont(font_name, 7.5)
    pdf.saveState()
    pdf.translate(7.2, (boxes[0][1] + boxes[0][3]) / 2.0)
    pdf.rotate(90)
    y_label = "MSD (Å²)"
    y_label_width = pdfmetrics.stringWidth(y_label, font_name, 7.5)
    pdf.drawString(-y_label_width / 2.0, 0.0, y_label)
    pdf.restoreState()

    x_label = "Time (ps)"
    x_label_width = pdfmetrics.stringWidth(x_label, font_name, 7.5)
    pdf.drawString((width_pt - x_label_width) / 2.0, 3.0, x_label)
    pdf.showPage()
    pdf.save()


def render_eps(png_path: Path, out_path: Path) -> None:
    """Write an EPS file using the already-rendered figure."""
    with Image.open(png_path) as image:
        image.convert("RGB").save(out_path, format="EPS", dpi=(DEFAULT_DPI, DEFAULT_DPI))


def render_matplotlib_pdf(
    panels: Sequence[Sequence[SeriesData]], out_path: Path
) -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": FONT_FAMILY, "font.size": 7})
    fig, axes = plt.subplots(1, 2, figsize=(FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN))
    x_max, x_ticks = x_axis(max(max(series.time_ps) for panel in panels for series in panel))
    y_max, y_ticks = nice_y_axis(max(max(series.msd) for panel in panels for series in panel))
    draw_order = {"Bulk Fe": 0, "Bulk B": 1, "Interface Fe": 2}

    for index, (axis, panel) in enumerate(zip(axes, panels)):
        for series in sorted(panel, key=lambda item: draw_order.get(item.spec.label, 3)):
            axis.plot(
                series.time_ps,
                series.msd,
                color=series.spec.color,
                linewidth=0.65,
                label=series.spec.label,
            )
        axis.set_xlim(0, x_max)
        axis.set_ylim(0, y_max)
        axis.set_xticks(x_ticks)
        axis.set_yticks(y_ticks)
        axis.tick_params(direction="out", length=2.5, width=0.55, pad=2)
        for spine in axis.spines.values():
            spine.set_linewidth(0.65)
            spine.set_color(BLACK)
        if index == 0:
            axis.set_ylabel("MSD (Å²)", labelpad=4)
        else:
            axis.set_yticklabels([])
        axis.set_xlabel("Time (ps)" if index == 0 else "")
        axis.legend(
            loc="upper left",
            frameon=True,
            facecolor="white",
            framealpha=0.9,
            edgecolor="none",
            handlelength=1.5,
            handletextpad=0.5,
            borderpad=0.3,
            labelspacing=0.5,
            fontsize=7,
        )

    fig.subplots_adjust(left=0.16, right=0.98, bottom=0.22, top=0.98, wspace=0.22)
    fig.savefig(out_path, dpi=600, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    if args.dpi <= 0:
        raise ValueError("--dpi must be positive")
    panels = load_panels(args.input_dir)
    args.out_stem.parent.mkdir(parents=True, exist_ok=True)
    png_path = args.out_stem.with_suffix(".png")
    svg_path = args.out_stem.with_suffix(".svg")
    pdf_path = args.out_stem.with_suffix(".pdf")
    eps_path = args.out_stem.with_suffix(".eps")
    render_png(panels, png_path, args.dpi)
    render_svg(panels, svg_path)
    render_matplotlib_pdf(panels, pdf_path)
    render_eps(png_path, eps_path)
    print(f"Wrote {png_path}")
    print(f"Wrote {svg_path}")
    print(f"Wrote {pdf_path}")
    print(f"Wrote {eps_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
