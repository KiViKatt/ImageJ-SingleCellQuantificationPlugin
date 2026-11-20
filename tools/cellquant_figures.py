# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Ivan Khimach

"""Cell Quant Studio Universal - figures from result tables.

Turns the CSV tables written by Cell Quant Studio into SVG figures for a
manuscript or a lab meeting. It runs with any Python 3 (3.6 or newer) and needs
no extra packages: open the SVG files in a web browser, Inkscape, Illustrator,
PowerPoint or Word. Each figure is written together with a CSV of exactly the
numbers it shows.

This is a desktop tool, not part of the Fiji application: run it with your
normal Python, not with Fiji's Jython.

Usage (from a terminal):

    python tools/cellquant_figures.py validation  <folder or validation_summary.csv> -o figures
    python tools/cellquant_figures.py depth       <folder or depth_profile.csv>      -o figures
    python tools/cellquant_figures.py spectrum    <folder or threshold_spectrum.csv> -o figures
    python tools/cellquant_figures.py robustness  <CellQuant_Robustness_... folder>  -o figures

A folder is searched (with its subfolders) for the table the figure needs, so
you can point it at a batch folder or at a folder holding many run folders.

Groups: images are coloured by group. By default the group is the start of the
image file name up to the first "-", "_" or space (so "Control-org1..." is
group "Control" and "Mutant-org1..." is group "Mutant"). Change the rule with
--group-regex; its first bracketed part is the group, for example
--group-regex "^(Control|Mutant)".

The source is kept plain ASCII with no encoding declaration on purpose.
"""

import argparse
import csv
import math
import os
import re
import sys

# -----------------------------------------------------------------------------
# Colours (validated default categorical palette, light mode) and text inks
# -----------------------------------------------------------------------------

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SURFACE = "#ffffff"
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_3 = "#8a8983"
GRID = "#e6e5e0"
FONT = "Helvetica, Arial, sans-serif"

DEFAULT_GROUP_REGEX = r"^([^-_ ]+)"


# -----------------------------------------------------------------------------
# Reading tables
# -----------------------------------------------------------------------------

def find_tables(path, names):
    """Return the CSV files called one of names, path being a file or a folder."""
    if os.path.isfile(path):
        return [path]
    found = []
    for root, dirs, files in os.walk(path):
        dirs.sort()
        for name in sorted(files):
            if name in names:
                found.append(os.path.join(root, name))
    if not found:
        raise SystemExit("No " + " or ".join(names) + " found in " + path)
    # A batch folder holds the combined batch_<name> table and, in each run
    # folder inside it, the same rows again. Use the combined table only.
    batch_dirs = [os.path.dirname(p) for p in found if os.path.basename(p).startswith("batch_")]
    kept = []
    for p in found:
        if not os.path.basename(p).startswith("batch_"):
            inside = [d for d in batch_dirs if os.path.abspath(p).startswith(os.path.abspath(d) + os.sep)]
            if inside:
                continue
        kept.append(p)
    return kept


def read_rows(paths):
    rows = []
    for path in paths:
        with open(path, newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                row["_file"] = path
                rows.append(row)
    return rows


def number(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(result) or math.isinf(result):
        return None
    return result


def group_of(text, regex):
    match = re.search(regex, text or "")
    if match and match.groups():
        return match.group(1)
    if match:
        return match.group(0)
    return "all"


def group_colours(groups):
    ordered = sorted(set(groups))
    if len(ordered) > len(SERIES):
        raise SystemExit(
            "More than %d groups (%s). Use --group-regex to define fewer groups." % (len(SERIES), ", ".join(ordered))
        )
    return dict((name, SERIES[index]) for index, name in enumerate(ordered))


def write_table(path, headers, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        for row in rows:
            writer.writerow([row.get(h, "") for h in headers])


# -----------------------------------------------------------------------------
# A very small SVG plotting layer
# -----------------------------------------------------------------------------

def esc(text):
    return (str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;"))


def nice_ticks(low, high, count=5):
    if low == high:
        low, high = low - 1.0, high + 1.0
    span = high - low
    raw = span / max(1, count)
    power = 10 ** math.floor(math.log10(raw))
    for step in [1, 2, 2.5, 5, 10]:
        if raw <= step * power:
            step = step * power
            break
    # Whole multiples of step from the one at or below low to the one at or
    # above high (counted as integers so rounding cannot drop the last tick).
    first = int(math.floor(low / step + 1e-9))
    last = int(math.ceil(high / step - 1e-9))
    return [round(k * step, 10) for k in range(first, last + 1)]


def tick_text(value):
    if value == 0:
        return "0"
    if abs(value) >= 10000 or abs(value) < 0.01:
        return "%.3g" % value
    if float(value).is_integer():
        return "%d" % value
    return ("%.3f" % value).rstrip("0").rstrip(".")


class Panel(object):
    """One plotting area with linear x and y scales inside an Svg."""

    def __init__(self, svg, x, y, width, height, xlim, ylim, xlabel="", ylabel="", title=""):
        self.svg = svg
        self.x0, self.y0, self.w, self.h = x, y, width, height
        self.xt = nice_ticks(xlim[0], xlim[1])
        self.yt = nice_ticks(ylim[0], ylim[1])
        self.xmin, self.xmax = min(self.xt[0], xlim[0]), max(self.xt[-1], xlim[1])
        self.ymin, self.ymax = min(self.yt[0], ylim[0]), max(self.yt[-1], ylim[1])
        self.xlabel, self.ylabel, self.title = xlabel, ylabel, title

    def sx(self, value):
        return self.x0 + (value - self.xmin) / float(self.xmax - self.xmin or 1) * self.w

    def sy(self, value):
        return self.y0 + self.h - (value - self.ymin) / float(self.ymax - self.ymin or 1) * self.h

    def axes(self, show_x_ticks=True):
        s = self.svg
        for t in self.yt:
            y = self.sy(t)
            s.line(self.x0, y, self.x0 + self.w, y, GRID, 1)
            s.text(self.x0 - 8, y + 4, tick_text(t), 11, INK_2, "end")
        if show_x_ticks:
            for t in self.xt:
                x = self.sx(t)
                s.line(x, self.y0 + self.h, x, self.y0 + self.h + 4, INK_3, 1)
                s.text(x, self.y0 + self.h + 18, tick_text(t), 11, INK_2, "middle")
        s.line(self.x0, self.y0 + self.h, self.x0 + self.w, self.y0 + self.h, INK_3, 1)
        if self.xlabel:
            s.text(self.x0 + self.w / 2.0, self.y0 + self.h + 38, self.xlabel, 12, INK_2, "middle")
        if self.ylabel:
            s.text(self.x0 - 52, self.y0 + self.h / 2.0, self.ylabel, 12, INK_2, "middle",
                   rotate=-90)
        if self.title:
            s.text(self.x0, self.y0 - 12, self.title, 13, INK, "start", bold=True)

    def polyline(self, points, colour, width=2, opacity=1.0, dash=None, tip=None):
        pts = [(self.sx(x), self.sy(y)) for x, y in points if x is not None and y is not None]
        if len(pts) >= 2:
            self.svg.polyline(pts, colour, width, opacity, dash, tip)

    def dot(self, x, y, colour, radius=4, tip=None, hollow=False):
        if x is None or y is None:
            return
        self.svg.circle(self.sx(x), self.sy(y), radius, colour, tip, hollow)

    def hline(self, y, colour=INK_3, dash="4 3", label=None):
        if y is None:
            return
        self.svg.line(self.x0, self.sy(y), self.x0 + self.w, self.sy(y), colour, 1, dash)
        if label:
            self.svg.text(self.x0 + self.w - 4, self.sy(y) - 5, label, 11, INK_2, "end")


class Svg(object):
    def __init__(self, width, height, title, description=""):
        self.width, self.height = width, height
        self.parts = []
        self.title, self.description = title, description

    def line(self, x1, y1, x2, y2, colour, width=1, dash=None):
        extra = ' stroke-dasharray="%s"' % dash if dash else ""
        self.parts.append('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="%s"%s/>'
                          % (x1, y1, x2, y2, colour, width, extra))

    def polyline(self, pts, colour, width, opacity, dash, tip):
        extra = ' stroke-dasharray="%s"' % dash if dash else ""
        body = " ".join("%.2f,%.2f" % p for p in pts)
        title = "<title>%s</title>" % esc(tip) if tip else ""
        self.parts.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="%s" stroke-opacity="%s" '
                          'stroke-linejoin="round" stroke-linecap="round"%s>%s</polyline>'
                          % (body, colour, width, opacity, extra, title))

    def circle(self, x, y, radius, colour, tip=None, hollow=False):
        title = "<title>%s</title>" % esc(tip) if tip else ""
        fill = SURFACE if hollow else colour
        self.parts.append('<circle cx="%.2f" cy="%.2f" r="%s" fill="%s" stroke="%s" stroke-width="2">%s</circle>'
                          % (x, y, radius, fill, colour if hollow else SURFACE, title))

    def rect(self, x, y, width, height, colour):
        self.parts.append('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" rx="2" fill="%s"/>'
                          % (x, y, width, height, colour))

    def text(self, x, y, text, size, colour, anchor="start", bold=False, rotate=None):
        weight = ' font-weight="600"' if bold else ""
        transform = ' transform="rotate(%s %.2f %.2f)"' % (rotate, x, y) if rotate else ""
        self.parts.append('<text x="%.2f" y="%.2f" font-size="%s" fill="%s" text-anchor="%s"%s%s>%s</text>'
                          % (x, y, size, colour, anchor, weight, transform, esc(text)))

    def legend(self, x, y, items):
        """items: list of (label, colour). Drawn in one row."""
        for label, colour in items:
            self.rect(x, y - 9, 14, 10, colour)
            self.text(x + 19, y, label, 12, INK)
            x += 30 + 7 * len(label)

    def save(self, path):
        head = ('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d" '
                'font-family="%s" role="img">' % (self.width, self.height, self.width, self.height, FONT))
        meta = "<title>%s</title><desc>%s</desc>" % (esc(self.title), esc(self.description))
        background = '<rect width="100%%" height="100%%" fill="%s"/>' % SURFACE
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(head + meta + background + "".join(self.parts) + "</svg>\n")


def limits(values, pad=0.05, include_zero=False, floor=None):
    """Axis range around values. floor: a value the range never goes below
    (0 for counts and percentages, which cannot be negative)."""
    values = [v for v in values if v is not None]
    if not values:
        return (0.0, 1.0)
    low, high = min(values), max(values)
    if include_zero:
        low, high = min(low, 0.0), max(high, 0.0)
    span = (high - low) or abs(high) or 1.0
    low, high = low - span * pad, high + span * pad
    if floor is not None and low < floor:
        low = floor
    return (low, high)


# -----------------------------------------------------------------------------
# Validation: detected against manual counts, and Bland-Altman
# -----------------------------------------------------------------------------

def figure_validation(args):
    rows = read_rows(find_tables(args.input, ["validation_summary.csv"]))
    points = []
    for row in rows:
        if row.get("Parent_ROI_Name") == "ALL_ROIS_POOLED":
            continue
        manual, detected = number(row.get("Manual_Count")), number(row.get("Detected_Count"))
        if manual is None or detected is None:
            continue
        points.append({
            "Source_Image": row.get("Source_Image", ""), "Parent_ROI_Name": row.get("Parent_ROI_Name", ""),
            "Channel_Name": row.get("Channel_Name", ""), "Manual_Count": manual, "Detected_Count": detected,
            "Mean_Count": (manual + detected) / 2.0, "Difference": detected - manual,
            "F1": row.get("F1", ""), "Group": group_of(row.get("Source_Image", ""), args.group_regex)
        })
    if len(points) < 2:
        raise SystemExit("Need at least two validated fields (rows of validation_summary.csv); found %d." % len(points))
    diffs = [p["Difference"] for p in points]
    n = len(diffs)
    bias = sum(diffs) / n
    sd = math.sqrt(sum((d - bias) ** 2 for d in diffs) / (n - 1))
    low, high = bias - 1.96 * sd, bias + 1.96 * sd
    colours = group_colours([p["Group"] for p in points])

    svg = Svg(900, 460, args.title or "Detected against manual counts",
              "Left: detected count against manual count per validated field, with the identity line. "
              "Right: Bland-Altman plot of detected minus manual against their mean, with bias and 95%% limits "
              "of agreement (n = %d fields)." % n)
    lim = limits([p["Manual_Count"] for p in points] + [p["Detected_Count"] for p in points], include_zero=True, floor=0.0)
    left = Panel(svg, 90, 70, 320, 300, lim, lim, "Manual count", "Detected count", "A  Agreement")
    left.axes()
    left.polyline([(lim[0], lim[0]), (lim[1], lim[1])], INK_3, 1, 1.0, "4 3")
    right = Panel(svg, 550, 70, 320, 300, limits([p["Mean_Count"] for p in points]),
                  limits(diffs + [low, high], include_zero=True),
                  "Mean of detected and manual", "Detected minus manual", "B  Bland-Altman")
    right.axes()
    right.hline(0.0, GRID, None)
    right.hline(bias, INK_2, "6 3", "bias %s" % tick_text(round(bias, 2)))
    right.hline(low, INK_3, "3 3", "-1.96 SD %s" % tick_text(round(low, 2)))
    right.hline(high, INK_3, "3 3", "+1.96 SD %s" % tick_text(round(high, 2)))
    for p in points:
        tip = "%s / %s: manual %d, detected %d" % (p["Source_Image"], p["Parent_ROI_Name"], p["Manual_Count"], p["Detected_Count"])
        left.dot(p["Manual_Count"], p["Detected_Count"], colours[p["Group"]], 4, tip)
        right.dot(p["Mean_Count"], p["Difference"], colours[p["Group"]], 4, tip)
    if len(colours) > 1:
        svg.legend(90, 30, sorted(colours.items()))
    svg.text(90, 440, "n = %d fields; bias %.2f; limits of agreement %.2f to %.2f" % (n, bias, low, high), 12, INK_2)
    base = os.path.join(args.output, "validation_agreement")
    svg.save(base + ".svg")
    write_table(base + ".csv", ["Group", "Source_Image", "Parent_ROI_Name", "Channel_Name", "Manual_Count",
                                "Detected_Count", "Mean_Count", "Difference", "F1"], points)
    print("Bias %.3f, SD %.3f, limits of agreement %.3f to %.3f, n = %d" % (bias, sd, low, high, n))
    return [base + ".svg"]


# -----------------------------------------------------------------------------
# Depth profile
# -----------------------------------------------------------------------------

def depth_metric(rows, wanted):
    """The measure to plot for one channel: the one asked for when this
    channel has it, otherwise the most informative one it has."""
    if wanted and any(number(r.get(wanted)) is not None for r in rows):
        return wanted
    for name in ["Object_Density_Per_mm2", "Object_Density_Per_Area", "Positive_Area_Fraction_Percent", "Mean_Intensity"]:
        if any(number(r.get(name)) is not None for r in rows):
            return name
    return "Mean_Intensity"


METRIC_LABELS = {
    "Object_Density_Per_mm2": "Objects per mm2",
    "Object_Density_Per_Area": "Objects per unit area",
    "Positive_Area_Fraction_Percent": "Positive area (%)",
    "Mean_Intensity": "Mean intensity",
    "Object_Count": "Objects",
}


def figure_depth(args):
    rows = read_rows(find_tables(args.input, ["depth_profile.csv", "batch_depth_profile.csv"]))
    rows = [r for r in rows if number(r.get("Bin_Index")) is not None]
    if args.frame:
        rows = [r for r in rows if r.get("Frame") == str(args.frame)]
    channels = []
    for r in rows:
        key = (r.get("Channel", ""), r.get("Channel_Name", ""))
        if key not in channels:
            channels.append(key)
    if not channels:
        raise SystemExit("The depth tables hold no bins.")
    unit = next((r.get("Depth_Unit") for r in rows if r.get("Depth_Unit")), "")
    panel_h, gap = 240, 110
    svg = Svg(820, 90 + len(channels) * (panel_h + gap), args.title or "Depth profile from the outline",
              "Per channel: the chosen measure against depth from the parent ROI outline. Thin lines: one per "
              "image and parent ROI. Thick lines: group means over those profiles.")
    plotted = []
    all_groups = set()
    files = []
    for index, (channel, name) in enumerate(channels):
        crow = [r for r in rows if r.get("Channel") == channel]
        metric = depth_metric(crow, args.metric)
        profiles = {}
        for r in crow:
            x = (number(r.get("Relative_Depth_Start")), number(r.get("Relative_Depth_End"))) if args.relative else None
            if args.relative:
                x = (x[0] + x[1]) / 2.0 if None not in x else None
            else:
                x = number(r.get("Depth_Mid"))
            y = number(r.get(metric))
            if x is None or y is None:
                continue
            key = (r.get("Source_Image", ""), r.get("Parent_ROI_Name", ""), r.get("Frame", ""))
            profiles.setdefault(key, []).append((int(float(r["Bin_Index"])), x, y))
            plotted.append({"Channel": channel, "Channel_Name": name, "Source_Image": key[0], "Parent_ROI_Name": key[1],
                            "Frame": key[2], "Group": group_of(key[0], args.group_regex), "Depth": x,
                            "Measure": metric, "Value": y})
        for key in profiles:
            all_groups.add(group_of(key[0], args.group_regex))
        files.append((index, channel, name, metric, profiles))
    colours = group_colours(all_groups or ["all"])
    for index, channel, name, metric, profiles in files:
        xs = [p[1] for prof in profiles.values() for p in prof]
        ys = [p[2] for prof in profiles.values() for p in prof]
        panel = Panel(svg, 100, 80 + index * (panel_h + gap), 680, panel_h, (0.0, limits(xs, 0.02)[1]),
                      limits(ys, include_zero=True, floor=0.0 if ys and min(ys) >= 0 else None),
                      "Relative depth (0 = outline, 1 = deepest point)" if args.relative else "Depth from the outline (%s)" % unit,
                      METRIC_LABELS.get(metric, metric), "C%s %s" % (channel, name))
        panel.axes()
        sums = {}
        for key, prof in sorted(profiles.items()):
            group = group_of(key[0], args.group_regex)
            prof.sort()
            panel.polyline([(p[1], p[2]) for p in prof], colours[group], 1, 0.35, None, "%s / %s" % (key[0], key[1]))
            for bin_index, x, y in prof:
                entry = sums.setdefault((group, bin_index), [0.0, 0.0, 0])
                entry[0] += x
                entry[1] += y
                entry[2] += 1
        for group in sorted(colours):
            means = sorted((b, e[0] / e[2], e[1] / e[2]) for (g, b), e in sums.items() if g == group)
            panel.polyline([(m[1], m[2]) for m in means], colours[group], 2.5, 1.0, None, "%s mean" % group)
    if len(colours) > 1:
        svg.legend(100, 36, sorted(colours.items()))
    base = os.path.join(args.output, "depth_profile" + ("_relative" if args.relative else ""))
    svg.save(base + ".svg")
    write_table(base + ".csv", ["Group", "Source_Image", "Parent_ROI_Name", "Frame", "Channel", "Channel_Name",
                                "Depth", "Measure", "Value"], plotted)
    return [base + ".svg"]


# -----------------------------------------------------------------------------
# Threshold spectrum
# -----------------------------------------------------------------------------

def figure_spectrum(args):
    rows = read_rows(find_tables(args.input, ["threshold_spectrum.csv", "batch_threshold_spectrum.csv"]))
    if args.frame:
        rows = [r for r in rows if r.get("Frame") == str(args.frame)]
    channels = []
    for r in rows:
        key = (r.get("Channel", ""), r.get("Channel_Name", ""))
        if key not in channels:
            channels.append(key)
    if not channels:
        raise SystemExit("The spectrum tables hold no rows.")
    panel_h, gap = 240, 110
    svg = Svg(820, 90 + len(channels) * (panel_h + gap), args.title or "Coverage across thresholds",
              "Per channel: percent of ROI pixels at or above each threshold (at or below for light-background "
              "channels). One line per image and parent ROI; the dot marks the threshold used in the run.")
    groups = set(group_of(r.get("Source_Image", ""), args.group_regex) for r in rows)
    colours = group_colours(groups)
    plotted = []
    for index, (channel, name) in enumerate(channels):
        crow = [r for r in rows if r.get("Channel") == channel]
        curves, configured = {}, {}
        for r in crow:
            key = (r.get("Source_Image", ""), r.get("Parent_ROI_Name", ""), r.get("Frame", ""))
            x, y = number(r.get("Threshold")), number(r.get("Coverage_Percent"))
            if x is None or y is None:
                continue
            if r.get("Level_Index") == "0":
                configured[key] = (x, y)
            else:
                curves.setdefault(key, []).append((x, y))
            plotted.append({"Channel": channel, "Channel_Name": name, "Source_Image": key[0], "Parent_ROI_Name": key[1],
                            "Frame": key[2], "Group": group_of(key[0], args.group_regex),
                            "Level_Index": r.get("Level_Index"), "Threshold": x, "Coverage_Percent": y})
        xs = [p[0] for c in curves.values() for p in c] + [v[0] for v in configured.values()]
        panel = Panel(svg, 100, 80 + index * (panel_h + gap), 680, panel_h, limits(xs, 0.02), (0.0, 100.0),
                      "Threshold (prepared-image intensity)", "Coverage (%)", "C%s %s" % (channel, name))
        panel.axes()
        for key, curve in sorted(curves.items()):
            group = group_of(key[0], args.group_regex)
            panel.polyline(sorted(curve), colours[group], 1.5, 0.6, None, "%s / %s" % (key[0], key[1]))
        for key, (x, y) in sorted(configured.items()):
            group = group_of(key[0], args.group_regex)
            panel.dot(x, y, colours[group], 4, "%s / %s: threshold %s, coverage %.2f%%" % (key[0], key[1], tick_text(x), y))
    if len(colours) > 1:
        svg.legend(100, 36, sorted(colours.items()))
    base = os.path.join(args.output, "threshold_spectrum")
    svg.save(base + ".svg")
    write_table(base + ".csv", ["Group", "Source_Image", "Parent_ROI_Name", "Frame", "Channel", "Channel_Name",
                                "Level_Index", "Threshold", "Coverage_Percent"], plotted)
    return [base + ".svg"]


# -----------------------------------------------------------------------------
# Robustness: specification curve
# -----------------------------------------------------------------------------

def figure_robustness(args):
    spec_path = find_tables(args.input, ["robustness_specifications.csv"])
    if len(spec_path) > 1:
        raise SystemExit("Several robustness runs found; point at one CellQuant_Robustness_... folder.")
    specs = read_rows(spec_path)
    folder = os.path.dirname(spec_path[0])
    measure = "Hedges_g" if args.metric in ("g", "Hedges_g") else "Difference_B_Minus_A"
    rows = []
    for r in specs:
        y = number(r.get(measure))
        rank = number(r.get("Curve_Rank"))
        if y is None or rank is None:
            continue
        rows.append(dict(r, _y=y, _rank=rank))
    if not rows:
        raise SystemExit("No specification has a value for " + measure + ".")
    rows.sort(key=lambda r: r["_rank"])
    for position, r in enumerate(rows):
        r["_x"] = position + 1
    parameters = []
    parameter_rows = []
    parameter_file = os.path.join(folder, "robustness_specification_parameters.csv")
    if os.path.isfile(parameter_file):
        parameter_rows = [p for p in read_rows([parameter_file]) if p.get("Applies", "true").lower() in ("true", "1", "yes")]
        for p in parameter_rows:
            key = (p.get("Parameter"), p.get("Value"))
            if key not in parameters:
                parameters.append(key)
        parameters.sort(key=lambda k: (k[0], number(k[1]) if number(k[1]) is not None else 0, k[1]))
    first = rows[0]
    group_a, group_b = first.get("Group_A", "A"), first.get("Group_B", "B")
    row_h = 16
    matrix_h = row_h * len(parameters)
    svg = Svg(900, 150 + 300 + (matrix_h + 40 if parameters else 0),
              args.title or "Specification curve: %s minus %s" % (group_b, group_a),
              "Each column is one combination of analysis settings, sorted by effect. Top: the group difference "
              "(%s) with its 95%% confidence interval; the configured settings are marked. Bottom: which setting "
              "values each combination used." % measure)
    ys = [r["_y"] for r in rows]
    if measure == "Difference_B_Minus_A":
        ys += [number(r.get("Difference_CI95_Low")) for r in rows] + [number(r.get("Difference_CI95_High")) for r in rows]
    metric_name = first.get("Metric", "")
    unit = first.get("Metric_Unit", "")
    ylabel = ("Hedges g" if measure == "Hedges_g" else "%s minus %s" % (group_b, group_a))
    panel = Panel(svg, 230, 80, 630, 300, (0.5, len(rows) + 0.5), limits(ys, include_zero=True),
                  "", ylabel, "%s%s, C%s %s" % (metric_name, " (%s)" % unit if unit else "", first.get("Channel", ""),
                                                first.get("Channel_Name", "")))
    panel.axes(show_x_ticks=False)
    panel.hline(0.0, INK_3, None)
    same = SERIES[0]
    other = SERIES[1]
    for r in rows:
        colour = same if r.get("Same_Sign_As_Reference", "").lower() in ("true", "1", "yes") else other
        if measure == "Difference_B_Minus_A":
            lo, hi = number(r.get("Difference_CI95_Low")), number(r.get("Difference_CI95_High"))
            if lo is not None and hi is not None:
                svg.line(panel.sx(r["_x"]), panel.sy(lo), panel.sx(r["_x"]), panel.sy(hi), colour, 1.5)
        reference = r.get("Is_Reference", "").lower() in ("true", "1", "yes")
        panel.dot(r["_x"], r["_y"], colour, 5 if reference else 3.5,
                  "%s: %s = %s, Welch p %s" % (r.get("Specification", ""), measure, tick_text(round(r["_y"], 4)),
                                                r.get("Welch_P", "")), hollow=reference)
        if reference:
            svg.text(panel.sx(r["_x"]), panel.sy(r["_y"]) - 10, "configured", 11, INK_2, "middle")
    svg.legend(230, 36, [("same direction as configured", same), ("opposite direction", other)])
    if parameters:
        top = 80 + 300 + 50
        uses = set((p.get("Specification_ID"), p.get("Parameter"), p.get("Value")) for p in parameter_rows)
        for index, (name, value) in enumerate(parameters):
            y = top + index * row_h
            svg.text(220, y + 4, "%s = %s" % (name, value), 11, INK_2, "end")
            svg.line(panel.x0, y, panel.x0 + panel.w, y, GRID, 1)
            for r in rows:
                if (r.get("Specification_ID"), name, value) in uses:
                    svg.circle(panel.sx(r["_x"]), y, 3, INK_2)
    summary_file = os.path.join(folder, "robustness_summary.csv")
    if os.path.isfile(summary_file):
        summary = read_rows([summary_file])[0]
        line = "%s of %s combinations give the same direction as the configured settings" % (
            summary.get("N_Same_Sign_As_Reference", "?"), summary.get("N_Specifications_Compared", "?"))
        svg.text(230, svg.height - 16, line + " (descriptive; combinations are not independent).", 12, INK_2)
    base = os.path.join(args.output, "specification_curve")
    svg.save(base + ".svg")
    write_table(base + ".csv", ["Curve_Rank", "Specification_ID", "Specification", "Is_Reference", "Group_A", "Group_B",
                                "Mean_A", "Mean_B", "Difference_B_Minus_A", "Difference_CI95_Low", "Difference_CI95_High",
                                "Hedges_g", "Welch_P", "Mann_Whitney_P", "Same_Sign_As_Reference"], rows)
    return [base + ".svg"]


# -----------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(description="Make SVG figures from Cell Quant Studio result tables.")
    parser.add_argument("figure", choices=["validation", "depth", "spectrum", "robustness"])
    parser.add_argument("input", help="a result table, or a folder searched for it")
    parser.add_argument("-o", "--output", default="figures", help="folder for the figures (default: figures)")
    parser.add_argument("--group-regex", default=DEFAULT_GROUP_REGEX,
                        help="rule giving each image its group from its name (first bracketed part)")
    parser.add_argument("--metric", default=None,
                        help="depth: a column such as Positive_Area_Fraction_Percent; robustness: g or difference")
    parser.add_argument("--relative", action="store_true", help="depth: use relative depth (0 to 1)")
    parser.add_argument("--frame", type=int, default=None, help="only this frame")
    parser.add_argument("--title", default=None, help="figure title")
    args = parser.parse_args(argv)
    if not os.path.isdir(args.output):
        os.makedirs(args.output)
    makers = {"validation": figure_validation, "depth": figure_depth, "spectrum": figure_spectrum,
              "robustness": figure_robustness}
    for path in makers[args.figure](args):
        print("Wrote " + path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
