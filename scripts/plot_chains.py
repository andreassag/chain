#!/usr/bin/env python3
"""
plot_chains.py - Generate standalone SVG synteny dot-plot from PAF output.
"""

import sys
import xml.etree.ElementTree as ET

def parse_paf(filepath):
    records = []
    max_q = 0
    max_r = 0
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
            q_name = parts[0]
            q_len = int(parts[1])
            q_start = int(parts[2])
            q_end = int(parts[3])
            strand = parts[4]
            r_name = parts[5]
            r_len = int(parts[6])
            r_start = int(parts[7])
            r_end = int(parts[8])
            circ = "false"
            for tag in parts[12:]:
                if tag.startswith("circ:Z:"):
                    circ = tag[7:]

            max_q = max(max_q, q_len)
            max_r = max(max_r, r_len)
            records.append({
                'q_name': q_name, 'q_len': q_len, 'q_start': q_start, 'q_end': q_end,
                'strand': strand, 'r_name': r_name, 'r_len': r_len, 'r_start': r_start,
                'r_end': r_end, 'circ': circ
            })
    return records, max_q, max_r

def generate_svg(records, max_q, max_r, out_path):
    width = 800
    height = 800
    margin = 80

    plot_w = width - 2 * margin
    plot_h = height - 2 * margin

    svg = ET.Element('svg', {
        'xmlns': 'http://www.w3.org/2000/svg',
        'width': str(width),
        'height': str(height),
        'viewBox': f'0 0 {width} {height}',
        'style': 'background-color: #ffffff; font-family: sans-serif;'
    })

    # Title
    title = ET.SubElement(svg, 'text', {
        'x': str(width // 2), 'y': '35',
        'text-anchor': 'middle', 'font-size': '18', 'font-weight': 'bold', 'fill': '#24292f'
    })
    title.text = "Synteny Dot-Plot (Minimizer Co-Linear Chains)"

    # Axes
    ET.SubElement(svg, 'rect', {
        'x': str(margin), 'y': str(margin),
        'width': str(plot_w), 'height': str(plot_h),
        'fill': '#f6f8fa', 'stroke': '#d0d7de', 'stroke-width': '1.5'
    })

    # Axis Labels
    x_label = ET.SubElement(svg, 'text', {
        'x': str(margin + plot_w // 2), 'y': str(height - 25),
        'text-anchor': 'middle', 'font-size': '14', 'fill': '#57606a'
    })
    x_label.text = "Reference Coordinate (bp)"

    y_label = ET.SubElement(svg, 'text', {
        'x': str(25), 'y': str(margin + plot_h // 2),
        'text-anchor': 'middle', 'font-size': '14', 'fill': '#57606a',
        'transform': f'rotate(-90, 25, {margin + plot_h // 2})'
    })
    y_label.text = "Query Coordinate (bp)"

    if max_r == 0 or max_q == 0:
        with open(out_path, 'wb') as f:
            f.write(ET.tostring(svg))
        return

    # Draw chain segments
    for rec in records:
        x1 = margin + (rec['r_start'] / max_r) * plot_w
        x2 = margin + (rec['r_end'] / max_r) * plot_w
        if rec['strand'] == '+':
            y1 = margin + plot_h - (rec['q_start'] / max_q) * plot_h
            y2 = margin + plot_h - (rec['q_end'] / max_q) * plot_h
            color = '#0969da' # Blue for forward
        else:
            y1 = margin + plot_h - (rec['q_end'] / max_q) * plot_h
            y2 = margin + plot_h - (rec['q_start'] / max_q) * plot_h
            color = '#cf222e' # Red for reverse

        if rec['circ'] == 'true':
            color = '#8250df' # Purple for circular wrapping chains

        ET.SubElement(svg, 'line', {
            'x1': f'{x1:.1f}', 'y1': f'{y1:.1f}',
            'x2': f'{x2:.1f}', 'y2': f'{y2:.1f}',
            'stroke': color, 'stroke-width': '3', 'stroke-linecap': 'round'
        })

    with open(out_path, 'wb') as f:
        f.write(ET.tostring(svg))
    print(f"Rendered dot-plot to {out_path}")

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: plot_chains.py <input.paf> <output.svg>")
        sys.exit(1)
    recs, mq, mr = parse_paf(sys.argv[1])
    generate_svg(recs, mq, mr, sys.argv[2])
