"""Digitize the unmodified 720x599 2006 Figure 1 (see audit report for URL).

Run from repository root; outputs build/audit/map_pixels.json for reconstruct.py.
No sequences, matrix coefficients, or reference energies enter this extraction.
"""
from PIL import Image
import json
im = Image.open('build/audit/native_maps.jpg').convert('L')
assert im.size == (720, 599)
out = {}
for name, x0, x1 in [('A', 54, 341), ('B', 428, 714)]:
    pairs = []
    for i in range(1, 49):
        for j in range(1, 49):
            x = round(x0 + i * (x1-x0)/50)
            y = round(544-j*(544-258)/50)
            if im.getpixel((x,y)) < 100:
                pairs.append((i-1,j-1))
    out[name] = pairs
    assert len(pairs) == 114 and len(set(pairs)) == 114
    assert all((j,i) in pairs and abs(i-j)>1 for i,j in pairs)
    print(name, len(pairs), all((j,i) in pairs for i,j in pairs))
with open('build/audit/map_pixels.json','w') as f:
    json.dump(out,f)
