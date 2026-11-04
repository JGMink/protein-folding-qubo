"""Embed the source contact graphs in the full grid, then score the result.

Requires NetworkX; run digitize.py first from repository root. Coordinates
are determined solely by graph isomorphism, before consulting any energies.
"""
import sys,json
sys.path[:0]=['.', 'build/pylibs']
import networkx as nx
from qubo_bench.mj85_benchmarks import definition,evaluate_coordinates
maps=json.load(open('build/audit/map_pixels.json'))
grid=nx.grid_graph(dim=[3,4,4])
out={}
for name,pairs in maps.items():
    contacts=sorted((i,j) for i,j in pairs if i<j)
    g=nx.Graph(); g.add_nodes_from(range(48)); g.add_edges_from(contacts)
    g.add_edges_from((i,i+1) for i in range(47))
    mapping=next(nx.isomorphism.GraphMatcher(g,grid).isomorphisms_iter())
    coords=[mapping[i] for i in range(48)]
    out[name]={'contacts':contacts,'coordinates':coords}
    for k in (range(1,4) if name=='A' else range(4,7)):
        key=f'irback2025_mj{k}'
        print(key,evaluate_coordinates(key,coords)['energy'])
json.dump(out,open('build/audit/native_folds.json','w'),indent=2)
