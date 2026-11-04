"""Recheck source-derived native witnesses and saved pilot states; standard library only.

Run from repository root: python experiments/irback2025_audit/verify.py
"""
import itertools
import json
from pathlib import Path
import random
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from qubo_bench.mj85_benchmarks import BENCHMARKS, definition
from qubo_bench.model import build, energy

HERE = Path(__file__).resolve().parent
SCAN = json.loads((HERE / 'matrix_scan_transcription.json').read_text())
NATIVE = json.loads((HERE / 'native_folds.json').read_text())
ALPHABET = SCAN['alphabet']
TRIANGLE = SCAN['lower_triangle_hundredths']
SITES = list(itertools.product(range(4), range(4), range(3)))


def coefficient(a, b):
    i, j = ALPHABET.index(a), ALPHABET.index(b)
    return TRIANGLE[max(i, j)][min(i, j)]


def direct(sequence, occupied, penalties):
    """Equations (3)-(6), integer counts, fresh-scan coefficients, geometric distances.

    occupied is a sequence of sets of coordinates, including arbitrary invalid rows.
    Does not use the builder's contact matrix, edges, or validator.
    """
    e1 = sum((len(row)-1)**2 for row in occupied)
    e2 = sum(len(a & b) for i,a in enumerate(occupied) for b in occupied[i+1:])
    e3 = sum(sum(sum((x-y)**2 for x,y in zip(u,v)) > 1 for u in a for v in b)
             for a,b in zip(occupied, occupied[1:]))
    contacts = [(i,j) for i,a in enumerate(occupied) for j in range(i+2,len(sequence))
                for u in a for v in occupied[j]
                if sum(abs(x-y) for x,y in zip(u,v)) == 1]
    physical = sum(coefficient(sequence[i],sequence[j]) for i,j in contacts)/100
    valid = e1 == e2 == e3 == 0
    return dict(interaction_energy=physical, penalties=[e1,e2,e3], valid=valid,
                physical_energy=physical if valid else None,
                qubo_energy=physical+sum(a*b for a,b in zip(penalties,(e1,e2,e3))),
                contact_count=len(contacts), contacts=contacts)


def bits_for(coords, metadata):
    sites=[tuple(p) for p in metadata['coordinates']]
    return ''.join('1' if p == tuple(coords[i]) else '0'
                   for i in range(48) for p in sites)


def score_bits(sequence, bits, metadata):
    coords=[tuple(p) for p in metadata['coordinates']]
    occupied=[{p for u,p in enumerate(coords) if bits[i*48+u]=='1'} for i in range(48)]
    return direct(sequence,occupied,metadata['penalties'])


def main():
    resource=json.loads((ROOT/'qubo_bench/mj85_table_vi.json').read_text())
    assert resource['alphabet']==ALPHABET
    assert resource['lower_triangle_hundredths']==TRIANGLE
    witnesses=[]
    rng=random.Random(923)
    for identifier,(topology,sequence,ref) in BENCHMARKS.items():
        native=NATIVE[topology]
        for penalties in ((1.5,2,2),(3,4,4)):
            q,n,offset,meta=build(sequence,(4,4,3),penalties,'MJ85_TableVI','x-fastest')
            bits=bits_for(native['coordinates'],meta)
            s=score_bits(sequence,bits,meta)
            assert s['valid'] and s['penalties']==[0,0,0] and s['contact_count']==57
            assert sorted(s['contacts'])==[tuple(pair) for pair in native['contacts']]
            assert round(s['physical_energy']*100)==ref
            assert abs(energy(q,offset,bits)-s['physical_energy'])<1e-8
            witnesses.append(dict(benchmark=identifier,penalty_weights=penalties,
                                  physical=s['physical_energy'],qubo=energy(q,offset,bits)))
        # Arbitrary invalid states test the unconstrained polynomial too.
        for density in (0,0.02,0.5,1):
            bits=''.join('1' if rng.random()<density else '0' for _ in range(n))
            s=score_bits(sequence,bits,meta)
            assert abs(energy(q,offset,bits)-s['qubo_energy'])<1e-6
    summary=[]
    for path in sorted(HERE.glob('mj*_*.jsonl')):
        identifier='irback2025_'+path.name.split('_')[0]
        item=definition(identifier)
        q,n,offset,meta=build(item['sequence'],(4,4,3),(1.5,2,2),'MJ85_TableVI','x-fastest')
        for record in map(json.loads,path.read_text().splitlines()):
            row=dict(file=path.name,seed=record['seed'],runtime_seconds=record['runtime_seconds'],
                     attempted_updates=record['attempted_updates'],accepted_updates=record['accepted_updates'])
            for prefix in ('best','final'):
                s=score_bits(item['sequence'],record[prefix+'_bits'],meta)
                assert abs(s['qubo_energy']-record[prefix+'_energy'])<1e-7
                assert abs(energy(q,offset,record[prefix+'_bits'])-s['qubo_energy'])<1e-7
                row[prefix]={k:v for k,v in s.items() if k!='contacts'}
                row[prefix]['matches_reference']=s['valid'] and round(s['physical_energy']*100)==round(item['reference']['energy']*100)
            row['stages']=[]
            for point in record['progress']:
                if 'bits' not in point: continue
                s=score_bits(item['sequence'],point['bits'],meta)
                assert abs(s['qubo_energy']-point['current_energy'])<1e-7
                row['stages'].append(dict(sweep=point['sweep'],beta=point['beta'],
                    best_energy=point['best_energy'],**{k:v for k,v in s.items() if k!='contacts'}))
            summary.append(row)
    result=dict(matrix_entries_checked=210,witness_checks=witnesses,invalid_state_checks=24,trials=summary)
    (HERE/'verified_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='trials'},indent=2))
    for row in summary:
        print(row['file'],row['seed'],row['best'],row['final'],row['runtime_seconds'])


if __name__=='__main__': main()
