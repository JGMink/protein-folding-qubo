"""Optional C++ diagnostics on native 2,304-bit assignments (run in Linux).

python3 experiments/irback2025_audit/check_cpp.py
"""
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from qubo_bench.io import read_qubo
from qubo_bench.model import energy
from verify import NATIVE,bits_for

exe='qubo_annealing/qubo_anneal'
for short,topology in [('mj1','A'),('mj5','B')]:
    path=f'build/audit/{short}.qubo'
    q,n,offset,meta=read_qubo(path)
    bits=bits_for(NATIVE[topology]['coordinates'],meta)
    def call(*args):
        return json.loads(subprocess.check_output([exe,path,*map(str,args)],text=True))
    assert abs(call('--energy',bits)['energy']-energy(q,offset,bits))<1e-8
    for index in (0,47,48,1023,2047,2303):
        flipped=bits[:index]+str(1-int(bits[index]))+bits[index+1:]
        delta=energy(q,offset,flipped)-energy(q,offset,bits)
        assert abs(call('--delta',bits,index)['delta']-delta)<1e-8
    print(short, 'native C++ energy and six bit deltas match full Python QUBO')
