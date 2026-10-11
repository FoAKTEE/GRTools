"""Data-only FD25 setup check; no author module is imported or executed.

Sources: ref-paper/arXiv-2506.19917/src/Sections/4_results.tex:61,125;
ref-code/callumf8-linear_ecc_planet/README.md:3-9. Limits below are the
scoped profile-characterisation ceilings of solid fd25-oracle-conventions,
not a new computed-response versus author-data comparison.
"""
from bisect import bisect_right
import hashlib
import json
import math
from pathlib import Path
import sys


def check_setup_oracle(root):
    discs=(("1.5","1.0"),("1.5","0.0"),("0.5","0.0"))
    eccentricities=("0.01","0.06","0.12","0.3")
    # Per-case bounds preserve the known physical jump and unexplained features.
    bounds={"1.5":(.004,.020,.016,.007),"0.5":(.060,.060,.050,.022)}
    rows=[]; hashes={}
    for p,q in discs:
        for j,e in enumerate(eccentricities):
            directory=root/f"p_{p}_q_{q}"/f"e_{e}"
            arrays=[]
            for name in ("rgrid.out","dTdr.out","amf.out"):
                file=directory/name; raw=file.read_bytes()
                hashes[str(file.relative_to(root))]=hashlib.sha256(raw).hexdigest()
                arrays.append([float(line) for line in raw.decode('ascii').splitlines() if line.strip()])
            r,t,f=arrays
            assert len(r)==len(t)==len(f)==10000
            assert all(math.isfinite(v) for ar in arrays for v in ar)
            assert all(b>a for a,b in zip(r,r[1:]))
            grid_error=max(abs(v/(.05*100**(i/9999))-1) for i,v in enumerate(r))
            assert grid_error<5e-7
            # Lander addition (lane w17-g43 review): the files print 19 digits, but every node is the
            # 7-significant-digit rounding of the formula, so 5e-7 is the rounding ceiling, not a margin.
            assert all(float('%.6e'%(.05*100**(i/9999)))==v for i,v in enumerate(r))
            assert f[0]>0 and f[-1]>0
            cumulative=[0.]
            for i in range(1,len(r)):
                cumulative.append(cumulative[-1]+.5*(t[i]+t[i-1])*(r[i]-r[i-1]))
            row=dict(p=float(p),q=float(q),e=float(e),grid_relative_error=grid_error,
                     positive_edge_flux=True,integrated_normalized_disc_torque=cumulative[-1])
            if q=='0.0':
                residual=max(abs(v-f[0]-c) for v,c in zip(f,cumulative))
                assert residual<bounds[p][j]
                row.update(balance_residual=residual,balance_bound=bounds[p][j])
            else:
                spreads=[]
                # Match the solid convention certificate's interpolation locations.
                for endpoints in ((.1,.2),(3.,4.5)):
                    v=[]
                    for x in endpoints:
                        k=bisect_right(r,x)-1
                        y=f[k]+(f[k+1]-f[k])*(x-r[k])/(r[k+1]-r[k])
                        v.append(x*y)
                    spreads.append((max(v)-min(v))/max(abs(z) for z in v))
                assert max(spreads)<.012
                row.update(thermal_flux_spread=spreads,thermal_flux_bound=.012)
            rows.append(row)
    return dict(cases=rows,input_sha256=hashes,
        conventions="dT/dr in FJ0/a; FJ in FJ0; torque on disc; arrays already normalized",
        limitation="Absolute normalization is grounded in the paper and the solid convention node; a common scale factor cannot be inferred from balance alone. No rates are supplied by these 12 profiles.")


if __name__=='__main__':
    result=check_setup_oracle(Path(sys.argv[1]))
    Path(sys.argv[2]).write_text(json.dumps(result,indent=2)+'\n')
    print('FD25 DATA-ONLY SETUP PASS: 12 profiles, 36 inputs hashed; derived quantities only')
    print('max relative log-grid residual',max(r['grid_relative_error'] for r in result['cases']))
