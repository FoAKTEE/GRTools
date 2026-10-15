"""Verify registered finite windows and render candidate products (no admission)."""
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np
from plot_first_results import plot_case
import matplotlib.pyplot as plt

OBS=['torque_disc','power_disc','pdot','edot']

def read_modes(path):
    with path.open() as stream:
        return [{k:int(v) if k in ['m','k','nodes'] else v=='true' if k=='conditional_secular' else float(v) for k,v in row.items()} for row in csv.DictReader(stream)]

def decode_matrix(encoded):
    return np.asarray(encoded['data'],dtype=float).reshape(encoded['size'],order='F')

def plot_ladder(directory, result, summary, ladder):
    fig,axes=plt.subplots(4,3,figsize=(10,8),sharex='col',layout='constrained')
    names={'m':'m cutoff','negative_k':'negative k extent','positive_k':'positive k extent'}
    for j,(name,label) in enumerate(names.items()):
        entries=[row for row in ladder if row['axis']==name]
        for i,key in enumerate(OBS):
            axis=axes[i,j];total=summary['totals'][key]
            values=[row[key]/total for row in entries] if total else [np.nan]*len(entries)
            axis.plot([row['edge'] for row in entries],values,'o-',lw=1,ms=3)
            reference=summary.get('criterion_references',{}).get(key,abs(total))
            if total and (not result.get('ridge_registration') or key in OBS[:2]):
                width=.01*reference/abs(total)
                axis.axhspan(1-width,1+width,color='0.8',alpha=.5)
            axis.axhline(1,color='0.4',lw=.6,ls='--')
            axis.grid(alpha=.2)
            if j==0:axis.set_ylabel(key+' / final')
            if i==3:axis.set_xlabel(label)
            if i==0:axis.set_title(name.replace('_',' '))
    status='TRUNCATION PASS' if summary['criterion_met'] else 'UNCONVERGED'
    fig.suptitle(f"{status}: independent signed sums; m > 0, secular omitted\n"
                 f"chi={result['chi']:+.1f}, p/M=20, e={result['e']:.1f}; M={result['m_max']}, "
                 f"k=[{result['k_min']},{result['k_plus']}]\n"
                 'Other cutoffs held fixed; ±1% reference: |T|, gross |P| for ridge window')
    fig.savefig(directory/'convergence.png',dpi=140)
    plt.close(fig)

def rates(rows,jac,mu):
    torque=math.fsum(r['torque_disc'] for r in rows)
    power=math.fsum(r['power_disc'] for r in rows)
    pe=jac @ [-power/mu,-torque/mu]
    return dict(torque_disc=torque,power_disc=power,pdot=float(pe[0]),edot=float(pe[1]),Edot_orbit=-power,Ldot_orbit=-torque)

def fitted_tail(shells):
    """Conditional geometric extrapolation of three decreasing absolute shell sums."""
    if len(shells)!=3 or not all(math.isfinite(x) and x>0 for x in shells):
        return dict(valid=False,estimate=None,reason='need three positive finite shells',shells=shells)
    if not shells[0]>shells[1]>shells[2]:
        return dict(valid=False,estimate=None,reason='shells are not strictly decreasing',shells=shells)
    y=np.log(shells);slope,intercept=np.polyfit(np.arange(3),y,1)
    residual=float(np.sum((y-(intercept+slope*np.arange(3)))**2))
    r2=1-residual/float(np.sum((y-np.mean(y))**2))
    ratio=max(math.exp(slope),shells[-1]/shells[-2])
    valid=bool(slope<0 and r2>=.9 and ratio<1)
    return dict(valid=valid,estimate=float(shells[-1]*ratio/(1-ratio)) if valid else None,
                slope=float(slope),r2=r2,ratio=ratio,shells=shells,interpretation='fitted absolute tail; no rigorous bound')

def ridge_window(m, registration):
    w=registration['window']
    return (min(-w['negative_floor'],math.floor(w['k_apo']*m)-w['padding']),
            max(w['positive_floor'],math.ceil(w['k_peri']*m)+w['padding']))

def band_diagnostics(main, result):
    """Measured shells and descriptive fits; never add fitted tails to the totals."""
    M=result['m_max'];shells=[];fits={}
    soft=.018*result.get('p',20)/(1-result['e']**2)
    for band in ('apocentre','pericentre'):
        radius=result.get('p',20)/(1+(result['e'] if band=='pericentre' else -result['e']))
        for m in range(1,M+1):
            rows=[r for r in main if r['m']==m and (r['k']<=0 if band=='apocentre' else r['k']>0)]
            edge=min(r['k'] for r in rows) if band=='apocentre' else max(r['k'] for r in rows)
            row=dict(m=m,band=band,k_edge=edge,modes=len(rows))
            for key in OBS[:2]:
                peak=max(abs(r[key]) for r in rows)
                gross=math.fsum(abs(r[key]) for r in rows)
                row[key]=math.fsum(r[key] for r in rows)
                row['gross_'+key]=gross
                row['edge_peak_'+key]=sum(abs(r[key]) for r in rows if r['k']==edge)/peak if peak else 0.
                row['edge4_gross_'+key]=math.fsum(abs(r[key]) for r in rows if abs(r['k']-edge)<4)/gross if gross else 0.
            shells.append(row)
        for label,lo in [('last30',max(1,M-29)),('post_onset',math.floor(radius/soft)+1)]:
            rows=[r for r in shells if r['band']==band and lo<=r['m']<=M]
            for key in OBS[:2]:
                name=f'{band}_{key}_{label}'
                fit=dict(range=[lo,M],samples=len(rows),expected_efold=radius/(2*soft),valid=False)
                values=np.array([r[key] for r in rows]);ms=np.array([r['m'] for r in rows])
                if len(rows)>=3 and np.all(values!=0) and len(set(np.sign(values)))==1:
                    slope,intercept=np.polyfit(ms,np.log(np.abs(values)),1)
                    residual=np.log(np.abs(values))-(intercept+slope*ms)
                    variance=float(np.sum((np.log(np.abs(values))-np.mean(np.log(np.abs(values))))**2))
                    fit.update(slope=float(slope),intercept=float(intercept),sign=float(np.sign(values[0])),
                        efold_m=float(-1/slope) if slope<0 else None,
                        r2=1-float(np.sum(residual**2))/variance if variance>0 else None,valid=bool(slope<0))
                else:
                    fit['reason']='need at least three nonzero same-sign shells; post-onset data may be unavailable'
                fits[name]=fit
    return dict(shells=shells,fits=fits,interpretation='descriptive measurement; no extrapolated tail in finite totals')

def assess(modes,result):
    pairs=[(r['m'],r['k']) for r in modes]
    M,Km,Kp=result['m_max'],-result['k_min'],result['k_plus']
    ridge=result.get('ridge_registration')
    if ridge:
        excluded={tuple(v) for v in ridge['named_exclusions']}
        bounds=[(m,ridge_window(m,ridge)) for m in range(1,M+1)]
        assert all(-Km<=lo<=hi<=Kp for m,(lo,hi) in bounds),'envelope misses ridge'
        expected={(m,k) for m,(lo,hi) in bounds for k in range(lo,hi+1)}-excluded
    else:
        expected={(m,k) for m in range(1,M+1) for k in range(-Km,Kp+1)}
    assert len(pairs)==len(set(pairs)) and set(pairs)==expected,'incomplete or duplicate registered window'
    for row in modes:
        assert row['m']>=1 and math.isfinite(row['torque_disc']) and math.isfinite(row['power_disc'])
        frequency=row['m']*result['Omega_phi']+row['k']*result['Omega_r']
        assert math.isclose(row['omega'],frequency,rel_tol=1e-13,abs_tol=1e-16)
        secular=(row['m'],row['k'])==(1,-1) or abs(frequency)<=64*np.finfo(float).eps*max(abs(row['m']*result['Omega_phi']),abs(row['k']*result['Omega_r']))
        assert row['conditional_secular']==secular
        assert math.isclose(row['power_disc'],row['omega']/row['m']*row['torque_disc'],rel_tol=1e-9,abs_tol=1e-35)
    jac=decode_matrix(result['ladder'][-1]['jacobian'])
    main=[r for r in modes if not r['conditional_secular']]
    totals=rates(main,jac,result['mu'])
    references={key:abs(totals[key]) for key in OBS}
    if ridge:
        references['power_disc']=math.fsum(abs(r['power_disc']) for r in main)
    for key in OBS:
        assert math.isclose(totals[key],result['ladder'][-1][key],rel_tol=1e-12,abs_tol=1e-35)
    sectors={}
    for name,rows in [('paper',[r for r in main if r['m']+r['k']>0]),('negative_frequency',[r for r in main if r['omega']<0]),('secular',[r for r in modes if r['conditional_secular']]),('excluded_integer',[r for r in main if r['m']+r['k']<=0])]:
        values=rates(rows,jac,result['mu']);sectors[name]=dict(modes=len(rows),**values,
            torque_share_of_nonsecular=values['torque_disc']/totals['torque_disc'] if totals['torque_disc'] else None,
            power_share_of_nonsecular=values['power_disc']/totals['power_disc'] if totals['power_disc'] else None)
    axes={}; ladder=[]
    for axis,bound,width in [('m',M,10 if M>=30 else 5),('negative_k',Km,10 if Km>=30 else 5),('positive_k',Kp,5)]:
        def inside(row,edge):
            return row['m']<=edge if axis=='m' else row['k']>=-edge if axis=='negative_k' else row['k']<=edge
        entries=[]
        for edge in sorted({0,*range(bound%width,bound+1,width)}):
            entry=dict(axis=axis,edge=edge,**rates([r for r in main if inside(r,edge)],jac,result['mu']))
            entries.append(entry);ladder.append(entry)
        if not entries or entries[-1]['edge']!=bound or len(entries)<4:
            axes[axis]=dict(passed=False,reason='three complete equal-width shells unavailable');continue
        changes={key:[abs(entries[-i][key]-entries[-i-1][key])/references[key] if references[key] else None for i in [1,2]] for key in OBS}
        tails={}
        for key in OBS[:2]:
            shells=[math.fsum(abs(r[key]) for r in main if inside(r,edge) and not inside(r,edge-width)) for edge in [bound-2*width,bound-width,bound]]
            tail=fitted_tail(shells);tail['fraction_of_total']=tail['estimate']/references[key] if tail['estimate'] is not None and references[key] else None;tails[key]=tail
        passed=all(v is not None and v<.01 for vs in changes.values() for v in vs) and all(t['valid'] and t['fraction_of_total'] is not None and t['fraction_of_total']<.01 for t in tails.values())
        axes[axis]=dict(passed=passed,width=width,last_two_changes=changes,tails=tails)
    integrity=all(r['residual']<=1e-6 and r['bc_residual']<=1e-6 for r in main)
    summary=dict(totals=totals,sectors=sectors,axes=axes,solver_integrity=integrity,
        criterion_met=integrity and all(v['passed'] for v in axes.values()) and Km>M*result['Omega_phi']/result['Omega_r'],
        omega_zero_boundary_contained=Km>M*result['Omega_phi']/result['Omega_r'],modes=len(main),
        verification_checks=7*len(modes)+len(OBS)+1,criterion_references=references)
    if ridge:
        summary.update(criterion_met=False,tail_adopted=False,
            coverage_guard='registered ridge padding; measured edge fall-off, not omega=0 coverage',
            convergence_hole='finite M stage only; no adopted remaining-tail bound',
            bands=band_diagnostics(main,result))
    return summary,ladder

def main():
    parser=argparse.ArgumentParser();parser.add_argument('evidence',type=Path);parser.add_argument('--no-plot',action='store_true');args=parser.parse_args()
    summaries=[]
    for path in sorted(args.evidence.glob('a*/result.json')):
        result=json.loads(path.read_text())
        modes=read_modes(path.with_name('modes.csv'))
        summary,ladder=assess(modes,result)
        summary.update(case=path.parent.name,chi=result['chi'],e=result['e'],closure=result['closure'],domain=result['domain'],m_max=result['m_max'],k_min=result['k_min'],k_plus=result['k_plus'])
        previous=args.evidence.parent.parent/'gr-first-results'/path.parent.name/'result.json'
        if previous.exists():
            old=json.loads(previous.read_text())['ladder'][-1]
            summary['previous_truncated']={k:old[k] for k in OBS}
            jac=decode_matrix(result['ladder'][-1]['jacobian'])
            previous_k=12 if result['e']<.3 else 10
            rows=[r for r in modes if not r['conditional_secular']]
            summary['comparison_decomposition']={
                'retained_old_sector':rates([r for r in rows if r['m']<=old['m_max'] and abs(r['k'])<=previous_k and r['omega']>0],jac,result['mu']),
                'negative_frequency_in_old_rectangle':rates([r for r in rows if r['m']<=old['m_max'] and abs(r['k'])<=previous_k and r['omega']<0],jac,result['mu']),
                'added_k_shells':rates([r for r in rows if r['m']<=old['m_max'] and abs(r['k'])>previous_k],jac,result['mu']),
                'added_m_shells':rates([r for r in rows if r['m']>old['m_max']],jac,result['mu'])}
            for key in OBS:
                assert math.isclose(math.fsum(v[key] for v in summary['comparison_decomposition'].values()),summary['totals'][key],rel_tol=1e-12,abs_tol=1e-35)
            summary['change_from_previous']={k:(summary['totals'][k]-old[k])/abs(old[k]) for k in OBS}
        circular=args.evidence.parent.parent/'gr-first-results'/path.parent.name.replace('e0.2','e0.0').replace('e0.5','e0.0').replace('_li','_ad')/'result.json'
        if circular.exists():
            old=json.loads(circular.read_text())['ladder'][-1]
            summary['circular_m20_unconverged']={k:old[k] for k in OBS[:2]}
        if not args.no_plot:
            summary['figures']=plot_case(path.parent,result,status_label='TRUNCATION PASS' if summary['criterion_met'] else 'UNCONVERGED',spectrum_label='signed-frequency spectrum')
            plot_ladder(path.parent,result,summary,ladder)
        path.with_name('convergence.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
        if 'bands' in summary:
            shells=summary['bands']['shells']
            with path.with_name('band_shells.csv').open('w') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(shells[0]));writer.writeheader();writer.writerows(shells)
        with path.with_name('independent_ladder.csv').open('w') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(ladder[0]));writer.writeheader();writer.writerows(ladder)
        summaries.append(summary);print('VERIFY',path.parent.name,'PASS',summary['verification_checks'],'checks; truncation=',summary['criterion_met'],flush=True)
    (args.evidence/'case_summaries.json').write_text(json.dumps(summaries,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':
    main()
