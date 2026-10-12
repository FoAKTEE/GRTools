"""Render candidate evidence and independently check modal sums; no ledger writes."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import struct

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def wake(path):
    with path.open('rb') as stream:
        nr, nphi = struct.unpack('qq', stream.read(16))
        r = np.fromfile(stream, dtype='<f8', count=nr)
        phi = np.fromfile(stream, dtype='<f8', count=nphi)
        field = np.fromfile(stream, dtype='<f8').reshape((nr, nphi), order='F')
    return r, phi, field


def truncation(result):
    rows = result['ladder']
    best = max(rows, key=lambda row: (row['m_max'], row['k_max']))
    prior_m = [row for row in rows if row['k_max'] == best['k_max'] and row['m_max'] < best['m_max']]
    prior_k = [row for row in rows if row['m_max'] == best['m_max'] and row['k_max'] < best['k_max']]
    trend = {}
    for axis, prior in [('m', prior_m), ('k', prior_k)]:
        previous = max(prior, key=lambda row: row[axis + '_max']) if prior else None
        trend[axis] = {}
        for quantity in ['torque_disc', 'power_disc', 'pdot', 'edot']:
            if previous is None or best[quantity] is None:
                change = None
            else:
                numerator = abs(best[quantity] - previous[quantity])
                change = numerator / abs(best[quantity]) if best[quantity] else (0.0 if numerator == 0 else None)
            trend[axis][quantity] = change
    required = [v for values in trend.values() for v in values.values()]
    passed = all(v is not None and v < .01 for v in required)
    return best, trend, passed


def plot_case(directory, result):
    profile = np.genfromtxt(directory / 'profiles.csv', delimiter=',', names=True)
    modes = np.genfromtxt(directory / 'modes.csv', delimiter=',', names=True)
    closure_label='adiabatic' if result['closure']=='Adiabatic' else 'LI option A'
    label = rf"$\chi={result['chi']:+.1f},\ p/M=20,\ e={result['e']:.1f}$" + ' — ' + closure_label
    scale = result['epsilon_ref'] * result['mu'] ** 2
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 6.5), sharex=True, layout='constrained')
    for axis, column, ylabel in zip(axes,
            ['amf_disc', 'torque_density_disc', 'power_density_disc'],
            [r'$F_J/(\epsilon_{ref}\mu^2)$', r'$(dT/dr)/(\epsilon_{ref}\mu^2)$', r'$(dP/dr)/(\epsilon_{ref}\mu^2)$']):
        axis.plot(profile['r'], profile[column] / scale, lw=.9)
        axis.axhline(0, color='0.6', lw=.5)
        for r in [20/(1+result['e']), 20/(1-result['e'])]:
            axis.axvline(r, color='0.6', lw=.6, ls='--')
        axis.set_ylabel(ylabel)
        axis.grid(alpha=.15)
    axes[-1].set_xlabel(r'$r/M$')
    fig.suptitle('PRELIMINARY: disc profiles; (1,−1) omitted\n' + label)
    fig.savefig(directory / 'profiles.png', dpi=140)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), layout='constrained')
    for axis, column, title in zip(axes, ['torque_disc', 'power_disc'], ['|torque|', '|power|']):
        spectrum = np.full((2*result['k_max']+1, result['m_max']), np.nan)
        for row in np.atleast_1d(modes):
            spectrum[int(row['k'])+result['k_max'], int(row['m'])-1] = abs(row[column]) / scale
        spectrum = np.ma.masked_where(~np.isfinite(spectrum) | (spectrum <= 0), spectrum)
        image = axis.imshow(np.ma.log10(spectrum), origin='lower', aspect='auto',
                            extent=[.5,result['m_max']+.5,-result['k_max']-.5,result['k_max']+.5], cmap='viridis')
        axis.set(xlabel='m', ylabel='k', title=title + r' / $(\epsilon_{ref}\mu^2)$')
        if result['e']:
            axis.plot(1,-1,'rx',ms=6)
        fig.colorbar(image,ax=axis,label='log10 amplitude')
    fig.suptitle('PRELIMINARY: positive-frequency spectrum; red × = conditional\n' + label)
    fig.savefig(directory / 'spectra.png', dpi=140)
    plt.close(fig)

    bankdir = Path(result['source_bank']).parent
    r, phi, peri = wake(bankdir / 'wake_peri.bin')
    ra, phia, apo = wake(bankdir / 'wake_apo.bin')
    assert np.array_equal(r,ra) and np.array_equal(phi,phia)
    assert np.isfinite(peri).all() and np.isfinite(apo).all()
    # One explicit robust limit across both phases, without changing numerical data.
    limit = float(np.quantile(np.abs(np.concatenate([peri.ravel(),apo.ravel()])), .995))
    if limit == 0:
        limit = 1.0
    X, Y = r[:,None]*np.cos(phi), r[:,None]*np.sin(phi)
    for name, field, radius, angle in [('peri',peri,result['wake']['r_peri'],0.0),
                                      ('apo',apo,result['wake']['r_apo'],result['wake']['phi_apo'])]:
        fig, axis = plt.subplots(figsize=(6,5.5), layout='constrained')
        image = axis.pcolormesh(X,Y,field,shading='gouraud',cmap='RdBu_r',vmin=-limit,vmax=limit,rasterized=True)
        axis.plot(radius*np.cos(angle),radius*np.sin(angle),'ko',ms=4)
        axis.add_patch(plt.Circle((0,0),result['domain'][0],color='white'))
        axis.set(aspect='equal',xlabel=r'$x/M$',ylabel=r'$y/M$')
        phase_name=('pericentre' if name=='peri' else 'apocentre') if result['e'] else (
            'circular phase 0' if name=='peri' else 'circular half radial-limit period')
        axis.set_title('PRELIMINARY: '+phase_name+'\n'+label)
        fig.colorbar(image,ax=axis,label=r'$\delta\epsilon/\epsilon$ (Eulerian)' + '\nshared 99.5% colour limits (clipped)')
        fig.savefig(directory / f'wake_{name}.png',dpi=140)
        plt.close(fig)
    return {'wake_color_limit':limit,'wake_max_abs':float(max(np.max(np.abs(peri)),np.max(np.abs(apo)))),
            'color_limit_rule':'joint 99.5 percentile of absolute fields; clipped display only'}


def verify_case(directory, result):
    with (directory/'modes.csv').open() as stream:
        modes = list(csv.DictReader(stream))
    expected={(m,k) for m in range(1,result['m_max']+1)
              for k in range(-result['k_max'],result['k_max']+1)
              if m*result['Omega_phi']+k*result['Omega_r']>0}
    observed=[(int(row['m']),int(row['k'])) for row in modes]
    assert len(observed)==len(set(observed)) and set(observed)==expected, 'Incomplete/duplicate mode rectangle'
    assert all((row['conditional_secular']=='true') == ((int(row['m']),int(row['k']))==(1,-1)) for row in modes)
    checks = 0
    for row in result['ladder']:
        selected = [m for m in modes if int(m['m'])<=row['m_max'] and abs(int(m['k']))<=row['k_max']
                    and m['conditional_secular']=='false']
        assert len(selected) == row['modes']
        for key in ['torque_disc','power_disc']:
            total = math.fsum(float(m[key]) for m in selected)
            assert math.isclose(total,row[key],rel_tol=5e-13,abs_tol=1e-40), (directory,key,total,row[key])
            checks += 1
        assert row['Edot_orbit'] == -row['power_disc']
        assert row['Ldot_orbit'] == -row['torque_disc']
        assert math.isclose(row['Edot_specific'],row['Edot_orbit']/result['mu'],rel_tol=2e-15)
        assert math.isclose(row['Ldot_specific'],row['Ldot_orbit']/result['mu'],rel_tol=2e-15)
        if row['pdot'] is not None:
            J = np.array(row['jacobian']['data']).reshape(row['jacobian']['size'],order='F')
            mapped = J @ [row['Edot_specific'],row['Ldot_specific']]
            assert np.allclose(mapped,[row['pdot'],row['edot']],rtol=2e-12,atol=1e-35)
        checks += 5
    for mode in modes:
        assert float(mode['omega']) > 0
        expected = int(mode['m'])*result['Omega_phi']+int(mode['k'])*result['Omega_r']
        assert math.isclose(float(mode['omega']),expected,rel_tol=1e-14)
        assert math.isfinite(float(mode['torque_disc'])) and math.isfinite(float(mode['power_disc']))
        checks += 3
    bankhash=hashlib.file_digest(Path(result['source_bank']).open('rb'),'sha256').hexdigest()
    assert bankhash == result['source_bank_sha256']
    assert result['p'] > result['separatrix']
    print(f'VERIFY {directory.name}: PASS {checks+2} checks',flush=True)
    return checks+2


def main():
    p=argparse.ArgumentParser()
    p.add_argument('evidence',type=Path)
    p.add_argument('--verify-only',action='store_true')
    a=p.parse_args()
    assert a.evidence.is_dir() and (a.evidence/'PREREG.json').is_file(), 'Missing evidence/preregistration'
    cases=[]
    checks=0
    commands=[json.loads(line) for line in (a.evidence/'commands.jsonl').read_text().splitlines()]
    for path in sorted(a.evidence.glob('a*_e*_*ad/result.json')) + sorted(a.evidence.glob('a*_e*_*li/result.json')):
        result=json.loads(path.read_text())
        checks+=verify_case(path.parent,result)
        best,trend,passed=truncation(result)
        figures = {} if a.verify_only else plot_case(path.parent,result)
        process_cost=[command['cpu_seconds'] for command in commands
                      if Path(command['output']).name.startswith('case_'+result['case'])
                      and 'first_results.jl case ' in command['command']]
        cases.append(dict(case=result['case'],chi=result['chi'],e=result['e'],closure=result['closure'],
            totals=best,last_shell_relative=trend,truncation_pass=passed,
            torque_changes_sign_on_ladder=(min(row['torque_disc'] for row in result['ladder'])<0<
                                           max(row['torque_disc'] for row in result['ladder'])),
            cpu_hours=(sum(process_cost) if process_cost else result['cpu_seconds'])/3600,
            cpu_cost_scope='all recorded case processes, including failed/stopped original runs' if process_cost else 'driver clock (process record unavailable)',
            last_run_bank_build_or_load_cpu_hours=result['bank_cpu_seconds']/3600,
            omitted_static=result['excluded_static_count'],omitted_retrograde=result['excluded_retrograde_count'],
            secular=result['secular'],max_residual=result['max_residual'],
            max_bc_residual=result['max_bc_residual'],domain=result['domain'],
            original_domain=result.get('original_domain',result['domain']),
            domain_amendment=result.get('domain_amendment'),figure_dir=str(path.parent),**figures))
    expected_cases={(chi,e,'Adiabatic') for chi in (-.9,0.,.9) for e in (0.,.2,.5)}
    observed_cases={(case['chi'],case['e'],case['closure']) for case in cases}
    assert expected_cases<=observed_cases, f'Missing requested cases: {expected_cases-observed_cases}'
    for case in cases:
        circular=next(c for c in cases if c['chi']==case['chi'] and c['e']==0 and c['closure']=='Adiabatic')
        case['circular_comparison']={key+'_ratio':case['totals'][key]/circular['totals'][key]
                                    for key in ('torque_disc','power_disc')}
        case['circular_comparison']['same_radial_domain']=case['domain']==circular['domain']
        if case['closure']=='LocallyIsothermalOptionA':
            adiabatic=next(c for c in cases if c['chi']==case['chi'] and c['e']==case['e'] and c['closure']=='Adiabatic')
            case['adiabatic_comparison']={key+'_relative_difference':abs(case['totals'][key]-adiabatic['totals'][key])/abs(adiabatic['totals'][key])
                                         for key in ('torque_disc','power_disc','pdot','edot')}
    if a.verify_only:
        print(f'VERIFICATION PASS: {checks} checks across {len(cases)} cases')
        return
    summary=dict(status='PRELIMINARY candidate evidence; reviewer admission pending',cases=cases,
        verification_checks=checks,units='G=M=c=1; mu=1e-5; epsilon_ref=1e-8; rates per coordinate time M',
        wake_coordinates='Boyer-Lindquist polar projection: x=r*cos(phi), y=r*sin(phi)',
        model_assumptions='Softened equal-time leading singular-field forcing; pSigma=1.5, h_ref=.06, b=.3, qTemp=0, epsilon_p=1e-5',
        scaling='T,P scale as epsilon_ref*mu^2; specific rates and pdot,edot as epsilon_ref*mu',
        scope='m>=1 and omega>0 only; (1,-1) conditional and excluded from main totals and wakes',
        holes=['gr-radial-convergence and gr-mode-convergence remain open',
               'positive-frequency restriction omits retrograde/static and all m=0 power',
               'conditional secular response requires its separate certificate',
               'circular orbital Jacobian is degenerate; pdot and edot not evaluated',
               'candidate radial rule transferred to p=20,e=.2/.5, not certified there',
               'requested M40/60 and larger k rectangles deferred under the preregistered cost cap',
               'no complete physical feedback or astrophysical density normalization claimed',
               'circular e=0 boundary-convention regression remains open at pinned HEAD',
               'wake radial sampling and profile display interpolation not convergence certified',
               'five primary eccentric cases and both LI cases use amended common domains; omitted radial intervals and boundary sensitivity unmeasured',
               'AD/LI comparisons use different boundary prescriptions; a matched-boundary comparison remains open'],
        decisions=json.loads((a.evidence/'CAMPAIGN.json').read_text()) if (a.evidence/'CAMPAIGN.json').exists() else {},
        domain_amendment=json.loads((a.evidence/'DOMAIN_AMENDMENT.json').read_text()) if (a.evidence/'DOMAIN_AMENDMENT.json').exists() else None,
        optional_decisions=json.loads((a.evidence/'OPTIONAL_DECISION.json').read_text()) if (a.evidence/'OPTIONAL_DECISION.json').exists() else None,
        requested_cases=9,independent_review='separate reviewer/lander; not performed by this lane')
    (a.evidence/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
    with (a.evidence/'conditional_secular.csv').open('w') as stream:
        writer=csv.writer(stream)
        keys=['torque_disc','power_disc','Edot_orbit','Ldot_orbit','pdot','edot']
        writer.writerow(['case','closure','chi','e','m','k',*keys])
        for case in cases:
            for row in case['secular']:
                writer.writerow([case['case'],case['closure'],case['chi'],case['e'],row['m'],row['k'],
                                 *(row[key] for key in keys)])
    lines=['# PRELIMINARY eccentric Kerr–disc results','',
        'Candidate evidence only; no ledger admission. '+summary['units']+'. '+summary['scaling']+'.',
        '',summary['model_assumptions']+'. '+summary['scope']+'. The reference radius is p/(1−e²), so each e uses its stated reference disc.',
        '', '| closure | chi | e | m,k max | T into disc | P into disc | dp/dt | de/dt | last m/k torque change | CPU h |',
        '|:---|---:|---:|:---:|---:|---:|---:|---:|:---:|---:|']
    fmt=lambda x: '—' if x is None else f'{x:.5e}'
    for case in cases:
        t=case['totals']; change=case['last_shell_relative']
        pct=lambda x: '—' if x is None else f'{100*x:.1f}%'
        closure='AD' if case['closure']=='Adiabatic' else 'LI-A'
        case_link=f"[{case['chi']:+.1f}]({Path(case['figure_dir']).name}/profiles.png)"
        lines.append(f"| {closure} | {case_link} | {case['e']:.1f} | {t['m_max']},{t['k_max']} | {fmt(t['torque_disc'])} | {fmt(t['power_disc'])} | {fmt(t['pdot'])} | {fmt(t['edot'])} | {pct(change['m']['torque_disc'])}/{pct(change['k']['torque_disc'])} | {case['cpu_hours']:.3f} |")
    passed_cases=sum(case['truncation_pass'] for case in cases)
    nonlinear=[case['case'] for case in cases if case['wake_max_abs']>1]
    sign_changes=[case['case'] for case in cases if case['torque_changes_sign_on_ladder']]
    lines += ['', f'{passed_cases}/{len(cases)} cases meet the preregistered truncation criterion. ' +
        ('Sampled wakes exceed |delta epsilon / epsilon|=1 for '+', '.join(nonlinear)+
         '; the linear-response interpretation fails at those sampled points.' if nonlinear else
         'All displayed-grid wake samples have |delta epsilon / epsilon|<=1; this is not a linearity certificate.'),
        '', 'Torque changes sign across measured ladder entries for '+', '.join(sign_changes)+'.'
        if sign_changes else '',
        '', 'Orbit totals are dE/dt=−P and dL/dt=−T. Each case has `ladder.csv` (all totals/rates), '
        '`modes.csv`, `profiles.csv`, `profiles.png`, `spectra.png`, `wake_peri.png`, and `wake_apo.png`. '
        'Totals use native solver grids; profile CSVs and wake images are sampled displays. '
        'Both wake phases share one colour scale and use the Boyer–Lindquist polar projection x=r cos(phi), y=r sin(phi). Separately conditional (1,−1) transfers are in '
        '[conditional_secular.csv](conditional_secular.csv) and are excluded from these totals. '
        'The circular cases supply torque/power comparisons; '
        'the eccentric Jacobian intentionally rejects e=0. Circular snapshots use reference phases '
        'separated by half the limiting radial period, since a circle has no distinct peri/apo.', '',
        'Truncation passes only if both last measured m/k shells change every rate and total by <1%; '
        'all measured changes, including failures, are in SUMMARY.json. This is not convergence certification.', '',
        '**Open obligations:** '+ '; '.join(summary['holes'])+'.', '',
        f'Independent arithmetic verifier: {checks} checks passed. Source banks and raw wake matrices are in '
        '`/tmp/chandra/emri-ecc/banks/w18-g49/`; hashes and commands are in RUN.json.']
    (a.evidence/'SUMMARY.md').write_text('\n'.join(lines)+'\n')
    print(f'RENDERED {len(cases)} cases; VERIFICATION PASS {checks} checks')

if __name__=='__main__':
    main()
