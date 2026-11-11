"""Validate the B3 v3 workflow with synthetic data; no AoU or PLINK execution.

Runs the notebook's actual sampling, GWAS parsing/selection, ALT harmonization,
MR and conditional covariance code. An independent OLS adapter supplies PLINK-
shaped synthetic reports/exports. This does NOT verify a real PLINK executable.
Run with a Python containing numpy, pandas and scipy.
"""
import ast
import hashlib
import json
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import linalg, sparse, stats


ROOT = Path(__file__).resolve().parents[1]
NB = ROOT / 'notebooks/AGEMR_B3_IVW_robust_MR_age_mismatch_v3.ipynb'


def validate():
    nb = json.loads(NB.read_text())
    code = {i: ''.join(c['source']) for i, c in enumerate(nb['cells']) if c['cell_type'] == 'code'}
    for i, source in code.items():
        compile(source, f'cell_{i}', 'exec')
    assert all(not c.get('outputs') for c in nb['cells'])
    g = dict(np=np, pd=pd, stats=stats, linalg=linalg, sparse=sparse,
             hashlib=hashlib, Path=Path, json=json, warnings=warnings,
             display=lambda *args: None, gc=__import__('gc'))
    tree = ast.parse(code[1])
    functions = [n for n in tree.body if isinstance(n, ast.FunctionDef)]
    exec(compile(ast.Module(body=functions, type_ignores=[]), '<setup helpers>', 'exec'), g)
    g.update(SEED=20261028, MIN_ROLE_N=100, AGE_BIN_YEARS=5,
             N_TARGET_REQUESTED=700, N_SOURCE_REQUESTED=700, N_REPEATS=3,
             N_BOOT=400, IV_P_THRESHOLD=1e-5, MAX_UNION_IV=128, THREADS=1,
             MEMORY_MB=1024, UNRELATED_KEEP=None, PLINK='synthetic-adapter', PLINK_VERSION='synthetic adapter',
             METHODS=('IVW', 'dIVW', 'weighted_median', 'MR_Egger'), OUTCOMES=('SBP', 'DBP'),
             SCENARIOS=[('matched', 0.), ('younger_moderate', -.5), ('younger_strong', -1.),
                        ('older_moderate', .5), ('older_strong', 1.)],
             KEY=['#CHROM', 'POS', 'ID', 'REF', 'ALT'], CONTRAST_FAMILY=32)
    rng = np.random.default_rng(702)
    n, k = 6500, 7
    ids = pd.Index([f'i{i:05d}' for i in range(n)], name='IID')
    age = rng.uniform(25, 80, n)
    C = pd.DataFrame({'AGE10': (age-50)/10, 'AGE2': ((age-50)/10)**2,
                      'SEX': rng.integers(0, 2, n), 'PC1': rng.normal(size=n),
                      'EXACT_DUPLICATE': (age-50)/5}, index=ids)
    G = rng.binomial(2, np.array([.25, .8, .35, .3, .4, .2, .4]), size=(n, k)).astype(float)
    y = (24 + .4*C.AGE10.to_numpy() + G[:, 0]*1.2 - G[:, 1]*1.1
         + G[:, 2]*(.10+.08*(age-50)) + .85*G[:, 4]
         + .8*G[:, 6] + rng.normal(size=n))
    target_mask = np.array([g['role'](iid) == 'target' for iid in ids])
    G[target_mask, 4] = 0  # Selected source IV can fail outcome harmonization.
    G[rng.random(G.shape) < .015] = np.nan
    # Phenotypes precede genotype missingness; no imputation in association fits.
    sbp = 95 + .5*y + rng.normal(0, 3, n)
    dbp = 65 + .25*y + rng.normal(0, 2, n)
    panel = pd.DataFrame({'#CHROM': ['1']*k, 'POS': np.arange(1, k+1)*1000,
                          'ID': [f's{j}' for j in range(k)], 'REF': ['A']*k, 'ALT': ['G']*k})
    dosage = pd.DataFrame(G, index=ids, columns=panel.ID)
    g.update(pheno={'BMI': pd.Series(y, index=ids, name='BMI'),
                    'SBP': pd.Series(sbp, index=ids, name='SBP'),
                    'DBP': pd.Series(dbp, index=ids, name='DBP')},
             covar={t: C.copy() for t in ('BMI', 'SBP', 'DBP')},
             age={t: pd.Series(age, index=ids) for t in ('BMI', 'SBP', 'DBP')},
             panel=panel, pfiles={1: Path('synthetic_chr1')})

    def independent_ols(yv, cv, gv):
        ok = np.isfinite(gv)
        X = np.column_stack([np.ones(ok.sum()), np.asarray(cv)[ok], gv[ok]])
        # lstsq/pinv with rank deficiency; GWAS adapter receives rank-reduced covariates.
        fit = np.linalg.lstsq(X, np.asarray(yv)[ok], rcond=None)[0]
        resid = np.asarray(yv)[ok]-X@fit
        rank = np.linalg.matrix_rank(X); df = len(resid)-rank
        se = np.sqrt((resid@resid/df)*np.linalg.pinv(X.T@X)[-1, -1])
        hc = np.linalg.pinv(X.T@X)@X.T
        hc_se = np.sqrt(np.sum((hc[-1]*resid)**2)*len(resid)/df)
        return fit[-1], se, 2*stats.t.sf(abs(fit[-1]/se), df), ok.sum(), hc_se

    commands = []

    def adapter(command, prefix):
        command = list(map(str, command)); commands.append(command)
        assert not set(command) & {'--clump', '--maf', '--mac'}
        def arg(flag): return command[command.index(flag)+1]
        keep = pd.read_csv(arg('--keep'), sep='\t')['#IID'].astype(str)
        selected = pd.read_csv(arg('--extract'), header=None)[0].astype(str)
        prefix = Path(prefix); prefix.parent.mkdir(parents=True, exist_ok=True)
        if '--export' in command:
            raw = pd.DataFrame({'FID': 0, 'IID': keep.to_numpy()})
            for snp in selected:
                raw[f'{snp}_G'] = dosage.loc[keep, snp].to_numpy()
            raw.to_csv(str(prefix)+'.raw', sep='\t', index=False, na_rep='NA')
            return
        cv = pd.read_csv(arg('--covar'), sep='\t').set_index('#IID').loc[keep]
        ph = pd.read_csv(arg('--pheno'), sep='\t').set_index('#IID').loc[keep, 'BMI']
        rows = []
        for snp in selected:
            b, se, p, obs, _ = independent_ols(ph, cv, dosage.loc[keep, snp].to_numpy())
            row = panel.set_index('ID').loc[snp].to_dict()
            # Mimic REF being the effect allele for one source SNP.
            flip = snp == 's1'
            row.update(ID=snp, A1='A' if flip else 'G', TEST='ADD', BETA=-b if flip else b,
                       SE=se, P=p, OBS_CT=obs, ERRCODE='.')
            rows.append(row)
        pd.DataFrame(rows).to_csv(str(prefix)+'.BMI.glm.linear', sep='\t', index=False, float_format='%.9g')

    with tempfile.TemporaryDirectory(prefix='agemr-b3-validation-') as tmp:
        out = Path(tmp)
        for folder in ('qc', 'design', 'source_GWAS', 'genotypes', 'associations', 'MR', 'figures'):
            (out/folder).mkdir()
        g.update(OUT=out, run_cmd=adapter)
        exec(compile(code[3], '<sampling>', 'exec'), g)
        # Target remains exactly the v2 target for the same inputs and seed.
        legacy_target = ids[target_mask]
        expected_target = pd.Index(g['stable_rng']('target').choice(
            legacy_target, min(g['N_TARGET_REQUESTED'], len(legacy_target)), replace=False)).sort_values()
        assert g['target_ids'].equals(expected_target)
        assert set(g['roles']) == {'target', 'source'}
        for s in g['source_sets']:
            assert not set(s['ids']) & set(g['target_ids'])
            assert len(s['ids']) == g['N_SOURCE']
        quotas = g['age_quotas']
        assert np.allclose(quotas.groupby('sid').realized_source_probability.sum(), 1)
        age_means = g['sample_qc'].groupby('scenario').age_mean.mean()
        assert age_means.younger_strong < age_means.younger_moderate < age_means.matched < age_means.older_moderate < age_means.older_strong
        exec(compile(code[5], '<GWAS selection>', 'exec'), g)
        assert g['selection'].N_selected.nunique() > 1, 'Fixture must exercise changing tool sets'
        for s in g['source_sets']:
            x = g['source_gwas'][s['sid']]
            assert set(s['iv_ids']) == set(x.loc[x.P.lt(g['IV_P_THRESHOLD']) & x.valid_GWAS, 'ID'])
        # Boundary checks: strict P, invalid P, errors and no implicit F threshold.
        boundary = g['source_gwas'][g['source_sets'][0]['sid']].iloc[:4].copy()
        boundary['ERRCODE'] = '.'; boundary['P'] = [1e-5, 1e-6, np.nan, 1e-12]
        boundary['F_classic'] = 0.; boundary.iloc[3, boundary.columns.get_loc('ERRCODE')] = 'ERR'
        _, chosen = g['select_source_instruments'](boundary, g['IV_P_THRESHOLD'])
        assert chosen.tolist() == [False, True, False, False]
        exec(compile(code[7], '<associations>', 'exec'), g)
        a = next(a for a in g['source_assoc'] if 's1' in set(a.ID))
        assert a.loc[a.ID.eq('s1'), 'orientation'].iloc[0] == -1
        for j, (s, a) in enumerate(zip(g['source_sets'], g['source_assoc'])):
            for row in a.itertuples(index=False):
                b, se, _, obs, hc = independent_ols(g['pheno']['BMI'].loc[s['ids']], C.loc[s['ids']], dosage.loc[s['ids'], row.ID].to_numpy())
                assert np.isclose(row.BETA_ALT, b, rtol=1e-8)
                assert np.isclose(row.SE_classic, se, rtol=1e-7)
                assert np.isclose(row.SE_HC1, hc, rtol=1e-7)
        covariance = g['SigmaX']
        assert np.linalg.eigvalsh(covariance).min() > -1e-10
        assert covariance.shape[0] == sum(len(a) for a in g['source_assoc'])
        # Independently reconstruct a cross-source HC1 influence covariance.
        vectors = []
        column_indices = []
        for j in (0, 3):
            s, a = g['source_sets'][j], g['source_assoc'][j]
            loc = int(np.flatnonzero(a.ID.eq('s0'))[0])
            gv = dosage.loc[s['ids'], 's0'].to_numpy(); ok = np.isfinite(gv)
            X = np.column_stack([np.ones(ok.sum()), C.loc[s['ids']].to_numpy()[ok], gv[ok]])
            inverse = np.linalg.pinv(X)
            residual = g['pheno']['BMI'].loc[s['ids']].to_numpy()[ok]-X@(inverse@g['pheno']['BMI'].loc[s['ids']].to_numpy()[ok])
            influence = inverse[-1]*residual*np.sqrt(ok.sum()/(ok.sum()-np.linalg.matrix_rank(X)))
            vector = pd.Series(0., index=ids)
            vector.loc[s['ids'][ok]] = influence*a.orientation.iloc[loc]
            vectors.append(vector.to_numpy())
            column_indices.append(int(g['offsets'][j])+loc)
        assert np.isclose(covariance[tuple(column_indices)], vectors[0]@vectors[1], rtol=1e-7)
        exec(compile(code[9], '<MR>', 'exec'), g)
        mr = g['mr']; assert len(mr) == 2*15*4
        assert g['contrasts'].repeat.eq('mean').sum() == 32
        assert (pd.read_csv(out/'qc/target_harmonization.tsv', sep='\t').N_invalid_target > 0).any()
        # Direct check of point estimates with each source's actual selected set.
        for s, a in zip(g['source_sets'], g['source_assoc']):
            target = g['target_assoc']['SBP'].set_index('ID')
            ok = target.loc[a.ID, 'status'].eq('ok').to_numpy()
            xx = a.BETA.to_numpy()[ok]
            yy = target.loc[a.ID, 'BETA'].to_numpy()[ok]*a.orientation.to_numpy()[ok]
            sy = target.loc[a.ID, 'SE_classic'].to_numpy()[ok]
            expected = np.sum(xx*yy/sy**2)/np.sum(xx**2/sy**2)
            actual = mr.loc[mr.sid.eq(s['sid']) & mr.outcome.eq('SBP') & mr.method.eq('IVW'), 'estimate'].iloc[0]
            assert np.isclose(actual, expected)
            sx = a.SE_classic.to_numpy()[ok]
            weights = sy**-2
            denominator = np.dot(weights, xx**2 - sx**2)
            expected_divw = np.dot(weights, xx*yy)/denominator
            expected_se = np.sqrt(np.dot(weights, xx**2) + expected_divw**2 *
                                  np.dot(weights**2, sx**2*(xx**2+sx**2)))/denominator
            divw = mr.loc[mr.sid.eq(s['sid']) & mr.outcome.eq('SBP') & mr.method.eq('dIVW')].iloc[0]
            assert np.isclose(divw.estimate, expected_divw)
            assert np.isclose(divw.SE_conventional, expected_se)
        # Empty/one/two/three SNPs and empty covariance must have explicit outcomes.
        for count in (0, 1, 2, 3):
            rows, _ = g['conventional_mr'](np.arange(1, count+1.), np.ones(count)*.1,
                                          np.arange(1, count+1.)*.5, np.ones(count)*.2)
            assert np.isfinite(rows['IVW']['estimate']) == (count > 0)
            assert np.isfinite(rows['dIVW']['estimate']) == (count > 0)
            assert np.isfinite(rows['weighted_median']['estimate']) == (count >= 3)
            assert np.isfinite(rows['MR_Egger']['estimate']) == (count >= 3)
        assert g['covariance_root'](np.empty((0, 0))).shape == (0, 0)
        weak_rows, _ = g['conventional_mr'](np.array([.1]), np.array([.1]),
                                             np.array([.2]), np.array([.1]))
        assert weak_rows['dIVW']['status'] == 'unstable_dIVW_denominator'
        assert np.isnan(weak_rows['dIVW']['estimate'])
        draws = g['mr_points'](np.array([[.2], [.05]]), np.array([[.1], [.1]]),
                              np.array([.2]), np.array([.1]))
        assert np.isclose(draws['dIVW'][0], 2/3) and np.isnan(draws['dIVW'][1])
        # A no-hit source alongside available sources must not shift covariance offsets.
        def one_no_hit(command, prefix):
            adapter(command, prefix)
            result = Path(str(prefix)+'.BMI.glm.linear')
            if Path(prefix).parent.name == 'r01_matched' and result.is_file():
                x = pd.read_csv(result, sep='\t'); x['P'] = 1.
                x.to_csv(result, sep='\t', index=False)
        g['run_cmd'] = one_no_hit
        for block in (5, 7, 9):
            exec(compile(code[block], f'<one no-hit source {block}>', 'exec'), g)
        assert g['mr'].loc[g['mr'].sid.eq('r01_matched'), 'estimate'].isna().all()
        assert g['mr'].loc[~g['mr'].sid.eq('r01_matched'), 'estimate'].notna().any()
        assert g['mean_mr'].loc[g['mean_mr'].scenario.eq('matched'), 'estimate_mean_repeats'].isna().all()
        # No-hit full path: all source GWAS reports become non-significant.
        normal_adapter = adapter
        def no_hits(command, prefix):
            normal_adapter(command, prefix)
            result = Path(str(prefix)+'.BMI.glm.linear')
            if result.is_file():
                x = pd.read_csv(result, sep='\t'); x['P'] = 1.
                x.to_csv(result, sep='\t', index=False)
        g['run_cmd'] = no_hits
        for block in (5, 7, 9):
            exec(compile(code[block], f'<no hits {block}>', 'exec'), g)
        assert g['K'] == 0 and g['mr'].estimate.isna().all()
        assert g['mean_mr'].estimate_mean_repeats.isna().all()
        assert not any('--clump' in command for command in commands)
    print('PASS: synthetic 15-source/2-outcome workflow, variable P-selected IVs, fixed target, ALT flips,')
    print('missing genotypes, target-invalid SNPs, independent OLS/HC1, IVW and dIVW checks, covariance,')
    print('strict P threshold, 0/1/2/3 tools, and full no-hit workflow. Real PLINK/AoU not executed.')


if __name__ == '__main__':
    validate()
