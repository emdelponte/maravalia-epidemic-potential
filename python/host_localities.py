import pandas as pd, numpy as np
from scipy import stats
R=pd.read_csv('epi_by_site_scenario.csv'); B=R[R.country=='BR'].pivot(index='site',columns='scen',values='EPI')
B['key']=B.S0_orig_BR.round(2)
H=pd.read_csv('../tables/brazilian_sites_epidemic_potential.csv')
H['key']=H.EPI_inf.round(2)
# nearest match on original EPI
H['cell']=[B.index[np.argmin(abs(B.S0_orig_BR.values-v))] for v in H.EPI_inf]
H['match_err']=[np.min(abs(B.S0_orig_BR.values-v)) for v in H.EPI_inf]
H=H.join(B.drop(columns='key'),on='cell')
print('max match error',H.match_err.max().round(3),'| unique cells',H.cell.nunique())
for sc in ['S0_orig_BR','S2_RH90_cap24','S4_RH90_P1_cap24']:
    x=H[sc]; a=np.log10(H['Infested Area (ha)'])
    r=stats.pearsonr(a,x); rho=stats.spearmanr(a,x)
    g=H.groupby('cell').agg(A=('Infested Area (ha)',lambda v:np.log10(v).mean()),E=(sc,'first')); rc=stats.spearmanr(g.A,g.E)
    print(f'{sc}: host-locality EPI mean {x.mean():.1f} range {x.min():.1f}-{x.max():.1f}; <20: {(x<20).sum()}/73; r(logArea)={r[0]:.3f} P={r[1]:.4f}; rho={rho[0]:.3f} P={rho[1]:.4f}; cell-level rho={rc[0]:.3f} P={rc[1]:.3f} (n={len(g)})')
print(H.sort_values('S2_RH90_cap24')[['Municipality','Infested Area (ha)','S0_orig_BR','S2_RH90_cap24','S4_RH90_P1_cap24']].round(1).to_string(index=False))
H.to_csv('host_localities_epi_scenarios.csv',index=False)
M=pd.read_csv('epi_monthly_scenario.csv'); mm=M.groupby(['country','scen','month']).EPI.mean().unstack().round(2)
print(mm.to_string())
