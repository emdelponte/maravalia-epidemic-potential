# Recompute EPI_inf under alternative, consistent wetness definitions (review, 2026-09-28)
import glob, numpy as np, pandas as pd, pyreadr
from scipy import stats
To,s,g,lam,k=21.259719067938295,2.2797747024024932,0.2983713151869794,8.366570952838913,8.928832606137053
Wb=12.5
def I(T,W,cap=None):
    Ws=W if cap is None else np.minimum(W,cap)
    sw=np.maximum(1.0,s+g*(Ws-Wb))
    out=np.exp(-.5*((T-To)/sw)**2)*(1-np.exp(-(W/lam)**k))
    return np.where(W<6,0,out)
def episodes(wet,T,split=None):
    w=wet.astype(np.int8); d=np.diff(np.r_[0,w,0]); st=np.where(d==1)[0]; en=np.where(d==-1)[0]
    if split:
        S,E=[],[]
        for a,b in zip(st,en):
            for x in range(a,b,split): S.append(x);E.append(min(x+split,b))
        st,en=np.array(S,int),np.array(E,int)
    cs=np.r_[0,np.cumsum(T)]; dur=en-st; mT=(cs[en]-cs[st])/np.maximum(dur,1)
    return st,dur,mT
SC={ # name: (wet rule, cap, split)
 'S0_orig_BR':(lambda d:(d.RH2M>=90)|(d.PRECTOTCORR>0),None,None),
 'S1_RH90':(lambda d:d.RH2M>=90,None,None),
 'S2_RH90_cap24':(lambda d:d.RH2M>=90,24,None),
 'S3_RH90_split24':(lambda d:d.RH2M>=90,24,24),
 'S4_RH90_P1_cap24':(lambda d:(d.RH2M>=90)|(d.PRECTOTCORR>=1),24,None),
}
rows=[];mon=[]
for country,pat,idc in [('BR','../nasapower_maravalia/hourly/*.rds','cell_id'),('AU','../nasapower_australia/hourly/*.rds','site_name')]:
    for f in sorted(glob.glob(pat)):
        d=pyreadr.read_r(f)[None]; T=d.T2M.values; site=d[idc].iloc[0]
        for sc,(rule,cap,split) in SC.items():
            wet=rule(d).values; st,dur,mT=episodes(wet,T,split); ip=I(mT,dur,cap)
            ny=d.YEAR.nunique()
            rows.append(dict(country=country,site=site,scen=sc,wet_frac=wet.mean(),eps_per_yr=len(dur)/ny,
                 mean_dur=dur.mean(),share_dur_gt24=(dur>24).mean(),EPI=ip.sum()/ny,
                 EPI_from_gt24=ip[dur>24].sum()/max(ip.sum(),1e-9)))
            m=d.MO.values[st]
            for mo in range(1,13): mon.append(dict(country=country,site=site,scen=sc,month=mo,EPI=ip[m==mo].sum()/ny))
R=pd.DataFrame(rows);M=pd.DataFrame(mon)
R.to_csv('epi_by_site_scenario.csv',index=False);M.to_csv('epi_monthly_scenario.csv',index=False)
summ=R.groupby(['country','scen']).agg(wet_frac=('wet_frac','mean'),eps_yr=('eps_per_yr','mean'),mean_dur=('mean_dur','mean'),
   share_gt24=('share_dur_gt24','mean'),EPI_from_gt24=('EPI_from_gt24','mean'),EPI_mean=('EPI','mean'),EPI_min=('EPI','min'),EPI_max=('EPI','max')).round(3)
print(summ.to_string())
mm=M.groupby(['country','scen','month']).EPI.mean().unstack()
print('\nMay-Jul share (BR):');print((mm.loc['BR'][[5,6,7]].sum(1)/mm.loc['BR'].sum(1)).round(3).to_string())
print('\nJan-Mar share (AU):');print((mm.loc['AU'][[1,2,3]].sum(1)/mm.loc['AU'].sum(1)).round(3).to_string())
print('\nBR peak months:');print(mm.loc['BR'].idxmax(1).to_string())
hi=['McLeod River','Laura Lea','Lakefield NP','Charters Towers','Inkerman']
print('\nAustralia per site:');au=R[R.country=='AU'].pivot(index='site',columns='scen',values='EPI').round(1);print(au.sort_values('S2_RH90_cap24',ascending=False).to_string())
for sc in SC:
    a=R[(R.country=='AU')&(R.scen==sc)]; h=a[a.site.isin(hi)].EPI; l=a[~a.site.isin(hi)].EPI
    t=stats.ttest_ind(h,l,equal_var=False); u=stats.mannwhitneyu(h,l)
    b=R[(R.country=='BR')&(R.scen==sc)].EPI
    print(f'{sc}: Welch P={t.pvalue:.3f} MWU P={u.pvalue:.3f} | BR range {b.min():.1f}-{b.max():.1f} vs AU high-impact {h.min():.1f}-{h.max():.1f}; BR cells above AU median {np.median(a.EPI):.1f}: {(b>np.median(a.EPI)).mean():.0%}')
