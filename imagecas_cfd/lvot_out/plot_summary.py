import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=pd.read_csv("lvot_out/lvot_curtain_results.csv"); f=d.curtain_end_mm.notna(); s=d[f]
SURF,INK,INK2,GRID='#fcfcfb','#0b0b0b','#52514e','#e4e3df'
B1,B2,B3='#86b6ef','#3987e5','#184f95'; BLUE,ORANGE,GRAY='#2a78d6','#eb6834','#b4b2ab'
plt.rcParams.update({'font.size':10,'axes.edgecolor':GRAY,'axes.labelcolor':INK2,'xtick.color':INK2,'ytick.color':INK2,
 'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.color':GRID,'grid.linewidth':0.6,'axes.axisbelow':True,
 'axes.facecolor':SURF,'figure.facecolor':SURF,'axes.titlesize':11,'axes.titleweight':'bold','axes.titlelocation':'left','text.color':INK})
fig,ax=plt.subplots(2,3,figsize=(15,8.6))
# A histogram
a=ax[0,0]; a.hist(s.curtain_end_mm,bins=np.arange(2.5,12.75,0.5),color=BLUE,edgecolor=SURF,linewidth=2)
m=s.curtain_end_mm.median(); a.axvline(m,color=INK,lw=1,ls='--'); a.text(9.3,a.get_ylim()[1]*0.93,f'median {m:.1f} mm\nIQR {s.curtain_end_mm.quantile(.25):.1f}–{s.curtain_end_mm.quantile(.75):.1f}',va='top',fontsize=9,color=INK)
a.set(title=f'A  Curtain end depth (n={f.sum()})',xlabel='mm below aortic valve plane',ylabel='cases'); a.grid(axis='x',visible=False)
# B outcomes
a=ax[0,1]; r=d.curtain_rule
cats=[('5 planes',r.eq('5 planes')),('4 planes',r.eq('4 planes')),('3 planes',r.eq('3 planes')),('2 planes agree',r.eq('2 planes agree')),
 ('2 planes disagree',r.str.startswith('not found (2')),('1 plane',r.str.startswith('not found (1')),('0 planes',r.str.startswith('not found (0'))]
lab=[c[0] for c in cats]; n=[c[1].sum() for c in cats]; col=[BLUE]*4+[GRAY]*3
y=np.arange(len(cats))[::-1]; a.barh(y,n,color=col,height=0.7)
for yi,ni in zip(y,n): a.text(ni+4,yi,f'{ni}',va='center',fontsize=9,color=INK2)
a.set_yticks(y,lab); a.grid(axis='y',visible=False); a.set_xlabel('cases')
a.set_title(f'B  Rule outcome: found {f.sum()} / {len(d)} ({f.mean():.0%})')
from matplotlib.patches import Patch; a.legend(handles=[Patch(color=BLUE,label='found'),Patch(color=GRAY,label='not found')],loc='lower right',frameon=False,fontsize=9)
# C found rate by LV/myo
a=ax[0,2]; bins=[0,.6,.8,1,1.2,5]; labs=['<0.6\n(systolic flag)','0.6–0.8','0.8–1.0','1.0–1.2','>1.2']
g=d.groupby(pd.cut(d.LV_myo_ratio,bins),observed=False).curtain_end_mm; rate=g.apply(lambda x:x.notna().mean()*100); cnt=g.size()
x=np.arange(5); a.bar(x,rate,color=[ORANGE]+[BLUE]*4,width=0.65)
for xi,(rt,c) in enumerate(zip(rate,cnt)): a.text(xi,rt+2,f'{rt:.0f}%\nn={c}',ha='center',fontsize=8.5,color=INK2)
a.set_xticks(x,labs,fontsize=8.5); a.set(ylim=(0,100),ylabel='% with curtain end found',xlabel='LV / myocardium volume ratio',title='C  Detection vs cardiac phase proxy'); a.grid(axis='x',visible=False)
# D-F by level (found cases)
lev=['2mm','mid','curtain_end']; xl=['2 mm','mid','curtain\nend']
def box(a,cols,title,yl):
    bp=a.boxplot([s[c] for c in cols],positions=range(3),widths=0.5,patch_artist=True,showfliers=False,medianprops=dict(color=INK,lw=1.5),whiskerprops=dict(color=GRAY),capprops=dict(color=GRAY))
    for p,c in zip(bp['boxes'],[B1,B2,B3]): p.set(facecolor=c,edgecolor=c)
    for i,c in enumerate(cols): a.text(i+0.3,s[c].median(),f'{s[c].median():.2f}' if s[c].median()<5 else f'{s[c].median():.0f}',va='center',fontsize=9,color=INK)
    a.set_xticks(range(3),xl); a.set(title=title,ylabel=yl); a.grid(axis='x',visible=False)
box(ax[1,0],[f'ellipticity_{l}' for l in lev],f'D  Ellipticity by depth (n={len(s)})','ellipticity (√ eigenvalue ratio)')
box(ax[1,2],[f'area_mm2_{l}' for l in lev],f'F  Cross-sectional area by depth (n={len(s)})','area (mm²)')
a=ax[1,1]; xs=np.arange(3)
for key,c,name,dx in [('major',BLUE,'major axis',-0.04),('minor',ORANGE,'minor axis',0.04)]:
    q=np.array([s[f'{key}_mm_{l}'].quantile([.25,.5,.75]).values for l in lev])
    a.fill_between(xs,q[:,0],q[:,2],color=c,alpha=0.15,lw=0); a.plot(xs,q[:,1],color=c,lw=2,marker='o',ms=8,mec=SURF,mew=2)
    a.text(2.08,q[2,1],f'{name}\n{q[2,1]:.1f} mm',color=INK,va='center',fontsize=9)
    a.text(-0.08,q[0,1],f'{q[0,1]:.1f}',color=INK2,va='center',ha='right',fontsize=9)
a.set_xticks(xs,xl); a.set_xlim(-0.4,2.6); a.set(title='E  Major vs minor axis (median, IQR band)',ylabel='mm'); a.grid(axis='x',visible=False)
fig.suptitle('LVOT on ImageCAS (DTU labels): end of LVOT = end of aorto-mitral curtain',x=0.01,ha='left',fontsize=13,fontweight='bold')
fig.text(0.01,0.935,f'{len(d)} cases (90, 141 excluded). Panels D–F: cases with a curtain end found (n={len(s)}); levels are 2 mm, half of curtain-end depth, and curtain end.',fontsize=9.5,color=INK2)
plt.tight_layout(rect=(0,0,1,0.93)); plt.savefig('lvot_out/lvot_curtain_summary.png',dpi=150)
