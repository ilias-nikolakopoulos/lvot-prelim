"""Short-axis CT + label QC strips (depths -3, 2, 5, 8, 12 mm) for case screening.
Run from the imagecas folder:  uv run --python 3.12 --with nibabel --with scipy --with numpy --with matplotlib python cfd/qc/shortaxis_qc.py <out.png> <id,id,...>
Used: systolic_shortaxis_1 = 13,33,146,166,253,364,370; systolic_shortaxis_2 = 388,539,576,586,592,808,825;
diastolic_shortaxis_1 = 142,479,29,587,467,663; diastolic_shortaxis_2 = 167,978"""
import sys, os, glob, numpy as np, nibabel as nib, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import ndimage as ndi
sys.path.insert(0, os.getcwd()); import lvot_curtain as lc
ids=[int(x) for x in sys.argv[2].split(',')]; outp=sys.argv[1]
where={int(os.path.basename(p).split('.')[0]):os.path.dirname(p) for p in glob.glob("ImageCAS_*/*/*.img.nii.gz")}
lab="labels_DTU/ImageCAS-STACOM2025-02-10-2025/segmentations"; deps=[-3,2,5,8,12]
f,axs=plt.subplots(len(ids),len(deps),figsize=(2.4*len(deps),2.5*len(ids)))
g=np.arange(-20,20.01,0.25); U,V=np.meshgrid(g,g,indexing='ij')
for r,cid in enumerate(ids):
    ln=nib.load(f"{lab}/{cid}.img.nii.gz"); L=np.asarray(ln.dataobj).astype(np.uint8); A=ln.affine; Ai=np.linalg.inv(A)
    I=np.asarray(nib.load(f"{where[cid]}/{cid}.img.nii.gz").dataobj).astype(np.float32)
    c,n,e1,e2=lc.frame(L,A); bl=np.median(I[L==3]); my=np.median(I[L==1])
    sizes=np.bincount(ndi.label(L==3)[0].ravel())[1:]
    for k,dep in enumerate(deps):
        pts=(c-dep*n)+U[...,None]*e1+V[...,None]*e2; ijk=(Ai[:3,:3]@pts.reshape(-1,3).T).T+Ai[:3,3]
        im=ndi.map_coordinates(I,ijk.T,order=1).reshape(U.shape); lb=ndi.map_coordinates(L,ijk.T,order=0).reshape(U.shape)
        a=axs[r,k]; a.imshow(im.T,cmap='gray',vmin=my-150,vmax=bl+150,origin='lower',extent=[-20,20,-20,20])
        for val,col in [(3,'cyan'),(6,'yellow'),(2,'magenta'),(1,'orange')]:
            if (lb==val).any(): a.contour(g,g,(lb==val).T.astype(float),[0.5],colors=col,linewidths=0.7)
        a.set_xticks([]); a.set_yticks([]); a.set_title(f'{cid}  {dep:+d} mm'+(f'  LVcc={sizes.tolist()[:3]}' if k==0 and len(sizes)>1 else ''),fontsize=7)
    print(cid,flush=True)
f.text(0.01,0.995,'cyan LV, yellow aorta, magenta LA, orange myo; x=e1 (toward LA, AP), y=e2; depth below valve plane',fontsize=8,va='top')
plt.tight_layout(rect=(0,0,1,0.99)); plt.savefig(outp,dpi=70)
