import nibabel as nib, numpy as np, json, sys, os
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
LEVELS=[0,2,4,6,8,10,12]
def analyse(path, save_qc=False):
    img=nib.load(path); a=np.asarray(img.dataobj).astype(np.uint8); A=img.affine.copy(); vx=np.abs(np.linalg.det(A[:3,:3]))
    keep=np.isin(a,[1,2,3,6]); bb=np.argwhere(keep); lo=np.maximum(bb.min(0)-20,0); hi=np.minimum(bb.max(0)+21,a.shape)
    a=a[lo[0]:hi[0],lo[1]:hi[1],lo[2]:hi[2]].copy(); T=np.eye(4); T[:3,3]=lo; A=A@T
    LV=a==3; AO=a==6; LA=a==2
    out={'id':os.path.basename(path).split('.')[0],'LV_ml':LV.sum()*vx/1000,'LA_ml':LA.sum()*vx/1000,'AO_ml':AO.sum()*vx/1000,'myo_ml':(a==1).sum()*vx/1000}
    # LV-aorta interface
    aod=ndi.binary_dilation(AO,structure=ndi.generate_binary_structure(3,1))
    itf=np.argwhere(LV&aod)
    if len(itf)<20: out['error']='no LV-aorta interface'; return out
    P=(A[:3,:3]@itf.T).T+A[:3,3]
    c=P.mean(0); w,v=np.linalg.eigh(np.cov((P-c).T)); n=v[:,0]
    # orient n from LV into aorta
    aop=(A[:3,:3]@np.argwhere(AO).T).T+A[:3,3]
    near=aop[np.linalg.norm(aop-c,axis=1)<15]
    if (near.mean(0)-c)@n<0: n=-n
    out.update(itf_area_mm2=float(len(itf)*vx**(2/3)), itf_planarity_sd=float(np.sqrt(w[0])), n=n.tolist(), c=c.tolist())
    e1=np.cross(n,[1,0,0]); 
    if np.linalg.norm(e1)<0.1: e1=np.cross(n,[0,1,0])
    e1/=np.linalg.norm(e1); e2=np.cross(n,e1)
    Ainv=np.linalg.inv(A); h=0.25; g=np.arange(-30,30+h,h); U,V=np.meshgrid(g,g,indexing='ij')
    LAd=ndi.distance_transform_edt(~LA,sampling=img.header.get_zooms()[:3]).astype(np.float32)  # mm to LA
    qc={}
    for d in LEVELS:
        o=c-d*n
        pts=o[None,None,:]+U[...,None]*e1+V[...,None]*e2
        ijk=(Ainv[:3,:3]@pts.reshape(-1,3).T).T+Ainv[:3,3]
        m=ndi.map_coordinates(LV.astype(np.uint8),ijk.T,order=0,mode='constant').reshape(U.shape)>0
        lab,nl=ndi.label(m)
        if nl==0: out[f'area_{d}']=np.nan; continue
        # component closest to plane centre
        cen=np.array(m.shape)//2; dist=ndi.distance_transform_edt(lab==0); 
        ids=lab[tuple(cen)] if lab[tuple(cen)]>0 else None
        if ids is None:
            cc=np.argwhere(lab>0); k=np.argmin(((cc-cen)**2).sum(1)); ids=lab[tuple(cc[k])]
        R=lab==ids; area=R.sum()*h*h
        q=np.argwhere(R)*h; q=q-q.mean(0); ev=np.sort(np.linalg.eigvalsh(np.cov(q.T)))[::-1]
        # contour fraction adjacent to LA (within 2 mm) -> mitral/inflow border flag
        edge=R&~ndi.binary_erosion(R); ep=np.argwhere(edge)
        epw=o+((ep[:,0]*h-30)[:,None]*e1+(ep[:,1]*h-30)[:,None]*e2)
        eijk=(Ainv[:3,:3]@epw.T).T+Ainv[:3,3]
        dLA=ndi.map_coordinates(LAd,eijk.T,order=1,mode='nearest')
        out[f'area_{d}']=float(area); out[f'ellip_{d}']=float(np.sqrt(ev[0]/ev[1])); out[f'Deq_{d}']=float(2*np.sqrt(area/np.pi))
        out[f'major_{d}']=float(4*np.sqrt(ev[0])); out[f'minor_{d}']=float(4*np.sqrt(ev[1])); out[f'LAborder_{d}']=float((dLA<2.5).mean())
        if save_qc: qc[d]=R
    if save_qc: return out,qc
    return out
if __name__=='__main__':
    res=[]
    for p in sys.argv[1:]:
        try: res.append(analyse(p))
        except Exception as e: res.append({'id':os.path.basename(p).split('.')[0],'error':repr(e)})
        print(res[-1]['id'],end=' ',flush=True)
    json.dump(res,open('lvot_cas_results.json','w'))
