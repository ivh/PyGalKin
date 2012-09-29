"""
plot.py
"""

from pylab import *
from PyGalKin import *
import PyGalKin.db as DB
import pylab as P # to make P not self-reference to plot.py
import matplotlib
import matplotlib.mlab as mlab
from matplotlib import rcParams, colors
from mpl_toolkits.axes_grid import AxesGrid
import aplpy

from DjCigale.table import models as M

#matplotlib.use('Agg')
#matplotlib.rc('text', usetex = True)

######### SAURON COLORMAP ##############
#
# Original format:
#x = [1.0, 43.5, 86.0, 86.0+20, 128.5-10, 128.5, 128.5+10, 171.0-20, 171.0, 213.5, 256.0]
#red =   [0.01, 0.0, 0.4,  0.5, 0.3, 0.0, 0.7, 1.0, 1.0,  1.0, 0.9]
#green = [0.01, 0.0, 0.85, 1.0, 1.0, 0.9, 1.0, 1.0, 0.85, 0.0, 0.9]
#blue =  [0.01, 1.0, 1.0,  1.0, 0.7, 0.0, 0.0, 0.0, 0.0,  0.0, 0.9]

_sauron_data = {
    'red': ((0.0,0.01,0.01),
            (0.169921875,0.0,0.0),
            (0.3359375,0.4,0.4),
            (0.4140625,0.5,0.5),
            (0.462890625,0.3,0.3),
            (0.501953125,0.0,0.0),
            (0.541015625,0.7,0.7),
            (0.58984375,1.0,1.0),
            (0.66796875,1.0,1.0),
            (0.833984375,1.0,1.0),
            (1.0,0.9,0.9)),
    'green':((0.0,0.01,0.01),
             (0.169921875,0.0,0.0),
             (0.3359375,0.85,0.85),
             (0.4140625,1.0,1.0),
             (0.462890625,1.0,1.0),
             (0.501953125,0.9,0.9),
             (0.541015625,1.0,1.0),
             (0.58984375,1.0,1.0),
             (0.66796875,0.85,0.85),
             (0.833984375,0.0,0.0),
             (1.0,0.9,0.9)),
    'blue':((0.0,0.01,0.01),
            (0.169921875,1.0,1.0),
            (0.3359375,1.0,1.0),
            (0.4140625,1.0,1.0),
            (0.462890625,0.7,0.7),
            (0.501953125,0.0,0.0),
            (0.541015625,0.0,0.0),
            (0.58984375,0.0,0.0),
            (0.66796875,0.0,0.0),
            (0.833984375,0.0,0.0),
            (1.0,0.9,0.9))
    }
_sauron_inv_data = {
    'red': ((0.0,1.0,1.0),
            (0.169921875,1.0,1.0),
            (0.3359375,1.0,1.0),
            (0.4140625,0.4,0.4),
            (0.462890625,0.7,0.7),
            (0.501953125,0.0,0.0),
            (0.541015625,0.3,0.3),
            (0.58984375,0.5,0.5),
            (0.66796875,0.0,0.0),
            (0.833984375,1.0,1.0),
            (1.0,0.01,0.01)),
    'green':((0.0,1.0,1.0),
             (0.169921875,0.0,0.0),
             (0.3359375,0.85,0.85),
             (0.4140625,1.0,1.0),
             (0.462890625,1.0,1.0),
             (0.501953125,0.9,0.9),
             (0.541015625,1.0,1.0),
             (0.58984375,1.0,1.0),
             (0.66796875,0.85,0.85),
             (0.833984375,0.0,0.0),
             (1.0,0.01,0.01)),
    'blue':((0.0,1.0,1.0),
            (0.169921875,0.0,0.0),
            (0.3359375,0.0,0.0),
            (0.4140625,0.0,0.0),
            (0.462890625,0.0,0.0),
            (0.501953125,0.0,0.0),
            (0.541015625,0.7,0.7),
            (0.58984375,1.0,1.0),
            (0.66796875,1.0,1.0),
            (0.833984375,1.0,1.0),
            (1.0,0.01,0.01))
    }

sauron=colors.LinearSegmentedColormap('sauron', _sauron_data)
sauron_inv=colors.LinearSegmentedColormap('sauron_inv', _sauron_inv_data)
register_cmap(cmap=sauron)
register_cmap(cmap=sauron_inv)
rcParams['image.cmap']='sauron'

def imshow(*args, **kwargs):
    kwargs.setdefault('interpolation','nearest')
    kwargs.setdefault('origin','lower')
    kwargs.setdefault('cmap','sauron')
    args = list(args)
    args[0] = N.transpose(args[0])
    P.imshow(*args, **kwargs)

def plot(*args, **kwargs):
    kwargs.setdefault('drawstyle','steps-mid')
    P.plot(*args, **kwargs)

def noticks():
  P.setp(P.gca(),'xticks',[])
  P.setp(P.gca(),'yticks',[])

def plotscale(p,left=0.65,right=0.9,frombottom=0.1,kpc=False):
  axi=N.array(P.gca().axis())
  y=axi[2]+frombottom*(axi[3]-axi[2])
  x1=axi[0]+left*(axi[1]-axi[0])
  x2=axi[0]+right*(axi[1]-axi[0])
  P.plot([x1,x2],[y,y],'k-')
  arcs=(x2-x1)*p['echelle']
  if not kpc:
      kpc=arcsec2kpc(arcs,p['vr0'])
  P.text(x1-0.01,y+0.1,r'$%.1f^{\prime\prime}\approx %.1fkpc$'%(arcs,kpc),fontsize=9)


def setYaxis_pc(inarr):
  """ set the plotting-axes to parsec """

  nticks=inarr.ny()*inarr.scale()/1000
  nticks=nticks.sum()
  factor=1.0

  while nticks > 10:
    nticks/=2
    factor/=2

  P.setp(P.gca(),'yticks',str(N.arange(nticks)/inarr.scale()*1000/factor))
  P.setp(P.gca(),'yticklabels',str(N.arange(nticks)/factor))
  P.ylabel('[kpc]',font)

def setXaxis_pc(inarr):
  """ set the plotting-axes to parsec """

  nticks=inarr.nx()*inarr.scale()/1000
  nticks=nticks[0]
  factor=1.

  while nticks > 10:
    nticks/=2
    factor/=2

  P.setp(P.gca(),'xticks',N.arange(nticks)/inarr.scale()*1000/factor)
  P.setp(P.gca(),'xticklabels',N.arange(nticks)/factor)
  P.xlabel('[kpc]',font)


def showsum(data,vmin=1E5,vmax=2E6,range='cat',Z=1.002912,typ='sum'):
    if typ=='sum': dat=selsum(data,range=range,axis=2)
    elif typ=='aver': dat=selav(data,range=range,axis=2)
    P.imshow(N.transpose(dat),origin='lower',interpolation='nearest',vmin=vmin,vmax=vmax)


def plotspec(data,region=None,plotlines=False,Z=1.0206,style=False,linestyle='steps',vminmax=None):
    from PyArgus import SpecLen,Lamb0,Step
    if style:
        P.plot((N.arange(SpecLen)*Step)+Lamb0,data,style,linestyle=linestyle)
    else:
        P.plot((N.arange(SpecLen)*Step)+Lamb0,data,linestyle=linestyle)

    if plotlines == True:
        CaTz= CaT * Z
        Paschenz=Paschen * Z
        EmissionLinesz=EmissionLines * Z
        for i in N.arange(len(CaTz)): P.plot([CaTz[i],CaTz[i]],P.axis()[2:],'r')
        for i in N.arange(len(Paschenz)): P.plot([Paschenz[i],Paschenz[i]],P.axis()[2:],'k')
        for i in N.arange(len(EmissionLinesz)): P.plot([EmissionLinesz[i],EmissionLinesz[i]],P.axis()[2:],'b')

    if region == 'cat': # legacy
        region=[8470,8700]
    if region != None:
        relevant=data[lamb2pix(region[0]*Z,Lamb0,Step):lamb2pix(region[1]*Z,Lamb0,Step)]
        if vminmax=='sigbased':
            vmax=relevant.mean() + 4*relevant.std()
            vmin=relevant.mean() - 4*relevant.std()
        else: vmin,vmax=relevant.min(),relevant.max()
        P.axis([float(region[0]*Z),float(region[1]*Z),vmin,vmax])

    else: pass


def fillederrorplot(x,y,e1,e2=None,f='r--',c='r',alpha=0.5,label=None):
        if e2==None: e2=e1/2.0; e1=e2
        P.plot(x,y,f,linewidth=2,label=label)
        P.plot(x,y-e1,c+'-',linewidth=0.5)
        P.plot(x,y+e2,c+'-',linewidth=0.5)
        P.fill(N.concatenate((x,x[::-1])),N.concatenate((y+e2,(y-e1)[::-1])),c,alpha=alpha)




##################################

def velSigmaPlot():
    gs1=M.Galax.objects.filter(sample=1,maxvel__isnull=False)
    gs2=M.Galax.objects.filter(sample=2,maxvel__isnull=False)
    P.plot([9.,59.],[9.,59.],'k--')
    v,s,inc=zip(*gs1.values_list('maxvel','sigma_cent','incl'))
    v=array(v)/2
    P.plot(s*sq3,v,'ok',label=r'$v_{max,1}$',mfc='w')
    P.plot(s*sq3,v/N.sin(N.radians(inc)),'kD',label='$v_{max,i,1}$',mfc='w')
    v,s,inc=zip(*gs2.values_list('maxvel','sigma_cent','incl'))
    v=array(v)/2
    P.plot(s*sq3,v,'ok',label=r'$v_{max,2}$',mfc='k')
    P.plot(s*sq3,v/N.sin(N.radians(inc)),'kD',label='$v_{max,i,2}$',mfc='k')

    P.ylabel(r'$v_{max}\quad (km\,s^{-1})$')
    P.xlabel(r'$\sigma^\ast_{c}\quad (km\,s^{-1})$')
    P.grid()
    P.axis((9,60,0,185))
    P.legend(loc='upper left')

def rotMassRandPlot():
    gs1=M.Galax.objects.filter(sample=1,mass_p2p__isnull=False)
    gs2=M.Galax.objects.filter(sample=2,mass_p2p__isnull=False)
    P.loglog([3E7,2E10],[3E7,2E10],'k--',label='$\mathrm{unity}$')
    ms,mv,inc=zip(*gs1.values_list('mass_sig','mass_p2p','incl'))
    P.loglog(ms,mv,'ok',mfc='w',label='$sample\, 1$')
    ms,mv,inc=zip(*gs2.values_list('mass_sig','mass_p2p','incl'))
    P.loglog(ms,mv,'ok',mfc='k',label='$sample\, 2$')
    P.ylabel(r'$M_{rot}\quad (M_\odot)$')
    P.xlabel(r'$M_{disp}\quad (M_\odot)$')
    P.axis((6E7,3E10,1E6,1E11))
    P.grid(True)
    P.legend(loc='lower right')

def masshistPlot():
    mfr1=array([g.mass_p2p/g.mass_sig for g in M.Galax.objects.filter(sample=1) ])
    mfr2=array([g.mass_p2p/g.mass_sig for g in M.Galax.objects.filter(sample=2) ])
    P.hist(log10(mfr1),bins=arange(-2.5,2.5,0.5),color='w',alpha=0.5,lw=2,hatch='--')
    P.hist(log10(mfr2),bins=arange(-2.5,2.5,0.5),color='w',alpha=0.5,lw=2,hatch='//')
    P.ylabel('$N$')
    P.xlabel('$\log (M_{rot}/M_{disp})$')

def photMassRandPlot():
    gs=M.Galax.objects.all()

    fus = [dynMassDisk, dynMassSphere]
    rs = ['reff','rholm','h25','h27']
    fac = [1,1,1.68,1.68]
    fmts = ['o','D']
    cols = ['b','r','g','y']
    lab1 = ['disk + ','sphere + ']
    for i,r in enumerate(rs):
        for j,fu in enumerate(fus):
            for g in gs:
                rad=getattr(g,r) or 0
                rad*=fac[i]
                mass = fu(arcsec2kpc(rad,g.vsys),g.sigma_cent or 0) or None
                P.loglog(g.mass_phot,mass,cols[i]+fmts[j])


    P.xlabel(r'$M_{phot}\quad (M_\odot)$')
    P.ylabel(r'$M_{dispersion}\quad (M_\odot)$')
    P.grid()

def sigmasPlot():
    gs1=M.Galax.objects.filter(sample=1)
    gs2=M.Galax.objects.filter(sample=2)
    P.plot([7,65],[7,65],'k--',label='$\mathrm{unity}$')
    s1,s2,s3=zip(*gs1.values_list('sigma_cent','sigma_weighmean','sigma_collap'))
    P.plot(s1,s2,'o',label=r'$\sigma_{i,1}$',mfc='w')
    P.plot(s1,s3,'D',label=r'$\langle\sigma\rangle_{i,1}$',mfc='w')
    s1,s2,s3=zip(*gs2.values_list('sigma_cent','sigma_weighmean','sigma_collap'))
    P.plot(s1,s2,'o',label=r'$\sigma_{i,2}$',mfc='k')
    P.plot(s1,s3,'D',label=r'$\langle\sigma\rangle_{i,2}$',mfc='k')
    P.xlabel(r'$\sigma_{c} (km\,s^{-1})$')
    P.ylabel(r'$\sigma\, (km\,s^{-1})$')
    P.legend(loc='upper left')
    P.axis((5,70,5,70))
    P.grid()

def LsigmaPlot():
    gs1=M.Galax.objects.filter(sample=1)
    gs2=M.Galax.objects.filter(sample=2)
    #sig,mb=zip(*gs1.values_list('sigma_cent','mb'))
    HaLum,sig=zip(*[(g.p.get('HaLum'),g.sigma_cent) for g in gs1])
    HaLum=N.log10(N.array(HaLum).astype('Float64'))
    sig=N.log10(N.array(sig) *sq3)
    P.plot(sig,HaLum,'o',label=r'$\sigma_{1}',mfc='w')

    HaLum,sig=zip(*[(g.p.get('HaLum'),g.sigma_cent) for g in gs2])
    HaLum=N.log10(N.array(HaLum).astype('Float64'))
    sig=N.log10(N.array(sig)*sq3)
    P.plot(sig,HaLum,'o',label=r'$\sigma_{2}',mfc='k')

#    P.grid()
#    P.axis((5E,-20.7,0,65))
    P.ylabel(r'$\log L(H\alpha)\,\,(erg\,s^{-1})$')
    P.xlabel(r'$\log \sigma^\ast_c\, (km\,s^{-1})$')

def MsigmaPlot():
    gs1=M.Galax.objects.filter(sample=1)
    gs2=M.Galax.objects.filter(sample=2)
    sig,mb=zip(*gs1.values_list('sigma_cent','mb'))
    sig=N.log10(N.array(sig) *sq3)
    P.plot(mb,sig,'o',label=r'$\sigma_{1}',mfc='w')

    sig,mb=zip(*gs2.values_list('sigma_cent','mb'))
    sig=N.log10(N.array(sig) *sq3)
    P.plot(mb,sig,'o',label=r'$\sigma_{2}',mfc='k')

    P.grid()
    P.axis((-13.2,-20.7,1.21,2.09))
    P.xlabel(r'$M_B$')
    P.ylabel(r'$\log \sigma_c^\ast\, (km\,s^{-1})$')

def mu_h_mergPlot():
    gs1=M.Galax.objects.filter(sample=1)
    gs2=M.Galax.objects.filter(sample=2)

    h1=N.array([arcsec2kpc(g.h27 or g.h25 or nan,g.vsys)[0] or nan for g in gs1])
    mu1=N.array([g.p.get('mu0host') or nan for g in gs1])
    merg1=N.array([g.doub+g.irrvf for g in gs1])

    P.scatter(log10(h1),array(mu1),s=(merg1+1)*40,c=merg1,cmap=cm.bone_r,alpha=0.8,lw=1.6,marker='D')

    h2=N.array([arcsec2kpc(g.h27 or g.h25 or nan,g.vsys)[0] or nan for g in gs2])
    mu2=N.array([g.p['mu0host'] or nan for g in gs2])
    merg2=N.array([g.doub+g.irrvf for g in gs2])

    P.scatter(log10(h2),array(mu2),s=(merg2+1)*40,c=merg2,cmap=cm.bone_r,alpha=0.8,lw=1.6,marker='o')
    P.axis((-0.75,0.9,25.8,18.1))

    P.xlabel(r'$log_{10}(h_r)$')
    P.ylabel(r'$\mu_{0,host}$')

def M_burst_mergPlot():
    gs1=M.Galax.objects.filter(sample=1)
    gs2=M.Galax.objects.filter(sample=2)

    burst1=N.array([g.p.get('bfrac27') or g.p.get('bfrac25') or nan for g in gs1])
    mb1=N.array([g.mb or nan for g in gs1])
    merg1=N.array([g.doub+g.irrvf for g in gs1])

    P.scatter(mb1,burst1,s=(merg1+1)*40,c=merg1,cmap=cm.bone_r,alpha=0.8,lw=1.4,marker='D')

    burst2=N.array([g.p.get('bfrac27') or g.p.get('bfrac25') or nan for g in gs2])
    mb2=N.array([g.mb or nan for g in gs2])
    merg2=N.array([g.doub+g.irrvf for g in gs2])

    P.scatter(mb2,burst2,s=(merg2+1)*40,c=merg2,cmap=cm.bone_r,alpha=0.8,lw=1.4,marker='o')

    P.axis((-13,-20.5,0,1))
    P.xlabel(r'$M_B$')
    P.ylabel(r'$\mathrm{burst\, fraction}$')

def EWsigmaPlot():
    gs=M.Galax.objects.filter(sample=1)
    merg=N.array([g.doub+g.irrvf for g in gs])
    ew=N.array([g.p.get('HaEW') for g in gs],dtype='Float64')
    sig=N.array([g.sigma_cent for g in gs]) * sq3
    P.scatter(N.log10(sig),N.log10(ew),s=(merg+1)*40,c=merg,cmap=cm.bone_r,alpha=0.8,lw=1.4,marker='D')

    gs=M.Galax.objects.filter(sample=2)
    merg=N.array([g.doub+g.irrvf for g in gs])
    ew=N.array([g.p.get('HaEW') for g in gs],dtype='Float64')
    sig=N.array([g.sigma_cent for g in gs]) * sq3
    P.scatter(N.log10(sig),N.log10(ew),s=(merg+1)*40,c=merg,cmap=cm.bone_r,alpha=0.8,lw=1.4,marker='o')

    P.xlabel(r'$\log \sigma_c^\ast_c\, (km\,s^{-1})$')
    P.ylabel(r'$\log \mathrm{EW}(H\alpha)$')

def L_EWPlot():
    gs=M.Galax.objects.filter(sample=1)
    merg=N.array([g.doub+g.irrvf for g in gs])
    ew=N.array([g.p.get('HaEW') for g in gs],dtype='Float64')
    LHa=N.array([g.p.get('HaLum') for g in gs],dtype='Float64')
    P.scatter(N.log10(ew),N.log10(LHa),s=(merg+1)*40,c=merg,cmap=cm.bone_r,alpha=0.8,lw=1.4,marker='D')

    gs=M.Galax.objects.filter(sample=2)
    merg=N.array([g.doub+g.irrvf for g in gs])
    ew=N.array([g.p.get('HaEW') for g in gs],dtype='Float64')
    LHa=N.array([g.p.get('HaLum') for g in gs],dtype='Float64')
    P.scatter(N.log10(ew),N.log10(LHa),s=(merg+1)*40,c=merg,cmap=cm.bone_r,alpha=0.8,lw=1.4,marker='o')

    P.ylabel(r'$\log L(H\alpha)\,\,(erg\,s^{-1})$')
    P.xlabel(r'$\log \mathrm{EW}(H\alpha)$')

def EW_burstPlot():
    gs=M.Galax.objects.filter(sample=1)
    merg=N.array([g.doub+g.irrvf for g in gs])
    ew=N.array([g.p.get('HaEW') for g in gs],dtype='Float64')
    burst=N.array([g.p.get('bfrac27') or g.p.get('bfrac25') or nan for g in gs])
    P.scatter(N.log10(ew),burst,s=(merg+1)*40,c=merg,cmap=cm.bone_r,alpha=0.8,lw=1.4,marker='D')

    gs=M.Galax.objects.filter(sample=2)
    merg=N.array([g.doub+g.irrvf for g in gs])
    ew=N.array([g.p.get('HaEW') for g in gs],dtype='Float64')
    burst=N.array([g.p.get('bfrac27') or g.p.get('bfrac25') or nan for g in gs])
    P.scatter(N.log10(ew),burst,s=(merg+1)*40,c=merg,cmap=cm.bone_r,alpha=0.8,lw=1.4,marker='o')

    P.ylabel('$\mathrm{burst\, fraction}$')
    P.xlabel(r'$\log \mathrm{EW}(H\alpha)$')

def L_vsigPlot():
    gs=M.Galax.objects.filter(sample=1)
    vsig=N.array([g.maxvel/2/g.sigma_cent/sq3 for g in gs])
    merg=N.array([g.doub+g.irrvf for g in gs])
    LHa=N.array([g.p.get('HaLum') for g in gs],dtype='Float64')
    P.scatter(vsig,N.log10(LHa),s=(merg+1)*40,c=merg,cmap=cm.bone_r,alpha=0.8,lw=1.4,marker='D')

    gs=M.Galax.objects.filter(sample=2)
    vsig=N.array([g.maxvel/2/g.sigma_cent/sq3 for g in gs])
    merg=N.array([g.doub+g.irrvf for g in gs])
    LHa=N.array([g.p.get('HaLum') for g in gs],dtype='Float64')
    P.scatter(vsig,N.log10(LHa),s=(merg+1)*40,c=merg,cmap=cm.bone_r,alpha=0.8,lw=1.4,marker='o')

    P.xlabel(r'$v_{max} / \sigma^\ast_c$')
    P.ylabel(r'$\log L(H\alpha)\,\,(erg\,s^{-1})$')

def EW_vsigPlot():
    gs=M.Galax.objects.filter(sample=1)
    vsig=N.array([g.maxvel/2/g.sigma_cent/sq3 for g in gs])
    merg=N.array([g.doub+g.irrvf for g in gs])
    ew=N.array([g.p.get('HaEW') for g in gs],dtype='Float64')
    P.scatter(vsig,N.log10(ew),s=(merg+1)*40,c=merg,cmap=cm.bone_r,alpha=0.8,lw=1.4,marker='D')

    gs=M.Galax.objects.filter(sample=2)
    vsig=N.array([g.maxvel/2/g.sigma_cent/sq3 for g in gs])
    merg=N.array([g.doub+g.irrvf for g in gs])
    ew=N.array([g.p.get('HaEW') for g in gs],dtype='Float64')
    P.scatter(vsig,N.log10(ew),s=(merg+1)*40,c=merg,cmap=cm.bone_r,alpha=0.8,lw=1.4,marker='o')

    P.xlabel(r'$v_{max} / \sigma_c^\ast$')
    P.ylabel(r'$\log \mathrm{EW}(H\alpha)$')

def classifPlot():
    P.axis((-13.1,-21,-0.2,4.4))
    P.yticks(range(5),['$%d$'%i for i in range(5)])
    gs=M.Galax.objects.filter(sample=1,mb__isnull=False)
    for g in gs:
        P.text(g.mb,g.doub+g.irrvf+0.1,'$\mathrm{%s}$'%g.p.get('morphclass'),
            bbox=dict(fc='black', alpha=0.1, ec='white'))

    gs=M.Galax.objects.filter(sample=2)
    for g in gs:
        P.text(g.mb,g.doub+g.irrvf-0.1,'$\mathrm{%s}$'%g.p.get('morphclass'))

    P.xlabel(r'$M_B$')
    P.ylabel(r'$P_v\, +\, P_d$')


def Plot():
    gs=M.Galax.objects.filter(sample=1)
    vsig=N.array([g.maxvel/2/g.sigma_cent/sq3 for g in gs])
    merg=N.array([g.doub+g.irrvf for g in gs])
    ew=N.array([g.p.get('HaEW') for g in gs],dtype='Float64')
    P.scatter(vsig,N.log10(ew),s=(merg+1)*40,c=merg,cmap=cm.bone_r,alpha=0.8,lw=1.4,marker='D')

    gs=M.Galax.objects.filter(sample=2)
    vsig=N.array([g.maxvel/2/g.sigma_cent/sq3 for g in gs])
    merg=N.array([g.doub+g.irrvf for g in gs])
    ew=N.array([g.p.get('HaEW') for g in gs],dtype='Float64')
    P.scatter(vsig,N.log10(ew),s=(merg+1)*40,c=merg,cmap=cm.bone_r,alpha=0.8,lw=1.4,marker='o')

    P.xlabel(r'$v_{max} / \sigma_c^\ast$')
    P.ylabel(r'$\log \mathrm{EW}(H\alpha)$')



##################################

def plotbadpix(color='w'):
    """ for the Argus FOV """

    P.plot([0,1],[0,1],color)
    P.plot([1,2],[0,1],color)
    P.plot([20,21],[0,1],color)
    P.plot([21,22],[0,1],color)
    P.plot([0,1],[13,14],color)
    P.plot([1,2],[13,14],color)
    P.plot([20,21],[13,14],color)
    P.plot([21,22],[13,14],color)
    P.plot([3,4],[4,5],color)
    P.plot([20,21],[8,9],color)
    P.plot([20,21],[9,10],color)
    P.plot([20,21],[10,11],color)

    P.plot([0,1],[1,0],color)
    P.plot([1,2],[1,0],color)
    P.plot([20,21],[1,0],color)
    P.plot([21,22],[1,0],color)
    P.plot([0,1],[14,13],color)
    P.plot([1,2],[14,13],color)
    P.plot([20,21],[14,13],color)
    P.plot([21,22],[14,13],color)
    P.plot([3,4],[5,4],color)
    P.plot([20,21],[9,8],color)
    P.plot([20,21],[10,9],color)
    P.plot([20,21],[11,10],color)


def plotonpix(data3,data2=None,layout='argus',line='k-',range=None):
    if not data2: data2=N.sum(data3,axis=2)
    if layout=='argus': x=22;y=14
    f=P.figure()
    P.imshow(N.transpose(data2))
    corn=P.gca().figbox.corners()
    sx,sy=corn[0,0],corn[0,1]
    lx,ly=corn[-1,0]-sx,corn[-1,1]-sy
    xl,yl=lx/x,ly/y
    xs,ys=sx,sy
    for i in N.arange(x):
        for j in N.arange(y):
            ax=P.axes([xs,ys,xl,yl],frameon=False)
            P.setp(ax, xticks=[], yticks=[])
            P.plot(data3[j,i,:],line)
            ys+=yl
        ys=sy
        xs+=xl

