"""
tool.py

general toolbox

"""

from PyGalKin import *

# SHORTCUTS
tab='\t'
nl='\n'
bs='\\'
sp=' '
null='\0'


# CONSTANTS
#
# units:
#        velocity: km/s
#        wavelength: Angstrom
lambHA=N.array([6562.7797852000003],'Float64')
sol=N.array([299792.458],'Float64')
c=sol
H0=N.array([70.],'Float64')
Grav=N.array([6.6726E-11*1.989E30/1000**3],'Float64') ### in solar masses and km
pc=N.array([3.086E13],'Float64') ## in km
sq3 = N.sqrt(3.0)
## Roberts values
## [S III] 9068.8
## Pa 10   9014.910
## Pa 11   8862.783
## Pa 12   8750.473
## Pa 13   8665.018
## Pa 14   8598.392
## Pa 15   8545.382
## Pa 16   8502.483
## Pa 17   8467.253
## O I     8446.455
## Pa 18   8437.955

# Paschen wavelengths
#                 Pa19         18  ....
Paschen=N.array([8413.317, 8437.955, 8467.253, 8502.483, 8545.382, 8598.392, 8665.018, 8750.473, 8862.783, 9014.910, 9229.014])
def PaLamb(number): return Paschen[19-number]


#PschenStrengths     9/10  10/10  11/10  12/10    ...            19/10
PaschStren=N.array([1.3812, 1.0, 0.7830, 0.6131, 0.4801, 0.3759, 0.2943, 0.2477, 0.2084, 0.1754, 0.1476])

PaschStren=PaschStren[::-1]

#                       SIII   OI  ClII  FeII
EmissionLines=N.array([9068.6,8446,8579,8617])
CaT=N.array([8498., 8542., 8662.])
Sulfur=9068.87



#
# BASE CLASSES
#
class numpdict(N.ndarray):
    def __new__(subtype, data, p=None, dtype=None, copy=False):
        # Make sure we are working with an array, and copy the data if requested
        subarr = N.array(data, dtype=dtype, copy=copy)

        # Transform 'subarr' from an ndarray to our new subclass.
        subarr = subarr.view(subtype)

        # Use the specified 'p' parameter if given
        if p is not None:
            subarr.p = p
        # Otherwise, use data's p attribute if it exists
        elif hasattr(data, 'p'):
            subarr.p = data.p

        # Finally, we must return the newly created object:
        return subarr

    def __array_finalize__(self,obj):
        if not hasattr(self, 'p'):
          self.p = getattr(obj, 'p', {})

    def __repr__(self):
        desc="""array(data=%(data)s,p=%(p)s)"""
        return desc % {'data': str(self), 'p':self.p }

class spec(numpdict):
    def __new__(subtype, data, p=None, dtype=None, copy=False):
        subarr=numpdict.__new__(subtype, data, p, dtype, copy)
        if p.has_key('l0'):
            subarr.l0=p['l0']
            if p.has_key('dl'):
                subarr.dl=p['dl']
                subarr.wave=N.arange(subarr.shape[-1])*subarr.dl + subarr.l0
        return subarr



# physical functions
def dis(arr1,arr2):
    """ returns the distance between two points"""
    #arr=(arr2-arr1)**2
    #return N.sqrt(arr.sum())
    return mlab.dist(arr1,arr2)

def app2abs(dist,m):
    """ convert apparent to absolute magnitudes, takes distance in Mpc"""
    return (-5*N.log10(dist*1000000.))+5+m

def balmer(m):
    """ caculate m'th balmer line"""
    return hydrogen(2,m+2)

def hydrogen(n,m):
    """ calculate rydberg wavelengths of hydrogen"""
    m,n=float(m),float(n)
    return (n**2 * m**2 / (m**2 - n**2)) / 1.0967758E7

def massKepler(r,v):
    """ returns the mass from keplers law: M(<r)=r*(v**2)*G
    input in pc and km/s
    mind the factors 1/2 !!
    """
    return ((v/2.)**2)*(r/2.)*pc/Grav

def dynMassDisk(r,sigma):
    'r in kpc, sigma in km/s, returns solar masses'
    return 7.9E5 * r * sigma**2

def dynMassSphere(r,sigma):
    'r in kpc, sigma in km/s, returns solar masses'
    return 1.1E6 * r * sigma**2

def lamb2vel(l,rest=lambHA):
    """ converts a wavelength in A wrt HA into a radial velocity """
    return ((l/rest)-1)*sol

def flux2mag(f):
    return -2.5*(N.log10(f))

def mag2flux(m):
    return 10**(m/2.5)

def z2vel(z):
    return z*sol

def vel2z(v):
    return v/sol

def lamb2freq(l):
    return 1.0E3*sol/l

def freq2lamb(f):
    return 1.0E3*sol/f

def Ghz2micron(f):
    return freq2lamb(f)*1.0E-3

def micron2Ghz(l):
    return lamb2freq(l)*1.0E-3

def arcsec2rad(arcsec):
    return N.radians(arcsec/3600.0)

def arcsec2kpc(arcsec=1.0,vsys=1000):
    return vsys/H0*1E3*arcsec2rad(arcsec)

def kpc2arcsec(kpc=1.0,vsys=1000):
    return kpc/arcsec2kpc(vsys=vsys)

def hubbledist(v):
    return v/H0*1000

def units(have,want,number=''):
    """
    uses the external procram "units" to convert units :-)
    """
    out=commands.getoutput('units -q -s ' + str(number) + have + ' ' + want + '| head -1 | cut -d " " -f2')
    return N.array([float(out)])

def lamb2pix(data,Lamb0,Step):
    if type(data) == type(1) or type(data) == type(1.0):
        return int(N.around((data-Lamb0)/Step).astype('Int32'))
    else: return N.around((data-Lamb0)/Step).astype('Int32')

def dlamb2vel(data,lamb0):
    return data/lamb0*c

def pix2lamb(data,Lamb0,Step):
    return (data*Step)+Lamb0

def pix2vel(data,lamb0,Step):
    return z2vel(((data*Step)+lamb0 )/ lamb0)

def pix2relvel(data,lamb0,Step):
    return data*Step/lamb0*c

def vel2lamb(data,lamb0):
    return vel2z(data) * lamb0

def vel2z(vel):
    return ((vel/c)+1)

def z2vel(z):
    return (z-1)*c

def vel2dis(v):
    return v/H0

def scalefromvarc(arcsecperpix,v):
    """ returns pc/pix"""
    return arcsecperpix/3600/360*2*pi * vel2dis(v)*1E6

def isconstant(data):
    return S.std(data)==0.0

def relz(v):
    return sqrt((1+(v/c))/(1-(v/c)))-1


# formatting
def deg2hms(deg):
    x = deg / 360.0 * 24.0
    h = int( x //1 )
    x = (x%1) * 60
    m = int( x //1 )
    s = (x%1) * 60
    return h,m,Decimal('%.1f'%round(s,1))

def deg2dms(deg):
    if deg < 0: f = -1
    else: f = 1
    deg *= f
    d = int(deg //1)
    deg = (deg%1)*60
    m = int(deg//1)
    s = int(round((deg%1)*60))
    return f*d,m,s


# Handy general functions
def getXY(data):
    i=N.indices((data.shape[0],data.shape[1]))
    return N.ravel(i[0]),N.ravel(i[1])

def shift(vec,i):
    """ Shift a vector.
    Usage: new_vec = shift(vec, i)
    vec:  The vector to be shifted
    i:  The steps to shift the vector with
    """
    n= vec.size
    i %= n
    return N.concatenate((vec[n-i:n],vec[0:n-i]))

def calcpeak(inarr,n):
    """ calculate the barycenter-velociy of the n highest pixels"""
    sorted,args=N.sort(inarr),N.argsort(inarr)
    erg = N.sum(sorted[-n:] * args[-n:]) / N.sum(sorted[-n:])
    return erg

def firstmoment(inarr):
    """ calculates the 1st moment from an array """
    return N.sum(N.arange(inarr.size)*inarr)/N.sum(inarr)

def secondmoment(inarr):
    """ calculates the 2nd moment from an array xx"""
    return N.sqrt(N.abs(N.sum((N.arange(inarr.size)-firstmoment(inarr))**2*inarr)/N.sum(inarr)))

def fwhm(inarr):
    """ returns the full widh at half maximum"""
    return 2.3548200450309493*secondmoment(inarr)


def doforeachpoint(data, function, *args, **keywords):
    """Apply a function (whcih takes a 1d-vector) to all values of x and
    y of a 3D-matrix. The output will have a z-dimension equal to the
    length of the output from the 'function'.

    Usage: new_arr = doforeachpoint(arr, function, arguments)
    data: The 3D-array input array

    """
    data=data.copy()
    x,y,z=data.shape
    xy=x*y
    data.shape=(xy,z)

    erg=None
    for i in N.arange(xy):
        tmp=function(data[i,:], *args, **keywords)
        #print i,tmp
        if not hasattr(tmp,'__len__'): tmp=N.array([tmp]);
        if erg==None:
            try:
                erg=N.zeros((xy,len(tmp)),dtype='float64')
                #print 'first value'
            except: continue

        erg[i,:]=tmp

    erg.shape=(x,y,-1)
    if erg.shape[2]==1: erg.shape=(x,y)
    return erg


def selective_sum(data,range='cat',Z=1.002912,axis=2):
    if range=='cat': zmin,zmax=lamb2pix(N.array([8470,8700])*Z,Lamb0,Step)
    else: zmin,zmax=0,data.shape[-1]
    print data.shape
    return N.sum(data[:,:,zmin:zmax],axis)
selsum=selective_sum

def selective_average(data,range='cat',Z=1.002912,axis=2):
    if range=='cat': zmin,zmax=lamb2pix(N.array([8470,8700])*Z,Lamb0,Step)
    else: zmin,zmax=0,data.shape[-1]

    if len(data.shape)==3: return N.average(data[:,:,zmin:zmax],axis)
    elif len(data.shape)==1: return N.average(data[zmin:zmax])
selav=selective_average

def pa2vec(pa):
    return N.array([-N.sin(radians(pa)),N.cos(radians(pa))])

def posvel(vf,dyncen,pa):
    vec=pa2vec(pa)
    x,y=N.indices(vf.shape)
    pos=N.inner(vec,N.transpose([dyncen[0]-x,dyncen[1]-y])).reshape(-1)
    vel=N.ma.array(vf,mask=vf.mask).reshape(-1)
    return pos, vel

def m2masks(angmap,pa,wedge):
    pa2=pa+pi
    mask1=angmap < pa-wedge
    if pa-wedge<0:
        mask1=mask_or(mask1, (angmap < 2*pi+pa-wedge) & (angmap>pi))
        mask1=mask_or(mask1, (angmap > pa+wedge) & (angmap<pi))
    else:
        mask1=mask_or(mask1, angmap > pa+wedge)
    mask2=angmap > pa2+wedge
    if pa2+wedge>2*pi:
        mask2=mask_or(mask2, (angmap > pa2-2*pi+wedge) & (angmap<pi))
        mask2=mask_or(mask2, (angmap < pa2-wedge) & (angmap>pi))
    else:
        mask2=mask_or(mask2, angmap < pa2-wedge)
    return mask1,mask2

def Angmap(x,y,pa):
    angmap=N.arctan(x/y)*(-1)
    angmap=N.where(y<0.0,angmap+pi,angmap)
    angmap=N.where(angmap<0.0,angmap+(2*pi),angmap)
    angmap-=pa
    angmap=N.where(angmap<0.0,angmap+(2*pi),angmap)
    angmap=N.where(N.isnan(angmap),pi,angmap)
    return angmap

def Dismap(x,y,pa,incl):
    vec=N.array([-N.sin(pa),N.cos(pa)])
    perp=N.array([N.cos(pa),N.sin(pa)])
    perp/=N.cos(incl)
    d1=(x*vec[0])+(y*vec[1])
    d2=(x*perp[0])+(y*perp[1])
    return N.sqrt(d1**2 + d2**2)


def binRC(rin,vin,rbin=1.0):
    n=N.ceil(rin.max()/rbin)
    R=(N.arange(n))*rbin
    V=N.zeros_like(R)
    S=N.zeros_like(R)
    for i,r in enumerate(R):
        vt=masked_where((rin<r)|(rin>r+rbin),vin)
        V[i]=vt.mean()
        S[i]=vt.std()
    return R+(rbin/2.0),V,S

def rotcur(vf,cen,pa,wedge,incl):
    """ calculate a rotation curve from a VF"""
    while pa < 0.0: pa+=180.0
    while pa >= 180.0: pa-=180.0
    pa,wedge,incl=map(N.radians,(pa,wedge,incl))
    x,y=getXY(vf)
    x.shape=vf.shape
    y.shape=vf.shape
    x=x.astype('f') - cen[0]
    y=y.astype('f') - cen[1]

    angmap=Angmap(x,y,pa)
    dismap=Dismap(x,y,pa,incl)

    vf=vf.copy()
    offset=vf[cen[0],cen[1]]
    vf-=offset
    vf=vf/N.abs(N.cos(angmap))/N.sin(incl)

    #mask1,mask2=m2masks(angmap,pa,wedge)
    mask1=(angmap>wedge) & (angmap<2*pi-wedge)
    mask2=mask_or(angmap<pi-wedge,angmap>pi+wedge,copy=True)
    r1=masked_array(dismap,mask1).flatten()
    r2=masked_array(dismap,mask2).flatten()
    v1=masked_array(vf,mask1).flatten()
    v2=masked_array(vf,mask2).flatten()

    #return vf,masked_array(N.cos(angmap),mask1&mask2)
    return r1,r2,v1+offset,v2+offset


#########################
####  Cross correlation
#########################

def xcorr(galaxy,star,filtgal=False,filtstar=None,range=N.array([700,1300]),baryN=15,plot=False,offset=50):

    Lamb0,Step,contSubtr=A.Lamb0,A.Step,A.contSubtr
    wavecal=pix2lamb(range,Lamb0,Step)

    if filtstar != None: gaussian_filter1d(star,filtstar)
    origshape=galaxy.shape
    galaxy=galaxy.copy()
    star=star.copy()
    star-=contSubtr(star,order=1)
    star=star[range[0]:range[1]]
    star,x=log_rebin(star,wavecal)

    if len(galaxy.shape) == 3:
        galaxy.shape=(origshape[0]*origshape[1],origshape[2])

    pos=N.zeros(galaxy.shape[0],'Float32')
    wid=N.zeros(galaxy.shape[0],'Float32')
    bary=N.zeros(galaxy.shape[0],'Float32')
    secmom=N.zeros(galaxy.shape[0],'Float32')
    cont=N.zeros(galaxy.shape[0],'Float32')
    amp=N.zeros(galaxy.shape[0],'Float32')
    h3=N.zeros(galaxy.shape[0],'Float32')
    h4=N.zeros(galaxy.shape[0],'Float32')

    fitl=(len(star)/2)-80+offset
    fitr=(len(star)/2)+80+offset
    print fitl,fitr
    for i in N.arange(len(pos)):
        if isconstant(galaxy[i]): continue
        gal=galaxy[i,range[0]:range[1]]
        gal-=contSubtr(gal,order=1)
        gal,x1=log_rebin(gal,wavecal)
        print x1[1]-x1[0]
        if filtgal:
            gal.ndim=1
            gal=bandfilt(gal)

        xc=Sig.correlate(gal,star,'same')
        if plot: P.clf(); P.plot(xc); sleep(3);P.clf()
        if plot: sleep(4);P.clf()
        xc=xc[fitl:fitr]
        fit=G.fitgaussh34(xc+1.0,err=1/xc,plot=plot,prin=True)
        #fit=fit2gauss(xc+1.0,plot=True)
        if fit == -1: print "somethings wrong!"
        else:
            cont[i],pos[i],amp[i],wid[i],h3[i],h4[i]=fit.params
            #cont,pos[i],amp,wid[i],x4,x5,x6=fit.params

        bary[i]=calcpeak(xc,baryN)
        xc=N.array(xc)
        secmom[i]=secondmoment(xc)


        if plot: P.plot([pos[i],bary[i]],[amp[i]+cont[i],amp[i]+cont[i]],'ro')
        P.draw()

    pos=pos-((fitr-fitl)/2.0)+offset
    bary=bary-((fitr-fitl)/2.0)+offset
    #vel=((exp(x[0])/exp(x[diff]))-1)*3E5
    if len(origshape) == 3:
        pos.shape=(origshape[0],origshape[1])
        wid.shape=(origshape[0],origshape[1])
        bary.shape=(origshape[0],origshape[1])
        secmom.shape=(origshape[0],origshape[1])
        cont.shape=(origshape[0],origshape[1])
        amp.shape=(origshape[0],origshape[1])
        h3.shape=(origshape[0],origshape[1])
        h4.shape=(origshape[0],origshape[1])

    return pos,bary,wid,secmom,cont,amp,h3,h4


def offset2vel(data,calib=4.93750657338e-05):
    """I *think* this is simply to apply a previously determined calibration"""
    return (N.exp(data*calib)-1) * c


def sigmacal(star,plot=False):
    sigmain=N.arange(0,20,1.0,'Float32')
    sigmaout=sigmain.copy()*0.0
    wavecal=[8000,8000+(Step*len(sigmain))]
    star,x=log_rebin(star,wavecal)
    for i in N.arange(len(sigmain)):
        st,x1=log_rebin(gaussian_filter1d(star,sigmain[i]),wavecal)
        xc=Sig.correlate(star,st,'same')
        xc=xc[280:-280]
        if plot: P.clf()
        fit=fitgauss(xc+1.0,plot=plot)
        print fit.params[3]
        sigmaout[i]=fit.params[3]
    return dlamb2vel(sigmain*Step,CaT[1]),sigmaout


def applysigcal(data,cal):
    mi=cal.x.min()
    ma=cal.x.max()

    dat=N.where(data <= ma, data, data*0.0 + mi)
    dat=N.where(dat >= mi, dat, data*0.0 +mi)
    cadat=cal(dat.flat)
    cadat.shape=dat.shape
    return cadat


#########################
####  IDL Wrappers
#########################

def log_rebin(spec,lamRange=None):
    """ wrapper for IDL's log_rebin"""

    # make a new IDL session
    idl=IDL()

    # give the variables to IDL
    idl.put('spec',spec)
    idl.put('lamRange',lamRange)

    #construct the IDL command and execute it
    idlcommand='LOG_REBIN, lamRange, spec, specNew, logLam, VELSCALE=velScale'
    idl.eval(idlcommand)

    # get the result
    specNew=N.array(idl.get('specNew'))
    logLam=N.array(idl.get('logLam'))

    return specNew,logLam


def ppxf():
    """ wrapper for ppxf in IDL"""
    #PPXF, star, galaxy, noise, velScale, start, sol, $
    #;       BESTFIT=bestFit, BIAS=bias, /CLEAN, DEGREE=degree, ERROR=error, $
    #;       GOODPIXELS=goodPixels, MDEGREE=mdegree, MOMENTS=moments, $
    #;       /OVERSAMPLE, /PLOT, /QUIET, VSYST=vsyst, WEIGHTS=weights



def voronoi2dbinning(data,Noise=False,targetSN=20,plot=True,quiet=True,returnall=False):
    """ wrapper to do voronoi binning
     CAREFUL: treats 2d-data as two spatial dimensions
    """
    origshape=data.shape
    if len(data.shape) == 3 and Noise==False:
        X,Y=getXY(data)
        data.shape=(origshape[0]*origshape[1],origshape[2])
        Signal=data.mean(axis=1)
        Noise=data.std(axis=1)
        data.shape=origshape
    elif len(data.shape) == 3 and Noise!=False:
        X,Y=getXY(data)
        data.shape=(origshape[0]*origshape[1],origshape[2])
        Signal=data.mean(axis=1)
        Noise=N.resize(Noise,Signal.shape)
        data.shape=origshape
    elif len(data.shape) == 2:
        Signal=N.ravel(data)
        if len(Noise) != len(Signal): Noise=N.resize(Noise,Signal.shape)
        X,Y=getXY(data)
    elif len(data.shape) == 1 and Noise!=False:
        Signal=data
        if len(Noise) != len(Signal): Noise=N.resize(Noise,Signal.shape)
        X,Y=getXY(data)
    else:
        print "must have a noise level for non-spectral data"
        return -1

    #print Signal.shape,Noise.shape,X.shape,Y.shape
    print max(Signal),S.average(Noise),X[N.argmax(Signal)],Y[N.argmax(Signal)]
    # make a new IDL session
    idl=IDL()

    # give the variables to IDL
    idl.put('X',X)
    idl.put('Y',Y)
    idl.put('Signal',Signal)
    idl.put('Noise',Noise)
    idl.put('targetSN',targetSN)

    #construct the IDL command
    idlcommand='VORONOI_2D_BINNING, X, Y, Signal, Noise, targetSN, BinNumber, xBin, yBin, xBar, yBar, SN, nPixels'
    if plot: idlcommand+=', /PLOT'
    if quiet: idlcommand+=', /QUIET'

    # run the command and save the plot

    try:
        idl.eval('set_plot,\'ps\'')
        idl.eval(idlcommand)
        idl.eval('device,/close')
    except:
        print "something went wron while running idl"
        return -1

    # collect the output
    BinNumber=N.array(idl.get('BinNumber'))
    xBin=N.array(idl.get('xBin'))
    yBin=N.array(idl.get('yBin'))
    xBar=N.array(idl.get('xBar'))
    yBar=N.array(idl.get('yBar'))
    SN=N.array(idl.get('SN'))
    nPixels=N.array(idl.get('nPixels'))

    if returnall: return BinNumber, xBin, yBin, xBar, yBar, SN, nPixels
    else: return BinNumber

def avbins2(data,BinNumber):
    """ average the data accordng to BinNumber, but return a 1d-vector instead of the same shape as data"""
    n=BinNumber.max()+1
    data=data.flatten()
    result=N.zeros(n,dtype='Float32')
    num=N.zeros(n,dtype='Int32')
    sig=N.zeros(n,dtype='Float32')
    for i in N.arange(n):
        result[i]=N.sum(N.where(BinNumber==i,data,0.0))
        num[i]=N.sum(N.where(BinNumber==i,1,0))
        sig[i]=masked_where(BinNumber!=i,data).std()
    return result/num,sig,num

def avbins3(data,BinNumber):
    """ average the data accordng to BinNumber, but return a 1d-vector instead of the same shape as data"""
    n=BinNumber.max()+1
    os=data.shape
    data=data.copy()
    data.shape=(os[0]*os[1],os[2])
    result=N.zeros((n,os[2]),dtype='Float32')
    num=N.zeros(n,dtype='Int32')
    for i in N.arange(n):
        num[i]=N.sum(N.where(BinNumber==i,1,0))
        result[i,:]=N.sum(data[N.where(BinNumber==i),:],axis=1) / num[i]
    return result,num

def spreadbins2(data,BinNumber,shape=None):
    if not shape: shape=(N.sqrt(BinNumber.shape).astype('i'),N.sqrt(BinNumber.shape).astype('i'))
    result=N.zeros(shape,dtype='Float32').flatten()
    for i,bin in enumerate(BinNumber):
        result[i]=data[bin]
    result.shape=shape
    return result



def average_bins3(data,BinNumber):
    """BinNumber is of length Npix and contains for each pix the bin-number that it belongs to"""
    orig=data.copy()
    Nbins=max(BinNumber)+1
    data=N.reshape(data,(orig.shape[0]*orig.shape[1],orig.shape[2]))

    BinValues=N.zeros((Nbins,orig.shape[2]),'Float32')
    counter=N.zeros((Nbins,))

    for i in N.arange(len(BinNumber)):
        BinValues[BinNumber[i],:] += data[i,:]
        counter[BinNumber[i]] += 1

    for i in N.arange(len(BinNumber)):
        data[i,:]=BinValues[BinNumber[i]] / counter[BinNumber[i]]

    data.shape=orig.shape
    return data

def average_bins2(data,BinNumber,prin=False):
    """BinNumber is of length Npix and contains for each pix the bin-number that it belongs to"""
    origshape=data.shape

    data=N.ravel(data.copy())
    sig=N.zeros_like(data)
    num=N.zeros_like(data)

    BinValues,BinSigmas,BinNum=binvalues(data,BinNumber)

    for i in N.arange(len(BinNumber)):
        data[i]=BinValues[BinNumber[i]]
        sig[i]=BinSigmas[BinNumber[i]]
        num[i]=BinNum[BinNumber[i]]

    data.shape=origshape
    sig.shape=origshape
    num.shape=origshape
    return data,sig,num


def binvalues(data,BinNumber):
    Nbins=max(BinNumber)+1
    BinValues=[N.array([])]*Nbins
    #print Nbins, BinValues.shape,data.shape
    for i,bin in enumerate(BinNumber):
        BinValues[bin]= N.hstack((BinValues[bin],data[i]))

    return map(N.mean,BinValues),map(N.std,BinValues),map(N.size,BinValues)

def rad_profile(data,xbin,ybin,xcen,ycen,BinNumber):
    BinValues=binvalues(data,BinNumber)
    diff=N.sqrt(((xbin-xcen)**2)+((ybin-ycen)**2))
    P.plot(diff,BinValues,'x')


def bandfilt(data):

    lpcf = 0.2
    lpsf = 0.25
    hpcf = 0.7
    hpsf = 0.6

    Rp = 2
    Rs = 20
    #print [lpcf,hpcf],[lpsf,hpsf],Rp,Rs
    #return lhpfilt(data,params=[[lpcf,hpcf],[lpsf,hpsf],Rp,Rs])
    return lhpfilt(data,params=[lpcf,lpsf,Rp,Rs])

def lhpfilt(data,params=[0.006,0.01,0,20]):
    """
    wp, ws -- Passband and stopband edge frequencies, normalized from 0
    to 1 (1 corresponds to pi radians / sample).  For example:
                   Lowpass:   wp = 0.2,          ws = 0.3
                   Highpass:  wp = 0.3,          ws = 0.2
                   Bandpass:  wp = [0.2, 0.5],   ws = [0.1, 0.6]
                   Bandstop:  wp = [0.1, 0.6],   ws = [0.2, 0.5]
      gpass -- The maximum loss in the passband (dB).
      gstop -- The minimum attenuation in the stopband (dB).
      """
    data.ndim=1
    wp,ws,gpass,gstop=params
    [n,Wn] = Sig.buttord(wp,ws,gpass,gstop)
    [b,a] = Sig.butter(n,Wn)
    return filtfilt(b,a,data)

def hpfilt(data):
    pass

#########################
####  HELPER FUNCTIONS
#########################

def new_wcs(pix,coord,delta,naxis=2,equinox=2000.0):
    wcs = pywcs.WCS(naxis=naxis)
    wcs.wcs.crpix = pix
    wcs.wcs.cdelt = N.array(delta) / 3600.0
    wcs.wcs.crval = coord
    wcs.wcs.ctype = ["RA---AIR", "DEC--AIR"]
    wcs.wcs.set_pv([(2, 1, 45.0)])
    wcs.wcs.equinox = equinox
    return wcs

def smooth_gauss(data,sigma):
    gauss=Sig.gaussian(10*sigma,sigma)
    return Sig.convolve(data,gauss/N.sum(gauss),mode='same')

def smooth_box(data,npix):
    return Sig.convolve(data,N.ones(npix)/npix,mode='same')

def fourier_CC(data,templ):
    return Sig.correlate(fft(data),fft(templ),mode='same')

def combinecubes(cubes,method='median'):
    origshape=cubes[0].shape
    bigcube=N.array([])
    for cube in cubes:
        bigcube=N.concatenate((bigcube,N.ravel(cube)))
    bigcube.shape=(len(cubes),N.product(origshape))
    return N.reshape(S.median(bigcube,axis=0),origshape)

def medianspec(data):
    """  """
    if len(data.shape) == 2:
        medi=N.median(data)
    elif len(data.shape) == 3:
        medi=N.reshape(data,(data.shape[0]*data.shape[1],data.shape[2]))
        medi=N.median(medi)
    else: medi=data

    return medi



def degrade_old(data,factor=4.25):
    oldlen=data.shape[-1]
    newlen=int(N.floor(oldlen/factor))
    degr=N.zeros(newlen,'Float32')
    for i in N.arange(newlen):
        lower=int(N.ceil(i*factor))
        upper=int(N.floor((i+1)*factor))-1
        if i%2==0: split=upper+1
        else: split=lower-1
        degr[i]=N.sum(data[lower:upper+1])+ (data[split]/2.0)

    return degr/factor

def degrade(data,factor=4.25,quadratic=False):
    extfactor=1
    while (factor*extfactor)%1 != 0:
        extfactor+=1
    #print extfactor
    oldlen=data.shape[-1]
    newlen=int(N.floor(oldlen/factor))
    ldata=N.resize(data,(extfactor,oldlen))
    ldata=N.transpose(ldata).flat
    degr=N.zeros(newlen,'Float32')
    fac=int(factor*extfactor)
    for i in N.arange(newlen):
        #print len(ldata[i*fac:(i+1)*fac])
        if quadratic:
            degr[i]=N.sqrt(N.sum((ldata[i*fac:(i+1)*fac])**2))/N.sqrt(fac)
        else:
            degr[i]=N.sum(ldata[i*fac:(i+1)*fac])/fac

    return degr

def degradeall(data,factor=4.25,quadratic=False):
    origshape=data.shape
    if len(data.shape) == 3:
        data.shape=(origshape[0]*origshape[1],origshape[2])

    npix=data.shape[0]
    newlen=int(N.floor(data.shape[-1]/factor))
    degrad=N.zeros((npix,newlen),'Float32')
    for i in N.arange(npix):
        degrad[i]=degrade(data[i,:],factor,quadratic=quadratic)

    #print origshape,data.shape
    data.shape=origshape
    if len(data.shape) == 3: degrad.shape=(origshape[0],origshape[1],newlen)
    return degrad

def sortbins(data,error,wave,start,binwidth=0.85,end=False,log=False):
    origshape=data.shape
    if len(data.shape) == 3:
        data.shape=(origshape[0]*origshape[1],origshape[2])
        error.shape=(origshape[0]*origshape[1],origshape[2])
        wave.shape=(origshape[0]*origshape[1],origshape[2])
    if start < wave[:,0].max():
        print "setting start to"+str(wave[:,0].max())
        start=wave[:,0].max()
    if not end: end=wave[:,-1].min()
    if end > wave[:,-1].min():
        print "setting end to"+str(wave[:,-1].min())
        send=wave[:,-1].min()
    leng=int((end-start)/binwidth)
    end=start+(leng*binwidth)
    print start,end,binwidth,leng

    dat=N.zeros((data.shape[0],leng),'Float32')
    err=dat.copy()
    count=dat.copy()
    for i in N.arange(data.shape[0]):
        bins=((wave-start)/binwidth).astype('Int32')
        for j in N.arange(data.shape[1]):
            if (bins[i,j] >= 0) and (bins[i,j] <leng):
                #print i,j,bins.shape,bins[i,j]
                dat[i,bins[i,j]] += data[i,j]
                err[i,bins[i,j]] += error[i,j]
                count[i,bins[i,j]] += 1.0
        #print dat[i,:],count[i,:]
    dat /= count
    err /= count
    err /= N.sqrt(count)

    data.shape=origshape
    error.shape=origshape
    wave.shape=origshape
    return dat,err


def intdegrade(data,n,method=N.average):
    """ decrease resolution of an image by an integer number"""
    nx,ny=data.shape
    x,y=N.arange(nx/n),N.arange(ny/n)
    erg=N.zeros((nx/n,ny/n),dtype='f')
    for i in x:
        for j in y:
            #print i,j,data[i*n:(i+1)*n,j*n:(j+1)*n]
            erg[i,j]=method(data[i*n:(i+1)*n,j*n:(j+1)*n])
    return erg

def intdegradespec(data,n,method=N.average):
    """ decrease resolution of spectra by an integer number"""
    nx,ny,nz=data.shape
    x=N.arange(nx)
    y=N.arange(ny)
    z=N.arange(nz//n)
    erg=N.zeros((nx,ny,nz//n),dtype='f')
    print erg.shape
    for i in x:
        for j in y:
            for k in z:
                erg[i,j,k]=method(data[i,j,k*n:(k+1)*n])
                #print method(data[i,j,k*n:(k+1)*n]),erg[i,j,k]
    return erg


#####################################
# COOKBOOK and other functions from various sources
#####################################

def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = mgrid[-size:size+1, -sizey:sizey+1]
    g = exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()

def blur_image(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    improc = signal.convolve(im,g, mode='valid')
    return(improc)



def sortout(inarr,banned=0):
    """
    gives back three vectors: the values in inarr that are !=
    value and the corresponding x and y coordinates. this is old and
    ugly, use masked arrays instead!
    """
    nx,ny=inarr.nx(),inarr.ny()
    xaxis=N.array([])
    yaxis=N.array([])
    val=N.array([])
    for x in N.arange(nx):
        for y in N.arange(ny):
            if not inarr[x,y] == banned:
                #print x,y
                xaxis=N.concatenate((xaxis,x))
                yaxis=N.concatenate((yaxis,y))
                val=N.concatenate((val,inarr[x,y]))
    return xaxis,yaxis,val
