#!/usr/bin/env python
"""
PyArgus.py

Stuff to handle ARGUS IFU data

"""

from PyGalKin import *

#####################################
#### Wave-cal and dimensions
#####################################

#Lamb0=8182.43
#SpecLen=1407
#Step=0.85

#Lamb0=8183.425
#SpecLen=1357
#Step=0.85

Lamb0=8183.213
SpecLen=2715
Step=0.425

#Lamb0=8182.
#Step=0.19996649916247891
#SpecLen=5980

#SpecLenOrg=5980
SpecLenOrg=2715
#SpecLenOrg=1357


dimX=22
dimY=14

#skyregion=N.array([8750,8800])
#skyregion=N.array([8860,8910])
skyregion=N.array([8810,8876]) #used for He 2-10

#haro 11
#skyregion=N.array([8865,8925])



########################
## CONSTRUCTING THE CUBE
########################

def image2cube(data,tablefile='/home/tom/projekte/PyGalKin/argus-fibres.txt'):
    """allows both a filename and a 2d-array as input. the latter has to be flipped already """

    if type(data) == type(''):
        data=read_fits(data)
        data=data[:,::-1]
    elif type(data) == type(N.array([])):
        pass
    else:
        print('unknown type of input')
        return -1

    if data.shape[1]==311:
        havesimcal=False
        missing=4
    elif data.shape[1]==316:
        havesimcal=True
        missing=4
    elif data.shape[0]==317:
        havesimcal=True
        missing=3
        data=N.transpose(data)
    else:
        print('unknown type of input')
        return -1

    cube=N.zeros((dimX,dimY,SpecLenOrg),'Float32')
    sky=N.array([],'Float32')
    if havesimcal: simcal=N.array([],'Float32')
    
    file=open(tablefile,'r')

    # two header lines and missing spectra
    file.readline()
    file.readline()
    print(str(missing) + ' missing spectra')
    for i in N.arange(missing): file.readline()

    #print data.shape,cube.shape
    
    for line in file.readlines():
        line=line.split()
        index=int(line[1])-(missing+1)
        
        if 'Sky' in line[4]:
            sky=N.concatenate((sky,data[:,index]))
           
        elif 'Calibration' in line[4]:
            if havesimcal: simcal=N.concatenate((simcal,data[:,index]))
            else: missing+=1
        else:
            x,y=int(line[-3])-1,int(line[-2])-1
            #print x,y,index
            cube[x,y,:]=data[:,index]
        
    file.close()
    sky.shape=(sky.size/SpecLenOrg,SpecLenOrg)
    #badpixels(cube)
    if havesimcal:
        simcal.shape=(simcal.size/SpecLenOrg,SpecLenOrg)
        return cube,sky,simcal
    else: return cube,sky

def badpixels(data, value=0.0):
    """ sets the known bad spectra in a cube to value"""

    if len(data.shape)==3:
        data[0,0,:]=value
        data[1,0,:]=value
        data[20,0,:]=value
        data[21,0,:]=value
        data[0,13,:]=value
        data[1,13,:]=value
        data[20,13,:]=value
        data[21,13,:]=value
        data[3,4,:]=value
        data[20,8,:]=value
        data[20,9,:]=value
        data[20,10,:]=value
    elif len(data.shape)==2:
        data[0,0]=value
        data[1,0]=value
        data[20,0]=value
        data[21,0]=value
        data[0,13]=value
        data[1,13]=value
        data[20,13]=value
        data[21,13]=value
        data[3,4]=value
        data[20,8]=value
        data[20,9]=value
        data[20,10]=value


#####################
## SUBRTACTING STUFF
#####################
def skysub(data,sky,factor=1.9,region=skyregion):
    """ wants data in 2d or 3d, sky is first medianned to 1d, then grown"""
    shape=data.shape
    if len(data.shape) == 3:
        data.shape=(shape[0]*shape[1],shape[2])
    sky=medianspec(sky)
    #sky=N.resize(sky,data.shape)
    factor=skyfit(data,sky,region,quiet=False)
    dataSS=data.copy()
    for i in N.arange(data.shape[0]):
        dataSS[i]=data[i]-(factor[i]*sky)
    data.shape=shape
    dataSS.shape=shape
    return dataSS

def skyfit(data,sky,region=skyregion,quiet=True):
    factor=N.zeros(data.shape[0],'Float32')
    region=lamb2pix(region,Lamb0,Step)
    parinfo=[]
    for i in range(2):
        parinfo.append({'value':1.0, 'fixed':0, 'limited':[0,0],'limits':[0.0, 0.0], 'step':0.0})

    for i in N.arange(data.shape[0]):
        sdata=data[i,region[0]:region[1]]
        ssky=sky[region[0]:region[1]]
        #print sdata.shape,ssky.shape
        fa={'data':sdata,'sky':ssky}
        fit=mpfit(skyfunc,functkw=fa,parinfo=parinfo,maxiter=200,quiet=quiet)
        print(fit.status)
        factor[i]=fit.params[1]
    return factor
    

def skyfunc(p, fjac=None, data=None, sky=None, returnmodel=False):
    model= p[0] + (p[1]*sky)
    if returnmodel==True:
        return model
    else:
        status = 0
        return([status, (data-model)])

def contSubtr(data,order=6,sigmaclip=1.0,plot=False):
    if len(data.shape)==1: return contFit(data,order=order,sigmaclip=sigmaclip,plot=plot)
    origshape=data.shape
    if len(data.shape) == 3:
        data.shape=(origshape[0]*origshape[1],origshape[2])

    contSub=N.zeros(data.shape,'Float32')
    for i in N.arange(data.shape[0]):
        contSub[i,:]=data[i,:]-contFit(data[i,:],order=order,sigmaclip=sigmaclip,plot=plot)
        #print str(i)+' done'

    data.shape=origshape
    contSub.shape=origshape
    return contSub
    
def contFit(data,order=6,sigmaclip=1.0,plot=False):
    data=N.where(N.isnan(data),0.0,data)
    x=N.arange(len(data))
    try: poly=P.polyfit(x,data,order)
    except: print(data)
    #print poly
    subtr=data-P.polyval(poly,x)
    flagged=N.where(N.abs(subtr) > (sigmaclip*N.std(subtr)),0,subtr)
    corrpoly=P.polyfit(x,flagged,order)
    finalfit=P.polyval(poly,x)+P.polyval(corrpoly,x)
    if plot:
        P.plot(data)
        P.plot(flagged)
        P.plot(finalfit)
        P.plot(data-finalfit)
    return finalfit


    

##################################
### PASCHEN AND OTHER LINE FITTING
##################################

def PaModel(p, fjac=None, x=None, y=None, err=None, returnmodel=False):
    
    model=N.zeros(len(x),'Float32')+1.0
    print(p)
    PaNumbers=N.array([10,11,12,13,14,15,16,17])
    PaLa=PaLamb(PaNumbers)*p[0]
    for i in N.arange(len(PaLa)):
        para=[0.0,PaLa[i],p[-8+i],p[1]]
        #para=[0.0,lamb2pix(PaLa[i],Lamb0,Step),p[-8+i],p[1]]
        model+=G.gauss(para,x=x,returnmodel=True)

    if returnmodel==True:
        return model
    else:
        status = 0
        
        return([status, (y-model)/err])
    
def fitAllPaschen(data,guessV=None,plot=False,prin=False,quiet=True):
    PaNumbers=[9,10,11,12,14,17]
    parinfo=[]
    parinfo.append({'value':vel2z(guessV), 'fixed':0, 'limited':[0,0],'limits':[0.0, 0.0], 'step':0.0})
    parinfo.append({'value':5.0, 'fixed':0, 'limited':[0,0],'limits':[0.0, 0.0], 'step':0.0})
    for i in PaNumbers:
        parinfo.append({'value':0.003, 'fixed':0, 'limited':[0,0],'limits':[0.0, 0.0], 'step':0.0})

    x=N.arange(len(data))
    err=x*0.0 + 0.001
    
    fa = {'x':x, 'y':data, 'err':err}
    #print parinfo
    try:
        fit=mpfit(PaModel,functkw=fa,parinfo=parinfo,maxiter=200,quiet=quiet)
    except OverflowError:
        return -1

    print('fitAllPaschen status: ',fit.status)
    
    if plot==True:
        P.plot(data,'r')
        P.plot(PaModel(fit.params,x=x,returnmodel=True),'b')
    if prin==True:
        print(fit.niter,fit.params,fit.status)
    
    return fit.params

def fitAllPaschen_old(data,err,velRange=None,guessV=None,PaNumbers=[10,11,12,14,17],parinfo=None,plot=False,prin=False,quiet=True):
    relevant=N.array([],dtype='Float32')
    relerr=N.array([],dtype='Float32')
    once=False
    for p in PaNumbers:
        p=PaLamb(p)
        Left,Right= vel2lamb(guessV-(velRange/2.),p),vel2lamb(guessV+(velRange/2.),p)
        Left,Right=int(lamb2pix(Left,Lamb0,Step)),int(lamb2pix(Right,Lamb0,Step))
        if not once:
            pixels=Right-Left-1
            once=True
        #print Left,Right, pixels
        rel=data[Left:Left+pixels]
        #print rel
        rele=err[Left:Left+pixels]
        
        #rel-=min(rel)
        relevant=N.concatenate((relevant,rel))
        relerr=N.concatenate((relerr,rele))

    nlines=len(PaNumbers)
    if parinfo==None:
        parinfo=[]
        parinfo.append({'value':pixels*0.5, 'fixed':0, 'limited':[0,0],'limits':[0.0, float(pixels)], 'step':0.0})
        parinfo.append({'value':pixels*0.05, 'fixed':0, 'limited':[0,0],'limits':[0.0, pixels*0.5], 'step':0.0})
        for i in N.arange(nlines):
            #print i,pixels,relevant.size,relevant[i*pixels:(i+1)*pixels]
            parinfo.append({'value':min(relevant[i*pixels:(i+1)*pixels]), 'fixed':0, 'limited':[0,0],'limits':[min(relevant[i*pixels:(i+1)*pixels]), max(relevant[i*pixels:(i+1)*pixels])], 'step':0.0})
            parinfo.append({'value':max(relevant[i*pixels:(i+1)*pixels])-min(relevant), 'fixed':0, 'limited':[0,0],'limits':[0.0, max(relevant[i*pixels:(i+1)*pixels])*1.2], 'step':0.0})

    x=N.arange(len(relevant))
    
    fa = {'x':x, 'y':relevant, 'err':relerr, 'n':nlines}
    
    try:
        fit=mpfit(funcAllPaschen_old,functkw=fa,parinfo=parinfo,maxiter=200,quiet=quiet,gtol=1E-5)
    except OverflowError:
        return -1

    print('fitAllPaschen status: ',fit.status)
    
    if plot==True:
        P.plot(relevant,'r')
        P.plot(funcAllPaschen_old(fit.params,x=N.arange(len(relevant)),n=nlines,returnmodel=True),'b')
    if prin==True:
        print(fit.niter,fit.params,fit.status)
    
    return fit.params

def funcAllPaschen_old(p, fjac=None, x=None, y=None, err=None, n=None,returnmodel=False):
    model=N.zeros(len(x),'Float32')
    pixels=len(x)/n
    
    for i in N.arange(n):
        #print x[i*pixels:(i+1)*pixels]
        model[i*pixels:(i+1)*pixels] += p[(2*i)+3]*N.exp( -1* ((x[i*pixels:(i+1)*pixels]-(p[0]+(i*pixels)))**2) / (2*(p[1]**2))) + p[(2*i)+2]

    if returnmodel==True:
        return model
    else:
        status = 0
        return([status, (y-model)/err])


def findLine(data,type='single',velRange=None,guessV=None,restlamb=Sulfur,parinfo=None,plot=False,prin=False,quiet=True):
    
    Left= vel2lamb(guessV-(velRange/2.),restlamb)
    Right= vel2lamb(guessV+(velRange/2.),restlamb)
    Left,Right=int(lamb2pix(Left,Lamb0,Step)),int(lamb2pix(Right,Lamb0,Step))
    
    relevant=data[Left:Right]
    print(relevant,Left,Right,restlamb)
    if type=='single':
        fit=G.fitgauss(relevant,parinfo=parinfo,plot=plot,prin=prin,quiet=quiet)
    elif type=='double':
        fit=G.fit2gauss(relevant,parinfo=parinfo,plot=plot,prin=prin,quiet=quiet)
        
    elif type=='h34':
        fit=G.fitgaussh34(relevant,parinfo=parinfo,plot=plot,prin=prin,quiet=quiet)
    else:
        print('Unknown type of fit')
        return -1

    Z=pix2lamb(fit.params[1]+Left,Lamb0,Step) / restlamb
    return fit,Z
                            

def emissionVF(data,velRange=None,guessV=None,restlamb=Sulfur,type='single',plot=False,parinfo=None):
    origshape=data.shape
    if len(data.shape) == 3:
        data.shape=(origshape[0]*origshape[1],origshape[2])

    if type=='single':
        allparams=N.zeros((data.shape[0],4),'Float32')
    elif type=='double':
        allparams=N.zeros((data.shape[0],7),'Float32')
    elif type=='h34':
        allparams=N.zeros((data.shape[0],6),'Float32')
    else:
        print('Unknown type of fit')
        return -1


    for i in N.arange(data.shape[0]):
        
        results=findLine(data[i,:],restlamb=restlamb,velRange=velRange,guessV=guessV,type=type,plot=plot,parinfo=parinfo)
        allparams[i,:]=results
        
    data.shape=origshape
    allparams.shape=(origshape[0],origshape[1],-1)
    
    #print data.shape,EmVF.shape
    #P.matshow(EmVF)
    return allparams




def createPa(paschenparam,Z,type,PaNumb,D1=0.0,D2=0.0):
    Pasch=Paschen * Z

    # don't subtract continuum
    paschenparam[0]=0.0

    x=N.arange(SpecLen)
    SynthSpec=N.zeros(SpecLen,'Float32')

    Stren=PaschStren / PaschStren[19-PaNumb]
    #print Stren
    for i in N.arange(len(Paschen)):
        para=paschenparam.copy()
        para[2]*=Stren[i]
        para[1]=lamb2pix(Paschen[i]*Z,Lamb0,Step)+D1
        if type=='double':
            para[5]*=Stren[i]
            para[4]=lamb2pix(Paschen[i]*Z,Lamb0,Step)+D2
            SynthSpec+=G.twogauss(para,x=x,returnmodel=True)
        else:
            SynthSpec+=G.gauss(para,x=x,returnmodel=True)

    return SynthSpec


def createPaschen(data,type='double',velRange=None,guessV=None,plot=False,plotfit=False,PaNumb=10):
    fitresults=findLine(data,type=type,velRange=velRange,guessV=guessV,restlamb=PaLamb(PaNumb),plot=plotfit)

    if fitresults==-1:
        return N.zeros(SpecLen,'Float32')
    else:
        fit,Z=fitresults
        fit=fit.params
        print(fit,Z)

    SynthSpec=createPa(fit,Z,D1=0.0,D2=0.0,type=type,PaNumb=PaNumb)
    
    if plot:    
        P.plotspec(SynthSpec)
        P.plotspec(data)
        P.plotspec(data-SynthSpec,Z=Z,region='cat',plotlines=False)
    
    return SynthSpec

def createPaschenSul(data,velRange=None,guessV=None,plot=False,plotfit=False,PaNumb=10):

    
    fitresults=findLine(data,velRange=velRange,guessV=guessV,plot=plotfit)
    if fitresults==-1:
        return N.zeros(SpecLen,'Float32')
    else:
        fit,Z=fitresults
            
    Pasch=Paschen * Z

    parinfo=[]
    for i in range(7):
        parinfo.append({'value':0.0, 'fixed':0, 'limited':[0,0],'limits':[0.0, 0.0], 'step':0.0})
    parinfo[0]['value']=fitpara[0]
    parinfo[1]['value']=D1
    parinfo[1]['fixed']=1
    parinfo[2]['value']=(max(data)-min(data))/2
    parinfo[3]['value']=fitpara[3]
    parinfo[3]['fixed']=1
    parinfo[4]['value']=D2
    parinfo[4]['fixed']=1
    try:
        relampl=fitpara[5]/fitpara[2]
        parinfo[5]['tied'] = str(relampl)+'*p[2]'
    except ZeroDivisionError:
        parinfo[5]['fixed'] = 1
    

    parinfo[6]['value']=fitpara[6]
    parinfo[6]['fixed']=1
    
    print("second fit")
    fitresults=findLine(data,velRange=velRange,guessV=z2vel(Z),restlamb=PaLamb(PaNumb),parinfo=parinfo,plot=plotfit)
    if fitresults==-1:
        return N.zeros(SpecLen,'Float32')
    else:
        Z,paschenparam,D1,D2=fitresults

        
    SynthSpec=createPa(paschenparam,Z,D1=0.0,D2=0.0,double=True,PaNumb=PaNumb)
    
    if plot:    
        P.plotspec(SynthSpec)
        P.plotspec(data)
        P.plotspec(data-SynthSpec,Z=Z,region='cat',plotlines=True)
    
    return SynthSpec

def subtrPaschen(data,velRange=None,guessV=None,PaNumb=9,fromSul=True,double=True):
    origshape=data.shape
    if len(data.shape) == 3:
        data.shape=(origshape[0]*origshape[1],origshape[2])

    subtracted=N.zeros(data.shape,'Float32')
    if fromSul:
        for i in N.arange(data.shape[0]):
            subtracted[i,:]=data[i,:]-createPaschenSul(data[i,:],velRange=velRange,guessV=guessV,PaNumb=PaNumb)
    else:
        for i in N.arange(data.shape[0]):
            subtracted[i,:]=data[i,:]-createPaschen(data[i,:],velRange=velRange,guessV=guessV,PaNumb=PaNumb,double=double)
    data.shape=origshape
    subtracted.shape=origshape
    return subtracted


## In case this file gets executed...
if __name__ == '__main__':
    demo()
        

def demo():
    print("This file defines some functions. It is not meant to be executed. Import it instead!")

