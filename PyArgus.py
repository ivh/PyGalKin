#!/usr/bin/env python

#
# Some functions to handle ARGUS IFU data
#

from time import sleep
from os.path import exists
from tool import *
import Gauss as G
from plot import plotspec

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
#skyregion=N.array([8810,8876]) used for He 2-10

#haro 11
skyregion=N.array([8865,8925])



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
        print 'unknown type of input'
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
        print 'unknown type of input'
        return -1

    cube=N.zeros((dimX,dimY,SpecLenOrg),'Float32')
    sky=N.array([],'Float32')
    if havesimcal: simcal=N.array([],'Float32')
    
    file=open(tablefile,'r')

    # two header lines and missing spectra
    file.readline()
    file.readline()
    print str(missing) + ' missing spectra'
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
        print fit.status
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
        print str(i)+' done'

    data.shape=origshape
    contSub.shape=origshape
    return contSub
    
def contFit(data,order=6,sigmaclip=1.0,plot=False):

    x=N.arange(len(data))
    poly=P.polyfit(x,data,order)
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
    print p
    PaNumbers=N.array([10,11,12,14,17])
    PaLa=PaLamb(PaNumbers)*p[0]
    for i in N.arange(len(PaLa)):
        para=[0.0,lamb2pix(PaLa[i]),p[-6+i],p[1]]
        model+=gauss(para,x=x,returnmodel=True)

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

    print 'fitAllPaschen status: ',fit.status
    
    if plot==True:
        P.plot(data,'r')
        P.plot(PaModel(fit.params,x=x,returnmodel=True),'b')
    if prin==True:
        print fit.niter,fit.params,fit.status
    
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

    print 'fitAllPaschen status: ',fit.status
    
    if plot==True:
        P.plot(relevant,'r')
        P.plot(funcAllPaschen_old(fit.params,x=N.arange(len(relevant)),n=nlines,returnmodel=True),'b')
    if prin==True:
        print fit.niter,fit.params,fit.status
    
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


def findLine(data,double=True,velRange=None,guessV=None,restlamb=Sulfur,parinfo=None,plot=False,prin=False,quiet=True):
    
    Left= vel2lamb(guessV-(velRange/2.),restlamb)
    Right= vel2lamb(guessV+(velRange/2.),restlamb)
    Left,Right=int(lamb2pix(Left,Lamb0,Step)),int(lamb2pix(Right,Lamb0,Step))
    
    relevant=data[Left:Right]
    print relevant,Left,Right,restlamb
    if double:
        fit=G.fit2gauss(relevant,parinfo=parinfo,plot=plot,prin=prin,quiet=quiet)
        if fit==-1:
            print "fit went wrong!"
            return fit
        PeakPos=N.argmax(G.twogauss(fit.params,x=N.arange(len(relevant)),returnmodel=True))
    else:
        fit=G.fitgauss(relevant,parinfo=parinfo,plot=plot,prin=prin,quiet=quiet)
        if fit==-1: return fit
        PeakPos=N.argmax(G.gauss(fit.params,x=N.arange(len(relevant)),returnmodel=True))

    if fit.status != 1: print "findLine status:",fit.status
    #Z=pix2lamb(PeakPos+Left) / restlamb
    Z=pix2lamb(fit.params[1]+Left,Lamb0,Step) / restlamb
    #print fit.params,Z,Left
    D1=fit.params[1]-PeakPos
    if double: D2=fit.params[4]-PeakPos
    else: D2=PeakPos
    return Z,fit.params,D1,D2
                                

def emissionVF(data,velRange=None,guessV=None,restlamb=Sulfur,double=False,plot=False,parinfo=None):
    origshape=data.shape
    if len(data.shape) == 3:
        data.shape=(origshape[0]*origshape[1],origshape[2])

    EmVF=N.zeros(data.shape[0],'Float32')
    Cont=N.zeros(data.shape[0],'Float32')
    Ampl=N.zeros(data.shape[0],'Float32')
    Width=N.zeros(data.shape[0],'Float32')
    allparams=N.zeros((data.shape[0],7),'Float32')
    for i in N.arange(len(EmVF)):
        
        results=findLine(data[i,:],restlamb=restlamb,velRange=velRange,guessV=guessV,double=double,plot=plot,parinfo=parinfo)
        if results==-1:
            Z,params,D1,D2=0.0,N.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0]),0.0,0.0
        else:
            Z,params,D1,D2=results

        #print params.shape
        allparams[i,:]=params
        EmVF[i]=Z
        Cont[i]=params[0]
        if len(params)==4:
            Ampl[i]=params[2]
            Width[i]=params[3]
        else:
            Ampl[i]=params[2]+params[5]
        #print i,EmVF[i]
    
    data.shape=origshape
    EmVF=z2vel(EmVF)
    Width=pix2relvel(Width,restlamb,Step)
    #print EmVF.shape, origshape
    EmVF.shape=(origshape[0],origshape[1])
    Ampl.shape=(origshape[0],origshape[1])
    Cont.shape=(origshape[0],origshape[1])
    Width.shape=(origshape[0],origshape[1])
    allparams.shape=(origshape[0],origshape[1],-1)
    
    #print data.shape,EmVF.shape
    #P.matshow(EmVF)
    #return EmVF,Width,Ampl,Cont
    return allparams



def interpasch(data,error,velRange=None,guessV=None,PaNumb=9,prefix='intPa'):
    #sub=data.copy()
    #for i in range(data.shape[0]):
    #    for j in range(data.shape[1]):
    #        dat=data[i,j,:].copy()
    #        err=error[i,j,:].copy()
    #        if isconstant(dat):
    #            print 'skipping '+str(i)+' '+str(j)
    #            sub[i,j,:]=data[i,j,:]
    #            continue
    #        print 'Cuttent pixel: %s %s' % (i,j)
    inter=interactplot(data,error,prefix=prefix,velRange=velRange,guessV=guessV,PaNumb=PaNumb)
    P.show()
    sub[i,j,:]=inter.data - inter.shiftscaled()

    return sub
            

class interactplot:
    def __init__(self,data,error,velRange,guessV,prefix='intPa',PaNumb=9 ):
        self.odata=data
        self.oerror=error
        self.velRange=velRange
        self.guessV=guessV
        self.prefix=prefix
        self.i=-1
        self.j=0
        self.tiltfac=0.0
        self.flag=0
        self.x=data.shape[0]
        self.y=data.shape[1]
        #print self.x,self.y
        
        #self.double=double
        self.PaNumb=PaNumb
        self.osn=self.odata/self.oerror
        self.step=0.1
        self.fact=1.0
        self.shift=0

        if exists(prefix+'.pick'): self.subtracted=load(prefix+'.pick')
        else: self.subtracted=N.zeros(data.shape,'Float32')

        self.currdata()
        self.file=open(prefix+'.dat','a')
        self.fig=P.figure(1)#,figsize=(14,10))
        
        canvas = self.fig.canvas
        canvas.mpl_connect('key_press_event', self.key_press_callback)
        canvas.mpl_connect('button_press_event', self.button_press_callback)
        self.canvas = canvas

        #first one
        self.nextone()
        
    def nextone(self):
        if self.i+1 < self.x: self.i +=1
        elif self.j+1 < self.y:
            self.i=0
            self.j +=1
        else:
            print "DONE!"
            self.quit()
            
        self.currdata()
        if isconstant(self.data):
            print "skipping",self.i,self.j
            self.nextone()
        elif not isconstant(self.subtr):
            print "This one has already been subtracted."
            get=raw_input('Redo (y/N): ')
            if get != 'y': self.nextone()
            else: self.startover()
        else:
            self.startover()

    def update(self):
        self.currdata()
        self.measurePa()
        self.makesynt()
        self.plot()

    def currdata(self):
        self.data=self.odata[self.i,self.j,:]
        self.error=self.oerror[self.i,self.j,:]
        self.sn=self.osn[self.i,self.j]
        self.subtr=self.subtracted[self.i,self.j,:]

    def save(self):
        dump(self.subtracted,self.prefix+'.pick')
        self.file.write('%s %s %s %s %s %s %s\n'%(self.i,self.j,self.PaNumb,self.fact,self.shift,self.tiltfac, self.flag))
        
    def choosepix(self):
        self.i=int(raw_input('i: '))
        self.j=int(raw_input('j: '))
        self.update()

    def makesynt(self):
        if self.PaNumb == 8:
            self.osynt=createPaschenSul(self.data+1,velRange=self.velRange,guessV=self.guessV,plotfit=False)
        else:
            self.osynt=createPaschen(self.data+1,double=True,velRange=self.velRange,guessV=self.guessV,PaNumb=self.PaNumb,plotfit=False)
        self.synt=self.osynt.copy()
 
        
    def measurePa(self):
        params=fitAllPaschen_old(self.data,self.error,velRange=self.velRange,guessV=self.guessV,plot=False,prin=False)
        #print params
        self.paparams=params
        

    def shiftscaled(self):
        return shift(self.synt,self.shift)*self.fact * self.tilt()

    def tilt(self):
        return N.ones(len(self.data),'Float32') + (self.tiltfac *N.arange(-1.,1.,2.0/len(self.data)))
        
    def accept(self):
        self.subtracted[self.i,self.j,:]=self.data - self.shiftscaled()
        self.save()
        self.nextone()

    def quit(self):
        self.file.close()
        P.close(self.fig)

    def reject(self):
        self.file.write('%s %s %s\n'%(self.i,self.j,'R'))
        self.nextone()
    
    def startover(self):
        #self.fact=1.0
        #self.shift=0
        self.flag=0
        #self.tiltfac=0.0
        self.update()
        

    def chooseline(self,key):
        self.PaNumb=int(key)
        if self.PaNumb < 8: self.PaNumb +=10
        self.update()

    def smooth(self):
        pass

    def toggleflag(self):
        if self.flag == 0: self.flag = 1
        else: self.flag = 0
        self.accept()
        
    def key_press_callback(self,event):
        
        if event.key == '+': self.fact += self.step
        elif event.key == '-': self.fact -= self.step 
        elif event.key == 'l': self.shift -= 1
        elif event.key == 'r': self.shift += 1
        elif event.key == 's': self.smooth()
        elif event.key == 'o': self.startover()
        elif event.key == 'a': self.accept()
        elif event.key == 'x': self.reject()
        elif event.key == 'c': self.choosepix()
        elif event.key == 'q': self.quit()
        elif event.key == 'u': self.update()
        elif event.key == 'b': self.toggleflag()
        elif event.key == 'm': self.tiltfac += 0.1
        elif event.key == 'n': self.tiltfac -= 0.1
        
        
        elif event.key in '0123456789': self.chooseline(event.key)
        else: print "Unknown key pressed, doing nothing"
        self.plot()
        
    def plot(self):
        self.fig.clf()
        #print 'currently at %s %s'%(self.i,self.j)
        # Around CaT
        ax=P.axes([0.02,0.68,0.70,0.27])
        P.setp(ax,xticks=[], yticks=[])
        plotspec(self.shiftscaled(),style='-b')
        plotspec(self.data,style='-k')
        plotspec(self.data-self.shiftscaled(),style='-r',vminmax='sigbased',Z=vel2z(self.guessV),plotlines=True)
        
        P.title('CaT and Pa 13, 14, 15, 16')
        
        # SIII
        ax=P.axes([0.74,0.68,0.23,0.27])
        self.plotaroundline(Sulfur)
        P.setp(ax,xticks=[], yticks=[])
        P.title('S[III]')

        ## Pa 9
        #ax=P.axes([0.02,0.35,0.23,0.27])
        #self.plotaroundline(PaLamb(9))
        #P.setp(ax,xticks=[], yticks=[])
        #P.title('Pa 9')
        ## Pa 10
        ax=P.axes([0.26,0.35,0.23,0.27])
        self.plotaroundline(PaLamb(10))
        P.setp(ax,xticks=[], yticks=[])
        P.title('Pa 10')
        ## Pa 11
        ax=P.axes([0.50,0.35,0.23,0.27])
        self.plotaroundline(PaLamb(11))
        P.setp(ax,xticks=[], yticks=[])
        P.title('Pa 11')
        ## Pa 12
        ax=P.axes([0.74,0.35,0.23,0.27])
        self.plotaroundline(PaLamb(12))
        P.setp(ax,xticks=[], yticks=[])
        P.title('Pa 12')
        
        ## Pa 17
        ax=P.axes([0.02,0.02,0.23,0.27])
        self.plotaroundline(PaLamb(17))
        P.setp(ax,xticks=[], yticks=[])
        P.title('Pa 17')
        
        ## PaStren Ratio
        ax=P.axes([0.28,0.02,0.23,0.27])
        lines=N.array([10,11,12,14,17])
        meas=N.zeros(len(lines),'Float32')
        j=0
        for i in N.arange(5)*2 + 3:
            meas[j]=self.paparams[i]
            j+=1
            
        ratio=meas / PaschStren[19-lines]
        P.plot(lines[::-1],ratio/ratio[-2],'bo')
        P.setp(ax,xticks=[9,10,11,12,14,17])
        P.title('Pa Strength Ratio')

        ## values
        ax=P.axes([0.86,0.02,0.10,0.27])
        
        P.text(0.1,0.9,'S/N: '+str(self.sn.mean()),transform = ax.transAxes)
        P.text(0.1,0.8,'X: %s  Y: %s'%(self.i,self.j),transform = ax.transAxes)
        P.text(0.1,0.7,'PaNumb: %s'%(self.PaNumb,),transform = ax.transAxes)
        P.text(0.1,0.6,'Tilt %s'%(self.tiltfac,),transform = ax.transAxes)
        P.text(0.1,0.5,'Fact: %s'%(self.fact),transform = ax.transAxes)
        P.text(0.1,0.4,'Shift: %s'%(self.shift),transform = ax.transAxes)
        P.text(0.1,0.3,'Flag: %s'%(self.flag),transform = ax.transAxes)
        
        P.setp(ax,xticks=[], yticks=[])
        P.title('Some Values')
        
        self.canvas.draw()

    def plotaroundline(self,lamb):
        region=[vel2lamb(-self.velRange /2.,lamb),vel2lamb(self.velRange /2.,lamb)]
        #print lamb,region, self.velRange, self.guessV
        plotspec(self.shiftscaled(),region=region,style='r',linestyle='steps')
        plotspec(self.data,region=region,style='k',linestyle='steps')

    def button_press_callback(self,event):
        pass
        
def createPa(paschenparam,Z,double,PaNumb,D1=0.0,D2=0.0):
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
        if double:
            para[5]*=Stren[i]
            para[4]=lamb2pix(Paschen[i]*Z,Lamb0,Step)+D2
            SynthSpec+=G.twogauss(para,x=x,returnmodel=True)
        else:
            SynthSpec+=G.gauss(para,x=x,returnmodel=True)

    return SynthSpec


def createPaschen(data,double=True,velRange=None,guessV=None,plot=False,plotfit=False,PaNumb=9):
    fitresults=findLine(data,double=double,velRange=velRange,guessV=guessV,restlamb=PaLamb(PaNumb),plot=plotfit)
    #print fitresults
    if fitresults==-1:
        return N.zeros(SpecLen,'Float32')
    else:
        Z,paschenparam,D1,D2=fitresults

    SynthSpec=createPa(paschenparam,Z,D1=0.0,D2=0.0,double=double,PaNumb=PaNumb)
    
    if plot:    
        plotspec(SynthSpec)
        plotspec(data)
        plotspec(data-SynthSpec,Z=Z,region='cat',plotlines=False)
    
    return SynthSpec

def createPaschenSul(data,velRange=None,guessV=None,plot=False,plotfit=False,PaNumb=9):

    
    fitresults=findLine(data,velRange=velRange,guessV=guessV,plot=plotfit)
    if fitresults==-1:
        return N.zeros(SpecLen,'Float32')
    else:
        Z,fitpara,D1,D2=fitresults
            
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

    fitresults=findLine(data,velRange=velRange,guessV=z2vel(Z),restlamb=PaLamb(PaNumb),parinfo=parinfo,plot=plotfit)
    if fitresults==-1:
        return N.zeros(SpecLen,'Float32')
    else:
        Z,paschenparam,D1,D2=fitresults

        
    SynthSpec=createPa(paschenparam,Z,D1=0.0,D2=0.0,double=True,PaNumb=PaNumb)
    
    if plot:    
        plotspec(SynthSpec)
        plotspec(data)
        plotspec(data-SynthSpec,Z=Z,region='cat',plotlines=True)
    
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
    print "This file defines some functions. It is not meant to be executed. Import it instead!"

