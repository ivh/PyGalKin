"""
fit.py

contains functions to fit and other functions that do it.

"""

from PyGalKin import *

## wrappers before fit:
def shiftnfit(data,fitfunc):
    s=int(data.size/2-data.argmax())
    dat=shift(data,s)
    fit=fitfunc(dat)
    return N.concatenate((fit,N.array([s])))


## FUNCTIONS TO FIT

def gauss(p, fjac=None, x=None, y=None, err=None, returnmodel=False):
    """p0=cont p1=ampl p2=center p3=sigma """
    model = p[0] + (p[2] * N.exp( -1* ((x-p[1])**2) / (2*(p[3]**2)) ) ) 

    #nomin=(x-p[2])**2
    #denom=(p[3]**2) * 2
    #expo=N.exp(-1*nomin/denom)
    #model=expo * p[1]
    #model+=p[0]

    # Non-negative status value means MPFIT should continue, negative means
    # stop the calculation.
    #print p
    #P.plot(model)
    if returnmodel==True:
        return model
    else:
        status = 0
        return([status, (y-model)/err])
 

def twogauss(p, fjac=None, x=None, y=None, err=None, returnmodel=False):
    """p0=cont p1=ampl p2=center p3=sigma """
    model = p[0] + (p[2] * N.exp( -1* ((x-p[1])**2) / (2*(p[3]**2)) ) )  + (p[5] * N.exp( -1* ((x-p[4])**2) / (2*(p[6]**2)) ) )

    #P.plot(model)
    if returnmodel==True:
        return model
    else:
        status = 0
        return([status, (y-model)/err])
 

def gaussh34(p, fjac=None, x=None, y=None, err=None, returnmodel=False):
    """p0=cont p1=x0 p2=ampl p3=sigma p4=h3 p5=h4"""
    X=(x-p[1])/p[3]
    h3=((2*N.sqrt(2)*(X**3)) - (3*N.sqrt(2)*X))/N.sqrt(6)
    h4=((4*(X**4)) - (12*(X**2)) + 3) / N.sqrt(24)
    G= N.exp( -1*(X**2) /2) * p[2] 
    model = p[0]+ G*(1 + (p[4] * h3) + (p[5]*h4))

    #print model.shape,x.shape,y.shape
    if returnmodel==True:
        return model
    else:
        status = 0
        return([status, (y-model)/err]) # for mpfit
        #return y-model

def gauss_overlap(p, fjac=None, x=None, y=None, err=None):
  """This is used by gauss_from_array(). It returnes a status flag and
      the weighted error of a gaussian fit to a data set.
      """
  # Calculate the exponents first
  arr_len = len(x)
  limit = -700
  #print p
  temp1 = N.maximum( -((x-p[1])**2)/(2*p[3]**2) , limit)
  temp1_right = N.maximum( -(((x+arr_len)-p[1])**2)/(2*p[3]**2) , limit)
  temp1_left = N.maximum( -(((x-arr_len)-p[1])**2)/(2*p[3]**2) , limit)
    
  # Compute the double gauss
  arr = p[0] + p[2]*N.exp(temp1)
  arr_left = p[2]*N.exp(temp1_left)
  arr_right = p[2]*N.exp(temp1_right)
  
  model_y = arr + arr_right + arr_left
  # We are not using this at the moment, 0 means successful computation
  status=0
  #P.plot(y,'b')
  #P.plot(err*(y-model_y),'r')
  # Return the status and the weighted error in the fit
  return [status, err*(y-model_y)]

def twogauss_overlap(p, fjac=None, x=None, y=None, err=None):
  """This is used by gauss_from_array(). It returnes a status flag and
      the weighted error of a gaussian fit to a data set.
      """
  # Calculate the exponents first
  arr_len = len(x)
  
  temp1 = -((x-p[1])**2)/(2*p[3]**2)
  temp2 = -((x-p[4])**2)/(2*p[6]**2)
  temp1_right = -(((x+arr_len)-p[1])**2)/(2*p[3]**2)
  temp2_right = -(((x+arr_len)-p[4])**2)/(2*p[6]**2)
  temp1_left = -(((x-arr_len)-p[1])**2)/(2*p[3]**2) 
  temp2_left = -(((x-arr_len)-p[4])**2)/(2*p[6]**2) 
    
  # Compute the double gauss
  arr = p[0] + p[2]*N.exp(temp1) + p[5]*N.exp(temp2)
  arr_left = p[2]*N.exp(temp1_left) + p[5]*N.exp(temp2_left)
  arr_right = p[2]*N.exp(temp1_right) + p[5]*N.exp(temp2_right)
  
  model_y = arr + arr_right + arr_left
  # We are not using this at the moment, 0 means successful computation
  status=0
  
  # Return the status and the weighted error in the fit
  return [status, err*(y-model_y)]

  


### FUNCTIONS THAT DO THE FITTING WORK

def fitgauss(data,err=None,parinfo=None,prin=False,plot=False,quiet=True):
    if isconstant(data):
        return -1

    data=data.astype('Float64')
    #data-=min(data)
    x=N.arange(len(data),dtype='Float64')
    #err=N.zeros(len(data))+1
    if err==None: err=1/N.sqrt(data)
    
    fa = {'x':x, 'y':data, 'err':err}

    if parinfo==None:
        parinfo=[]
        for i in range(4):
            parinfo.append({'value':0.0, 'fixed':0, 'limited':[0,0],'limits':[0.0, 0.0], 'step':0.0})

#        parinfo[0]['value']=min(data)
        parinfo[0]['value']=1.0
        parinfo[0]['fixed']=0
        parinfo[0]['limited']=[0,0]
        parinfo[0]['limits']=[min(data),max(data)]
        parinfo[1]['value']=N.argmax(data)
        parinfo[1]['limited']=[0,0]
        parinfo[1]['limits']=[0.0,len(data)]
        parinfo[2]['value']=(max(data)-min(data))
        parinfo[2]['limited']=[0,0]
        parinfo[2]['limits']=[0.0,max(data)]
        parinfo[3]['value']=len(data)/16.
        parinfo[3]['limited']=[0,0]
        parinfo[3]['limits']=[0.0,len(data)/2.]
        

    try:
        fit=mpfit(gauss,functkw=fa,parinfo=parinfo,maxiter=200,quiet=quiet)
    except OverflowError:
        return -1
    
    if plot==True:
        P.clf()
        P.plot(data,'r',linestyle='steps')
        P.plot(gauss(fit.params,x=N.arange(len(data)),returnmodel=True),'b')
    if prin==True:
        print fit.niter,fit.params,fit.status
    
    return fit


def fit2gauss(data,parinfo=None,plot=False,prin=False,quiet=True,fitfunc=None,x=None,):
    if isconstant(data):
        return -1

    data=data.astype('Float64')
    #data-=min(data)
    #err=N.zeros(len(data))+1
    err=1/N.sqrt(data)
    
    if x == None: x=N.arange(len(data),dtype='Float64')
    
    fa = {'x':x, 'y':data, 'err':err}

    if parinfo==None:
        parinfo=[]
        for i in range(7):
            parinfo.append({'value':0.0, 'fixed':0, 'limited':[0,0],'limits':[0.0, 0.0], 'step':0.0})

        parinfo[0]['value']=min(data)
        parinfo[0]['limited']=[1,1]
        parinfo[0]['limits']=[min(data),max(data)]
        parinfo[1]['value']=N.argmax(data)
        parinfo[1]['limited']=[1,1]
        parinfo[1]['limits']=[0.0,len(data)]
        parinfo[2]['value']=(max(data)-min(data))
        parinfo[2]['limited']=[1,1]
        parinfo[2]['limits']=[0.0,max(data)]
        parinfo[3]['value']=len(data)/6.
        parinfo[3]['limited']=[1,1]
        parinfo[3]['limits']=[0.0,len(data)/2.]
        parinfo[4]['value']=N.argmax(data)
        parinfo[4]['limited']=[1,1]
        parinfo[4]['limits']=[0.0,len(data)]
        parinfo[5]['value']=0.0
        parinfo[5]['limited']=[1,1]
        parinfo[5]['limits']=[0,max(data)]
        parinfo[6]['value']=len(data)/2.
        parinfo[6]['limited']=[1,1]
        parinfo[6]['limits']=[0.0,len(data)/2.]
    else:
        parinfo[1]['value']+=len(data)/2
        parinfo[4]['value']+=len(data)/2
      

    print data,x,err,parinfo
    if fitfunc==None: fitfunc=twogauss
    try:
        fit=mpfit(fitfunc,functkw=fa,parinfo=parinfo,maxiter=200,quiet=quiet)
    except OverflowError:
        return -1
        
    if plot==True:
        P.clf()
        P.plot(data,'r',linestyle='steps')
        P.plot(twogauss(fit.params,x=N.arange(len(data)),returnmodel=True),'g')
        P.plot(gauss(fit.params,x=N.arange(len(data)),returnmodel=True),'b')
        P.plot(gauss(fit.params[[0,4,5,6]],x=N.arange(len(data)),returnmodel=True),'r')
    if prin==True:
        print fit.niter,fit.params,fit.status
    return fit



def fitgaussh34(data,err=None,parinfo=None,prin=False,plot=False,quiet=True):
    if isconstant(data):
        return -1
    
    data=data.astype('Float64')
    x=N.arange(len(data),dtype='Float64')
    if err==None: err=1/N.sqrt(data)
    
    fa = {'x':x, 'y':data, 'err':err}
    #plot(x,data)
    #print type(x),type(data)

    if parinfo==None:
        parinfo=[]
        for i in range(6):
            parinfo.append({'value':0.0, 'fixed':0, 'limited':[0,0],'limits':[0.0, 0.0], 'step':0.0})
        parinfo[0]['value']=min(data)
        parinfo[0]['limited']=[1,1]
        parinfo[0]['limits']=[min(data),max(data)]
        #parinfo[0]['value']=1.0
        #parinfo[0]['fixed']=1
        
        parinfo[1]['value']=N.argmax(data)
        parinfo[1]['limited']=[1,1]
        parinfo[1]['limits']=[0.0,len(data)]

        parinfo[2]['value']=(max(data)-min(data))
        parinfo[2]['limited']=[1,1]
        parinfo[2]['limits']=[0.0,max(data)]

        parinfo[3]['value']=len(data)/16.
        parinfo[3]['limited']=[0,0]
        #parinfo[3]['limits']=[0.0,len(data)/2.]

        parinfo[4]['value']=len(data)/16.
        parinfo[4]['limited']=[0,0]
        #parinfo[4]['limits']=[0.0,len(data)/2.]

        parinfo[5]['value']=len(data)/16.
        parinfo[5]['limited']=[0,0]
        #parinfo[5]['limits']=[0.0,len(data)/2.]



    #print x.shape,data.shape,err.shape
    try:
        fit=mpfit(gaussh34,functkw=fa,parinfo=parinfo,maxiter=200,quiet=quiet)
    except OverflowError:
        fit.params=array([0,0,0,0,0,0],dtype='f')
    
    if plot==True:
        P.clf()
        P.plot(data,'r',linestyle='steps')
        P.plot(gaussh34(fit.params,x=N.arange(len(data)),returnmodel=True),'b')
    if prin==True:
        print fit.niter,fit.params,fit.status
    
    return fit.params

