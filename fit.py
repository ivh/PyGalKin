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

def errfu(p, x, y, fu, err=None):
    if err: return (y - fu(p,x)) / err
    else: return y- fu(p, x)

def errfu3(p, x, y, fu, err=None, returnmodel=False):
    """
    replicates the data threefold
    """
    n = len(x)
    p = list(p)
    sl = slice(n//2 -1, -n//2 +1)
    x3 = N.concatenate((x-n,x,x+n))[sl]
    model = N.zeros(len(x3),dtype='f')
    model += fu([p[0]/3., p[1]-n]+p[2:],x3)
    model += fu([p[0]/3., p[1]]+p[2:],x3)
    model += fu([p[0]/3., p[1]+n]+p[2:],x3)
    if returnmodel: return model
    y3 = N.concatenate((y,)*3)[sl]
    if err: return (y3 - model) / err
    else: return y3 - model

def startpara(y):
    return [y.min(), y.argmax().astype('f'), y.max()-y.min(), len(y)/10.0, 0.0, 0.0]

def gauss(p, x):
    return p[0] + (p[2] * N.exp( -1* ((x-p[1])**2) / (2*(p[3]**2)) ) )

def twogauss(p, x):
    return p[0] + (p[2] * N.exp( -1* ((x-p[1])**2) / (2*(p[3]**2)) ) ) \
                + (p[5] * N.exp( -1* ((x-p[4])**2) / (2*(p[6]**2)) ) )

def gaussh34(p, x):
    """p0=cont p1=x0 p2=ampl p3=sigma p4=h3 p5=h4"""
    X = (x-p[1])/p[3]
    h3 = ((2*N.sqrt(2)*(X**3)) - (3*N.sqrt(2)*X))/N.sqrt(6)
    h4 = ((4*(X**4)) - (12*(X**2)) + 3) / N.sqrt(24)
    G = N.exp( -1*(X**2) /2) * p[2]
    return p[0]+ G*(1 + (p[4] * h3) + (p[5]*h4))

def dofit(y, x=None, err=None, fu=gaussh34, errfu=errfu, fullout=False, prin=False, p=None):
    if isconstant(y): return -1
    if not x: x=N.arange(len(y),dtype='float64')
    y = y.astype('float64')

    if not p:
        p = startpara(y)
    if fullout or prin:
        p,cov,infd,msg,status = \
            LS(errfu, p, args=(x,y,fu,err), full_output=True)
        if prin: print msg
        return p,cov,infd,msg,status
    else:
        p, status = LS(errfu, p, args=(x,y,fu,err))
        if status not in [1,2,3,4]: return [0.0]*len(p)
        else: return p
