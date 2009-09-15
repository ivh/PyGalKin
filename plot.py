"""
plot.py
"""

from matplotlib import rcParams, colors
LUTSIZE = rcParams['image.lut']

from tool import *

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

sauron=colors.LinearSegmentedColormap('bone  ', _sauron_data, LUTSIZE)
sauron_inv=colors.LinearSegmentedColormap('bone  ', _sauron_inv_data, LUTSIZE)

def imshow(X, cmap=sauron, norm=None, aspect=None, interpolation='nearest', alpha=1.0, vmin=None, vmax=None, origin='lower', extent=None):
    P.imshow(N.transpose(X),cmap=cmap, norm=norm, aspect=aspect, interpolation=interpolation, alpha=alpha, vmin=vmin, vmax=vmax, origin=origin, extent=extent)



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
      kpc=(x2-x1)*arcsec2kpc(arcs,p['vr0'])/1.E3
  P.text(x1-0.01,y+0.01,r'$%.1f^{\prime\prime}\approx %.1fkpc$'%(arcs,kpc),fontsize=9)
  

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



# MATPLOTLIB EVENT HANDLERS
#
def on_click(event):
  """Handle click events in plotting windows. Not to be called from the command 
      line. It saves the integer coordinates for the clicks with button 1 in a 
      temporary file. Button 3 closes the figure.
      """
  file=open('/tmp/MPclick.dat','a')

  if (event.button == 1):
    if (event.inaxes):
      print event.xdata, event.ydata
      file.write('%d %d\n' % (event.xdata, event.ydata))
  elif (event.button == 2):
    pass
  elif (event.button == 3):
    P.close()

  file.close()

def on_click_float(event):
  """Handle click events in plotting windows. Not to be called from the command 
      line. It saves the float coordinates for the clicks with button 1 in a 
      temporary file. Button 3 closes the figure.
      """
  file=open('/tmp/MPclick.dat','a')

  if (event.button == 1):
    if event.inaxes:
      print event.xdata, event.ydata
      file.write('%f %f\n' % (event.xdata, event.ydata))
  elif (event.button == 2):
    pass
  elif (event.button == 3):
    P.close()

  file.close()
#
# END: MATPLOTLIB EVENT HANDLERS


class inspectdata():
    def __init__(self,data):
        self.data=data
        self.shape=self.data.shape
        if len(self.data.shape)>2:
            self.data.shape=(self.shape[0]*self.shape[1],self.shape[2])
        self.nx=self.data.shape[0]
        self.nz=self.data.shape[1]

        self.ngauss=N.zeros(self.nx,'i')+1
        self.cur=-1

        self.gauss1=N.zeros(self.nx,'i')

        self.fig1=P.figure()
        self.fig1.canvas.mpl_connect('key_press_event',self.keyhandler)     
        self.next()
        
        
    def init1(self):
        self.fig1.clf()
        self.ax1=self.fig1.add_subplot(1,1,1)
        self.ax1.set_xlabel('chanel')
        self.ax1.set_ylabel('intensity')
        self.ax1.set_title('Spectrum %d out of %d'%(self.cur,self.nx))

    def next(self):
        self.cur+=1
        self.init1()
        self.ax1.plot(self.data[self.cur,:],color='k',linestyle='steps')
        self.fig1.canvas.draw()
    
    def fit(self):
        pass

    def keyhandler(self,event):
        if event.key=='s':
            print "Pick the point"
            self.clickconn=self.fig1.canvas.mpl_connect('pick_event',self.plotpicked)
        if event.key=='1': self.ngauss[self.cur]=1;self.next()
        if event.key=='2': self.ngauss[self.cur]=2;self.next()
        if event.key=='3': self.ngauss[self.cur]=3;self.next()
        if event.key=='b': self.ngauss[self.cur]=0;self.next()
        if event.key=='q': self.fig1.canvas.mpl_disconnect(self.clickconn)
 

class rotcur_int():
    def __init__(self,vf,cen=[],pa=90.0,wedge=10.0,incl=45.0,rbin=1.0):
        if cen==[]: cen=N.array(vf.shape)/2.0
        self.vf=vf
        self.cen=N.array(cen)
        self.pa=pa
        self.wedge=wedge
        self.incl=incl
        self.rbin=rbin

        self.fig=P.figure(1)
        self.canvas=self.fig.canvas
        self.ax1=self.fig.add_subplot(2,2,1)
        self.ax2=self.fig.add_subplot(2,1,2)
        self.ax3=self.fig.add_subplot(2,2,2)
        
        self.canvas.mpl_connect('key_press_event',self.keyhandler)
        self.canvas.mpl_connect('button_press_event',self.mousehandler)
        
        self.flip=True
        self.showindiv=True
        
        self.update()

    def update(self):
        self.ax1.clear()
        self.ax2.clear()
        self.ax3.clear()
        self.calc()
        self.plot()
        self.canvas.draw()

    def plot(self):
        self.ax1.imshow(N.transpose(self.vf),interpolation='nearest',origin='lower',cmap=sauron)
        vec=N.array([-N.sin(N.radians(self.pa)),N.cos(N.radians(self.pa))])
        len=self.vf.shape[0]/2
        x=[self.cen[0]-len*vec[0],self.cen[0]+len*vec[0]]
        y=[self.cen[1]-len*vec[1],self.cen[1]+len*vec[1]]
        self.ax1.plot(x,y,'-k')
        vec=N.array([-N.sin(N.radians(self.pa - self.wedge)),N.cos(N.radians(self.pa - self.wedge))])
        x=[self.cen[0]-len*vec[0],self.cen[0]+len*vec[0]]
        y=[self.cen[1]-len*vec[1],self.cen[1]+len*vec[1]]
        self.ax1.plot(x,y,'--k')
        vec=N.array([-N.sin(N.radians(self.pa+self.wedge)),N.cos(N.radians(self.pa+self.wedge))])
        x=[self.cen[0]-len*vec[0],self.cen[0]+len*vec[0]]
        y=[self.cen[1]-len*vec[1],self.cen[1]+len*vec[1]]
        self.ax1.plot(x,y,'--k')
        self.ax1.axis([-0.5,self.vf.shape[0]-0.5,-0.5,self.vf.shape[1]-0.5])

        if self.showindiv:
            self.ax2.plot(self.r1,self.v1,'r1')
            self.ax2.plot(self.r2,self.v2,'g1')
        self.ax2.errorbar(self.rr,self.vr,self.sr,fmt='go')
        self.ax2.errorbar(self.rl,self.vl,self.sl,fmt='ro')
        self.ax2.grid(alpha=0.2,ls='--')

        self.ax3.text(0.05,0.9,'cen: %.1f %.1f'%(self.cen[0],self.cen[1]))
        self.ax3.text(0.05,0.8,'pa: %d'%self.pa)
        self.ax3.text(0.05,0.7,'wedge: %d'%self.wedge)
        self.ax3.text(0.05,0.6,'incl: %d'%self.incl)
        self.ax3.text(0.05,0.5,'rbin: %.1f'%self.rbin)
        self.ax3.text(0.05,0.4,'vsys: %.1f'%self.vsys)
        self.ax3.set_xticks([])
        self.ax3.set_yticks([])
        

    def calc(self):
        self.r1,self.r2,self.v1,self.v2=\
            rotcur(self.vf,self.cen,self.pa,\
                       self.wedge,self.incl)
        self.rl,self.vl,self.sl=binRC(self.r1,self.v1,self.rbin)
        self.rr,self.vr,self.sr=binRC(self.r2,self.v2,self.rbin)
        if self.flip:
            self.rl*=-1
            self.r1*=-1
        else:
            self.vl*=-1
            self.v1*=-1
        
        self.vsys=self.vf[self.cen[0],self.cen[1]]

    def keyhandler(self,event):
        if event.key=='c': self.pa+=2; self.update()
        elif event.key=='v': self.pa-=2; self.update()
        elif event.key=='i': self.incl+=5; self.update()
        elif event.key=='o': self.incl-=5; self.update()
        elif event.key=='w': self.wedge+=5; self.update()
        elif event.key=='e': self.wedge-=5; self.update()
        elif event.key=='b': self.rbin*=1.5; self.update()
        elif event.key=='n': self.rbin/=1.5; self.update()
        elif event.key=='l': self.flip=not self.flip; self.update()
        elif event.key=='k': self.showindiv=not self.showindiv; self.update()
        elif event.key=='left': self.cen[0]-=1; self.update()
        elif event.key=='right': self.cen[0]+=1; self.update()
        elif event.key=='up': self.cen[1]+=1; self.update()
        elif event.key=='down': self.cen[1]-=1; self.update()
        else: pass#print event.key

     

    def mousehandler(self,event):
        self.cen=N.array([event.xdata,event.ydata])
        self.update()

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
