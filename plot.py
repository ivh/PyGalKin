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

sauron=colors.LinearSegmentedColormap('bone  ', _sauron_data, LUTSIZE)


def imshow(X, cmap=sauron, norm=None, aspect=None, interpolation='nearest', alpha=1.0, vmin=None, vmax=None, origin='lower', extent=None):
    P.imshow(N.transpose(X),cmap=cmap, norm=norm, aspect=aspect, interpolation=interpolation, alpha=alpha, vmin=vmin, vmax=vmax, origin=origin, extent=extent)


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
        if vminmax=='sigbased': vmin,vmax=([-1,1]*0.2*relevant.std())+relevant.mean(); print 'huhuuu'
        else: vmin,vmax=relevant.min(),relevant.max()
        P.axis([region[0]*Z,region[1]*Z,vmin,vmax])
    
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
  
