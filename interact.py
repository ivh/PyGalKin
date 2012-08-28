"""
interact.py

A few classes that allow interaction with data
using matplotlibs key- and mouse-bindings.

"""


from PyGalKin import *

class measure(object):
    def __init__(self,g,img,vf,sig,fit,cube):
        self.g = g
        self.img = img
        self.vf = vf
        self.sig = sig
        self.fit = fit
        self.cube=cube
        self.fig = P.figure(figsize=(12,9))
        self.canvas = self.fig.canvas
        self.ax1 = self.fig.add_subplot(231)
        self.ax2 = self.fig.add_subplot(232, sharex=self.ax1, sharey=self.ax1)
        self.ax3 = self.fig.add_subplot(233, sharex=self.ax1, sharey=self.ax1)
        self.ax4 = self.fig.add_subplot(234, sharex=self.ax1, sharey=self.ax1)
        self.ax5 = self.fig.add_subplot(235, sharex=self.ax1, sharey=self.ax1)
        self.axspec = self.fig.add_subplot(4,6,18)
        self.axpv = self.fig.add_subplot(4,3,12)
        self.fig.subplots_adjust(wspace=0.01,hspace=0.01,left=0.01,right=0.99,top=0.99,bottom=0.01)
        self.imaxes = [self.ax1, self.ax2, self.ax3, self.ax4, self.ax5]
        for ax in self.imaxes: ax.set_xticks(()); ax.set_yticks(())
        self.canvas.mpl_connect('key_press_event', self.key)
        self.canvas.mpl_connect('button_press_event', self.button)

        self.avreg = 5
        self.cx, self.cy = (self.g.p['relsize'],)*2
        self.px = 0
        self.py = 0
        self.pixdiff2radius= 1

        self.ax1.imshow(N.transpose(self.img))
        self.ax2.imshow(N.transpose(self.vf))
        self.ax3.imshow(N.transpose(self.sig))
        self.ax4.imshow(N.transpose(self.fit[:,:,4]),vmin=-0.5,vmax=0.5)
        self.ax5.imshow(N.transpose(self.fit[:,:,5]),vmin=-0.5,vmax=0.5)
        self.draw_spec()
        self.draw_pvdiag()
        self.draw_markers()

        if self.g.p.get('axis'):
            self.ax1.axis(self.g.p.get('axis'))

        self.canvas.set_window_title('%s'%self.g.name)
        P.P.show()

    def getAvRange(self,previous=False):
        if previous: cx,cy=self.px,self.py
        else: cx,cy=self.cx,self.cy
        av = self.avreg //2
        extra = self.avreg %2
        xav = slice(cx - av, cx + av + extra)
        yav = slice(cy - av, cy + av + extra)
        return xav,yav

    def meas_sigma(self):
        xav,yav=self.getAvRange()
        self.g.sigma_x, self.g.sigma_y = self.cx, self.cy
        self.g.sigma_cent = N.mean(self.sig[xav,yav])
        self.g.mass_sig = dynMassSphere((self.g.h27 or self.g.h25)*1.68,self.g.sigma_cent)
        print "Measured sigma: %s   -> mass=%.1e)"%\
            (self.g.sigma_cent,self.g.mass_sig)
        self.g.save()

    def meas_vels(self):
        self.g.maxvel_x1, self.g.maxvel_y1 = self.px,self.py
        self.g.maxvel_x2, self.g.maxvel_y2 = self.cx,self.cy
        xav,yav=self.getAvRange()
        v1=N.mean(self.vf[xav,yav])
        xav,yav=self.getAvRange(previous=True)
        v2=N.mean(self.vf[xav,yav])
        self.g.maxvel=N.abs(v1-v2)

        scale = scalefromvarc(self.g.p['echelle'],self.g.vsys) / 1000.
        npix = dis(N.array([self.cx,self.cy]),N.array([self.px,self.py]))
        self.g.mass_p2p = massKepler(npix * scale, self.g.maxvel)
        print "|v1-v2|=%.1f, r=%.1f kpc -> mass=%.1e"%\
            (self.g.maxvel, npix*scale, self.g.mass_p2p)
        self.g.save()

    def add_poi(self):
        if not 'poi' in self.g.p:
            self.g.p['poi'] = []
        self.g.p['poi'].append((self.cx,self.cy))
        print "adding %s to pois."%(self.g.p['poi'][-1],)
        self.g.save()

    def button(self,event):
        if event.inaxes not in self.imaxes: return
        x,y=int(round(event.xdata)),int(round(event.ydata))
        if event.button==1:
            self.px,self.py=self.cx,self.cy
            self.cx,self.cy=x,y
            print 'x=%d, y=%d, v=%.1f, s=%.1f'%(x,y,self.vf[x,y],self.sig[x,y])
            self.draw_spec()

    def draw_spec(self):
        self.axspec.clear()
        spec=self.cube[self.cx,self.cy,:]
        self.axspec.plot(spec,'b',linestyle='steps')
        fit=F.gaussh34(self.fit[self.cx,self.cy], N.arange(len(spec)))
#                """p0=cont p1=x0 p2=ampl p3=sigma p4=h3 p5=h4
        self.axspec.plot(fit,'g-',linestyle='steps')
        self.canvas.draw()


    def draw_pvdiag(self):
        self.axpv.clear()
        pos,vel=posvel(self.vf,self.g.p['dyncen'],self.g.p['pa'])
        r1,r2,v1,v2=rotcur(self.vf,self.g.p['dyncen'],self.g.p['pa'],self.g.p['wedge'],self.g.incl)
        self.axpv.plot(pos,vel-self.g.vsys,'.y',alpha=0.6)
        self.axpv.plot(-r1,v1,'sr')
        self.axpv.plot(r2,v2,'sb')
        self.canvas.draw()

    def draw_markers(self):
        self.ax2.plot(*self.g.dyncen,color='k',marker='o')
        self.ax2.plot([self.g.maxvel_x1,self.g.maxvel_x2],
            [self.g.maxvel_y1,self.g.maxvel_y2],'Dk')

        vec=pa2vec(self.g.pa)*self.g.p['relsize']
        ends=N.transpose(N.vstack((self.g.dyncen+vec,self.g.dyncen-vec)))
        self.ax2.plot(ends[0],ends[1],'k-')

        self.ax3.plot(self.g.sigma_x,self.g.sigma_y,'ok')
        self.canvas.draw()

    def key(self,event):
        if event.key=='S':
            self.meas_sigma()
            self.draw_markers()
        elif event.key=='v':
            self.meas_vels()
        elif event.key=='i':
            self.add_poi()
        elif event.key=='r':
            self.g.p['poi'] = []
            self.g.save()
        elif event.key=='w':
            self.g.p['wedge'] = input('new wedge: ')
            self.g.save()
        elif event.key=='a':
            self.g.pa = input('new PA: ')
            self.g.save()
            self.draw_markers()
        elif event.key=='c':
            self.g.p['dyncen']=[self.cx,self.cy]
            self.g.save()
            self.draw_pvdiag()
            self.draw_markers()
        elif event.key in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
            self.avreg = int(event.key)
            print 'Averaging region set to %s pixels.'%self.avreg
        elif event.key=='z':
            self.g.p['axis']=self.ax2.axis()
            self.g.save()
        elif event.key=='q':
            sys.exit()

#############################################
class doublecomp:
    def __init__(self,fitresults=None,masks=None,cube=None,z=None,vmin=None,vmax=None,\
                     clipcube=True,extra=1.1,restl=Sulfur,savefile='Haro11-multicomp.pick'):
        self.cx=12
        self.cy=9
        self.extra=extra
        self.restl=restl
        self.savefile=savefile
        self.vmin=vmin
        self.vmax=vmax
        self.l0=vel2lamb(vmin/extra,restl)
        self.l1=vel2lamb(vmax*extra,restl)

        if cube!=None:
            if clipcube:
                p0=lamb2pix(self.l0,cube.l0,cube.dl)-1
                p1=lamb2pix(self.l1,cube.l0,cube.dl)
                cube=cube[:,:,p0:p1]

            self.z=cube.shape[-1]
            cont=C.contfrommin(cube,10)
            cont.shape=(cube.shape[0],cube.shape[1],1)
            self.cube=cube - cont

        elif z: self.z=z
        else: self.z=100

        self.zv=N.arange(self.z,dtype='float64')/self.z
        self.zv*=self.vmax*extra - self.vmin/extra
        self.zv+=self.vmin/extra

        if fitresults != None:
            a1,v1,s1,a2,v2,s2,a3,v3,s3=fitresults
            #self.a1,self.v1,self.s1=masked_array(a1),masked_array(v1),masked_array(s1)
            #self.a2,self.v2,self.s2=masked_array(a2),masked_array(v2),masked_array(s2)
            #self.a3,self.v3,self.s3=masked_array(a3),masked_array(v3),masked_array(s3)
            self.a1,self.v1,self.s1=a1,v1,s1
            self.a2,self.v2,self.s2=a2,v2,s2
            self.a3,self.v3,self.s3=a3,v3,s3
        else:
            z=N.zeros((self.cube.shape[0],self.cube.shape[1]),dtype='float64')
            self.a1,self.v1,self.s1=masked_array(z,copy=True),masked_array(z,copy=True),masked_array(z,copy=True)
            self.a2,self.v2,self.s2=masked_array(z,copy=True),masked_array(z,copy=True),masked_array(z,copy=True)
            self.a3,self.v3,self.s3=masked_array(z,copy=True),masked_array(z,copy=True),masked_array(z,copy=True)
            #self.all2()

        if masks==None:
            self.mask1=N.zeros_like(self.a1).astype('bool')
            self.mask2=N.zeros_like(self.a2).astype('bool')
            self.mask3=N.zeros_like(self.a2).astype('bool')
        else:
            self.mask1=masks[0]
            self.mask2=masks[1]
            self.mask3=masks[2]

        #self.sortbyamp()

        self.fig=P.figure(1)
        self.fig.clf()
        self.canvas=self.fig.canvas
        self.ax1=self.fig.add_subplot(3,3,1)
        self.ax2=self.fig.add_subplot(3,3,2)
        self.ax3=self.fig.add_subplot(3,3,4)
        self.ax4=self.fig.add_subplot(3,3,5)
        self.ax5=self.fig.add_subplot(3,3,7)
        self.ax6=self.fig.add_subplot(3,3,8)
        self.ax7=self.fig.add_subplot(2,3,3)

        self.canvas.mpl_connect('key_press_event', self.key)
        self.canvas.mpl_connect('button_press_event', self.button)

        self.setregion()
        self.applymask()
        self.newgausses()
        self.plotvfs()
        self.plotcurspec()
        #self.canvas.draw()
        self.update()

    def button(self,event):
        if event.inaxes in [self.ax1,self.ax2,self.ax3,self.ax4,self.ax5,self.ax6]:
            x,y=int(round(event.xdata)),int(round(event.ydata))
            if event.button==1:
                self.cx,self.cy=x,y
                self.update()
            elif event.button==2:
                self.setregion(x,y)
                self.update()
            elif event.button==3:
                self.setregion(x,y)
                self.update()
        elif event.inaxes in [self.ax7]:
            x,y=event.xdata,event.ydata
            if event.button==1:
                self.vcl1,self.acl1=x,y
                self.parinfo[1]['value']=x
                self.parinfo[2]['value']=y
            elif event.button==2:
                self.vcl2,self.acl2=x,y
                self.parinfo[4]['value']=x
                self.parinfo[5]['value']=y
            elif event.button==3:
                self.vcl3,self.acl3=x,y
                self.parinfo[7]['value']=x
                self.parinfo[8]['value']=y

    def key(self,event):
        print event.key
        if event.key=='x':
            IO.dump((self.mask1,self.mask2,self.mask3,self.a1,self.v1,self.s1,self.a2,self.v2,self.s2,self.a3,self.v3,self.s3),self.savefile)
        elif event.key=='w':
            self.switch(); self.update()
        elif event.key=='e':
            self.switch2(); self.update()
        elif event.key=='b':
            self.maskregion(1); self.update()
        elif event.key=='n':
            self.maskregion(2); self.update()
        elif event.key=='m':
            self.maskregion(3); self.update()
        elif event.key=='t':
            self.unmaskregion(1); self.update()
        elif event.key=='y':
            self.unmaskregion(2); self.update()
        elif event.key=='u':
            self.unmaskregion(3); self.update()
        elif event.key=='i':
            self.cy+=1; self.update()
        elif event.key=='k':
            self.cy-=1; self.update()
        elif event.key=='j':
            self.cx-=1; self.update()
        elif event.key=='l':
            self.cx+=1; self.update()
        elif event.key=='1':
            self.fit1(); self.unmaskregion(); self.maskregion(2); self.maskregion(3); self.update()
        elif event.key=='2':
            self.fit2(); self.unmaskregion(); self.maskregion(3); self.update()
        elif event.key=='3':
            self.fit3(); self.unmaskregion(); self.update()
        elif event.key=='8':
            self.parinfo[1]['fixed']=1
            self.parinfo[4]['fixed']=1

    def fit1(self,x=None,y=None):
        if not x: x=self.cx
        if not y: y=self.cy
        data=N.array(self.cube[x,y,:])+1
        fit=F.fitgauss(data,x=self.zv,parinfo=self.parinfo)
        if fit != -1:
            c,v,a,s,c,c,c,c,c,c=fit.params
        else: a,v,s=0,0,0
        self.a1[x,y]=a
        self.v1[x,y]=v
        self.s1[x,y]=s
        self.a2[x,y]=0
        self.v2[x,y]=0
        self.s2[x,y]=0
        print('Fitted x:%d y:%d with a:%e v:%e s:%e'%(x,y,a,v,s))

    def fit2(self,x=None,y=None):
        if not x: x=self.cx
        if not y: y=self.cy
        data=N.array(self.cube[x,y,:])+1
        fit=F.fit2gauss(data,x=self.zv,parinfo=self.parinfo,prin=True)
        if fit != -1:
            c,v1,a1,s1,v2,a2,s2,c,c,c=fit.params
        else: v1,a1,s1,v2,a2,s2=(0.0,)*6
        self.a1[x,y]=a1
        self.v1[x,y]=v1
        self.s1[x,y]=s1
        self.a2[x,y]=a2
        self.v2[x,y]=v2
        self.s2[x,y]=s2
        print('Fitted x:%d y:%d with a1:%e v1:%e s1:%e a2:%e v2:%e s2:%e'%(x,y,a1,v1,s1,a2,v2,s2))

    def fit3(self,x=None,y=None):
        if not x: x=self.cx
        if not y: y=self.cy
        data=N.array(self.cube[x,y,:])+1
        fit=F.fit3gauss(data,x=self.zv,parinfo=self.parinfo,prin=True)
        if fit != -1:
            c,v1,a1,s1,v2,a2,s2,v3,a3,s3=fit.params
        else: v1,a1,s1,v2,a2,s2,v3,a3,s3=(0.0,)*9
        self.a1[x,y]=a1
        self.v1[x,y]=v1
        self.s1[x,y]=s1
        self.a2[x,y]=a2
        self.v2[x,y]=v2
        self.s2[x,y]=s2
        self.a3[x,y]=a3
        self.v3[x,y]=v3
        self.s3[x,y]=s3
        print('Fitted x:%d y:%d with a1:%e v1:%e s1:%e a2:%e v2:%e s2:%e a3:%e v3:%e s3:%e'%(x,y,a1,v1,s1,a2,v2,s2,a3,v3,s3))


    def make_parinfo(self,x,data):
        parinfo=[]
        for i in range(10):
            parinfo.append({'value':0.0, 'fixed':0, 'limited':[1,0],'limits':[0.0, 0.0], 'step':0.0})
        parinfo[0]['value']=data.min()
        parinfo[1]['value']=x[data.argmax()]
        parinfo[2]['value']=data.max()
        parinfo[3]['value']=(x.max()-x.min())/10.
        parinfo[4]['value']=x[data.argmax()]
        parinfo[5]['value']=data.max()/2
        parinfo[6]['value']=(x.max()-x.min())/10.
        parinfo[7]['value']=x[data.argmax()]
        parinfo[8]['value']=data.max()/2
        parinfo[9]['value']=(x.max()-x.min())/10.

        #parinfo[1]['fixed']=1
        #parinfo[4]['fixed']=1
        return parinfo

    def all1(self):
        for i in N.arange(self.cube.shape[0]):
            for j in N.arange(self.cube.shape[1]):
                self.fit1(i,j)

    def all2(self):
        for i in N.arange(self.cube.shape[0]):
            for j in N.arange(self.cube.shape[1]):
                self.fit2(i,j)

    def setregion(self,x=None,y=None):
        if not x: x=self.cx
        if not y: y=self.cy
        self.rx=[x,self.cx]
        self.ry=[y,self.cy]
        self.rx.sort()
        self.ry.sort()

    def maskregion(self,what=4):
        if what==1:
            self.mask1[self.rx[0]:self.rx[1]+1,self.ry[0]:self.ry[1]+1]=True
        elif what==2:
            self.mask2[self.rx[0]:self.rx[1]+1,self.ry[0]:self.ry[1]+1]=True
        elif what==3:
            self.mask3[self.rx[0]:self.rx[1]+1,self.ry[0]:self.ry[1]+1]=True

    def unmaskregion(self,what=4):
        if what==1:
            self.mask1[self.rx[0]:self.rx[1]+1,self.ry[0]:self.ry[1]+1]=False
        elif what==2:
            self.mask2[self.rx[0]:self.rx[1]+1,self.ry[0]:self.ry[1]+1]=False
        elif what==3:
            self.mask2[self.rx[0]:self.rx[1]+1,self.ry[0]:self.ry[1]+1]=False
        else:
            for mask in [self.mask1,self.mask2,self.mask3]:
                mask[self.rx[0]:self.rx[1]+1,self.ry[0]:self.ry[1]+1]=False

    def sortbyamp(self):
        self.switchbymask(self.a1<self.a2)


    def switchbymask(self,mask):
        self.a1,self.a2=N.where(mask,self.a2,self.a1),N.where(mask,self.a1,self.a2)
        self.v1,self.v2=N.where(mask,self.v2,self.v1),N.where(mask,self.v1,self.v2)
        self.s1,self.s2=N.where(mask,self.s2,self.s1),N.where(mask,self.s1,self.s2)


    def switch(self):
        sl=slice(self.rx[0],self.rx[1]+1),slice(self.ry[0],self.ry[1]+1)
        self.a1[sl],self.a2[sl]= self.a2[sl].copy(),self.a1[sl].copy()
        self.v1[sl],self.v2[sl]= self.v2[sl].copy(),self.v1[sl].copy()
        self.s1[sl],self.s2[sl]= self.s2[sl].copy(),self.s1[sl].copy()
        self.mask1[sl],self.mask2[sl]= self.mask2[sl].copy(),self.mask1[sl].copy()

    def switch2(self):
        sl=slice(self.rx[0],self.rx[1]+1),slice(self.ry[0],self.ry[1]+1)
        self.a3[sl],self.a2[sl]= self.a2[sl].copy(),self.a3[sl].copy()
        self.v3[sl],self.v2[sl]= self.v2[sl].copy(),self.v3[sl].copy()
        self.s3[sl],self.s2[sl]= self.s2[sl].copy(),self.s3[sl].copy()
        self.mask3[sl],self.mask2[sl]= self.mask2[sl].copy(),self.mask3[sl].copy()


    def plotcurspec(self,x=0,y=0):
        self.ax7.set_title('Current pixel')
        if hasattr(self,'cube'):
            cub=self.cube[self.cx,self.cy,:]
            self.meas,=self.ax7.plot(self.zv,cub,'k-',linestyle='steps',label='meas')
        self.gauss1,=self.ax7.plot(self.zv,self.g1,'b-',label='1')
        self.gauss2,=self.ax7.plot(self.zv,self.g2,'g-',label='2')
        self.gauss3,=self.ax7.plot(self.zv,self.g3,'y-',label='3')
        self.gauss,=self.ax7.plot(self.zv,self.g,'r-',label='sum')
        self.ax7.set_yticks([])
        self.ax7.set_xbound(self.vmin/self.extra,self.vmax*self.extra)
        #self.ax7.legend(loc=2)

    def plotcross(self):
        crx,cry=[self.cx-0.5,self.cx+0.5],[self.cy-0.5,self.cy+0.5]
        csx,csy=[self.cx-0.5,self.cx+0.5],[self.cy+0.5,self.cy-0.5]
        self.cross1a,=self.ax1.plot(crx,cry,'w-')
        self.cross1b,=self.ax1.plot(csx,csy,'k-')
        self.cross2a,=self.ax2.plot(crx,cry,'w-')
        self.cross2b,=self.ax2.plot(csx,csy,'k-')
        self.cross3a,=self.ax3.plot(crx,cry,'w-')
        self.cross3b,=self.ax3.plot(csx,csy,'k-')
        self.cross4a,=self.ax4.plot(crx,cry,'w-')
        self.cross4b,=self.ax4.plot(csx,csy,'k-')
        self.cross5a,=self.ax5.plot(crx,cry,'w-')
        self.cross5b,=self.ax5.plot(csx,csy,'k-')
        self.cross6a,=self.ax6.plot(crx,cry,'w-')
        self.cross6b,=self.ax6.plot(csx,csy,'k-')


    def plotvfs(self):
        self.vf1=self.ax1.imshow(N.transpose(self.v1),interpolation='nearest',vmin=self.vmin,vmax=self.vmax,cmap=P.sauron)
        self.vf2=self.ax2.imshow(N.transpose(self.v2),interpolation='nearest',vmin=self.vmin,vmax=self.vmax,cmap=P.sauron)
        self.amp1=self.ax3.imshow(N.transpose(self.a1),interpolation='nearest',vmin=0.0,vmax=self.a1.max(),cmap=P.sauron_inv)
        self.amp2=self.ax4.imshow(N.transpose(self.a2),interpolation='nearest',vmin=0.0,vmax=self.a1.max(),cmap=P.sauron_inv)
        self.sig1=self.ax5.imshow(N.transpose(self.s1),interpolation='nearest',vmin=self.s1.min(),vmax=self.s1.max(),cmap=P.sauron_inv)
        self.sig2=self.ax6.imshow(N.transpose(self.s2),interpolation='nearest',vmin=self.s1.min(),vmax=self.s1.max(),cmap=P.sauron_inv)
        ax=self.ax2.axis()
        self.plotcross()
        self.ax1.axis(ax)
        self.ax2.axis(ax)
        self.ax3.axis(ax)
        self.ax4.axis(ax)
        self.ax5.axis(ax)
        self.ax6.axis(ax)
        self.ax1.set_title('VF 1')
        self.ax2.set_title('VF 2')
        self.ax3.set_title('Amp 1')
        self.ax4.set_title('Amp 2')
        self.ax5.set_title('Sig 1')
        self.ax6.set_title('Sig 2')
        self.ax1.set_yticks([]);self.ax1.set_xticks([])
        self.ax2.set_yticks([]);self.ax2.set_xticks([])
        self.ax3.set_yticks([]);self.ax3.set_xticks([])
        self.ax4.set_yticks([]);self.ax4.set_xticks([])
        self.ax5.set_yticks([]);self.ax5.set_xticks([])
        self.ax6.set_yticks([]);self.ax6.set_xticks([])

    def newgausses(self):
        self.g1=F.gauss([0.0,self.v1[self.cx,self.cy],self.a1[self.cx,self.cy],self.s1[self.cx,self.cy]],x=self.zv, returnmodel=True)
        self.g2=F.gauss([0.0,self.v2[self.cx,self.cy],self.a2[self.cx,self.cy],self.s2[self.cx,self.cy]],x=self.zv, returnmodel=True)
        self.g3=F.gauss([0.0,self.v3[self.cx,self.cy],self.a3[self.cx,self.cy],self.s3[self.cx,self.cy]],x=self.zv, returnmodel=True)
        self.g=self.g1+self.g2+self.g3
        self.g1-=self.g.min()
        self.g2-=self.g.min()
        self.g3-=self.g.min()
        self.g-=self.g.min()


    def applymask(self):
        self.a1.mask=self.mask1
        self.v1.mask=self.mask1
        self.s1.mask=self.mask1
        self.a2.mask=self.mask2
        self.v2.mask=self.mask2
        self.s2.mask=self.mask2
        self.a3.mask=self.mask3
        self.v3.mask=self.mask3
        self.s3.mask=self.mask3

    def updateplots(self):
        crx,cry=[self.cx-0.5,self.cx+0.5],[self.cy-0.5,self.cy+0.5]
        csx,csy=[self.cx-0.5,self.cx+0.5],[self.cy+0.5,self.cy-0.5]
        self.cross1a.set_data(crx,cry)
        self.cross1b.set_data(csx,csy)
        self.cross2a.set_data(crx,cry)
        self.cross2b.set_data(csx,csy)
        self.cross3a.set_data(crx,cry)
        self.cross3b.set_data(csx,csy)
        self.cross4a.set_data(crx,cry)
        self.cross4b.set_data(csx,csy)
        self.cross5a.set_data(crx,cry)
        self.cross5b.set_data(csx,csy)
        self.cross6a.set_data(crx,cry)
        self.cross6b.set_data(csx,csy)
        self.vf1.set_data(N.transpose(self.v1))
        self.vf2.set_data(N.transpose(self.v2))
        self.amp1.set_data(N.transpose(self.a1))
        self.amp2.set_data(N.transpose(self.a2))
        self.sig1.set_data(N.transpose(self.s1))
        self.sig2.set_data(N.transpose(self.s2))
        meas=self.cube[self.cx,self.cy]
        self.meas.set_ydata(meas)
        self.gauss1.set_ydata(masked_array(self.g1))
        self.gauss2.set_ydata(masked_array(self.g2))
        self.gauss3.set_ydata(masked_array(self.g3))
        self.gauss.set_ydata(masked_array(self.g))
        self.ax7.set_ylim(meas.min()/1.15,meas.max()*1.05)


    def update(self):
        self.setregion()
        self.applymask()
        self.parinfo=self.make_parinfo(self.zv,self.cube[self.cx,self.cy,:]+1)
        self.newgausses()
        self.updateplots()
        self.canvas.draw()


############################################
class PaSub:
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

        self.PaNumb=PaNumb
        self.osn=self.odata/self.oerror
        self.step=0.1
        self.fact=1.0
        self.shift=0

        if path.exists(prefix+'.pick'): self.subtracted=load(prefix+'.pick')
        else: self.subtracted=N.zeros(data.shape,'Float32')

        self.currdata()
        self.file=open(prefix+'.dat','w')

        self.fig=P.figure(1)

        self.canvas = self.fig.canvas
        self.canvas.mpl_connect('key_press_event', self.key_press_callback)
        self.canvas.mpl_connect('button_press_event', self.button_press_callback)

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
            print "using SIII"
            self.osynt=A.createPaschenSul(self.data+1,velRange=self.velRange,guessV=self.guessV,plotfit=False)
        else:
            print "using PaNumb %d"%self.PaNumb
            self.osynt=A.createPaschen(self.data+1,type='single',velRange=self.velRange,guessV=self.guessV,PaNumb=self.PaNumb,plotfit=False)
        self.synt=self.osynt.copy()


    def measurePa(self):
        params=A.fitAllPaschen_old(self.data,self.error,velRange=self.velRange,guessV=self.guessV,plot=False,prin=False)
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
        close(self.fig)
        sys.exit()

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
        ax=P.axes([0.02,0.68,0.94,0.27])
        plotspec(self.shiftscaled()+1,style='-b')
        plotspec(self.data,style='-k')
        plotspec(self.data-self.shiftscaled(),style='-r',vminmax='sigbased',Z=vel2z(self.guessV),plotlines=True,region='cat')
        P.setp(ax,xticks=[], yticks=[])
        P.title('CaT and Pa 13, 14, 15, 16')

        # SIII
        #ax=P.axes([0.74,0.68,0.23,0.27])
        ax=P.axes([0.02,0.35,0.23,0.27])
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
        P.plot(self.shiftscaled())
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
        ax=P.axes([0.38,0.02,0.23,0.27])
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
        P.setp(ax,xticks=[], yticks=[])
        P.title('Pa 17')

        ## PaStren Ratio
        ax=P.axes([0.38,0.02,0.23,0.27])
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
        ax=P.axes([0.76,0.02,0.20,0.27])

        #text(0.1,0.9,'S/N: '+str(self.sn.mean()),transform = ax.transAxes)
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



###############################################
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



##########################################
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





#############################

class vf_inter:
	def __init__(self, data, pixels=7, avg_square=3):

		# Basic Values
		self.avg_square_size = avg_square
		self.pixels = pixels
		self.spectra_size = (0.12,0.08)
		self.spectra_dist = (0.01,0.01)
		self.spectra_coord = (0.015,0.46)
		self.spectra_break = 6
		self.cenx = 0
		self.ceny = 0

		self.trusted = []
		self.trustgraph = []
		self.windownr = 1

		self.vel_pars = None
		self.flux_pars = None
		self.spectral_pars = None
		self.residual = None
		self.residual_sum = None
		self.residual_vel = None
		self.field = None
		self.modelled_flux = None
		self.radii = None
		self.synth_spectrum = None



		# import and process data cube
		self.odata = data.copy()
		self.cdata = data.cliparoundcenter()
		self.cdata = self.cdata
		self.cdata.p = self.odata.p
		self.sdata = G.sum(self.cdata)

		# import or calculate velocity field
		self.vel_file = self.odata.p['objname'] + '_peakvf_' + `self.pixels` + 'pixels.pick'
		if os.path.exists(self.vel_file):
			self.ini_velf = G.load(self.vel_file)
			print 'imported old velocity field'
		else:
			print 'recalculating velocity field ...,'
			self.ini_velf = G.peakvel(self.cdata,self.pixels)
			G.dump(self.ini_velf, self.vel_file)
			print 'done!\n'

		# setup main figure
		self.fig=P.figure(1)
		canvas = self.fig.canvas
		canvas.mpl_connect('key_press_event', self.key_press_callback)
		canvas.mpl_connect('button_press_event', self.button_press_callback)
		self.canvas = canvas
		self.plot()


	def plot(self):
		self.fig.clf()

        # Total Flux
		self.axflux=P.axes([0.00,0.59,0.34,0.34])
		axfluxresmin, axfluxresmax = self.nicegraph(self.sdata)
		P.setp(self.axflux,xticks=[], yticks=[])
		P.imshow(self.sdata,interpolation='nearest',vmin=axfluxresmin,vmax=axfluxresmax)
		P.title('Total Flux')
		P.colorbar()

		# Velocity Field
		self.axvelf=P.axes([0.34,0.59,0.34,0.34])
		axvelresmin, axvelresmax = self.nicegraph(self.ini_velf)
		P.setp(self.axvelf,xticks=[], yticks=[])
		P.imshow(self.ini_velf, interpolation='nearest',vmin=axvelresmin, vmax=axvelresmax)
		P.title('Velocity Field')
		P.colorbar()

		# Local Spectrum
		self.axspec=P.axes([0.015,0.46,0.12,0.08],xticks=[], yticks=[])


	def key_press_callback(self,event):
        # key bindings
		if event.key == 'a': self.fitvellinear()
		elif event.key == 'n': self.del_all_trusted()
		elif event.key == 's': self.fitfluxexp(modeltype='expfit_func')
		elif event.key == 'x': self.fitfluxexp(modeltype='powerlawfit_func')
		elif event.key == 'p': self.printfit()
		elif event.key == 'd': self.modelspectra_from_trusted()
		elif event.key == 'c': self.modelgauss_from_trusted()
		elif event.key == 'v': self.get_modelled_spectrum()
		elif event.key == 'q': self.close_all()
		elif event.key == 'm': self.open_menu()
		elif event.key == 'u': self.dump()
		elif event.key == 'i': self.load()
		else: print "Unknown key pressed:", event.key, '@ ('+`event.x`+','+`event.y`+ ") ... doing nothing"


	def button_press_callback(self,event):
		# mouse bindings, depending on axes
		if event.inaxes in [self.axflux, self.axvelf]:
			(x,y) = self.coord(event)
			if event.button == 1:
				self.plotat((x,y))
			elif event.button == 3:
				self.add_to_trusted((x,y))
			elif event.button == 2:
				self.set_center((x,y))

		elif event.inaxes in self.trustgraph:
			if event.button == 3:
				self.drop_trusted(event.inaxes)

		else: print '.'


	def res_key_press_callback(self, event):
		# key bindings for result window
		if   event.key == 'a': self.res_log_flux()
		elif event.key == 'y': self.res_log_flux(log=0)
		elif event.key == 's': self.res_log_vel()
		elif event.key == 'x': self.res_log_vel(log=0)
		else: print 'unused key:', event.key


	def res_button_press_callback(self, event):
		print 'unused button:', event.button


	def open_menu(self):
		# opens menu window to easily access models, modify graph range, etc.
		app = wx.PySimpleApp()
		name = self.odata.p['objname'] + ' Menu Window'
		frame = MenuWindow(None, -1, name)
		app.MainLoop()
		print 'huhu'
		###########################
		###########################
		###########################
		###########################


	def plotat(self, xy):
		# refresh local spectrum
		P.figure(1)
		self.axspec.cla()
		self.axspec.plot(self.cdata[xy[0],xy[1],:], 'r', linewidth=2)
		P.setp(self.axspec, xticks=[], yticks=[], title='Local')
		self.currdat = self.cdata[xy[0],xy[1],:]
		self.canvas.draw()


	def add_to_trusted(self, xy):
		# add coordinates, spectrum and velocity to trusted data
		P.figure(1)
		temp = N.zeros((self.cdata.shape[2],))
		cnt = 0.0
		dev = self.avg_square_size / 2
		for i in range(self.avg_square_size):
			for j in range(self.avg_square_size):
				di = xy[0]+i-dev
				dj = xy[1]+j-dev
				if di < self.cdata.shape[0] and dj < self.cdata.shape[1] and di >= 0 and dj >= 0:
					tmp = self.cdata[di,dj,:]
					temp += tmp
					cnt += 1.0
				else:
				    print 'Point at the border chosen! Continuing anyways.'
		temp = temp/cnt
#		veli = 'na'#G.calcpeak(temp, 7)
		newlistitem = [temp, xy, self.ini_velf[xy[0],xy[1]], self.sdata[xy[0],xy[1]]]
		self.trusted.append(newlistitem)
		print 'added point @', xy, 'to trusted points'
		self.update_trusted()


	def drop_trusted(self, axis):
		# remove trusted data
		P.figure(1)
		i = self.trustgraph.index(axis)
		print 'point @', self.trusted[i][1], 'being removed'
		self.trusted = self.trusted[:i] + self.trusted[i+1:]
		self.update_trusted(delf=i)


	def del_all_trusted(self):
		# remove all trusted data-points
		P.figure(1)
		for i in range(len(self.trusted)):
			self.trusted = self.trusted[:-1]
			self.update_trusted(delf=i)
			self.canvas.draw()
		print 'all trusted points removed'


	def update_trusted(self, delf=-1):
		# update graphs of trusted data
		P.figure(1)
		trustlen = len(self.trusted)

		# add new graphs
		if len(self.trustgraph) == trustlen - 1:
			self.trustgraph.append(P.axes(self.newaxis(trustlen-1, cutat = self.spectra_break), xticks=[],yticks=[]))
			P.plot(self.trusted[trustlen-1][0])

		# delete old ones
		elif len(self.trustgraph) == trustlen + 1 and delf >= 0:
			P.delaxes(self.trustgraph[-1])
			self.trustgraph = self.trustgraph[:-1]
			for i in range(delf,len(self.trusted)):
				self.trustgraph[i].cla()
				P.setp(self.trustgraph[i], xticks=[],yticks=[])
				self.trustgraph[i].plot(self.trusted[i][0])

		# notify errors
		elif len(self.trustgraph) == trustlen:
			print 'plotter idling ...'
		else:
			print 'something strange this way comes ...'

		# updating stuff
		self.canvas.draw()



	def fitvellinear(self):
		# perform a linear fit to trusted points, giving modelled VF as output

		P.figure(1)
		try: P.delaxes(self.axlinf)
		except: pass
		else: print 'old axis killed'

		print 'velocity fit started ...'
		# get and process input data
		ini_values = []
		ini_coordx = []
		ini_coordy = []
		ini_velocs = []
		self.data_shape = self.cdata.shape
		for i in self.trusted:
			ini_values.append(i[0])
			ini_coordx.append(i[1][0]-self.cdata.shape[0]/2-self.cenx)
			ini_coordy.append(i[1][1]-self.cdata.shape[1]/2-self.ceny)#-self.data_shape[1]/2-self.ceny)
			ini_velocs.append(i[2])
		ini_values = N.array(ini_values)
		ini_coordx = N.array(ini_coordx)
		ini_coordy = N.array(ini_coordy)
		ini_velocs = N.array(ini_velocs)
		print 'got values\n'

		# initialize mpfit fitting
		parinfo=[]
		parinfo.append({'value':self.odata.p['pa'], 'fixed':0, 'limited':[0,0],'limits':[0.0, 0.0], 'step':0.0, 'parname':'pa'})		# position angle
		parinfo.append({'value':1.0, 'fixed':0, 'limited':[0,0],'limits':[0.0, 0.0], 'step':0.1, 'parname':'vel gradient'})		# gradient
		parinfo.append({'value':self.odata.p['vr0'], 'fixed':0, 'limited':[0,0],'limits':[0.0, 0.0], 'step':1.0, 'parname':'vel gradient'})		# system

		err = 0
		functkw = {'x':ini_coordx, 'y':ini_coordy, 'z':ini_velocs, 'err':err}

		linfitres=mpfit.mpfit(linfit_func,parinfo=parinfo,functkw=functkw,quiet=False)

		# show results & store them
		print '\n\njuchu!'
		print linfitres.params
		print 'status', linfitres.status, '\n'
		print 'pa', (360-(linfitres.params[0]%360))
		self.vel_pars = {'model':'lin', 'posang':linfitres.params[0], 'gradient':linfitres.params[1], 'system':linfitres.params[2]}

		# checking the results
		seen=linfit_func(linfitres.params, None,  x=ini_coordx, y=ini_coordy, z=None, err=None, returnmodel=True)
		print '\ntrue velocities:\n', ini_velocs
		print '\nmodel velocities:\n', seen

		# visual output
		all_coordx = N.array(range(-self.cdata.shape[0]/2-self.cenx,self.cdata.shape[0]/2-self.cenx)*self.cdata.shape[1])
		all_coordx = all_coordx[::-1]
		all_coordy = N.array(range(-self.cdata.shape[1]/2-self.ceny,self.cdata.shape[1]/2-self.ceny)*self.cdata.shape[0])
		all_coordy.sort()
		self.field = linfit_func(linfitres.params, None,  x=all_coordx, y=all_coordy, z=None, err=None, returnmodel=True)
		self.field.shape = ((self.data_shape[0],self.data_shape[1]))
		self.field = N.transpose(self.field)
		self.axlinf=P.axes([0.68,0.02,0.32,0.32],xticks=[], yticks=[], title='Velocity Model')
		P.imshow(self.field, interpolation='nearest')

		self.canvas.draw()
		print 'done!\n'


	def fitfluxexp(self, modeltype, contnr=3):
		# fit exponential law to flux of trusted points

		P.figure(1)
		# throw away your television
		try: P.delaxes(self.axfluxmod)
		except: pass
		else: print 'old axis killed'

		print 'doing: ', modeltype, '...'
		if modeltype == 'expfit_func':
			modeltype = expfit_func
			modelt = 'exp'
		elif modeltype == 'powerlawfit_func':
			modeltype = powerlawfit_func
			modelt = 'powerlaw'

		# get & process values
		fluxes = N.array([])
		radii = []
		ini_coordx = []
		ini_coordy = []

		for i in self.trusted:
			ini_coordx.append(i[1][0]-self.cdata.shape[0]/2-self.cenx)
			ini_coordy.append(i[1][1]-self.cdata.shape[1]/2-self.ceny)
			tmp = i[0].copy()
			tmp.sort()
			tmp = tmp[:contnr]
			cont = tmp.sum()
			cont /= contnr
			flux = i[3] - self.cdata.shape[2] * cont
			fluxes = N.concatenate((fluxes, N.array([flux])))
			x = (i[1][0]-self.cdata.shape[0]/2-self.cenx)
			y = (i[1][1]-self.cdata.shape[1]/2-self.ceny)
			r = M.sqrt(x**2+y**2)
			radii.append(r)
		self.radii = N.array(radii)
		fluxes = N.log(fluxes)


		# initialize mpfit
		ini_coordx = N.array(ini_coordx)
		ini_coordy = N.array(ini_coordy)
		parinfo=[]
		parinfo.append({'value':10.0, 'fixed':0, 'limited':[1,0],'limits':[0.0, 0.0], 'step':0.0, 'parname':'amp'})		# amp
		parinfo.append({'value':-0.1, 'fixed':0, 'limited':[0,0],'limits':[0.0, 0.0], 'step':0.0, 'parname':'slope'})	# slope
		parinfo.append({'value':10.0, 'fixed':0, 'limited':[0,0],'limits':[0.0, 0.0], 'step':0.0, 'parname':'cutoff'})	# cutoff

		err = N.sqrt(fluxes)
		functkw = {'x':ini_coordx, 'y':ini_coordy, 'z':fluxes, 'err':err}

		expfitres=mpfit.mpfit(modeltype,parinfo=parinfo,functkw=functkw,quiet=False)

		# store data
		print expfitres.params
		print 'radii:',self.radii
		self.flux_pars = {'model':modelt, 'amp':expfitres.params[0], 'slope':expfitres.params[1], 'cutoff':expfitres.params[2]}

		# nice output
		self.axfluxmod=P.axes([0.72,0.75,0.27,0.20], yticks=[], title='Flux Model')
		self.axfluxmod.plot(self.radii, M.e**fluxes,'bo')
		self.axfluxmod.semilogy()

		# get fitted field
		maxrad = max(self.radii)
		modelled_flux = N.array([])
		l = N.array(range(1,int(round(maxrad))+1))
		self.modelled_flux = M.e **(modeltype(expfitres.params, None,  x=l, y=N.zeros((len(l),)), z=None, err=None, returnmodel=True))
		self.axfluxmod.plot(l, self.modelled_flux, 'k')

		self.canvas.draw()
		print 'done!\n'


	def modelspectra_from_trusted(self):
		# get typical spectral shape from trusted points

		P.figure(1)
		# die,die,die!
		try: P.delaxes(self.axspecmod)
		except: pass
		else: print 'old axis killed'

		# construct model spectrum
		try:
			vel_base = self.vel_pars['system']
		except TypeError:
			vel_base = self.odata.vel1st()
			print 'taking velocity from par file'
		else:
			print 'got fitted systemic velocity'

		vel_step = self.cdata.fsr()
		channr = self.cdata.shape[2]
		midchan = channr/2
		sspectra = N.zeros(channr)
		sspectra = N.array(sspectra, dtype='float32')
		for i in self.trusted:
			vel = i[2]
			chan = int(round((i[2]-vel_base)/vel_step))
			sspec = G.PyCigale.shift(i[0], chan)
			sspectra += sspec
		norm = N.sum(sspectra)
		sspectra /= norm

		# store fit in plot
		self.spectral_pars = {'model':'empirical', 'spectrum':sspectra}

		# plot it
		self.axspecmod=P.axes([0.70,0.47,0.29,0.20], xticks=[], yticks=[], title='Model Spectrum')
		self.axspecmod.plot(self.spectral_pars['spectrum'], 'r', linewidth=2)
		self.canvas.draw()
		print 'model spectrum done\n'


	def modelgauss_from_trusted(self):
		# assume gaussian shape and perform fit in width, amplitude and noise
		amplis = []
		widths = []
		widerr = []
		centrs = []
		noises = []
		coords = []

		parinfo=[]
		parinfo.append({'value':1.0, 'fixed':0, 'limited':[1,0],'limits':[0.0, 0.0], 'step':0.0, 'parname':'amp'})		# amplitude
		parinfo.append({'value':1.0, 'fixed':0, 'limited':[1,1],'limits':[0.0, self.cdata.shape[2]], 'step':0.0, 'parname':'width'})	# width
		parinfo.append({'value':self.cdata.shape[2]/2, 'fixed':0, 'limited':[1,1],'limits':[0.0, self.cdata.shape[2]], 'step':0.0, 'parname':'center'})	# center
		parinfo.append({'value':1.0, 'fixed':0, 'limited':[1,0],'limits':[0.0, 0.0], 'step':0.0, 'parname':'noise'})	# noise
	    # fit each and every spectrum
		for i in self.trusted:
			y = N.array(i[0])
			x = N.arange(len(y))
			err = N.sqrt(y)
			functkw = {'x':x, 'y':y, 'err':err}

			gaussfitres=mpfit.mpfit(gauss_fit_func,parinfo=parinfo,functkw=functkw,quiet=True)

			pp = gaussfitres.params
			amplis.append(pp[0])
			widths.append(pp[1])
			widerr.append(gaussfitres.perror[1])
			centrs.append(pp[2])
			noises.append(pp[3])

			dx = i[1][0] - self.cdata.shape[0] - self.cenx
			dy = i[1][1] - self.cdata.shape[1] - self.ceny
			dr = M.sqrt(dx**2+dy**2)
			coords.append(dr)

			#print '\n\n'

		# convert to arrays
		amplis = N.array(amplis)
		widths = N.array(widths)
		widerr = N.array(widerr)
		centrs = N.array(centrs)
		noises = N.array(noises)
		coords = N.array(coords)

		# show all results
		print '\n'
		print amplis,'\n'
		print widths,'\n'
		print centrs,'\n'
		print noises,'\n'
		P.figure()
		P.plot(coords, widths, 'bo')
		#P.show()
		self.windownr += 1

		# linear fit to width
		parinfo=[]
		parinfo.append({'value':1.0, 'fixed':0, 'limited':[1,0],'limits':[0.0, 0.0], 'step':0.0, 'parname':'offset'})		# initial
		parinfo.append({'value':1.0, 'fixed':0, 'limited':[0,0],'limits':[0.0, 0.0], 'step':0.0, 'parname':'slope'})	# slope
		functkw = {'x':coords, 'y':widths, 'err':widerr}
		widthfit=mpfit.mpfit(lin_fit_func,parinfo=parinfo,functkw=functkw,quiet=False)

		# show fit
		offset = widthfit.params[0]
		slope =  widthfit.params[1]
		xi = N.arange(min(coords),max(coords))
		yi = offset + xi * slope
		P.plot(xi,yi,'k')

		####################
		####################
		####################
		####################


	def get_modelled_spectrum(self):
		# use all of the results from fitting procedures to form model spectrum

		# check for performed fits, else abort
		if self.vel_pars == None:
			print '!!! No velocity model yet. Please perform modelling first (a) !!!'
			return
		elif self.flux_pars == None:
			print '!!! No flux model yet. Please perform modelling first (s,x) !!!'
			return
		elif self.spectral_pars == None:
			print '!!! No spectral model yet. Please perform modelling first (d) !!!'
			return

		# get velocities from vel model, rescale with flux model, give obtained spectral shape
		print 'alright, starting modelling with current fits ...'

		# get vel model
		if self.vel_pars['model'] == 'lin':
			vel_model = linfit_func
			vel_pars = [self.vel_pars['posang'], self.vel_pars['gradient'], self.vel_pars['system']]

		# get flux model
		if self.flux_pars['model'] == 'exp':
			flux_model = expfit_func
			flux_pars = [self.flux_pars['amp'], self.flux_pars['slope']]
		elif self.flux_pars['model'] == 'powerlaw':
			flux_model = powerlawfit_func
			flux_pars = [self.flux_pars['amp'], self.flux_pars['slope'], self.flux_pars['cutoff']]

		# get spectral model
		if self.spectral_pars['model'] == 'empirical':
			spectral_model = self.spectral_pars['spectrum']

		# get spectral model
		start = N.array(spectral_model)
		start.shape = (len(spectral_model),1,1)
		inter = N.repeat(start, self.cdata.shape[0], axis=1)
		model = N.repeat(inter, self.cdata.shape[1], axis=2)
		model = N.transpose(model, (1,2,0))

		# build coordinate arrays
		wdt = self.cdata.shape[0]/2
		hgt = self.cdata.shape[1]/2
		all_coordx = N.array(range(-wdt,wdt)*hgt*2)
		all_coordx = all_coordx[::-1]
		all_coordy = all_coordx.copy()
		all_coordy.sort()

		# get velocity model
		velmod = vel_model(vel_pars, None,  x=all_coordx, y=all_coordy, z=None, err=None, returnmodel=True)
		velmod.shape = (2*wdt,2*hgt)
		velmod = N.transpose(velmod)

		# get flux model
		fluxmod = N.exp(flux_model(flux_pars, None,  x=all_coordx, y=all_coordy, z=None, err=None, returnmodel=True))
		fluxmod.shape = (2*wdt,2*hgt)
		fluxmod = N.transpose(fluxmod)
		self.fluxmod = fluxmod

		# get unshifted spectral model
		intermediate_model = fluxmod[:,:,N.newaxis] * model

		# get indexshifts
		channr = self.cdata.shape[2]
		midchan = channr/2
		vel_base = self.vel_pars['system']
		vel_step = self.cdata.fsr()/channr
		shiftarray = N.around((velmod-vel_base)/vel_step)
		shiftarray.dtype = 'Int32'

		# final model
		fin_mod = intermediate_model.copy()
		for i in range(fin_mod.shape[0]):
			for j in range(fin_mod.shape[1]):
				fin_mod[i,j,:] = G.PyCigale.shift(fin_mod[i,j,:], int(shiftarray[i,j]))
		print 'done!\n'

		print 'creating residuals ...'
		self.synth_spectrum = G.PyCigale.array(fin_mod)
		self.synth_spectrum.p = self.cdata.p
		self.synth_spectrum.p['vr0'] = self.vel_pars['system']
		self.residual = self.cdata - self.synth_spectrum
		self.residual_sum = G.sum(self.residual)
		self.residual_vel = G.peakvel(self.residual, 7)

		# show it!!
		self.result = P.figure()
		self.rescanvas = self.result.canvas
		P.show()

		# Total Flux
		axfluxresmin, axfluxresmax = self.nicegraph(self.residual_sum)
		self.axfluxres=P.axes([0.00,0.59,0.36,0.36], xticks=[], yticks=[],title='Total Flux (Residual)')
		self.axfluxres.imshow(self.residual_sum, interpolation='nearest', vmin=axfluxresmin,vmax=axfluxresmax)
		self.canvas.draw()
		P.colorbar()

		# Velocity Field
		axvelresmin, axvelresmax = self.nicegraph(self.residual_vel)
		self.axvelres=P.axes([0.36,0.59,0.36,0.36],xticks=[], yticks=[],title='Velocity Field (Residual)')
		self.axvelres.imshow(self.residual_vel, interpolation='nearest', vmin=axvelresmin,vmax=axvelresmax)
		self.canvas.draw()
		P.colorbar()

#		P.show()
		self.windownr += 1
		self.rescanvas.mpl_connect('key_press_event', self.res_key_press_callback)
		self.rescanvas.mpl_connect('button_press_event', self.res_button_press_callback)
		print 'done!\n'


	def res_log_flux(self,log=1):
		# log flux plot
		if log:
			temp = N.log(N.where(self.residual_sum>0.001,self.residual_sum,0.001))
		else:
			temp = self.residual_sum
		self.axfluxres.cla()

		minmaxflux = N.array(temp)
		minmaxflux.shape = (-1,)
		minmaxflux = N.sort(minmaxflux)
		takeflux = int(round(len(minmaxflux)*0.01))
		axfluxresmax = max(int(round(minmaxflux[-takeflux])),1)
		axfluxresmin = max(int(round(minmaxflux[takeflux])), 0)
		if log:
			P.setp(self.axfluxres, xticks=[], yticks=[], title='Flux (residual,log)')
		else:
			P.setp(self.axfluxres, xticks=[], yticks=[], title='Flux (residual)')
		self.axfluxres.imshow(temp, interpolation='nearest', vmin=axfluxresmin,vmax=axfluxresmax)
		self.rescanvas.draw()


	def res_log_vel(self,log=1):
		# log vel plot
		if log:
			temp = N.log(N.where(self.residual_vel>0.001,self.residual_vel,0.001))
		else:
			temp = self.residual_vel
		self.axvelres.cla()
		print temp

		minmaxflux = N.array(temp)
		minmaxflux.shape = (-1,)
		minmaxflux = N.sort(minmaxflux)
		takeflux = int(round(len(minmaxflux)*0.01))
		axfluxresmax = max(minmaxflux[-takeflux],1)
		axfluxresmin = max(minmaxflux[takeflux],0)
		print axfluxresmax, axfluxresmin
		if log:
			P.setp(self.axvelres, xticks=[], yticks=[], title='Flux (residual,log)')
		else:
			P.setp(self.axvelres, xticks=[], yticks=[], title='Flux (residual)')
		self.axvelres.imshow(temp, interpolation='nearest', vmin=axfluxresmin,vmax=axfluxresmax)
		self.rescanvas.draw()


	def dump(self):
		# dump all relevant data
		dumper = {}
		dumper['vel_pars'] = self.vel_pars
		dumper['vel_field'] = self.field
		dumper['flux_pars'] = self.flux_pars
		dumper['modelled_flux'] = self.modelled_flux
		dumper['radii'] = self.radii
		dumper['spectral_pars'] = self.spectral_pars
		dumper['synth_spectrum'] = self.synth_spectrum
		dumper['residual'] = self.residual
		dumper['residual_sum'] = self.residual_sum
		dumper['residual_vel'] = self.residual_vel
		dumper['trusted'] = self.trusted

		dumpname = raw_input('Filename for dump? ')#self.odata.p['objname'] + '_interactvf_dumpfile.pick'
		G.dump(dumper, dumpname)
		print 'all dumped'


	def load(self):
		# load the dumped data
		#dumpname = self.odata.p['objname'] + '_interactvf_dumpfile.pick'
		dumpname=raw_input('Filename for load? ')
		if os.path.exists(dumpname):
			dumper = G.load(dumpname)
			print 'dump imported\n'

			print 'rebuilding data structure ...'
			self.vel_pars = dumper['vel_pars']
			self.field = dumper['vel_field']
			self.flux_pars = dumper['flux_pars']
			self.modelled_flux = dumper['modelled_flux']
			self.radii = dumper['radii']
			self.spectral_pars = dumper['spectral_pars']
			self.synth_spectrum = dumper['synth_spectrum']
			self.residual = dumper['residual']
			self.residual_sum = dumper['residual_sum']
			self.residual_vel = dumper['residual_vel']
			self.trusted = dumper['trusted']
			print 'done!\n'


		else:
			print 'no dumpfile found'


	def set_center(self,xy):
		# select new center for the galaxy
		self.cenx = xy[0] - self.cdata.shape[0]/2
		self.ceny = xy[1] - self.cdata.shape[1]/2
		print 'new center selected @', xy


	def printfit(self):
		# print parameters from fitting
		print '\nParameters for velocity field model:\n', self.vel_pars
		print '\nParameters for radial flux model:\n', self.flux_pars
		print '\nParameters for spectral model:\n', self.spectral_pars
		print '\n'


	def nicegraph(self, values, cut=0.02):
		# asjust upper & lower limit for graphs
		minmax = N.array(values)
		minmax.shape = (-1,)
		minmax = N.sort(minmax)
		take = int(round(len(minmax)*cut))
		tmax = max(int(round(minmax[-take])),1)
		tmin = max(int(round(minmax[take])), 0)
		return (tmin, tmax)


	def close_all(self):
		# close all opened windows
		for i in range(self.windownr):
			try: P.close(self.windownr-i)
			except: pass
		print 'interactvf closed. goodbye!\n\n'


	def coord(self, event):
		# extract proper coordinates of mouse event
		y = int(round(event.xdata))
		x = int(round(event.ydata))
		self.xy = (x,y)
		return self.xy


	def newaxis(self, i, cutat=4):
		# iterate window positions for trusted data
		row = (i+1)/cutat
		col = (i+1)%cutat
		coordadd = (row*(self.spectra_size[0] + self.spectra_dist[0]), -col*(self.spectra_size[1] + self.spectra_dist[1]))
		coords = (self.spectra_coord[0] + coordadd[0], self.spectra_coord[1] + coordadd[1])
		return [coords[0], coords[1], self.spectra_size[0], self.spectra_size[1]]


def radii_from_position(coords, pars):
	# calc radii from x,y coordinates
	pa = pars['pa']*M.pi/180.0
	x = coords[0]
	y = coords[1]
	cenx = pars['centr_offset_x']
	ceny = pars['centr_offset_y']
	radii = []

	for i in range(len(x)):
		sx = (x[i] - cenx)
		sy = (y[i] - ceny)
		if sy == 0:
			if sx >= 0: ang = M.pi/2
			else: ang = -M.pi/2
		else: ang = M.atan(sx/sy)
		t = sx * M.sin(-pa) + sy * M.cos(pa)
		o = sx * M.cos(pa)  + sy * M.sin(pa)
		dx = cenx + t * M.sin(-pa) - x[i]
		dy = cenx + t * M.cos(pa)  - y[i]
		if o >= 0: sgn = +1
		else: sgn = -1
		r = sgn*M.sqrt(dx**2 + dy**2)
		radii.append(r)

	radii = N.array(radii)
	return radii


#######################################
# GENERIC MATPLOTLIB EVENT HANDLERS
def on_click(event):
  """Handle click events in plotting windows. Not to be called from
      the command line. It saves the integer coordinates for the
      clicks with button 1 in a temporary file. Button 3 closes the
      figure.
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
  """Handle click events in plotting windows. Not to be called from
      the command line. It saves the float coordinates for the
      clicks with button 1 in a temporary file. Button 3 closes
      the figure.
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
