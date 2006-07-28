"""
 PyCigale
 
 Thomas Marquart July 2004
 
 This is the major file for my Velocity Field and ADHOC package.
 
 Set the PYTHONPATH environment varable to the directory
 that contains the directory containing this file or put it into
 your equivalent to /usr/local/lib/python2.3/site-packages/

 Start a python shell.
 Import this package with "import PyCigale.PyCigale as PyC"
 Now all functions and classes defined by this package are available
 via PyC.foo
   
 This file defines two classes:
  * ADP:
       * a class for contianing all important parameters, for a list see a .adp
         file or look into InOutput.py
       * it can write and read .adp files with its methods
       * it has some more attributes apart from the ones in .adp files

  * adhoc:
       * can hold 2D and 3D data
       * understands and can write the AD3 and AD2 formats with headers
       * uses matplotlib for visualisation
       * it is derived from "numarray", a very powerful and fast numeric library
         therefore the arrays are very comfortable and have loads of built-in
         methods, accessible via mydata.foo
       * it has reimplemented some tasks that ADHOC in principle also can do for
         convenience and speed: masking, summing, cutting images etc. more will come
       * it is also derived from the ADP class. this makes (for example) the
         object scan wavelength accessible via mydata.xlbneb
       * to read a ad3 or ad2 file use data=PyC.fromAD('filename')
 
 In addition to the classes, this file contains various global functions that are useful
 to apply to adhoc-type instances.

 The functions containing functions for Input and Output (to files) are sourced out to the
 the file InOutput.py in the same directory as this file

 The same is true for the functions that generate .adm macros (MacroGen.py) and the ones that
 do (not yet very advanced) model velocity fields (ModelVF.py)

 Additional documentation can be found in the doc string (accessible via PyC.foo.__doc__) of each
 function or method.
 
"""


# IMPORTS
#
import numarray as N
import numarray.numarraycore as _n
import os
import copy
import commands
import matplotlib
import matplotlib.pylab as MP
#matplotlib.use('GTK')
import Numeric
import mpfit
from math import pi,e,radians
#
# END IMPORTS


# THE OTHER FILES
#
from MacroGen import *
from InOutput import *
from ModelVF import *
from Gauss import *
#
# END THE OTHER FILES

# THE GLOBAL PARAMETER DICTIONARY
#import PyGalKin
#__PAR__=PyGalKin.__PAR__



# CONSTANTS
#
# units:
#        velocity: km/s
#        wavelength: Angstrom
lambHA=N.array([6562.7797852000003],'Float32')
sol=N.array([299792.458],'Float32')
H0=N.array([72.],'Float32')
G=N.array([6.6726E-11*1.989E30/1000**3],'Float32') ### in solar masses and km
pc=N.array([3.086E13],'Float32') ## in km
#
# END CONSTANTS


# SHORTCUTS
#
tab='\t'
nl='\n'
bs='\\'
sp=' '
null='\0'

font = {'fontname'   : 'Courier',
        'color'      : 'k',
        #'fontweight' : 'bold',
        'fontsize'   : 20}
#
# END SHORTCUTS


# CLASS: ADHOC
#
class adhoc(_n.NumArray, object):
    """The class for CIGALE datacubes and images
       see the head of the file PyCigale.py for a general description
    """
    
    #    def __init__(self, n, a):
    #        ''' n provides the length of each dimension,
    #        a is the constant value to be plugged.
    #        '''
    #        arr= N.numarraycore.array(sequence= n * n * [a], shape= (n, n))
    #        self.__setstate__(arr.__getstate__())
    
    #def __init__(self,filename):

    #    data,p=fromAD(filename)
    #    self=data
    #    self.p=p


    def __repr__(self):
        " Return printable representation of instance."
        className= self.__class__.__name__
        className= className.zfill(5).replace('0', ' ')
        arr= self.copy()
        arr.__class__= _n.NumArray
        rep= className + _n.NumArray.__repr__(arr)[5:]
        return rep

    def __str__(self):
        " Return a pretty printed string of the instance."
        stri= self.copy()
        stri.__class__= _n.NumArray
        return _n.NumArray.__str__(stri)
 
    # Wrappers for numarray-methods to conserve the dictionary p
    #
    def __add__(self, a):
        tmp = _n.NumArray.__add__(self, a)
        if hasattr(self,'p'): tmp.p = self.p.copy()
        return tmp
    
    def __sub__(self, a):
        tmp = _n.NumArray.__sub__(self, a)
        if hasattr(self,'p'): tmp.p = self.p.copy()
        return tmp
    
    def __mul__(self, a):
        tmp = _n.NumArray.__mul__(self, a)
        if hasattr(self,'p'): tmp.p = self.p.copy()
        return tmp
    
    def __div__(self, a):
        tmp = _n.NumArray.__div__(self, a)
        if hasattr(self,'p'): tmp.p = self.p.copy()
        return tmp

    def __radd__(self, a):
        tmp = _n.NumArray.__radd__(self, a)
        if hasattr(a,'p'): tmp.p = a.p.copy()
        elif hasattr(self,'p'): tmp.p = self.p.copy()
        return tmp
    
    def __rsub__(self, a):
        tmp = _n.NumArray.__rsub__(self, a)
        if hasattr(a,'p'): tmp.p = a.p.copy()
        elif hasattr(self,'p'): tmp.p = self.p.copy()
        return tmp
    
    def __rmul__(self, a):
        tmp = _n.NumArray.__rmul__(self, a)
        if hasattr(a,'p'): tmp.p = a.p.copy()
        elif hasattr(self,'p'): tmp.p = self.p.copy()
        return tmp
    
    def __rdiv__(self, a):
        tmp = _n.NumArray.__rdiv__(self, a)
        if hasattr(a,'p'): tmp.p = a.p.copy()
        elif hasattr(self,'p'): tmp.p = self.p.copy()
        return tmp

    def __pow__(self, a):
        tmp = _n.NumArray.__pow__(self, a)
        if hasattr(self,'p'): tmp.p = self.p.copy()
        return tmp
    
    def __getitem__(self,i):
       tmp=_n.NumArray.__getitem__(self,i)
       if isinstance(tmp,float): tmp=array([tmp])
       if hasattr(self,'p'): tmp.p = self.p.copy()
       return tmp
    
    def sum(self,axis=2):
        tmp=N.sum(self,axis)
        if hasattr(self,'p'): tmp.p = self.p.copy()
        return tmp

     #END: Wrappers for numarray-methods to conserve the dictionary p
    
    
    def cenx(self):
      return self.p['cen'][0]
    
    def ceny(self):
      return self.p['cen'][1]
    
    def fsr(self):
      """ calculate the free spectral range from self.xil and self.xlp
          uses lamb2vel
      """
      return lamb2vel(self.p['xlbneb'] + (self.p['xil'] / 2)) - lamb2vel(self.p['xlbneb'] - (self.p['xil'] / 2))
    
    def vel1st(self):
        """return the velocity of the first channel """
        #return lamb2vel(self.p['xlbneb']) - (self.fsr() / 2.)
        return lamb2vel(self.p['xl1'])
        
    def get_copy(self):
      """ returns a copy of the instance """
      temp = copy.copy(self)
      temp.p = copy.deepcopy(self.p)
      return temp
    
    def toADP(self,filename=None):
      """ write parameters to an ADP file
          if no filename is given it uses the one which was read before!
      """
      toADP(self,filename)
    
    def fromADP(self,filename):
      """ read parameters from an ADP file """
      fromADP(self,filename)
    
    def toPAR(self,filename='par'):
      """ write parameters to an ADP file
          if no filename is given it uses the one which was read before!
      """
      toPAR(self,filename)
    
    def fromPAR(self,filename='par'):
      """ read parameters from an ADP file """
      fromPAR(self,filename)
    
    def distance(self):
      """ returns the distance (Mpc) to the object using Hubble's law"""
      return self.p['vr0']/H0
    
    def scale(self):
      """returns the physical scale for the object (pc/pix) """
      return self.p['echelle']/3600/360*2*pi * self.distance()*1E6
    
    def M(self,m=None):
      """ returns absolute magnitute"""
      if (m==None): 
        m=self.p['mB']
      return app2abs(self.distance(),m)
    
    def fixcenterofrings(self):
      """ write a first guess for the center of rings for phase computation"""
      
      if (self.p['p'] == 793):
        self.p['xc']=302.5
        self.p['yc']=248.0
      elif (self.p['p'] == 1938):
        self.p['xc']=302.5
        self.p['yc']=248.0
      else: 
        pass
      
      self.toADP()
    
    def getsaltzer(self,filename='/home/tom/data/CIGALE-2004/saltzer.dat'):
      """ reads the data from saltzer etal."""
      
      for line in open(filename):
        if (self.p['objname'] in line):
          self.p['mB']=float(line.split()[1])
          self.p['mBe']=float(line.split()[2])
          self.p['BV']=float(line.split()[3])
          self.p['BVe']=float(line.split()[4])
          self.p['WHb']=float(line.split()[5])
    
    def toAD(self,filename):
      """write to an adhoc format file
      recognises it its a 2d or 3d file
      """
      toAD(self,filename)
    
    def secondmoment(self):
      """wrapper for function secondmoment"""
      return secondmoment(self)
    
    def firstmoment(self):
      """wrapper for function firstmoment"""
      return firstmoment(self)
    
    def cliparoundcenter(self, relsize=None):
      """clips spatially the data around .cen with the size
      given and returns the clipped array (no inplace clipping)
      """
      if (relsize==None): 
        relsize=self.p['relsize']
      
      outarr = self[self.cenx()-relsize:self.cenx()+relsize,self.ceny()-relsize:self.ceny()+relsize]
      outarr.p=self.p.copy()
      outarr.p['cen']=N.array([relsize,relsize],shape=(2,))
      outarr.p['dyncen']=outarr.p['cen']+(self.p['dyncen']-self.p['cen'])
      return outarr
    
    def maskwith(self,Mask,lower=None,upper=None,value=0):
      """ wrapper for mask-funktion
          works inplace!
      """
      mask(self,Mask,lower,upper,value)
    
    def show(self,upper=None,lower=None,title='',savefig='/tmp/show.png'):
      """ plotting a 2D image or a spectrum"""
      if (lower == None):
        lower=self.min()
      if (upper == None):
        upper=self.max()
      
      if (self.shape.__len__() == 3):
        print '3d plots not implemented'
      elif (self.shape.__len__() == 2):
        MP.figure(num = 1, figsize=(8, 8), dpi=80, facecolor='w', edgecolor='k', frameon=True) 
        MP.imshow(N.swapaxes(self,0,1),vmin=lower, vmax=upper, interpolation='nearest', origin='lower', aspect='preserve',extent=(0,self.nx()-1,0,self.ny()-1))
        MP.savefig(savefig)
        MP.show()
      elif (self.shape.__len__() == 1):
        MP.plot(self,'g-')
        MP.plot(self,'ro')
        MP.savefig(savefig)
        MP.show()
      else:
        pass
    
    def showmap(self,lower=None,upper=None,savefig=None):
      """plotting an image with kpc axes"""
      
      if ('Mono' in self.p['imagename']):
        type='Monochromatic H-alpha Flux'
        if (savefig==None): 
          savefig='/tmp/plots/'+self.p['objname']+'_Mono.png'
        if (lower == None):
          lower=self.p['monocuts'][0]
        if (upper == None):
          upper=self.p['monocuts'][1]
      elif ('Cont' in self.p['imagename']):
        type='Continuum Flux'
        if (savefig==None): 
          savefig='/tmp/plots/'+self.p['objname']+'_Cont.png'
        if (lower == None):
          lower=self.p['contcuts'][0]
        if (upper == None):
          upper=self.p['contcuts'][1]        
      else:
        type='unknown'
        savefig='/tmp/plots/unknown.png'
        
      MP.figure(num=1, figsize=(8.14, 8), dpi=80, facecolor='w', edgecolor='k')
      MP.imshow(N.swapaxes(self,0,1), vmin=lower, vmax=upper, interpolation='nearest', origin='lower', aspect='preserve')
      setXaxis_pc(self)
      setYaxis_pc(self)
      MP.title(self.p['objname'] + ' - ' + type,font)
      MP.axis([0,self.nx()-1,0,self.ny()-1])
      MP.savefig(savefig)
      MP.connect('button_press_event', on_click)
      MP.show()
    
    def showVF(self, lower=None, upper=None, plotaxis=True, slitname=None, gray=False, savefig=None):
      """Plots a velocity field with kpc axes.
          
          Usage: vf.showVF(lower,upper,plotaxis,slitname,gray,savefig)
          
          lower: The lowest value to plot
          upper: The highest value to plot
          plotaxis: True or False, plots the major axis defined by the pa and centrum
            default is True
          slitname: Plots the slit 'slitname.slit'
          gray: False or True, plot grayscale, default is False
          savefig: Save the plot to 'savefig'
          """
      if (lower == None):
        lower=self.p['velcuts'][0]
      if (upper == None):
        upper=self.p['velcuts'][1]
      
      if (savefig==None): 
        savefig='/tmp/plots/'+self.p['objname']+'_VF.png'
      
      MP.figure(num=1, figsize=(8.14, 8), dpi=80, facecolor='w', edgecolor='k')
      
      if (gray==True):
        MP.gray()
      else:
        MP.jet()
      
      MP.imshow(N.swapaxes(self,0,1), vmin=lower, vmax=upper, interpolation='nearest', origin='lower', aspect='preserve')
      MP.colorbar()
      
      # Plot the pa-axis
      if (plotaxis==True):
        vec=N.array([-N.sin(radians(self.p['pa'])),N.cos(radians(self.p['pa']))])
        p1=self.p['dyncen']+(vec*self.p['relsize']*1.5)
        p2=self.p['dyncen']-(vec*self.p['relsize']*1.5)
        MP.plot([p1[0],p2[0]],[p1[1],p2[1]],'w-',linewidth=3)
      
      # Plot the slit
      if (slitname!=None):
        pars = read_slit(slitname)
        angle=pars['angle']
        center=self.p['cen'] + pars['offset']
        d = pars['slitwidth']/self.p['echelle']
        
        if (pars['slit_angle_dep']==1): 
          D = d/2 + 1/2 + N.sin(2*radians(angle))*(N.sqrt(2)-1)/2
        else:
          D = d/2 + 1/2
        
        vec=N.array([-N.sin(radians(angle)),N.cos(radians(angle))])
        perp=N.array([N.cos(radians(angle)),N.sin(radians(angle))])
        
        p1 = center + D*perp + vec*self.p['relsize']*1.5
        p2 = center + D*perp - vec*self.p['relsize']*1.5
        p3 = center - D*perp + vec*self.p['relsize']*1.5
        p4 = center - D*perp - vec*self.p['relsize']*1.5
        
        MP.plot([p1[0],p2[0]],[p1[1],p2[1]] , 'w-', linewidth=2)
        MP.plot([p3[0],p4[0]],[p3[1],p4[1]] , 'w-', linewidth=2)
      
      # Set the axis-scales, title and plot the figure
      setXaxis_pc(self)
      setYaxis_pc(self)
      MP.title(self.p['objname'] + ' - ' + 'Radial Velocity',font)
      MP.axis([0,self.nx()-1,0,self.ny()-1])
      MP.savefig(savefig)
      MP.connect('button_press_event', on_click)
      MP.show()
    
    def PVdiag(self):
      """Calculate the position-velocity diagram from a velocity field.
          
          Usage: pos,vel = vf.PVdiag()
          
          vf: The velocity field to calculate from
          pos: A vector of positions
          vel: A vector of corresponding velocities
          """
      vec=N.array([-N.sin(radians(self.p['pa'])),N.cos(radians(self.p['pa']))])
      
      pos=N.zeros(0,'Float32')
      vel=N.zeros(0,'Float32')
      
      for i in N.arange(self.nx()):
        for j in N.arange(self.ny()):
          if self[i,j] != 0 :
            pos_temp=N.innerproduct(vec,N.array([self.p['dyncen'][0]-i,self.p['dyncen'][1]-j]))
            pos=N.concatenate((pos, pos_temp*self.scale()/1000))
            vel=N.concatenate((vel, self[i,j]))
      
      return pos,vel
    
    def showRC(self, pos_and_vel , upper=None, lower=None, savefig=None):
      """Plot a rotation curve or position velocity diagram.
          
          Usage: vf.showRC([pos,vel], upper, lower, savefig)
          
          pos: A vector of positions
          vel: A vector of corresponding velocities
          lower: The lowest value to plot
          upper: The highest value to plot
          savefig: Save the plot to 'savefig', default is None
          """
      if (savefig==None): 
        savefig='/tmp/plots/'+self.p['objname']+'_PV.png'
      
      pos=pos_and_vel[0]
      vel=pos_and_vel[1]
      
      if (lower == None):
        lower=self.p['velcuts'][0]
      if (upper == None):
        upper=self.p['velcuts'][1]
      
      MP.plot(pos,vel,'+')
      MP.title(self.p['objname'] + ' - ' + 'PV-diagram',font)
      MP.ylabel('Radial Velocity [km/s]',font)
      MP.xlabel('[kpc]',font)
      MP.axis([pos.min(),pos.max(),lower,upper])
      MP.savefig(savefig)
      MP.connect('button_press_event', on_click)
      MP.show()
    
    def RCslit(self, slitname='1', outfile=None):
      """Calculate the position-velocity diagram along a slit from a velocity field.
          
          Usage: pos,vel = vf.RCslit(slitname, outfile)
          
          vf: The velocity field to calculate from
          slitname: slitname.slit is used as .slit-file, default is '1'
          outfile: The result is written to outfile.rc, default is None
          pos: A vector of positions
          vel: A vector of corresponding velocities
          """
      # Read slit data
      pars=read_slit(slitname)
      angle=pars['angle']
      center=self.p['cen'] + pars['offset']
      d = pars['slitwidth']/self.p['echelle']
      
      if (pars['slit_angle_dep']==1): 
        D = d/2 + 1/2 + N.sin(2*radians(angle))*(N.sqrt(2)-1)/2
      else: 
        D = d/2 + 1/2
      
      # Vectors needed for calculations
      vec=N.array([-N.sin(radians(angle)),N.cos(radians(angle))])
      perp=N.array([N.cos(radians(angle)),N.sin(radians(angle))])
      
      # Arrays for the data
      pos=N.zeros(0,'Float32')
      vel=N.zeros(0,'Float32')
      err=N.zeros(0,'Float32')
      
      for i in N.arange(self.nx()):
        for j in N.arange(self.ny()):
          if (self[i,j] != 0):
            # The points must be within the slit
            dist=N.innerproduct(perp,N.array([center[0]-i,center[1]-j]))
            if (abs(dist) <= D):
              pos_temp=N.innerproduct(vec,N.array([center[0]-i,center[1]-j]))
              pos=N.concatenate((pos, pos_temp))
              vel=N.concatenate((vel, self[i,j]))
              err=N.concatenate((err, self.fsr()/self.p['lz']))
      
      # Write rc-outfile if a filename is given
      if (outfile != None):
        write_rc(outfile, pos, vel, err)
      
      return pos,vel
    
    def WidthMap(self):
      """ giving back a 2d map with the fwhm of the 3d cube"""
      return doforeachpoint(self,fwhm)/self.nz()*self.fsr()
    
    def ndim(self):
      """returns the dimension of the cube """
      return self.shape.__len__()
    
    def nx(self):
      """returns the size of the x-dimension """
      return self.shape[0]
    
    def ny(self):
      """returns the size of the y-dimension """
      if self.ndim() > 1:
        return self.shape[1]
      else:
        return 1
    
    def nz(self):
      """returns the size of the z-dimension """
      if self.ndim() > 2:
        return self.shape[2]
      else:
        return 1
#
# END: CLASS: ADHOC


# WRAPPERS
#
def fromfile(*args, **keys):
  a = N.fromfile(*args, **keys)
  a.__class__ = adhoc
  #a.p = {}
  return a

def array(*args, **keys):
  a = N.array(*args, **keys)
  a.__class__ = adhoc
  a.p = {}
  return a

def stripadhoc(arr):
  tmp=N.zeros(arr.shape)
  tmp[:]=arr[:]
  return tmp

def sum(arr,axis=2):
  tmp = N.sum(arr,axis)
  if hasattr(arr,'p'): tmp.p = arr.p.copy()
  return tmp
#
# END: WRAPPERS


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
    MP.close()

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
    MP.close()

  file.close()
#
# END: MATPLOTLIB EVENT HANDLERS


# TOOLS
#
def shift_spectra(data, i=0):
  """ Shift all spectra in an array and correctly adds the velocity difference
      to p['vr_offset'].
      
      Usage: new_arr = shift(arr, i)
      
      arr:  The array to be shifted
      i:  The steps to shift each vector with
      new_arr: The array shifted i steps
      """
  temp = doforeachpoint(data, shift, i)
  # Offset in km/s
  temp.p['vr_offset'] = temp.p['vr_offset'] - (i*data.fsr()/data.p['lz'])
  return temp

def VF_gauss(data, second=False):
  """Returns the velocity of the first or second peak from an adhoc-object in
      gaussian form.
      
      Usage: vf = VF_gauss(obj_gauss, second)
      
      vf: The returned velocity field
      obj_gauss: The adhoc-object in gaussian form
      second: False or True. If False the first peak is raturned, if True the
        second peak is returned.
      """
  
  # Create output array
  temp = N.resize(data, (data.nx(), data.ny()))
  
  # Check that the object is in gaussian form
  if (data.p['is_gauss']==True):
    if (second == False):
      # Compute the velocity for each point
      temp[:,:] = data[:,:,1]/data.p['lz']*data.fsr() + lamb2vel(data.p['xlbneb']-0.5*data.p['xil']) + data.p['vr_offset']
    elif (second == True):
      # Compute the velocity for each point
      temp[:,:] = data[:,:,4]/data.p['lz']*data.fsr() + lamb2vel(data.p['xlbneb']-0.5*data.p['xil']) + data.p['vr_offset']
    return temp
  else:
    print 'Object is not in gaussian form!'
  
def VF_firstmom(data, slitname=None):
  """ giving back a 2d map with the first moment of the 3d cube"""
  
  temp = N.resize(data, (data.nx(), data.ny()))
  
  temp[:,:] = doforeachpoint(data, firstmoment)/data.p['lz']*data.fsr() + lamb2vel(data.p['xlbneb']-0.5*data.p['xil']) + data.p['vr_offset']
   
  return temp

def shift(vec,i):
  """ Shift a vector.
      Usage: new_vec = shift(vec, i)
      
      vec:  The vector to be shifted
      i:  The steps to shift the vector with
      new_vec: The vector shifted i steps
      """
  temp = vec.get_copy()
  n= temp.size()
  i %= n
  temp[:] = N.concatenate((temp[n-i:n],temp[0:n-i]))
  return temp
  
def sortout(inarr,banned=0):
    """ give back three vectors: the values in inarr that are != value and the corresponding x and y coordinates """
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

def setYaxis_pc(inarr):
  """ set the plotting-axes to parsec """
  
  nticks=inarr.ny()*inarr.scale()/1000
  nticks=nticks.sum()
  factor=1.0
  
  while nticks > 10:
    nticks/=2
    factor/=2

  MP.setp(MP.gca(),'yticks',str(N.arange(nticks)/inarr.scale()*1000/factor))
  MP.setp(MP.gca(),'yticklabels',str(N.arange(nticks)/factor))
  MP.ylabel('[kpc]',font)

def setXaxis_pc(inarr):
  """ set the plotting-axes to parsec """
  
  nticks=inarr.nx()*inarr.scale()/1000
  nticks=nticks.sum()
  factor=1.

  while nticks > 10:
    nticks/=2
    factor/=2
  
  MP.setp(MP.gca(),'xticks',N.arange(nticks)/inarr.scale()*1000/factor)
  MP.setp(MP.gca(),'xticklabels',N.arange(nticks)/factor)
  MP.xlabel('[kpc]',font)

def dis(arr1,arr2):
  """ returns the distance between two points""" 
  arr=(arr2-arr1)**2
  return N.sqrt(arr.sum())

def app2abs(dist,m):
  """ convert apparent to absolute magnitudes, takes distance in Mpc"""
  return (-5*N.log10(dist*1000000.))+5+m

def chan2absvel(inarr):
    """ convert channels into absolute velocity"""
    #return lamb2vel(((inarr - (inarr.p['lz']/2.)) * inarr.p['xil'] /inarr.p['lz']) + inarr.p['xlbneb'] )
    return chan2relvel(inarr) + inarr.vel1st() + inarr.helioc()
    
def chan2relvel(inarr):
    """ convert channels into relative velocity"""
    #return lamb2vel((inarr * inarr.p['xil'] /inarr.p['lz'])+inarr.p['xlbneb']) - lamb2vel(inarr.p['xlbneb'])
    return inarr

def peakvel(inarr,n=3):
    """ calculate velocity field from the peak"""

    erg = doforeachpoint(inarr,calcpeak,n)
    erg.p=inarr.p.copy()
    return chan2absvel(erg)

def calcpeak(inarr,n):
    """ calculate the barycenter-velociy of the n highest pixels"""
    sorted,args=N.sort(inarr),N.argsort(inarr)
    erg = N.sum(sorted[-n:] * args[-n:]) / N.sum(sorted[-n:])
    return erg

def firstmoment(inarr):
  """ calculates the 1st moment from an array """
  shif=-inarr.argmax()+(inarr.size()/2)
  inarr[:]=shift(inarr,shif)
  return ((inarr * N.arange(inarr.size())).sum() / float(inarr.sum())) - shif

def secondmoment(inarr):
  """ calculates the 2nd moment from an array xx"""
  inarr=shift(inarr,-inarr.argmax()+(inarr.size()/2))
  return (((N.arange(inarr.size()) - firstmoment(inarr))**2) * inarr).sum() / float(inarr.sum())

def thirdmoment(inarr):
  """ calculates the 2nd moment from an array xx"""
  inarr=shift(inarr,-inarr.argmax()+(inarr.size()/2))
  return (((N.arange(inarr.size()) - firstmoment(inarr))**3) * inarr).sum() / float(inarr.sum())

def fwhm(inarr):
  """ returns the full widh at half maximum"""
  return 2.3548200450309493*N.sqrt(secondmoment(inarr))

def contfrommin(inarr,n):
    """ create a continuum map by averaging the n lowest channels"""
    return doforeachpoint(inarr,getavmin,n)
    
def getavmin(inarr,n):
    """ average the n lowest value in an array"""
    return N.sort(inarr)[0:n].mean()

def doforeachpoint(data, function, *args_extra):
  """Apply a function to all points, one at a time, in an array. The output
      will have a z-dimension equal to the length of the output from the
      'function'.
      
      Usage: new_arr = doforeachpoint(arr, function, *args_extra)
      
      arr: The 2D- or 3D-array
      function: The function to apply to each point in arr
      *args_extra: Arguments to be sent to 'function'
      new_arr: The output array
      """
  temp = data.get_copy()
  temp2 = data.get_copy()
  dim = data.nx()*data.ny()
  temp.setshape((dim, data.nz()))
  
  # First point
  args = (temp[0,:],) + args_extra
  point_zero = N.array(apply(function, args))
  
  # Resize to the length of the return from the function
  temp2 = N.resize(temp2, (dim, point_zero.size()))
  
  # Set first point
  temp2[0] = point_zero
  
  # Loop through the rest of the points
  for i in N.arange(1,dim):
    args = (temp[i,:],) + args_extra
    temp2[i] = N.array(apply(function, args))
  
  # This is ugly, but I can't find anything else that works
  if (point_zero.size() == 1):
    temp2.setshape((data.nx(), data.ny()))
  else:
    temp2.setshape((data.nx(), data.ny(), point_zero.size()))
  
  return temp2

def mask(bemasked,mask,lower=None,upper=None,value=0):
  """
  masks one array with another (numarray repeats the mask
  if it does not match in size)

  attention: Inplace masking, you may want to give a .copy()

  conditions say where the data are replaced by "value"
  """
  if (lower==None): 
    lower=bemasked.p['minmask']

  if ((lower != None) and (upper != None)):
    condition=N.where(mask > upper,1,0) | N.where(mask < lower,1,0)
  elif (upper != None):
    condition=N.where(mask > upper,1,0)
  elif (lower != None):
    condition=N.where(mask < lower,1,0)
  else:
    condition=0
    
  N.putmask(bemasked,condition,value)
  
def resample_vf(data, m):
  temp = data.get_copy()
  
  new_len_x = int(data.nx()/m)
  new_len_y = int(data.ny()/m)
  
  temp = temp.resize((new_len_x,new_len_y))
  for x in (N.arange(new_len_x)):
    for y in (N.arange(new_len_y)):
      temp[x,y] = data[m*x:m*x+m,m*y:m*y+m].mean()

  temp.p['dyncen'][0] = int(temp.p['dyncen'][0]/m)
  temp.p['dyncen'][1] = int(temp.p['dyncen'][1]/m)
  
  return temp
  
def LPfilter_vf(data):
  temp = data.get_copy()

  len_x = data.nx()
  len_y = data.ny()
  
  mvf = N.ones((3,3))/9.0
  temp[0,0] = (mvf*data[0:3,0:3]).sum()
  temp[0,len_y-1] = (mvf*data[0:3,len_y-3:len_y]).sum()
  temp[len_x-1,0] = (mvf*data[len_x-3:len_x,0:3]).sum()
  temp[len_x-1,len_y-1] = (mvf*data[len_x-3:len_x,len_y-3:len_y]).sum()
  
  mvf = N.ones((4,3))/12.0
  temp[1,0] = (mvf*data[0:4,0:3]).sum()
  temp[1,len_y-1] = (mvf*data[0:4,len_y-3:len_y]).sum()
  temp[len_x-2,0] = (mvf*data[len_x-4:len_x,0:3]).sum()
  temp[len_x-2,len_y-1] = (mvf*data[len_x-4:len_x,len_y-3:len_y]).sum()
  
  mvf = N.ones((3,4))/12.0
  temp[0,1] = (mvf*data[0:3,0:4]).sum()
  temp[0,len_y-2] = (mvf*data[0:3,len_y-4:len_y]).sum()
  temp[len_x-1,1] = (mvf*data[len_x-3:len_x,0:4]).sum()
  temp[len_x-1,len_y-2] = (mvf*data[len_x-3:len_x,len_y-4:len_y]).sum()
  
  mvf = N.ones((4,4))/16.0
  temp[1,1] = (mvf*data[0:4,0:4]).sum()
  temp[1,len_y-2] = (mvf*data[0:4,len_y-4:len_y]).sum()
  temp[len_x-2,1] = (mvf*data[len_x-4:len_x,0:4]).sum()
  temp[len_x-2,len_y-2] = (mvf*data[len_x-4:len_x,len_y-4:len_y]).sum()
  
  mvf = N.ones((5,3))/15.0
  for x in (N.arange(len_x-4)+2):
    temp[x,0] = (mvf*data[x-2:x+3,0:3]).sum()
    temp[x,len_y-1] = (mvf*data[x-2:x+3,len_y-3:len_y]).sum()
  
  mvf = N.ones((5,4))/20.0
  for x in (N.arange(len_x-4)+2):
    temp[x,1] = (mvf*data[x-2:x+3,0:4]).sum()
    temp[x,len_y-2] = (mvf*data[x-2:x+3,len_y-4:len_y]).sum()
  
  mvf = N.ones((3,5))/15.0
  for y in (N.arange(len_y-4)+2):
    temp[0,y] = (mvf*data[0:3,y-2:y+3]).sum()
    temp[len_x-1,y] = (mvf*data[len_x-3:len_x,y-2:y+3]).sum()
  
  mvf = N.ones((4,5))/20.0
  for y in (N.arange(len_y-4)+2):
    temp[1,y] = (mvf*data[0:4,y-2:y+3]).sum()
    temp[len_x-2,y] = (mvf*data[len_x-4:len_x,y-2:y+3]).sum()
  
  mvf = N.ones((5,5))/25.0
  for x in (N.arange(len_x-4)+2):
    for y in (N.arange(len_y-4)+2):
      temp[x,y] = (mvf*data[x-2:x+3,y-2:y+3]).sum()
  
  return temp


def fixsaturation():
    pass

#
# END: TOOLS


# PHYSICAL FUNCTIONS
#
def balmer(m):
  """ caculate m'th balmer line"""
  return hydrogen(2,m+2)
  
def hydrogen(n,m):
  """ calculate rydberg wavelengths of hydrogen"""
  m,n=float(m),float(n)
  return (n**2 * m**2 / (m**2 - n**2)) / 1.09677583E7

def massKepler(r,v):
  """ returns the mass from keplers law: M(<r)=r*(v**2)*G
   input in pc and km/s
   mind the factors 1/2 !!
  """
  return ((v/2)**2)*(r/2)*pc/G

def lamb2vel(l):
  """ converts a wavelength in A wrt HA into a radial velocity """
  return ((l/lambHA)-1)*sol

def flux2mag(f):
  return -2.5*(M.log10(f))

def mag2flux(m):
  return M.pow(10,((m)/2.5))

def absMag(app,dist):
  return (app + 5 - (5*M.log10(dist)))

def z2vel(z):
  return z*sol

def vel2z(v):
  return v/sol

def lamb2freq(l):
  return sol/l

def freq2lamb(f):
  return sol/f

def units(have,want,number=''):
  """
  uses the external procram "units" to convert units :-)
  """
  out=commands.getoutput('units -q -s ' + str(number) + have + ' ' + want + '| head -1 | cut -d " " -f2')
  return N.array([float(out)])
#
# END: PHYSICAL FUNCTIONS


# IN CASE SOMEONE TRIES TO EXECUTE THIS FILE
#
def demo():
    print "This file is the main file of a package. It is not meant to be executed. Import it instead!"

if __name__ == '__main__':
    demo()
#
# END: IN CASE SOMEONE TRIES TO EXECUTE THIS FILE
