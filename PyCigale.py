"""
PyCigale.py

"""

from PyGalKin import *

class adhoc(numpdict):
    def __new__(subtype, data, p=None, dtype=None, copy=False):
        subarr=numpdict.__new__(subtype, data, p, dtype, copy)
        subarr.calcstuff()
        return subarr

    def cenx(self):
        return self.p['cen'][0]
    
    def ceny(self):
      return self.p['cen'][1]

    def mymax(self):
        if self.ndim != 3: return -1
        else:
            erg=N.zeros(self.shape[:-1])
            for i in N.arange(erg.shape[0]):
                for j in N.arange(erg.shape[1]):
                    erg[i,j]=self[i,j,:].max()
        return erg
    
    def mymin(self):
        if self.ndim != 3: return -1
        else:
            erg=N.zeros(self.shape[:-1])
            for i in N.arange(erg.shape[0]):
                for j in N.arange(erg.shape[1]):
                    erg[i,j]=self[i,j,:].min()
        return erg

    def calcstuff(self):
        try: self.p['fsr']=lamb2vel(self.p['xl1'] + self.p['xil']) - lamb2vel(self.p['xl1'])
        except: pass
        

    def fsr(self):
      """ calculate the free spectral range from self.xil and self.xlp
          uses lamb2vel
      """
      return self.p['fsr']
    
    def vel1st(self):
        """return the velocity of the first channel """
        #return lamb2vel(self.p['xlbneb']) - (self.fsr() / 2.)
        return lamb2vel(self.p['xl1'])

    def wcal(self):
        """ returns vector of length lz with the velocities for each channel """
        if  self.nz() != 1:
            return (N.arange(self.nz(),dtype='f')/self.nz()*self.fsr())+self.p['xl1']
        
        
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
      return vel2dis(self.p['vr0'])
  
    def helioc(self):
      """ return the value for helocentric correction"""
      return self.p['corrv']
    
    def scale(self):
      """returns the physical scale for the object (pc/pix) """
      return scalefromvarc(self.p['echelle'],self.p['vr0'])
    
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
      outarr.p['cen']=N.array([relsize,relsize])
      outarr.p['dyncen']=outarr.p['cen']+(self.p['dyncen']-self.p['cen'])
      return outarr
    
    def mask(self,cond):
      """ wrapper for mask-funktion
          works inplace!
      """
      self=masked_where(cond,self)
    
   
    def PVdiag(self):
      """Calculate the position-velocity diagram from a velocity field.
          
          Usage: pos,vel = vf.PVdiag()
          
          vf: The velocity field to calculate from
          pos: A vector of positions
          vel: A vector of corresponding velocities
          """

      pos,vel=posvel(self,pa=self.p['pa'],dyncen=self.p['dyncen'])
      return pos*self.scale()/1000

    
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
      perp=N.array([-N.sin(radians(angle)),N.cos(radians(angle))])
      vec=-N.array([N.cos(radians(angle)),N.sin(radians(angle))])
      
      # Arrays for the data
      pos=N.zeros(1000,'Float32')
      vel=N.zeros(1000,'Float32')
      err=N.zeros(1000,'Float32')
      count = 0
      for i in N.arange(self.nx()):
        for j in N.arange(self.ny()):
          if (self[i,j] != 0):
            # The points must be within the slit
            dist=N.inner(perp,N.array([center[0]-i,center[1]-j]))
            if (abs(dist) <= D):
              pos[count]=N.inner(vec,N.array([center[0]-i,center[1]-j]))
              vel[count]=self[i,j]
              err[count]=self.fsr()/self.nz()
              count+=1
      
      # Write rc-outfile if a filename is given
      if (outfile != None):
        write_rc(outfile, pos, vel, err)
      
      return pos[:count],vel[:count]
    
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
  return N.fromfile(*args, **keys).view(adhoc)
  
def fromstring(*args, **keys):
  return N.fromstring(*args, **keys).view(adhoc)
  

def array(*args, **keys):
  return N.array(*args, **keys).view(adhoc)

def stripadhoc(arr):
  return arr.view(N.ndarray)

def sum(arr,axis=2, dtype=None, out=None):
  return N.sum(arr,axis=axis, dtype=dtype, out=out)

#
# END: WRAPPERS



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
  temp = adhoc(temp,data.p)
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

def chan2absvel(inarr):
    """ convert channels into absolute velocity.."""
    return lamb2vel( inarr.p['xl1'] + (inarr*inarr.p['xil']/inarr.p['lz']) )
    
def chan2relvel(inarr):
    """ convert channels into relative velocity"""
    #return lamb2vel((inarr * inarr.p['xil'] /inarr.p['lz'])+inarr.p['xlbneb']) - lamb2vel(inarr.p['xlbneb'])
    return lamb2vel(lambHA + (inarr * inarr.p['xil'] /float(inarr.p['lz'])))

def peakvel(inarr,n=3):
    """ calculate velocity field from the peak"""
    erg = doforeachpoint(inarr,calcpeak,n)
    return chan2absvel(erg)

def contfrommin(inarr,n=5):
    """ create a continuum map by averaging the n lowest channels"""
    return doforeachpoint(inarr,getavmin,n)
    
def getavmin(inarr,n):
    """ average the n lowest value in an array"""
    return N.sort(inarr)[0:n].mean()

def monofromcont(inarr,n=5):
    csum=sum(inarr,axis=2)
    cont=contfrommin(inarr,n)
    return csum-(cont*inarr.nz())


def resample_vf(data, m):
  temp = data.copy()
  
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
  temp = data.copy()

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



#
# END: TOOLS




# IN CASE SOMEONE TRIES TO EXECUTE THIS FILE
#
def demo():
    print "This file is the main file of a package. It is not meant to be executed. Import it instead!"

if __name__ == '__main__':
    demo()
#
# END: IN CASE SOMEONE TRIES TO EXECUTE THIS FILE
