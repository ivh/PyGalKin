"""
 InOutput.py

 For documentation see the major file PyCigale.py

 This file contains the file-in/output functions for the PyCigale package

"""
from PyGalKin import *

## SHORTCUTS
tab='\t'
nl='\n'
bs='\\'
sp=' '
null='\0'

def dump(data,filename):
  file=open(filename,'w')
  pickle.dump(data,file)
  file.close()

def load(filename):
  file=open(filename,'r')
  data=pickle.load(file)
  file.close()
  return data

def fromPAR(p={},filename='par'):
  """ reading my own parameter file with relevant parameters """

  p['parname']=filename
  for line in open(filename):
    line=line.split()
    if line[0] == 'cen':  p['cen'] = int(line[1]),int(line[2])
    elif line[0] == 'relsize':  p['relsize']=int(line[1])
    elif line[0] == 'pa':  p['pa']=float(line[1])
    elif line[0] == 'dyncen':  p['dyncen']= int(line[1]),int(line[2])
    elif line[0] == 'monocuts':  p['monocuts']= int(line[1]),int(line[2])
    elif line[0] == 'contcuts':  p['contcuts']= int(line[1]),int(line[2])
    elif line[0] == 'velcuts':  p['velcuts']= int(line[1]),int(line[2])
    elif line[0] == 'sigcuts':  p['sigcuts']= int(line[1]),int(line[2])
    elif line[0] == 'minmask':  p['minmask']=float(line[1])
    elif line[0] == 'incl':  p['incl']=float(line[1])
    elif line[0] == 'wedge':  p['wedge']=float(line[1])
    elif line[0] == 'gid':  p['gid']= int(line[1])
    elif line[0] == 'coor': p['coor'] = float(line[1]),float(line[2])
  return p

def toPAR(inarr,filename=None):
  """ writing my own parameter file with relevant parameters """
  if filename == None:
    if inarr.parname == None:
      print("I know no filename")
      return 1
    else:
      filename=inarr.parname

  file=open(filename,'w')

  file.write('cen'+sp+str(inarr.p['cen'][0])+sp+str(inarr.p['cen'][1])+nl)
  file.write('relsize'+sp+str(inarr.p['relsize'])+nl)
  file.write('pa'+sp+str(inarr.p['pa'])+nl)
  file.write('dyncen'+sp+str(inarr.p['dyncen'][0])+sp+str(inarr.p['dyncen'][1])+nl)
  file.write('monocuts'+sp+str(inarr.p['monocuts'][0])+sp+str(inarr.p['monocuts'][1])+nl)
  file.write('contcuts'+sp+str(inarr.p['contcuts'][0])+sp+str(inarr.p['contcuts'][1])+nl)
  file.write('velcuts'+sp+str(inarr.p['velcuts'][0])+sp+str(inarr.p['velcuts'][1])+nl)
  file.write('minmask'+sp+str(inarr.p['minmask'])+nl)

  file.close()

def fromAD(filename, readparams=True):
  """
  reading files from adhoc file format, setting the parameters right. returns adhoc-class
  """
  tmpF=N.fromfile(filename,dtype='Float32')
  tmpI=N.fromfile(filename,dtype='Int32')
  ndim=tmpI[-64]
  nx=tmpI[-61]
  ny=tmpI[-60]
  nz=tmpI[-59]

  if ndim == 2:
    shape=(nx,ny)
    data = N.fromfile(filename,dtype='Float32',count=nx*ny).reshape(shape)
    data = data[::-1,:]
  elif ndim ==3:
    shape=(nx,ny,nz)
    data = N.fromfile(filename,dtype='Float32',count=nx*ny*nz).reshape(shape)
    data = data[::-1,:,:]
  else:
    pass


  data=data.swapaxes(0,1)    # to be make origin at lower right

  # setting some parameters
  p={}
  p['cen']=nx/2,ny/2
  p['imagename']=filename

  # setting the Adhoc-parameters
  p['echelle']=float(tmpF[-58])
  p['xl1']=float(tmpF[-54])
  p['xil']=float(tmpF[-53])
  p['vr0']=float(tmpF[-52])
  p['p']=float(tmpF[-50])
  p['xlp']=float(tmpF[-49])
  p['xlneb']=float(tmpF[-48])
  p['v1']=float(tmpF[-47])
  p['interfr_kms']=float(tmpF[-46])

  if (readparams==True):
    fromADP(p,'par.adp')
    fromPAR(p,'par')
    #data.getsaltzer()

  data=N.ascontiguousarray(data)
  return C.adhoc(data,p=p)

def toAD(input,filename=None):
  """ writes ADHOC files. for file format see: adw_what.htm
  by default uses original filename
  """
  inarr=input.copy()

  if filename == None:
    if inarr.p['imagename'] == None:
      print("I know no filename")
      return 1
    else:
      filename=inarr.p['imagename']

  file=open(filename,'w')
  inarr=inarr.swapaxes(0,1)   # flipping axes again to be consistent with ADHOC
  inarr=inarr[::-1,:,:]
  inarr.astype('Float32').tofile(file)                         # the data
  N.array([inarr.ndim()],'Int32').tofile(file) # number of dimensios
  N.array([0],'Int32').tofile(file)          # id1 what is this?
  N.array([0],'Int32').tofile(file)          # id2
  N.array([inarr.nx()],'Int32').tofile(file)   # x-dim
  N.array([inarr.ny()],'Int32').tofile(file)   # y-dim
  N.array([inarr.nz()],'Int32').tofile(file)   # z-dim
  N.array([inarr.p['echelle']],'Float32').tofile(file)  # pixel scale
  N.array([0],'Int32').tofile(file)      # ix0
  N.array([0],'Int32').tofile(file)      # iy0
  N.array([1],'Float32').tofile(file)      # zoom

  if inarr.ndim() == 2:
    print("Writing 2D file")
    N.array([0],'Int32').tofile(file)
    N.array([0],'Float32').tofile(file)
    N.array([1],'Float32').tofile(file)
    N.array([256],'Int32').tofile(file)
    N.array([1],'Int32').tofile(file)
    N.arange(17,dtype='Int32').tofile(file) # unsused 68 bytes

  elif inarr.ndim() == 3:
    print("Writing 3D file")
    N.array([inarr.p['xl1']],'Float32').tofile(file)
    N.array([inarr.p['xil']],'Float32').tofile(file)
    N.array([inarr.p['vr0']],'Float32').tofile(file)
    N.array([0],'Float32').tofile(file) # velocity correction
    N.array([inarr.p['p']],'Float32').tofile(file)
    N.array([inarr.p['xlp']],'Float32').tofile(file)
    N.array([inarr.p['xlneb']],'Float32').tofile(file)
    N.array([inarr.p['v1']],'Float32').tofile(file)
    N.array([inarr.p['interfr_kms']],'Float32').tofile(file)
    N.arange(13,dtype='Int32').tofile(file) # unsused 52 bytes

  else:
    pass

  N.arange(32,dtype='Int32').tofile(file) # the 128 byte of comment

  file.close()

def fromADP(p,filename):
  """ read an ADP file into an ADP class
  """

  p['paraname']=filename

  for line in open(filename).readlines():

    line=line.split()
    if line.__len__() <= 1: p['objname']=line[0]
    elif line[1] == 'lx':  p['lx']=int(line[0])
    elif line[1] == 'ly': p['ly']=int(line[0])
    elif line[1] == 'lz': p['lz']=int(line[0])
    elif line[1] == 'ix0': p['ix0']=int(line[0])
    elif line[1] == 'iy0': p['iy0']=int(line[0])
    elif line[1] == 'iz0': p['iz0']=int(line[0])
    elif line[1] == 'modefp': p['modefp']=int(line[0])
    elif line[1] == 'p': p['p']=int(line[0])
    elif line[1] == 'xlp': p['xlp']=float(line[0])
    elif line[1] == 'xleta': p['xleta']=float(line[0])
    elif line[1] == 'xlbeta': p['xlbeta']=float(line[0])
    elif line[1] == 'xlneb': p['xlneb']=float(line[0])
    elif line[1] == 'xlbneb': p['xlbneb']=float(line[0])
    elif line[1] == 'xil': p['xil']=float(line[0])
    elif line[1] == 'xl1': p['xl1']=float(line[0])
    elif line[1] == 'vr0': p['vr0']=float(line[0])
    elif line[1] == 'corrv': p['corrv']=float(line[0])
    elif line[1] == 'echelle': p['echelle']=float(line[0])
    elif line[1] == 'xc': p['xc']=float(line[0])
    elif line[1] == 'yc': p['yc']=float(line[0])
    elif line[1] == 'ellipt': p['ellipt']=float(line[0])
    elif line[1] == 'mode_c': p['mode_c']=int(line[0])
    elif line[1] == 'mode_e': p['mode_e']=int(line[0])
    elif line[1] == 'min_a': p['min_a']=float(line[0])
    elif line[1] == 'max_a': p['max_a']=float(line[0])
    elif line[1] == 'min_r': p['min_r']=float(line[0])
    elif line[1] == 'max_r': p['max_r']=float(line[0])
    else: pass

def toADP(inADP,filename=None):
  """ write adp-class to adp-file
  uses original filename bu default
  """

  if filename == None:
    if inADP.p['paraname'] == None:
      print("I know no filename")
      return 1
    else:
      filename=inADP.p['paraname']

  file=open(filename,'w')

  file.write(inADP.p['objname'] +nl)
  file.write(tab + str(inADP.p['lx']) + tab + 'lx' + tab + 'dimension X' +nl)
  file.write(tab + str(inADP.p['ly']) + tab + 'ly' + tab + 'dimension Y' +nl)
  file.write(tab + str(inADP.p['lz']) + tab + 'lz' + tab + 'dimension Z' +nl)
  file.write(tab + str(inADP.p['ix0']) + tab + 'ix0' + tab + 'corner X' +nl)
  file.write(tab + str(inADP.p['iy0']) + tab + 'iy0' + tab + 'corner Y' +nl)
  file.write(tab + str(inADP.p['iz0']) + tab + 'iz0' + tab + 'corner Z' +nl)
  file.write(tab + str(inADP.p['modefp']) + tab + 'modefp' + tab + '1=PF,2=tabectro' +nl)
  file.write(tab + str(inADP.p['p']) + tab + 'p' + tab + 'PF interference order' +nl)
  file.write(tab + str(inADP.p['xlp']) + tab + 'xlp' + tab + 'lmd reference p' +nl)
  file.write(tab + str(inADP.p['xleta']) + tab + 'xleta' + tab + 'lmd calibration line' +nl)
  file.write(tab + str(inADP.p['xlbeta']) + tab + 'xlbeta' + tab + 'lmd calibration scan' +nl)
  file.write(tab + str(inADP.p['xlneb']) + tab + 'xlneb' + tab + 'lmd zero object' +nl)
  file.write(tab + str(inADP.p['xlbneb']) + tab + 'xlbneb' + tab + 'lmd object scan' +nl)
  file.write(tab + str(inADP.p['xil']) + tab + 'xil' + tab + 'tabectral length/interfringe in Angstroems' +nl)
  file.write(tab + str(inADP.p['xl1']) + tab + 'xl1' + tab + 'lmd channel 1' +nl)
  file.write(tab + str(inADP.p['vr0']) + tab + 'vr0' + tab + 'mean RV of the object in km/s' +nl)
  file.write(tab + str(inADP.p['corrv']) + tab + 'corrv' + tab + 'heliocentric correction to add to measures in km/s' +nl)
  file.write(tab + str(inADP.p['echelle']) + tab + 'echelle' + tab + 'scale ("/pix)' +nl)
  file.write(tab + str(inADP.p['xc']) + tab + 'xc' + tab + 'center rings X' +nl)
  file.write(tab + str(inADP.p['yc']) + tab + 'yc' + tab + 'center rings Y' +nl)
  file.write(tab + str(inADP.p['ellipt']) + tab + 'ellipt' + tab + 'X/Y pixel ratio' +nl)
  file.write(tab + str(inADP.p['mode_c']) + tab + 'mode_c' + tab + 'mode center/ellips of rings (0=normal; 1=enforced)' +nl)
  file.write(tab + str(inADP.p['mode_e']) + tab + 'mode_e' + tab + 'force ellips to 1 if CCD (0=normal; 1=enforced)' +nl)
  file.write(tab + str(inADP.p['min_a']) + tab + 'min_a' + tab + 'min angle for phase' +nl)
  file.write(tab + str(inADP.p['max_a']) + tab + 'max_a' + tab + 'max angle for phase' +nl)
  file.write(tab + str(inADP.p['min_r']) + tab + 'min_r' + tab + 'min radius for phase' +nl)
  file.write(tab + str(inADP.p['max_r']) + tab + 'max_r' + tab + 'max radius for phase' +nl)

  file.write("""     0       alpha hours
     0       alpha minutes
   0.0000000 alpha seconds
     0       delta degrees
     0       delta minutes
   0.0000000 alpha seconds
   0.0000000 reference year
   0.0000000 time minutes
     0       time hours
     0       time day
     0       time month
     0       time year
     0       dummy
     0       dummy
     0       dummy
     0       dummy
     0       dummy
   0.0000000 dummy
   0.0000000 dummy
   0.0000000 dummy
   0.0000000 dummy
   0.0000000 dummy
     0       numcom  last command number executed
""")

  file.close()

def read_fits (file):
  """
  gives back the data from a fits file
  """
  ima1=pyfits.open(file)
  ima=ima1[0]
  im=ima.data
  return im

def write_fits (data, name):
  """
  writes a matrix into a fits file
  no special headers so far
  """
  fitsfile=pyfits.HDUList()
  primary=pyfits.PrimaryHDU()
  primary.data=data
  fitsfile.append(primary)
  fitsfile.writeto(name)

def read_model(filename):
  """Reads filename.model and returns a model_list that can be used with
      functions in ModelVF.

      Usage: model_list = read_model(filename)
      """
  pars={}
  models=[]
  file=open(filename + '.model','r')

  # Read parameters and models
  for line in file:
    line=line.split()

    if len(line) == 0: pass
    elif line[0]=='model_sys': models += [str(line[1])]
    elif line[0]=='model_rot': models += [str(line[1])]
    elif line[0]=='model_exp': models += [str(line[1])]

    elif line[0]=='dim': pars['dim']=int(line[1])
    elif line[0]=='a_scale': pars['a_scale']=[float(line[1]),int(line[2])]
    elif line[0]=='v_system': pars['v_system']=[float(line[1]),int(line[2])]
    elif line[0]=='inclination': pars['inclination']=[float(line[1]),int(line[2])]
    elif line[0]=='pa': pars['pa']=[float(line[1]),int(line[2])]
    elif line[0]=='v_expansion': pars['v_expansion']=[float(line[1]),int(line[2])]
    elif line[0]=='v_max': pars['v_max']=[float(line[1]),int(line[2])]
    elif line[0]=='exp_max': pars['exp_max']=[float(line[1]),int(line[2])]
    elif line[0]=='r_max': pars['r_max']=[float(line[1]),int(line[2])]
    elif line[0]=='x0': pars['centr_offset_x']=[float(line[1]),int(line[2])]
    elif line[0]=='y0': pars['centr_offset_y']=[float(line[1]),int(line[2])]

  file.close()

  # Create the model_list
  model_list = []
  for i in range(len(models)):
    model_list += [[models[i], pars]]

  return model_list

def write_model(filename, model_list):
  """Writes the model_list to filename.model in the .model-format. Note: Only
      one model of each type (system, rotation and expansion) is supported.

      Warning! If the file exists it will be overwritten!

      Usage: write_model(filename, model_list)
      """

  # Open the file for writing
  file=open(filename + '.model', 'w')

  # Write the models in the list
  for i in range(len(model_list)):
    if (model_list[i][0] == 'system'): file.write('model_sys'+tab+model_list[i][0]+nl)
    elif (model_list[i][0] == 'expansion'): file.write('model_exp'+tab+model_list[i][0]+nl)
    elif (model_list[i][0] == 'linear'): file.write('model_rot'+tab+model_list[i][0]+nl)
    elif (model_list[i][0] == 'disk'): file.write('model_rot'+tab+model_list[i][0]+nl)
    elif (model_list[i][0] == 'kepler'): file.write('model_rot'+tab+model_list[i][0]+nl)
    elif (model_list[i][0] == 'pure_kepler'): file.write('model_rot'+tab+model_list[i][0]+nl)
    elif (model_list[i][0] == 'expansion'): file.write('model_rot'+tab+model_list[i][0]+nl)

  # Write the model-parameters
  file.write(nl)
  file.write('dim'+tab+tab+tab+str(model_list[0][1]['dim'])+nl)
  file.write('v_system'+tab+str(model_list[0][1]['v_system'][0])+tab+tab+str(model_list[0][1]['v_system'][1])+nl)
  file.write('a_scale'+tab+tab+str(model_list[0][1]['a_scale'][0])+tab+tab+str(model_list[0][1]['a_scale'][1])+nl)
  file.write('inclination'+tab+str(model_list[0][1]['inclination'][0])+tab+tab+str(model_list[0][1]['inclination'][1])+nl)
  file.write('pa'+tab+tab+tab+str(model_list[0][1]['pa'][0])+tab+tab+str(model_list[0][1]['pa'][1])+nl)
  file.write('v_expansion'+tab+str(model_list[0][1]['v_expansion'][0])+tab+tab+str(model_list[0][1]['v_expansion'][1])+nl)
  file.write('v_max'+tab+tab+str(model_list[0][1]['v_max'][0])+tab+tab+str(model_list[0][1]['v_max'][1])+nl)
  file.write('exp_max'+tab+tab+str(model_list[0][1]['exp_max'][0])+tab+tab+str(model_list[0][1]['exp_max'][1])+nl)
  file.write('r_max'+tab+tab+str(model_list[0][1]['r_max'][0])+tab+tab+str(model_list[0][1]['r_max'][1])+nl)
  file.write('x0'+tab+tab+tab+str(model_list[0][1]['centr_offset_x'][0])+tab+tab+str(model_list[0][1]['centr_offset_x'][1])+nl)
  file.write('y0'+tab+tab+tab+str(model_list[0][1]['centr_offset_y'][0])+tab+tab+str(model_list[0][1]['centr_offset_y'][1])+nl)

  file.close()

def read_slit (filename):
  """Read slit parameters from a filename.slit. Returns a dictionary with the
      parameters.

      Usage: pars = read_slit(filename)
  """
  pars={}
  file=open(filename + '.slit','r')

  for line in file:
    line=line.split()
    if len(line) == 0: pass
    elif line[0]=='slitwidth': pars['slitwidth']=float(line[1])
    elif line[0]=='angle': pars['angle']=float(line[1])
    elif line[0]=='offset': pars['offset']=N.array([float(line[1]), float(line[2])])
    elif line[0] == 'slit_angle_dep': pars['slit_angle_dep']=int(line[1])

  file.close()
  return pars

def write_rc (filename,position, velocity, error):
  """Writes an rc output-file (from RCslit or gauss_slit).

      Usage: write_rc(filename, pos, vel, err)

      pos, vel and err are arrays with the positions, velocities and errors.
  """
  outfile=open(filename + '.rc','w')

  for i in range(len(position)):
    outfile.write('%g %g %g \n' % (position[i], velocity[i],error[i]))

  outfile.close()



#### NOT WORKING RIGHT NOW! LEGACY!!!
def ADTtoADP(outname,objname,ADT):
    """ make a parameter file, using values from an ADT file
        NOT WORKING RIGHT NOW! LEGACY!!!
    """

    outfile=open(outname,'w')
    adt=open(ADT+'.ADT','r')

    outfile.write(objname+nl)
    adt.readline()
    if adt.readline().split()[1] != ADT: print("observation numer mismatch")
    adt.readline()
    if adt.readline().split()[2] != objname: print("object name mismatch")
    pf=adt.readline().split()[3]
    filter=adt.readline().split()[2]
    adt.readline()
    tmp=adt.readline().split()
    nx,ny=int(tmp[5]),int(tmp[6])
    nz=int(adt.readline().split()[4])
    adt.readline()
    exppchan=float(adt.readline().split()[5])
    scanlamb=float(adt.readline().split()[2])
    queen=float(adt.readline().split()[2])

    ## extracting the number of cycles
    tmp=adt.readlines()
    tmp.reverse()
    for i in range(5):
        if ' cy=' in tmp[i]:
            numbercycles=int(tmp[i].split()[5].split('=')[1])
            break


    # GETTING INFO FROM .ADT FILE FINISHED
    adt.close()

    print(pf,filter,nx,ny,nz,exppchan,scanlamb,queen)

    outfile.write(tab+str(nx)+tab+'lx'+tab+'dimension X'+nl)
    outfile.write(tab+str(ny)+tab+'lx'+tab+'dimension X'+nl)
    outfile.write(tab+str(nz)+tab+'lx'+tab+'dimension X'+nl)
    outfile.write(tab+"""1       ix0     corner X
        1       iy0     corner Y
        1       iz0     corner Z
        1       modefp  1=PF,2=spectro""" + nl)
    if pf == 'OM798': outfile.write(tab+'793'+tab+'p'+tab+'PF interference order'+nl)
    elif pf == 'OM1938': outfile.write(tab+'1983'+tab+'p'+tab+'PF interference order'+nl)
    else: outfile.write(tab+'0'+tab+'p'+tab+'PF interference order'+nl)
    outfile.write(tab+"""6562.7797852 xlp     lmd reference p
    6598.9531250 xleta   lmd calibration line
    6598.9531250 xlbeta  lmd calibration scan
    6562.7797852 xlneb   lmd zero object"""+nl)
    outfile.write(tab+str(scanlamb)+tab+'xlbneb'+tab+'lmd object scan')
    bulk="""
       8.3223772 xil     spectral length/interfringe in Angstroems
    6583.0000000 xl1     lmd channel 1
    1099.0000000 vr0     mean RV of the object in km/s
       0.0000000 corrv   heliocentric correction to add to measures in km/s
       0.4170000 echelle scale (\"/pix)
        256.0000 xc      center rings X
        256.0000 yc      center rings Y
       1.0000000 ellipt  X/Y pixel ratio
         0       mode_c  mode center/ellips of rings (0=normal; 1=enforced)
         1       mode_e  force ellips to 1 if CCD (0=normal; 1=enforced)
          0.0000 min_a   min angle for phase
        360.0000 max_a   max angle for phase
          0.0000 min_r   min radius for phase
      99999.0000 max_r   max radius for phase
         0       alpha hours
         0       alpha minutes
       0.0000000 alpha seconds
         0       delta degrees
         0       delta minutes
       0.0000000 alpha seconds
       0.0000000 reference year
       0.0000000 time minutes
         0       time hours
         0       time day
         0       time month
         0       time year
         0       dummy
         0       dummy
         0       dummy
         0       dummy
         0       dummy
       0.0000000 dummy
       0.0000000 dummy
       0.0000000 dummy
       0.0000000 dummy
       0.0000000 dummy
         0       numcom  last command number executed
    """

    outfile.write(bulk)
    outfile.close()

    return numbercycles



