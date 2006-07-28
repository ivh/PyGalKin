# IMPORTS
#
import numarray as N
import os
import copy
import commands
import matplotlib
import matplotlib.pylab as MP
#matplotlib.use('GTK')
import Numeric
from mpfit.mpfit import mpfit
from math import pi,e,radians
#
# END IMPORTS

from InOutput import *

# CONSTANTS
#
# units:
#        velocity: km/s
#        wavelength: Angstrom
lambHA=N.array([6562.7797852000003],'Float32')
sol=N.array([299792.458],'Float32')
H0=N.array([75.],'Float32')
G=N.array([6.6726E-11*1.989E30/1000**3],'Float32') ### in solar masses and km
pc=N.array([3.086E13],'Float32') ## in km
#
# END CONSTANTS

def lamb2vel(l):
  """Convert a wavelength in A wrt HA into a radial velocity."""
  return ((l/lambHA)-1)*sol
  
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

def on_click_float(event):
  """Handle click events in plotting windows. Not to be called from the command line.
      It saves the coordinates for the clicks with button 1 in a temporary file. Button 3 closes the figure.
      """

  file=open('/tmp/MPclick.dat','a')

  if (event.button == 1):
    if event.inaxes:
      #print event.xdata, event.ydata
      file.write('%f %f\n' % (event.xdata, event.ydata))
  elif (event.button == 2):
    pass
  elif (event.button == 3):
    MP.close()

  file.close()
  
def find_peaks(data):
  """Finds statistically significant peaks in a data set. If two peaks are found,
      the position of the peaks are returned in a list. A value of -1 in the list
      means that no peak were found.
      
      Usage:  peaks = find_peaks(data)
      
      data: A spectra to find peaks in.
      peaks: A list of length two with the two largest peaks.
      """
  # Temporary list to store all peaks in
  temp=[]
  # The list to be returned, set to a "no peak list" to begin with
  peaks=[-1,-1]
  # The limit used to decide if a peak is statistically significant or not
  limit = 0.1
  
  # Find all peaks in the form _-^-_
  for i in N.arange(data.size()):
    if ((shift(data,-1)[i] < data[i]) and (shift(data,1)[i] < data[i]) and (shift(data, -2)[i] < shift(data, -1)[i]) and (shift(data, 2)[i] < shift(data, 1)[i])):
      temp = temp + [i]
    elif ( (shift(data,-1)[i]  > data[i]) and (shift(data,1)[i]  > data[i]) and (shift(data, -3)[i] < shift(data, -2)[i]) and (shift(data, 3)[i] < shift(data, 2)[i]) and (shift(data, 2)[i] < shift(data, 1)[i]) and (shift(data, -2)[i] < shift(data, -1)[i])):
      temp = temp + [i]
  
  # Bubble sort, but only two times (enough to move the two biggest peaks to the end of the list)
  n = len(temp)
  if (n==0):
    return peaks
  elif (n==1):
    if ((data[temp[0]] - data.min()) > (limit*N.sqrt(data[temp[0]]))):
      peaks[0]=temp[0]
    return peaks
  else:
    for i in range(2):
      n=n-1
      for j in range(n):
        if(data[temp[j]]>data[temp[j+1]]):
          a=temp[j]
          temp[j]=temp[j+1]
          temp[j+1]=a
          
    
    n = len(temp)
    
    # The (peak - minimum) must be bigger than limit * sqrt(peak)
    if ((data[temp[n-1]] - data.min()) > (limit*N.sqrt(data[temp[n-1]]))):
      peaks[0] = temp[n-1]
    
    if ((data[temp[n-2]] - data.min()) > (limit*N.sqrt(data[temp[n-2]]))):
      peaks[1] = temp[n-2]
    
    return peaks

def count_double_single(data):
  """Test how many single, double and empty fits there are in a gauss fitted image. The result is printed to the terminal.
      Usage: test_double_single(data)
      
      data: A gauss fitted adhoc-object.
      """
  temp=data.get_copy()
  dim = temp.nx()*temp.ny()
  
  temp.setshape((dim,14))
  
  empty=0
  double=0
  single=0
  
  for i in range(dim):
    if ((temp[i][5] == 0) and (temp[i][2] == 0)):
      empty = empty + 1
    elif (temp[i][5] == 0):
      single = single +1
    else:
      double = double +1
  
  print 'Total points '+str(single+double+empty)
  print 'Single fits: '+str(single)
  print 'Double fits: '+str(double)
  print 'Empty fits: '+str(empty)
  
def gauss_to_spectra(data, array_length=None, first_peak=False, second_peak=False):
  """Converts all points from gaussian parameters to data arrays.
      
      Input:
        data = The object to convert
        array_length = The desired length of the data arrays
        
      Output:
        The object in array form
        """
  # If the data is not in gauss form, return a copy
  if (data.p['is_gauss']==False):
    print 'Object is not in gaussian form!'
    return data.get_copy()
  
  # Default to length lz
  if (array_length == None):
    array_length = data.p['lz']
  
  # Temporary variables (get_copy() ensures that properties of the object is kept)
  temp = data.get_copy()
  temp2 = data.get_copy()
  
  # Dimension and lengths
  dim = temp.nx()*temp.ny()
  length_x = temp.nx()
  length_y = temp.ny()
  length_z = temp.nz()
  
  # Reshape for convenience
  temp.setshape((dim,length_z))
  
  # Resize the arrays to the needed length
  temp2 = N.resize(temp2, (dim, array_length))
  
  # Convert the gauss parameters to arrays of desired length
  if (first_peak==True):
    for i in range(dim):
      params = [0,0,0,1,0,0,1]
      params[0] = temp[i][0]
      params[1] = temp[i][1]
      params[2] = temp[i][2]
      params[3] = temp[i][3]
      
      temp2[i] = gauss_to_array(params, array_length)
  elif (second_peak == True):
    for i in range(dim):
      params = [0,0,0,1,0,0,1]
      params[4] = temp[i][4]
      params[5] = temp[i][5]
      params[6] = temp[i][6]
      
      temp2[i] = gauss_to_array(params, array_length)
  else:
    for i in range(dim):
      params = [0,0,0,1,0,0,1]
      params[0] = temp[i][0]
      params[1] = temp[i][1]
      params[2] = temp[i][2]
      params[3] = temp[i][3]
      params[4] = temp[i][4]
      params[5] = temp[i][5]
      params[6] = temp[i][6]
      
      temp2[i] = gauss_to_array(params, array_length)
  
  # Object is no longer in gaussian form
  temp2.p['is_gauss']=False
  # Reshape
  temp2.setshape((length_x,length_y,array_length))
    
  return temp2

def spectra_to_gauss(data, double=False, try_double=False, do_shift=False, limit=5, manual=None):
  """
      Converts all points from data arrays to gaussian parameters. 
      Can fit a single or a double gauss to the data (the 'double' parameter). 
      If double is True and two peaks are found, a double fit is performed, otherwise a single fit is used.
      The 'try_double' parameter can be used to force a double fit, even if only one peak is found.
      You can select rectangular areas where you manually chooses the initial values for the peaks
      , this is done by using the 'manual' parameter.
      Only converts arrays with the max higher than the limit parameter,
      if it is lower it is set to zero in all points.
      
      Input:
        data = The object to convert
        double = True or False (default is False)
        try_double = True or False
        limit = The max value in the array must be above this to be converted
        manual = A list in the form [[(x0,y0),(x1,y1)],[(x2,y2),(x3,y3)],...].
          This set defines a rectangle with corners x0,y0 and x1,y1 (all points between x0-x1 and y0-y1),
          and a rectangle with corners x2,y2 and x3,y3. You can choose as many rectangles as you like.
        
      Output:
        The object in gauss form
        """
  # If the object is already in gauss form, return a copy
  if (data.p['is_gauss']==True):
    print 'Object already is in gaussian form!'
    return data.get_copy()
  
  # Temporary variables (copy() ensures that properties of the object is kept)
  temp = data.get_copy()
  temp2 = data.get_copy()
  
  # Dimension and lengths
  dim=temp.nx()*temp.ny()
  length_x = temp.nx()
  length_y = temp.ny()
  length_z = temp.nz()
  
  # Reshape for convenience
  temp.setshape((dim,length_z))
  # Resize the arrays to the needed length
  temp2 = N.resize(temp2, (dim,14))

  manual_area = N.zeros((length_x, length_y))
  if (manual != None):
    for i in range(len(manual)):
      manual_area[manual[i][0][1]:manual[i][1][1], manual[i][0][0]:manual[i][1][0]] = 1
    print manual_area
  manual_area.setshape((dim))

  # Only mark points with the maximum value greater than 'limit' for computation
  points=[]
  manual_points=[]
  nopoints=[]
  for i in range(dim):
    if (i%1000 == 0): print i
    if (temp[i].max() > limit) and (manual_area[i] == 1):
      manual_points = manual_points + [i,]
    elif (temp[i].max() > limit):
      points = points + [i,]
    else:
      nopoints = nopoints + [i,]
      
  
  print 'Interactive points: '+str(len(manual_points))
  
  # Compute the interactive fits
  if (len(manual_points) > 0):
    count=0
    for i in manual_points:
      count+=1
      if (count%100 == 0): print str(i)+' '+str(count)
      params, error = gauss_from_array_interactive(temp[i])
      temp2[i] = N.concatenate((params, error))
  
  
  print 'Automatic points: '+str(len(points))
  
  # Compute the gauss fits and store them in temp2
  if (len(points) > 0):
    count=0
    for i in points:
      count+=1
      #if (count%100 == 0): 
      print str(i)+' '+str(count)
      params, error = gauss_from_array(temp[i], double=double, try_double=try_double, do_shift=do_shift)
      temp2[i] = N.concatenate((params, error))
      #print 'Params: '+str(params)
      #print 'Error: '+str(error)

  # Set all other points to a zero-gauss
  if (len(nopoints) > 0):
    params = [0,0,0,1,0,0,1]
    error = [0,0,0,0,0,0,0]
    temp2[nopoints] = N.concatenate((params, error))

  # Object is in gaussian form
  temp2.p['is_gauss']=True
  # Reshape
  temp2.setshape((length_x,length_y,14))

  return temp2

def gauss_to_array(p, array_length=None):
  """Converts one set of gaussian parameters to a data array.
      
      Input:
        p = The parameters to convert
        array_length = The desired length of the data arrays
        
      Output:
        The data array computed from the gaussian parameters
        """
  
  # Create the "x-values"
  x = N.arange(float(array_length))

  # Calculate the exponents first
  arr_len = len(x)
  limit = -700

  temp1 = N.maximum( -((x-p[1])**2)/(2*p[3]**2) , limit)
  temp2 = N.maximum( -((x-p[4])**2)/(2*p[6]**2) , limit)
  temp1_right = N.maximum( -(((x+arr_len)-p[1])**2)/(2*p[3]**2) , limit)
  temp2_right = N.maximum( -(((x+arr_len)-p[4])**2)/(2*p[6]**2) , limit)
  temp1_left = N.maximum( -(((x-arr_len)-p[1])**2)/(2*p[3]**2) , limit)
  temp2_left = N.maximum( -(((x-arr_len)-p[4])**2)/(2*p[6]**2) , limit)

  # Compute the double gauss
  arr = p[0] + p[2]*Numeric.exp(temp1) + p[5]*Numeric.exp(temp2)
  arr_left = p[2]*Numeric.exp(temp1_left) + p[5]*Numeric.exp(temp2_left)
  arr_right = p[2]*Numeric.exp(temp1_right) + p[5]*Numeric.exp(temp2_right)
  
  model_y = arr + arr_left + arr_right
  
  # Reuse the array, set it to the computed double gauss and return it
  x[:] = model_y[:]
  
  return x

def gauss_from_array_interactive(arr):
  """Converts a data array to a set of gaussian parameters.
      You click in the data set where you want peaks to be fitted (two
      peaks maximum).
      
      If no peaks are selected, a set of parameters that generates zero 
      in all points is returned.
      
      If one distinct peak is found, a single gauss fit is performed, unless
      'try_double' is True, then a double fit is performed anyway.
      
      If 'double' is True and either two distinct peaks are found or 'try_double'
      is True, a double fit is performed (if 'double is False, 'try_double' has no
      effect).
      
      Always returns an array of 7 parameters: 
        y0, x01, A1, sigma1, x02, A2, sigma2
      If a single fit is performed the last 3 is set to 0, 0, 1.
      
      Input:
        arr = The array to be fitted
        double = True or False (default is False)
        try_double = True or False (default is False)
        
      Output:
        The gaussian parameters computed from the data array
        """
  # View the spectra and read clicks
  try:
    os.remove('/tmp/MPclick.dat')
  except:
    pass
  
  MP.plot(arr, 'bo')
  MP.connect('button_press_event', on_click_float)
  MP.show()
  
  # Read the choosen points
  points = []
  for line in open('/tmp/MPclick.dat', 'r').readlines():
    line = line.split()
    points += [[float(line[0]), float(line[1])]]
  
  # If no points are choosen
  if (len(points)==0):
    params = [0, 0, 0, 1, 0, 0, 1]
    params_error_red = [1, 1, 1, 1, 1, 1, 1]
    return params, params_error_red
  
  # Initial values
  if (arr.min()<0):
    minshift = arr.min()
    arr=arr-minshift
  else:
    minshift=0
  
  y0=arr.min()
  
  x01 = points[0][0]
  A1 = points[0][1]-y0
  sigma1 = 3.0
  
  if (len(points)>1):
    x02 = points[1][0]
    A2 = points[1][1]-y0
    sigma2 = 3.0
  
  # Setup the parameter list
  if (len(points)>1):
    number_of_pars = 7
  else:
    number_of_pars = 4
  
  parinfo = []
  for i in range(number_of_pars):
    parinfo.append({'value':0.0, 'fixed':0, 'limited':[0,0],'limits':[0.0, 0.0], 'step':0.0})
  
  parinfo[0]['value'] = y0
  parinfo[0]['limited'][0] = 1
  parinfo[0]['limits'][0] = 0.0 #float(N.minimum(0, arr.min()))
  parinfo[0]['limited'][1] = 1
  
  arr_sorted = N.sort(arr)
  arr_sorted_mean = arr_sorted[0:5].mean()
  arr_sorted_stddev = arr_sorted[0:5].stddev()
  y0_upper_limit = arr_sorted_mean + arr_sorted_stddev
  
  parinfo[0]['limits'][1] = arr.min() + 0.5*(arr.max()-arr.min()) #N.maximum(0, y0_upper_limit)
  
  parinfo[1]['value'] = x01
  parinfo[1]['limited'][0] = 1
  parinfo[1]['limits'][0] = float(x01-1)
  parinfo[1]['limited'][1] = 1
  parinfo[1]['limits'][1] = float(x01+1)
  
  parinfo[2]['value'] = A1
  parinfo[2]['limited'][0] = 1
  parinfo[2]['limits'][0] = float(A1-1)
  parinfo[2]['limited'][1] = 1
  parinfo[2]['limits'][1] = float(A1+1)
  
  parinfo[3]['value'] = sigma1
  parinfo[3]['limited'][0] = 1
  parinfo[3]['limits'][0] = 0.5
  parinfo[3]['limited'][1] = 1
  parinfo[3]['limits'][1] = float(arr.size())/2.0
  
  if (len(points)>1):
    # If this is a double fit, set initial values for the second gauss
    parinfo[4]['value'] = x02
    parinfo[4]['limited'][0] = 1
    parinfo[4]['limits'][0] = float(x02-1)
    parinfo[4]['limited'][1] = 1
    parinfo[4]['limits'][1] = float(x02+1)
    
    parinfo[5]['value'] = A2
    parinfo[5]['limited'][0] = 1
    parinfo[5]['limits'][0] = float(A2-1)
    parinfo[5]['limited'][1] = 1
    parinfo[5]['limits'][1] = float(A2+1)
  
    parinfo[6]['value'] = sigma2
    parinfo[6]['limited'][0] = 1
    parinfo[6]['limits'][0] = 0.5
    parinfo[6]['limited'][1] = 1
    parinfo[6]['limits'][1] = float(arr.size())/2.0
  
  # Setup arrays with the data to fit
  x = N.arange(float(arr.size())).astype('Float64')
  y = N.zeros(arr.shape, 'Float32')
  y[:] = arr[:]
  err = N.ones(arr.shape, 'Float32')
  
  functkw = {'x':x, 'y':y, 'err':err}
  
  # Fit parameters
  if (len(points)>1):
    m = mpfit(gauss_func_double, parinfo=parinfo, functkw=functkw, quiet=1)
    params = m.params
    params_error = m.perror
  else:
    m = mpfit(gauss_func_single, parinfo=parinfo, functkw=functkw, quiet=1)
    params = N.concatenate((m.params,[0,0,1]))
    params_error = N.concatenate((m.perror,[0,0,0]))

  # Reduced chi-square error
  dof = len(x) - len(m.params)
  params_error_red = params_error * N.sqrt(m.fnorm / dof)

  # Set the highest peaks parameters to the left
  if (params[2] < params[5]):
    params = N.concatenate((params[0], params[4:7], params[1:4]))
  
  # If one peak is very weak compared to the other one, 
  # or if the peaks are located close to each other (within half the sigma
  # of the stronger peak), make a single gauss fit.
  if (len(points)>1):
    # Set limits to use for QA of the fit
    limit_err = 500
    limit_amp = 0.05
    limit_width = 1.0
    
    if (params[5] < limit_amp*params[2]):
      params, params_error_red = gauss_from_array(arr, double=False)
    elif (N.abs(params[1]-params[4]) < N.maximum(limit_width*params[3], limit_width*params[6])):
      params, params_error_red = gauss_from_array(arr, double=False)
    else:
      # Compute the sum of squared errors
      temp=(err*(gauss_to_array(params, len(y))-y))**2
      
      # If the error is big, do a single gauss fit
      if (temp.sum() > limit_err):
        params, params_error_red = gauss_from_array(arr, double=False)
  
  params[0]+=minshift
  return params, params_error_red

def gauss_from_array(arr, double=False, try_double=False, do_shift=False):
  """Converts a data array to a set of gaussian parameters.
      Can fit a single or a double gauss to the data (the double parameter).
      
      If no distinct peak is found in the data set, a set of parameters that
      generates zero in all points is returned.
      
      If one distinct peak is found, a single gauss fit is performed, unless
      'try_double' is True, then a double fit is performed anyway.
      
      If 'double' is True and either two distinct peaks are found or 'try_double'
      is True, a double fit is performed (if 'double is False, 'try_double' has no
      effect).
      
      Always returns an array of 7 parameters: 
        y0, x01, A1, sigma1, x02, A2, sigma2
      If a single fit is performed the last 3 is set to 0, 0, 1.
      
      Input:
        arr = The array to be fitted
        double = True or False (default is False)
        try_double = True or False (default is False)
        
      Output:
        The gaussian parameters computed from the data array
        """
  if (double == True):
    # Find distinct peaks
    peaks = find_peaks(arr)
    
    # If no peaks are found return 0-gauss
    if ((peaks[0] == -1) and (peaks[1] == -1)):
      params = [0, 0, 0, 1, 0, 0, 1]
      params_error_red = [1, 1, 1, 1, 1, 1, 1]
      return params, params_error_red
    
    # If one peak is found, use single gauss fit
    if ((peaks[0] != -1) and (peaks[1] == -1)):
      if (try_double==False):
        params, params_error_red = gauss_from_array(arr, double=False)
        return params, params_error_red
      else:
        peaks[1] = peaks[0]
    
    # Else if two peaks are found, set one gauss in each peak
    if (arr.min()<0):
      minshift = arr.min()
      arr=arr-minshift
    else:
      minshift=0
    
    y0=arr.min()
    
    x01=peaks[0]
    A1 = arr[x01]-y0
    sigma1 = 3.0
    
    x02=peaks[1]
    A2 = arr[x02]-y0
    sigma2 = 3.0

  else:
    # Set initial values for single fit
    if (arr.min()<0):
      minshift = arr.min()
      arr=arr-minshift
    else:
      minshift=0
    
    y0=arr.min()
    
    x01 = arr.argmax()
    A1 = arr[x01]-y0
    sigma1 = 3.0

  #print str([y0, x01, A1, sigma1, x02, A2, sigma2])
  
  # Shift the spectra so the strongest peak is in the middle of the array
  # Maybe we can remove this as it doesn't affect the fit at all?
  if (do_shift == True):
    if (double==True):
      shift_i = int(arr.size()/2-peaks[0])
    else:
      shift_i = int(arr.size()/2-arr.argmax())
    
    arr = P.shift(arr, shift_i)
    x01 = x01 + shift_i
    if (double==True):
      x02 = x02 + shift_i

  # Set limits on the parameters
  if (double==True):
    number_of_pars = 7
  else:
    number_of_pars = 4
  
  parinfo = []
  for i in range(number_of_pars):
    parinfo.append({'value':0.0, 'fixed':0, 'limited':[0,0],'limits':[0.0, 0.0], 'step':0.0})
  
  parinfo[0]['value'] = y0
  parinfo[0]['limited'][0] = 1
  parinfo[0]['limits'][0] = 0.0 #float(N.minimum(0, arr.min()))
  parinfo[0]['limited'][1] = 1

  arr_sorted = N.sort(arr)
  arr_sorted_mean = arr_sorted[0:5].mean()
  arr_sorted_stddev = arr_sorted[0:5].stddev()
  y0_upper_limit = arr_sorted_mean + arr_sorted_stddev

  parinfo[0]['limits'][1] = arr.min() + 0.5*(arr.max()-arr.min()) #N.maximum(0, y0_upper_limit)
  
  parinfo[1]['value'] = x01
  parinfo[1]['limited'][0] = 1
  parinfo[1]['limits'][0] = 0.0 #-float(arr.size()-1)
  parinfo[1]['limited'][1] = 1
  parinfo[1]['limits'][1] = float(arr.size()-1) #2 * float(arr.size()-1)
  
  parinfo[2]['value'] = A1
  parinfo[2]['limited'][0] = 1
  parinfo[2]['limits'][0] = 0.0
  parinfo[2]['limited'][1] = 1
  parinfo[2]['limits'][1] = float(1.5*arr.max())
  
  parinfo[3]['value'] = sigma1
  parinfo[3]['limited'][0] = 1
  parinfo[3]['limits'][0] = 0.5
  parinfo[3]['limited'][1] = 1
  parinfo[3]['limits'][1] = float(arr.size())/2.0
  
  if (double==True):
    # If this is a double fit, set initial values for the second gauss
    parinfo[4]['value'] = x02
    parinfo[4]['limited'][0] = 1
    parinfo[4]['limits'][0] = 0.0 #-float(arr.size()-1)
    parinfo[4]['limited'][1] = 1
    parinfo[4]['limits'][1] = float(arr.size()-1) #2 * float(arr.size()-1)
    
    parinfo[5]['value'] = A2
    parinfo[5]['limited'][0] = 1
    parinfo[5]['limits'][0] = 0.0
    parinfo[5]['limited'][1] = 1
    parinfo[5]['limits'][1] = float(1.5*arr.max())
  
    parinfo[6]['value'] = sigma2
    parinfo[6]['limited'][0] = 1
    parinfo[6]['limits'][0] = 0.5
    parinfo[6]['limited'][1] = 1
    parinfo[6]['limits'][1] = float(arr.size())/2.0
  
  # Setup arrays with the data to fit
  x = N.arange(float(arr.size())).astype('Float64')
  y = N.zeros(arr.shape, 'Float32')
  y[:] = arr[:]
  err = N.ones(arr.shape, 'Float32')
  
  functkw = {'x':x, 'y':y, 'err':err}
  
  # Fit parameters
  no_converge = False
  if (double==True):
    m = mpfit(gauss_func_double, parinfo=parinfo, functkw=functkw, quiet=1)
    params = m.params
    if (m.perror == None):
      no_converge = True
      params_error = N.array([0,0,0,0,0,0,0])
    else:
      params_error = m.perror
  else:
    m = mpfit(gauss_func_single, parinfo=parinfo, functkw=functkw, quiet=1)
    params = N.concatenate((m.params,[0,0,1]))
    if (m.perror == None):
      params = N.array([0,0,0,1,0,0,1])
      params_error = N.array([0,0,0,0,0,0,0])
    else:
      params_error = N.concatenate((m.perror,[0,0,0]))
  
  # Reduced chi-square error
  dof = len(x) - len(m.params)
  params_error_red = params_error * N.sqrt(m.fnorm / dof)

  # If the peaks were shifted, shift back
  # Maybe we can remove this alltogether?
  if (do_shift == True):
    params[1] = params[1] - shift_i
    if (double==True):
      params[4] = params[4] - shift_i
  
  # Set the highest peaks parameters to the left
  if (params[2] < params[5]):
    params = N.concatenate((params[0], params[4:7], params[1:4]))
  
  # If one peak is very weak compared to the other one, 
  # or if the peaks are located close to each other (within half the sigma
  # of the stronger peak), make a single gauss fit.
  if (double == True):
    # Set limits to use for QA of the fit
    limit_err = 500
    limit_amp = 0.05
    limit_width = 1.0
    
    if (no_converge == True):
      params, params_error_red = gauss_from_array(arr, double=False)
    elif (params[5] < limit_amp*params[2]):
      params, params_error_red = gauss_from_array(arr, double=False)
    elif (N.abs(params[1]-params[4]) < N.maximum(limit_width*params[3], limit_width*params[6])):
      params, params_error_red = gauss_from_array(arr, double=False)
    else:
      # Compute the sum of squared errors
      temp=(err*(gauss_to_array(params, len(y))-y))**2
      
      # If the error is big, do a single gauss fit
      if (temp.sum() > limit_err):
        params, params_error_red = gauss_from_array(arr, double=False)
  
  params[0]+=minshift
  return params, params_error_red

def gauss_func_single(p, fjac=None, x=None, y=None, err=None):
  """This is used by gauss_from_array(). It returnes a status flag and
      the weighted error of a gaussian fit to a data set.
      """
  # Calculate the exponents first
  arr_len = len(x)
  limit = -700
  
  temp1 = N.maximum( -((x-p[1])**2)/(2*p[3]**2) , limit)
  temp1_right = N.maximum( -(((x+arr_len)-p[1])**2)/(2*p[3]**2) , limit)
  temp1_left = N.maximum( -(((x-arr_len)-p[1])**2)/(2*p[3]**2) , limit)
    
  # Compute the double gauss
  arr = p[0] + p[2]*Numeric.exp(temp1)
  arr_left = p[2]*Numeric.exp(temp1_left)
  arr_right = p[2]*Numeric.exp(temp1_right)
  
  model_y = arr + arr_right + arr_left
  # We are not using this at the moment, 0 means successful computation
  status=0
  
  # Return the status and the weighted error in the fit
  return [status, err*(y-model_y)]

def gauss_func_double(p, fjac=None, x=None, y=None, err=None):
  """This is used by gauss_from_array(). It returnes a status flag and
      the weighted error of a gaussian fit to a data set.
      """
  # Calculate the exponents first
  arr_len = len(x)
  limit = -700
  
  temp1 = N.maximum( -((x-p[1])**2)/(2*p[3]**2) , limit)
  temp2 = N.maximum( -((x-p[4])**2)/(2*p[6]**2) , limit)
  temp1_right = N.maximum( -(((x+arr_len)-p[1])**2)/(2*p[3]**2) , limit)
  temp2_right = N.maximum( -(((x+arr_len)-p[4])**2)/(2*p[6]**2) , limit)
  temp1_left = N.maximum( -(((x-arr_len)-p[1])**2)/(2*p[3]**2) , limit)
  temp2_left = N.maximum( -(((x-arr_len)-p[4])**2)/(2*p[6]**2) , limit)
    
  # Compute the double gauss
  arr = p[0] + p[2]*Numeric.exp(temp1) + p[5]*Numeric.exp(temp2)
  arr_left = p[2]*Numeric.exp(temp1_left) + p[5]*Numeric.exp(temp2_left)
  arr_right = p[2]*Numeric.exp(temp1_right) + p[5]*Numeric.exp(temp2_right)
  
  model_y = arr + arr_right + arr_left
  # We are not using this at the moment, 0 means successful computation
  status=0
  
  # Return the status and the weighted error in the fit
  return [status, err*(y-model_y)]

  
def gauss_slit(data, slitname='1', outfile=None, limit=5):
  """Interactive gauss fits along a slit.
      
      Input:
      data = The 3D adhoc-object to be fitted
      slitname = The slit-file to be used
      outfile = If given, the result is written to a rc-file with the name 
        given, default is None.
      limit = The minimum peak-height to consider useful data to fit.
      
      Output:
      pos,vel = Two lists with position and velocity for all points.
      """
  
  # Read slit data
  pars=read_slit(slitname)
  angle=pars['angle']
  center=data.p['cen'] + pars['offset']
  d = pars['slitwidth']/data.p['echelle']
  
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

  # Use interactive gauss fit on all points within the slit
  for i in N.arange(data.nx()):
    for j in N.arange(data.ny()):
      # The peak must be large enough
      if (data[i,j].max() > limit):
        # The points must be within the slit
        dist=N.innerproduct(perp,N.array([center[0]-i,center[1]-j]))
        if (abs(dist) <= D):
          # Gauss fit
          pos_temp=N.innerproduct(vec,N.array([center[0]-i,center[1]-j]))
          par_temp=gauss_from_array_interactive(data[i,j])
          
          # Check if any peaks are found, if so, add them to the list
          if (par_temp[0][1] != 0):
            vel_temp=par_temp[0][1]/data.p['lz']*data.fsr() + lamb2vel(data.p['xlbneb']-0.5*data.p['xil']) + data.p['vr_offset']
            vel=N.concatenate((vel, vel_temp))
            pos=N.concatenate((pos, pos_temp))
            err=N.concatenate((err, data.fsr()/data.p['lz']))
          
          if (par_temp[0][4] != 0):
            vel_temp=par_temp[0][4]/data.p['lz']*data.fsr() + lamb2vel(data.p['xlbneb']-0.5*data.p['xil']) + data.p['vr_offset']
            vel=N.concatenate((vel, vel_temp))
            pos=N.concatenate((pos, pos_temp))
            err=N.concatenate((err, data.fsr()/data.p['lz']))

  # Write rc-outfile if a filename is given
  if (outfile != None):
    write_rc(outfile, pos, vel, err)

  return pos,vel
