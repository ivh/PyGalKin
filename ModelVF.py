"""
 ModelVF.py

 For a general description of the whole package and how to import,
 see the head of the file PyCigale.py
 
 Usage:
 PyC.model('disk','output.fits',256,100,10,60,70)

 Models implemented so far:
 * linear: gives a velocity field with a simple linear gradient
 * shell: velocity field of a shell that is beeing looked at
   ( it turned out to be linear as well :-))
 * disk: a thin disk in DM-potential

 See the function descriptions as well.

 This file needs better integration into the package PyCigale!!
"""

import numarray as N
import Scientific.Functions.LeastSquares as LS
import math as M
import os
import InOutput as IO
import PyCigale as PyCigale
import matplotlib
import matplotlib.pylab as MP
#matplotlib.use('GTK')
import mpfit

###
import pickle
import pylab as P
import PyGalKin as G
from matplotlib.numerix.ma import masked_where 
###




class interactvf:
	def __init__(self, data, pixels=7, avg_square=3):
	
		# Basic Values
		self.avg_square_size = avg_square
		self.pixels = pixels
		self.spectra_size = (0.15,0.10)
		self.spectra_dist = (0.01,0.01)
		self.spectra_coord = (0.015,0.44)
		self.spectra_break = 4
		self.cenx = 0
		self.ceny = 0
			
		self.trusted = []
		self.trustgraph = []
#		self.trustlength = 0

		self.vel_pars = None
		self.flux_pars = None
		self.spectral_pars = None
		
		# import and process data cube
		self.odata = data.get_copy()
		self.cdata = data.cliparoundcenter()
		self.cdata = N.transpose(self.cdata,axes=(1,0,2))
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
#		P.title(self.odata.p['objname'])
		canvas = self.fig.canvas
		canvas.mpl_connect('key_press_event', self.key_press_callback)
		canvas.mpl_connect('button_press_event', self.button_press_callback)
		self.canvas = canvas
		self.plot()
		
        
	def plot(self):
		self.fig.clf()
        
        # Total Flux
		self.axflux=P.axes([-0.03,0.59,0.36,0.36])
		minmaxflux = N.array(self.sdata)
		minmaxflux.shape = (-1,)
		minmaxflux = N.sort(minmaxflux)
		takeflux = int(round(len(minmaxflux)*0.01))
		axfluxresmax = max(int(round(minmaxflux[-takeflux])),1)
		axfluxresmin = max(int(round(minmaxflux[takeflux])), 0)
		P.setp(self.axflux,xticks=[], yticks=[])
		P.imshow(N.transpose(self.sdata),interpolation='nearest',vmin=axfluxresmin,vmax=axfluxresmax)
		P.title('Total Flux')
		
		# Velocity Field
		self.axvelf=P.axes([0.27,0.59,0.36,0.36])
		minmaxvel = N.array(self.ini_velf)
		minmaxvel.shape = (-1,)
		minmaxvel = N.sort(minmaxvel)
		takevel = int(round(len(minmaxvel)*0.01))
		axvelresmax = max(int(round(minmaxvel[-takevel])),1)
		axvelresmin = max(int(round(minmaxvel[takevel])), 0)
		P.setp(self.axvelf,xticks=[], yticks=[])
		P.imshow(N.transpose(self.ini_velf), interpolation='nearest',vmin=axvelresmin, vmax=axvelresmax)
		P.title('Velocity Field')
		
		# Local Spectrum
		self.axspec=P.axes([0.015,0.44,0.15,0.10],xticks=[], yticks=[])
#		P.title('Local Spectrum')
		
		
	def key_press_callback(self,event):
        # key bindings
		if event.key == 'q': self.fitvellinear()
		elif event.key == 'v': self.del_all_trusted()
		elif event.key == 'a': self.fitfluxexp(modeltype='expfit_func')
		elif event.key == 's': self.fitfluxexp(modeltype='powerlawfit_func')
		elif event.key == 'p': self.printfit()
		elif event.key == 'y': self.modelspectra_from_trusted()
		elif event.key == 'b': self.get_modelled_spectrum()
		elif event.key == 'control': P.close(1); print 'interactvf closed. goodbye!'
#		elif event.key == 't': print event.x, event.y, event.inaxes
		else: print "Unknown key pressed:", event.key, '@ ('+`event.x`+','+`event.y`+ ") ... doing nothing"
#		self.plot()


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
		
		
	def plotat(self, xy):
		# refresh local spectrum
#		P.delaxes(self.axspec)
#		self.axspec=P.axes([0.62,0.59,0.36,0.36])
		self.axspec.cla()
		self.axspec.plot(self.cdata[xy[0],xy[1],:], 'r', linewidth=2)
		P.setp(self.axspec, xticks=[], yticks=[], title='Local Spectrum')
#		P.title(self.axspec, 'Local Spectrum')
		self.currdat = self.cdata[xy[0],xy[1],:]
		self.canvas.draw()
		
		
	def add_to_trusted(self, xy):
		# add coordinates, spectrum and velocity to trusted data
		temp = N.zeros((self.cdata.shape[2],))
		cnt = 0.0
		dev = self.avg_square_size / 2
		for i in range(self.avg_square_size):
			for j in range(self.avg_square_size):
				try:
					tmp = self.cdata[xy[0]+i-dev,xy[1]+j-dev,:]
#					self.move_to_peak(tmp)
					temp += tmp
				except: print 'Point at the border chosen! Continuing anyways.'
				else: cnt += 1.0
		temp = temp/cnt
		veli = G.calcpeak(temp, 7)
		newlistitem = [temp, xy, self.ini_velf[xy[0],xy[1]], self.sdata[xy[0],xy[1]]]
		self.trusted.append(newlistitem)
		print 'added point @', xy, 'to trusted points'
		self.update_trusted()
		
		
	def drop_trusted(self, axis):
		# remove trusted data
		i = self.trustgraph.index(axis)
		print 'point @', self.trusted[i][1], 'being removed'
		self.trusted = self.trusted[:i] + self.trusted[i+1:]
		self.update_trusted(delf=i)
		
		
	def del_all_trusted(self):
		# remove all trusted data-points
		for i in range(len(self.trusted)):
			self.trusted = self.trusted[:-1]
			self.update_trusted(delf=i)
			self.canvas.draw()
		print 'all trusted points removed'
			
		
	def update_trusted(self, delf=-1):
		# update graphs of trusted data				
		trustlen = len(self.trusted)
		
		# add new graphs		
		if len(self.trustgraph) == trustlen - 1:
			self.trustgraph.append(P.axes(self.newaxis(trustlen-1, cutat = self.spectra_break), xticks=[],yticks=[]))
#			tmplen = len(self.trustgraph)-1
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
		all_coordy = N.array(range(-self.cdata.shape[1]/2-self.ceny,self.cdata.shape[1]/2-self.ceny)*self.cdata.shape[0])
		all_coordy.sort()
		self.field = linfit_func(linfitres.params, None,  x=all_coordx, y=all_coordy, z=None, err=None, returnmodel=True)
		self.field.shape = ((self.data_shape[0],self.data_shape[1]))
		self.axlinf=P.axes([0.61,0.03,0.36,0.36])
		P.setp(self.axlinf,xticks=[], yticks=[])
		P.imshow(self.field, interpolation='nearest')
		P.title('Velocity Model')
		P.colorbar()
		
		print 'done!\n'
		
		# redrawings
		self.canvas.draw()
		
			
	def fitfluxexp(self, modeltype, contnr=3):
		# fit exponential law to flux of trusted points
		
		# throw away your television
		try: P.delaxes(self.axfluxmod)
		except: pass
		else: print 'old axis killed'


		print 'doing:', modeltype, '...'
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
			fluxes = N.concatenate((fluxes, flux))
			x = (i[1][0]-self.cdata.shape[0]/2-self.cenx)
			y = (i[1][1]-self.cdata.shape[1]/2-self.ceny)
			r = M.sqrt(x**2+y**2)
			radii.append(r)
		radii = N.array(radii)
		fluxes = N.log(fluxes)																		####################
		
		ini_coordx = N.array(ini_coordx)
		ini_coordy = N.array(ini_coordy)
		
		
		# initialize mpfit
		parinfo=[]
		parinfo.append({'value':10.0, 'fixed':0, 'limited':[1,0],'limits':[0.0, 0.0], 'step':0.0, 'parname':'max'})		# position angle
		parinfo.append({'value':-0.1, 'fixed':0, 'limited':[0,0],'limits':[0.0, 0.0], 'step':0.0, 'parname':'slope'})	# gradient

		err = N.sqrt(fluxes)
		functkw = {'x':ini_coordx, 'y':ini_coordy, 'z':fluxes, 'err':err}

		expfitres=mpfit.mpfit(modeltype,parinfo=parinfo,functkw=functkw,quiet=False)
		
		# store data
		print expfitres.params
		print 'radii:',radii
		self.flux_pars = {'model':modelt, 'amp':expfitres.params[0], 'slope':expfitres.params[1]}
		
		# nice output
		self.axfluxmod=P.axes([0.62,0.75,0.36,0.20], yticks=[], title='Flux Model')
		self.axfluxmod.plot(radii, M.e**fluxes,'bo')
		self.axfluxmod.semilogy()
		
		# get fitted field
		maxrad = max(radii)
		modelled_flux = N.array([])
		l = N.array(range(1,int(round(maxrad))+1))
		modelled_flux = M.e **(modeltype(expfitres.params, None,  x=l, y=N.zeros((len(l),)), z=None, err=None, returnmodel=True))
		self.axfluxmod.plot(l, modelled_flux, 'k')

		print 'done!\n'

		# redrawings
		self.canvas.draw()
			
			
			
	def modelspectra_from_trusted(self):
		# get typical spectral shape from trusted points

		# die,die,die!
		try: P.delaxes(self.axspecmod)
		except: pass
		else: print 'old axis killed'

		# construct model spectrum
		try:
			vel_base = self.vel_pars['system']
		except TypeError:
			vel_base = self.odata.vel1st()
			print 'taking velocity from par'
		else:
			print 'got fitted systemic velocity'
			
		vel_step = self.cdata.fsr()
		channr = self.cdata.shape[2]
		midchan = channr/2
		sspectra = N.zeros(channr)
		sspectra = N.array(sspectra, type=N.Float32)
		for i in self.trusted:
			vel = i[2]
			chan = int(round((i[2]-vel_base)/vel_step))
			sspec = G.PyCigale.shift(i[0], chan)
			sspectra += sspec
		norm = N.sum(sspectra)
		sspectra /= norm
		
		# store fit in plot
		self.spectral_pars = {'model':'empirical', 'spectrum':sspectra}
		print 'model spectrum done\n'
		
		# plot it
		self.axspecmod=P.axes([0.62,0.47,0.36,0.20], xticks=[], yticks=[], title='Model Spectrum')
		self.axspecmod.plot(self.spectral_pars['spectrum'], 'r', linewidth=2)

		# redrawings
		self.canvas.draw()
	
	
	def get_modelled_spectrum(self):
		# use all of the results from fitting procedures to form model spectrum
		
		# check for performed fits, else abort
		if self.vel_pars == None:
			print '!!! No velocity model yet. Please perform modelling first (q) !!!'
			return
		elif self.flux_pars == None:
			print '!!! No flux model yet. Please perform modelling first (a,s) !!!'
			return
		elif self.spectral_pars == None:
			print '!!! No spectral model yet. Please perform modelling first (y) !!!'
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
			flux_pars = [self.flux_pars['amp'], self.flux_pars['slope']]
		
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
		all_coordy = all_coordx.copy()
		all_coordy.sort()
		
		# get velocity model
		velmod = vel_model(vel_pars, None,  x=all_coordx, y=all_coordy, z=None, err=None, returnmodel=True)
		velmod.shape = (2*wdt,2*hgt)
		
		# get flux model
		fluxmod = N.exp(flux_model(flux_pars, None,  x=all_coordx, y=all_coordy, z=None, err=None, returnmodel=True))
		fluxmod.shape = (2*wdt,2*hgt)
		self.fluxmod = fluxmod
		
		# get unshifted spectral model
		intermediate_model = fluxmod[:,:,N.NewAxis] * model
		
		# get indexshifts
		channr = self.cdata.shape[2]
		midchan = channr/2
		vel_base = self.vel_pars['system']
		vel_step = self.cdata.fsr()/channr
		shiftarray = N.around((velmod-vel_base)/vel_step)
		shiftarray.type = 'Int32'
		
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
#		self.synth_total = G.sum(self.synth_spectrum)		
#		self.synth_vel = G.peakvel(self.synth_spectrum, 7)
		self.residual = self.cdata - self.synth_spectrum
		self.residual_sum = G.sum(self.residual)
		self.residual_vel = G.peakvel(self.residual, 7)
		
		# show it!!
#		self.result = P.figure()
#		self.rescanvas = self.result.canvas

		# Total Flux
		try: 							################
			P.delaxes(self.axflux)		################
			P.delaxes(self.axfluxres)	################
		except: pass					################
		minmaxflux = N.array(self.residual_sum)
		minmaxflux.shape = (-1,)
		minmaxflux = N.sort(minmaxflux)
		takeflux = int(round(len(minmaxflux)*0.01))
		axfluxresmax = max(int(round(minmaxflux[-takeflux])),1)
		axfluxresmin = max(int(round(minmaxflux[takeflux])), 0)
		self.axfluxres=P.axes([-0.03,0.59,0.36,0.36], xticks=[], yticks=[],title='Total Flux (Residual)')
		self.axfluxres.imshow(self.residual_sum, interpolation='nearest', vmin=axfluxresmin,vmax=axfluxresmax)
		
		# Velocity Field
		try: 							################
			P.delaxes(self.axvelf)		################
			P.delaxes(self.axvelf1)		################
		except: pass					################
		minmaxvel = N.array(self.residual_vel)
		minmaxvel.shape = (-1,)
		minmaxvel = N.sort(minmaxvel)
		takevel = int(round(len(minmaxvel)*0.01))
		axvelresmax = max(int(round(minmaxvel[-takevel])),1)
		axvelresmin = max(int(round(minmaxvel[takevel])), 0)
		self.axvelf1=P.axes([0.29,0.59,0.36,0.36],xticks=[], yticks=[],title='Velocity Field (Residual)')
		self.axvelf1.imshow(self.residual_vel, interpolation='nearest', vmin=axvelresmin,vmax=axvelresmax)

		self.canvas.draw()
		
		print 'done!\n'
		
		
	def set_center(self,xy):
		# select new center for the galaxy
		self.cenx = xy[0] - self.cdata.shape[0]/2
		self.ceny = xy[1] - self.cdata.shape[1]/2
		print 'new center selected @', xy
		

	def printfit(self):
		print '\nParameters for velocity field model:\n', self.vel_pars
		print '\nParameters for radial flux model:\n', self.flux_pars
		print '\nParameters for spectral model:\n', self.spectral_pars, '\n'
	
	
	def coord(self, event):
		# extract proper coordinates of mouse event
		x = int(round(event.xdata))
		y = int(round(event.ydata))
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
			if sx >= 0:
				ang = M.pi/2
			else:
				ang = -M.pi/2
		else:
			ang = M.atan(sx/sy)

		t = sx * M.sin(-pa) + sy * M.cos(pa)
		o = sx * M.cos(pa) + sy * M.sin(pa)
		dx = cenx + t * M.sin(-pa) - x[i]
		dy = cenx + t * M.cos(pa) - y[i]
		if o >= 0:
			sgn = -1
		else:
			sgn = +1
		
		r = sgn*M.sqrt(dx**2 + dy**2)
		radii.append(r)
		
	radii = N.array(radii)
	return radii


def interact_model_linear(r, pars):
	# calculate velocities from gradient and radius
	arr = r[:].copy()
	velgrad = pars['gradient']
	sysvel = pars['system']
	arr = arr * velgrad + sysvel
	return arr


def galmod_fitlin(x, y, p):
	# returns linear vel field from coordinates x,y and fitting parameters p
	pars = {}
	pars['pa'] = p[0]
	pars['gradient'] = p[1]
	pars['system'] = p[2]
	pars['centr_offset_x'] = 0
	pars['centr_offset_y'] = 0
	
	radii = radii_from_position((x,y),pars)
	vels = interact_model_linear(radii, pars)
	return vels
	


def linfit_func(p, fjac, x=None, y=None, z=None, err=None, returnmodel=False):
	# linear galaxy rotation modelling, conform to mpfit
	model = galmod_fitlin(x, y, p)
	if returnmodel: return model
	else: return([0, z-model])
		
		
def expfit_func(p, fjac, x=None, y=None, z=None, err=None, returnmodel=False):						####################
	# exponential galaxy flux modelling, conform to mpfit
	model = p[0]+(p[1]*N.sqrt(x**2+y**2))															####################
	if returnmodel: return model
	else: return([0, z-model])
		
		
def powerlawfit_func(p, fjac, x=None, y=None, z=None, err=None, returnmodel=False):		
	# exponential galaxy flux modelling, conform to mpfit
	model = p[0]+p[1]*N.log(N.sqrt(x**2+y**2))
	if returnmodel: return model
	else: return([0, z-model])
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
##### trashy stuff ######
##	
	def move_to_peak(selfmo, row):
		
		parinfo=[]
		parinfo.append({'value':1.0, 'fixed':0, 'limited':[1,0],'limits':[0.0, 0.0], 'step':0.0})
		parinfo.append({'value':1.0, 'fixed':0, 'limited':[1,0],'limits':[0.0, 0.0], 'step':0.0})
		parinfo.append({'value':24.0, 'fixed':0, 'limited':[1,1],'limits':[0.0, len(row)], 'step':0.0})
		
		x = range(len(row))
		y = row
		err = 1/abs(sqrt(y))
		functkw = {'x':x, 'y':y, 'err':err}
		
		gaussres=mpfit.mpfit(gauss_fit_func,parinfo=parinfo,functkw=functkw,quiet=False)
		print gaussres.status
		print gaussres.params
		
		## + a lot more			
	
def gauss_fit_func(p, fjac, x=None, y=None, err=None):
	# gaussian fitting function, conform to mpfit
	model = p[0]*exp(-p[1]*(x-P[2])**2)
	return([0, y-model])
##
##### end trashy stuff #####


#################################################################


def rotation_curve(rings, outfile=None):
  """Generates a position and velocity vector from a list of rings generated
      by tilted_ring_model(). It can also write the output to a .rc-file.
      
      Usage: pos,vel = rotation_curve(rings, outfile)
      
      rings: A list of rings generated by tilted_ring_model()
      outfile: Optional filename to write the data to outfile.rc (default is None)
      pos,vel: Vectors containing position and velocity along a rotation curve.
      """
  # Create empty lists
  vel=[]
  pos=[]
  err=[]
  
  # For each ring, read the radius (scaled to pc's) and velocity
  for i in range(len(rings)):
    pos.append(float(rings[i]['radius']*rings[i]['scale']))
    vel.appens(float(rings[i]['vel']))
    err.append(float(0.0))
  
  # Write the outfile if a name is given
  if (outfile != None):
    write_rc(outfile, pos, vel, err)
  
  return pos,vel
  
def tilted_ring_model(data, width):
  """Creates a tilted ring-model from a velocity field.
      
      Usage: system,rings = tilted_ring_model(data, width)
      
      data: The velocity field to model
      width: The width of each ring, this will set the number of rings
      system: The system velocity
      rings: A list of dictionaries (one for each ring) with parameters on the 
      form:
        [{'centr_x', 'centr_y', 'vel', 'inclination', 'pa', 'radius', 'width', 
        'dim', 'scale'},
        {'centr_x', 'centr_y', 'vel', 'inclination', 'pa', 'radius', 'width', 
        'dim', 'scale'}, ... ]
        
        centr_x,centr_y: The center of the ring
        vel: The velocity of the ring
        inclination: The inclination of the ring
        pa: The position angle of the ring
        radius: The inner radius of the ring
        width: The width of the ring
        dim: The dimension of the array for the ring (always a square array)
        scale: The parsec/pixel-scale
      """
  # We do not want to edit a link to the data
  temp = data.get_copy()
  
  # Understandable variables
  dim = temp.nx()
  centr_x = temp.p['dyncen'][0]
  centr_y = temp.p['dyncen'][1]
  vel = 50
  inclination = 45
  pa = temp.p['pa']
  minimum = temp.min()
  
  # Fit the system velocity, only one parameter
  parinfo = []
  parinfo.append({'value':0.0, 'fixed':0, 'limited':[0,0],'limits':[0.0, 0.0], 'step':0.0})
  parinfo[0]['value'] = temp.mean()
  
  functkw = {'x':temp, 'y':dim, 'err':minimum}
  temp1 = mpfit.mpfit(model_system_func, parinfo=parinfo, functkw=functkw, quiet=1)
  v_system = temp1.params[0]
  
  print v_system
  
  # Empty list for the rings
  pars = []
  
  # Compute the rings
  for i in range(int(0.5*dim/width)):
    radius = i*width
    print 'Ring: '+str(i+1)+', radius: '+str(radius)
    
    # Parameters to fit
    parinfo = []
    for j in range(5):
      parinfo.append({'value':0.0, 'fixed':0, 'limited':[0,0],'limits':[0.0, 0.0], 'step':0.0})
    
    parinfo[0]['value'] = pa
    parinfo[0]['limited'][0] = 0
    parinfo[0]['limits'][0] = 0
    parinfo[0]['limited'][1] = 0
    parinfo[0]['limits'][1] = 0
    parinfo[0]['fixed'] = 0
    
    parinfo[1]['value'] = vel
    parinfo[1]['limited'][0] = 0
    parinfo[1]['limits'][0] = 0
    parinfo[1]['limited'][1] = 0
    parinfo[1]['limits'][1] = 0
    parinfo[1]['fixed'] = 0
    
    parinfo[2]['value'] = inclination
    parinfo[2]['limited'][0] = 0
    parinfo[2]['limits'][0] = 0
    parinfo[2]['limited'][1] = 0
    parinfo[2]['limits'][1] = 0
    parinfo[2]['fixed'] = 0
    
    parinfo[3]['value'] = centr_x
    parinfo[3]['limited'][0] = 0
    parinfo[3]['limits'][0] = 0
    parinfo[3]['limited'][1] = 0
    parinfo[3]['limits'][1] = 0
    parinfo[3]['fixed'] = 0
    
    parinfo[4]['value'] = centr_y
    parinfo[4]['limited'][0] = 0
    parinfo[4]['limits'][0] = 0
    parinfo[4]['limited'][1] = 0
    parinfo[4]['limits'][1] = 0
    parinfo[4]['fixed'] = 0
    
    functkw = {'x':temp, 'y':[radius, width, dim, v_system, minimum], 'err':0}
    
    temp2 = mpfit.mpfit(model_tilted_ring_func, parinfo=parinfo, functkw=functkw, quiet=1)
    p = temp2.params
    
    # Add the ring to pars
    pars.append({'centr_x':p[3], 'centr_y':p[4], 'vel':p[1], 'inclination':p[2], 'pa':p[0], 'radius':radius, 'width':width, 'dim':dim, 'scale':temp.scale()})
    
    print str(pars[i])
  
  return [v_system, pars]

def model_system_func(p, fjac=None, x=None, y=None, err=None): 
  """Used by tilted_ring_model() to fit a system velocity to a velocity field.
      It will only consider points containing data, not empty points.
      
      Note: This is not supposed to be used from the command line.
      """
  # Set meaningful variable names
  data = x
  dim = y
  minimum = err
  
  # Array to store the result
  result = N.zeros((dim,dim), type='Float64')
  # Uniform velocity field
  vf_system = N.zeros((dim,dim), type='Float64') + p[0]
  # Points where there is data
  points = N.where(data > minimum)
  # Difference between model and data
  result[points] = vf_system[points] - data[points]
  result.setshape(dim*dim)
  
  status=0
  return [status, result]
  
def model_tilted_ring_func(p, fjac=None, x=None, y=None, err=None):
  """Used by tilted_ring_model() to fit a tilted ring to a velocity field.
      It will only consider points containing data (both in the velocity field 
      and the ring), not empty points.
      
      Note: This is not supposed to be used from the command line.
      """
  # Setup variables with understandable names
  data = x
  radius = y[0]
  width = y[1]
  dim = y[2]
  v_system = y[3]
  minimum = y[4]
  
  # Create a system velocity field from the known system velocity
  vf_system = N.zeros((dim,dim), type='Float64') + v_system

  # Create a ring and add the system velocity
  ring = create_ring(p, radius, width, dim)
  
  # Only consider points where the ring has some data and
  # where the 'data'  has data
  points1 = N.where(ring != 0)
  points2 = N.where(data > minimum)
  
  ring[points1] = ring[points1] + vf_system[points1]
  
  # The difference between the data and the model
  result_temp = N.zeros((dim,dim), type='Float64')
  result_temp[points1] = ring[points1] - data[points1]
  
  result = N.zeros((dim,dim), type='Float64')
  result[points2] = result_temp[points2]
  result.setshape(dim*dim)
  
  print str(p)
  
  status=0
  return [status, result]

def create_ring(p, radius, width, dim):
  """Creates a tilted ring as viewed from the telescope.
      
      Usage: ring_tr = create_ring(p, radius, width, dim)
      
      p: A list of parameters for the ring [pa,vel,inclination,centr_x,centr_y]
        centr_x, centr_y: The centre of the ring
        vel: The velocity
        inclination: The inclination of the ring
        pa: The position angle of the ring
      radius: The inner radius of the ring
      width: The width of the ring
      dim: The dimension of the array (the array is always a square)
      ring_tr: A 2D-array with the tilted ring
      """
  # Set understandable variable names
  centr_x = p[3]
  centr_y = p[4]
  vel = p[1]
  inclination = p[2]
  pa = p[0]
  
  # Convert the pa to rads and transform it to the correct definition
  pa = (-pa-90)*M.pi/180
  # Convert the inclination to rads
  inclination = inclination*M.pi/180

  # Arrays needed for the computation
  r = N.zeros((dim,dim), type='Float64')
  phi = N.zeros((dim,dim), type='Float64')
  ring_non_tr = N.zeros((dim,dim), type='Float64')
  ring_tr = N.zeros((dim,dim), type='Float64')
  
  # Create an array with the angle and one with the radius in each point
  # calculated from the ring-centre.
  for x in range(dim):
    for y in range(dim):
      temp_x = x-centr_x
      temp_y = y-centr_y
      
      phi[x,y] = N.arctan2(temp_y, temp_x)
      r[x,y] = N.sqrt(temp_x**2+temp_y**2)

  # Mark points within the inner ring as t2 and points within the outer ring as t1
  t1 = N.where(r < (radius+width))
  t2 = N.where(r < radius)
  
  # Create a non transformed ring
  ring_non_tr[t1] = vel*N.cos(phi[t1])*N.sin(inclination)
  ring_non_tr[t2] = 0.0
  
  # Transform the ring to the telescope-coordinates (from the galaxy coords.)
  for x in range(dim):
    for y in range(dim):
      temp_x = x-centr_x
      temp_y = y-centr_y
      
      x_corr = temp_x*N.cos(pa) - temp_y*N.sin(pa)
      y_corr = temp_x*N.sin(pa) + temp_y*N.cos(pa)
      y_corr = y_corr/N.cos(inclination)
      
      x_corr2 = int(x_corr+centr_x)
      y_corr2 = int(y_corr+centr_y)
      
      if (x_corr2 >= 0 and x_corr2 < dim) and (y_corr2 >= 0 and y_corr2 < dim):
        ring_tr[x,y] = ring_non_tr[x_corr2, y_corr2]
  
  return ring_tr


def fit_parameters(data, model_list, box_size=3):
  """Fits a model defined in model_list to the data. You interactivly choose
      boxes of size 'box_size' to use for comparing the model and the data.
      
      An average value will be calculated for each box you choose and
      compared with the average from the same box in the model. This
      difference is minimized.
      
      The output is a model_list with the best fit.
      
      Usage: model_list_out = fit_parameters(data, model_list_in, box_size)
      
      data: The velocity field
      model_list_in: A list of models to fit and the initial values
      box_size: The size of the boxes to average (a kind of low pass filter)
      model_list_out: A list of the models and the best fit parameters
      """
  # Ask the user for boxes
  boxes = get_boxes(data)
  # Average the boxes
  boxes_avg = get_boxes_avg(data, boxes, box_size)
  
  # Fixes some obscure bug, linking was used instead of copying
  model_list = N.copy.deepcopy(model_list)
  
  # Set limits on the parameters
  parinfo = []

  for i in range( len(model_list)*10 + 1 ):
    parinfo.append({'value':0.0, 'fixed':1, 'limited':[0,0],'limits':[0.0, 0.0], 'step':0.0})
  
  parinfo[0]['value'] = model_list[0][1]['dim']
  
  for i in range(len(model_list)):
    parinfo[i*10+1]['value'] = model_list[i][1]['pa']
    parinfo[i*10+1]['limited'][0] = 1
    parinfo[i*10+1]['limits'][0] = 0
    parinfo[i*10+1]['limited'][1] = 1
    parinfo[i*10+1]['limits'][1] = 360
    parinfo[i*10+1]['fixed'] = 0
    
    parinfo[i*10+2]['value'] = model_list[i][1]['inclination']
    parinfo[i*10+2]['limited'][0] = 1
    parinfo[i*10+2]['limits'][0] = 0
    parinfo[i*10+2]['limited'][1] = 1
    parinfo[i*10+2]['limits'][1] = 90
    parinfo[i*10+2]['fixed'] = 0
    
    parinfo[i*10+3]['value'] = model_list[i][1]['exp_max']
    parinfo[i*10+3]['limited'][0] = 1
    parinfo[i*10+3]['limits'][0] = 0
    parinfo[i*10+3]['limited'][1] = 1
    parinfo[i*10+3]['limits'][1] = model_list[0][1]['dim']
    parinfo[i*10+3]['fixed'] = 1
    
    parinfo[i*10+4]['value'] = model_list[i][1]['r_max']
    parinfo[i*10+4]['limited'][0] = 1
    parinfo[i*10+4]['limits'][0] = 0
    parinfo[i*10+4]['limited'][1] = 1
    parinfo[i*10+4]['limits'][1] = model_list[0][1]['dim']
    parinfo[i*10+4]['fixed'] = 1
    
    parinfo[i*10+5]['value'] = model_list[i][1]['v_max']
    parinfo[i*10+5]['limited'][0] = 1
    parinfo[i*10+5]['limits'][0] = 0
    parinfo[i*10+5]['limited'][1] = 1
    parinfo[i*10+5]['limits'][1] = 1000
    parinfo[i*10+5]['fixed'] = 0
    
    parinfo[i*10+6]['value'] = model_list[i][1]['v_system']
    parinfo[i*10+6]['limited'][0] = 1
    parinfo[i*10+6]['limits'][0] = 0
    parinfo[i*10+6]['limited'][1] = 1
    parinfo[i*10+6]['limits'][1] = 2000
    parinfo[i*10+6]['fixed'] = 0
    
    parinfo[i*10+7]['value'] = model_list[i][1]['v_expansion']
    parinfo[i*10+7]['limited'][0] = 1
    parinfo[i*10+7]['limits'][0] = 0
    parinfo[i*10+7]['limited'][1] = 1
    parinfo[i*10+7]['limits'][1] = 1000
    parinfo[i*10+7]['fixed'] = 1
    
    parinfo[i*10+8]['value'] = model_list[i][1]['a_scale']
    parinfo[i*10+8]['limited'][0] = 1
    parinfo[i*10+8]['limits'][0] = 0
    parinfo[i*10+8]['limited'][1] = 1
    parinfo[i*10+8]['limits'][1] = 100
    parinfo[i*10+8]['fixed'] = 0
    
    parinfo[i*10+9]['value'] = model_list[i][1]['centr_offset_x']
    parinfo[i*10+9]['limited'][0] = 1
    parinfo[i*10+9]['limits'][0] = -int(0.5*model_list[0][1]['dim'])
    parinfo[i*10+9]['limited'][1] = 1
    parinfo[i*10+9]['limits'][1] = int(0.5*model_list[0][1]['dim'])
    parinfo[i*10+9]['fixed'] = 1
    
    parinfo[i*10+10]['value'] = model_list[i][1]['centr_offset_y']
    parinfo[i*10+10]['limited'][0] = 1
    parinfo[i*10+10]['limits'][0] = -int(0.5*model_list[0][1]['dim'])
    parinfo[i*10+10]['limited'][1] = 1
    parinfo[i*10+10]['limits'][1] = int(0.5*model_list[0][1]['dim'])
    parinfo[i*10+10]['fixed'] = 1
  
  # Input to mpfit
  functkw = {'x':model_list, 'y':boxes_avg, 'err':[boxes, box_size]}
  
  # Fit parameters
  print model_list,boxes_avg
  m = mpfit.mpfit(model_func, parinfo=parinfo, functkw=functkw, quiet=0)
  p = m.params
  print p,m.status
  
  # Create output model_list
  for i in range(len(model_list)):
    model_list[i][1]['pa'] = p[1]
    model_list[i][1]['inclination'] = p[2]
    model_list[i][1]['exp_max'] = p[3]
    model_list[i][1]['r_max'] = p[4]
    model_list[i][1]['v_max'] = p[5]
    model_list[i][1]['v_system'] = p[6]
    model_list[i][1]['v_expansion'] = p[7]
    model_list[i][1]['a_scale'] = p[8]
    model_list[i][1]['centr_offset_x'] = p[9]
    model_list[i][1]['centr_offset_y'] = p[10]

  return model_list,parinfo

  
def model_func(p, fjac=None, x=None, y=None, err=None):
  """Only used by mpfit() in fit_parameters(). Creates a model using the
      parameters to test and returns the difference between the model
      and the data in the boxes given.
      """
  # Set meaningful names
  model_list = x
  boxes_avg = N.array(y)
  boxes = err[0]
  box_size = err[1]

  print p
  # Create a model_list
  for i in range(len(model_list)):
    model_list[i][1]['pa'][0] = p[1]
    model_list[i][1]['inclination'][0] = p[2]
    model_list[i][1]['exp_max'][0] = p[3]
    model_list[i][1]['r_max'][0] = p[4]
    model_list[i][1]['v_max'][0] = p[5]
    model_list[i][1]['v_system'][0] = p[6]
    model_list[i][1]['v_expansion'][0] = p[7]
    model_list[i][1]['a_scale'][0] = p[8]
    model_list[i][1]['centr_offset_x'][0] = p[9]
    model_list[i][1]['centr_offset_y'][0] = p[10]

  print model_list
  
  # Create the model
  model = create_vf(model_list)
  # Read the average in the boxes
  model_boxes_avg = N.array(get_boxes_avg(model, boxes, box_size))
  
  status = 0
  
  # Return status and model_boxes_avg-boxes_avg
  return [status, (model_boxes_avg - boxes_avg)]

  
def get_boxes(data):
  """Used by fit_parameters(). The velocity field is shown to the user that
      chooses boxes by clicking.
      """
  # Clean the click-file
  try:
    os.remove('/tmp/MPclick.dat-nono')
  except:
    pass
  
  # Show the velocity field and connect click-events to on_click_float
  lower=data.min()-1
  upper=data.max()+1
  font = {'fontname'   : 'Courier','color'      : 'k','fontsize'   : 20}

  MP.figure(num=1, figsize=(8.14, 8), dpi=80, facecolor='w', edgecolor='k')
  MP.imshow(N.swapaxes(data,0,1), vmin=lower, vmax=upper, interpolation='nearest', origin='lower', aspect='preserve')
  #MP.colorbar()
  #PyCigale.setXaxis_pc(data)
  #PyCigale.setYaxis_pc(data)
  MP.title(data.p['objname'] + ' - ' + 'Radial Velocity',font)
  MP.axis([0,data.nx()-1,0,data.ny()-1])
  MP.connect('button_press_event', PyCigale.on_click_float)
  MP.show()
  
  # Read the clicks and put the coordinates in the boxes-list
  boxes = []
  for line in open('/tmp/MPclick.dat', 'r').readlines():
    line = line.split()
    boxes += [[float(line[0]), float(line[1])]]
  
  return boxes
  

def get_boxes_avg(data, boxes, box_size):
  """Used by fit_parameters(). Computes the average value in each box defined
      by boxes and box_size.
      """
  boxes_avg = N.zeros(len(boxes),type='Float32')
  for i in range(len(boxes)):
    p1 = N.maximum(0,int(boxes[i][0]-box_size/2))
    p2 = N.minimum((data.nx()-1),int(boxes[i][0]+box_size/2))
    p3 = N.maximum(0,int(boxes[i][1]-box_size/2))
    p4 = N.minimum((data.ny()-1),int(boxes[i][1]+box_size/2))
    boxes_avg[i]= data[p1:p2,p3:p4].mean()
  
  return boxes_avg


def create_vf(model_list):
  vf = range(len(model_list))
  for i in range(len(model_list)):
    if (model_list[i][0] == 'system'):
      vf[i] = create_system_vf(model_list[i][1])
    elif (model_list[i][0] == 'linear'):
      vf[i] = create_rot_vf(model_linear, model_list[i][1])
    elif (model_list[i][0] == 'kepler'):
      vf[i] = create_rot_vf(model_kepler, model_list[i][1])
    elif (model_list[i][0] == 'pure_kepler'):
      vf[i] = create_rot_vf(model_pure_kepler, model_list[i][1])
    elif (model_list[i][0] == 'disk'):
      vf[i] = create_rot_vf(model_disk, model_list[i][1])
    elif (model_list[i][0] == 'expansion'):
      vf[i] = create_exp_vf(model_expansion, model_list[i][1])
  
  vf_total = 0
  for i in range(len(model_list)):
    vf_total = vf_total + vf[i]
  
  return PyCigale.array(vf_total)


def parameter_dict(dim=512, pa=0.0, inclination=0.0, exp_max=512, r_max=512, v_max=100, v_system=0, v_expansion=0, a_scale=1, centr_offset_x=0, centr_offset_y=0):
  """Returnes a dictionary with parameters for a model.
      
      Usage: pars = parameter_dict(dim,pa,inclination,exp_max,r_max,v_max,
        v_system,v_expansion,a_scale,centr_offset_x,centr_offset_y)
        
      pars: The dictionary with parameters
      dim..centr_offset_y: The parameters for model_list
      """
  pars = {}
  pars['dim'] = int(dim)
  pars['pa'] = pa
  pars['inclination'] = inclination
  pars['r_max'] = r_max
  pars['v_max'] = v_max
  pars['v_system'] = v_system
  pars['v_expansion'] = v_expansion
  pars['a_scale'] = a_scale
  pars['exp_max'] = exp_max
  pars['centr_offset_x'] = centr_offset_x
  pars['centr_offset_y'] = centr_offset_y

  return pars

  
def arguments_list(pars):
  """Produces two arrays, one for the radius in each point and one for the
      angle in each point.
      
      The parameters from a parameter dictionary is used to create a map 
      of the angle and radius for each point in the galaxy (in its own 
      coordinate system) viewed through the telescope.
      
      For each point in the velocity field it calculates what radius and angle
      that point corresponds to in the galaxy-coordinates (where pa
      defines the major axis).
      
      This method is from the Fortran-program.
      
      Usage r_and_phi = arguments_list(pars)
      
      pars: A parameter dictionary
      r_and_phi: [r, phi], Two 2D-arrays with the transformed coordinates.
      """
  # Coordinates for the galaxy centrum
  centr_x = pars['dim']/2 +pars['centr_offset_x']
  centr_y = pars['dim']/2 +pars['centr_offset_y']
  
  # The empty maps
  r = N.zeros((pars['dim'],pars['dim']), type='Float64')
  phi = N.zeros((pars['dim'],pars['dim']), type='Float64')
   
  # Transform each point
  for x in range(pars['dim']):
    for y in range(pars['dim']):
      temp_x = x-centr_x
      temp_y = y-centr_y
      
      pa = (-pars['pa']-90)*M.pi/180
      inclination = pars['inclination']*M.pi/180
      
      x_corr = temp_x*N.cos(pa) - temp_y*N.sin(pa)
      y_corr = temp_x*N.sin(pa) + temp_y*N.cos(pa)
      
      y_corr = y_corr*N.cos(inclination)
      
      phi[x,y] = N.arctan2(y_corr, x_corr)
      r[x,y] = N.sqrt(x_corr**2+y_corr**2)

  r_and_phi = [r, phi]
  
  return r_and_phi

  
def create_exp_vf(model, pars):
  """Creates a VF one of the model-functions (model_expansion is the only one 
      available atm.) with the parameter dictionary 'pars'.
      
      Usage: vf_exp = create_expansion_vf(model, pars)
      
      pars: The parameter dictionary to use
      vf_exp: A velocity field (note: this field is just a pure numarray, so no 
        p-list). If you want to keep the p-list from 'vf' you can do this:
        vf_exp = vf.get_copy()
        vf_exp[:,:] = create_expansion_vf(pars)
  """
  r_and_phi = arguments_list(pars)

  arguments = [r_and_phi[0], pars]
  v_model_exp = apply(model_expansion, arguments)

  vf = v_model_exp*N.sin(r_and_phi[1])*N.sin(pars['inclination'][0]*M.pi/180)
   
  return vf

  
def create_rot_vf(model, pars):
  """Creates a VF using one of the model-functions (model_disk, model_kepler, 
      model_linear, model_pure_kepler) with the parameter dictionary 'pars'.
      
      Usage: vf_rot = create_rot_vf(model, pars)
      
      model: One of the functions listed above
      pars: The parameter dictionary to use
      vf_rot: A velocity field (note: this field is just a pure numarray, so no 
        p-list). If you want to keep the p-list from 'vf' you can do this:
        vf_rot = vf.get_copy()
        vf_rot[:,:] = create_rot_vf(model, pars)
  """
  r_and_phi = arguments_list(pars)

  arguments = [r_and_phi[0], pars]
  # Call model-function
  v_model_rot = apply(model, arguments)

  # Projection to the rotation plane
  vf = v_model_rot*N.cos(r_and_phi[1])*N.sin(pars['inclination']*M.pi/180)
  
  return vf

  
def create_system_vf(pars):
  """Creates a VF with a system-velocity using the parameter dictionary 'pars'.
      
      Usage: vf_sys = create_system_vf(pars)
      
      pars: The parameter dictionary to use
      vf_sys: A velocity field (note: this field is just a pure numarray, so no 
        p-list). If you want to keep the p-list from 'vf' you can do this:
        vf_sys = vf.get_copy()
        vf_sys[:,:] = create_system_vf(pars)
  """
  #print pars['v_system']
  vf = N.zeros((pars['dim'],pars['dim']), type='Float64') + pars['v_system']
  return vf
  
  
def model_expansion(r, pars):
  """Creates an expansion velocity field. The velocity is v_max in the centrum
      and decreases as 1/r**2, stretched such that it is 0.1*v_max at the radius 
      exp_max.
      
      Usage: arr = model_expansion(r, pars)
      
      r: An array with radii in every point, preferably from argument_list()
      pars: The parameter dictionary to use
      arr: An array with velocities in every point
      """
  # Bugfix, we don't want to change r, just a copy of it
  arr = r.copy()
  
  # Set variables
  len_x = arr.shape[0]
  len_y = arr.shape[1]
  
  v_max = pars['v_expansion'][0]
  exp_max = pars['exp_max'][0]

  # Special case if exp_max=0
  if (exp_max == 0):
    arr[:,:] = 0
    return arr
  
  #Calculate the velocity in ervery point
  arr.setshape((len_x*len_y))

  arr = v_max/(( ((N.sqrt(10)-1)/exp_max)*arr+1)**2)
  
  arr.setshape((len_x, len_y)) 
  return arr

  
def model_disk(r, pars):
  """Creates an disk rotation velocity field. 
      
      Usage: arr = model_disk(r, pars)
      
      r: An array with radii in every point, preferably from argument_list()
      pars: The parameter dictionary to use
      arr: An array with velocities in every point
      """
  # Bugfix, we don't want to change r, just a copy of it
  arr = r.copy()
  
  # Set variables
  len_x = arr.shape[0]
  len_y = arr.shape[1]
  
  a = pars['a_scale']
  vm = pars['v_max']

  # Calculate the velocity in ervery point
  arr.setshape((len_x*len_y))
  
  t1 = N.where(arr==0)
  t2 = N.where(arr>0)
  
  arr[t1] = 0.0
  arr[t2] = N.sqrt(N.fabs(vm**2 * ( 1 - ( N.arctan2(arr[t2],a) * (a/arr[t2]) ) ) ) )
  
  arr.setshape((len_x, len_y)) 
  
  return arr
  
  
def model_kepler(r, pars):
  """Creates an kepler rotation velocity field. Linear up to v_max*r/r_max and 
      then decreasing as v_max*sqrt(r_max/r).
      
      Usage: arr = model_disk(r, pars)
      
      r: An array with radii in every point, preferably from argument_list()
      pars: The parameter dictionary to use
      arr: An array with velocities in every point
      """
  # Bugfix, we don't want to change r, just a copy of it
  arr = r.copy()
  
  # Set variables
  len_x = arr.shape[0]
  len_y = arr.shape[1]
  
  rm = pars['r_max'][0]
  vm = pars['v_max'][0]
  
  # Special case when r_max = 0
  if (rm == 0):
    arr[:,:] = 0
    return arr
  
  # Calculate the velocity in ervery point
  arr.setshape((len_x*len_y))
  
  t1 = N.where(arr<rm)
  t2 = N.where(arr>=rm)
  
  arr[t1] = vm*arr[t1]/rm
  arr[t2] = vm*N.sqrt(rm/arr[t2])
  
  arr.setshape((len_x, len_y))
  return arr
  

def model_pure_kepler(r, pars):
  """Creates an pure kepler rotation velocity field. Decreasing as 
      v_max*sqrt(r_max/r). The velocity is set to 0 in the centrum.
      
      Usage: arr = model_disk(r, pars)
      
      r: An array with radii in every point, preferably from argument_list()
      pars: The parameter dictionary to use
      arr: An array with velocities in every point
      """
  # Bugfix, we don't want to change r, just a copy of it
  arr = r.copy()
  
  # Set variables
  len_x = arr.shape[0]
  len_y = arr.shape[1]
  
  rm = pars['r_max'][0]
  vm = pars['v_max'][0]
  
  # Special case when r_max = 0
  if (rm == 0):
    arr[:,:] = 0
    return arr
  
  # Calculate the velocity in ervery point
  arr.setshape((len_x*len_y))
  
  t1 = N.where(arr<0.5)
  t2 = N.where(arr>=0.5)
  
  arr[t1] = 0.0
  arr[t2] = vm*N.sqrt(rm/arr[t2])
  
  arr.setshape((len_x, len_y))
  return arr
  
  
def model_linear(r, pars):
  """Creates an linear rotation velocity field. Increasing as v_max*r/r_max,
      up to r_max, then it is held fix (at v_max).
      
      Usage: arr = model_disk(r, pars)
      
      r: An array with radii in every point, preferably from argument_list()
      pars: The parameter dictionary to use
      arr: An array with velocities in every point
      """
  # Bugfix, we don't want to change r, just a copy of it
  arr = r.copy()
  
  # Set variables
  len_x = arr.shape[0]
  len_y = arr.shape[1]
  
  rm = pars['r_max']
  vm = pars['v_max']

  # Calculate the velocity in ervery point
  arr.setshape((len_x*len_y))
  
  t1 = N.where(arr<rm)
  t2 = N.where(arr>=rm)
  
  arr[t1] = vm*arr[t1]/rm
  arr[t2] = vm
  
  arr.setshape((len_x, len_y))
  return arr
 


 
###############################################
# The old stuff
###############################################

def model(model,name,*args):
  
  if model == 'linear':
    print "creating the linear field..."
    result,err=create_linear (args)
    if err != 0:
      print "error"
    else:  
      print "done."
    
    
  elif model == 'disk':
    print "creating the disk..."
    result,err=create_disk (args)
    if err != 0:
      print "error"
    else:  
      print "done."
    
  elif model == 'shell':
    print "creating the shell..."
    result,err=create_shell (args)
    if err != 0:
      print "error"
    else:  
      print "done."

  elif model == 'sphere':
    print "creating the sphere..."
    result,err=create_sphere (args)
    if err != 0:
      print "error"
    else:  
      print "done."
    
  else:
    print "unknown model"
    result,err=N.zeros((1,1)),1

  # write the file if all went fine
  if err == 0:
    os.popen('rm ' + name)
    IO.write_fits(result, name)

  else:
    print "no file written"



def create_linear (arg):
  """
  makes a velocity field with a simple gradient.
  rotation possible.
  """
  if len(arg) != 3:
    print "Wrong number of arguments"
    return N.zeros((1,1)),1
  else:
    dim=arg[0]             # dimension of the file.
    vel=float(arg[1])      # peak velocity at r=dim (it becomes bigger than
                           # that in the corners if you rotate of course).
    rot=M.radians(arg[2])  # rotation from top to left in degrees.
    
  vf = N.zeros((dim,dim), type=N.Float32)
  cent = dim / 2.
  diff = vel / cent
  normvec=[M.sin(rot),M.cos(rot)]
  
  for i in range(dim):
    icent=i-cent
    for j in range(dim):

      vf[i,j]= diff * N.dot(normvec,[icent,j-cent])

      
  return vf,0



def create_disk (arg):
  """
  creates a velocity field of a thin disk galaxy
  according to Sparke, Gallagher's
  "Galaxies in the Universe", equations 2.20 and 5.4 .
  Parameters taken from there as well.
  """
  if len(arg) != 5:
    print "Wrong number of arguments"
    return N.zeros((1,1)),1
  else:
    dim=arg[0]                # dimension of the file.
    vel=float(arg[1])         # max velocity
    a=float(arg[2])           # scale parameter
    incl=M.radians(arg[3])    # inclination
    rot=M.radians(arg[4])     # rotation from top to left in degrees.
    
  vf = N.zeros((dim,dim), type=N.Float32)
  cent = dim / 2.
  minaxis=[M.sin(rot),M.cos(rot)]
  majaxis=[M.cos(rot),-M.sin(rot)]


  for i in range(dim):
    icent=i-cent
    for j in range(dim):
      
      d_min=N.dot(majaxis,[icent,j-cent])
      d_maj=N.dot(minaxis,[icent,j-cent]) / N.cos(incl)
      
      r = N.hypot(d_min, d_maj)
      if r == 0.0:
        vf[i,j]=0

      else:  
        
        v_r = v_r_disk(vel,a,r)

        vf[i,j]= v_r * N.sin(incl) * (d_min / r) # note that (d_min / r) = cos Phi
        

  return vf,0    


def create_shell(arg):
  """
  Velocity field of a rotating (opt thick) ball, looking at the surface
  Parameter vel is the velocity IN the image plane at the center.
  R is the radius of the sphere which need not be the same as the
  image size.
  rotation means the rot-axis from top to left in degrees.
  """
  if len(arg) != 4:
    print "Wrong number of arguments"
    return N.zeros((1,1)),1
  else:
    dim=arg[0]                # dimension of the file.
    R=float(arg[1])           # Radius (see above)
    vel=float(arg[2])         # vel. (see above)
    rot=M.radians(arg[3])     # rotation from top to left in degrees.
    
  vf = N.zeros((dim,dim), type=N.Float32)
  cent = dim / 2.
  perpaxis=[M.sin(rot),M.cos(rot)]
  rotaxis=[M.cos(rot),-M.sin(rot)] 
  omega=vel / R


  for i in range(dim):
    icent=i-cent
    for j in range(dim):
      
      d_rot=N.dot(perpaxis,[icent,j-cent]) # distance from rotation axis
      d_cent=N.dot(rotaxis,[icent,j-cent]) # distance from center along rot-axis
      d_centS = d_cent / R * M.pi / 2.
      r = N.hypot(d_rot, d_cent)
      if r == 0:
        vf[i,j]= 0
      elif r >= R:
        vf[i,j]= 0
      else:
        vf[i,j]= d_rot * omega
        #r_ring= N.sqrt(R**2 - d_cent**2) # radius of the circle a height d_cent
        #alpha = N.arcsin(d_rot / r_ring)
        #vf[i,j]= N.sin(alpha) * (omega * r_ring)
        #N.sqrt((omega*d_rot)**2 * ( 1 - (d_rot / r_ring)**2 ) )
        #v_r_ring(r_ring * omega,d_rot,r)
        
  return vf,0    
  


def create_sphere(arg):
  """
  Velocity field of a rotating opt. thin ball, looking through it.
  This will be more tricky and probably needs more assumtions.
  """
  return N.zeros((1,1)),1


  

def v_r_disk(vel,a,r):
  """
  V(r) according to eq. 2.20 from Sparke, Gallagher
  """
  return  N.sqrt(N.fabs(vel**2 * ( 1 - ( N.arctan2(r,a) * (a/r) ) ) ) ) 



def get_box(vf,x):
  """
  gives the average in predefined boxes
  """

  if x == 1:
    aver=mean(mean(vf[20:30,100:110]))
  elif x==2:
    aver=mean(mean(vf[70:80,110:120]))
  elif x==3:
    aver=mean(mean(vf[95:105,95:105]))
  elif x==4:
    aver=mean(mean(vf[105:115,65:75]))
  elif x==5:
    aver=mean(mean(vf[95:105,25:35]))
  elif x==6:
    aver=mean(mean(vf[50:60,15:25]))
  elif x==7:
    aver=mean(mean(vf[20:30,40:50]))
  elif x==8:
    aver=mean(mean(vf[20:30,70:80]))
  elif x==9:
    aver=mean(mean(vf[20:30,100:110]))
  elif x==10:
    aver=mean(mean(vf[20:30,100:110]))
  else:
    aver=0

  return aver


def modelfu_linear(param,x):
  """
  gives back the model velocity at the desired points
  """

  dim=64.

  if x == 1:
    X=25
    Y=105
  elif x==2:
    X=75
    Y=115
  elif x==3:
    X=100
    Y=100
  elif x==4:
    X=110
    Y=70
  elif x==5:
    X=100
    Y=30
  elif x==6:
    X=55
    Y=20
  elif x==7:
    X=25
    Y=45
  elif x==8:
    X=25
    Y=75
  else:
    return 1

  rot=param[1] / 360 * 2 * M.pi
  #print param[0], rot, rot / 2 / M.pi * 360, X, Y
  return (param[0]/dim) * Nu.dot([Nu.sin(rot),Nu.cos(rot)],[X-dim,Y-dim])



def fit_linear(name):
  """
  testing a least square fit
  """



  rawdata=IO.read_fits(name)

  data= [ (1, get_box(rawdata,1)),
          (2, get_box(rawdata,2)),
          (3, get_box(rawdata,3)),
          (4, get_box(rawdata,4)),
          (5, get_box(rawdata,5)),
          (6, get_box(rawdata,6)),
          (7, get_box(rawdata,7)),
          (8, get_box(rawdata,8)) ]
  
  result,chi = LS.leastSquaresFit(modelfu_linear, (100,-30,), data)
  print result[0], M.fmod(result[1],360), chi



def fit_test():
  """
  testing a least square fit
  """

  def modelfu_test(param,t):
    return param[0]*t

  data= [ (1, 1),
          (2, 2),
          (3, 3),
          (4, 4),
          (5, 5),
          (6, 6),
          (7, 7),
          (8, 8) ]
  
  print LS.leastSquaresFit(modelfu_test, (2,), data)


