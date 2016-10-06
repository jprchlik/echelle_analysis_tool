import matplotlib
matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
#implement the deault mpl key bindings
from matplotlib.backend_bases import key_press_handler,MouseEvent
import tkMessageBox as box
import tkFileDialog as Tkf
import numpy as np
import pyfits


#import subprocess and glob to run interpolation routine for MARCS and MOOG 
import glob,subprocess

#reads in moog file
import read_moog

#Scipy for spline
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d

#scipy for integration
from scipy.integrate import simps

#Astropymodeling for gauss fit

from astropy.modeling import models, fitting

import matplotlib.pyplot as plt
import sys,os

#For Gaussian smoothing
from astropy.convolution import Gaussian1DKernel, convolve

if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk

from fancy_plot import *
#from astroML.plotting import setup_text_plots

#setup_text_plots(usetex=True,fontsize=28)

#dictionary to convert proton number into element symbol
edict = {}
edict['1'] = 'H'
edict['2'] = 'Hi'
edict['3'] = 'Li'
edict['4'] = 'Be'
edict['5'] = 'B'
edict['6'] = 'C'
edict['7'] = 'N'
edict['8'] = 'O'
edict['9'] = 'F'
edict['10'] = 'Ne'
edict['11'] = 'Na'
edict['12'] = 'Mg'
edict['13'] = 'Al'
edict['14'] = 'Si'
edict['15'] = 'P'
edict['16'] = 'S'
edict['17'] = 'Cl'
edict['18'] = 'Ar'
edict['19'] = 'K'
edict['20'] = 'Ca'
edict['21'] = 'Sc'
edict['22'] = 'Ti'
edict['23'] = 'V'
edict['24'] = 'Cr'
edict['25'] = 'Mn'
edict['26'] = 'Fe'
edict['27'] = 'Co'
edict['28'] = 'Ni'
edict['29'] = 'Cu'
edict['30'] = 'Zn'
edict['31'] = 'Ga'
edict['32'] = 'Ge'
edict['33'] = 'As'
edict['34'] = 'Se'
edict['35'] = 'Br'
edict['36'] = 'Kr'
edict['37'] = 'Rb'
edict['38'] = 'Sr'
edict['39'] = 'Y'
edict['40'] = 'Zr'
edict['41'] = 'Nb'
edict['42'] = 'Mo'
edict['43'] = 'Tc'
edict['44'] = 'Ru'
edict['45'] = 'Rh'
edict['46'] = 'Pd'
edict['47'] = 'Ag'
edict['48'] = 'Cd'
edict['49'] = 'In'
edict['50'] = 'Sn'
edict['51'] = 'Sb'
edict['52'] = 'Te'
edict['53'] = 'I'
edict['54'] = 'Xe'
edict['55'] = 'Cs'
edict['56'] = 'Ba'
edict['57'] = 'La'
edict['58'] = 'Ce'
edict['59'] = 'Pr'
edict['60'] = 'Nd'
edict['61'] = 'Pm'
edict['62'] = 'Sm'
edict['63'] = 'Eu'
edict['64'] = 'Gd'
edict['65'] = 'Tb'
edict['66'] = 'Dy'
edict['67'] = 'Ho'
edict['68'] = 'Er'
edict['69'] = 'Tm'
edict['70'] = 'Yb'
edict['71'] = 'Lu'
edict['72'] = 'Hf'
edict['73'] = 'Ta'
edict['74'] = 'W'
edict['75'] = 'Re'
edict['76'] = 'Os'
edict['77'] = 'Ir'
edict['78'] = 'Pt'
edict['79'] = 'Au'
edict['80'] = 'Hg'
edict['81'] = 'Tl'
edict['82'] = 'Pb'
edict['83'] = 'Bi'
edict['84'] = 'Bo'
edict['85'] = 'At'
edict['86'] = 'Rn'
edict['87'] = 'Tr'
edict['88'] = 'Ra'
edict['89'] = 'Ac'
edict['90'] = 'Th'
edict['91'] = 'Pa'
edict['92'] = 'U'
edict['93'] = 'Np'
edict['94'] = 'Pu'
edict['95'] = 'Am'
edict['96'] = 'Cm'



#dictionary for ionization level
idict = {}
idict['0'] = 'I'
idict['1'] = 'II'
idict['2'] = 'III'
idict['3'] = 'IV'
idict['4'] = 'V'
idict['5'] = 'VI'
idict['6'] = 'VII'




#Convert from air to vacuum wavelengths
def air_vac(wav0):
    wav0 = wav0/(1.0 + 2.735182*10**(-4.0) + (131.4182/wav0**2.) + (2.76249*10.**8./wav0**4.))
    return wav0



#correct wavelengths for rv
def rv_correct(lam0,rv):

    cs = 2.99792458e5 #km/s
    lam = lam0*(np.sqrt((1.+(rv/cs))/(1.-(rv/cs))))

    return lam
#Pop-up to show abundance trends
class abun_popup:
    global root
    def __init__(self,parent,figcol,figrow):
#        top = self.top = Tk.Toplevel(parent)

#Create the figure
        self.abunfig,self.abunax = plt.subplots(nrows=figrow,ncols=figcol,figsize=(8,10))
#        self.abunfig.subplots_adjust(hspace=.001,wspace=.001)
#Create window for the plot
        self.abuncanvas = FigureCanvasTkAgg(self.abunfig,master=parent)
        self.abuncanvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1.0)
        self.abuncanvas.show()
#Display matplotlib widgets
        self.toolbar = NavigationToolbar2TkAgg(self.abuncanvas,parent)
        self.toolbar.update()
        self.abuncanvas._tkcanvas.pack(side=Tk.TOP,fill=Tk.BOTH,expand=1)
    def destroy(self):
        self.parent.destroy()


#Line EW measurement edit popup
class ew_popup:
    global root
    def __init__(self,parent):
        top = self.top = Tk.Toplevel(parent)
#        self.myLabel = Tk.Label(top,text='EW Edit')
#        self.myLabel.pack()
        self.gButton = Tk.Button(top,text='Gauss Fit',command=self.sendg)
        self.gButton.pack(side=Tk.LEFT)
        self.gButton = Tk.Button(top,text='Set Gauss',command=self.sendsg)
        self.gButton.pack(side=Tk.LEFT)
        self.gButton = Tk.Button(top,text='Simpson Rule Integration',command=self.sends)
        self.gButton.pack(side=Tk.LEFT)
        self.gButton = Tk.Button(top,text='Voigt Fit',command=self.sendv)
        self.gButton.pack(side=Tk.LEFT)
        self.gButton = Tk.Button(top,text='Delete',command=self.sendd)
        self.gButton.pack(side=Tk.LEFT)
        self.gButton = Tk.Button(top,text='Cancel',command=self.sendc)
        self.gButton.pack(side=Tk.LEFT)

    def sendd(self):
        self.editparm = 'd'
        self.top.destroy()
    def sendg(self):
        self.editparm = 'g'
        self.top.destroy()
    def sendsg(self):
        self.editparm = 'sg'
        self.top.destroy()
    def sendc(self):
        self.editparm = 'c'
        self.top.destroy()
    def sendv(self):
        self.editparm = 'v'
        self.top.destroy()
    def sends(self):
        self.editparm = 's'
        self.top.destroy()

#Class for asking for order number to plot first
#No longer used switched to text box
class order_popup:
    global root
    def __init__(self,parent):
        top = self.top = Tk.Toplevel(parent)
        self.myLabel = Tk.Label(top,text='Order')
        self.myLabel.pack()
        self.myEntryBox = Tk.Entry(top)
        self.myEntryBox.pack()
        self.mySubmitButton = Tk.Button(top,text='Submit',command=self.send)
        self.mySubmitButton.pack()

    def send(self):
        self.order = self.myEntryBox.get()
        self.top.destroy()

#Class for asking for fits dimension with echelle data
#Only does this if echelle data not in zeroth dimension 
class dimension_popup:
    global root
    def __init__(self,parent):
        top = self.top = Tk.Toplevel(parent)
        self.myLabel = Tk.Label(top,text='Fits Dimension')
        self.myLabel.pack()
        self.myEntryBox = Tk.Entry(top)
        self.myEntryBox.pack()
        self.mySubmitButton = Tk.Button(top,text='Submit',command=self.send)
        self.mySubmitButton.pack()

    def send(self):
        self.order = self.myEntryBox.get()
        self.top.destroy()




#Main gui create and modification loop
class gui_c(Tk.Frame):

    def __init__(self,parent):
        Tk.Frame.__init__(self,parent,background='white')
        self.parent = parent
#Default setting for finding the fits dimension
        self.fitsd = 0
#Intialize clicker
        self.click = 0
#Initial gauss find tolerance
        self.findtol = .2
#set default order
        self.order = 1
#radial velocity corrections to line list
#        self.rv = -269.16
#        self.hrv = -2.14
#radial velocity corrections to line list
        self.rv = 0.0 
        self.hrv = 0.0
#looking fore continuum
        self.seeking = False
#line list input
        self.linelistin = False
#EW have been measured
        self.ewmeas = False
#Order has been measured
        self.oewmeas = False
#Edit EW  measurement
        self.ewedit = False
#Set the line min and half max by hand
        self.setew = False
#maximum line width in one direction for gauss fit
        self.maxlinewidth = 2.
#Make sure file is saved before you destroy data
        self.saved = True
#Variable for simpsons integration
        self.sint = False
#line with in plot
        self.linewidth = 1
#Smoothing FWHM
        self.smofwhm = 0.0

#IS an ESO spectrum
        self.eso_spec = False

#Stellar parameters for Moog and Marcs 
        self.teff = 0.
        self.feh = 0.
        self.logg = 0.
        self.vturb = 0.

#check to see if marcs and moog alias are defined
#first initialize variables
        self.findmoog = False
        self.findmarcs = False
        self.modlistin = False
        self.modeldir = ''
        self.marcs_alias = 'interp_seis_marcs'
        self.marcs_alias = 'interp_marcs'
        self.checkalias()
#Do not let event work unless model has been interpolated
        self.marcsin = False

#set the maximum number of continuum points for spline interpolation
        self.contpointmax = 10


#arrays for writing information about the continuum
        self.acounts = []
        self.arootms = []
        self.apoints = []
        self.aorders = []
        self.awavels = []




#Ask whether models are parallel or spherical
#Initial equal parallel 
        self.ps = 'p'
#Ask whether models are alpha poor
#Initial equal normal galatic evolution
        self.alpha = 'a+0.40'

#Start the creation of the window and GUI
        self.centerWindow()
        self.FigureWindow()
        self.initUI()


#Check if marcs and moog alias are defined in .cshr
    def checkalias(self):
        checker = subprocess.Popen('csh -i -c alias',shell=True,stdout=subprocess.PIPE)
        for line in iter(checker.stdout.readline,''):
            splitline = line.split('\t')
            alias = splitline[0]
            if alias == 'moog':
#throw moog flag if moog alias is found
                self.findmoog = True
                self.moogloc = splitline[1].replace('\n','')
#throw moog flag if interp_marcs alias is found
            if alias == self.marcs_alias:
                self.findmarcs = True
                self.marcsloc = splitline[1].replace('\n','')

#Create area and window for figure
    def FigureWindow(self):
#set the information based on screen size
        x =  self.parent.winfo_screenwidth()
        y =  self.parent.winfo_screenheight()

        aratio = float(x)/float(y)
#Create the figure
        self.f,self.a = plt.subplots(figsize=(8*aratio,8*aratio*.5))
#Create window for the plot
        self.canvas = FigureCanvasTkAgg(self.f,master=self)
#Draw the plot
        self.canvas.draw()
#Turn on matplotlib widgets
        self.canvas.get_tk_widget().pack(side=Tk.TOP,fill=Tk.BOTH,expand=1)
#Display matplotlib widgets
        self.toolbar = NavigationToolbar2TkAgg(self.canvas,self)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=Tk.TOP,fill=Tk.BOTH,expand=1)
#Connect mpl to mouse clicking
        self.f.canvas.mpl_connect('button_press_event',self.on_click_event)
#Connect mpl to mouse clicking
        self.f.canvas.mpl_connect('key_press_event',self.on_key_event)

#create button for initialization of Continuum normalization
        specframe = Tk.Frame(self)
#Check box for alpha poor Models
        self.checkvaral = Tk.IntVar()
        self.C3 = Tk.Checkbutton(master=specframe,text="Alpha Poor",variable=self.checkvaral,onvalue=1,offvalue=0,command=self.alpha_mod)
        self.C3.pack(side=Tk.LEFT)
#Check box for Spherical Models
        self.checkvarsp = Tk.IntVar()
        self.C2 = Tk.Checkbutton(master=specframe,text="Spherical Mod",variable=self.checkvarsp,onvalue=1,offvalue=0,command=self.sp_mod)
        self.C2.pack(side=Tk.LEFT)
#Check box for vaccuum correction
        self.checkvar = Tk.IntVar()
        self.C1 = Tk.Checkbutton(master=specframe,text="Vaccum",variable=self.checkvar,onvalue=1,offvalue=0,command=self.vac_mod)
        self.C1.pack(side=Tk.LEFT)
#create button to go up an order
        upbutton = Tk.Button(master=specframe,text='Increase Order',command=self.increaseorder)
        upbutton.pack(side=Tk.LEFT)
#create button to go down an order
        downbutton = Tk.Button(master=specframe,text='Decrease Order',command=self.decreaseorder)
        downbutton.pack(side=Tk.LEFT)
#Put in same frame as other spec information
        downbutton = Tk.Button(master=specframe,text='Normalize Continuum',command=self.click_cont)
        downbutton.pack(side=Tk.LEFT)
#Create button to ask for Moog formated line list
        mooggrab = Tk.Button(master=specframe,text='Linelist',command=self.line_grab)
        mooggrab.pack(side=Tk.LEFT)
        specframe.pack(side=Tk.TOP)
#Create button to ask for Moog formated line list
        moogsave = Tk.Button(master=specframe,text='Save EW',command=self.line_save)
        moogsave.pack(side=Tk.RIGHT)
#Crate button to reset plot limits without clearing data
        plotrbutton = Tk.Button(master=specframe,text='Reset Limits',command=self.reset_limits)
        plotrbutton.pack(side=Tk.RIGHT)
#Crate button to reset plot limits without clearing data
        openbutton = Tk.Button(master=specframe,text='Open File',command=self.open_file)
        openbutton.pack(side=Tk.RIGHT)
#Add to run MOOG 
        self.moogb = Tk.Button(master=specframe,text='MOOG',command=self.startmoog)
        self.moogb.pack(side=Tk.RIGHT) 
#Add to run MARCS interpolation 
        self.marcb = Tk.Button(master=specframe,text='Interp. Phot.',command=self.marcinterp)
        self.marcb.pack(side=Tk.RIGHT) 
#Pack the spec frame together

#        frame2 = Tk.Frame(self)
        specframe.pack(side=Tk.TOP)



#Call MOOG
    def startmoog(self):
        if self.findmoog:
            self.marcinterp()
            if self.marcsin:
                try:
                    self.in_moog_lines()
                    self.write_abfind()
                    callmoog = subprocess.Popen(['csh -i -c moog'],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,close_fds=True,shell=True)
                    stdout = callmoog.communicate(input=b'y\nabfind_pipe.par\n'+'\n'*200)
                    self.plot_moog()
#Error because linelist has not been saved
                except AttributeError:
                    self.error = 107
                    self.onError()
#Error and say no model has been interpolated
            else:
                self.error = 108
                self.onError()
#Error because no moog alias is defined               
        else:
            self.error = 101
            self.onError()
#Start MARCS interpolation
    def marcinterp(self):
#Open the marcs model file to find the location of marcs models to be interpolated
        if self.findmarcs:
            marcsfile = open(self.marcsloc,'r')
            self.marcstext = marcsfile.readlines()
            marcsfile.close()
#Ue the list from before to find models
            if self.modlistin:
                self.select_marcs_models()
                findmarcs = subprocess.call("csh -i -c "+self.marcs_alias,shell=True)
#If the file list has not already been read in
            else:
                self.modeldir = self.marcstext[20].replace('\n','').split(' = ')[1]
                print self.modeldir
                self.find_marcs_models()
                self.select_marcs_models()
                findmarcs = subprocess.call("csh -i -c "+self.marcs_alias,shell=True)
#Format the marcs interpolated model atmosphere for moog
            self.marcs2moog()
#Tell the program a Marcs model has been read in
            self.marcsin = True
#If Marcs alias is not set throw an error when trying to interpolate
        else:
            self.error = 100
            self.onError()
            return

#Find list of marcs models
    def find_marcs_models(self):
        modellist = glob.glob(self.modeldir+'/'+self.ps+'*.mod')
        try:
            modellist.remove(self.modeldir+'/sun.mod')
        except ValueError:
            print 'No Sun in list'
        self.modelteffs = []
        self.modelfehs = []
        self.modelloggs = [] 
        for i in modellist:
            cut = i.replace(self.modeldir+'/','')
            cut = cut.replace('.mod','')
            self.modelteffs.append(cut[1:5])
            self.modelloggs.append(cut[7:11])
            self.modelfehs.append(cut[25:30])
#Convert the list into numpy arrays so the the bookends to the input parameters can be found
        self.modelteffs = np.unique(np.array(self.modelteffs))
        self.modelloggs = np.unique(np.array(self.modelloggs))
        self.modelfehs = np.unique(np.array(self.modelfehs))

#Selects 8 models to input into the interpolation file
    def select_marcs_models(self):
        fminfeh = float(self.feh)-self.modelfehs.astype(float)
        fminlogg = float(self.logg)-self.modelloggs.astype(float)
        fminteff = float(self.teff)-self.modelteffs.astype(float)

#find the bottom bracket values 
        tfehmin, = np.where(fminfeh >= 0.)
        tloggmin, = np.where(fminlogg >= 0.)
        tteffmin, = np.where(fminteff >= 0.)
#Return to main loop if value is out of range of parameters
        if tfehmin.size == 0:
            self.error = 110
            self.onError()
        if tloggmin.size == 0:
            self.error = 111
            self.onError()
        if tteffmin.size == 0:
            self.error = 112
            self.onError()

        ttfehmin, = np.where(fminfeh[tfehmin] == fminfeh[tfehmin].min())
        ttloggmin, = np.where(fminlogg[tloggmin] == fminlogg[tloggmin].min())
        ttteffmin, = np.where(fminteff[tteffmin] == fminteff[tteffmin].min())
#putvalues into variables
        self.fehminval = self.modelfehs[tfehmin][ttfehmin][0]
        self.loggminval = self.modelloggs[tloggmin][ttloggmin][0]
        self.teffminval = self.modelteffs[tteffmin][ttteffmin][0]
#find the top bracket values
        fmaxfeh = -fminfeh
        fmaxlogg = -fminlogg
        fmaxteff = -fminteff
        tfehmax, = np.where(fmaxfeh >= 0.)
        tloggmax, = np.where(fmaxlogg >= 0.)
        tteffmax, = np.where(fmaxteff >= 0.)
#Return to main loop if value is out of range of parameters
        if tfehmax.size == 0:
            self.error = 113
            self.onError()
        if tloggmax.size == 0:
            self.error = 114
            self.onError()
        if tteffmax.size == 0:
            self.error = 115
            self.onError()
        ttfehmax, = np.where(fmaxfeh[tfehmax] == fmaxfeh[tfehmax].min())
        ttloggmax, = np.where(fmaxlogg[tloggmax] == fmaxlogg[tloggmax].min())
        ttteffmax, = np.where(fmaxteff[tteffmax] == fmaxteff[tteffmax].min())

        self.fehmaxval = self.modelfehs[tfehmax][ttfehmax][0]
        self.loggmaxval = self.modelloggs[tloggmax][ttloggmax][0]
        self.teffmaxval = self.modelteffs[tteffmax][ttteffmax][0]

#create model list
        print self.alpha
        print self.modeldir+'/'+self.ps+self.teffminval+'_g'+self.loggminval+'_*z'+self.fehminval+'_'+self.alpha
        self.m1 = 'set model1 = '+glob.glob(self.modeldir+'/'+self.ps+self.teffminval+'_g'+self.loggminval+'_*z'+self.fehminval+'_'+self.alpha+'*.mod')[0].replace(self.modeldir+'/','')+'\n'
        self.m2 = 'set model2 = '+glob.glob(self.modeldir+'/'+self.ps+self.teffminval+'_g'+self.loggminval+'_*z'+self.fehmaxval+'_'+self.alpha+'*.mod')[0].replace(self.modeldir+'/','')+'\n'
        self.m3 = 'set model3 = '+glob.glob(self.modeldir+'/'+self.ps+self.teffminval+'_g'+self.loggmaxval+'_*z'+self.fehminval+'_'+self.alpha+'*.mod')[0].replace(self.modeldir+'/','')+'\n'
        self.m4 = 'set model4 = '+glob.glob(self.modeldir+'/'+self.ps+self.teffminval+'_g'+self.loggmaxval+'_*z'+self.fehmaxval+'_'+self.alpha+'*.mod')[0].replace(self.modeldir+'/','')+'\n'
        self.m5 = 'set model5 = '+glob.glob(self.modeldir+'/'+self.ps+self.teffmaxval+'_g'+self.loggminval+'_*z'+self.fehminval+'_'+self.alpha+'*.mod')[0].replace(self.modeldir+'/','')+'\n'
        self.m6 = 'set model6 = '+glob.glob(self.modeldir+'/'+self.ps+self.teffmaxval+'_g'+self.loggminval+'_*z'+self.fehmaxval+'_'+self.alpha+'*.mod')[0].replace(self.modeldir+'/','')+'\n'
        self.m7 = 'set model7 = '+glob.glob(self.modeldir+'/'+self.ps+self.teffmaxval+'_g'+self.loggmaxval+'_*z'+self.fehminval+'_'+self.alpha+'*.mod')[0].replace(self.modeldir+'/','')+'\n'
        self.m8 = 'set model8 = '+glob.glob(self.modeldir+'/'+self.ps+self.teffmaxval+'_g'+self.loggmaxval+'_*z'+self.fehmaxval+'_'+self.alpha+'*.mod')[0].replace(self.modeldir+'/','')+'\n'
#write models into new file
        self.marcstext[27] = self.m1
        self.marcstext[28] = self.m2
        self.marcstext[29] = self.m3
        self.marcstext[30] = self.m4
        self.marcstext[31] = self.m5
        self.marcstext[32] = self.m6
        self.marcstext[33] = self.m7
        self.marcstext[34] = self.m8
#write stellar parameters into new file         
        self.marcstext[37] = 'foreach Tref   ( '+str(int(float(self.teff)))+' )\n'
        self.marcstext[38] = 'foreach loggref   ( '+str(round(float(self.logg),2))+' )\n'
        self.marcstext[39] = 'foreach zref   ( '+str(round(float(self.feh),2))+' )\n'
#overwrite previous .com file
        marcsfile = open(self.marcsloc,'w')
        for i in self.marcstext:
            marcsfile.write(i)
        marcsfile.close()
        
         

#Resets the limits of the plot
    def reset_limits(self):
        self.a.set_xlim([self.wav0.min(),self.wav0.max()])
        maxer = np.max(self.spec0[self.spec0.size/2-.1*self.spec0.size:self.spec0.size/2+.1*self.spec0.size])
#Find the middle vaule in array to set limit if it is finite
        if np.isfinite(maxer):
            self.a.set_ylim([0.,maxer*.1+maxer])
        else:
            self.a.set_ylim([0.,np.max(self.spec0[np.isfinite(self.spec0)])])
        self.canvas.draw()

#Save Updated EW
    def line_save(self):
        self.linelistout = Tkf.asksaveasfile(mode='w')
        linedata = open(self.linelistout.name+'_linedata.dat','w')
        linedata.write('Measured by echelle analysis tool\n')
        linedata.write('    Wave       Ion       RHM       BHM      Depth     Cont       sigmad       sigmaHM    sigmadHM\n')
        if self.linelistout is None:
            return
#Added header so EW can be read back in
        self.linelistout.write('EW Measured by echelle analysis tool\n')
        if self.checkvar.get() == 0:
#Write out line ew in MOOG format
            for i,j in enumerate(self.line):
                self.linelistout.write('{0:^10.2f}{1:^10}{2:^10.2f}{3:^10.4f}{5}{4:^10.2f}\n'.format(self.line[i],str(self.lineel[i]),self.lineep[i],self.linegf[i],self.lineew[i],' '*20))
                linedata.write('{0:15.2f}{1:15.1f}{2:15.5f}{3:15.5f}{4:15.3f}{5:15.4f}{6:15.6f}{7:15.6f}{8:15.6f}\n'.format(self.line[i],float(self.lineel[i]),self.linerhm[i],self.linebhm[i],self.linedep[i],self.linecot[i],self.linesa1[i],self.linesa3[i],self.linesa13[i]))
#Always write out in air wavelengths
        if self.checkvar.get() == 1:
#Write out line ew in MOOG format
            for i,j in enumerate(self.aline):
                self.linelistout.write('{0:^10.2f}{1:^10}{2:^10.2f}{3:^10.4f}{5}{4:^10.2f}\n'.format(self.aline[i],str(self.lineel[i]),self.lineep[i],self.linegf[i],self.lineew[i],' '*20))
                linedata.write('{0:15.2f}{1:15.1f}{2:15.5f}{3:15.5f}{4:15.3f}{5:15.4f}{6:15.6f}{7:15.6f}{8:15.6f}\n'.format(self.aline[i],float(self.lineel[i]),self.linerhm[i],self.linebhm[i],self.linedep[i],self.linecot[i],self.linesa1[i],self.linesa3[i],self.linesa13[i]))
#close the file
        self.linelistout.close()
        linedata.close() 

#Removes flag to warn user of unsave data
        self.saved = True

#Grab Moog line list
    def line_grab(self):
        m = 0
#If the values already exist set the line ew measure to True
        self.ewmeas = True
        while m == 0:
            self.linelist = self.get_file()
            try:
                fil = open(self.linelist,'r')
                fils = fil.readlines()[1:]
                line = []
                lineel = []
                lineels = []
                lineep = []
                linegf = []
                lineew = []
                for i in fils:
                    line.append(float(i[:10]))
                    lineel.append(i[10:20].replace(' ',''))
                    eltext = i[10:20].replace(' ','').split('.')
                    ion = idict[eltext[1]]
                    el = edict[eltext[0]]
                    lineels.append(el+ion)
                    lineep.append(float(i[20:30]))
                    linegf.append(float(i[30:40]))
#Try to grab EW from line list file
                    if len(i) >= 60:
                        lineew.append(round(float(i[60:80]),1))
#Append value -9999.9 if EW have not been measured
                    else:
#If the any values don't exist set the line ew measure to False
                        self.ewmeas = False
                        lineew.append(round(-9999.9,1))
                
                m =1    
            except ValueError:
                print linegf
                print lineew
                self.error = 7
                self.onError()


        self.line = np.array(line)
        self.lineel = np.array(lineel)
        self.lineels = np.array(lineels)
        self.lineep = np.array(lineep)
        self.linegf = np.array(linegf)
        self.lineew = np.array(lineew)
        try:
            ld =  np.loadtxt(self.linelist+'_linedata.dat',skiprows=2)
            self.linebhm = ld[:,3] 
            self.linerhm = ld[:,2]
            self.linedep = ld[:,4] 
            self.linecot = ld[:,5] 
            self.linesa1 = ld[:,6]
            self.linesa3 = ld[:,7]
            self.linesa13 = ld[:,8]
        except IOError:
            self.linebhm =  np.zeros(self.lineew.size)-9999.9
            self.linerhm =  np.zeros(self.lineew.size)-9999.9
            self.linedep =  np.zeros(self.lineew.size)-9999.9
            self.linecot =  np.zeros(self.lineew.size)-9999.9
            self.linesa1 =  np.zeros(self.lineew.size)-9999.9
            self.linesa3 =  np.zeros(self.lineew.size)-9999.9
            self.linesa13 =  np.zeros(self.lineew.size)-9999.9
        self.linelistin = True
        self.text_loc()
        self.plot_lines()        

#Find location for labeling lines
    def text_loc(self):
#If the spectrum is normalized just point linelist at normal y values
        if np.median(self.spec0) < 10:
            self.elloc = 1.01+np.zeros(self.line.size)
            self.ewloc = 1.05+np.zeros(self.line.size)
#Else try half heartly to fit the continuum and scale
        else:
            forfit, = np.where(np.abs(self.spec0-np.median(self.spec0)) <= np.std(self.spec0)*3.)
            input = rv_correct(self.line,self.rv-self.hrv)
            specloc = np.interp(input,self.wav0[forfit],self.spec0[forfit])
            self.elloc = specloc*.1+specloc
            self.ewloc = specloc*.15+specloc

#Plot lines from the Moog line list
    def plot_lines(self):
        maxw = max(self.wav0)
        minw = min(self.wav0)
        findlines = rv_correct(self.line,self.rv-self.hrv)
        use, = np.where((findlines >= minw) & (findlines <= maxw))
#locations of lines in order
#location on spectrum
        self.olines = findlines[use]
#rest wavelengths of lines in order
        self.orestlines = self.line[use]

#An array containing the text values to be updated when a manual measurement is applied
        self.ewtext = []

        for i in use:
            self.a.text(rv_correct(self.line[i],self.rv-self.hrv),self.elloc[i],str(self.lineels[i]),rotation='vertical',fontsize=18)
            self.ewtext.append(self.a.text(rv_correct(self.line[i],self.rv-self.hrv),self.ewloc[i],str(self.lineew[i]),rotation='vertical',fontsize=16))
        self.canvas.draw()

#Do wavelength conversion of linelist to vaccum
    def vac_mod(self):
#release cursor from entry box and back to the figure
#needs to be done otherwise key strokes will not work
        self.f.canvas._tkcanvas.focus_set()
        if ((self.checkvar.get() == 1) & (self.linelistin)):
            self.aline = self.line.copy()
            self.line = air_vac(self.line)
            self.spec_plot()
        else:
#Uncheck box if linelist is not read in
            if (self.linelistin) == False:
                self.error = 8
                self.onError()
#Uncorrect if mistakenly pressed
            if (self.checkvar.get() == 0):
                self.line = self.aline
                self.spec_plot()
            self.C1.deselect()

#Check to see if asking for alpha poor models
    def alpha_mod(self):
        self.f.canvas._tkcanvas.focus_set()
        if (self.checkvaral.get() == 1):
            self.alpha = 'a+0.00'
        else:
            self.alpha = 'a+0.40'
            self.C3.deselect()


#Check whether models are spherical or plane-parallel
    def sp_mod(self):
#release cursor from entry box and back to the figure
#needs to be done otherwise key strokes will not work
        self.f.canvas._tkcanvas.focus_set()
        if (self.checkvarsp.get() == 1):
            self.ps = 's'
        else:
            self.ps = 'p'
            self.C2.deselect()
        
        

#Basic click event
    def on_click_event(self,click):
#Click envents for continuum selection
#Make sure you click inside the plot
        try:
            test = click.ydata-0.
            test = click.xdata-0.
            if self.seeking:
    #In case of error do not continue to try to ask about EW editing
                self.ewdit = False
                self.contx.append(click.xdata)
                self.conty.append(click.ydata)
                self.pcontx.append(click.x)
                self.pconty.append(click.y)
                self.a.scatter(click.xdata,click.ydata,marker='x',color='red',zorder=8,s=35)
                self.canvas.draw()
                if len(self.pconty) >= self.contpointmax:
                    self.seeking = False
                    self.cont_norm()
                
            if self.ewedit :
                self.contx.append(click.xdata)
                self.conty.append(click.ydata)
                self.mangaussfitlist.append(self.a.scatter(click.xdata,click.ydata,marker='x',color='red',zorder=8,s=35))
    #Commands to inform the user what is available during line fitting            
                try:
                    self.clickcommands.remove()
                except ValueError:
                    print 'Command not yet defined'
                if len(self.contx) == 1:
                    self.clickcommands = self.a.text(self.olines[self.editewlineindex],self.ymarker,'Click Line Center')
                if len(self.contx) == 2:
                    self.clickcommands = self.a.text(self.olines[self.editewlineindex],self.ymarker,'Left FWHM')
                if len(self.contx) == 3:
                    self.clickcommands = self.a.text(self.olines[self.editewlineindex],self.ymarker,'Right FWHM')
                if len(self.contx) == 4:
                    self.clickcommands = self.a.text(self.olines[self.editewlineindex],self.ymarker,'Right Continumm')
                self.canvas.draw()
    #Once all data points have been ented 
                if len(self.contx) == 5:
                    self.calc_gauss()
                    self.check_ew()
    #Remove the X markers from the gaussian fit
                    for i in self.mangaussfitlist:
                        i.set_visible(False)
                    self.canvas.draw()
    #Do simpsons rule integration
            if self.sint:
                self.contx.append(click.xdata)
                self.conty.append(click.ydata)
    
                self.intlist.append(self.a.scatter(click.xdata,click.ydata,marker='x',color='red',zorder=8,s=35))
                self.clickcommands.remove()
                if len(self.contx) == 1:
                    self.clickcommands = self.a.text(self.olines[self.editewlineindex],self.ymarker,'Right Side')
    #Do the integration
                if len(self.contx) == 2:
                    range, = np.where((self.wav0 >= np.min(self.contx)) & (self.wav0 <= np.max(self.contx)))
                    self.ewupdate = 1000.*simps(1.-self.spec0[range],x=self.wav0[range],even='avg')
                    self.canvas.draw()
    #check EW
                    self.check_ew()
                    self.sint = False
                    try:
                        self.fitlist[self.editewlineindex][0].remove()
                    except ValueError:
                        print 'line already removed'
                    self.fitlist[self.editewlineindex] = self.a.plot(self.wav0[range],self.spec0[range],color='red',zorder=8)
    #Remove the X markers from the gaussian fit
                    for i in self.intlist:
                        i.set_visible(False)
                self.canvas.draw()
#Throw error if clicked outside the plot
        except TypeError:
            self.error = 20
            self.onError()


#Function to check if the EW should be saved
    def check_ew(self): 
        if box.askyesno('Save','Should the New EW be added? \n '+str(round(self.ewupdate,2))):
            close = np.abs(self.orestlines[self.editewlineindex]-self.line)
            update, = np.where(close == close.min())
            self.lineew[update] = round(float(self.ewupdate),2)
            if self.sint:
                self.linebhm[update] = -9999.9
                self.linerhm[update] = -9999.9
                self.linedep[update] = -9999.9
                self.linesa1[update] = -9999.9
                self.linesa3[update] = -9999.9
                self.linesa13[update] = -9999.9
            
            else: 
                self.linebhm[update] = round(float(self.bhmupdate),8)
                self.linerhm[update] = round(float(self.rhmupdate),8)
                self.linedep[update] = round(float(self.depupdate),8)
                self.linecot[update] = round(float(self.cotupdate),8)
                self.linesa1[update] = round(float(self.sa1update),8)
                self.linesa3[update] = round(float(self.sa3update),8)
                self.linesa13[update] = round(float(self.sa13update),8)
#Plot the new EW
            self.ewtext[self.editewlineindex].remove()
            self.ewtext[self.editewlineindex] = self.a.text(rv_correct(self.line[update],self.rv-self.hrv),self.ewloc[update],str(self.lineew[update][0]),color='red',rotation='vertical')
            self.canvas.draw()
#After the ew has been calculated remind use to save if they quit 
            self.saved = False
        else:
            try:
                self.fitlist[self.editewlineindex][0].remove()
            except ValueError:
                print 'already removed'
#After the all the points have been clicked drop editing and calculate gaussian
        self.ewedit = False
#        self.spec_plot()
          
#Select points for continuum normalization
    def click_cont(self):
        if self.seeking == False:
            self.seeking = True
            self.contx  = []
            self.conty = []
#Pixel value so points can be deleted
            self.pcontx  = []
            self.pconty = []
            self.clickcommands = self.a.text(np.median(self.wav0),0.0,'Click on the Continuum then Click q to quit and fit \n d deletes bad points')
            self.canvas.draw()

#Actually normalize the continuum
    def cont_norm(self):
#First a second order legendre polynomail
#Changed to spline based on Gray Analysis of stellar photospheres
#        lcoef =  np.polynomial.legendre.legfit(self.contx,self.conty,len(self.contx)-1)

        if len(self.contx) > 3:
#sort the continuum in x
            sorted = np.argsort(np.array(self.contx))
#Spline interpolation of continuum points
#            fspline = interp1d(np.array(self.contx)[sorted],np.array(self.conty)[sorted],bounds_error=False)
#extroplate spline if not in fit (ext=0)
#iterate smooth of fit until  (s<= len(x))(Default)
#k degree of smoothing spline (default k is a free parameter used by s)
            fspline = UnivariateSpline(np.array(self.contx)[sorted],np.array(self.conty)[sorted],k=3,s=0,ext=0)


#Find the rms of the fit
            rms = np.sqrt(fspline.get_residual())
            rms = np.sqrt(np.sum(np.power(np.array(self.conty)[sorted]/fspline(np.array(self.contx)[sorted])-1.,2)))
            add = '_rms_cont'
#Put fit values into arrays to be added to file if accepted
            self.rms_vals[self.order-1,1] = np.min(self.wav0)
            self.rms_vals[self.order-1,2] = np.max(self.wav0)
            self.rms_vals[self.order-1,3] = rms
            self.rms_vals[self.order-1,4] = len(self.conty)

            addcont = '_cont_placement'

            
        

#Calc the continuum using spline
#        self.cont = np.polynomial.legendre.legval(self.wav0,lcoef)
            self.cont = fspline(self.wav0)
            contplot = self.a.plot(self.wav0,self.cont,color='red',zorder=16)
            self.canvas.draw()
#If the continuum normalization is good apply it        
            if box.askyesno("Continuum",'Is the Continuum Normalization Good?'):
#Add Continuum Point Selection uncertainty
               waves = np.array(self.contx)[sorted]
               contfile = open(self.filename.__str__().replace('.fits','')+addcont,'w')
               contfile.write('Order   Wavelength    ContCounts    Contrms     Npoints\n')
               for ll in waves:
                   finder, = np.where(np.abs(self.wav0-ll) < 0.5)
                   self.acounts.append(np.sum(self.spec0[finder]))
                   self.arootms.append(np.sqrt(np.sum(((self.spec0[finder]/self.cont[finder])-1.)**2)))
                   self.apoints.append(finder.size)
                   self.aorders.append(self.order)
                   self.awavels.append(ll)
               for ll in range(len(self.aorders)):
                   contfile.write('{0:^6d} {1:^10.2f} {2:^10.0f} {3:^10.6} {4:^10d}\n'.format(self.aorders[ll],self.awavels[ll],self.acounts[ll],self.arootms[ll],self.apoints[ll]))
               contfile.close()

                   

               try:
#2-D echelle format read in
                   if len(self.spec.data.shape) == 2:
                       self.spec.data[self.order-1] = self.spec0/self.cont
#1-D echelle format read in
                   if len(self.spec.data.shape) == 1:
                       self.spec.data = self.spec0/self.cont
#For ESO spectrum
               except AttributeError:
                   if self.eso_spec:
#                       self.spec[0][1] = self.spec0/self.cont   
                       self.spec0 = self.spec0/self.cont




#Write RMS information to file if you are satisfied with values
               rmsfile = open(self.filename.__str__()+add,'w')
               rmsfile.write(' Order   Startwav  Endwav    rmsfit     npoints\n')
               for i in np.arange(self.ordermax):
                   rmsfile.write('{4:^10.0f} {0:^10.2f} {1:^10.2f} {2:^10.8f} {3:^10.0f}\n'.format(self.rms_vals[i,1],self.rms_vals[i,2],self.rms_vals[i,3],self.rms_vals[i,4],i+1))
               rmsfile.close()


               self.spec_set()
               self.spec_plot()
            else:
#Ask to keep points in array if continuum normalization is improper
               if box.askyesno('Continuum','Keep Trying?'):
                   self.seeking = True
#remove pervious continuum
                   contplot[0].remove()        
               else:
                   self.spec_plot()
        else:
#Error if not enough points for spline fit
             self.error = 15
             self.seeking = True
             self.onError()
            


#Defind key strokes
    def on_key_event(self,event):

#Exits continuum normalization and fits
        if event.key == 'q':
            self.seeking = False
            self.cont_norm()

#Deletes continuum normalization point
        if event.key == 'd':
            if ((self.seeking) & (self.ewedit == False)) :
                subx = np.array(self.pcontx)-event.x
                suby = np.array(self.pconty)-event.y
                close = np.sqrt(subx**2.+suby**2.)
                delcol, = np.where(close == np.min(close))
                delcol = delcol[0]
#Delete index with closest pixel value to where d was pressed
                pointx = self.contx[delcol]
                pointy = self.conty[delcol]
                del self.conty[delcol]
                del self.contx[delcol] 
                del self.pconty[delcol]
                del self.pcontx[delcol]
                self.a.scatter(pointx,pointy,marker='x',color='blue',zorder=15,s=35)
                self.canvas.draw()

#key event for editing EW Measurement
        if event.key == 'e':
#If the EW for a given measure have been measured
#            if (((self.oewmeas[self.order-1]) | (self.ewmeas)) & (self.seeking == False)):
            if (self.seeking == False):
                subx = np.array(self.olines)-event.xdata
#Find the nearest point in wavelength space
                self.editewlineindex, = np.where(np.abs(subx) == np.min(np.abs(subx)))
#If two points are the same distance away error
                if self.editewlineindex.size > 1:
                    self.error = 12
                else:
                    self.edit_ew()
            else:
                self.error = 11
                self.onError()

#key event for interpolating Marcs Models
        if event.key == 'i':
            self.marcinterp()
#key event for calculation of Moog abundances
        if event.key == 'm':
            self.startmoog()


#Another way to trigger continuum normalization
        if event.key == 'c':
            self.click_cont()
#Keystroke to reset limits
        if event.key == 'r':
            self.reset_limits()

#Another way to add and subtract orders to be similar to IRAF
        if event.key == 'j':
            self.decreaseorder()
#Another way to add and subtract orders to be similar to IRAF
        if event.key == 'k':
            self.increaseorder()
#Automattically Fit lines
        if event.key == 'a':
            self.auto_gauss_fit()
#smooth by number in the smoothing box
        if event.key == 's':
            self.g_smooth()



#Edit EW measurement
    def edit_ew(self):
#set EW range to focus on line
        line = self.olines[self.editewlineindex]
        self.a.set_xlim([line-5.,line+5.])
        findmin, = np.where((self.wav0 >= line-5.) & (self.wav0 <= line+5.) & (np.isfinite(self.spec0)))
        ymin = np.min(self.spec0[findmin])
        self.a.set_ylim([ymin,1.2])
        self.ymarker = np.abs(np.diff(self.a.get_ylim()))*.1+np.min(self.a.get_ylim())
        self.canvas.draw()

        choice = ew_popup(root)
        root.wait_window(choice.top)
        edit = choice.editparm
#set gauss parameters
        if edit == 'sg':
            self.mangaussfitlist = []
            self.setew = True
            self.man_gauss_fit()
#manual gauss fit
        if edit == 'g':
#points for the gaussian fit
            self.mangaussfitlist = []
            self.man_gauss_fit()
#replace line EW with -9999.9
        if edit == 'd':
            minfind = np.abs(self.orestlines[self.editewlineindex]-self.line)
            remove, = np.where(minfind == minfind.min())
            self.lineew[remove] = -9999.9
            self.linecot[remove] = -9999.9
            self.linesa1[remove] = -9999.9
            self.linesa13[remove] = -9999.9
            self.linesa3[remove] = -9999.9
            self.linedep[remove] = -9999.9
            self.linerhm[remove] = -9999.9
            self.linebhm[remove] = -9999.9
#Plot the new EW
	    self.ewtext[self.editewlineindex].remove()
	    self.ewtext[self.editewlineindex] = self.a.text(rv_correct(self.line[remove],self.rv-self.hrv),self.ewloc[remove],str(self.lineew[remove][0]),color='red',rotation='vertical',fontsize=18)
#Remove poor fits
            try:
                self.fitlist[self.editewlineindex][0].remove()
            except ValueError:
                print 'line already removed'
            self.canvas.draw()
            
#Return to normal part of spectrum
#            self.reset_limits()
#Do straight integration
        if edit == 's':
            self.intlist = []
            self.man_int()

#Do not edit if c go straight to spec_plot
        if edit == 'c':
            print 'Do Nothing'
#            self.reset_limits()
#Do a Voigt fit Not yet implimented
        if edit == 'v':
            self.error = 9
            self.onError()
#Return to normal part of spectrum
#            self.reset_limits()


#Calculate a Voigt from manual fit
    def calc_voigt(self):
#First element is the continuum
        cont = (self.conty[0]+self.conty[-1])/2.
#Second element is line center
        depth = cont-self.conty[1]
#Width of feature (estimate of simga)
        width = np.abs(self.contx[2]-self.contx[3])/2.355

#Put human input into a first model
        v_model = models.Gaussian1D(amplitude=depth,mean=self.contx[1],stddev=width)+models.Lorentz1D(aplitude=depth,x_0=self.contx[1],fwhm=.1)
        fit_v = fitting.LevMarLSQFitter()
#Use input parameters to find LSQ fit to 3 sigma of line
        use, = np.where((self.wav0 >= g_model.mean-3*(width)) & (self.wav0 <= g_model.mean+3.*(width)))
        voigt_out = fit_v(v_model,self.wav0[use],cont-self.spec0[use])
       


#Calculate a Gaussian from manual fit
    def calc_gauss(self):
#First element is the continuum
        cont = (self.conty[0]+self.conty[-1])/2.
#Second element is line center
        depth = cont-self.conty[1]
#Width of feature (estimate of simga)
        width = np.abs(self.contx[2]-self.contx[3])/2.355

#Put human input into a first model
        g_model = models.Gaussian1D(amplitude=depth,mean=self.contx[1],stddev=width)
        fit_g = fitting.LevMarLSQFitter()
#Use input parameters to find LSQ fit to 3 sigma of line
        use, = np.where((self.wav0 >= g_model.mean-2*(width)) & (self.wav0 <= g_model.mean+2.*(width)))
        gauss_out = fit_g(g_model,self.wav0[use],cont-self.spec0[use])
#Get the covariance of the fit
        sa1 = np.sqrt(fit_g.fit_info['param_cov'][0,0])
        sa13 = np.sqrt(np.abs(fit_g.fit_info['param_cov'][0,2]))
        sa3 = np.sqrt(fit_g.fit_info['param_cov'][2,2])
#If you are manually setting the gaussian use the input values. (i.e. ignore the fit gaussian model)
#Then reset output paramater to fit 
        if self.setew:
            gauss_out = g_model
#if you are setting the FWHM manually set the varaince to a value perscribed by McWilliams 1995
            sa3 = 0.16
            self.setew = False
#Read EW in mA
        ew = 1000.*(gauss_out.amplitude.value*np.sqrt(2.*np.pi)*(np.abs(gauss_out.stddev.value)))
        rhm = 2.355*gauss_out.stddev.value/2.
        bhm = -rhm
        dep = gauss_out.amplitude.value
             

#Plot found EW
        sig3 = np.linspace(self.contx[1]-4.*gauss_out.stddev.value,self.contx[1]+4.*gauss_out.stddev.value,100)
#Put gaussian plots in a list so they can be deleted if unsatisfactory
#First remove old plot
#need the zero since for whatever reason plots are different than text
        try:
            self.fitlist[self.editewlineindex][0].remove()
        except ValueError:
            print 'No previous measure'
#Now plot new gaussian in place
        self.fitlist[self.editewlineindex] = self.a.plot(np.array(sig3),cont-gauss_out(sig3),color='red',zorder=3)
#Found gaussian mistake no longer need to use half width
#        findhalf = np.abs(gauss_out(sig3)-.5*depth)
#        fhalfmax, = np.where(findhalf == findhalf.min())
#        halfmax = sig3[fhalfmax[0]]-gauss_out.mean.value
        
        self.a.text(self.contx[1],.5,str(round(ew,1)),color='red',rotation='vertical')
        self.canvas.draw()
# the EW to pass to the plotting routine
        self.ewupdate = round(ew,1)
        self.rhmupdate = round(rhm,2)
        self.bhmupdate = round(bhm,2)
        self.depupdate = round(dep,2)
        self.cotupdate = round(cont,2)
        self.sa1update = round(sa1,2)
        self.sa13update = round(sa13,2)
        self.sa3update = round(sa3,2)

#Simpson's Rule integration
    def man_int(self):
        self.sint = True
        self.contx = []
        self.conty = []
        self.clickcommands = self.a.text(self.olines[self.editewlineindex],self.ymarker,'Click Left intregation point')
        self.canvas.draw()

#Command to manually fit Gaussian
    def man_gauss_fit(self):
        self.ewedit = True
        self.contx  = []
        self.conty = []
        self.clickcommands = self.a.text(self.olines[self.editewlineindex],self.ymarker,'Click Left Continuum')
        self.canvas.draw()

#Command to fit gaussians to lines in given order
    def auto_gauss_fit(self):
        if self.linelistin:
#set EW measured flag for order to be true
            self.oewmeas[self.order-1] = True

#List for gaussian plots
            self.fitlist = []

#Use the min and maximum wavlength values in order to search for lines
            maxw = max(self.wav0)
            minw = min(self.wav0)
            findlines = rv_correct(self.line,self.rv-self.hrv)
            use, = np.where((findlines >= minw) & (findlines <= maxw))
#array for newly calculated ew,hm,and depth
            newew = [] 
            bluehmlist = []
            redhmlist = []
            depthlist = []
            contlist = [] 
            sa1list = [] 
            sa13list = [] 
            sa3list = [] 
#holds plot values for temp. created EW values   
            tempew = []
            for i in use:
#            self.a.set_xlim([findlines[i]-5.,findlines[i]+5.])
#look for minimum value with some tolerance value of line center (initally setting to .1)
                lookbot, = np.where(np.abs(self.wav0-findlines[i]) <= self.findtol)
                if lookbot.size > 0:
                    minspec = self.spec0[lookbot].min()
                    lookbot, = np.where(self.spec0 == minspec)
#Only use line if minimum is found
                if lookbot.size > 0:
                    minwav = self.wav0[np.where(self.spec0 == minspec)][0]
#Find the line depth and half max
                    depthspec = np.float(1.-minspec)
                    halfspec = 1.-(depthspec*.5)
     
#find the blue line side information
#If more than one minimum pick the first
                    bluefind = lookbot[0]
                    listwav = []
                    listspec = []
                    distance = 0.
#look for the half max unless it is larger than 2 Angstroms in width
                 
                    while ((self.spec0[bluefind] <= halfspec) & (distance <= self.maxlinewidth) & (bluefind >= 0)):
                        bluefind = bluefind-1
                        listwav.append(np.float(self.wav0[bluefind]))
                        listspec.append(1.-self.spec0[bluefind])
                        distance = minwav-self.wav0[bluefind]
#reverse the blue wavelengths so they are in the correct order
                    listwav.reverse()
                    listspec.reverse() 
        
#append the center of the line to the list
                    listwav.append(np.float(minwav))
                    listspec.append(depthspec)
#find the red line side information            
                    redfind = lookbot[0]
                    distance = 0.
#look for the half max unless it is larger than 2 Angstroms in width
                    while ((self.spec0[redfind] <= halfspec) & (distance <= self.maxlinewidth) & (redfind < self.spec0.size-2)):
                        redfind = redfind+1
                        listwav.append(np.float(self.wav0[redfind]))
                        listspec.append(1.-self.spec0[redfind])
                        distance = self.wav0[redfind]-minwav
                   
#Only include lines with greater than 3 pixels in the line
#and the sigma is greater than 1 pixel seperation in wavelength
                    redhm = np.mean(listwav[-2:])
                    bluehm = np.mean(listwav[:2])
                    

#Put continuum in list to be writen as well as initial variables for the uncertainty in the fit
                    contlist.append(1.00)
                    sa1 = -9999.9
                    sa13 = -9999.9
                    sa3 = -9999.9

#convet half width estimate into sigma estimate
                    sigmaest = (redhm-bluehm)/(2.355)
#                    if ((len(listspec) > 2) & (sigmaest >= np.diff(self.wav0).min())):
                    if ((len(listspec) > 2)):
#Fit a gaussian model to EW
                        g_model = models.Gaussian1D(amplitude=depthspec,mean=minwav,stddev=sigmaest)
                        fit_g = fitting.LevMarLSQFitter()
#Use input parameters to find LSQ fit to 2 sigma of line
                        useg, = np.where((self.wav0 >= g_model.mean-2.*(sigmaest)) & (self.wav0 <= g_model.mean+2.*(sigmaest)))
#Testing to see if I can set gauss with out fitting better 11/11/15
                        gauss_out = g_model
#If there are not enough points to fit a gaussian skip the line
                        try:
                            gauss_out = fit_g(g_model,self.wav0[useg],1.-self.spec0[useg])
#                        gauss_out = fit_g(g_model,np.array(listwav),np.array(listspec))
                            sig3 = np.linspace(minwav-4.*gauss_out.stddev.value,minwav+4.*gauss_out.stddev.value,100)

#Use Half max relation to calc. EW
#Don't need to do this already fixed problem with EW integration
#                        findhalf = np.abs(gauss_out(sig3)-.5*depthspec)
#                        fhalfmax, = np.where(findhalf == findhalf.min())
#                        halfmax = sig3[fhalfmax[0]]-gauss_out.mean.value
#                        ew = 1000.*gauss_out.amplitude.value*1.06447*2.*np.abs(halfmax)
#Read EW in mA
                            ew = 1000.*gauss_out.amplitude.value*np.sqrt(2.*np.pi)*(np.abs(gauss_out.stddev.value))
                            depth = gauss_out.amplitude.value
                            minwav = gauss_out.mean.value
                            redhm = gauss_out.stddev*2.355/2.+minwav
                            bluehm = -gauss_out.stddev*2.355/2.+minwav
                            sa1 = np.sqrt(fit_g.fit_info['param_cov'][0,0])
                            sa13 = np.sqrt(np.abs(fit_g.fit_info['param_cov'][2,0]))
                            sa31 = np.sqrt(np.abs(fit_g.fit_info['param_cov'][0,2]))
                            sa3 = np.sqrt(fit_g.fit_info['param_cov'][2,2])
                        except TypeError:
                            ew = -9999.9
                        if (ew > 9.9) & (ew < 200.):
#Append new EW to list
                            newew.append(round(ew,1))
#Plot found EW 
                            self.fitlist.append(self.a.plot(np.array(sig3),1.-gauss_out(sig3),'--',color='teal',zorder=3,linewidth=4))
#                            self.fitlist.append(self.a.plot(np.array(sig3),1.-gauss_out(sig3),'--',color='white',zorder=4,linewidth=2))
#Update temp EW measurements
                            tempew.append(self.a.text(findlines[i],.6,str(round(ew,1)),color='red',rotation='vertical'))
                        else:
                            newew.append(round(-9999.9,1))
                            self.fitlist.append(self.a.plot(np.arange(2),np.arange(2),color='red',zorder=3))
#need the zero since for whatever reason plots are different than text
                            self.fitlist[-1][0].remove()
                            tempew.append(self.a.text(findlines[i],.5,str(-9999.9),color='red',rotation='vertical'))
                    else:
                        newew.append(-9999.9)
#Create the plot so everything remains in the same order but immediately drop the plot from the window
                        self.fitlist.append(self.a.plot(np.arange(2),np.arange(2),color='red',zorder=3))
#need the zero since for whatever reason plots are different than text
                        self.fitlist[-1][0].remove()
                        tempew.append(self.a.text(findlines[i],.2,str(-9999.9),color='red',rotation='vertical'))
#                    print listspec
#                    print listwav
#                    print 'No line found'
#draw found EW
                else:
                    newew.append(-9999.9)
#Create the plot so everything remains in the same order but immediately drop the plot from the window
                    self.fitlist.append(self.a.plot(np.arange(2),np.arange(2),color='red',zorder=3))
#need the zero since for whatever reason plots are different than text
                    self.fitlist[-1][0].remove()
                    tempew.append(self.a.text(findlines[i],.4,str(-9999.9),color='red',rotation='vertical'))
#Write the hm width to array
                redhmlist.append(redhm-minwav)
                bluehmlist.append(bluehm-minwav)
                depthlist.append(depthspec)
                sa1list.append(sa1)
                sa13list.append(sa13)
                sa3list.append(sa3)


            self.canvas.draw()

#Ask if Line EW should be updated
            if box.askyesno('Update','Should EW be Updated?'):
#Update lines in order
                for j,i in enumerate(use):
                    self.lineew[i] = newew[j]
                    self.linebhm[i] = bluehmlist[j]
                    self.linerhm[i] = redhmlist[j]
                    self.linedep[i] = depthlist[j]
                    self.linecot[i] = contlist[j]
                    self.linesa1[i] = sa1list[j]
                    self.linesa3[i] = sa3list[j]
                    self.linesa13[i] = sa13list[j]
#Plot the new EW
                    self.ewtext[j].remove()
                    self.ewtext[j] = self.a.text(findlines[i],self.ewloc[i],str(self.lineew[i]),color='red',rotation='vertical')
                self.canvas.draw()
#Remind the used the changes have not been saved if they try to exit
                self.saved = False
#                self.spec_plot()

        else:
            self.error = 8
            self.onError()
            self.canvas.draw()
#Remove red temp. EW measurements
        for i in tempew:
            i.remove()
      
#Command to increase the order to plot new spectrum
    def increaseorder(self):
        self.order = self.order+1
        if self.order > self.ordermax:
            self.order = 1 
        self.sorder.set(str(int(self.order)))
        self.a.clear()
        self.spec_set()
        self.spec_plot()

#Command to decrease order to plot new spectrum
    def decreaseorder(self):
        self.order = self.order-1
        if self.order < 1:
            self.order = self.ordermax
        self.sorder.set(str(int(self.order)))
        self.a.clear()
        self.spec_set()
        self.spec_plot()

#Create window in center of screen
    def centerWindow(self):
        w = 2000
        h = 1200
        sw = self.parent.winfo_screenwidth()
        sh = self.parent.winfo_screenheight()
       
        x = (sw-w)/2
        y = (sh-h)/2
        self.parent.geometry('%dx%d+%d+%d' % (w,h,x,y))

#Intialize the GUI
    def initUI(self):
	self.parent.title("Echelle Spectrum Analysis")


#create frame for plotting
        frame = Tk.Frame(self,relief=Tk.RAISED,borderwidth=1)
        frame.pack(fill=Tk.BOTH,expand=1)
 
        self.pack(fill=Tk.BOTH,expand=1)

#set up okay and quit buttons
        quitButton = Tk.Button(self,text="Quit",command=self.onExit)
        quitButton.pack(side=Tk.RIGHT,padx=5,pady=5)
        saveButton = Tk.Button(self,text="Save Fits",command=self.save_fits)
        saveButton.pack(side=Tk.RIGHT)

#set up velocity information to be added to line list
        rvText = Tk.StringVar()
        rvText.set("Radial Velocity (km/s)")
        rvDir = Tk.Label(self,textvariable=rvText,height=4)
        rvDir.pack(side=Tk.LEFT)
#Add so Velocity can be updated
        rv = Tk.StringVar()
        rv.set(str(round(self.rv,2)))
        self.rvval = Tk.Entry(self,textvariable=rv,width=10)
        self.rvval.bind("<Return>",self.rv_callback)
        self.rvval.pack(side=Tk.LEFT,padx=5,pady=5)
   
#set up heliocentric velocity information to be added to line list
        hrvText = Tk.StringVar()
        hrvText.set("Heliocentric Correction (km/s)")
        hrvDir = Tk.Label(self,textvariable=hrvText,height=4)
        hrvDir.pack(side=Tk.LEFT)
#Add so Velocity can be updated
        hrv = Tk.StringVar()
        hrv.set(str(round(self.hrv,2)))
        self.hrvval = Tk.Entry(self,textvariable=hrv,width=10)
        self.hrvval.bind("<Return>",self.rv_callback)
        self.hrvval.pack(side=Tk.LEFT,padx=5,pady=5)
   
#set up effective temperature
        teffText = Tk.StringVar()
        teffText.set("Teff (K)")
        teffDir = Tk.Label(self,textvariable=teffText,height=4)
        teffDir.pack(side=Tk.LEFT)
#Add so effective temperature can be updated
        teff = Tk.StringVar()
        teff.set(str(round(self.teff,2)))
        self.teffval = Tk.Entry(self,textvariable=teff,width=10)
        self.teffval.bind("<Return>",self.stellar_param)
        self.teffval.pack(side=Tk.LEFT,padx=5,pady=5)
   
#set up [Fe/H]
        fehText = Tk.StringVar()
        fehText.set("[Fe/H]")
        fehDir = Tk.Label(self,textvariable=fehText,height=4)
        fehDir.pack(side=Tk.LEFT)
#Add so [Fe/H] can be updated
        feh = Tk.StringVar()
        feh.set(str(round(self.feh,2)))
        self.fehval = Tk.Entry(self,textvariable=feh,width=10)
        self.fehval.bind("<Return>",self.stellar_param)
        self.fehval.pack(side=Tk.LEFT,padx=5,pady=5)
#set up log(g)
        loggText = Tk.StringVar()
        loggText.set("log(g)")
        loggDir = Tk.Label(self,textvariable=loggText,height=4)
        loggDir.pack(side=Tk.LEFT)
#Add so log(g) can be updated
        logg = Tk.StringVar()
        logg.set(str(round(self.logg,2)))
        self.loggval = Tk.Entry(self,textvariable=logg,width=10)
        self.loggval.bind("<Return>",self.stellar_param)
        self.loggval.pack(side=Tk.LEFT,padx=5,pady=5)
   
#set up microturbulent velocity
        vturbText = Tk.StringVar()
        vturbText.set("Vturb (km/s)")
        vturbDir = Tk.Label(self,textvariable=vturbText,height=4)
        vturbDir.pack(side=Tk.LEFT)
#Add so microturbulent velocity can be updated
        vturb = Tk.StringVar()
        vturb.set(str(round(self.vturb,2)))
        self.vturbval = Tk.Entry(self,textvariable=vturb,width=10)
        self.vturbval.bind("<Return>",self.stellar_param)
        self.vturbval.pack(side=Tk.LEFT,padx=5,pady=5)
#Add for smoothing 
        smoText = Tk.StringVar()
        smoText.set('s(FWHM)')
        smoDir = Tk.Label(self,textvariable=smoText,height=4)
        smoDir.pack(side=Tk.LEFT)
#Add so smoothing can be updated
        smo = Tk.StringVar()
        smo.set('{0:3.2f}'.format(self.smofwhm))
        self.sval = Tk.Entry(self,textvariable=smo,width=10)
        self.sval.bind("<Return>",self.smooth)
        self.sval.pack(side=Tk.LEFT,padx=5,pady=5)
       


#set up order number
        orderText = Tk.StringVar()
        orderText.set("Order")
        orderDir = Tk.Label(self,textvariable=orderText,height=4)
        orderDir.pack(side=Tk.LEFT)
#Add so order number can be updated
        self.sorder = Tk.StringVar()
        self.sorder.set(str(int(self.order)))
        self.orderval = Tk.Entry(self,textvariable=self.sorder,width=10)
        self.orderval.bind("<Return>",self.on_order_box)



        self.orderval.pack(side=Tk.LEFT,padx=5,pady=5)
   
   
   

#set up popup menu
        self.menu = Tk.Menu(self.parent,tearoff=0)
        self.menu.add_command(label="Beep",command=self.bell())
        self.menu.add_command(label="Exit",command=self.onExit)
        self.parent.bind("<Button-3>",self.showMenu)
        self.pack()

#set up Submenu
        menubar = Tk.Menu(self.parent)
        self.parent.config(menu=menubar)

        fileMenu = Tk.Menu(menubar)
        subMenu = Tk.Menu(fileMenu)
        subMenu.add_command(label="New File",command=self.open_file)
        fileMenu.add_cascade(label="Import",menu=subMenu,underline=0)

#create another item in menu
        fileMenu.add_separator()

        fileMenu.add_command(label='Exit',underline=0,command=self.onExit)
        menubar.add_cascade(label="File",underline=0,menu=fileMenu)

        self.open_file()



    def open_file(self):
#ask to open fits file
            self.filename = self.get_file()
    
#open the fits file
            self.fitfile = self.check_spec() 
#number of dimension of fits file
            self.dmax = len(self.fitfile)
    
            self.order = 1

            m = 0
#Set up parameters for POP ESO stuff           
            if self.filename.__str__().split('.')[-1] == 'tfits':
                self.ordermax = 1
                m = 1
                self.fitsd = 1
                self.eso_spec = True
#use the fits dimension with the spectrum in it (default is 0)
            while m == 0:
                try:
                    if self.eso_spec:
                        self.fitsd = 1
                        m = 1
                    else:
                        self.ordermax = np.shape(self.fitfile[self.fitsd])[0]
                        m = 1
                except IndexError:
#ask for fits dimension with ordered spectrum inside
                    self.error = 2
                    self.onError()
                    self.fitsd = self.get_fitsd()
                    m = 1
    
#set up boolean array with whether order has been measured or not
            self.oewmeas = np.ones(self.ordermax,dtype=bool)
#Set the default pyfits information to be extracted
            self.spec = self.fitfile[self.fitsd]

#Set up file to report continuum fit rms
            self.setup_rms()


    
#starts plotting the fits file
#First ask for order
#Changed to text box
#        self.order = self.on_order()
#set standard header 
            self.standard_header()
    
#Set up the infromation need to plot the spectrum
            self.spec_set()
#Draw said spectrum
            self.spec_plot()


#Look for file containing rms per order
    def setup_rms(self):
#setup file addendum
        add = '_rms_cont'
#create numpy array for the rms values 
#check to see if file exists 
        try:
            self.rms_vals = np.loadtxt(self.filename.__str__()+add,skiprows=1)
#If file doesn't exist create new array
        except IOError:
            self.rms_vals = np.zeros((self.ordermax,5))-9999.9



#RV Callback corrects line list for radial velocities
    def rv_callback(self,event):
#release cursor from entry box and back to the figure
#needs to be done otherwise key strokes will not work
        self.f.canvas._tkcanvas.focus_set()
        if self.linelistin:
            self.rv = float(self.rvval.get())
            self.hrv = float(self.hrvval.get())
            self.spec_plot()
            self.text_loc()
            self.plot_lines()
        else:
            self.error = 8
            self.onError() 

#Set Gaussian smoothing of spectrum
    def smooth(self,event):
        self.f.canvas._tkcanvas.focus_set()
        try:
            self.smofwhm = float(self.sval.get())
            self.g_smooth()
        except ValueError:
            self.error = 10
            self.onError()

#Start Gaussian smoothing
    def g_smooth(self):
       wavpix = np.mean(np.diff(self.wav0))
       g = Gaussian1DKernel(self.smofwhm/2.355/wavpix)
       self.spec0 = 1.-convolve(1.-self.spec0,g,boundary='extent')
       self.spec_plot()


    def stellar_param(self,event):
#release cursor from entry box and back to the figure
#needs to be done otherwise key strokes will not work
        self.f.canvas._tkcanvas.focus_set()
        try:
            self.teff = float(self.teffval.get())
            self.feh = float(self.fehval.get())
            self.logg = float(self.loggval.get())
            self.vturb = float(self.vturbval.get())
#error if not floats
        except ValueError:
            self.error = 10
            self.onError()

#save fits file
    def save_fits(self):
        self.fileout = Tkf.asksaveasfilename(filetypes=(('FITS File','*.fits'),('FITS File','*.FITS'),('FITS File','*.fit'),('FITS File','*.FIT')))
#cannot overwrite current file or core will dump
        if self.fileout == self.fitfile.filename():
            self.error = 14
            self.onError()
        else:
            try:
                self.spec.writeto(self.fileout,output_verify='fix+ignore')
            except IOError:
                if box.askyesno("Over write file?","Do you want to over write "+self.fileout+"?"):
                    self.spec.writeto(self.fileout,clobber=True,output_verify='fix+ignore')

#Function to get the dimension of 
    def get_fitsd(self):
        m = 0
        while m == 0:
            try:
                inputd = dimension_popup(root)
                root.wait_window(inputd.top)
                dim = int(inputd.order)
                    
                if ((dim > 0) & (dim <= self.dmax)):
                    m = 1
                else:
#Error order is out of range
                    self.error = 5
                    self.onError()
            except ValueError:
#IF order is work ESO assumed ESO formatted spectrum
                if inputd.order.upper() == 'ESO':
                   self.get_ESO()
                   dim = 0
                   m = 1

                else:
#Error order is not an integer
                    self.error = 4
                    self.onError()

        return dim



#Function for retrieving order from popup
    def on_order(self):
        m = 0
        while m == 0:
            try:
                inputO = order_popup(root)
                root.wait_window(inputO.top)
                order = int(inputO.order)
                if ((order > 0) & (order <= self.ordermax)):
                    m = 1
                else:
#Error order is out of range
                    self.error = 3
                    self.onError()
            except ValueError:
#Error order is not an integer
                self.error = 4
                self.onError()

        return order

#Function for retrieving order from text box
    def on_order_box(self,event):
#release cursor from entry box and back to the figure
#needs to be done otherwise key strokes will not work
        self.f.canvas._tkcanvas.focus_set()
        m = 0
        while m == 0:
            try:
                order = self.orderval.get()
                order = int(order)
                if ((order > 0) & (order <= self.ordermax)):
                    m = 1
                    self.order = order
                    self.spec_set()
                    self.spec_plot()
                else:
#Error order is out of range
                    self.error = 3
                    self.onError()
            except ValueError:
#Error order is not an integer
                self.error = 4
                self.onError()
            

#Shows a menu with right click
    def showMenu(self,e):
        self.menu.post(e.x_root,e.y_root)

#Exits the program
    def onExit(self):
        if self.saved == False:
            if box.askyesno('Warning','You have unsaved EW Data!\n Are you sure you want to Exit?'):
                self.quit()
            else:
                self.line_save()
        else:
           self.quit()

#Tells Why Order information is incorrect
    def onError(self):
        if self.error == 1:
            box.showerror("Error","File Not Found")
        if self.error == 2:
            box.showerror("Error","Fits Dimension Does not contain proper ordered Spectrum")
        if self.error == 5:
            box.showerror("Error","Fits Dimension is Out of range (0-"+str(self.dmax)+")")
        if self.error == 3:
            box.showerror("Error","Order Value is out of Range (1-"+str(self.ordermax)+')')
        if self.error == 4:
            box.showerror("Error","Value Must be an Integer")
        if self.error == 6:
            box.showerror("Error","File is not in Fits Format")
        if self.error == 7:
            box.showerror("Error","File is not formated for Moog")
        if self.error == 8:
            box.showerror("Error","No line list has been read")
        if self.error == 9:
            box.showerror("Error","Voigt Function is Currently in Developement")
        if self.error == 10:
            box.showerror("Error","Value Must be Float")
        if self.error == 11:
            box.showerror("Error","EW Have not been Measured (press a)")
        if self.error == 12:
            box.showerror("Error","Could not determine line desired for editing \n Please try again")
        if self.error == 14:
            box.showerror("Error","Cannot overwrite currently opened fits file \n Please choose a different file name")
        if self.error == 15:
            box.showerror("Error","Not enough points in continuum \n Please try again")
        if self.error == 20:
            box.showerror("Error","Must Select Inside Plot Bounds")
        if self.error == 100:
            box.showerror("Error","interp_marcs alias has not been defined\n please create a "+self.marcs_alias+" alais in ~/.cshrc")
        if self.error == 101:
            box.showerror("Error","moog alias has not been defined\n please create a moog alais in ~/.cshrc")
        if self.error == 107:
            box.showerror("Error","You need to save line list before running moog")
        if self.error == 108:
            box.showerror("Error","No Marcs model has been interpolated")
        if self.error == 110:
            box.showerror("Error","[Fe/H] is out of range of available models (too low)")
        if self.error == 111:
            box.showerror("Error","log(g) is out of range of available models (too low)")
        if self.error == 112:
            box.showerror("Error","Teff is out of range of available models (too low)")
        if self.error == 113:
            box.showerror("Error","[Fe/H] is out of range of available models (too high)")
        if self.error == 114:
            box.showerror("Error","log(g) is out of range of available models (too high)")
        if self.error == 115:
            box.showerror("Error","Teff is out of range of available models (too high)")

#routine for opening file
#check to make sure file is found
    def get_file(self):
        filename = Tkf.askopenfilename()
        while os.path.isfile(filename) == False:
            self.error = 1
            self.onError()
            filename = Tkf.askopenfilename()
            
        return filename

#make sure file is a fits file
    def check_spec(self):
        ffile = 0
        while ffile == 0:
            try:
                fitfile = pyfits.open(self.filename)
                ffile = 1
            except IOError:
                self.error = 6
                self.onError()
                self.filename = self.get_file()
        return fitfile

#Uses IRAF WAT Headers
    def standard_header(self):
        text = ''
        for i in self.spec.header["WAT2*"].itervalues():
            add = i
#Correct for any WAT Values less than 68 Characters
            while len(add) < 68:
                add = add+' '
            text = text+add
#Create an array containing one order of information per array
        text = text.replace('spe c','spec').replace('sp ec','spec').replace('s pec','spec')
        breakt = text.replace('multipspec ','').split('spec')[2:]
            
#loop over all orders creating a standard header which will be able to be read
        for i,j in enumerate(breakt):
            split = j.split(' = ')
            astr = split[0]
            while len(astr) < 2:
                astr = '0'+astr
            wavd = split[1].split(' ')
       
       
            try:
                lam0 = np.float(wavd[3])
                dlam = np.float(wavd[4])
       
       
                self.spec.header['CRVL1_'+astr] = lam0
                self.spec.header['CDLT1_'+astr] = dlam
            except ValueError:
                print 'Trying any way' 

#set variables spectrum of a given order
    def spec_set(self):
        popeso = 0
        ostring = str(self.order)
        while len(ostring) < 2:
            ostring = '0'+ostring
        try:
            self.dlam = self.spec.header['CDLT1_'+ostring]
        except KeyError:
            try:
               self.dlam = self.spec.header['CDELT'+str(self.order)]
            except KeyError:
               popeso = 1
        try:
            self.lam0 = self.spec.header['CRVL1_'+ostring]
        except KeyError:
            try:
               self.lam0 = self.spec.header['CRVAL'+str(self.order)]
            except KeyError:
               popeso = 1
               
#Use orders if 2-D spectrum
        if popeso == 0:
            if len(self.spec.data.shape) == 2:
                self.wav0 = np.arange(self.spec.data[self.order-1].size).astype(float)*self.dlam+self.lam0
                self.spec0 = self.spec.data[self.order-1]
#Use data if 1-D spectrum
            if len(self.spec.data.shape) == 1:
                self.wav0 = np.arange(self.spec.data.size).astype(float)*self.dlam+self.lam0
               
                self.spec0 = self.spec.data
        try:
            APTEST = self.spec.header['CTYPE'+str(self.order)] == 'LOG-LINEAR'
            if APTEST:
                self.wav0 = 10.**(self.wav0)
        except KeyError:
            print 'ASSUMING LINEAR WAVELENGTH SCALE'
#Specail to read Pop ESO spectrum
#        else:
#             self.get_ESO()
#            self.wav0 = self.spec.data['Wave']
#            self.spec0 = self.spec.data['Flux']
#            self.eso_spec = True
            
#Use archieve ESO spectrum
    def get_ESO(self):
        self.fitsd = 1
        self.spec = self.fitfile[self.fitsd].data
        self.spec0 = self.spec[0][1]
        self.wav0  = self.spec[0][0]
        self.order = 0
        self.ordermax = 0
        self.eso_spec = True

#actually plot the spectrum
    def spec_plot(self):
        self.a.clear()
        self.a.plot(self.wav0,self.spec0,color='black',zorder=1,linewidth=self.linewidth)
        self.a.set_xlabel(r'Wavelength (\AA)',fontsize=24)
        self.a.set_ylabel(r'Counts',fontsize=24)
        self.a.set_xlim([min(self.wav0),max(self.wav0)])
        fancy_plot(self.a)
        self.canvas.draw()
        if self.linelistin:
            self.text_loc()
            self.plot_lines()
    
#Convert inperpolate MARCS to MOOG FORMAT
#a whole lot of formating has to happen so while this is long it is not that interesting
    def marcs2moog(self):
        self.moogmodel = str(int(float(self.teff)))+'g'+str(round(float(self.logg),2))+'z'+str(round(float(self.feh),2))+'.alt'
        marcsmoogout = open(self.moogmodel,'w')
        marcsmoogout.write('GENERIC\n')
        marcsmoogout.write('{0:5.0f}/{1:4.2f}/{2:5.2f}\n'.format(float(self.teff),float(self.logg),float(self.feh)))
        marcsinterpout = open(str(int(float(self.teff)))+'g'+str(round(float(self.logg),2))+'z'+str(round(float(self.feh),2))+'.interpol','r')
        textin = []
        j = 0
        tau = [] 
        t = []
        pg = []
        pe = []
        akross = []
        m = 0
        while m == 0:
            text = marcsinterpout.readline()
            if j == 0:
                addlater = text[19:24]
                j = j+1
            else:
                if len(text) <= 70:
                    tau.append(text[1:8])
                    t.append(text[9:17])
                    pe.append(text[19:26])
                    pg.append(text[28:35])
                    akross.append(text[46:60])         
                    j = j+1
                else:
                    marcsmoogout.write('ntau=      '+str(j-1)+'\n')
                    m = 1
        marcsmoogout.write('{0:11.3f}\n'.format(float(addlater)))
        for i in range(j-1):
            marcsmoogout.write('{0:^10.4f}{1:^10.2f}{2:^10.4f}{3:^10.4f}{4:^10.4f}{5:^10.3e}\n'.format(float(tau[i]),float(t[i]),float(pg[i]),float(pe[i]),1.000,float(akross[i])))

        marcsmoogout.write('{0:13.3e}\n'.format(float(self.vturb)))
        marcsmoogout.write('NATOMS    {0:2d}  {1:6.2f}\n'.format(0,float(self.feh)))
        marcsmoogout.write('NMOL       {0:2d}\n'.format(16))
        marcsmoogout.write('{0:8.1f}{1:8.1f}{2:8.1f}{3:8.1f}{4:8.1f}{5:8.1f}{6:8.1f}{7:8.1f}\n'.format(101.,106.,107.,108.,606.,607.,608.,707.))
        marcsmoogout.write('{0:8.1f}{1:8.1f}{2:8.1f}{3:8.1f}{4:8.1f}{5:8.1f}{6:8.1f}{7:8.1f}\n'.format(708.,808.,10108.,60808.,109.,6.1,7.1,8.1))
        marcsmoogout.close()



#Write out abfind_pipe.par to be called by MOOG
    def write_abfind(self):
        abfindf = open('abfind_pipe.par','w')
        abfindf.write('abfind\n')
        abfindf.write('terminal     x11\n')
        abfindf.write('standard_out stan_chell_out\n')
        abfindf.write('summary_out  sum_chell_out\n')
        abfindf.write("model_in       '"+self.moogmodel+"'\n")
        abfindf.write("lines_in       '"+self.mooglinelist.name.replace(os.getcwd()+'/','')+"'\n")
        abfindf.write("atmosphere    2\n")
        abfindf.write("molecules     2\n")
#add curve of growth limits
#        abfindf.write("coglimits\n")
#minimum maximum and step of COG
#        abfindf.write("  -7.0    -4.7    0.01\n")
#        abfindf.write("plot          2\n")
        abfindf.write("lines         2\n")
        abfindf.write("strong        0\n")
        abfindf.write("flux/int      0\n")
        abfindf.write("damping       2\n")
        abfindf.write("trudamp       1\n")
        abfindf.write("obspectrum    5\n")

#Format the line list to be input into moog
    def in_moog_lines(self):
        self.mooglinelist = open(self.linelistout.name+'.moog','w')
        linedata = open(self.linelistout.name+'_linedata_sorted.dat','w')
        linedata.write('Measured by echelle analysis tool\n')
        linedata.write('    Wave       Ion       RHM       BHM      Depth     Cont       sigmad       sigmaHM   sigmadHM\n')
#Added header so EW can be read back in
        self.mooglinelist.write('EW Measured by echelle analysis tool\n')
#Sort by element
        sortedel = np.argsort(self.lineel)
        if self.checkvar.get() == 0:
#Write out line ew in MOOG format
            for i,j in enumerate(self.line):
                if float(self.lineew[sortedel][i]) > 0.:
                    self.mooglinelist.write('{0:^10.2f}{1:^10}{2:^10.2f}{3:^10.4f}{5}{4:^10.2f}\n'.format(self.line[sortedel][i],str(self.lineel[sortedel][i]),self.lineep[sortedel][i],self.linegf[sortedel][i],self.lineew[sortedel][i],' '*20))
                    linedata.write('{0:15.2f}{1:15.1f}{2:15.5f}{3:15.5f}{4:15.3f}{5:15.4f}{6:15.6f}{7:15.6f}{8:15.6f}\n'.format(self.line[sortedel][i],float(self.lineel[sortedel][i]),self.linerhm[sortedel][i],self.linebhm[sortedel][i],self.linedep[sortedel][i],self.linecot[sortedel][i],self.linesa1[sortedel][i],self.linesa3[sortedel][i],self.linesa13[sortedel][i]))
#Always write out in air wavelengths
        if self.checkvar.get() == 1:
#Write out line ew in MOOG format
            for i,j in enumerate(self.aline):
                if float(self.lineew[sortedel][i]) > 0.:
                    self.mooglinelist.write('{0:^10.2f}{1:^10}{2:^10.2f}{3:^10.4f}{5}{4:^10.2f}\n'.format(self.aline[sortedel][i],str(self.lineel[sortedel][i]),self.lineep[sortedel][i],self.linegf[sortedel][i],self.lineew[sortedel][i],' '*20))
                    linedata.write('{0:15.2f}{1:15.1f}{2:15.5f}{3:15.5f}{4:15.3f}{5:10.4f}{6:15.6f}{7:15.6f}{8:15.6f}\n'.format(self.aline[sortedel][i],float(self.lineel[sortedel][i]),self.linerhm[sortedel][i],self.linebhm[sortedel][i],self.linedep[sortedel][i],self.linecot[sortedel][i],self.linesa1[sortedel][i],self.linesa3[sortedel][i],self.linesa13[sortedel][i]))
#close the file
        self.mooglinelist.close()
        linedata.close()
 
    def plot_moog(self):
#        indlines,stat = read_moog.read_moog(self.mooglinelist.name.replace(os.getcwd()+'/',''))
        indlines,stat = read_moog.read_moog('sum_chell_out')
#Find the most square figure you can create
        nelements = stat['el'].size
        figcol = np.sqrt(nelements)
        if figcol-int(figcol) > 0.0001:
            figcol=figcol+1
#Convert to an integer that will be passed
        figcol = int(figcol)

#Lop off a row if it will be empty
        if figcol**2-nelements >= figcol:
            figrow = figcol-1
        else:
            figrow = figcol

#kill previous windows if open
        try:
            self.root2.destroy()
            self.root3.destroy()
            self.root4.destroy()
        except AttributeError:
            print 'Do nothing'
#Create the new windows for the popups
        self.root2 = Tk.Toplevel()
        self.root3 = Tk.Toplevel()
        self.root4 = Tk.Toplevel()
        abun_ep = abun_popup(self.root2,figcol,figrow)
        abun_rw = abun_popup(self.root3,figcol,figrow)
        abun_wa = abun_popup(self.root4,figcol,figrow)

#loop through to plot
        m = 0
        for i in range(figcol):
            for j in range(figcol):
                if m <= stat['el'].size-1:
#Find element 
                    useel, = np.where(indlines['el'] == stat['el'][m])
#Plot abundance vs. Ionization level 
                    abun_ep.abunax[i,j].scatter(indlines['ep'][useel],indlines['abun'][useel],color='black')
                    abun_ep.abunax[i,j].plot([indlines['ep'][useel].min(),indlines['ep'][useel].max()],np.zeros(2)+float(stat['ave'][m]),'--',color='black')
                    if float(stat['epb'][m]) > -10.:
                        abun_ep.abunax[i,j].plot([indlines['ep'][useel].min(),indlines['ep'][useel].max()],float(stat['epm'][m])*np.array([indlines['ep'][useel].min(),
                                                  indlines['ep'][useel].max()])+float(stat['epb'][m]),'--',color='red')
                    abun_ep.abunax[i,j].set_xlabel('EP (eV)')
                    abun_ep.abunax[i,j].set_ylabel('Abundance')
                     
                    self.add_text(abun_ep.abunax[i,j],str(stat['el'][m]),str(stat['epm'][m]),str(stat['epb'][m]),str(stat['ave'][m]),str(stat['std'][m]))
                    
#Plot abundance vs. REW  
                    abun_rw.abunax[i,j].scatter(indlines['rew'][useel],indlines['abun'][useel],color='black')
                    abun_rw.abunax[i,j].plot([indlines['rew'][useel].min(),indlines['rew'][useel].max()],np.zeros(2)+float(stat['ave'][m]),'--',color='black')
                    if float(stat['rwb'][m]) > -10.:
                        abun_rw.abunax[i,j].plot([indlines['rew'][useel].min(),indlines['rew'][useel].max()],float(stat['rwm'][m])*np.array([indlines['rew'][useel].min(),
                                                  indlines['rew'][useel].max()])+float(stat['rwb'][m]),'--',color='red')
                    abun_rw.abunax[i,j].set_xlabel('R.EW')
                    abun_rw.abunax[i,j].set_ylabel('Abundance')
                    self.add_text(abun_rw.abunax[i,j],str(stat['el'][m]),str(stat['rwm'][m]),str(stat['rwb'][m]),str(stat['ave'][m]),str(stat['std'][m]))
              
#Plot abundance vs. wavelength  
                    abun_wa.abunax[i,j].scatter(indlines['wav'][useel],indlines['abun'][useel],color='black')
                    abun_wa.abunax[i,j].plot([indlines['wav'][useel].min(),indlines['wav'][useel].max()],np.zeros(2)+float(stat['ave'][m]),'--',color='black')
                    abun_wa.abunax[i,j].set_xlabel('Wave')
                    abun_wa.abunax[i,j].set_ylabel('Abundance')
                    self.add_text(abun_wa.abunax[i,j],str(stat['el'][m]),'NaN','NaN',str(stat['ave'][m]),str(stat['std'][m]))
              


#add 1 to loop
                    m = m+1

#Draw the plots
        abun_ep.abuncanvas.draw()
        abun_rw.abuncanvas.draw()
        abun_wa.abuncanvas.draw()
            
    
#Add element text to plot
    def add_text(self,ax,text,textm,textb,textave,textsig):
        xs = np.array(ax.get_xlim())
        ys = np.array(ax.get_ylim())
        
        xdiff = xs.max()-xs.min()
        ydiff = ys.max()-ys.min()
        eltext = text.replace(' ','').split('.')
        ion = idict[eltext[1]]
        el = edict[eltext[0]]
        ax.text(xdiff*.1+xs.min(),ydiff*.9+ys.min(),el+' '+ion+' = '+textave+'+/-'+textsig,color='black',fontsize=18)
        ax.text(xdiff*.1+xs.min(),ydiff*.1+ys.min(),textm+'x+'+textb,color='black')


def main():
#Main loop
    global root
    root = Tk.Tk()
    app = gui_c(root)
    root.mainloop()

#directory label
#L1 = Tk.Label(root,text="Echelle Directory")
#L1.pack(side='left')
#E1 = Tk.Entry(root,bd=5)
#E1.pack(side='right')

if __name__=="__main__":
#create root frame
    main()
