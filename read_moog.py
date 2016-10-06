import numpy as np
from astropy.table import Table
from astropy.io import ascii

#Asplund 2009 Photosphere Comp.
aphot = ascii.read('/wandajune_home/jakub/MyModules/asplund_sun.dat')

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


def read_moog(filein):

    filedat = open(filein,'r')

#stat array stuff
    elnam = []
    elave = []
    elstd = []
    elnum = []
    elrwm = []
    elepm = []
    elrwb = []
    elepb = []
    elrwr = []
    elepr = []
   
#big array stuff
    bigwav = []
    bigabu = []
    bigew = []
    bigrew = []
    bigep = []
    bigel = []
    biggf = []
    for j,i in enumerate(filedat):
        if j >= 3:
            try:
                bigwav.append(float(i[:10]))
                bigrew.append(float(i[49:56]))
                bigew.append(float(i[40:46]))
                bigabu.append(float(i[60:66]))
                bigep.append(float(i[23:29]))
                bigel.append(round(float(i[12:17]),1))
                biggf.append(float(i[31:40]))
            except ValueError:
                if i[:17] == 'average abundance':
                    elnam.append(bigel[-1])
                    elave.append(i[20:26])
                    elstd.append(i[49:54])
                    elnum.append(i[68:71])
                if i[:35] == ' No statistics done for E.P. trends':
                    elepm.append(-9999.9)
                    elepb.append(-9999.9)
                    elepr.append(-9999.9)
                if i[:35] == ' No statistics done for R.W. trends':
                    elrwm.append(-9999.9)
                    elrwb.append(-9999.9)
                    elrwr.append(-9999.9)
                if i[:16] == 'E.P. correlation':
                    elepm.append(float(i[26:34]))
                    elepb.append(float(i[49:55]))
                    elepr.append(float(i[71:]))
                if i[:16] == 'R.W. correlation':
                    elrwm.append(float(i[26:34]))
                    elrwb.append(float(i[49:55]))
                    elrwr.append(float(i[71:]))


#stats
    elnam = np.array(elnam)
    elave = np.array(elave)
    elnum = np.array(elnum)
    elstd = np.array(elstd)
    elrwm = np.array(elrwm)
    elrwb = np.array(elrwb)
    elrwr = np.array(elrwr)
    elepm = np.array(elepm)
    elepb = np.array(elepb)
    elepr = np.array(elepr)
#Ind. values
    bigwav = np.array(bigwav)
    bigabu = np.array(bigabu)
    bigew =  np.array(bigew)
    bigrew = np.array(bigrew)
    bigep =  np.array(bigep)
    bigel =  np.array(bigel)
    biggf =  np.array(biggf)

    bigarray = Table([bigwav,bigabu,bigew,bigrew,bigep,bigel,biggf],names=['wav','abun','ew','rew','ep','el','loggf'])
    statarray = Table([elnam,elave,elnum,elstd,elepm,elepb,elrwm,elrwb,elepr,elrwr],names=['el','ave','num','std','epm','epb','rwm','rwb','rep','rrw'])

    return bigarray, statarray

def read_synth(filein):

    p =0
    k =0
    filedat = open(filein,'r')

    flux = []
    for j,i in enumerate(filedat):
        if k == 1:
            fd = i.split(' ')
            for l,q in enumerate(fd):
              try:
                 fluxer = float(q) 
                 flux.append(fluxer)
              except ValueError:
                 continue
        if i[:3] == 'ALL':
            k = 0
            if p == 1:
                aflux = np.array(flux).ravel()
            if p > 1:
#                while len(flux) > aflux.shape[1]:
#                    flux = flux[:-1]
#                while len(flux) < aflux.shape[1]:
#                    flux.append(1.0)
                uflux = np.array(flux).ravel()
                aflux = np.vstack((aflux,np.array(flux).ravel()))
            flux = []


        if i[:7] == 'element':
            print 'do nothing'
        if i[:5] == 'MODEL':
            print 'do nothing'
        if i[:2] == '  ':
            print 'New model'
            wavestart = float(i[2:11])
            wavestep = float(i[28:33])
            k = 1
            p = p+1
    filedat.close()
    if p == 1:
        aflux = np.array(flux).ravel()
    if p > 1:
        while len(flux) > aflux.shape[1]:
            flux = flux[:-1]
        while len(flux) < aflux.shape[1]:
            flux.append(1.0)
        aflux = np.vstack((aflux,np.array(flux).ravel()))
    return aflux,wavestart,wavestep


def read_smooth(files):

    fil = open(files,'r')

    p = 0
    k = 0

    flux = []
    wave = []

    for j,i in enumerate(fil):
        if k == 1:
            try:
                wave.append(float(i[:10]))
                flux.append(float(i[14:]))
            except ValueError:
                print 'skipped'
                ll = 55
        if i[:3] == 'the':
            print 'end'
            k = 0
            if p == 1:
                aflux = np.array(flux)
                awave = np.array(wave)
            if p > 1:
                uflux = np.array(flux)
                uwave = np.array(wave)
                aflux = np.vstack((aflux,uflux))
                awave = np.vstack((awave,uwave)) 
        if i[:5] == 'start':
            print 'restart'
            flux = []
            wave = []
            k = 1
            p = p+1

    fil.close()
    if p == 1:
        aflux = np.array(flux)
        awave = np.array(wave)
    if p > 1:
        aflux = np.vstack((aflux,np.array(flux)))
        awave = np.vstack((awave,np.array(wave)))
    return awave,aflux


def read_ewfind(files,alld=False):

    fil = open(files,'r')

    p = 0
    k = 0
    wav = []
    ews = []
    els = []
    lgf = []
    eps = []
    seek = False
    for j,i in enumerate(fil):

        while seek:
            wav.append(float(i[:10]))
            ews.append(float(i[54:]))
            els.append(i[32:42])
            if alld:
                eps.append(float(i[15:23]))
                lgf.append(float(i[25:33]))


            seek = False
        if i[:10] == 'wavelength':
            seek = True

    if alld:
        ews, wav, els, eps, lgf = np.array(ews), np.array(wav), np.array(els), np.array(eps), np.array(lgf)
        return wav,els,eps,lgf,ews

    else:
        ews, wav, els = np.array(ews), np.array(wav), np.array(els)
        return wav,els,ews




