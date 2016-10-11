# echelle_analysis_tool
A GUI written in python (TkInter), which calculates equivalent widths (EWs) for echelle spectra.
The tools works on all HIRES spectra reduced by MAKEE and spectra reduced via IRAF (i.e. makes assumptions about IRAF wavelength keywords). 
However, the GUI can do more than EWs.
If you install MOOG and the MARCs model interpolater (and create aliases) you are able to perform a full spectral analysis (i.e. determine Teff, log g, metal content, and microturbulence). 


I tried to make the program have minimal dependancies, so a typical anconda python install will contain everything nessacary to run this program.
Therefore, you may just download this program and run it by typing python echelle_analysis_tool.
However, I doubt you will keep all your spectra in the same directory as this program, so I would create an alias in your .cshrc file (e.g. alias eat "python /pathtopythonfile/echelle_analysis_tool.py").
In fact I highly recommend you create this alias, since EAT saves files in whatever directory you run it. 

If you want to use MOOG you need to create two additional aliases. One for the MARCS model intepolator (e.g. alias interp_marcs /pathto/marcinterpolator). You can download the interpolator supported by EAT at http://marcs.astro.uu.se/software.php. The other alias is for moog (e.g. alias moog /pathto/moog). You can download MOOG at http://www.as.utexas.edu/~chris/moog.html. The read_moog module included in the program requires the most recent MOOG build (June 2014), which outputs in a different format than previous verisons. 
