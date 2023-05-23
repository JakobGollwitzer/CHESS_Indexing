#### import modules
import os
import sys

import numpy as np

import pyFAI
import fabio
import os
from nexusformat.nexus import *
from spec2nexus.spec import SpecDataFile
import inspect
import tqdm


### get input
print("All Directories end with a / at the end please.")


val = input("PONI Filepath?"+'\n')
poni_file = str(val)

val = input("Mask Filepath?"+'\n')
mask_file = "/home/jg2364/chess_oct2022/calibration/mask_16keV.edf"

mask_file = str(val)
val = input("Sample Name?"+'\n')

sample_name = str(val)
val = input("Path of directory to save in?"+'\n')

outpath = str(val)
val = input("Specfile path?"+'\n')

specfile = str(val)


val = input("enter list of scans to stack as e.g. 1,2,3,4"+'\n')
val = np.array(val.split(','))




### make the stack, partly written by J. Ruff
def make_stack(scan_num):
    ai=pyFAI.load(poni_file)
    imgmask=fabio.open(mask_file)

  

    #read from spec
    spexus = SpecDataFile(specfile)
    
    
    scanexus = spexus.getScan(scan_num)
    
    
    phi=np.asarray(scanexus.data['phi'])
    chi=float(scanexus.positioner['chi'])
    mu=float(scanexus.positioner['mu'])
    eta=float(scanexus.positioner['th'])
    icnorm=np.asarray(scanexus.data['ic2'])
    icnorm=icnorm/(np.average(icnorm))

    # read in geometry info from pyFAI
    pol=ai.twoThetaArray()*np.cos(ai.chiArray()+(np.pi*0.5))
    az=-ai.twoThetaArray()*np.sin(ai.chiArray()+(np.pi*0.5))
    psi=ai.chiArray()
    qmag=ai.qArray()
    wl=ai.wavelength #assumes that the calibrant and data energy is the same

    imgshape=imgmask.shape

    solidangle=ai.solidAngleArray()

    ##Define the nexus fields and object to be created

    phi=NXfield(phi, name='phi')
    xpixel=NXfield(range(0,imgshape[0]), name='x')
    ypixel=NXfield(range(0,imgshape[1]), name='y')
    pol=NXfield(pol, name='pol')
    az=NXfield(az, name='az')
    qmag=NXfield(qmag, name='qmag')
    psi=NXfield(psi, name='psi')
    icnorm=NXfield(icnorm, name='icnorm')
    solidangle=NXfield(solidangle, name='solidangle')

    W=NXroot()

    #W.scan_dir=NXentry()
    #W.scan_dir.name=scan_dir
    
    W.mask_file=NXentry()
    W.mask_file.name=mask_file

    W.psic=NXentry() #psi-circle angles from spec
    W.psic.eta=eta
    W.psic.mu=mu
    W.psic.chi=chi

    W.sample=NXentry() #sample info, read from spec / meta / command line
    W.sample.compound=NXfield(specfile,name="compound")
    W.sample.id=NXfield("",name="sampleid")
    W.sample.temperature=NXfield(300,name="temperature")


    W.geo=NXentry() #area detector geometry
    W.geo.pol=pol
    W.geo.az=az
    W.geo.psi=psi
    W.geo.qmag=qmag
    W.geo.wl=NXfield(wl,name='wavelength')

    W.norm=NXentry()
    W.norm.icnorm=icnorm
    W.norm.solidangle=solidangle

    W.powder=NXentry() #mimic powder pattern
    W.data=NXentry()
    counts=0
    W.data=NXdata(counts,(xpixel,ypixel,phi))
    
    outfile=sample_name+"_"+str(scan_num)+".nxs" 
    fout = outpath+outfile
    
    W.save(fout, 'w')
    print("wrote stack for scan: "+str(scan_num))



#### main loop, stack scan numbers 1 to 18, in 18 steps

for scan_num in range(0, len(val), 1):
    make_stack(int(val[scan_num]))
