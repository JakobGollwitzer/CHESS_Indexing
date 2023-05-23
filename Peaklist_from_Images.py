#### import modules
import os
import sys


import numpy as np
import ctypes as ct
import os
import sys
import signal
import time
from nexusformat.nexus import *
import tqdm
import fabio


#### get input

print("All Directories end with a / at the end please.")


val = input("Stack Path?"+'\n')
stack_path = str(val)

val = input("Image Directory ?"+'\n')
image_dir = str(val)


val = input("Filepath to write peaklist to?"+'\n')
peaklistout = str(val)

val = input("Threshold? Recommnded is 1e5"+'\n')
threshold = float(val)


stack=nxload(stack_path)



# Disable
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


#### peak finding function
def findpeaks(stack, image_dir, intens_thresh):
    dtype = np.dtype(np.uint16)
    scan_dir=image_dir#str(stack.scan_dir.name)
    #imgmask = str(stack.mask_file.name)
    #imgmask=fabio.open(imgmask)
    imgfiles = [filename for filename in sorted(os.listdir(scan_dir)) if filename.endswith('.cbf')]
    #print(len(imgfiles))
    #counter for number of files processed
    count=0

    
    
    #in this case, we are rocking theta (eta in psi-C notation). Need to change for other angle scans
    tth=0.0
    eta=stack.psic.eta.nxdata
    chi=stack.psic.chi.nxdata
    phi=stack.data.phi.nxdata
    phi = phi[:len(phi)-1]
    nu=0.0
    mu=0.0
    #print(len(phi))
    pol2=stack.geo.pol.nxdata
    az2=stack.geo.az.nxdata
    n=range(0,len(phi))
    listofpeaks =  [] 
    for i in tqdm.tqdm(range(0,len(phi))):

        if (stack.norm.icnorm[i]>0.00):
            framenorm=1.0/(stack.norm.icnorm[i]*stack.norm.solidangle) #normalize solid angle and ionchamber
            count=count+1
            #blockPrint()
            with HiddenPrints():
                imgcurrent=fabio.open(scan_dir+imgfiles[i]).data #open file 
            #enablePrint()
            peaks = imgcurrent>intens_thresh
            peaklist = np.asarray(np.where(peaks)).T
            
            for j in range(len(peaklist)):
                listofpeaks.append(np.array([peaklist[j][0], peaklist[j][1], i, imgcurrent[peaklist[j][0]][peaklist[j][1]]]))
    
    listofpeaks = np.array(listofpeaks)
    

    return listofpeaks


### find the peaks within a threshold and write them to a file

peaklist = findpeaks(stack, image_dir, threshold)
print("Number of Peaks in the Generated List:")
print(len(peaklist))



import os
with open(peaklistout, 'wb') as f:
    np.save(f, peaklist)

print("wrote peaklist")




