## import modules and environment
import os
import sys
sys.path.append("/home/jg2364/indexing/hklconv_new/")
os.chdir("/home/jg2364/indexing/hklconv_new/")

import numpy as np
import ctypes as ct
import os
import sys
import signal
import time
from nexusformat.nexus import *
import hkl
from scipy import *
import tqdm
import fabio

#some physical constants
PC=4.13566733e-15       # planks constant in si
c=2.99792458e8          # speed of light

nxsetmemory(2000)


stack_dir = "/data2/2023/chess_mar23/SiP/stacks/" 


#### This is the UB matrix obtained in UB matrix finding Jupyter Notebook
U = np.array([[-0.832761, -0.800328,  0.001698],
 [-0.800132 , 0.832608 , 0.023912],
 [ 0.017793 ,-0.016064 , 1.154748]])

image_dir = "/data2/2023/chess_mar23/SiP/"


outpath = "/data2/2023/chess_mar23/SiP/indexed_objects/"


#### Define HKL space dimensions and extent
H=np.arange(-4.8,4.8, 0.02) 
K=np.arange(-4.8,4.8, 0.02) 
L=np.arange(-4,5.5, 0.005)


#create the data storage arrays
data=np.zeros((len(H)*len(K)*len(L)),dtype=np.float32)
norm=np.zeros((len(H)*len(K)*len(L)),dtype=np.float32)
errors=np.zeros((len(H)*len(K)*len(L)),dtype=np.float32)


#### function initally written by J. Ruff, modified to use less RAM
def anglerock(H,K,L,stack, image_dir):
    dtype = np.dtype(np.uint16)
    
    imgmask = "/home/jg2364/chess_mar2023/calibration/mask_new.edf"
    imgmask=fabio.open(imgmask)
    scan_dir = image_dir
    imgfiles = [filename for filename in sorted(os.listdir(scan_dir)) if filename.endswith('cbf')]
    
    WL=stack.geo.wl*1e10#wavelength
    
    #counter for number of files processed
    count=0

    
    angs=np.linspace(0,0,6)
    
    tth=0.0
    eta=stack.psic.eta.nxdata
    chi=stack.psic.chi.nxdata
    phi=stack.data.phi.nxdata 
    phi = phi[:len(phi)-1]
    nu=0.0
    mu=0.0
    
    pol2=stack.geo.pol.nxdata
    az2=stack.geo.az.nxdata
    n=range(0,len(phi))

    for i in tqdm.tqdm(range(0,len(phi))):


        if (stack.norm.icnorm[i]>0.00):
            framenorm=1.0/(stack.norm.icnorm[i]*stack.norm.solidangle) #normalize solid angle and ionchamber
            count=count+1
            IN=hkl.Calc_HKL(pol2,az2,eta,mu,chi,phi[i],WL,U)
            imgcurrent=fabio.open(scan_dir+imgfiles[i]).data #open file 
            imgcurrent[imgmask.data>0.5]=-2
            hkl.HIST(IN,(np.array(imgcurrent)*framenorm).ravel(),1.0,H,K,L,data,norm,errors)


## dont print the output
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

# image directories
image_directory = [image_dir+'SiP_003/', image_dir+'SiP_004/', image_dir+'SiP_005/']
with HiddenPrints():
#### main loop, scan numbers 1 to 18, in 18 steps
    for scan_num in range(3, 6, 1):
        stack=nxload(stack_dir+'SiP_%d.nxs'%scan_num)
        W2=anglerock(H,K,L,stack, image_directory[scan_num-3])
        stack = 0


#### output the indexed object
outfile = "SiP_s1_3_15_1_hiRes"
fout=outpath+outfile+'.nxs'

dataout=data.clip(0.0)/norm.clip(0.9)

dataout=dataout.reshape(len(H),len(K),len(L))
H=H.astype('float32')
K=K.astype('float32')
L=L.astype('float32')
H=NXfield(H,name='H',long_name='H')
K=NXfield(K,name='K',long_name='K')
L=NXfield(L,name='L',long_name='L')
dataout=NXfield(dataout,name='counts',long_name='counts')

G=NXdata(dataout,(H,K,L))

while os.path.exists(fout):
    fout=fout[:-4]+"more.nxs"

G.save(fout)



