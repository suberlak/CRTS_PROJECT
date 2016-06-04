# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 16:53:08 2015

@author: astronomy

Does the same as sf_yumi.py  EXCEPT  fitting anything - it just saves 
the mag_difference_squared,  the mjd_difference, and the mag_diff_error 
 ( made by adding the corresponding errors in quadrature )

Based on sf_yumi,  I read-in files  in the directory, and
plot the SF  : 
- for each file, I calculate the delta_mag vs delta_time  for each pair of 
  measurements
- From all those files together I bin the data into bins of with time_bin , 
  and for each bin I calculate rms from the <delta_mag>: this is SF
- I plot the SF vs delta_time 

UPDATE : 
03/18/2015
The difference between sf_load.py   and sf_load_NEW.py , is that instead
of making one big 'master' file, I decided to switch into making more 
smaller files, that may be much easier to handle 

03/23/2015
Added saving error_ij  (err_m_i and err_m_j  added in quadrature) , removed
storing <mag> and <mag_err> per LC, since this data is already stored as one 
number per LC in SDSS-CRTS cross-matched catalog 

MAJOR UPDATE MAJOR UPDATE MAJOR UPDATE MAJOR UPDATE MAJOR UPDATE MAJOR UPDATE 
MAJOR UPDATE MAJOR UPDATE MAJOR UPDATE MAJOR UPDATE MAJOR UPDATE MAJOR UPDATE

11/25/2015 
I made this program  based on sf_load_NEW.py  , to make individual master
files per lightcurve, which should make selection easier - instead of selecting 
lines from master files that correspond to objects satisfying the cut, 
we select master files for objects that satisfy the cut 

02/03/2016
A small update in line  75 , changing the stellar input catalog from 
'../stars_CRTS_proc_err_w_good/'  to '../stars_CRTS_LC_err_w/'
and   in line 94, so that instead of making SF for all the stars, 
it will only calculate it for those that are missing between the two 
directories   ( '../stars_CRTS_proc_err_w_good/'  is the one whereby 
an incorrect filter was applied to  '../stars_CRTS_processing_err_w/' 
'../stars_CRTS_LC_err_w/'  is the same as '../stars_CRTS_processing_err_w/' ) 

06/03/2016 
A major update : moving an entire project to a new location, making a new 
directory structure to make everything more transparent.   


NOTE : 
the program can be called   sf_load_NEW.py   inDir  outDir  outRoot  objType

"""

import os
import numpy as np 
import sys

def update_progress(progress):
    ''' A simple function updating the time progress. 
    
    progress : a value (float or int) between 0 and 100 indicating 
               percentage progress 
    '''
    print '\r[%-10s] %0.2f%%' % ('#' * int(progress/10), progress),

# Read in the LC files : 
# if the parameters are provided , from the user input
#args = sys.argv
#if len(args) > 1 : 
#    inDir = args[1]
#    outDir = args[2]
#    outRoot = args[3]
#    objType = args[4]
#    
#    if objType == 'qso'  : start = 4 ; end = -8
#    if objType == 'stars' : start = 4 ; end = -4
#    

# Read in the LC files  if there is no user input 
#if len(args) == 1 :


# CHOOSE HERE WHETHER WE SHOULD CRUNCH  STELLAR OR QSO  
# LIGHTCURVES TO MAKE MASTER FILES ! 
choice = 'qso'
if choice == 'stars':
    inDir =  '../proc_LC_CRTS/stars/'
    start= 4
    end=  -8
    
if choice == 'qso' :
    inDir = '../proc_LC_CRTS/qso/'  
    start= 4
    end=  -4
    

# regardless of choice, outDir stays the same 
outDir = '../data_products/sf_file_per_LC/' +choice+'/'

# Check if the directory for output exists, and if not, make one 
if not os.path.exists(outDir): os.system('mkdir %s' % outDir) 

# Load LC files 
inFiles = os.listdir(inDir)  

count = 0
total = float(len(inFiles))

# Loop over lightcurves 
for i in range(len(inFiles)): 
   
    percent = 100*(count / total)
    if (count % 10) == 0 : # every tenth loop.. 
        update_progress(int(percent))
    count += 1

    file = str(inFiles[i])

    # load the mjd, flux, and flux_error  from the read file 
    mjd,flx4,err = np.loadtxt(inDir+'%s' % (file),usecols=(0,1,2),unpack=True)
    
    # check for any nan's...
    add = np.where(np.isnan(flx4) == True)[0]
   
    # Delete all rows with nan's - there should be none in files I provide.. 
    mjd = np.delete(mjd,add); flx4 = np.delete(flx4,add); err = np.delete(err,add)
        
    # Sort according to mjd's 
    ind = np.argsort(mjd)
    mjd = mjd[ind]; flx4 = flx4[ind]; err = err[ind]
    
    # Calculate tau, mag difference and err difference for each pair (k,j), 
    # where tau=t_j-t_k (j>k) 

    delflx2  = []; delflx = [];  delflxerr = []; tau = []
    for j in range(len(mjd)-1):
        for k in range(j+1):     
            tau.append(mjd[j+1]-mjd[k])  # j from 1 and k<j
            noise2 = err[k]**2+err[j+1]**2 
            delflx.append((flx4[k]-flx4[j+1]))
            delflx2.append((flx4[k]-flx4[j+1])**2.0)
            delflxerr.append((noise2)**0.5)  
            
    # Change lists to arrays..         
    tau = np.array(tau); delflx = np.array(delflx)
    delflx2 = np.array(delflx2);   delflxerr = np.array(delflxerr) 
    
    # sort according to tau
    int0 = np.argsort(tau)
    tau = tau[int0]; delflx = delflx[int0]
    delflx2 = delflx2[int0]; delflxerr = delflxerr[int0]
    
    # grab the object name
    obj_name = file[start:end]
    
    # Make arrays that repeat the LC stats 
    # ID_arr = np.array(len(tau)*[file[start:end]])
    #avg_mag_arr = np.ones_like(tau, dtype=float) * np.mean(flx4)
    #avg_err_arr = np.ones_like(tau,dtype=float) * np.mean(err)    
    #avg_err_hold   = np.append(avg_err_hold, avg_err_arr)
    #avg_mag_hold   = np.append(avg_mag_hold, avg_mag_arr)
 
    ##### FILE SAVING #######
    DATA = np.column_stack((delflx,tau, delflxerr))    
    outfile = outDir + 'SF_' + obj_name + '.txt'
    np.savetxt(outfile, DATA, delimiter =' ', fmt="%s") 

    #sys.stdout.write('\nProcessing LC %d out of %d \r' % (i, len(inFiles)))
    #sys.stdout.flush()
    #sys.stdout.write('\nSaving... %s \n'% outfile)
    #sys.stdout.flush()
       
print '\nAll in all, we saved %d SF files ' % (count)  
