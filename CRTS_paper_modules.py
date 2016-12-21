# -*- coding: iso-8859-1 -*-
'''
Some routines that are common to my paper-making programs (those that make 
Figures 2,3,4  : 

get_qso_catalog   : read in cross-matched CRTS-SDSS catalog of qso 
get_stars_catalog : read in cross-matched CRTS-SDSS catalog of stars
cut_qso : make selection of qso from the catalog 
cut_stars : make selection of stars from the catalog 

add_tau_delflx
read_xi_ei 


MAJOR UPDATE : 2016/06/03 : moved everything to a new directory structure , to

2016/06/07 : modified  read_xi_ei,  to allow it not to read in red stars  : 
             moved the red stars IDs argument to the very end,  as an optional 
             kw.  This helps  rework Fig.3 more efficiently. Fig.2, Fig.4 code
             affected, since the order of things as they are read in has to be 
             changed,  to inDirStars, good_ids_S_blue, inDirQSO, good_ids_QSO, 
                          good_ids_S_red

2016/06/14 : modified cut_qso , to allow for quasar redshift to be returned 
             for the chosen quasars. This helps to correct the timescale 
	     to the restframe. 
'''
import numpy as np
import os
import sys 

def update_progress(progress):
    ''' A simple function updating the time progress. 
    
    progress : a value (float or int) between 0 and 100 indicating 
               percentage progress 
    '''
    sys.stdout.write('\r[%-10s] %0.2f%%' % ('#' * int(progress/10), progress),)
    sys.stdout.flush()

def get_qso_catalog():
    File = '../data_products/CRTS_SDSS_cross_matched_qso_DB_QSO_catalog.txt'
    colnames = open(File,'r').read().splitlines()[0][1:].split()
    datatable = np.genfromtxt(File, dtype=str)
    qso_catalog = {}
    print('Zipping CRTS-SDSS quasars catalog from %s' % File)
    for label, column in zip(colnames, datatable.T):
        qso_catalog[label] = column
    
    print('Read in %d quasars from CRTS' % len(qso_catalog['redshift']))
    return  colnames, qso_catalog
    
def get_stars_catalog():
    File = '../data_products/CRTS_SDSS_cross_matched_stars_catalog.txt'
    colnames = open(File,'r').read().splitlines()[0][1:].split()
    datatable = np.genfromtxt(File)
    stars_catalog = {}
    print('zipping CRTS-SDSS stars catalog...')
    for label, column in zip(colnames, datatable.T):
        stars_catalog[label] = column
    
    print('Read in catalog for %d stars from CRTS '% len(stars_catalog['g_mMed']))
    return  colnames, stars_catalog

# Perform cuts 
def cut_qso(qso_cat=None, mMin=-9, mMax=19,   
            mErrMin = -9, mErrMax = 0.3,cut_mag='r', redshift = None, match_deg_rad = 1.0 / 3600):
    ''' A short  routine to select CRTS quasars according to desired parameters.
    
    Parameters:
    ----------
    qso_cat  - a CRTS-SDSS matched catalog, the output of get_qso_catalog()
    mMin / mMax - minimum / maximum  desired magnitude , of SDSS filter specified 
              by cut_mag  
    mErrMin - minimum / maximum  desired average CRTS lightcurve error, after 
              day-averaging 
    cut_mag - desired SDSS filter for the cut,  allowed values  are u,g,r,i,z. 
              Default setting is 'r' 
    match_deg_rad  - a maximum matching radius allowed  - the CRTS-SDSS matched 
              catalog stores the matching radius in degrees. By default, we 
              choose 1 arcsec , i.e. 1/3600 of a degree as an appropriate 
              matching radius for QSO.  For the CRTS sample, only 15 do not 
              have an SDSS counterpart within that range. 

    Returns:
    --------
    qso_id : the CRTS id's of the quasars accepted  through the selection criteria 
    '''
    print('Returning only QSO which had an SDSS counterpart within %f radians'%match_deg_rad)
    mask_rad = (qso_cat['m_ang_deg'].astype(float) < match_deg_rad)
    mask_mag = (qso_cat[cut_mag].astype(float) > mMin) * (qso_cat[cut_mag].astype(float) < mMax) 
    mask_err = (qso_cat['CRTS_avg_e'].astype(float) > mErrMin) * (qso_cat['CRTS_avg_e'].astype(float) < mErrMax)
    mask = mask_rad * mask_mag * mask_err 
    qso_id = qso_cat['CRTS_id'][mask]
    print('\n These cuts reduced the number of qso  in the sample from %d to %d' %
        (len(qso_cat['redshift']), len(qso_id)))
          

    if redshift is not None:
        print('Also returning quasar redshifts...')
        qso_redshift = np.array(qso_cat['redshift'][mask]).astype(float)
        return qso_id, qso_redshift 
    else:
        return  qso_id

def cut_stars(star_cat=None, mMin=-9, mMax=19, mErrMin = -9, 
              mErrMax = 0.3, gi_Min = -1, gi_Max=1 , cut_mag='r_mMed'):
    ''' A short  routine to select CRTS stars according to desired parameters.
    
    Parameters:
    ----------
    star_cat  - a CRTS-SDSS matched catalog, the output of get_star_catalog()
    mMin / mMax - minimum / maximum  desired magnitude , of SDSS filter specified 
              by cut_mag  
    mErrMin - minimum / maximum  desired average CRTS lightcurve error, after 
              day-averaging 
    gi_Min / gi_Max - a tool to select CRTS stars based on SDSS photometry : the 
              minimum / maximum value of the g-i color to allow 
    cut_mag - desired SDSS filter for the cut,  allowed values  are 'u_mMed',
              'g_mMed','r_mMed','i_mMed','z_mMed'.Default setting is 'r_mMed' 
  
    Returns:
    --------
    qso_id : the CRTS id's of the quasars accepted  through the selection criteria 
    '''
    mask_mag = (star_cat[cut_mag] > mMin) * (star_cat[cut_mag] < mMax) 
    mask_err = (star_cat['CRTS_Merr'] > mErrMin) * (star_cat['CRTS_Merr'] < mErrMax)
    SDSS_gi = star_cat['g_mMed'] - star_cat['i_mMed']
    mask_color = (SDSS_gi > gi_Min ) * (SDSS_gi < gi_Max)
    mask = mask_mag * mask_err * mask_color
    star_id_f = star_cat['crts_id'][mask]
 
    # convert floats to strings without comma and zeros
    star_id = np.array(["{:.0f}".format(name) for name in star_id_f])
    print('\n These cuts reduced the number of stars  in the sample from %d to %d'%
        (len(star_cat['CRTS_M']),  len(star_id)))
    return  star_id


# inside the main loop : get tau, delflx from a master file, either qso or star
def add_tau_delflx(File, inDir, data, z=None):
    ''' A function used by read_xi_ei() method to add delta_mag, delta_time, 
    and corresponding error err from an individual master file  
    to the list.

    Parameters:
    -----------
    File : a name of a 'master file' with delta_time, delta_mag, err. 
          It is assumed that the name coresponds to the  original_CRTS_name in 
          a following way:   out_original_CRTS_name.dat.txt

    inDir : a directory where the master file should be found 
    data : a storage array with four fields, passed on from read_xi_ei() 
    
    Returns:
    ---------
    delflx, tau, err, master_names : updated elements of data array, 
         to which we appended  the values of delflx, tau, err, and a 
         name of the object (master file)  multiplied by the number of 
         rows in the master file 
 
    '''
    # read in storage arrays
    delflx = data[0]  
    tau = data[1]
    err = data[2]
    master_acc_list = data[3]   
    
    # grab the object name 
    master_name = File[3:-4]
    
    # read in the i-th master file 
    master =  np.genfromtxt(inDir+File, dtype=str)
    
    # read in tau,  del_mag,  del_mag_err for quasars on the list 
    delflx = np.append(delflx, master[:,0].astype(float))
    if z is not None:
        tau = np.append(tau, master[:,1].astype(float) / (1.0+z))
    else:
        tau = np.append(tau, master[:,1].astype(float))
    err = np.append(err, master[:,2].astype(float))
    master_names  = np.append(master_acc_list, np.array(len(master[:,0])*[master_name]))
    
    return delflx, tau, err, master_names
    
def read_xi_ei(inDirStars, good_ids_S_blue, inDirQSO,
                 good_ids_QSO, good_ids_S_red=None, redshift=None):
    ''' A routine to read the delta_mag (xi), delta_time (tau), and error (ei) 
    for CRTS stars and quasars. Stars and quasar master files are read from 
    inDirStars and inDirQSO , and only those files are selected to be read in 
    that are on the list of good_ids_S_blue, good_ids_S_red  for blue 
    and red stars, and good_ids_QSO  for quasars.  
    
    '''           
    # use new names to make things code-compliant...      
    inDir_S = inDirStars
    good_ids_S_blue    = good_ids_S_blue
    if good_ids_S_red is not None:
       good_ids_S_red    = good_ids_S_red
    inDir_Q       = inDirQSO
      
    # Read the Stellar Master file names 
    masterFiles_S = os.listdir(inDir_S)  # need to shorten these names ... 
    masterFilesS1 = [name[3:-4] for name in masterFiles_S]
    
    # pick out Blue and Red stars, based on their names ... 
    good_masterSB = np.array(masterFiles_S)[np.in1d(masterFilesS1, good_ids_S_blue)]
    if good_ids_S_red is not None:
        good_masterSR = np.array(masterFiles_S)[np.in1d(masterFilesS1, good_ids_S_red)]
    
    # Read the QSO Master file names 
    #masterFiles_Q = os.listdir(inDir_Q)
    #masterFilesQ1 = [name[3:-4] for name in masterFiles_Q]
    #good_masterQ =  np.array(good_ids_QSO) #np.array(masterFiles_Q)[np.in1d(masterFilesQ1, good_ids_QSO)]
    good_masterQ = np.array(['SF_' +qso+'.txt' for qso in good_ids_QSO])

    # If no previous read-in xi, ei exists, initialize arrays    
    print('making new delflx, tau, xi arrays')
    delflx_S      = np.empty(0,dtype=float)
    tau_S         = np.empty(0,dtype=float)
    err_S         = np.empty(0,dtype=float)
    master_acc_list_S = np.empty(0, dtype=str)

    delflx_Q      = np.empty(0,dtype=float)
    tau_Q         = np.empty(0,dtype=float)
    err_Q         = np.empty(0,dtype=float)
    master_acc_list_Q = np.empty(0, dtype=str)

    # Initialize the data structures to which more and more delta_t and delta_mag
    # are addded from each consecutive master file 
    qso_data = [delflx_Q, tau_Q, err_Q, master_acc_list_Q] 
    star_data_blue = [delflx_S, tau_S, err_S, master_acc_list_S]
    
    if good_ids_S_red is not None:
        star_data_red  = [delflx_S, tau_S, err_S, master_acc_list_S]
    
    ### READ IN QUASARS ### 
    
    print('Reading in quasars...')
    if redshift is not None: 
      # correcting QSO delta_time to restframe
      print('\n')
      print('Correcting delta_time to restframe, t_rest = t_obs / (1+z)')
      c = 0
      for i in range(len(good_masterQ)): #  len(masterFiles_Q)
          z = redshift[i]
          File = good_masterQ[i]
          qso_data = add_tau_delflx(File,inDir_Q, qso_data, z)
          c += 1 
          if c % 5 == 0:
            progress = (100.0*c) / float(len(good_masterQ))
            update_progress(progress)
            #print('\r[%-10s] %0.2f%%' % ('#' * int(progress/10), progress),)
              #print('\r----- Already read %d%% of qso'%pers),
    else:
      # returning delta_time in observed frame 
      print('\n')
      print('Returning delta_time in observed frame, t_obs')
      c = 0
      for File in good_masterQ: #  len(masterFiles_Q)
          qso_data = add_tau_delflx(File,inDir_Q, qso_data)
          c += 1
          if c % 5 == 0:
            progress = (100.0*c) / float(len(good_masterQ))
            update_progress(progress)
            #print('\r[%-10s] %0.2f%%' % ('#' * int(progress/10), progress),)
              #update_progress(pers)
	            #print('\r----- Already read %d%% of qso'%pers),
    
    ### READ IN BLUE STARS ###
    
    print('\n')
    print('Reading in blue stars ...')
    c = 0                   
    for File in good_masterSB:    # [:len(good_masterQ)]
        star_data_blue = add_tau_delflx(File, inDir_S,star_data_blue)
        c += 1 
        if c % 5 == 0:
            progress = (100.0*c) / float(len(good_masterSB))
            update_progress(progress)

    ### READ IN RED STARS ###          
    
    print('\n')
    print('Reading in red stars ...')
    if good_ids_S_red is not None:
        c = 0                         
        for File in good_masterSR:   # [:len(good_masterQ)]
            star_data_red = add_tau_delflx(File, inDir_S, star_data_red)      
            c += 1               
            if c % 5 == 0:
                progress  = (100.0*c) / float(len(good_masterSR))
                update_progress(progress)         
                
    if good_ids_S_red is not None:
        return  qso_data, star_data_blue, star_data_red
    else: 
        return qso_data, star_data_blue

