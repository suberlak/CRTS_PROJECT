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


##############################
# CONVERT RA DEC TO DEGREES  #
##############################

def get_ra_dec_CRTS(ra_hms, dec_hms):
    """
    Extracting RA, DEC information from the QSO  name for the CRTS case 
    """
    
    dec_hms_split = np.empty(0, dtype=str)
    for i in range(0,len(dec_hms)):
        dec_hms_split = np.append(dec_hms_split,  dec_hms[i][0:3]+' '+dec_hms[i][3:5]+' '+dec_hms[i][5:9])
    
    
    ra_hms_split = np.empty(0, dtype=str)
    for i in range(0,len(ra_hms)):
        ra_hms_split = np.append(ra_hms_split, ra_hms[i][0:2]+' '+ ra_hms[i][2:4] + ' ' + ra_hms[i][4:9])
    return ra_hms_split, dec_hms_split
     
def HMS2deg(ra='', dec=''):
    """
    From http://www.bdnyc.org/2012/10/15/decimal-deg-to-hms/  
    Converting  ra and dec from h:m:s   and deg:m:s  to  degrees.decimal 
    I assume they are using ICRS coordinates 
    """
    RA, DEC, rs, ds = '', '', 1, 1
    if dec:
        D, M, S = [float(i) for i in dec.split()]
        if str(D)[0] == '-':
            ds, D = -1, abs(D)
        deg = D + (M/60) + (S/3600)
        DEC = '{0}'.format(deg*ds)
      
    if ra:
        H, M, S = [float(i) for i in ra.split()]
        if str(H)[0] == '-':
            rs, H = -1, abs(H)
        deg = (H*15) + (M/4) + (S/240)
        RA = '{0}'.format(deg*rs)
      
    if ra and dec:
        return (RA, DEC)
    else:
        return RA or DEC
    
    
def convert_to_deg(ra_split, dec_split):
    '''
    Converts ra and dec from h:m:s  extracted from the quasar name with 
    get_ra_dec_CRTS()  to  degrees.decimal , using HMS2deg() function
    '''
    ra_deg  =   np.empty(0, dtype=float)
    dec_deg =   np.empty(0, dtype=float)
    
    for i in range(0,len(dec_split)):
        dec_deg = np.append(dec_deg,float(HMS2deg(dec=dec_split[i])))
    
    
    for i in range(0,len(ra_split)):
        ra_deg = np.append(ra_deg,float(HMS2deg(ra=ra_split[i])))
    return ra_deg , dec_deg
    

    
def get_qso_catalog():
    File = '../data_products/CRTS_SDSS_catalogs/CRTS_SDSS_cross_matched_qso_DB_QSO_catalog.txt'
    colnames = open(File,'r').read().splitlines()[0][1:].split()
    datatable = np.genfromtxt(File, dtype=str)
    qso_catalog = {}
    print('Zipping CRTS-SDSS quasars catalog from %s' % File)
    for label, column in zip(colnames, datatable.T):
        qso_catalog[label] = column
    
    print('Read in %d quasars from CRTS' % len(qso_catalog['redshift']))
    return  colnames, qso_catalog
    
def get_stars_catalog():
    File = '../data_products/CRTS_SDSS_catalogs/CRTS_SDSS_cross_matched_stars_catalog.txt'
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
            mErrMin = -9, mErrMax = 0.3,cut_mag='r', redshift = None, match_radius_arcsec = 2.0 , survey='CRTS'):
    ''' A short  routine to select CRTS quasars according to desired parameters.
    
    Parameters:
    ----------
    qso_cat  - a CRTS-SDSS matched catalog, the output of get_qso_catalog()
    mMin / mMax - minimum / maximum  desired magnitude , of SDSS filter specified 
              by cut_mag  
    mErrMin / mErrMax - minimum / maximum  desired average CRTS/PTF lightcurve error, AFTER 
              day-averaging 
    cut_mag - desired SDSS filter for the cut,  allowed values  are u,g,r,i,z. 
              Default setting is 'r' 
    match_radius_arcsec  - a maximum matching radius  allowed   in arcsec - the CRTS-SDSS matched 
              catalog stores the matching radius in degrees. By default, we 
              choose 2 arcsec , i.e. 2/3600 of a degree as an appropriate 
              matching radius for QSO.  For the CRTS sample, only 15 do not 
              have an SDSS counterpart within that range. 

    Returns:
    --------
    qso_id : the CRTS id's of the quasars accepted  through the selection criteria 
    '''
    print('\n Returning only QSO with  an SDSS counterpart within %f arcsec'%match_radius_arcsec)
    
    if survey  == 'CRTS' : 
        match_radius_degrees = match_radius_arcsec / 3600.0 
        mask_rad = (qso_cat['m_ang_deg'].astype(float) < match_radius_degrees )
        mask_mag = (qso_cat[cut_mag].astype(float) > mMin) * (qso_cat[cut_mag].astype(float) < mMax) 
        mask_err = (qso_cat['CRTS_avg_e'].astype(float) > mErrMin) * (qso_cat['CRTS_avg_e'].astype(float) < mErrMax)
        mask = mask_rad * mask_mag * mask_err 
        qso_id = qso_cat['CRTS_id'][mask]
    
    if survey =='PTF':
        mask_rad = (qso_cat['match_radius_arcsec'] < match_radius_arcsec )
        mask_mag = (qso_cat[cut_mag] > mMin) * (qso_cat[cut_mag] < mMax)
        mask_err = (qso_cat['avg_day_err']> mErrMin) * (qso_cat['avg_day_err'] < mErrMax)
        mask = mask_rad * mask_mag * mask_err 
        qso_id = qso_cat['ra_sdss'][mask]


    print('\n These cuts reduced the number of qso  in the sample from %d to %d' %
        (len(qso_cat['redshift']), len(qso_id)))
          

    if redshift is not None:
        print('\n  Also returning quasar redshifts...')
        qso_redshift = np.array(qso_cat['redshift'][mask]).astype(float)
        return qso_id, qso_redshift 
    else:
        return  qso_id

def cut_stars(star_cat=None, mMin=-9, mMax=19, mErrMin = -9, 
              mErrMax = 0.3, gi_Min = -1, gi_Max=1 , cut_mag='r_mMed', survey='CRTS'):
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
    print('\nChoosing stars with  SDSS   %.2f<g-i<%.2f'%(gi_Min, gi_Max))
    mask_mag = (star_cat[cut_mag] > mMin) * (star_cat[cut_mag] < mMax) 
    SDSS_gi = star_cat['g_mMed'] - star_cat['i_mMed']
    mask_color = (SDSS_gi > gi_Min ) * (SDSS_gi < gi_Max)
    
    if survey == 'CRTS' : 
        mask_err = (star_cat['CRTS_Merr'] > mErrMin) * (star_cat['CRTS_Merr'] < mErrMax)
        mask = mask_mag * mask_err * mask_color
        star_id_f = star_cat['crts_id'][mask]
 
        # convert floats to strings without comma and zeros
        star_id = np.array(["{:.0f}".format(name) for name in star_id_f])
        print(' These cuts reduced the number of stars  in the sample from %d to %d'%
            (len(star_cat['CRTS_M']),  len(star_id)))

    if survey == 'PTF' : 
        mask_err = (star_cat['avg_day_err'] > mErrMin) * (star_cat['avg_day_err'] < mErrMax)
        mask = mask_mag * mask_err * mask_color
        star_id = star_cat['ra_sdss'][mask]
        print(' These cuts reduced the number of stars  in the sample from %d to %d'%
            (len(star_cat['ra_sdss']),  len(star_id)))

    return  star_id


def faster_read_xi_ei(inDirSF = None, good_ids = None , detailed = None ):
    ''' A trimmed down version of fast_read_xi_ei : there is no  need really
    to have all variables as input : it is conceptually easier to 
    call a simple function three times, than to make a function three 
    times as complicated.  Plus,  that way we don't necessarily have 
    to be reading qso, starsR , or starsR : it could be stars, 
    or just qso, or only starsR...   Using the same engine as fast_read_xi_ei(), 
    so speed of execution should be the same . 
   
    Parameters: 
    -----------
    inDirSF : a directory where we can find the 'master' files.  Whether it should point to the 'simple'
          or 'detailed' master files, is up to user. Examples:   
          inDirSF = '../data_products/sf_file_per_LC/stars/'
          inDirSF = '../data_products/sf_file_per_LC/qso/' 
          inDirSF = '../data_products/sf_file_per_LC/stars_detailed/'
          inDirSF = '../data_products/sf_file_per_LC/qso_detailed/'
          inDirSF = '../data_products/sf_file_per_LC_PTF/stars/'
          inDirSF = '../data_products/sf_file_per_LC_PTF/qso/' 
    good_ids : a list of object names for which we should read in the master files. 
        We assume that each   master file is called    'SF_' + name + '.txt'
    detailed =  if not None,  then apart from three standard columns in master files 
         (delta_time,  delta_mag,  error(delta_mag)),  we also expect   t1, t2,  m1, m2,  e1,  e2, 
          i.e. all the data taht was used to prepare master files...     

    Returns :
    ----------
    store : a dictionary with keys corresponding to nonzero either  xi, tau, ei, or  that plus 
          t1, t2,  m1, m2,  e1,  e2,  eg. 

          store = {'xi':xi_flat, 'tau':tau_flat, 'ei':ei_flat}
          or 
          store = {'xi':xi_flat, 'tau':tau_flat, 'ei':ei_flat, 't1':t1_flat, 
                   't2':t2_flat, 'm1': m1_flat, 'm2':m2_flat, 'e1':e1_flat, 'e2': e2_flat}

    '''
    # rename some variables to slim the code 
    ids, inDir =  good_ids, inDirSF
    print('\nReading in tau,xi,ei  for %d objects'%len(good_ids))
    print('\nUsing structure function master files from %s'%inDir)
    

    # initialize storage arrays... 
    tau_store = list(np.zeros(len(ids)))
    xi_store = list(np.zeros(len(ids)))
    ei_store = list(np.zeros(len(ids)))
    
    if detailed is not None : 
        print('\nReading in detailed SF : t1, t2, e1, e2, m1, m2')
        t1_store = list(np.zeros(len(ids)))
        t2_store = list(np.zeros(len(ids)))
        e1_store = list(np.zeros(len(ids)))
        e2_store = list(np.zeros(len(ids)))
        m1_store = list(np.zeros(len(ids)))
        m2_store = list(np.zeros(len(ids)))

    # loop over all master files .. 
    c = 0

    for i in range(len(ids)) : 
        File = 'SF_' + ids[i] +'.txt'
        address = inDir+File
        data =  np.loadtxt(address)

        xi_store[i]  = data[:,0]
        tau_store[i] = data[:,1]
        ei_store[i]  = data[:,2]

        if detailed is not None : 
            t1_store[i] = data[:,3]
            t2_store[i] = data[:,4]
            m1_store[i] = data[:,5]
            m2_store[i] = data[:,6]
            e1_store[i] = data[:,7]
            e2_store[i] = data[:,8]
      
        c += 1 
        if c % 5 == 0:
            progress = (100.0*c) / float(len(ids))
            update_progress(progress)

    # flatten the arrays 
    tau_flat = np.concatenate(tau_store)
    xi_flat  = np.concatenate(xi_store)
    ei_flat  = np.concatenate(ei_store)
    if detailed is not None : 
        t1_flat = np.concatenate(t1_store)
        t2_flat = np.concatenate(t2_store)
        m1_flat = np.concatenate(m1_store)
        m2_flat = np.concatenate(m2_store)
        e1_flat = np.concatenate(e1_store)
        e2_flat = np.concatenate(e2_store)

    # assign them to the output dic.  
    if detailed is None : 
        store = {'xi':xi_flat, 'tau':tau_flat, 'ei':ei_flat}
    elif detailed is not None : 
        store = {'xi':xi_flat, 'tau':tau_flat, 'ei':ei_flat, 't1':t1_flat, 
        't2':t2_flat, 'm1': m1_flat, 'm2':m2_flat, 'e1':e1_flat, 'e2': e2_flat}

    print('\nFinished reading all master files for the selected objects ...')
    return store 

def fast_read_xi_ei(inDirStars = '../data_products/sf_file_per_LC/stars/', 
                    inDirQSO= '../data_products/sf_file_per_LC/qso/',  
                    good_ids_S_blue = None , good_ids_S_red = None , good_ids_QSO = None, detailed = None ):
    ''' A function combining what read_xi_ei and add_tau_delflx  would do much slower. 
   
    Parameters: 
    -----------
    inDirStars : by default, 
    inDirQSO
    good_ids_S_blue
    good_ids_S_red
    good_ids_QSO
    
    Returns :
    ----------
    store : a dictionary with keys corresponding to nonzero object lists provided as input. 
    '''
    store = {}   
    
    good_ids_list = [good_ids_S_blue, good_ids_S_red, good_ids_QSO]
    obj_type = ['starsB', 'starsR', 'qso']
    inDir_list = [inDirStars,inDirStars,inDirQSO]


    for i in range(3) : 
        good_ids = good_ids_list[i]
        if good_ids is not None : 
            obj = obj_type[i]
            print('\nReading in tau,xi,ei for  %s'%obj)

            ids = good_ids
            inDir = inDir_list[i]

            # initialize storage arrays... 
            tau_store = list(np.zeros(len(ids)))
            xi_store = list(np.zeros(len(ids)))
            ei_store = list(np.zeros(len(ids)))
            
            if detailed is not None : 
                t1_store = list(np.zeros(len(ids)))
                t2_store = list(np.zeros(len(ids)))
                e1_store = list(np.zeros(len(ids)))
                e2_store = list(np.zeros(len(ids)))
                m1_store = list(np.zeros(len(ids)))
                m2_store = list(np.zeros(len(ids)))
            # loop over all master files .. 
            c = 0

            for i in range(len(ids)) : 
                File = 'SF_' + ids[i] +'.txt'
                address = inDir+File
                data =  np.loadtxt(address)

                xi_store[i]  = data[:,0]
                tau_store[i] = data[:,1]
                ei_store[i]  = data[:,2]

                if detailed is not None : 
                    t1_store[i] = data[:,3]
                    t2_store[i] = data[:,4]
                    m1_store[i] = data[:,5]
                    m2_store[i] = data[:,6]
                    e1_store[i] = data[:,7]
                    e2_store[i] = data[:,8]
              
                c += 1 
                if c % 5 == 0:
                    progress = (100.0*c) / float(len(ids))
                    update_progress(progress)

            # flatten the arrays 
            tau_flat = np.concatenate(tau_store)
            xi_flat  = np.concatenate(xi_store)
            ei_flat  = np.concatenate(ei_store)
            if detailed is not None : 
                t1_flat = np.concatenate(t1_store)
                t2_flat = np.concatenate(t2_store)
                m1_flat = np.concatenate(m1_store)
                m2_flat = np.concatenate(m2_store)
                e1_flat = np.concatenate(e1_store)
                e2_flat = np.concatenate(e2_store)

            # assign them to the output dic. 
            if detailed is None : 
                store[obj] = [xi_flat, tau_flat, ei_flat]
            elif detailed is not None : 
                store[obj] = [xi_flat, tau_flat, ei_flat, t1_flat, t2_flat, 
                              m1_flat, m2_flat, e1_flat, e2_flat]


    print('\nFinished reading all master files for the selected objects ...')
    return store 





# inside the main loop : get tau, delflx from a master file, either qso or star
def add_tau_delflx(File, inDir, data, z=None, tau_min = None, tau_max = None):
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
    z : if not None, then we correct for redshift, using provided z. It should 
         be a scalar corresponding to the object for which we are reading in 
         xi, ei, tau points 

    tau_min : if not None, a scalar corresponding to the min tau that we want to 
         keep
    tau_max : if not None, a scalar corresponding to the max tau that we want to 
         keep
    
    Returns:
    ---------
    delflx, tau, err, master_names : updated elements of data array, 
         to which we appended  the values of delflx, tau, err, and a 
         name of the object (master file)  multiplied by the number of 
         rows in the master file 
 
    '''
    # read in storage arrays
    store_xi  = data[0]  
    store_tau = data[1]
    store_ei  = data[2]
    store_object_list = data[3]   
    
    # grab the object name 
    object_name = File[3:-4]
    
    # read in the i-th master file 
    object_data =  np.genfromtxt(inDir+File, dtype=str)
    file_xi  = object_data[:,0].astype(float)
    file_tau = object_data[:,1].astype(float)
    file_ei  = object_data[:,2].astype(float)

    # introduce a mask, by default it accepts everything, 
    # since time difference is definitely bigger than 0 ... 
    m1 = 0 < file_tau
    m2 = 0 < file_tau 

    # take only xi,ei,tau corresponding to the desired range ...
    if tau_min is not None: 
        m1 = tau_min < file_tau 

    if tau_max is not None : 
        m2 = file_tau < tau_max 

    m = m1 * m2 # combine the two masks 

    # read in tau,  del_mag,  del_mag_err for quasars on the list 
    store_xi = np.append(store_xi, file_xi[m])
    if z is not None:
        store_tau = np.append(store_tau,  file_tau[m] / (1.0+z))
    else:
        store_tau = np.append(store_tau, file_tau[m])

    store_ei = np.append(store_ei, file_ei[m])
    store_object_list = np.append(store_object_list, np.array(len(file_xi[m])*[object_name]))
    
    return store_xi, store_tau, store_ei, store_object_list
    
def read_xi_ei(inDirStars, good_ids_S_blue, inDirQSO,
                 good_ids_QSO, good_ids_S_red=None, redshift=None,tau_min = None, 
                 tau_max = None):
    ''' A routine to read the delta_mag (xi), delta_time (tau), and error (ei) 
    for CRTS / PTF stars and quasars. Stars and quasar master files are read from 
    inDirStars and inDirQSO , and only those files are selected to be read in 
    that are on the list of good_ids_S_blue, good_ids_S_red  for blue 
    and red stars, and good_ids_QSO  for quasars.  

    Parameters 
    -----------
    inDirStars : dir with stellar master files 
    good_ids_S_blue : ids of stars within the cut
    inDirQSO : dir with QSO master files 
    good_ids_QSO : ids of QSO within a cut 
 
    Optional Parameters :
    ----------------------
    good_ids_S_red=None :  if specified, only red stars would be read 
    redshift=None : if not None,  then correcting QSOs for redshift ( expect an array of z )
    tau_min , tau_max : if not None,  scalars corresponding to the limits on 
         tau that we want to read-in from each master file 
    
    Returns:
    ---------
    qso_data, star_data_blue, star_data_red : 1D arrays with [ xi, tau, ei,  name ]
        per object type ... 

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
    masterFiles_Q = os.listdir(inDir_Q)
    masterFilesQ1 = [name[3:-4] for name in masterFiles_Q]
    good_masterQ =  np.array(masterFiles_Q)[np.in1d(masterFilesQ1, good_ids_QSO)]
    #good_masterQ = np.array(['SF_' +qso+'.txt' for qso in good_ids_QSO])

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
    print('\n')
    print('Reading in quasars...')
    if redshift is not None: 
      # correcting QSO delta_time to restframe
      print('Correcting delta_time to restframe, t_rest = t_obs / (1+z)')
      c = 0
      for i in range(len(good_masterQ)): #  len(masterFiles_Q)
          z = redshift[i]
          File = good_masterQ[i]
          if tau_min is not None : 
              qso_data = add_tau_delflx(File,inDir_Q, qso_data, z, tau_min=tau_min, tau_max=tau_max)
          else:
              qso_data = add_tau_delflx(File,inDir_Q, qso_data, z)
          c += 1 
          if c % 5 == 0:
            progress = (100.0*c) / float(len(good_masterQ))
            update_progress(progress)
            #print('\r[%-10s] %0.2f%%' % ('#' * int(progress/10), progress),)
              #print('\r----- Already read %d%% of qso'%pers),
    else:
      # returning delta_time in observed frame 
      print('Returning delta_time in observed frame, t_obs')
      c = 0
      for File in good_masterQ: #  len(masterFiles_Q)
          if tau_min is not None :
            qso_data = add_tau_delflx(File,inDir_Q, qso_data, tau_min=tau_min, tau_max=tau_max)
          else: 
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
        if tau_min is not None : 
          star_data_blue = add_tau_delflx(File, inDir_S,star_data_blue, tau_min=tau_min, tau_max=tau_max)
        else: 
          star_data_blue = add_tau_delflx(File, inDir_S,star_data_blue)
        c += 1 
        if c % 5 == 0:
            progress = (100.0*c) / float(len(good_masterSB))
            update_progress(progress)

    ### READ IN RED STARS ###          
    
    print('\n')
    if good_ids_S_red is not None:
        print('Reading in red stars ...')
        c = 0                         
        for File in good_masterSR:   # [:len(good_masterQ)]
            if tau_min is not None : 
                star_data_red = add_tau_delflx(File, inDir_S, star_data_red, tau_min=tau_min, tau_max=tau_max)    
            else: 
                star_data_red = add_tau_delflx(File, inDir_S, star_data_red) 
            c += 1               
            if c % 5 == 0:
                progress  = (100.0*c) / float(len(good_masterSR))
                update_progress(progress)         

    if good_ids_S_red is not None:
        return  qso_data, star_data_blue, star_data_red
    else: 
        return qso_data, star_data_blue

