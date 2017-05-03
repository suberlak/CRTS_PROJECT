# CRTS_PROJECT

Files needed to reproduce the analysis of CRTS errors. We compute Structure Functions for CRTS data: quasars and stars. We use quasars from the SDSS S82 quasar catalog, for which we obtain lightcurves from the CRTS DR2 database. We use stars  also  from the SDSS S82 standard stars catalog : lightcurves for 10 % randomly chosen  stars  are obtained from the CRTS DR2  (the full SDSS S82 catalog has over one  million objects).  

All dependencies (custom-written modules) are in CRTS_paper_modules.py

All SDSS information about the S82 standard stars is from http://www.astro.washington.edu/users/ivezic/sdss/catalogs/stripe82.html  

All SDSS information the S82 quasars is from DB_QSO_S82.dat from http://www.astro.washington.edu/users/ivezic/cmacleod/qso_dr7/Southern.html

All time-series data (light curves) are from CRTS, available locally at the UW (courtesy of B. Sesar). Another way to query for CRTS  QSO and standard stars is to use multi-object query http://nesssi.cacr.caltech.edu/cgi-bin/getmulticonedb_release2.cgi   , qiven the list of SDSS ra,dec . However, this query only supports up to 100 objects at a time. CRTS is not currently (May 2017) accessible via IRSA 	http://irsa.ipac.caltech.edu/applications/Gator/.  If you know of any direct way to efficiently query CRTS for multiple light curves, please contact  suberlak@uw.edu, and we will gladly update this page. 

Steps taken : 

1) Use A_Fig_1_QSO_CRTS_day_averaged_stats.ipynb and A_Fig_1_stars_CRTS_stats.ipynb
- for CRTS QSO and stellar light curves, we select only those that have more than 10  epochs. 
- for illustration and to describe the dataset, per raw lightcurve we calculate the light curve time span , number of rows , mean magnitude, mean error ('raw_timespan_obs', 'raw_lc_length', 'raw_mean_mag', 'raw_mean_err'). This information is stored as Fig_1_raw_CRTS_**_LC_stats.npy  file, where ** stands for 'QSO' or 'stars' 
- to improve statistics, we day-average the raw lightcurves. The weighted mean is calculated using inverse errors as weights, as described in the paper
- for the day-averaged lightcurves (called 'processed' in the code), we calculate  light curve time span, number of rows , number of observations per night that were averaged into one epoch, mean magnitude, mean error, median pairwise time difference ('proc_mjd_span', 'proc_lc_length', 'proc_mean_N_day', 'proc_mean_mag', 'proc_mean_err', 'proc_median_dt'). This information is stored as Fig_1_proc_CRTS_**_LC_stats.npy  file, where ** stands for 'QSO' or 'stars'
- the before- and after-  day-averaging aggregate statistics catalogs are combined into one to ease comparison : CRTS_**_raw_and_proc_combined_agg.dat file, where ** stands for 'QSO' or 'stars'
- using this catalog we make Fig 1 with statistics for day-averaged CRTS S82 QSO 

2) Use B_CRTS_SDSS_matching.ipynb  
- match the CRTS quasars and stars by ra,dec position to the SDSS counterparts (they were obtained based on SDSS ra,dec so there is a perfect match, but we were only provided lightcurves with ra,dec of each object, so this matching allows to do 1-to-1 mapping of CRTS to SDSS objects) 
- we use AstroPy SkyCoord to do 2D catalog matching in the most efficient manner 
- save the CRTS-SDSS catalogs to CRTS_SDSS_combined_**_catalog.dat file, where ** stands for 'QSO' or 'stars'
- (note :  B_sf_CRTS_SDSS_matching_NEW.py  is now deprecated )

3) Use C_sf_load_single_LC.py  
- this program makes 'master' files per object light curve,  containing   delta_time,  delta_magnitude,  combined_error 
- we also investigated in detail the dependence between cadence and structure of the pairwise covariance of magnitude measurements. To do that we extracted delta_time, delta_magnitude, combined_error,  t1, t1, m1, m2,  e1, e2  , using C_sf_load_single_LC_detail.ipynb 


4) Use D_Fig_2_CRTS_sel_r_cut.ipynb
- based on SDSS r-magnitude we select objects in a given class (stars or QSO), within the magnitude bin (17-18, 18-18.5, 18.5-19). Stars are split into 'red' and 'blue' samples  by SDSS g-i  color
- given the objectId in each sample, we read the delta_mag (aka xi) and combined_error (aka ei) from 'master' files 


