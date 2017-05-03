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
- given the objectId in each sample, we collect the delta_mag (aka xi) and combined_error (aka ei) from 'master' files 
- the interesting dimension is  delta_time.  We bin along  delta_time,  plotting the raw data, and then per each linearly spaced bin we find the following statistics for delta_magnitude points in the bin  :  the standard deviation, the robust standard deviation, the structure function, finding the intrinstic variability (sigma) using the AstroML methodology http://www.astroml.org/book_figures/chapter5/fig_posterior_gaussgauss.html
- four panels with these quantities for quasars in 18.5-19 magnitude range become Fig.2 

5) Use E_Fig_3_histogram_panels.ipynb  
- we gather all delta_magnitude  points in a given magnitude bin into broader delta_time bins :   
	0 < log(delta_time) < 1.7 ,  
	2.3 < log(delta_time) < 2.5 , 
	2.8 < log(delta_time) < 3.0
- we use compare blue stars and quasars plotting histograms of a  quantity  chi = delta_magnitude / combined_error . This is Fig.3 in the paper

6) Use F_Fig_4_SF_three_panels.ipynb
- we use correction factors for combined_errors derived from the robust standard deviation of chi quantity for blue stars  (where this quantity should be a chi-like distribution centered on 1 width 1 if error model was correct)
- we apply the blue stars - derived correction to quasars and stars,  and show that in each magnitude range stars drop to 0 in standard deviation (as expected by construction). 

7) Use G_Fig_5_SF_panels_redshift_corr.ipynb
- we do the same as before in Fig 3-4, but applying to quasars in restframe
- this does not alter the results - there is still no significant signal for quasars for timescales < 100 days . 

8) Use 01.2017_PTF_data_for_CRTS_project.ipynb 
- we obtain light curves for the same objects based on their SDSS ra, dec , querying the PTF catalog https://irsa.ipac.caltech.edu/frontpage/ 

9) PTF_*  code 
- we follow the same steps as above, adapting the code to data for stars and quasars from PTF survey 

10) Investigating a special object : PG1302-102.ipynb
- we investigate an object found to have an unusually high variability from Vaughan+2016 http://adsabs.harvard.edu/cgi-bin/bib_query?2016MNRAS.461.3145V 
- we obtain the lightcurve from CRTS DR2 http://nesssi.cacr.caltech.edu/cgi-bin/getcssconedbid_release2.cgi 
- we combine epochs within a time window of a given span, and calculate statistics such as weighted magnitude, standard deviation weighted  by errors, chi2dof, robust standard deviation
- we repeat the same statistics for other quasars from CRTS, as well as for nearby stars  


11) Use CRTS_stars_robust_sigmaG.ipynb 
- for all CRTS stars, we perform 2.5 sigma clipping, counting the number of outliers, and for the remaining part of the light curve, calculate the median  magnitude , median error, robust standard deviation , weighted mean magnitude 
- for all CRTS stars, after sigma-clipping, based on a quantity z = magnitude - mean_magnitude / error  we calculate chi2dof, and robust interquartile-based standard deviation of z/ 
- we plot median(error), sigmaG(magnitude), and sigmaG(z)  as a function of median(magnitude) - this becomes the Appendix 

12) Use D_Fig_2_CRTS_sel_r_cut.ipynb
- we investigate the wiggles 
- this may become a new notebook - at the moment it is part of D_Fig2 because we compare the results with Fig2.  


