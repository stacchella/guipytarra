c
c     variables are defined as double precision
c     arrays used for images are either real*4 or integer
c
      implicit none
      double precision mirror_area, integration_time, gain,
     *     time_since_previous, decay_rate,  read_noise,
     *     dark_mean, dark_sigma, ktc, voltage_offset
      double precision dark_mean_cv3, dark_sigma_cv3, gain_cv3,
     *     read_noise_cv3

      double precision effective_wl_nircam, width_nircam, 
     *     system_transmission, bandwidth, wavelength,
     *     photplam, photflam, stmag, abmag
      double precision ra_dithered, dec_dithered
      double precision scale
      double precision linearity_gain,  lincut, well_fact
      double precision zodiacal_scale_factor, bkg, background
      double precision mstar, alpha, phistar, p_evol, q_evol
      double precision zdist, expected, int_zdist, mag, counts,n_central
      double precision bright, faint, apm_bright, apm_faint, zmin, zmax,
     *     solid_angle
      double precision  h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
      double precision rcore_arcsec
      double precision pi, q
      double precision xc, yc, new_ra, new_dec, cosdec, osim_scale
      double precision dx, dy, pa_v3, ra0, dec0, x0, y0, angle, 
     *     dither_arc_sec, dtheta, dither_step, pa_degrees
      double precision xshift, yshift, xmag, ymag, xrot,
     *     yrot, xrms, yrms, oxshift, oyshift, oxmag, oymag,
     *     oxrot, oyrot, oxrms, oyrms
      double precision mag_array, int_lf, lumfun,m_tot, m_disk,
     *     m_bulge, mmax, lf_scale, inc_d
      double precision ra_stars, dec_stars, mag_stars, stellar_photons
      double precision ra_galaxies, dec_galaxies, magnitude,
     *     z, nsersic, ellipticity, re, theta, flux_ratio, id
      double precision total_time, tot
      double precision two_n, bn, n_profile, r_eff
c
      double precision integrated_psf
      double precision filters, filtpars
c
      double precision xdither, ydither
c
      real cr_matrix, cr_flux, cr_accum
c
c      real  
      double precision 
     *     equinox, crpix1, crpix2, crval1, crval2, cdelt1,cdelt2,
     *     cd1_1, cd1_2, cd2_1, cd2_2
      real    image, tframe, tgroup, accum,latent_image, gain_image, 
     *     base_image, dark_image, well_depth, bias, linearity
c
      integer  order
      integer   nxy, nxny
      integer   n_psf_x, n_psf_y, n_image_x, n_image_y
      integer   nnn, nz, max_stars, max_objects, npts, nstars,
     *     ngal, max_groups, max_order
      integer   over_sampling_rate, time_step, verbose, bkg_mode
      integer   frame,  bitpix, sca_id, nf, cr_mode,nf_used
      integer   ngroups, nframe, nskip, idither, ndither
      integer include_ktc, include_bg, include_cr, include_dark,
     *     include_latents, include_readnoise, include_non_linear,
     *     include_1_over_f, include_reference,
     *     include_stars, include_star_cluster, include_galaxies,
     *     include_cloned_galaxies, brain_dead_test
      integer cube, overlap
      integer option, use_filter, cat_filter, filters_in_cat
      integer   seed, i, i1, i2, indx, j, icat_f
      integer n_cr_levels, ncr
      integer number_primary, number_subpixel
      integer parallel, do_dither, first_dither, last_dither
c
      logical noiseless, psf_add, ipc_add
c
c     new (ISIM CV2) keywords for subarray
c
      integer nfilter_wl, nfilters, npar, nbands, nsub, ncomponents

      character filterid*20, temp*20, galaxy_catalogue*180, 
     *     star_catalogue*180, filename_param2*180, filter_id*5


      integer   colcornr, rowcornr, naxis1, naxis2
      character subarray*8
c
      character psf_file*180, read_patt*10, filename*180, text*10,
     *     clone_path*180, noise_name*180
      character object*20, partname*20, module*20
      character primary_dither*20, subpixel_dither*20, camera*20
c
      parameter (nfilter_wl = 25000, nbands=200, npar=30)
      parameter (nnn=2048,nz=160, ncr=21, nfilters=54)
      parameter (overlap = 50)
      parameter (max_stars=10000, max_objects=50000, nsub=4)
      parameter (max_order = 7)
      parameter (nxny = 2048*2048)
c
      dimension zdist(nnn), expected(nnn), int_zdist(nnn)
      dimension dx(1000), dy(1000),pa_v3(1000)
      dimension mag_array(nnn), lumfun(nnn), int_lf(nnn)
c
c     SCA parameters
c
      dimension dark_mean_cv3(10), dark_sigma_cv3(10), gain_cv3(10),
     *     read_noise_cv3(10)
      dimension dark_mean(10), dark_sigma(10), gain(10),
     *     read_noise(10)
c     
c     images
c
      dimension cube(nnn, nnn, overlap)
      dimension base_image(nnn,nnn)
      dimension accum(nnn,nnn), image(nnn,nnn), latent_image(nnn,nnn)
      dimension gain_image(nnn,nnn), dark_image(nnn,nnn,2),
     *     well_depth(nnn,nnn), linearity(nnn,nnn,max_order),
     &     bias(nnn,nnn)
c
c     filter parameters, average responses from filters
c
      dimension effective_wl_nircam(nfilters), width_nircam(nfilters), 
     *     system_transmission(nfilters)
      dimension bkg(nfilters), 
     &     cat_filter(nfilters),use_filter(nfilters)

      dimension filters(nbands, nfilter_wl), filtpars(nbands,npar), 
     *     filterid(nbands)
      dimension psf_file(nfilters), integrated_psf(nxny)
c
c     dithers
c
      dimension xdither(nnn), ydither(nnn)
c
c     catalogues
c
      dimension cr_matrix(ncr,ncr,10000),cr_flux(10000),cr_accum(10000)
      dimension ra_stars(max_stars), dec_stars(max_stars),
     *     mag_stars(max_stars ,nfilters), counts(max_stars)
c
      dimension xshift(10), yshift(10), xmag(10), ymag(10), xrot(10),
     *     yrot(10)
      dimension oxshift(10), oyshift(10), oxmag(10), oymag(10),
     *      oxrot(10), oyrot(10)
c
      dimension ra_galaxies(max_objects), dec_galaxies(max_objects), 
     *     z(max_objects),
     *     magnitude(max_objects,nfilters), ncomponents(max_objects),
     *     id(max_objects),
     *     nsersic(max_objects, nsub), ellipticity(max_objects, nsub),
     *     re(max_objects, nsub), theta(max_objects, nsub), 
     *     flux_ratio(max_objects, nsub),
     *     clone_path(max_objects)
c
c     images
c     
      common /history/ cube
      common /base/   base_image
      common /dark_/  dark_image
      common /gain_/  gain_image
      common /well_d/ well_depth, bias, linearity,
     *     linearity_gain,  lincut, well_fact, order
      common /latent/ latent_image
      common /images/ accum, image, n_image_x, n_image_y
      common /psf/ integrated_psf, n_psf_x, n_psf_y, nxy, 
     *     over_sampling_rate
c
      common /cr_list/ cr_matrix, cr_flux, cr_accum, n_cr_levels
c
c     catalogues
c
      common /stars/ ra_stars, dec_stars, mag_stars, nstars
      common /galaxy/ra_galaxies, dec_galaxies, z, magnitude, 
     *     nsersic, ellipticity, re, theta, flux_ratio, ncomponents
c
      common /transform/ xshift, yshift, xmag, ymag, xrot, yrot
      common /otransform/ oxshift, oyshift, oxmag, oymag, oxrot, oyrot
c
      common /parameters/ mirror_area, integration_time, gain,
     *     decay_rate, time_since_previous, read_noise, 
     *     dark_mean, dark_sigma, ktc, voltage_offset 
      common /throughput/ effective_wl_nircam, width_nircam,
     *     system_transmission
c
      common /filter/filters, filtpars, filterid
c
      common /gdif_par/ two_n
      common /sersic_par/ bn, n_profile, r_eff
c      common /cosmology/ h0, omega_m, omega_l, omega_k, omega_r, 
c     *     omega_nu, omega_cm, w, wprime
c      common /schechter_params/ mstar, alpha, phistar, p_evol, q_evol
c      common /sample_params/   bright, faint, apm_bright, apm_faint, 
c     *     zmin, zmax,solid_angle
c
      common /wcs/ equinox, crpix1, crpix2, crval1, crval2,
     *     cdelt1,cdelt2, cd1_1, cd1_2, cd2_1, cd2_2
c
c     centre coordinates for combined SCAs (OSIM coordinates)
c
c      data xc, yc /0.0d0, -315.6d0/
c      data osim_scale/1.59879d0/ 
      data xc, yc/0.0d0, 0.0d0/
      data osim_scale/60.d0/ 
c
c     ISIM CV3 values in e-/sec from K. Misselt
c      
      data dark_mean_cv3 /0.001d0, 0.003d0, 0.003d0, 0.003d0, 0.033d0,
     *     0.002d0, 0.001d0, 0.001d0, 0.001d0, 0.040d0/
      
      data dark_sigma_cv3/0.002d0, 0.005d0, 0.003d0, 0.002d0, 0.006d0,
     *     0.002d0, 0.002d0, 0.003d0, 0.004d0, 0.004d0/
c     e-
      data read_noise_cv3/11.3d0, 10.5d0, 10.2d0, 10.3d0, 8.9d0,
     *     11.5d0, 12.7d0, 11.3d0, 12.0d0, 10.5d0/
c
c     e-/ADU
c     gain_cv3 was measured using IDL's biweight_mean on the ISIMCV3
c     images, removing NaN, infinities and reference pixels
c     2016-07-05
c      data gain_cv3/2.077d0, 2.020d0, 2.166d0, 2.01d90, 1.845d0,
c     *     2.005d0, 2.440d0, 1.9381d0, 2.248d0, 1.7977d0/
      data gain_cv3/10*1.0d0/
c         
c     Cosmic rays
c     cr_mode:
c     0         -  Use ACS observed counts for rate and energy levels
c     1         -  Use M. Robberto models for quiet Sun
c     2         -  Use M. Robberto models for active Sun
c     3         -  Use M. Robberto models for solar flare
c
c     Zodiacal background
c
c     If using the SED this value seems to match 
c     the ST ETC low background
c      zodiacal_scale_factor = 0.858801d0
c     nominal value from MJR
c     zodiacal_scale_factor = 1.20d0
c     this gives a better match with MJR's values
c
c      zodiacal_scale_factor = 2.00302d0
c
c     bkg_mode ( 1= MJR;  2 = STSci ETC ; 3 = Calculate from SED, 
c     4 = read from input)
c     
c     constants
c
      pi = dacos(-1.0d0)
      q  = pi/180.d0
      nframe = 1
      noiseless = .false.
      psf_add   = .true.
      ipc_add   = .true.

c=======================================================================
c
c     Read dither position, the associated catalogue, detector
c     and filter for this scene. The catalogue must have been
c     already tailored to the detector area
c
      print *,' this should be run in batch mode :'
      print *,' guitarra < params_F356W_490_001.input'
      read(5,92) filename_param2
      print *,'other parameter file = ', filename_param2
      read(5,*,err=89) idither
 89   print *,'idither = ', idither
      read(5,*,err=90) filename
 90   print *, filename
      read(5,*,err=91) noise_name
 91   print *,noise_name
      read(5,*) ra0
      read(5,*) dec0
      read(5,*) new_ra
      read(5,*) new_dec
      read(5,*) dx(idither)
      read(5,*) dy(idither)
      read(5,*) sca_id
      read(5,*) indx
      read(5,*) icat_f
      read(5,92) star_catalogue
 92   format(a180)
      read(5,92) galaxy_catalogue
      read(5,*) filters_in_cat

    
c
c     Read parameters from parameter file 2
c
      call read_parameters(filename_param2, verbose, brain_dead_test,
     *     module, read_patt, ngroups, nframe,
     *     subarray, colcornr, rowcornr, naxis1, naxis2, 
     *     camera, primary_dither, subpixel_dither, 
     *     number_primary, number_subpixel,
     *     dither_arc_sec, ra0, dec0, pa_degrees, 
     *     include_ktc, include_dark, include_readnoise,
     *     include_reference,
     *     include_1_over_f, include_latents, include_non_linear,
     *     include_cr, cr_mode, include_bg, bkg_mode,
     *     zodiacal_scale_factor, 
     *     include_stars, nstars, star_catalogue,
     *     include_galaxies, ngal, include_cloned_galaxies, 
     *     galaxy_catalogue, nf, use_filter, cat_filter)
c
c=======================================================================
c
c     Set  cosmological parameters
c
c      call is_cosmology_defined
c
c     Instrument-related parameters
c
c     
c     NIRCam constants
c     
      object   ='Simulation'
      if(module .eq. '481') then
         sca_id = 481
      end if

      if(module .eq. '482') then
         sca_id = 482
      end if

      if(module .eq. '483') then
         sca_id = 483
      end if

      if(module .eq. '484') then
         sca_id = 484
      end if

      if(module .eq. '485') then
         sca_id = 485
      end if

      if(module .eq. '486') then
         sca_id = 486
      end if

      if(module .eq. '487') then
         sca_id = 487
      end if

      if(module .eq. '488') then
         sca_id = 488
      end if

      if(module .eq. '489') then
         sca_id = 489
      end if

      if(module .eq. '490') then
         sca_id = 490
      end if
c
c      if(sca_id .eq. 1) then
c         module   = 'modA'
c      end if
c      if(sca_id .eq. 2) then
c         module   = 'modB'
c      end if
      n_image_x = naxis1
      n_image_y = naxis2
c
c     SCA
c
      partname = 'none'
      partname = '16989'
c
      if(brain_dead_test.eq.1) then
         sca_id             = 0
         module             = 'none'
      end if
c     
c     mirror area from MJR (2014-May meeting)
c
      mirror_area = 25.37d0 * 1.0D4 
c     
c=======================================================================
c
      call set_params(read_patt, nframe, nskip, max_groups)
      max_groups = 160
      if(ngroups.gt.max_groups) ngroups  = max_groups
      if(max_groups.gt.nz) then
         print *,' maximum number of groups is ', nz
         stop
      end if
c     
c     integration time per read (is a function of the array size)
c     
      tframe   = (dble(naxis1)/2048.d0)*(dble(naxis2)/2048.d0)
c
c     readout time for CV3 is 10.73676 seconds
c
      tframe   = real(10.73676* tframe) ! single precision 
      integration_time = tframe
      tgroup   = tframe * (nframe+nskip)
c     
c     find total integration time given read mode, ngroups etc.
c
      ndither = number_primary * number_subpixel
      tot = total_time(nframe, nskip, ngroups, ndither,
     *     integration_time)
c
c=======================================================================
c
c     Instrument related parameters
c
c     Voltage offset (in e-) to add as a baseline to all pixels
c
      voltage_offset  = 500.d0
      voltage_offset  = 0.d0
c

c     kTC noise in e- (measured by K. Misselt)
c
      ktc = 30.d0               ! e-
      do i = 1, 10
         dark_mean(i)  = dark_mean_cv3(i)
         dark_sigma(i) = dark_sigma_cv3(i)
         gain(i)       = gain_cv3(i)
         read_noise(i) = read_noise_cv3(i)
      end do
c     
c     Latents:
c     decay rate comes Marcia Rieke's report for the EIDP OSIM SCAs
c     page 5-4: "latent images drop to < 1 % of the saturating signal
c     after the reset cycle finishes (e.g., by 31.8sec)"
c     If 1% is used: decay_time = 6.9053 s or
c     
      decay_rate      = 0.144817d0
c
c     Chad quotes 0.1% after 60 seconds. This corresponds to
c      decay_rate      = 0.11513d0
c     The value quoted in the report is
c      decay_rate      = 1.d0/55.12d0
c
c     WFC3 has  a latent decay rate of dark = 2.212 * 0.931**time + 0.446
c     this is the same as 2.212 * exp(time * ln(0.931)) + 0.446 
c     or                  2.212 * exp(-0.07150 * time) + 0.446
c                         2.212 * exp(-time/13.9868) + 0.446
c      decay_rate      = dlog(0.931d0)
c
c     reset time (time between end of exposure and start of next)
c     is about 31.8 s (M. Rieke EIDP_OSIM.pdf)
c
      time_since_previous = 31.8d0
c
c=======================================================================
c
c     read filter parameters
c
      call read_filter_parameters(nf_used, verbose)
c     
c     read list of fits filenames of point-spread-function
c
      call read_psf_list(psf_file)
c     
c     read zodiacal background
c
      call zodi(zodiacal_scale_factor, mirror_area, bkg, bkg_mode,
     &     nf_used, verbose)
c
c     read parameters that will be used to generate the WCS keywords
c
      call load_osim_transforms(verbose)
c
c=======================================================================
c
c     Initialise the random number generator
c     set seed to zero if the system clock will be used to seed
c     the random numbers. For repeatability set to 1
c
      seed = 0
      call zbqlini(seed)
c     
c**************************************************************************
c
c     Read source catalogues 
c
c**************************************************************************
c
      if(include_stars.gt.0) then 
         call read_star_catalogue(star_catalogue, nf_used, verbose)
         if(verbose.ge.2) print *,'read star catalogue'
      end if
c
      if(include_galaxies .eq. 1 .and. ngal .gt. 0) then 
         call read_fake_mag_cat(galaxy_catalogue, cat_filter, 
     &        filters_in_cat, nf, ngal)
      end if
      if(verbose.ge.2) then
         do i = 1, 10
            print * ,ra_galaxies(i), dec_galaxies(i), nf,
     &           (magnitude(i,j),j = 1,nf)           
         end do
      end if
c
c**************************************************************************
c
c     set more parameters
c
c**************************************************************************
c
      scale               =   0.0317d0
      if(sca_id .eq. 485 .or. sca_id .eq.490) then
         scale = 0.0648d0
      end if
 110  print *, idither, ra0, dec0, new_ra, new_dec, dx(idither),
     *     dy(idither), sca_id, indx, icat_f
      j                   = indx
      
      temp                = filterid(j)
      filter_id           = temp(1:5)
      wavelength          = filtpars(j,5) ! effective_wl_nircam(j)
      bandwidth           = filtpars(j,16) ! width_nircam(j)
      system_transmission = filtpars(j,8) 
      photplam            = filtpars(j,18)
      photflam            = filtpars(j,26)
      stmag               = filtpars(j,24)
      abmag               = filtpars(j,25)
      background          = bkg(j)
c     
c     use appropriate column of object catalogue
c
c      icat_f              = cat_filter(indx)
c     
c     use appropriate filters for SW/LW
c     check for inconsistencies
c
      if(scale.eq.0.0317d0.and.wavelength.gt.2.4d0) go to 999
      if(scale.eq.0.0648d0.and.wavelength.lt.2.4d0) go to 999
c     
      print 120, sca_id, j, filter_id, scale, wavelength,
     &     psf_file(j)
 120  format('main:       sca',2x,i3,' filter ',i4,2x,a5,
     &     ' scale',2x,f6.4,' wavelength ',f7.4,
     &     2x,a180)
      PRINT *,'PA_DEGREES ', PA_DEGREES, 'PAUSE'
c
c***********************************************************************
c     
c     create scene for SCA sca_id in filter filter_id at dither
c     position  idither
c
c***********************************************************************
c     
      call sca_image(idither, dx(idither), dy(idither),
     &     pa_degrees,
     &     filename, noise_name,
     *     sca_id, module, brain_dead_test,
     *     x0, y0, pa_v3, osim_scale, scale,
     *     include_ktc, include_dark, include_readnoise, 
     &     include_reference,
     *     include_1_over_f, include_latents, 
     *     include_non_linear, include_cr, cr_mode, include_bg,
     *     include_stars, include_galaxies, nstars, ngal,
     *     bitpix, ngroups, nframe, nskip, tframe, tgroup,
     *     object, subarray, colcornr, rowcornr, 
     *     naxis1, naxis2,
     *     filter_id, wavelength, bandwidth, 
     *     system_transmission,
     *     photplam, photflam, stmag, abmag,
     *     background, icat_f,j, 
     *     psf_file(j), over_sampling_rate,
     *     noiseless, psf_add, ipc_add, verbose)
      go to 1000
 999  continue
      print *,' in guitarra.f'
      print 998, j, filterid(j)
 998  format('filter number:', i3,2x,'id: ',a20)
      print *,'scale : ',scale,' wavelength: ',wavelength
      print *,'something is inconsistent; possible error:'
      print *,'read_filter_parameters not reading appropriate',
     &     ' list; pause'
      read(*,'(A)')
 1000 continue
      stop
      end
c
