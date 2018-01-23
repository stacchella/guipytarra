c
c     Make a sky image from catalogue
c
c     cnaw 2017-04-26
c     Steward Observatory, University of Arizona
c     
      subroutine sky_image(idither, ra_dithered, dec_dithered, 
     *     pa_degrees,
     *     filename, noise_name,
     *     sca_id, module, brain_dead_test, 
     *     xc, yc, pa_v3, osim_scale,scale,
     *     include_ktc, include_dark, include_readnoise, 
     *     include_reference,
     *     include_1_over_f, include_latents, include_non_linear,
     *     include_cr, cr_mode, include_bg,
     *     include_stars, include_galaxies, nstars, ngal,
     *     bitpix, ngroups, nframe, nskip, tframe, tgroup, object,
     *     subarray, colcornr, rowcornr, naxis1, naxis2,
     *     filter_id, wavelength, bandwidth, system_transmission,
     *     photplam, photflam, stmag, abmag,
     *     background, icat_f,filter_index, psf_file, 
     *     over_sampling_rate, noiseless, psf_add,
     *     ipc_add, verbose)
      
      implicit none
c
      double precision mirror_area, integration_time, gain,
     *     decay_rate, time_since_previous, read_noise, 
     *     dark_mean, dark_sigma, ktc, voltage_offset
      double precision wavelength, bandwidth, system_transmission, 
     *     background, scale, bias_value
      double precision  even_odd
      double precision x_sca, y_sca, ra_sca, dec_sca
      double precision xc, yc, pa_v3, q, posang, ra_dithered,
     *     dec_dithered, osim_scale, x_osim, y_osim,
     *     pa_degrees
      double precision linearity_gain,  lincut, well_fact
      double precision zero_point, tol, eps
c
c     images are either real*4 or integer
c
      real  image,  accum, latent_image, base_image, bias, well_depth,
     *     gain_image, linearity, variate, scratch
      double precision photplam, photflam, stmag, abmag
c
      integer max_order, order,over_sampling_rate
      integer zbqlpoi
      integer filter_index, icat_f
c
c     parameters
c
      real tframe, tgroup
c      real  
      double precision
     *     equinox, crpix1, crpix2, crval1, crval2, cdelt1,cdelt2,
     *     cd1_1, cd1_2, cd2_1, cd2_2
c
      integer colcornr, rowcornr, naxis1, naxis2
c
      integer frame,  bitpix, groupgap, sca_id, ibitpix,iunit
      integer time_step, verbose, naxes, n_image_x, n_image_y,idither,
     *     nskip, naxis, ngroups, nframe, job, indx
c
      logical ipc_add, psf_add, noiseless
      integer cr_mode, seed
      integer include_ktc, include_bg, include_cr, include_dark,
     *     include_latents, include_readnoise, include_non_linear,
     *     include_stars, include_galaxies, include_cloned_galaxies,
     &     brain_dead_test, include_1_over_f, include_reference

      integer i, j, k, loop, nlx, nly, level
      integer nnn, nstars, ngal, in_field
      character filename*120, latent_file*120, psf_file*120
      character object*20, partname*5, module*20, filter_id*5
      character noise_name*120
c     
      logical subarray
c
      parameter (nnn=2048, max_order=7)
c
      dimension dark_mean(10), dark_sigma(10), gain(10),
     *     read_noise(10),  even_odd(8)
c
c     images
c
      dimension base_image(nnn,nnn), naxes(3)
      dimension accum(nnn,nnn),image(nnn,nnn),latent_image(nnn,nnn),
     *     well_depth(nnn,nnn), linearity(nnn,nnn,max_order),
     *     bias(nnn,nnn), gain_image(nnn,nnn), scratch(nnn,nnn)
c
c     images
c
      common /gain_/ gain_image
      common /base/ base_image
      common /latent/ latent_image
      common /images/ accum, image, n_image_x, n_image_y
      common /scratch_/ scratch
      common /well_d/ well_depth, bias, linearity,
     *     linearity_gain,  lincut, well_fact, order
c
      common /wcs/ equinox, crpix1, crpix2, crval1, crval2,
     *     cdelt1,cdelt2, cd1_1, cd1_2, cd2_1, cd2_2
c
c     Parameters
c
      common /parameters/ mirror_area, integration_time, gain,
     *     decay_rate, time_since_previous, read_noise, 
     *     dark_mean, dark_sigma, ktc, voltage_offset 
c                    
      q   = dacos(-1.0d0)/180.d0
      tol = 1.d-10
      eps = 1.d-16
c
      if(verbose.gt.0) then
         print 10,sca_id, filter_index, filter_id, idither
 10      format('sca_image:  sca ', i4,' filter ', i4, 2x, a5,
     *        ' dither ',i4)
      end if
c
c     translate SCA id into partname. Good to keep track that the correct
c     files are being read
c
c      call get_sca_id(sca_id, partname)
c      indx = sca_id - 480
c     
cc     define FITS keywords
cc
cc     How to tag this to a previous image ?
cc
c      write(latent_file, 1120) filter_id, iabs(sca_id)
c 1120 format('latent_',a5,'_',i3.3,'.fits')
c
c     WCS keywords ; to take distortions into account need to
c     include the SIP keywords. For now assume all is plane.
c
      x_sca = 1024.d0
      y_sca = 1024.d0
      call wcs_keywords(sca_id, x_sca, y_sca, xc, yc, osim_scale,
     *     ra_dithered, dec_dithered,  pa_degrees,verbose)
      call osim_coords_from_sca(sca_id, x_sca, y_sca, x_osim, y_osim)
      call sca_to_ra_dec(sca_id, 
     *     ra_dithered, dec_dithered,
     *     ra_sca, dec_sca, pa_degrees, 
     *     xc, yc, osim_scale, x_sca, y_sca)
c
c
c     Read PSF
c
      if(verbose.gt.1) then
         PRINT *,'SCA_IMAGE:', sca_id, ra_dithered, dec_dithered, 
     *        x_sca, y_sca, ra_sca, dec_sca
         print 1130, psf_file
      end if
 1130 format('read psf ', a120)
      call read_psf(psf_file, verbose)
      if(verbose.gt.2) print *,'psf has been read'
c
c     This ensures that bitpix will be 16 and bzero=32K, bscale = 1
c     for unsigned integers
c
      bitpix    = 20
c
      groupgap =   nskip
      naxis    =   3
c
c     set image sizes for sub-array cases
c
      if(subarray .eqv. .true.) then
         naxes(1) = naxis1
         naxes(2) = naxis2
      else
         naxes(1) =  nnn
         naxes(2) =  nnn
      end if
c      naxes(3) =  ngroups
c      if(verbose.gt.1) print *, 'open_big_fits_cube for ', sca_id
c
c     read gain, bias and well-depth images; set the average and sigma values
c     for the darks.
c
c      call read_sca_calib(sca_id, verbose)
c
cc     
cc     if latents are being considered read the previous image
cc     which will be attenuated as time progresses.
cc
c      if(include_latents .eq.1 .and.idither .gt.1) then
c         call  read_fits(latent_file,latent_image, nlx, nly)
c         if(nlx .ne. n_image_x .or. nly .ne. n_image_y) then
c            print *,' size of latent image does not match:',
c     *           nlx, nly, 'image ', n_image_x, n_image_y
c            include_latents = 0
c            print *,' latents will be ignored'
c         end if
c      end if
cc     
cc     Read the cosmic ray distribution if set;
cc     The sca id is required to set the wavelength for the
cc     M. Robberto cosmic ray models
cc
c      if(include_cr .eq. 1) then
c         call cr_distribution(cr_mode, sca_id)
c      end if
c
cc
cc     Initialise image and 
cc     create baseline image (bias + ktc)
c     
c      if(include_ktc.eq.0) then
c         bias_value = 10.0d0
c      end if
c      if(include_ktc.eq.1) then
c         bias_value = ktc
c      end if
cc
cc     Brain-dead test
cc
c      if(brain_dead_test.eq.1) then
c         print *,' ktc = 1000 '
c         do j = 1, n_image_y
c            do i = 1, n_image_x
c               image(i,j)      =    0.0
c               base_image(i,j) = 1000.0
c            end do
c         end do
c      end if
cc
cc     Add ktc for random SCA (i.e., not a NIRCam one)
cc
c      if(noiseless .eqv. .true.) then
c         do j = 1, n_image_y
c            do i = 1, n_image_x
c               image(i,j)      = 0.0
c               base_image(i,j) = 0.0
c            end do
c         end do
c      else
c         if(include_ktc.eq.-1 .or.sca_id.eq.0) then 
c            do j = 1, n_image_y
c               do i = 1, n_image_x
c                  image(i,j)      = 0.0
c                  base_image(i,j) = 0.0
c               end do
c            end do
c         else
cc
cc     otherwise
cc
c      do j = 1, n_image_y
c         do i = 1, n_image_x
c            image(i,j)     = 0.0
c     
c     fix NaNs or negative values
c
c                  if(bias(i,j).ne.bias(i,j).or.bias(i,j).le. 0.0) then
c                     variate = 0.0 
c                  else
c                     variate = real(zbqlpoi(dble(bias(i,j))))
c                  end if
c                  base_image(i,j) = variate + voltage_offset ! units [ADU] 
c               end do
c            end do
c         end if
c      end if
c
      do j = 1, n_image_y
         do i = 1, n_image_x
            image(i,j) = 0.0
         end do
      end do
c
      print *, 'image has been initialised'
cc
cc=======================================================================
cc
cc     Initialise output image and write FITS keywords
cc
c      call open_big_fits_cube(filename, iunit, 
c     *     n_image_x, n_image_y, ngroups, bitpix,
c     *     naxis, naxes,
c     *     nframe, tframe, groupgap, tgroup, ngroups, nskip,
c     *     object, partname, sca_id, module, filter_id,
c     *     photplam, photflam, stmag, abmag,
c     *     subarray,  colcornr, rowcornr,naxis1, naxis2, job,
c     *     include_ktc, include_bg, include_cr, include_dark,
c     *     include_latents, include_readnoise, include_non_linear,
c     *     bias_value, read_noise, background, verbose)
c      if(verbose.gt.1) print *, 'open_big_fits_cube done'
c
c
c     Loop through the number of groups
c
c      do k = 1, ngroups
c         if(verbose.ge.1) then
c            print 1140,idither, k, ngroups,nframe,nskip
c 1140       format('dither, group, ngroups,nframe, nskip ', 7I6)
c         end if
cc
cc     loop through the cycles
cc
c         frame = 0 
         call clear_accum
c         do loop = 1, nframe+nskip
c            if(verbose .ge. 2) then
c               if(loop.le.nframe) then
c                  print *,'dither',idither,' group ', k , 
c     *                 ' of ',ngroups,', read ',loop,
c     *                 ' of  nframe+nskip =',nframe+nskip,
c     *                 ' nskip ',nskip
c               else
c                  print *,'dither',idither,' group ', k , 
c     *                 ' of ',ngroups,', skip ',loop-nframe,
c     *                 ' of  nframe+nskip =',nframe+nskip,
c     *                 ' nskip ',nskip
c               end if
c            end if
cc     
cc     Keep on accumulating signal whether the frame is read out or not
cc     (in the "image" matrix). The "accum" matrix will contain the
cc     scene to which readout noise is added.
cc
cc     add stars
cc     [e-]
         if(include_stars.eq.1 .and. nstars.gt.0) then
            if(verbose.ge.2) print *, 'sca_image: add_stars '
            call add_stars(ra_dithered, dec_dithered,pa_degrees,
     *           xc, yc, osim_scale, sca_id, filter_index, seed,
     *           subarray, colcornr, rowcornr, naxis1, naxis2,
     *           wavelength, bandwidth,system_transmission,
     *           mirror_area,integration_time, in_field, 
     *           noiseless, psf_add, ipc_add, verbose)
            if(verbose.ge.2) then
               print *, 'added ',in_field, ' stars of ', nstars
            end if
         end if
c     
c     add galaxies
c     [e-]
         if(include_galaxies .eq. 1 .and. ngal .gt. 0) then
            if(verbose.ge.2)print *, 'sca_image: add_galaxies',
     *           wavelength, bandwidth, system_transmission
            call add_modelled_galaxy(sca_id,
     *           ra_dithered, dec_dithered, pa_degrees,
     *           xc, yc, osim_scale, icat_f,
     *           ngal, scale,
     *           wavelength, bandwidth, system_transmission, 
     *           mirror_area, integration_time, seed, in_field,
     &           noiseless, psf_add, ipc_add,verbose)
            if(verbose.ge.2) then
               print *, 'sca_image: added', in_field,' galaxies of ',
     &              ngal
            end if
         end if
c     
c     add sky background [e-]
c     
         if(include_bg .eq. 1) then 
            if(verbose.ge.2) print *, 'going to add sky background',
     &           background,' e-/sec/pixel'
            call add_sky_background(background,
     *           subarray, colcornr, rowcornr, naxis1, naxis2,
     *           integration_time, verbose)
         end if
c     
c     add cosmic rays [e-]
c     
c         if(include_cr .eq. 1) then 
c            call add_modelled_cosmic_rays(n_image_x, n_image_y,
c     *           cr_mode, subarray, naxis1, naxis2, integration_time)
cc     *              subarray, colcornr, rowcornr, naxis1, naxis2)
c         end if
cc     
cc     add dark current  [e-]
cc
c            if(include_dark .eq. 1) then
c               if(verbose.ge.2) print *, 'add dark'
c               call add_dark(brain_dead_test,
c     *              subarray, colcornr, rowcornr, naxis1, naxis2,
c     *              dark_mean, dark_sigma, integration_time)
c            end if
cc     
cc     Add latent charge [e-]
cc     
c            if(include_latents .eq. 1) then 
c               time_step = (k-1) * (nframe+nskip) +loop
c               call add_latents(sca_id, time_step, integration_time,
c     *              decay_rate, time_since_previous)
c            end if
cc
cc     Add charge to reference pixels (and whole image) [e-]
cc
c            if(include_reference.eq.1) then
c               if(verbose.gt.1) print *,'add reference pixels'
c               call add_reference_pixels(read_noise, even_odd,
c     &              subarray,colcornr, rowcornr, naxis1, naxis2)
c            end if
cc
c            if(loop.le.nframe) then
cc
cc     add readnoise; however, readnoise should not be accumulated
cc     as it is a measurement error, not an additive error that is
cc     aggregated to the signal. This is added to the accumulated image,
cc     as this will correspond to a "measured" quantity. The
cc     same holds true for 1/F noise
cc     read_noise units [e-];    1/f noise  units [e-]
cc
c               if (include_readnoise .eq. 1) then 
c                  level = (k-1)*(nskip+nframe) + loop
c                  if(verbose.ge.2) then
c                     print *, 'sca_image: add read noise ',
c     &                  '(loop, k, level)',loop, k, level,
c     &                    read_noise(indx)
c                  end if
c                  call add_read_noise(brain_dead_test,read_noise, 
c     *                 subarray,colcornr, rowcornr, naxis1, naxis2)
c               end if
c               if(include_1_over_f.eq.1) then
cc                  write(noise_name,900) sca_id
cc 900              format('ng_hxrg_noise_',i3,'.fits')
c                  call add_one_over_f_noise(noise_name, level,
c     *                 subarray,colcornr, rowcornr, naxis1, naxis2)
c               end if
c               call coadd
c            end if
cc     
cc     if this is the last frame read in a group, calculate the average 
cc     
c            if(loop .eq. nframe) then
c               call divide_image(nframe)
cc     
cc     undo the linearity correction and convert into ADU
cc     
c               if(verbose.ge.2) then
c                  print *,'SCA_image : gain(indx)',indx, gain(indx)
c               end if
c               if(indx.gt.0) then
c                  call linearity_incorrect(include_non_linear,
c     *                 gain(indx), subarray,colcornr, rowcornr,
c     *                  naxis1, naxis2, verbose)
c               else
c                  do j = 1, naxis2
c                     do i = 1, naxis1
c                        scratch(i,j) = image(i,j)
c                     end do
c                  end do
c               end if
cc
cc     add baseline noise [ADU]
cc
c               call add_baseline(
c     &              subarray,colcornr, rowcornr, naxis1, naxis2)
cc
cc     write to FITS data cube
cc
c               if(verbose.ge.2) print *,'sca_image : write image plane'
c               call write_frame(filename, iunit, bitpix, k, 
c     *              naxis, naxes)
c            endif
cc     
cc  if this is the last read of a group, there are no more frames to skip
cc  exit all integrations
cc
c            if(loop .eq. nframe .and. k .eq. ngroups) go to 1000
c         end do                 ! close loop on nframe+ nskip
c      end do                    ! close loop on ngroups
cc
c 1000 continue
c
c      if(verbose .ge. 2) print *,'dither: exit main loop'
cc
c      call closefits(iunit)
c
c     write image with accumulated counts but without the read noise
c
      write(filename,1100) filter_id,iabs(sca_id),idither
 1100 format('sim_',a5,'_',i3.3,'_',i3.3,'.fits')
      ibitpix   = -32
      do j = 1, n_image_y
         do i = 1, n_image_x
            image(i,j) = image(i,j)/integration_time
         end do
      end do
      call write_float_2d_image(filename, image, n_image_x, n_image_y,
     *     ibitpix, nframe, tframe, nskip, tgroup, ngroups, object, 
     *     partname, sca_id, module, filter_id,
     *     subarray, colcornr, rowcornr, naxis1, naxis2, job)
c
c     store/update image for latents
c
c      if(include_latents .eq. 1) then
c         do i = 1, n_image_x
c            do j = 1, n_image_y
c               latent_image(i,j) = image(i,j)
c            end do
c         end do
c         ibitpix = -32
c         call write_float_2d_image(latent_file, latent_image,
c     *        n_image_x, n_image_y,
c     *        ibitpix, nframe, tframe, nskip, tgroup, ngroups, object, 
c     *        partname, sca_id,module, filter_id,
c     *        subarray, colcornr, rowcornr, naxis1, naxis2, job)
c      endif
      if(verbose.ge.1) then
         print * ,'end of dither', sca_id, idither, n_image_x, 
     *        n_image_y, ngroups
      end if
      return
      end
