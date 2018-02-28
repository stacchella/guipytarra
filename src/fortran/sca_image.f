c
c     Create data cube for an SCA at a given dither position.
c     Start by reading the bias and gain images for this SCA.
c     If this is the first dither, ignore reading the latent 
c     image.
c
c     cnaw 2015-01-27
c     Steward Observatory, University of Arizona
c
c     modified to make it compatible STScI data model
c     cnaw 2017-05-01
c     Steward Observatory, University of Arizona
c     
      subroutine sca_image(idither, ra_dithered, dec_dithered, 
     *     pa_degrees,
     *     filename, noise_name,
     *     sca_id, module, brain_dead_test, 
     *     xc, yc, osim_scale,scale,
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
      integer cube, overlap
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
      character subarray*(*)
c
      parameter (nnn=2048, max_order=7, overlap=30)
c
      dimension dark_mean(10), dark_sigma(10), gain(10),
     *     read_noise(10),  even_odd(8)
c
c     images
c
      dimension cube(nnn,nnn,overlap)
      dimension base_image(nnn,nnn), naxes(3)
      dimension accum(nnn,nnn),image(nnn,nnn),latent_image(nnn,nnn),
     *     well_depth(nnn,nnn), linearity(nnn,nnn,max_order),
     *     bias(nnn,nnn), gain_image(nnn,nnn), scratch(nnn,nnn)
c
c     images
c
      common /history/ cube
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
      job = 1
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
      call get_sca_id(sca_id, partname)
      indx = sca_id - 480
c     
c     define FITS keywords
c
c     How to tag this to a previous image ?
c
      write(latent_file, 1120) filter_id, iabs(sca_id)
 1120 format('latent_',a5,'_',i3.3,'.fits')
c
c     WCS keywords ; to take distortions into account need to
c     include the SIP keywords. For now assume all is plane.
c
      x_sca = 1024.d0
      y_sca = 1024.d0
c
      call wcs_keywords(sca_id, x_sca, y_sca, xc, yc, osim_scale,
     *     ra_dithered, dec_dithered,  pa_degrees,verbose)
      call osim_coords_from_sca(sca_id, x_sca, y_sca, x_osim, y_osim)
      call sca_to_ra_dec(sca_id, 
     *     ra_dithered, dec_dithered,
     *     ra_sca, dec_sca, pa_degrees, 
     *     xc, yc, osim_scale, x_sca, y_sca)
      print *,'crval1, crval2', crval1, crval2
      print *,'crpix1, crpix2', crpix1, crpix2
      print *,'ra_dithered, dec_dithered',ra_dithered, dec_dithered
      print *,'ra_sca, dec_sca',ra_sca, dec_sca
      print *,'pa_degrees, osim_scale',pa_degrees, osim_scale
      print *,'xc, yc', xc, yc
c      print *,'stop at sca_image'
c      stop
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
      if(subarray .eq. 'FULL') then
         naxes(1) =  nnn
         naxes(2) =  nnn
      else
         naxes(1) = naxis1
         naxes(2) = naxis2
      end if
      naxes(3) =  ngroups
      if(verbose.gt.1) print *, 'open_big_fits_cube for ', sca_id
c
c     read gain, bias and well-depth images; set the average and sigma values
c     for the darks.
c
      call read_sca_calib(sca_id, verbose)

c     
c     if latents are being considered read the previous image
c     which will be attenuated as time progresses.
c
      if(include_latents .eq.1 .and.idither .gt.1) then
         call  read_fits(latent_file,latent_image, nlx, nly)
         if(nlx .ne. n_image_x .or. nly .ne. n_image_y) then
            print *,' size of latent image does not match:',
     *           nlx, nly, 'image ', n_image_x, n_image_y
            include_latents = 0
            print *,' latents will be ignored'
         end if
      end if
c     
c     Read the cosmic ray distribution if set;
c     The sca id is required to set the wavelength for the
c     M. Robberto cosmic ray models
c
      if(include_cr .eq. 1) then
         call cr_distribution(cr_mode, sca_id)
      end if

c
c     Initialise image and 
c     create baseline image (bias + ktc)
c     
      if(include_ktc.eq.0 .and.noiseless .eqv. .false.) then
         bias_value = 10.0d0
      end if
      if(include_ktc.eq.1) then
         bias_value = ktc
      end if
c
c     Brain-dead test
c
      if(brain_dead_test.eq.1) then
         print *,' ktc = 1000 '
         do j = 1, n_image_y
            do i = 1, n_image_x
               image(i,j)      =    0.0
               base_image(i,j) = 1000.0
            end do
         end do
      end if
c
c     Add ktc for random SCA (i.e., not a NIRCam one)
c
      if(noiseless .eqv. .true.) then
         do j = 1, n_image_y
            do i = 1, n_image_x
               image(i,j)      = 0.0
               base_image(i,j) = 0.0
            end do
         end do
      else
         if(include_ktc.eq.-1 .or.sca_id.eq.0) then 
            do j = 1, n_image_y
               do i = 1, n_image_x
                  image(i,j)      = 0.0
                  base_image(i,j) = 0.0
               end do
            end do
         else
c
c     otherwise
c
            do j = 1, n_image_y
               do i = 1, n_image_x
                  image(i,j)     = 0.0
c
c     fix NaNs or negative values
c
                  if(bias(i,j).ne.bias(i,j).or.bias(i,j).le. 0.0) then
                     variate = 0.0 
                  else
                     variate = real(zbqlpoi(dble(bias(i,j))))
                  end if
                  base_image(i,j) = variate + voltage_offset ! units [ADU] 
               end do
            end do
         end if
      end if
c
      print *, 'image has been initialised'
c
c=======================================================================
c
c     Initialise output image and write FITS keywords
c
      call open_big_fits_cube(filename, iunit, 
     *     n_image_x, n_image_y, ngroups, bitpix,
     *     naxis, naxes,
     *     nframe, tframe, groupgap, tgroup, ngroups, nskip,
     *     ra_dithered, dec_dithered, pa_degrees,
     *     object, partname, sca_id, module, filter_id,
     *     photplam, photflam, stmag, abmag,
     *     subarray,  colcornr, rowcornr,naxis1, naxis2, job,
     *     include_ktc, include_bg, include_cr, include_dark,
     *     include_latents, include_readnoise, include_non_linear,
     *     bias_value, read_noise(indx), background, verbose)
      if(verbose.gt.1) print *, 'open_big_fits_cube done'

c
c     Loop through the number of groups
c
      do k = 1, ngroups
         if(verbose.ge.1) then
            print 1140,idither, k, ngroups,nframe,nskip
 1140       format('dither, group, ngroups,nframe, nskip ', 7I6)
         end if
c
c     loop through the cycles
c
         frame = 0 
         call clear_accum
         do loop = 1, nframe+nskip
            if(verbose .ge. 2) then
               if(loop.le.nframe) then
                  print *,'dither',idither,' group ', k , 
     *                 ' of ',ngroups,', read ',loop,
     *                 ' of  nframe+nskip =',nframe+nskip,
     *                 ' nskip ',nskip
               else
                  print *,'dither',idither,' group ', k , 
     *                 ' of ',ngroups,', skip ',loop-nframe,
     *                 ' of  nframe+nskip =',nframe+nskip,
     *                 ' nskip ',nskip
               end if
            end if
c     
c     Keep on accumulating signal whether the frame is read out or not
c     (in the "image" matrix). The "accum" matrix will contain the
c     scene to which readout noise is added.
c
c     add stars
c     [e-]
            if(include_stars.eq.1 .and. nstars.gt.0) then
               if(verbose.ge.2) print *, 'sca_image: add_stars '
               call add_stars(ra_dithered, dec_dithered,pa_degrees,
     *              xc, yc, osim_scale, sca_id, filter_index, seed,
     *              subarray, colcornr, rowcornr, naxis1, naxis2,
     *              wavelength, bandwidth,system_transmission,
     *              mirror_area,integration_time, in_field, 
     *              noiseless, psf_add, ipc_add, verbose)
               if(verbose.ge.2) then
                  print *, 'added ',in_field, ' stars of ', nstars
               end if
            end if
c     
c     add galaxies
c     [e-]
            if(include_galaxies .eq. 1 .and. ngal .gt. 0) then
               if(verbose.ge.2)print *, 'sca_image: add_galaxies',
     *              wavelength, bandwidth, system_transmission
               call add_modelled_galaxy(sca_id,
     *              ra_dithered, dec_dithered, pa_degrees,
     *              xc, yc, osim_scale, icat_f,
     *              ngal, scale,
     *              wavelength, bandwidth, system_transmission, 
     *              mirror_area, integration_time, seed, in_field,
     &              noiseless, psf_add, ipc_add,verbose)
               if(verbose.ge.2) then
                  print *, 'sca_image: added', in_field,' galaxies of ',
     &                 ngal
               end if
            end if
c     
c     add sky background [e-]
c     
            if(include_bg .eq. 1) then 
               if(verbose.ge.2) print *, 'going to add sky background',
     &              background,' e-/sec/pixel'
               call add_sky_background(background,
     *              subarray, colcornr, rowcornr, naxis1, naxis2,
     *              integration_time, noiseless, verbose)
            end if
c
c     add cosmic rays [e-]
c     
            if(include_cr .eq. 1) then 
               call add_modelled_cosmic_rays(n_image_x, n_image_y,
     *              cr_mode, subarray, naxis1, naxis2, integration_time,
     *              ipc_add, verbose)
c     *              subarray, colcornr, rowcornr, naxis1, naxis2)
            end if
c     
c     add dark current  [e-]
c
            if(include_dark .eq. 1) then
               if(verbose.ge.2) print *, 'add dark'
               call add_dark(brain_dead_test,
     *              subarray, colcornr, rowcornr, naxis1, naxis2,
     *              dark_mean, dark_sigma, integration_time)
            end if
c     
c     Add latent charge [e-]
c     
            if(include_latents .eq. 1) then 
               time_step = (k-1) * (nframe+nskip) +loop
               call add_latents(sca_id, time_step, integration_time,
     *              decay_rate, time_since_previous)
            end if
c
c     Add charge to reference pixels (and whole image) [e-]
c
            if(include_reference.eq.1) then
               if(verbose.gt.1) print *,'add reference pixels'
               call add_reference_pixels(read_noise, even_odd,
     &              subarray,colcornr, rowcornr, naxis1, naxis2)
            end if
c
            if(loop.le.nframe) then
c
c     add readnoise; however, readnoise should not be accumulated
c     as it is a measurement error, not an additive error that is
c     aggregated to the signal. This is added to the accumulated image,
c     as this will correspond to a "measured" quantity. The
c     same holds true for 1/F noise
c     read_noise units [e-];    1/f noise  units [e-]
c
               if (include_readnoise .eq. 1) then 
                  level = (k-1)*(nskip+nframe) + loop
                  if(verbose.ge.2) then
                     print *, 'sca_image: add read noise ',
     &                  '(loop, k, level)',loop, k, level,
     &                    read_noise(indx)
                  end if
                  call add_read_noise(brain_dead_test,read_noise, 
     *                 subarray,colcornr, rowcornr, naxis1, naxis2)
               end if
               if(include_1_over_f.eq.1) then
c                  write(noise_name,900) sca_id
c 900              format('ng_hxrg_noise_',i3,'.fits')
                  call add_one_over_f_noise(noise_name, level,
     *                 subarray,colcornr, rowcornr, naxis1, naxis2)
               end if
               call coadd
            end if
c     
c     if this is the last frame read in a group, calculate the average 
c     
            if(loop .eq. nframe) then
               call divide_image(nframe)
c     
c     undo the linearity correction and convert into ADU
c     
               if(verbose.ge.2) then
                  print *,'SCA_image : gain(indx)',indx, gain(indx)
               end if
               if(indx.gt.0) then
                  call linearity_incorrect(include_non_linear,
     *                 gain(indx), subarray,colcornr, rowcornr,
     *                  naxis1, naxis2, verbose)
               else
                  do j = 1, naxis2
                     do i = 1, naxis1
                        scratch(i,j) = image(i,j)
                     end do
                  end do
               end if
c
c     add baseline noise [ADU]
c
               call add_baseline(
     &              subarray,colcornr, rowcornr, naxis1, naxis2)
c
c     write to FITS data cube
c
               if(verbose.ge.2) print *,'sca_image : write image plane'
               call write_frame(filename, iunit, bitpix, k, 
     *              naxis, naxes)
            endif
c     
c  if this is the last read of a group, there are no more frames to skip
c  exit all integrations
c
            if(loop .eq. nframe .and. k .eq. ngroups) go to 1000
         end do                 ! close loop on nframe+ nskip
      end do                    ! close loop on ngroups
c
 1000 continue
c
      if(verbose .ge. 2) print *,'dither: exit main loop'
c
      call closefits(iunit)
c
c     write image with accumulated counts but without the read noise
c
      write(filename,1100) filter_id,iabs(sca_id),idither
 1100 format('sim_',a5,'_',i3.3,'_',i3.3,'.fits')
      ibitpix   = -32
      call write_float_2d_image(filename, image, n_image_x, n_image_y,
     *     ibitpix, nframe, tframe, nskip, tgroup, ngroups, object, 
     *     partname, sca_id, module, filter_id,
     *     subarray, colcornr, rowcornr, naxis1, naxis2, job)
c
c     store/update image for latents
c
      if(include_latents .eq. 1) then
         do i = 1, n_image_x
            do j = 1, n_image_y
               latent_image(i,j) = image(i,j)
            end do
         end do
         ibitpix = -32
         call write_float_2d_image(latent_file, latent_image,
     *        n_image_x, n_image_y,
     *        ibitpix, nframe, tframe, nskip, tgroup, ngroups, object, 
     *        partname, sca_id,module, filter_id,
     *        subarray, colcornr, rowcornr, naxis1, naxis2, job)
      endif
      if(verbose.ge.1) then
         print * ,'end of dither', sca_id, idither, n_image_x, 
     *        n_image_y, ngroups
      end if
c
c      write(filename,1200) filter_id,iabs(sca_id),idither
c 1200 format('history_',a5,'_',i3.3,'_',i3.3,'.fits')
c      ibitpix = 32
c      call write_int_3d_image(filename, cube, n_image_x, n_image_y,
c     *     ibitpix,
c     *     nframe, tframe, groupgap, tgroup, overlap, object, partname, 
c     *     partname, sca_id,module, filter_id,
c     *     subarray, colcornr, rowcornr, naxis1, naxis2, job)
c
      return
      end
