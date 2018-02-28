c
c----------------------------------------------------------------------     
c
      subroutine open_big_fits_cube(filename, iunit, nx, ny, nz,bitpix,
     *     naxis, naxes,
     *     nframe, tframe, groupgap, tgroup, ngroup, nskip,
     *     ra, dec, pa_degrees,
     *     object, partname, sca_id, module, filter, 
     *     photplam, photflam, stmag, abmag,
     *     subarray, colcornr, rowcornr, naxis1, naxis2,job,
     *     include_ktc, include_bg, include_cr, include_dark,
     *     include_latents, include_readnoise, include_non_linear,
     *     bias_value, readnoise, background, verbose)
c
      implicit none
      character telescop*20, instrume*20, filter*5,
     *     module*20, partname*5,comment*40, object*20      
      real image, tframe, tgroup
      integer status, bitpix, naxes, naxis, pcount, gcount, block,
     *     groupgap, sca_id, job, verbose
      integer include_ktc, include_bg, include_cr, include_dark,
     *     include_latents, include_readnoise, include_non_linear
      double precision bias_value, readnoise, background, exptime
      double precision total_time
      double precision photplam, photflam, stmag, abmag

c      double precision  bzero, bscale
      integer iunit, nx, ny, nz,  nframe, ngroup, nnn, nskip
      logical simple,extend
      character subarray*8
      logical dhas_subarray
      integer  colcornr, rowcornr, naxis1, naxis2
      character filename*120
c
      double precision sec, ut, jday, mjd
      double precision ra, dec, pa_degrees
      double precision expstart, expmid, expend, effexptm
      integer ih, im, month, day, year
      character  date_obs*10, time_obs*12, full_date*23,
     &      string*40, pw_filter*8, fw_filter*8, card*80
c
      parameter (nnn=2048)
c
      dimension image(nnn,nnn),naxes(3)
c
c     define the required primary array keywords
c
      if(verbose .ge.1) print *,'entered open_big_fits_cube'
      simple   = .true.
      extend   = .true.
      naxis    = 2
      if(subarray .eq. 'FULL') then 
         naxes(1) = nnn
         naxes(2) = nnn
         dhas_subarray = .False.
      else
         naxes(1) = naxes(1)
         naxes(2) = naxes(2)
         dhas_subarray = .True.
      end if
      if(nz.gt.1) then
         naxis    =  3
         naxes(3) = nz
      end if
c     
C     Get an unused Logical Unit Number to use to open the FITS file
c
      status = 0
      call ftgiou(iunit, status)
      if (status .gt. 0) call printerror(status)
      if(verbose.ge.2) print *,'open_big_fits_cube iunit', iunit
c      status = 0
c      call ftgiou(iunit, status)
c
c     delete previous version of the file, if it exists
c
      if(verbose.ge.2) print *,'open_big_fits_cube ftopen'
      call ftopen(iunit, filename, 1, block, status)
      if (status .eq. 0) then
         call ftdelt(iunit, status)
      else
c     clear the error message stack
         call ftcmsg
      end if
      status = 0
c
c     create the fits file
c
      if(verbose.ge.2) print *,'open_big_fits_cube ftinit'
      call ftinit(iunit, filename, 1, status)
      if (status .gt. 0) then 
         call printerror(status)
         print *,'pause: enter return to continue'
         read(*,'(A)')
      endif
      status =  0
c     
      if(verbose.ge.2) print *,'open_big_fits_cube ftphr'
      call ftphpr(iunit,simple, bitpix, naxis, naxes, 
     & 0,1,extend,status)
      if (status .gt. 0) then
         print *, 'simple,bitpix,naxis'
         call printerror(status)
      end if
      status =  0
c
c      call ftpkyj(iunit,'BSCALE',bscale,comment,status)
c       if (status .gt. 0) then
c          call printerror(status)
c          print *, 'BSCALE'
c       end if
c
      status =  0
c      call ftpkyj(iunit,'BZERO',bzero,comment,status)
c       if (status .gt. 0) then
c          call printerror(status)
c          print *, 'BZERO'
c       end if
      status =  0
c
c
      if(verbose.ge.2) print *,'open_big_fits_cube ftpkyj INC_KTC'
      comment = 'include KTC F(0) T(1)'
      call ftpkyj(iunit,'INC_KTC',include_ktc,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'INC_KTC'
          status = 0
       end if
c
      if(verbose.ge.2) print *,'open_big_fits_cube ftpkyj INC_RON'
      comment = 'include readnoise F(0) T(1)'
      call ftpkyj(iunit,'INC_RON',include_readnoise,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'INC_NRON'
          status = 0
       end if
c
      if(verbose.ge.2) print *,'open_big_fits_cube ftpkyj INC_BKG'
      comment = 'include background F(0) T(1)'
      call ftpkyj(iunit,'INC_BKG',include_bg,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'INC_BKG'
          status = 0
       end if
c
      if(verbose.ge.2) print *,'open_big_fits_cube ftpkyj INC_CR'
      comment = 'include Cosmic Rays F(0) T(1)'
      call ftpkyj(iunit,'INC_CR',include_dark,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'INC_CR'
          status = 0
       end if
c
      if(verbose.ge.2) print *,'open_big_fits_cube ftpkyj INC_DARK'
      comment = 'include darks F(0) T(1)'
      call ftpkyj(iunit,'INC_DARK',include_dark,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'INC_DARK'
          status = 0
       end if
c
      if(verbose.ge.2) print *,'open_big_fits_cube ftpkyj INC_LAT'
      comment = 'include latents F(0) T(1)'
      call ftpkyj(iunit,'INC_LAT',include_latents,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'INC_LAT'
          status = 0
       end if
c
      if(verbose.ge.2) print *,'open_big_fits_cube ftpkyj INC_NLIN'
      comment = 'include non-linearity F(0) T(1)'
      call ftpkyj(iunit,'INC_NLIN',include_non_linear,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'INC_NLIN'
          status = 0
       end if
c
      if(verbose.ge.2) print *,'open_big_fits_cube ftpkye KTC',
     &      bias_value
      comment = 'KTC value (e-)'
      call ftpkye(iunit,'KTC',real(bias_value),-7,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'KTC'
          status = 0
       end if
c
      if(verbose.ge.2) print *,'open_big_fits_cube ftpkye readnoise'
      comment = 'Readout noise (e-)'
      call ftpkye(iunit,'RDNOISE',real(readnoise),-7,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'RDNOISE'
          status = 0
       end if
c
      if(verbose.ge.2) print *,'open_big_fits_cube ftpkye background'
      comment = 'background (e-/sec/pixel)'
      call ftpkye(iunit,'BKG',real(background),-7,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'BKG'
          status = 0
       end if
c
c     exposure time for this ramp
c
      exptime = total_time(nframe, nskip, ngroup, 1,
     *     dble(tframe))
      if(exptime.le.0) then
         print *,'open_big_fits_cube: total_time.f', exptime
         print *,'negative time ?'
         stop
      endif
      if(verbose.ge.2) print *,'open_big_fits_cube ftpkye exptime',
     &     exptime
      comment = 'Exposure time'
      call ftpkye(iunit,'EXPTIME',real(exptime),-7,
     *     comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'EXPTIME'
      end if
      status =  0
c
c---------------------------
c
      call get_date_time(date_obs, time_obs)
      read(date_obs, 130) year, month, day
 130  format(i4,1x,i2,1x,i2)
      write(full_date,110) date_obs, time_obs
 110  format(a10,'T',a12)
      read(time_obs,120) ih, im, sec
 120  format(i2,1x,i2,1x,f6.3)
c
c     calculate UT and other parameters
c
      call julian_day(year, month, day, ut, jday, mjd)
c
      ut = ih + im/60.d0 + sec/3600.d0
      expstart = mjd + ut/24.d0
c
c----------------------------
c
c     coordinates
c
      comment = 'Target RA at mid time of exposure'
      call ftpkyd(iunit,'TARG_RA',ra,-15,
     *     comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'TARG_RA'
      end if
      status =  0
c
      comment = 'Target DEC at mid time of exposure'
      call ftpkyd(iunit,'TARG_DEC',dec,-15,
     *     comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'TARG_DEC'
      end if
      status =  0
c      
      comment = 'Position angle of V3-axis of JWST'
      call ftpkyd(iunit,'PA_V3', pa_degrees,-7,
     *     comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'PA_V3'
      end if
      status =  0
      
c
c     write more header keywords
c
      call write_nircam_keywords(iunit,nframe, tframe, groupgap,
     *      tgroup, ngroup, object, partname, sca_id,module, filter, 
     *     dhas_subarray, colcornr, rowcornr,job)
c      
      comment = 'Pivot wavelength Angstroms'
      call ftpkye(iunit,'PHOTPLAM',real(photplam*1.d04),7,
     *     comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'PHOTPLAM'
      end if
      status =  0
c
      comment = 'Flux for 1e s-1 in erg cm**-2 s-1 A-1'
      call ftpkye(iunit,'PHOTFLAM',real(photflam),-7,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'PHOTFLAM'
      end if
      status =  0
c
      comment = 'STMAG zeropoint'
      call ftpkye(iunit,'STMAG',real(stmag),-6,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'STMAG'
      end if
      status =  0
c
      comment = 'ABMAG zeropoint'
      call ftpkye(iunit,'ABMAG',real(abmag),-6,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'ABMAG'
      end if
      status =  0
c
      return
      end
