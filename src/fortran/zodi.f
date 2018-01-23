c
c-----------------------------------------------------------------------
c
      subroutine zodi(zodiacal_scale_factor, mirror_area, bkg, mode, nf,
     &     verbose)
c
c     Average zodiacal light contribution coming from Marcia's table
c     (only for W filters though) or from STScI ETC for low or 
c     average values. Assume gain = 1.0 when converting from e-
c     into photons
c
      implicit none
      double precision bkg, bkg_stsci, bkg_mjr, bkg_calc,
     *     zodiacal_scale_factor
      double precision effective_wl_nircam, width_nircam, response
      double precision pi, scale,dlam_over_lam, hplanck, pixel_area,
     *     stsci, f_nu, ratio, low, average, cee, f_lam, e_flux, 
     *     arc_sec_per_radian, mirror_area
      integer ii, nfilters, mode, npar, indx, nf, verbose
      character path_guitarra*100
c
      parameter(nfilters=54, npar=30)
c
      dimension bkg(nfilters), bkg_stsci(nfilters), bkg_mjr(nfilters),
     *     bkg_calc(npar)
      dimension effective_wl_nircam(nfilters), width_nircam(nfilters),
     *     response(nfilters)
c
      common /throughput/ effective_wl_nircam, width_nircam, response

c
c     MJR values for the background in photons/sec/pixel
c
      data bkg_mjr/    0.09d0, 0.09d0, 0.18d0, 0.18d0, 0.15d0, 0.15d0,
     *     0.00d0, 0.0d0, 
     *     0.16d0, 0.16d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0,
     *     0.00d0, 0.00d0, 0.00d0, 0.00d0,
     *     0.18d0, 0.18d0,
     *     0.00d0, 0.00d0, 0.00d0, 0.00d0,
     *     0.00d0, 0.00d0, 0.61d0, 0.61d0, 
     *     0.00d0, 0.00d0, 0.00d0, 0.00d0,0.00d0, 0.00d0, 
     *     0.43d0, 0.43d0, 
     *     0.00d0, 0.00d0, 0.00d0, 0.00d0,0.00d0, 0.00d0, 0.00d0,0.00d0,
     *     1.20d0, 1.2d0,
     *     0.00d0, 0.00d0, 0.00d0,0.00d0,0.00d0, 0.00d0, 0.00d0,0.00d0/
c
c     constants
c
      hplanck     = 6.626068963d-27 ! erg s 
      cee         = 2.99792458d10 ! cm/s
      pi          = dacos(-1.0d0)
      arc_sec_per_radian = 180.d0 * 3600.d0/pi 
      pixel_area  = 1.8d-4*1.8d-4 ! pixel area in cm**2
c
c     calculate the background from SED
c
      call zodi_background(nf,mirror_area, 
     *     zodiacal_scale_factor, bkg_calc, verbose)
      if(verbose.gt.0) then
         print 10
 10      format(4x,'Wavelength',5x,'STScI', 5x,'MJR model',5x,
     &        'ratio',8x,'f_nu', 12x,'f_lam',/,
     &        17x,'e-/sec/pix',2x,'e-/sec/pix',18x,'Jy',9x,
     &        'erg/(cm**2 sec A)')
      end if
c     
c     STSCI ETC background in e-/sec/arc_sec**2
c
      call getenv('GUITARRA_HOME',path_guitarra)
      open(1,file=path_guitarra(1:len_trim(path_guitarra))
     +            //'data/etc_stsci_sky_bkg_2015_02_12.dat')
c     
c     This is necessary as there are A/B filters, while
c     ETC uses only a generic filter
c
      read(1,*)
      indx = 0
      do ii = 1, 27
         read(1,20) low, average 
 20      format(6x,2(1x, f7.3))
         indx = indx + 1
         bkg_stsci(indx) =  average
         if(nf.eq.54) then
            indx = indx + 1
            bkg_stsci(indx) =  average
         end if
      end do
c
      do ii = 1, nf
         if(effective_wl_nircam(ii).gt.2.2d0) then
            scale = 0.0648d0
         else
            scale = 0.0317d0
         end if
c
c     make this conversion so ETC is in same units as MJR's (e-/sec/pixel)
c
         stsci = bkg_stsci(ii) * scale*scale
         if(mode .eq. 1) then
            e_flux = bkg_mjr(ii)
         end if
         if(mode .eq. 2) then
            e_flux = stsci
         end if
         if(mode.eq.3) then
            e_flux  =bkg_calc(ii)
         end if
c     convert from e-/sec/pixel to  e-/sec/sr/
         e_flux  = e_flux * (scale/arc_sec_per_radian)**2 
         e_flux  = e_flux / pixel_area     ! e-/sec/sr/cm**2
         e_flux  = e_flux/width_nircam(ii) ! e-/sec/sr/cm**2/micron
c
         f_lam   = e_flux * hplanck * cee/effective_wl_nircam(ii)
c
c     since
c     f_nu = (wl**2/cee) *f_lam = e_flux * hplanck * cee/wl * (wl**2/cee)
c
         f_nu     = e_flux * hplanck * effective_wl_nircam(ii)
         f_nu     = f_nu * 1.d23
         ratio = 0.d0
c         if(bkg_mjr(ii) .ne.0.0d0) ratio = stsci/bkg_mjr(ii)
         ratio = stsci/bkg_calc(ii)
c         bkg(ii)    = bkg_mjr(ii)
         bkg(ii)    = bkg_calc(ii)
         if(verbose.gt.0) then
            print 130, effective_wl_nircam(ii), stsci, bkg(ii), 
     *           ratio,f_nu, f_lam
 130        format(4(2x,f10.3),2(2x,1pe16.10))
         end if
      end do
      return
      end
