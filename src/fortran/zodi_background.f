c-------------------------------------------------------------------------
c
c     Calculate the background due to zodiacal light
c     using expressions from M. Rieke, 2005 jwst_calc_003894
c
c     2015-02-16
c
      subroutine zodi_background(nf, mirror_area, zodiacal_scale_factor,
     *     results, verbose)
      implicit none
      double precision mirror_area, zodiacal_scale_factor, results
      double precision pi, cee, h, boltzmann, hc_over_k, two_hcc,
     *     aperture, arc_sec_per_radian, dwl, dwl_cm, wl_cm, wl0, wlf, 
     *     wl_nominal, pixel, pixel_sr, wavelength, throughput,
     *     photons, photon_flux, electron_flux, electrons, events,
     *     factor, eflux, e_photon, per_aperture, eflux_total, f_nu,
     *     f_lambda
      double precision e_scatter, t_scatter, tau_thermal, t_thermal,
     *     to_si_per_micron, exp1, denom1, term1_cgs, term1_si,
     *     exp2, denom2, term2_cgs, term2_si, sb, bkg
      double precision filters, filtpars
c
      integer nf, nnn, mmm, lll, npar, npts, index, k, verbose

      character filterid*20, temp*20
c
      parameter (nnn=25000,mmm=200,lll=1000,npar=30)
c
      dimension results(npar)
c
      dimension filters(2,nnn,mmm),filtpars(mmm,npar), filterid(mmm)
      common /filter/filters, filtpars, filterid
c     
c     CGS units
c
      pi  = dacos(-1.00d0)
      cee = 2.99792458d10       ! cm/s
      h   = 6.626068963d-27     ! erg s 
      boltzmann = 1.38065042d-16 ! erg/K
c
      hc_over_k = h*cee/boltzmann
      two_hcc   = 2.d0 * h * cee*cee
c
      aperture = 2.5d0
      arc_sec_per_radian = 180.d0 * 3600.d0/pi 
c
c     QE and instrument efficiencies are included in the
c     filter throughput curves
c
      e_scatter = 3.95d-14      ! ?
      t_scatter = 5301.d0       ! K
c
      tau_thermal = 2.79d-8     ! ?
      t_thermal   = 282.d0      ! K
c
c     set the zodiacal light level to 20% higher than minimum
c     (Rieke, 2005 jwst_calc_003894.pdf)
c      zodiacal_scale_factor = 1.20d0
c
c     this seems to match the ST ETC low background
c      zodiacal_scale_factor = 0.858801d0
c
c     this gives a better match with MJR's values
c      zodiacal_scale_factor = 2.00302d0
      
      npts       = 2000
      wavelength = 0.5d0
      dwl        = (5.5d0 - 0.5d0)/npts
      dwl        = (5.8d0 - 0.5d0)/npts
      dwl_cm     = dwl * 1.d-04
c     converts from ergs/cm2/cm/sec to W/m2/micron:
      to_si_per_micron = 1.0d-7 
c
      if(verbose.gt.0) then
         print 111, aperture
 111     format('zodi_background: aperture radius', f6.2,
     &    /,11x,'wavelength  e-/pixel/sec events/aperture/sec',
     &        2x,'e-/arcsec**2/sec   MJy/sr     W/(m**2 um sr)')
      end if
      do index = 1, nf
c     
c     find wavelength limits of filter passband (measurements in microns)
c     
         wl0 = filtpars(index, 1)
         wlf = filtpars(index, 2)
         dwl = (wlf - wl0)/(lll)
         wl_nominal = filtpars(index,10)
         npts = nint((wlf - wl0)/dwl)
c
         pixel = 0.0317d0
         if(wl_nominal .gt. 2.2) pixel = 0.0648d0
         pixel_sr = (pixel/arc_sec_per_radian)**2
c
c     calculate the background sed within filter boundaries
c     
         f_lambda     = 0.0d0
         f_nu         = 0.0d0
         photon_flux  = 0.0d0
         do k = 1, lll
            factor = 1.0d0
            if(k.eq.1 .or.k.eq.npts) factor = 0.5d0
            wavelength =  filters(1, k, index)
            throughput =  filters(2, k, index) 
            wl_cm = wavelength * 1.d-04 ! um --> cm
c     
c     Zodiacal light : scattering term
c     
            exp1 = hc_over_k / (wl_cm * t_scatter)
            denom1 = dexp(exp1) - 1.d0
            term1_cgs = e_scatter * two_hcc /(denom1*wl_cm**5)
            term1_cgs =  term1_cgs * zodiacal_scale_factor 
c     
            term1_si  = term1_cgs * to_si_per_micron
c     
c     Zodiacal light : thermal emission
c     
            exp2 = hc_over_k / (wl_cm * t_thermal)
            denom2 = dexp(exp2) - 1.d0
            term2_cgs = tau_thermal * two_hcc /(denom2*wl_cm**5)
            term2_cgs =  term2_cgs * zodiacal_scale_factor 
c     
            term2_si  = term2_cgs * to_si_per_micron
c     
c     scattered background : f_nu --> f_lambda
c
            sb = 1.3d-19 * dwl_cm * cee/(wl_cm**2)
            sb = sb /21.d0
c     
c     total background
c     
            bkg    = term1_cgs + term2_cgs + sb
c     expressed in photon flux
            e_photon         = h * cee / wl_cm
            photons          = bkg/e_photon
            photon_flux = photon_flux + photons*throughput*factor
c     background incident on the telescope (prior to optics)
            f_nu     = f_nu     + (bkg * wl_cm**2) / cee 
            f_lambda = f_lambda + bkg * to_si_per_micron
c           print *, wavelength, e_photon, photons, bkg
c           PRINT *, k,wavelength,wl_cm, denom1, denom2, t_thermal
         end do
c
c     this will be the total number of photo-electrons integrated
c     over the filter band pass per pixel, taking into account 
c     the instrumental efficiency and telescope area
c
         electron_flux = photon_flux * dwl * 1.d-04 
         eflux_total   = electron_flux * mirror_area
         eflux_total   = eflux_total/arc_sec_per_radian**2
         events  = electron_flux * pixel_sr * mirror_area
         per_aperture =  events * pi * aperture**2
         results(index) = events
         temp = filterid(index)
c
c     Flux incident on the mirror
c
         f_lambda = f_lambda*dwl/(wlf-wl0) !* 1.d-04
         f_nu = f_nu*(dwl/(wlf-wl0))*1.d17! erg/(sec cm**2 hz sr) --> MJy/sr
         if(verbose.gt.0) then
            print 20, temp(1:6), wl_nominal, events , per_aperture,
     *           eflux_total,f_nu, f_lambda,electron_flux
 20         format(2x,a6,2x,f10.6,3x,f10.6,5x,f10.3,9x,f10.3,
     &           3x,4(2x,e13.5))
         end if
      end do
      return
      end
