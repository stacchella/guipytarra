c
c-----------------------------------------------------------------------
c
c     Calculate the mean energy flux through a photon-counting system
c     following  Bessell & Murphy 2012, PASP, 124, 140
c     This combines expressions (A12) for average f_nu, (A29) for
c     AB magnitude
c
c     Input spectrum must be in units of erg/sec/cm**2/A
c     Filters must already be rebinned
c     output average flux will be in units of
c     [f_nu] = [erg/( sec Hz cm**2)]
c     [f_lambda] = [erg/(sec A cm**2)]
c
c     2015-06-15
c     cnaw@as.arizona.edu
c
      subroutine sedflux(indx, z, wlsed, sed, nsed, f_lambda, f_nu)
      implicit none
      double precision wlsed, sed, sed_at_wl, effective_wl,h,
     *     wl0, dwl, resp, wavelength, numerator, denominator,cee,
     *     hplanck, f_lambda, f_nu, wl_a, ab_mag, zmag, total_flux,
     *     z, zplus1, rest_wl
      double precision previous, wldif, dl
      integer indx, nsed, iii, i, k, mmm, nnn, lll, npar
      double precision filters, filtpars
      character filterid*20
c
      parameter (nnn=25000,mmm=200,lll=1000,npar=30)
c
      dimension wlsed(nnn), sed(nnn)
      dimension filters(mmm,nnn),filtpars(mmm,npar), filterid(mmm)
      common /filter/filters, filtpars, filterid
      cee     = 2.99792458d18       !  c in A/s
      hplanck = 6.626068963d-27 ! erg s 
      wl0     = filtpars(indx,1)
      dwl     = filtpars(indx,11)
      numerator   = 0.0d0
      f_nu        = 0.0d0
      f_lambda    = 0.0d0
      zplus1      = z + 1.0d0
      previous    = 0.0d0
      do k = 1, lll
         h = 1.0d0
         if(k.eq.1 .or.k.eq.lll) h = 0.5d0
         wavelength=wl0 + dble(k-1) * dwl
         rest_wl   = wavelength/zplus1
         if(rest_wl  .ge. wlsed(1) 
     *        .and. rest_wl .le. wlsed(nsed)) then
            call linear_interpolation(nsed, wlsed, 
     *           sed, rest_wl, sed_at_wl, nnn)

            if(previous .ne.0.d0) then
               wldif    = rest_wl - previous
               previous = rest_wl
            else
               wldif    = dwl
               previous = rest_wl
            end if
c
            resp        = filters(indx, k) 
! convert wavelength from microns to Angstroms
            wl_a        = rest_wl * 1.d04
            numerator   = numerator + h*resp*sed_at_wl*wl_a
            f_nu        = f_nu     +  h * resp / wl_a
            f_lambda    = f_lambda +  h * resp * wl_a
c            if(mod(k,500).eq.0) print *,k, wl_a, f_nu
         end if
      end do
      if(f_nu.ne.0.0d0) then
         f_nu     = numerator/(cee*f_nu)
      end if
      if(f_lambda .ne.0.0d0) then
         f_lambda = numerator/f_lambda
      end if
c     
c     attenuate because of distance 
c     f_nu  =  f_nu / ( 1.d0*dl(z)**2) where 1.d10 = [1 Mpc/10pc]**2
c     f_lam = f_lam / ( 1.d0*dl(z)**2)
c
c     attenuate because of wavelength stretching
c
c     f_nu  =  f_nu * (1+z)
c     f_lam = f_lam / (1+z)
c
      if(z.gt.0.0d0) then
         f_nu     = f_nu   * zplus1 /(1.d10*(dl(z))**2)
         f_lambda = f_lambda/(zplus1* 1.d10*(dl(z))**2)
      end if
c     
      return
      end
c
c----------------------------------------------------------------------
c
c     Integrate flux through a filter
c
c     2014 04 30: While using wl_nominal (Reach) gives the best match to
c     results reported by Rieke et al. (2008) for 2MASS, IRAC and MIPS
c     filters, when calculating effective or average fluxes for a flat
c     F_NU input SED, the resulting average value differs from the
c     input spectrum. The only way the results will match is using
c     the Fukugita et al. procedure (for CCDFLUX) and the effective-
c     wavelength. 
c     2015-02-13 This had been noted by Bessell & Murphy 2012, PASP,124,140
c     appendix a.2.1, where they define this as the pivot wavelength
c
c      subroutine sedflux(index, z, wlsed, sed, nsed, f_lambda, f_nu)
c      implicit double precision (a-h,o-z)
c      double precision mag_vega, mag_ab, nominal
c      character filterid*20
c      parameter (nnn=25000,mmm=200,lll=1000,npar=30)
c      dimension wlsed(nnn), sed(nnn), deriv(nnn)
c      dimension filters(mmm,nnn),filtpars(mmm,npar), filterid(mmm)
c      common /filter/filters, filtpars, filterid
cc
c      cee = 2.99792458d18       !  c in A/s
cc
cc     Spline interpolate SED
cc     units of wavelength are microns = 10**4 A
cc     units of flux are (f_lambda) = erg/(s cm**2 A)
cc                       (f_nu)     = erg/(s cm**-2 Hz) 
cc
c      call spline(wlsed, sed, nsed,0.0d0,0.0d0, deriv)
cc
c      zplus1  = z + 1.0d0
cc     
c      wl0     = filtpars(index,1)
c      dwl     = filtpars(index,11)
c      offset  = filtpars(index,13)
c      nominal = filtpars(index,9)
c      
c      filter     = 0.0d0
c      ccdflux    = 0.0d0
c      ccdresp    = 0.0d0
c      rieke_resp = 0.0d0
c      effective  = 0.0d0
c      previous   = 0.0d0
c      flux       = 0.0d0
c      do k = 1, lll
c         h = 1.0d0
c         if(k.eq.1 .or.k.eq.lll) h = 0.5d0
c         wavelength=wl0 + dble(k-1) * dwl
c         rest_wl   = wavelength/zplus1
c         if(rest_wl .lt. wlsed(1) 
c     *        .or. rest_wl .gt. wlsed(nsed)) go to 1000
cc
cc     Find object flux at rest wavelength corresponding to this
cc     observed wavelength and redshift
cc
c         call splint(wlsed, sed, deriv, nsed,rest_wl, sed_at_rest_wl)
cc     
cc     Filter response
cc
c         throughput  = filters(index, k) 
cc
cc     Because of the sampling rate, keep track of the last 
cc     wavelength that was used, otherwise the total flux will
cc     be made too large
cc
c         if(previous .ne.0.d0) then
c            wldif    = rest_wl - previous
c            previous = rest_wl
c         else
c            wldif    = dwl
c            previous = rest_wl
c         end if
cc         
cc     These come from Fukugita et al. 1995 PASP 107, 945 in the
cc     case of a filter convolved with the detector QE, and doing the
cc     following conversions:
cc     numerator (ccdflux)
cc     dnu * f_nu * s(nu)/nu = dlambda * f_lambda * s(lambda) * lambda/cee
cc     and denominator (ccdresp)
cc     dnu * s(nu)/ nu  = dlambda * (cee/lambda**2) * (lambda/cee) *s(lambda)
c 
c         ccdflux      = ccdflux + 
c     *        throughput * sed_at_rest_wl * rest_wl * h * wldif
c         ccdresp   = ccdresp   + (throughput / rest_wl) * h * wldif
c         filter    = filter    + throughput * h 
cc     using rest wavelengths
c         rieke_resp = rieke_resp + rest_wl   * throughput * wldif * h
cc     using observed wavelengths
c         effective  = effective + wavelength * throughput * h
c         flux       = flux + sed_at_rest_wl  * throughput * wldif * h
c      end do
cc
cc
cc     Since ccdresp uses rest wavelengths, wl_nominal will not
cc     the same as nominal
cc      
c      filter      = filter * dwl
c      wl_nominal  = filter/ccdresp
cc
cc     effective flux
cc
c      effective        = effective * dwl
c      wl_effective     = effective/filter
c      f_effective      = flux/filter
cc      print *,' wl_nominal ', wl_nominal
cc       print *,'sedflux: effective      ', wl_effective f_effective
cc
cc     Flux integrated using CCDs (Fukugita et al. 1995); 
cc     units should be erg/(sec cm2 micron)
cc
cc      ccdflux   = ccdflux/cee 
c      pivot     = dsqrt(effective/ccdresp)
c      f_lambda  = ccdflux/ccdresp
c      fis_flambda  = ccdflux/filter/wl_nominal
c      print *,'sedflux FIS               ', f_lambda, fis_flambda,pivot
cc
cc     f_lambda according to Rieke et al. 2008
cc
c      f_lambda   = ccdflux / rieke_resp
cc      print *,'sedflux Rieke et al. 2008 ', f_lambda
cc
cc     f_nu is in units of  erg/(s cm**-2 Hz) 
cc     the 1.d-04 factor is due to 10**4 (A/micron)/10**8(A/cm)
cc
cc      print *,'sedflux : different f_nu estimators'
cc      f_nu       = fis_flambda *  (wl_effective*1.d4)**2/cee
cc      print *,'sedflux : f_nu 1', fis_flambda, wl_effective, f_nu
cc     This is the expression that reproduces a flat f_nu spectrum
c      f_nu       = fis_flambda *  (wl_nominal*1.d4)**2/cee
cc      print *,'sedflux : f_nu 2 ', fis_flambda,wl_nominal,  f_nu
cc
cc      f_nu       = f_effective *  (wl_effective*1.d4)**2/cee
cc      print *,'sedflux : f_nu 3', f_effective, wl_effective, f_nu
cc      f_nu       = f_effective *  (wl_nominal*1.d4)**2/cee
cc      print *,'sedflux : f_nu 4',f_effective, wl_nominal, f_nu
ccc
cc      f_nu       = f_lambda    *  (wl_effective*1.d04)**2/cee
cc      print *,'sedflux : f_nu 5', f_lambda, wl_effective, f_nu
cc      f_nu       = f_lambda    *  (wl_nominal*1.d4)**2/cee
cc      print *,'sedflux : f_nu 6', f_lambda , wl_nominal, f_nu
c      f_nu_jy    = f_nu * 1.0d23
c      print 200, filterid(index), z, nominal, wl_nominal,f_lambda,
c     *     f_nu_jy
c 200  format(a20,1x,3(1x,f8.4),2(1x,1pe13.6),' jy sedflux')
c      return
c 1000 f_lambda  =  1.d-200
c      f_nu      =  1.d-200
c      print *, 'sedflux: wavelength range of spectrum too small'
c      print *, 'sedflux:filter ',index, wl0,wl0+lll*dwl
c      print *, 'sedflux:sed    ', wlsed(1), wlsed(nsed)
c      return
c      end

