c
c-----------------------------------------------------------------------
c
c     Rebin and attenuate spectrum to a redshift of z
c
c     Use expression (6) of  Hogg et al. 2002, ArXiV 0210394 :
c     Luminosity(emitted) = 4 * pi * r**2 * (1+z) * f_lam(obs)
c
c     If one assumes model spectra are "absolute", i.e., measured at
c     10 pc, the attenuation can be approximated as
c
c     attenuation = 10.**{-0.4*[distance_modulus(z) + 2.5d0*dlog10(1+z)]}
c     
      subroutine redshift_spectrum(z_0, z,wl, flux, nwl, wlz, fluxz,
     *     nwlz, debug)
      implicit double precision (a-h,o-z)
      integer debug
      parameter (nnn=25000)
      dimension wl(nnn), flux(nnn), wlz(nnn), fluxz(nnn)
c
c      pi = dacos(-1.0d0)
c      r_mpc        = 3.0856775807d24 ! cm
c      parsec       = 3.0856775807d18 ! cm
c
      zplus1 = 1.d0 + z
c
      corr  = distance_modulus(z) + 2.5d0 *dlog10(zplus1)
      attenuation  = 10.d0**(-0.4d0*corr)
c
      nwlz = 0 
      do i =1, nwl
         wlz(i)   = wl(i) * zplus1
         fluxz(i) = flux(i) * attenuation 
         nwlz     = nwlz + 1
      end do
      if(debug .eq.1) then
         print 100,z_0,z,attenuation, -2.5d0*dlog10(attenuation)
 100  format('redshift_spectrum :z_0, z, attenuation ',
     *        10(1x,e15.6))
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
c
c     This is used for cloning:
c
c     Rebin and attenuate spectrum to a redshift of z
c
c     Use expression (6) of  Hogg et al. 2002, ArXiV 0210394 :
c     f_lam (emitted) = 4 * pi * r**2 * (1+z) * f_lam(obs)
c     
      subroutine redshift_spectrum_h(z_0,z,wl, flux, nwl, wlz, fluxz,
     *     nwlz, debug)
      implicit double precision (a-h,o-z)
      double precision mpc_in_cm
      integer debug
      parameter (nnn=25000)
      dimension wl(nnn), flux(nnn), wlz(nnn), fluxz(nnn)
c
      pi = dacos(-1.0d0)
      mpc_in_cm    = 3.0856775807d24 ! cm
c      parsec       = 3.08567802d18 ! cm
c      ten_pc       = parsec *10.d0
c
      zplus1        = 1.0 + z
      if(z.eq.0.0d0) then
         r_mpc = 1.d-05             ! set at 10 pc
      else
         r_mpc = dl(z)
      end if
      r_cm  = r_mpc * mpc_in_cm
c
      attenuation = 4.d0 * pi * r_cm * r_cm * zplus1
      attenuation = 1.d0/attenuation
      print *, 'attenuation',attenuation
c
      nwlz = 0 
      do i =1, nwl
         wlz(i)   = wl(i) * zplus1
         fluxz(i) = flux(i) * attenuation 
         amag0    = -2.5d0*dlog10(flux(i))
         amag1    = -2.5d0*dlog10(fluxz(i))
         if(mod(i,100).eq. 1) print 90, wl(i), flux(i),fluxz(i), 
     *        amag0, amag1
 90      format(f9.3,2(1x,e15.3), 2(1x,f8.3))
         nwlz     = nwlz + 1
      end do
      if(debug .eq.1) then
         print 100,z,attenuation,-2.5d0*dlog10(attenuation)!,bandwidth
 100  format('redshift_spectrum :z, attenuation, att. in mag',
     *        10(1x,f9.5))
      end if
      return
      end
c      z_ratio       = ((1.d0+z_0)/zplus1)**3
c      attenuation   = z_ratio   ! * dlambda_ratio 
c      if(z.eq.0.0d0) then
c         r = 1.d-05             ! set at 10 pc
c      else
c         r             = dl(z)
c      end if
c     commenting this gives consistent results with Bouwens et al.
c     otherwise the galaxies become exceedingly luminous
c         fluxz(i) = flux(i) * attenuation 
c      bandwidth   = -2.5d0*dlog10((1.d0+z_0)/zplus1)
