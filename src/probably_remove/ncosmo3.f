c
c***********************************************************************
c
      subroutine is_cosmology_defined
c
c     Verify if cosmology parameters are set and select mode to calculate:
c     i.e., if reproducing the plots (and values) obtained using 
c     Hogg, D. W. 1999, astro-ph/990511v4 from 2000-12-16 (mode 1)
c     or from Ned Wright's Cosmology calculator (modes 2 and 3), described
c     in  Wright, E. L. 2006, PASP, 118, 1711
c
c     options
c     1. lambda  = 0, no neutrinos       (e.g., Hogg 1999)
c     2. lambda != 0, no neutrinos       (e.g., Hogg 1999)
c     3. massless neutrinos (original calculator)
c     4. neutrino masses + equation of state
c
c     In option 1 omega_r = 0,  omega_l = 0 omega_cm = omega_m 
c     In option 2 omega_r = 0               omega_cm = omega_m 
c     In option 3 omega_r = 4.165d-05/h**2, omega_cm = omega_m
c     in option 4 omega_r = 2.473d-05/h**2, omega_m = omega_nu + omega_cm
c     where omega_m = omega_cm + omega_nu
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-03-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      integer mode
      double precision tnow, t0, n_eff, eps, h
      double precision omega_neutrino
      double precision  h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
      parameter(eps = 1.d-16)
      common /cosmology/ h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
c
      if(h0 .le. 0.0d0 .or. omega_m .le.0) then
         print *, '  '
         print *, 'need to define cosmology parameters !'
         print *, 'h0, omega_m, omega_l', h0,omega_m, omega_l
         print *, '  '
         stop
      end if
c
      h = h0/100.d0
c
c     from Planck collaboration 2015, Astro-ph/1502.01589v2
c
      n_eff = 3.15              
c
      if(omega_r .eq.0.d0 .and. omega_nu .eq.0.0d0) mode = 1
      if(omega_r .ne.0.d0 .and. omega_nu .eq.0.0d0) mode = 2
      if(omega_r .ne.0.d0 .and. omega_nu .ne.0.0d0) mode = 3
c
c     Uncomment if you really want to use the values of Ned Wright's calculator
c      if(mode.eq.3) mode = 4
c
c     mode 1 reproduces the Hogg (2000) plots for Einstein-de Sitter,
c     omega_m = 0.05, omega_l = 0.0 and omega_l=0.2, omega_l = 0.8 
c     models; mode = 4 uses Ned Wright's values of his Advanced calculator
c
c     
      if(mode.eq.1) then
         omega_cm = omega_m
         w        = -1.0d0
         wprime   =  0.0d0
         omega_k  = 1 - omega_m - omega_l
         if(dabs(omega_k).lt.eps) omega_k = 0.0d0
         return
      end if
c
c     mode 2 reproduces E. L. Wright's (2006) default calculator
c     http://www.astro.ucla.edu/~wright/CC.html
c     using 3 massless neutrinos and T0 = 2.72528 K in addition to
c     values of omega_m and omega_l
c
      if(mode.eq.2) then
         omega_cm = omega_m
c     omega_nu =  0.0d0
c     omega_r  = 4.165d-05/(h*h)
c     w        = -1.0d0
c     wprime   =  0.0d0
         omega_k  = 1 - omega_m - omega_l - omega_r
         if(dabs(omega_k).lt.eps) omega_k = 0.0d0
         return
      end if
c
c     mode 3 follows E. L. Wright's (2006) advanced cosmology calculator
c     http://www.astro.ucla.edu/~wright/ACC.html
c     which includes the equation of energy and neutrino masses
c     The paper quotes a value of 2.778d-05 for the CMB radiation density, while
c     the calculator uses 2.477d-05. The calculator uses Tnow = 2.72528d0
c     The value for the CMB radiation density quoted in the 
c     Review of Particle Physics 2014 edition 
c     (Olive, K. A. et al. 2014 Chin. Phys. C, 38, 090001) is
c     Omega_r = 2.473e-5 (T/2.7255)**4 h**−2
c     The most recent value of Tnow = 2.7255d0 +/- 0.0006 K by 
c     Fixsen 2009, ApJ, 707, 916
c     
      if(mode.eq.3) then
         omega_cm = omega_m - omega_nu
         tnow     = 2.7255d0
         t0       = tnow
         omega_r  = 2.473d-05/(h*h)*(t0/tnow)**4 
         omega_k  = 1 - omega_m - omega_l - omega_r
      end if
c
c     force to use Ned Wright's values. However the differences are 
c     insignificant relative to the Fixsen + Review of Particle Physics values:
c     ~ 0.02% for the co-moving volume and 2.2d-05% for the luminosity distance
c     To make this run you have to uncomment the if(mode.eq.3) line above
c
      if(mode.eq.4) then
         omega_cm = omega_m - omega_nu
         tnow     = 2.72528d0 ! value used by Wright (2006)
         t0       = tnow
         omega_r  = 2.477d-05/(h*h)*(t0/tnow)**4 ! value used by Wright (2006)
         omega_k  = 1 - omega_m - omega_l - omega_r
      end if
      return
      end
c
c***********************************************************************
c
      double precision function omega_neutrino(z)
c
c     Calculate the contribution to the matter density due to the
c     Cosmic Neutrino Background. This is a translation of
c     E. L. Wright's Advanced Cosmology Calculator at
c     http://www.astro.ucla.edu/~wright/ACC.html
c     and described in  Wright, E. L. 2006, PASP, 118, 1711
c
c     Mangano et al. (2005) :
c     omega_nu = rho_nu/rho_critical = 3 * m_nu/(93.14*h**2 eV) (eq. 18)
c     contribution of neutrinos to present energy density is n_eff:
c     rho_rad = [ 1 + (7/8) * (4/11)**(4/3) * n_eff] * rho_photons (eq. 1)
c     where rho_photons is the present energy density of photons coming
c     from CMB temperature
c
      implicit none
      double precision m_e, m_mu, m_tau, cee, h, m_rel, sum, beta, a, 
     *     tnow, t0, m_nu, nu_temp, momentum, result, nurho, one_ev, z
      integer i
      dimension m_nu(3)
      double precision  h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
      common /cosmology/ h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
c
      h           = h0/100.d0
      cee         = 299792.458d0
      m_e         = 0.001d0     ! Neutrino mass in eV/cee**2
      m_mu        = 0.009d0     ! eV/cee**2
      m_tau       = 0.049d0     ! ev/cee**2
      tnow        = 2.72528d0   ! K
      t0          = tnow        ! should this be a variable ?
      one_ev      = 11604.505d0
c
      m_nu(1)     = m_e
      m_nu(2)     = m_mu
      m_nu(3)     = m_tau
c
c     The  present day neutrino temperature relative to the photon 
c     temperature given by the CMB (Wright 2006, Section 5):
c
      nu_temp = t0 * (4.d0/11.d0)**0.333333333333d0
c
c     Mass of a neutrino "that is just now relativistic"
c
      momentum = 3.151d0 * nu_temp 
      m_rel    = (momentum/one_ev) * (t0/tnow)
c     The values printed below match those on paper (Section 5)
c      print *,' nu_temp, momentum, m_rel',nu_temp, momentum, m_rel
c
c     While the paper quotes a value of 93.14 eV for the normalisation (eq.18),
c     the ACC itself uses different values for the electron and mu and  tau
c     neutrinos (93.64, 93.90, 93.90); 
c     Hannestad & Madsen 1995,  quote (92.55, 92.80, 92.80), for T=2.736K
c     Mangano et al. (2005) hep-ph/0506164 use 93.14; a normalisation value
c     of 93.04 eV for the three species is used in the 
c     Review of Particle Physics 2014 edition:
c     Olive, K. A. et al. 2014 Chin. Phys. C, 38, 090001
c
      a    = 1.d0/(1.d0+z)
      sum  = 0.0d0
      omega_neutrino = 0.0d0
      do i = 1, 3
         result = nurho(m_rel, m_nu(i)*a)
c         if(i.eq.1) then
c            omega_neutrino = omega_neutrino + result * m_nu(i)/93.64d0
c         else
c            omega_neutrino = omega_neutrino + result * m_nu(i)/93.90d0
c         end if
         omega_neutrino = omega_neutrino + result * m_nu(i)/93.04d0
      end do
      omega_neutrino = omega_neutrino *(t0/tnow)**3
      omega_neutrino = omega_neutrino/(h*h)
      return
      end
c
c-----------------------------------------------------------------------
c
c     Calculate the neutrino density over rest mass density
c     Translation of N. Wright's AAC calculator and 
c     Wright, E. L. 2006, PASP, 118, 1711, eqs. 16 (and 18).
c
      double precision function nurho(mnurel, mnu)
      double precision y, beta, inv_beta,mnurel, mnu
      beta     = 1.842d0
      inv_beta = 1.d0/beta
      y        = 1.d0 + (mnurel/mnu)**beta
      nurho    = y**inv_beta
      return
      end
c
c
c-----------------------------------------------------------------------
c
c     Calculate the contribution to Omega due to dark energy
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c     
      double precision function rho_de(z)
      implicit none
      double precision a, expo1, expo2, z
      double precision  h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
      common /cosmology/ h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
c
      a       = 1.d0/(1.d0+z)
c
      expo1   = 3.d0+3.d0*w+ 6.d0*wprime
      expo1   = -expo1
      expo2   = 6.d0*wprime*(a -1.0d0)
      rho_de  = omega_l * a**expo1 * dexp(expo2)
      return
      end
c
c----------------------------------------------------------------------
c
      double precision function th_in_gyr(th_hubble)
c
c     Convert Hubble times into Gyr
c     constants from 
c     http://pdg.lbl.gov/2014/reviews/rpp2014-rev-astrophysical-constants.pdf
c     http://pdg.lbl.gov/2014/reviews/contents_sports.html
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision r_mpc, year,th_hubble 
      double precision  h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
      common /cosmology/ h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
      r_mpc     = 3.08567758149d19 ! km
      year      =  31556925.2 ! tropical year in seconds for 2011
      th_in_gyr = r_mpc/h0/year/1.d9 * th_hubble
      return
      end

c
c----------------------------------------------------------------------
c
      double precision function ez(z)
c
c     Hogg (2000) eq. 14 for the function E(z) which is 
c     "proportional to the time derivative of the logarithm of the scale factor"
c     and used to calculate the line of sight comoving distance (eq. 15)
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision h, eps, omega_tot, omega_nu_z, z
      double precision omega_neutrino, rho_de, omega_de
      parameter (eps=1.0d-12)
      double precision  h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
      common /cosmology/ h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
c
c      if(h0.eq.0.0d0) then
c         print *,'h0 = ',h0
c         print *,'need to define cosmology parameters !'
c         stop
c      end if
c      call is_cosmology_defined
c
c     fall back onto the Hogg (2000) or Wright (2006) CC.html calculation
c     for omega_nu = 0.0. Otherwise use the Wright (2006) advanced calculator
c     expression where omega_k is a function of redshift that takes into account
c     the changing omega_nu and dark energy values
c
      if(omega_nu.eq.0.0d0) then
         ez = omega_m * (1.d0+z)**3 + omega_k * (1.d0+z)**2 + omega_l +
     *        omega_r * (1.d0+z)**4
      else
         omega_nu_z = omega_neutrino(z) 
         omega_de   = rho_de(z)
         ez  = (omega_cm+omega_nu_z)*(1.d0+z)**3 + omega_r * (1.d0+z)**4
     *        + omega_de + omega_k * (1.d0+z)**2
      end if
      ez = dsqrt(ez)
c
      return
      end
c
c-----------------------------------------------------------------------
c
      double precision function ezinverse(z)
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision ez, z
      external ez
      ezinverse = 1.d0/ez(z)
c      print *, 'ezinverse ', z, ezinverse
      return
      end
c
c-----------------------------------------------------------------------
c
c     Line of sight comoving distance (Hogg eqs. 15 and 17)
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      double precision function dc(z)
      implicit none
      double precision ezinverse, cee, eps, dh, term1, term2,
     *     term3, zmin, zmax, result,HRVINT, fac, z, yy, dm_over_dh
      integer mfinal, m1
      parameter (cee=299792.458d0, eps =1.0d-16)
      double precision  h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
      common /cosmology/ h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
c
      external ezinverse
      m1     = 16 ! number of iterations for numerical integration
c
      if(h0.eq.0.0d0) then
         print *,'h0 = ',h0
         print *,'need to define cosmology parameters !'
         stop
      end if
c
      call is_cosmology_defined
c
      dh = cee/h0
c
      if(z.eq.0.d0) then
         dc = 0.d0
         return
      else
c
c     if omega_lambda = 0, there are  analytical solutions for D_m.
c     for omega_k = 0 , dc = dm;
c     for omega_k > 0 a solution for D_c is derived from  Eqs. (16) and (17):
c     arcsinh[(dm/dh) * sqrt(omega_k)] = sqrt(omega_k) *dc/dh
c     where
c     arcsinh(z) = ln(z+sqrt(1+z**2))
c
         if(omega_l .eq. 0.d0) then ! Eq. 17
            term1 =  2.d0 - omega_m * (1.d0 - z)
            term2 = (2.d0 - omega_m) * dsqrt(1.d0 + omega_m*z)
            term3 = (1.d0+z) * omega_m**2
            dm_over_dh = 2.d0 * (term1 - term2)/ term3
            if(omega_k.eq.0.0d0) then
               dc = dm_over_dh*dh
            else
               yy = dm_over_dh * dsqrt(omega_k)
               dc = dh * dlog(yy+dsqrt(1.d0+yy**2))/dsqrt(omega_k)
            end if
         else                   ! Eq. 15
            zmin = 0.0d0
            zmax = z
            result = hrvint(ezinverse,zmin,zmax,m1, eps, fac, mfinal)
c            print *, 'dc(z): zmin, zmax, result, fac',
c     8           zmin, zmax, result, fac
c            print *, 'dc(z)', h0, omega_m, omega_l, omega_k, omega_r, 
c     *     omega_nu, omega_cm, w, wprime
            dc = dh * result
         end if
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      double precision function dm(z)
c
c     Calculate the transverse comoving distance (eq. 16)
c     which is equivalent to the proper motion distance. This uses
c     eqs 15 or 17 (used in dc(z) above).
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision cee, eps, dc, z, x, dh, term1, term2, term3
      parameter (cee=299792.458d0, eps=1.0d-12)
      double precision  h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
      common /cosmology/ h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
c
      external dc
c
      dh = cee/h0
      if(omega_l .eq. 0.0d0) then
         term1 =  2.d0 - omega_m * (1.d0 - z)
         term2 = (2.d0 - omega_m) * dsqrt(1.d0 + omega_m*z)
         term3 = (1.d0+z) * omega_m**2
         dm    = dh * 2.d0 * (term1 - term2)/ term3
         return
      end if
c
      if(omega_k .eq. 0.0d0) then
         dm = dc(z)
         return
      end if
c
      if(omega_k .lt. 0.d0) then
         x  = dsqrt(dabs(omega_k)) * dc(z) / dh
         dm = dh * dsin(x)/dsqrt(dabs(omega_k))
         return
      end if
c  
      if(omega_k .gt. 0.d0) then
         x = dsqrt(omega_k) * dc(z) / dh
         dm = (dexp(-x) - dexp(x))/2.d0
         dm = dabs(dm)
         dm = dh * dm /dsqrt(omega_k)
         return
      end if
c
      return
      end
c
c-----------------------------------------------------------------------
c
      double precision function da(z)
c
c     Angular diameter distance Hogg (Eq. 18)
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision z, dm
      external dm
      da = dm(z)/(1.d0 + z)
      return
      end
c
c-----------------------------------------------------------------------
c
c
      double precision function dl(z)
c
c     Luminosity distance Hogg (eq. 21)
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision z, dm
      external dm
      dl = dm(z) * (1.0d0 + z)
      return 
      end

c-----------------------------------------------------------------------
c
      double precision function tlookback(z)
c
c     Lookback time Hogg (eq. 30), referencing 
c     Peebles, 1993, Principles of Physical Cosmology, Princeton Univ. Press,
c     pages 313-315
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      integer m1, mfinal
      double precision z, zplus1ezi, th, zmin, zmax, result, 
     *     hrvint, eps, fac
      double precision  h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
      common /cosmology/ h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
c     
      external zplus1ezi
c     
      call is_cosmology_defined
c
      eps = 1.0d-16 
      m1  = 16
c     
      th = 1.d0/h0
      if(z.eq.0.d0) then
         tlookback = 0.d0
         return
      else
         zmin = 0.0d0
         zmax = z
         tlookback = hrvint(zplus1ezi,zmin,zmax,m1, eps, fac, mfinal)
      end if
      tlookback = tlookback * th
      return
      end
c
c-----------------------------------------------------------------------
c
      double precision function zplus1ezi(z)
c
c     function used to determine the lookback time
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision z, ez
      external ez
      zplus1ezi = 1.d0/(ez(z)*(1.d0+z))
      return
      end
c-----------------------------------------------------------------------
c
c     
      double precision function tform(z)
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      integer m1, mfinal
      double precision z, th, zmax, zmin, result, hrvint, zplus1ezi,
     *     eps, fac
      double precision  h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
      common /cosmology/ h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
      external zplus1ezi
c
      call is_cosmology_defined
c
      m1 = 16
      eps = 1.0d-16
c
      th = 1.d0/h0
      zmax = 3000.0
      if(z.ge.zmax) then
         tform = 0.d0
         return
      else
         zmin = z
         zmax = zmax
         tform = hrvint(zplus1ezi,zmin,zmax,m1, eps, fac, mfinal)
      end if
      tform = tform * th
      return
      end
c
c-----------------------------------------------------------------------
c
      double precision function vol_elm(z,solid_angle)
c
c     Comoving Volume element Hogg (2000) eq. 28, referenced to
c     Weinberg 1972, Gravitation and Cosmolgy: Principles and Applications
c     of the General Theory of Relativity, John Wiley & Sons, New York, pg 486
c     and 
c     Peebles, 1993, Principles of Physical Cosmology, Princeton Univ. Press,
c     pages 331-333
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision z, solid_angle, cee, ez, da, dh, daz
      double precision  h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
      parameter (cee=299792.458d0)
      common /cosmology/ h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
      external ez, da
c     
      dh = cee/h0
      if(z.le.0.0d0) then
         vol_elm = 0.d0
         return
      end if
      daz = da(z)
      vol_elm = dh * (1.d0 + z) **2 * daz * daz * solid_angle/ez(z)
      return
      end
c
c-----------------------------------------------------------------------
c
      double precision function arcsinh(z)
      implicit none
      double precision z
      arcsinh = dlog(z+ dsqrt(1.d0+z**2))
      return
      end
c
c
c-----------------------------------------------------------------------
c
      double precision function zvolumen(solid_angle,zmin,zmax)
c
c     Find the volume in Mpc**3 between zmin and zmax
c     using Hogg (2000) eq (29), referencing 
c     Carroll, Press and Turner 1992, ARAA, 30, 311
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision solid_angle, zmin, zmax, cee, eps, dm,
     *     dh, dvol0, rm, rm_over_dh, term1, term3, term2, arcsinh,z
      double precision  h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
c
      parameter (cee=299792.458d0, eps=1.0d-12)
c
      common /cosmology/ h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
c
      if(zmax .eq. 0.0d0 .or. zmin.eq.zmax) then
         zvolumen = 0.0d0
         return
      end if
c
      dh = cee/h0
c
      if(zmin.eq.0.0d0. or.zmin.lt.eps) then 
         dvol0 = 0.0d0
      else
c     
c     the call to DM will update rho_de and omega_k in the case of
c     a general cosmology through calls to DC and EZ
c
         rm  = dm(zmin)
c     
         if(omega_k .eq. 0.0d0) then
            dvol0 = (solid_angle/3.0d0) * rm**3
         else
            rm_over_dh = rm/dh
            term1 = rm_over_dh*dsqrt(1.0d0 + omega_k * rm_over_dh**2)
            term3 = rm_over_dh*dsqrt(dabs(omega_k))
c     
            if(omega_k .gt. 0.0d0) then
               term2 = (1.d0/dsqrt(dabs(omega_k))) * arcsinh(term3)
            else
               term2 = (1.d0/dsqrt(dabs(omega_k))) * dasin(term3)
            end if
c     
            dvol0 = (solid_angle*dh**3)/(2.0d0*omega_k) * (term1-term2)
         end if
      end if
c
      rm  = dm(zmax)
c
      if(omega_k .eq. 0.0d0) then
         zvolumen = (solid_angle/3.0d0) * rm**3 - dvol0
         return
      end if
c
      rm_over_dh = rm/dh
      term1 = rm_over_dh*dsqrt(1.0d0 + omega_k * rm_over_dh**2)
      term3 = rm_over_dh*dsqrt(dabs(omega_k))
c     
      if(omega_k .gt. 0.0d0) then
         term2 = (1.d0/dsqrt(dabs(omega_k))) * arcsinh(term3)
      else
         term2 = (1.d0/dsqrt(dabs(omega_k))) * dasin(term3)
      end if
c
      zvolumen = ((solid_angle*dh**3)/(2.0d0*omega_k)) *(term1-term2)
      zvolumen = zvolumen - dvol0
c      print 10, zmax, zvolumen, rm/dh,defint, term3, omega_k
 10   format(10(1pe13.5))
      return
      end
c-----------------------------------------------------------------------
c
      double precision function distance_modulus(z)
c
c     Distance modulus (Eq. 25) but in Mpc
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision z, dl
      external dl
      if (z.eq.0.0d0) then
         distance_modulus = 0.0d0
         return
      end if
      distance_modulus = 5.0d0 * dlog10(dl(z)) + 25.d0
      return
      end
c
c-----------------------------------------------------------------------
c
c     
      double precision function intersection_probability(z)
c
c     dimensionless probability of intersecting objects along the
c     line of sight Hogg (2000) eq. 31
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-20
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision z, ez, dh, cee
      double precision  h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
      parameter (cee=299792.458d0)
      common /cosmology/ h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
c
      call is_cosmology_defined
      dh = cee/h0
      intersection_probability = dh * (1.d0+z)*(1.d0+z)/ez(z)
      return
      end
c
c-----------------------------------------------------------------------
c
      double precision function zdistmodn(distmod)
c
c     Obtain redshift from the distance modulus in magnitudes
c     using the bisection method of
c     R. Brent, 1973 in "Algorithms for Minimization without Derivatives", 
c     Prentice-Hall.
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision zeroin, tol, eps, zmin, zmax, dmod, distmod,
     *     distance_modulus
      common /func_param/ dmod
      external distance_modulus, distmodn
c
      dmod = distmod
      zmin = 1.0d-03
      zmax = 20.0d0

      tol = 1.0d-16
      eps = 1.0d-20
c
      zdistmodn=zeroin(zmin,zmax,distmodn, tol,eps)
      return
      end
c-----------------------------------------------------------------------
c
      double precision function distmodn(z)
c
c     Auxiliary function for zdistmodn
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision z, dmod, distance_modulus
      common /func_param/ dmod
      distmodn = distance_modulus(z) - dmod
c      print 100, distmodn
 100  format(1pe17.10)
      return
      end
c
c
c-----------------------------------------------------------------------
c     Find absolute magnitude only considering the distance modulus
c
      double precision function getnmabs(mapp,z)
c
c     Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c     Based on original routine written by David Reiss, UCSC, 1992
c
      implicit none
      double precision mapp, z, distance_modulus
      getnmabs = mapp - distance_modulus(z)
      return
      end
c
c-----------------------------------------------------------------------
c     Find the apparent magnitude by adding the distance modulus
c     and the kcorrection calculated elsewhere
c
      double precision function getnmappar(absm,z,kcorval)
c
c     Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c     Based on original routine written by David Reiss, UCSC, 1992
c
      implicit none
      double precision absm, z, kcorval, distance_modulus
      getnmappar = absm + distance_modulus(z) + kcorval
      return
      end
c
c-----------------------------------------------------------------------
c     Find projected distance in Mpc for an angle in radians
c
      double precision function getnproj(angle,z)
c
c     Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c     Based on original routine written by David Reiss, UCSC, 1992
c
      implicit none
      double precision da, angle, z
      external da
      getnproj = angle * da(z)
      return
      end
c-----------------------------------------------------------------------
c     Find angle on the sky for a given separation in Mpc
c
      double precision function getnzangle(projdist,z)
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
      implicit none
      double precision da, projdist, z
      external da
      getnzangle = projdist / da(z)
      return
      end
c
c-----------------------------------------------------------------------
c
c     Calculate the volume for a given set of concentric shells as
c     well as the total volume.
c
      subroutine getnvol(zmin, zmax,zzz,nbins,solid_angle, vtot, dvol)
c
c     zmin, zmax  - limits in redshift 
c     zzz         - array with external redshift  of shell
c     nbins       - length of zzz
c     solid_angle - solid angle in steradians
c     vt          - total volume
c     dvol        - volume of each shell correspoding to zzz
c
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-27
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision zmin, zmax,zzz, solid_angle, vtot, dvol, pi,
     *     zvolumen
      integer max, i, nbins
      parameter(max = 1000)
      dimension dvol(max),zzz(max)
      if(nbins.gt.max) then
         print *, 'ncosmo : getnvol nbins > max ', nbins, max
      end if
c
      pi=dacos(-1.0d0)
c
c     Find volume contained by each shell
c     
      dvol(1) = zvolumen(solid_angle,zmin,zzz(1)) 
c
      do i=2,nbins
         dvol(i)= zvolumen(solid_angle, zzz(i-1),zzz(i))
      enddo
c
c     total volume
c
      vtot= zvolumen(solid_angle, zmin, zmax)
      return
      end
c
c-----------------------------------------------------------------------
c
c     Find z as a function of vel (km/s)
c
      double precision function getz(v)
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision cee, v
      cee = 299792.458d0
      if(v.eq.0.0d0) then
         getz = 0.0d0
      else
         getz = ((1.d0 + v/cee)/(1.d0-v/cee))**0.5d0 - 1.d0
      end if
      return
      end
c-----------------------------------------------------------------------
c
c     Find velocity (km/s) as a function of z
c
      double precision function getvel(z)
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision z, cee
      cee = 299792.458d0
      if(z.eq.0.0d0) then
         getvel = 0.0d0
      else
         getvel = cee*((2.d0*z + z**2.d0)/(2.d0 + 2.d0*z + z**2.d0))
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
c     Correct Magnitudes at a given z, to convert from different
c     cosmologies. Not entirely tested !
c
      double precision function cormag2(mag,h_in,h_out,omega_in,
     *     omega_out,lambda_in, lambda_out,z)
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-27
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision h_in,h_out,omega_in,omega_out, z, 
     *     distance_modulus, distmod1, distmod2
      double precision mag, lambda_in, lambda_out
      double precision  h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
      common /cosmology/ h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
c
      h0       = h_in
      omega_m  = omega_in
      omega_l  = lambda_in
      distmod1 =  distance_modulus(z)
c
      h0       = h_out
      omega_m  = omega_out
      omega_l  = lambda_out
      distmod2 =  distance_modulus(z)
      cormag2  = mag + (distmod1-distmod2)
      return
      end
c
c
c-----------------------------------------------------------------------
c     Correct the densities
c
      double precision function phicor2(phi_in,h_in,h_out,omega_in,
     *     omega_out,lambda_in, lambda_out,z)
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision phi_in,h_in,h_out,omega_in,
     *     omega_out, z, dvdz1, dvdz2, vol_elm
      double precision lambda_in, lambda_out
c
      double precision  h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
      common /cosmology/ h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
c
      h0      = h_in
      omega_m = omega_in
      omega_l = lambda_in
      dvdz1   = vol_elm(z, 1.d0)
c
      h0      = h_out
      omega_m = omega_out
      omega_l = lambda_out
      dvdz2   = vol_elm(z, 1.d0)
c
c     This seems to do better (Comparing Poli et al. 2001 fits)
c
      phicor2  = phi_in * (dvdz1/dvdz2) 
      return
      end
c
c-----------------------------------------------------------------------
c     Correct the densities considering a shell; untested !
c
      double precision function phicor3(phi_in,h_in,h_out,omega_in,
     *     omega_out,lambda_in, lambda_out,z1, z2)
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision phi_in,h_in,h_out,omega_in, omega_out, z1,z2,
     *     v_in, zvolumen, v_out
      double precision lambda_in, lambda_out
      double precision  h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
      common /cosmology/ h0, omega_m, omega_l, omega_k, omega_r, 
     *     omega_nu, omega_cm, w, wprime
c
c      z = z2
      h0      = h_in
      omega_m = omega_in
      omega_l = lambda_in
      v_in    = zvolumen(1.0d0, z1 ,z2)
c
      h0      = h_out
      omega_m = omega_out
      omega_l = lambda_out
      v_out   = zvolumen(1.0d0, z1 ,z2)
      phicor3  = phi_in * (v_in/v_out) 
      return
      end
c
c----------------------------------------------------------------------
c
      double precision function zfromr(r)
c
c     Given a proper distance in Mpc, find the corresponding redhshift
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-27
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision r, rdist, tol, eps, zmin, zmax, zeroin, 
     *     get_z_from_r
      common /func_param/ rdist
      external get_z_from_r
      rdist = r
      tol = 1.0d-09
      eps = 1.0d-16
      zmin = 0.d0
      zmax = 16.d0
      zfromr=zeroin(zmin,zmax,get_z_from_r, tol,eps)
      return
      end
c     
c----------------------------------------------------------------------
c
      double precision function get_z_from_r(z)
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision z, dl, r, dm
      common /func_param/r
      external dm
      get_z_from_r = dl(z) - r
      return
      end
c
c----------------------------------------------------------------------
c
      double precision function z_from_tform(time)
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision time, t, tol, eps, zmin, zmax, zeroin
      common /func_param/ t
      external get_z_from_tform
c     Once this fraction of a Hubble time is reached z = 0
c     (the bissection method fails and z < 9.9d-05)
      if(time.ge.0.9635788d0) then
         z_from_tform = 0.d00
      else
         t = time
         tol  = 1.0d-8
         eps  = 1.0d-16
         zmin =  0.0d0
         zmax = 300.d0
         z_from_tform = zeroin(zmin,zmax,get_z_from_tform, tol,eps)
      end if
      return
      end
c     
c----------------------------------------------------------------------
c
      double precision function get_z_from_tform(z)
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision z, t, tform
      common /func_param/t
      external tform
      get_z_from_tform = tform(z) - t
c      print *,  z, tform(z), t
      return
      end
c     
c----------------------------------------------------------------------
c
      double precision function flux_emitted(f_nu, lambda, kcorr, z)
c
c     f_nu in Jy, lambda in microns, where lambda is the 
c     characteristic wavelength of the filter. Emitted flux
c     in solar luminosities
c     Constants from 
c     http://www.nist.gov/pml/div684/fcdc/codata.cfm
c     published in 2010, Reviews of Modern Physics, 84, 1527
c     Available at:
c     http://www.nist.gov/pml/div684/fcdc/upload/Wall-chart-3.pdf
c     and 
c     http://pdg.lbl.gov/2014/reviews/rpp2014-rev-astrophysical-constants.pdf
c     http://pdg.lbl.gov/2014/reviews/contents_sports.html

c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-03-19
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision f_nu, z, r_mpc, cee, four_pi, constant, dist_mpc,
     *     dl
      double precision lsol,jy_erg, micron_cm, lambda, kcorr
c   
      lsol    = 3.828d33        ! solar luminosity in erg/s
      r_mpc   = 3.08567758149d24 ! Mpc --> cm
      cee     = 2.99792458d10   ! cm/s
      four_pi = 4.0d0 * dacos(-1.0d0)
c     
      jy_erg    = 1.0d-23
      micron_cm = 1.0d-04
c
      constant = four_pi * cee * jy_erg * r_mpc**2
      constant = constant/(micron_cm*lsol)
c
      dist_mpc     = dl(z)
      flux_emitted = f_nu * dist_mpc * dist_mpc / lambda
      flux_emitted = constant * flux_emitted
      return
      end
c     
c----------------------------------------------------------------------
c
      double precision function flux_observed(lum, lambda, z, kcorr)
c
c     f_nu in Jy, lambda in microns, where lambda is the 
c     characteristic wavelength of the filter. Emitted flux
c     in solar luminosities
c
c
c     Constants from 
c     http://www.nist.gov/pml/div684/fcdc/codata.cfm
c     published in 2010, Reviews of Modern Physics, 84, 1527
c     Available at:
c     http://www.nist.gov/pml/div684/fcdc/upload/Wall-chart-3.pdf
c     and
c     http://pdg.lbl.gov/2014/reviews/rpp2014-rev-astrophysical-constants.pdf
c     http://pdg.lbl.gov/2014/reviews/contents_sports.html
c
c
c     (C) Copyright 2015 Christopher Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-19
c     2015-03-18
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision z, r_mpc, cee, four_pi, constant, dist_mpc, dl
      double precision lsol,jy_erg, micron_cm, lambda,lum,kcorr
c   
      lsol    = 3.828d33        ! solar luminosity in erg/s
      r_mpc   = 3.08567758149d24 ! Mpc --> cm
      cee     = 2.99792458d10   ! cm/s
      four_pi = 4.0d0 * dacos(-1.0d0)
c     
      jy_erg    = 1.0d-23
      micron_cm = 1.0d-04
c
      constant = four_pi * cee * jy_erg * r_mpc**2
      constant = (micron_cm*lsol)/constant
c
      dist_mpc      = dl(z)
      flux_observed = lum * lambda/(dist_mpc * dist_mpc)
      flux_observed = constant * flux_observed
c      print *, 'ncosmo:flux_observed ', constant, lum, lambda,z
      return
      end
c
