c
c-----------------------------------------------------------------------
c
c     Create a galaxy catalogue with positions, magnitude, B/T etc.
c     following a Schechter function
c
      subroutine create_galaxy_catalogue(option,ngal, ra0, dec0, pa_v3,
     *     verbose, seed)
      implicit none
      double precision ra0, dec0, pa_v3
      double precision ra, dec, z
      double precision magnitude, nsersic, ellipticity, re, theta, 
     *     flux_ratio
      double precision ran, pi, dm, mmax, dtheta
      double precision zbqlu01, ab_mag_to_photon_flux, distance_modulus
      double precision mstar, alpha, phistar, p_evol, q_evol
      double precision bright, faint, apm_bright, apm_faint, 
     *     zmin, zmax,solid_angle
      double precision lf_scale, mag_array, lumfun, int_lf
      double precision z_scale, zdist, expected, int_zdist
c
      double precision v2_min, v2_max, v3_min, v3_max, osim_scale,
     *     x_osim, y_osim

      double precision x0, y0, x_sca, y_sca, xc, yc, xg, yg,
     *     ra_c, dec_c
 
c       double precision effective_wl_nircam, width_nircam, 
c     *     system_transmission
c
      integer option, ngal, verbose, seed
      integer nbins, nzdist, i, j, max_objects, nnn, nfilters
      integer ii, nsub, ncomponents, indx, nc
c
      parameter (max_objects= 50000,nnn=2048, nfilters=54, nsub=4)

c      dimension effective_wl_nircam(nfilters), width_nircam(nfilters), 
c     *     system_transmission(nfilters)
      dimension mag_array(nnn), int_lf(nnn), lumfun(nnn)
      dimension zdist(nnn), int_zdist(nnn), expected(nnn)
      dimension ra(max_objects), dec(max_objects), z(max_objects),
     *     magnitude(max_objects,nfilters), ncomponents(max_objects),
     *     nsersic(max_objects, nsub), ellipticity(max_objects, nsub),
     *     re(max_objects, nsub), theta(max_objects, nsub), 
     *     flux_ratio(max_objects, nsub)
c
      common /schechter_params/ mstar, alpha, phistar, p_evol, q_evol
      common /sample_params/   bright, faint, apm_bright, apm_faint, 
     *     zmin, zmax,solid_angle
      common /galaxy/ra, dec, z, magnitude, nsersic, ellipticity, re,
     *     theta, flux_ratio, ncomponents
c
      data xc, yc /0.0d0, -315.6d0/
      data osim_scale/1.59879d0/ ! arc sec average for X and Y
c
      pi      = dacos(-1.0d0)
c     
c-----------------------------------------------------------------------
c     Distribution of magnitudes
c     
      nbins = 200
      call simple_schechter_lf(mag_array,lumfun,int_lf, nbins)
      lf_scale = int_lf(nbins)
c     
c     Calculate the expected redshift distribution
c     
      nzdist = (zmax-zmin)/0.01d0
      call zdistribution (zdist,expected,int_zdist,nzdist)
      print *, 'z', zdist(nzdist),int_zdist(nzdist), nzdist
      z_scale = int_zdist(nzdist)
c
c-----------------------------------------------------------------------
c
c     Initialise variables to assign RA, Dec coordinates
c
      if(verbose.eq.1) then
         print * , 'option ', option
         print *, 'sca id (481-490), or module (0 = both, 1=A, 2=B)'
      end if
c     
c     for this configuration, find range in V2, V3
c
      call sca_boundaries(option, v2_min, v2_max, v3_min, v3_max)
      print *, option, v2_min, v2_max, v3_min, v3_max 
c     
c     this provides coordinates of the SCA centre
c     
      x0 = (v2_min + v2_max)/2.0d0
      y0 = (v3_min + v3_max)/2.0d0
c     
c     find ra_c dec_c corresponding to x_sca, y_sca
c     
      if(option.ge.481 .and. option.le.490) then
         x_sca  = 1024.d0
         y_sca  = 1024.d0
         call sca_to_ra_dec(option,  ra0, dec0,ra_c, dec_c,
     *        pa_v3, xc, yc, osim_scale, x_sca, y_sca)
         call osim_coords_from_sca(option, x_sca, y_sca, x_osim, y_osim)
      end if
c     
c     for mod A use 485 
c
      if(option .eq.1) then
         x_sca  = 1024.d0
         y_sca  = 1024.d0
         call sca_to_ra_dec(485,  ra0, dec0, ra_c, dec_c,
     *        pa_v3, xc, yc, osim_scale, x_sca, y_sca)
         call osim_coords_from_sca(485, x_sca, y_sca, x_osim, y_osim)
      end if
c
c     for mod B use 490 
c     
      if(option .eq.2) then
         x_sca  = 1024.d0
         y_sca  = 1024.d0
         call sca_to_ra_dec(490, ra0, dec0, ra_c, dec_c,
     *        pa_v3, xc, yc, osim_scale, x_sca, y_sca)
         call osim_coords_from_sca(495, x_sca, y_sca, x_osim, y_osim)
      end if
c     
c     otherwise, use XC, YC
c
      if(option .eq.0) then
         ra_c  = ra0
         dec_c = dec0
      end if
      
      print 90,x0, y0, ra0, dec0, ra_c, dec_c
 90   format(' X0, Y0',2(2x,f8.2),'  ra0, dec0 ',2(2x,f16.10),
     *     '  ra_c, dec_c ',2(2x,f16.10))
      if(option .ne.0) then
         print *, option,pa_v3, xc, yc, osim_scale, x_sca, y_sca,
     *        x_osim, y_osim
      end if
c
c-----------------------------------------------------------------------
c
c     Create list of galaxies
c     
      open(1,file= 'galaxy.cat')
      do i = 1, ngal
c
c     find a random redshift given the z distribution
c
         ran  = zbqlu01(seed) * z_scale
         call linear_interpolation(nzdist, int_zdist, zdist, ran,
     *        z(i), nnn)
         dm      = distance_modulus(z(i))
c    
c     Get a random magnitude. Finding the magnitude is modulated by 
c     the redshift distribution, thus restrict the search to the absolute
c     magnitude range at z(k)
c     
         mmax       = apm_faint  - dm
         call find_index(nbins, mag_array, mmax, ii, nnn)
         lf_scale  = int_lf(ii)
         ran  = zbqlu01(seed) * lf_scale
         call linear_interpolation(nbins, int_lf,mag_array, ran,
     *        magnitude(i,1), nnn)
         magnitude(i,1) = magnitude(i,1) + dm
         magnitude(i,1) = 17.0d0
c     
c     set structural parameters
c     
         nc               = 0
c
c         nc              = nc + 1
c         flux_ratio(i,nc) = 0.6d0
c         nsersic(i,nc)    = 0.5d0
c         re(i,nc)         = 0.2d0
c         theta(i,nc)      = 0.0d0
c
         nc               = nc + 1
         flux_ratio(i,nc)  = 0.2d0 + zbqlu01(seed) * 0.8d0
         nsersic(i,nc)    = 4.d0
         re(i,nc)         = 2.d0 + zbqlu01(seed) * 4.d0
         theta(i,nc)      = zbqlu01(seed)* 360.d0
         if(flux_ratio(i,nc) .eq.1.0d0) then
            ncomponents(i) = 1
         else
            nc            = nc + 1
            nsersic(i,nc)     = 1.d0
            flux_ratio(i, nc) = 1.d0 - flux_ratio(i,nc-1)
            re(i,nc)          = re(i,nc-1) + zbqlu01(seed)* 2.d0
            dtheta            = -10.d0 + zbqlu01(seed)*20.d0
            theta(i,nc)       = theta(i,nc) +dtheta
            ncomponents(i)       = nc
         end if
c     
c     set some coordinates
c     
         xg           = 5.d0 + zbqlu01(seed) * 2039.d0
         yg           = 5.d0 + zbqlu01(seed) * 2039.d0
         call sca_to_ra_dec(option, ra0, dec0, ra(i), dec(i),
     *        pa_v3, xc, yc, osim_scale, xg, yg)
         call ra_dec_xy_osim(ra0, dec0, ra(i), dec(i),
     *        pa_v3, xc, yc, osim_scale, x_osim, y_osim)
         write(1,110) i, ra(i), dec(i), x_osim, y_osim,xg, yg,
     *        magnitude(i,1), nsersic(i,1), re(i,1), theta(i,1),
     *        flux_ratio(i,1), z(i)
 110     format(i5,2(2x, f16.12), 4(1x,f8.3), 29(2x,f8.3))
      enddo
      close(1)
      return
      end
