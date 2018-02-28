      subroutine read_parameters(filename_param2, verbose,
     *     brain_dead_test, module, mode, ngroups, nframe,
     *     subarray, colcornr, rowcornr, naxis1, naxis2,
     *     camera, primary_dither, subpixel_dither, 
     *     number_primary, number_subpixel, 
     *     dither_arc_sec, ra0, dec0, pa_degrees, 
     *     include_ktc, include_dark, include_readnoise, 
     &     include_reference,
     *     include_1_over_f, include_latents, include_non_linear,
     *     include_cr, cr_mode, include_bg, bkg_mode,
     *     zodiacal_scale_factor,
     *     include_stars, nstars, star_catalogue,
     *     include_galaxies, ngal, include_cloned_galaxies, 
     *     galaxy_catalogue,
c     *     h0, omega_m, omega_l, omega_r, omega_nu, w, wprime,
c     *     phistar, mstar, alpha, q_evol, p_evol, zmin, zmax,
c     *     bright, faint, apm_bright, apm_faint, solid_angle,
     *     nf, use_filter, cat_filter)
      
      implicit none
      character subarray*8
      character filename_param2*180
c
      character mode*10, galaxy_catalogue*180, star_catalogue*180
      character primary_dither*20, subpixel_dither*20, camera*20,
     &     module*20
      character path_guitarra*100
c
      integer number_primary, number_subpixel
      integer verbose, ngroups, nframe, ndither, nf,
     *     use_filter, i, brain_dead_test, cat_filter
      integer include_ktc, include_dark, include_readnoise,
     *     include_non_linear, include_latents,include_1_over_f, 
     &     include_reference,
     *     include_cr, include_bg, include_stars,
     *     include_galaxies, 
     *     include_cloned_galaxies
      integer cr_mode, bkg_mode
      integer colcornr, rowcornr, naxis1, naxis2
      integer nstars, ngal, junk
c
      double precision ra0, dec0, pa_degrees, dither_arc_sec
      double precision h0, omega_m, omega_l, omega_r, omega_nu,
     *     w, wprime
      double precision  bright, faint, apm_bright, apm_faint, 
     *     zmin, zmax,solid_angle
      double precision rcore_arcsec
      double precision mstar, alpha, phistar, p_evol, q_evol 
      double precision zodiacal_scale_factor, bkg
c
      dimension use_filter(54),cat_filter(54)
c
      open(1,file=filename_param2)
      read(1,10) verbose
 10   format(i12)
      print *, 'verbose = ', verbose
      read(1,10) brain_dead_test
      print *, 'brain_dead_test     ', brain_dead_test
c
      read(1,11) module
 11   format(a20)
      print *, 'module  = ', module
      read(1,15) mode
 15   format(a10) 
      print 15, mode
c
      read(1,10) ngroups
      print *, 'ngroups ', ngroups
c      read(1,10) nframe
c      print *, 'nframe  ', nframe
c
      read(1, 17) subarray
 17   format(a8)
      print 18, subarray
 18   format('subarray is ',a8)
      read(1,10) colcornr
      print *, 'colcornr ', colcornr
      read(1,10) rowcornr
      print *, 'rowcornr ', rowcornr
      read(1,10) naxis1
      print *, 'naxis1 ', naxis1
      read(1,10) naxis2
      print *, 'naxis2 ', naxis2
c
c     dither section
c
      read(1,10) number_primary
      print *, 'number of primary dithers', number_primary
      read(1,10) number_subpixel
      print *, 'number of sub-pixel dithers', number_subpixel
      read(1,19) camera
 19   format(a20)
      print 21, 'camera ', camera
 21   format(a20,2x, a20)
      read(1,19) primary_dither
      print 21, 'primary_dither', primary_dither
      read(1,19) subpixel_dither
      print 21, 'subpixel_dither', subpixel_dither
c      read(1,20) dither_arc_sec
 20   format(f12.7)
c      print *, 'dither_arc_sec ', dither_arc_sec
      read(1,20) ra0
      print *, 'ra0 ', ra0
      read(1,20) dec0
      print *, 'dec0 ', dec0
      read(1,20) pa_degrees
      print *, 'PA ', pa_degrees
c
      read(1,10) include_ktc
      print *,'include_ktc        ', include_ktc
      read(1,10) include_dark
      print *,'include_dark       ', include_dark
      read(1,10) include_readnoise
      print *,'include_readnoise  ', include_readnoise
      read(1,10) include_reference
      print *,'include_reference  ', include_reference
      read(1,10) include_non_linear
      print *,'include_non_linear ', include_non_linear
      read(1,10) include_latents
      print *,'include_latents    ', include_latents
      read(1,10) include_1_over_f
      print *,'include_1_over_f   ', include_1_over_f
c
      read(1,10) include_cr
      print *,'include_cr         ', include_cr
      read(1,10) cr_mode
      print *,'cr_mode            ', cr_mode
c
      read(1,10) include_bg 
      print *,'include_bg               ', include_bg
      read(1,10) bkg_mode
      print *,'bkg_mode                 ', bkg_mode
      read(1,20) zodiacal_scale_factor
      print *,'zodiacal_scale_factor    ',zodiacal_scale_factor
c
      read(1,10) include_stars
      print *, 'include_stars             ', include_stars
      read(1,10) nstars
      print *, 'nstars                    ', nstars
      read(1, 40) star_catalogue
      print 50, 'read star   catalogue    ', star_catalogue
c
      read(1,10) include_galaxies
      print *, 'include_galaxies          ', include_galaxies
      read(1,10) ngal
      print *, 'ngal                     ', ngal
c      read(1,10) include_cloned_galaxies
c      print *, 'include_cloned_galaxies  ', include_cloned_galaxies
      read(1, 40) galaxy_catalogue
 40   format(a80)
      print 50, 'read galaxy catalogue    ', galaxy_catalogue
 50   format(a30,2x,a80)
c     
c      read(1,20) h0
c      print *, 'h0                      ',h0
c      read(1,20) omega_m
c      print *, 'omega_m                 ',omega_m
c      read(1,20) omega_l
c      print *, 'omega_l                 ',omega_l
c      read(1,20) omega_r
c      print *, 'omega_l                 ',omega_r
c      read(1,20) omega_nu
c      print *, 'omega_l                 ',omega_nu
c      read(1,20) w
c      print *, 'w                       ',w
c      read(1,20) wprime
c      print *, 'wprime                  ',wprime
cc 
c      read(1,20) phistar
c      print *, 'phistar                 ',phistar
c      read(1,20) mstar
c      print *, 'mstar                   ',mstar
c      read(1,20) alpha
c      print *, 'alpha                   ',alpha
c      read(1,20) q_evol
c      print *, 'q_evol                  ',q_evol
c      read(1,20) p_evol
c      print *, 'p_evol                  ',p_evol
cc
c      read(1,20) zmin
c      print *, 'zmin                    ',zmin
c      read(1,20) zmax
c      print *, 'zmax                    ',zmax
cc
c      read(1,20) bright
c      print *, 'bright                  ',bright
c      read(1,20) faint
c      print *, 'faint                   ',faint
c      read(1,20) apm_bright
c      print *, 'apm_bright              ',apm_bright
c      read(1,20) apm_faint
c      print *, 'apm_faint               ',apm_faint
cc
c      read(1,20) solid_angle
c      print *, 'solid_angle             ', solid_angle
c
      read(1,10) nf
      print *, 'number of filters        ', nf
      do i = 1, nf
         read(1,10) use_filter(i)
         print *,'filter index             ', use_filter(i)
         read(1,10) junk
         print *,'catalogue filter index   ', junk
         cat_filter(i) = junk
c     read(1,10) cat_filter(i)
c         print *,'filter index             ', cat_filter(i)
      end do
      close(1)
      return
      end
