c
c     Create dither moves for a given set pattern selected from tables
c     computed by J. Anderson. This is valid for a single pointing
c     (i.e., not an entire survey)
c
      subroutine nircam_dithers(camera, primary_dither, 
     *     subpixel_dither, number_primary, number_subpixel,
     *     xdither, ydither, nd, nnn)
c
      implicit none
c
      double precision xdither, ydither
      integer nd, nnn, number_primary, number_subpixel
      character primary_dither*20, subpixel_dither*20, camera*20
c
      double precision smallx, smally, sca_x, sca_y, full3x,
     *     full3y, tight3x, tight3y, full6x, full6y,full9x, full9y,
     *     full15x, full15y, full21x,full21y, full45x,full45y,
     *     intrasca_small_x,intrasca_small_y,
     *     intrasca_medium_x,intrasca_medium_y,
     *     intrasca_large_x,intrasca_large_y,
     *     gen_x, gen_y
c
      double precision size1, size2,small_x, small_y
c
      dimension xdither(nnn), ydither(nnn), small_x(nnn), small_y(nnn)
c
      dimension smallx(9,9), smally(9,9), sca_x(16), sca_y(16), 
     *     full3x(3), full3y(3), tight3x(3),tight3y(3), 
     *     full6x(6), full6y(6), full9x(9), full9y(9),
     *     full15x(15), full15y(15), full21x(21),full21y(21),
     *     full45x(45),full45y(45),
     *     intrasca_small_x(25),intrasca_small_y(25),
     *     intrasca_medium_x(25),intrasca_medium_y(25),
     *     intrasca_large_x(25),intrasca_large_y(25),
     *     gen_x(64), gen_y(64)
c     
      print 45, camera
      print 45, primary_dither
      print 45, subpixel_dither
 45   format(a20)
c     
c     read J. Anderson's table of dithers. Returned values are in Arc Sec
c     
      call read_dithers(smallx, smally, sca_x, sca_y, full3x,
     *     full3y, tight3x, tight3y, full6x, full6y,full9x, full9y,
     *     full15x, full15y, full21x,full21y, full45x,full45y,
     *     intrasca_small_x,intrasca_small_y,
     *     intrasca_medium_x,intrasca_medium_y,
     *     intrasca_large_x,intrasca_large_y,
     *     gen_x, gen_y)
c     
c
c     Full 3
c     
      if(primary_dither.eq.'full3') then
         if(camera.eq.'short') then
            size1 = 1.d0
            size2 = 1.d0
         else
            size1 = 1.d0        ! this corresponds to pixel size for intrasca
            size2 = 2.d0
         end if
c     
         print *,'number of dithers: primary, sub-pixel',
     *        number_primary,number_subpixel
         call dither_arrays(full3x, full3y,
     *        number_primary,
     *        small_x, small_y, number_subpixel, size1, size2,
     *        xdither, ydither, nd)
c     
      end if
c     
c     Tight 3
c     
      if(primary_dither.eq.'tight3') then
         if(camera.eq.'short') then
            size1 = 1.d0
            size2 = 1.d0
         else
            size1 = 1.d0     
            size2 = 2.d0
         end if
c     
         print *,'number of dithers: primary, sub-pixel',
     *        number_primary,number_subpixel
         call dither_arrays(tight3x, tight3y,
     *        number_primary,
     *        small_x, small_y, number_subpixel, size1, size2,
     *        xdither, ydither, nd)
c     
      end if
c     
c     Full 6
c     
      if(primary_dither.eq.'full6') then
         if(camera.eq.'short') then
            size1 = 1.d0
            size2 = 1.d0
         else
            size1 = 1.d0
            size2 = 2.d0
         end if
c     
         print *,'number of dithers: primary, sub-pixel',
     *           number_primary,number_subpixel
         call dither_arrays(full6x, full6y,
     *        number_primary,
     *        small_x, small_y, number_subpixel, size1, size2,
     *        xdither, ydither, nd)
c     
      end if
c     
c     Full 9
c     
      if(primary_dither.eq.'full9') then
         if(camera.eq.'short') then
            size1 = 1.d0
            size2 = 1.d0
         else
            size1 = 1.d0
            size2 = 2.d0
         end if
c     
         print *,'number of dithers: primary, sub-pixel',
     *        number_primary,number_subpixel
         call dither_arrays(full9x, full9y,
     *        number_primary,
     *        small_x, small_y, number_subpixel, size1, size2,
     *        xdither, ydither, nd)
c     
      end if
c     
c     INTRA-MODULE
c     
      if(primary_dither.eq.'intramodule') then
         if(camera.eq.'short') then
            size1 = 1.d0
            size2 = 1.d0
         else
            size1 = 1.d0
            size2 = 2.d0
         end if
c
         if(number_primary.gt.16) number_primary = 16
         print *,'number of dithers: primary, sub-pixel',
     *        number_primary,number_subpixel
         call dither_arrays(sca_x, sca_y,
     *        number_primary,
     *        small_x, small_y, number_subpixel, size1, size2,
     *        xdither, ydither, nd)
c     
      end if
c     
c     INTRA-SCA
c     
      if(primary_dither.eq.'intrasca') then
         if(camera.eq.'short') then
            size1 = 1.d0
            size2 = 1.d0
         else
            size1 = 1.d0
            size2 = 2.d0
         end if
c     
         print *,'number of dithers: primary, sub-pixel',
     *        number_primary,number_subpixel
         if(subpixel_dither .eq.'small') then
            call dither_arrays(intrasca_small_x, intrasca_small_y,
     *           number_primary,
     *           small_x, small_y, number_subpixel, size1, size2,
     *           xdither, ydither, nd)
         end if
c     
         if(subpixel_dither .eq.'medium') then
            call dither_arrays(intrasca_medium_x, intrasca_medium_y,
     *           number_primary,
     *           small_x, small_y, number_subpixel, size1, size2,
     *           xdither, ydither, nd)
         end if
c     
         if(subpixel_dither .eq.'large') then
            call dither_arrays(intrasca_large_x, intrasca_large_y,
     *           number_primary,
     *           small_x, small_y, number_subpixel, size1, size2,
     *           xdither, ydither, nd)
         end if
c     
      end if
c     
      return
      end
