c     
c     Create a composite dither pattern by adding sub-pixel dithers
c     to a "major" dither.
c
c     cnaw@as.arizona.edu
c     2015-10-07
c
      subroutine dither_arrays(intrasca_x, intrasca_y, number_primary,
     *     small_x, small_y, number_subpixel, size1, size2,  
     *     xdither, ydither, nd)
      implicit none
      double precision xdither, ydither, intrasca_x, intrasca_y,
     *     small_x, small_y, size1, size2
      integer number_primary,number_subpixel, nd
      integer k, j
c
      dimension small_x(9,9), small_y(9,9), 
     *     intrasca_x(25), intrasca_y(25)
      dimension xdither(1000), ydither(1000)
c      print *,'dither_arrays: number of dithers: primary, sub-pixel',
c     *     number_primary,number_subpixel
c
      nd = 0
      if(number_subpixel .gt. 9) number_subpixel =  9
      if(number_primary  .gt.25) number_primary  = 25
      do k = 1, number_primary
         do j = 1, number_subpixel
            nd = nd + 1
            xdither(nd) = intrasca_x(k) *size1 + 
     *           small_x(number_subpixel, j) * size2
            ydither(nd) = intrasca_y(k) * size1 + 
     *           small_y(number_subpixel, j) * size2
            xdither(nd) = xdither(nd) 
            ydither(nd) = ydither(nd) 
            print 120, nd, intrasca_x(k),intrasca_y(k),
     *           small_x(number_subpixel,j),small_y(number_subpixel,j),
     *            xdither(nd), ydither(nd)
 120        format('NIRCam dithers ', i4, 6(1x,f12.6))
         end do
      end do
      return
      end

c
