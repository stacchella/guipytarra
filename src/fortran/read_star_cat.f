c
c     Read star catalogue with ID, RA, DEC, list of mags
c
      subroutine read_star_cat(star_catalogue, filters_in_cat)
c     
      implicit none
      integer filters_in_cat, max_stars, nstars, i, j, idstar, nfilters,
     &     sca_id
      integer colcornr, rowcornr, naxis1, naxis2
      integer ixmin, ixmax, iymin,iymax
      character subarray*8
      double precision x_stars, y_stars,  mag_stars
      double precision rra, ddec, array
      double precision ra_dithered, dec_dithered, pa_degrees, xc, yc, 
     &     osim_scale
      character star_catalogue*120
c
      parameter(max_stars=10000, nfilters = 54)
c
      dimension x_stars(max_stars), y_stars(max_stars), 
     *     mag_stars(max_stars, nfilters)
      dimension array(nfilters)
c
      common /stars/ x_stars, y_stars, mag_stars, nstars
c
      open(1,file=star_catalogue)
      nstars = 0
      do i = 1, max_stars
         read(1,110,end= 1000) idstar, rra, ddec,
     *        (array(j),j=1,filters_in_cat)
 110     format(i5,2(2x, f16.12), 2(1x,f8.3), 58(2x,f8.3))
c
         nstars = nstars + 1
         ra_stars(nstars) = rra
         dec_stars(nstars) = ddec
c         
         do j = 1, filters_in_cat
            mag_stars(nstars,j) = array(j)
         end do
c
 100     continue
      end do
      print 120, nstars, star_catalogue
 120  format('read ',i6,' objects from',/,a80)
 1000 close(1)
      return
      end
