      subroutine read_star_catalogue(star_catalogue, nfilters, verbose)
c     &     subarray, colcornr, rowcornr, naxis1, naxis2,
c     &     ra_dithered, dec_dithered, pa_degrees, xc, yc, osim_scale,
c     &     sca_id, old_style, verbose)
      implicit none
      integer nfilters, max_stars, nstars, i, j, indx, max_filters,
     &     sca_id, old_style, verbose
      integer colcornr, rowcornr, naxis1, naxis2
      integer ixmin, ixmax, iymin,iymax
      character subarray*8
      double precision x_stars, y_stars,  mag_stars
      double precision rra, ddec, x_sca, y_sca, array
      double precision ra_dithered, dec_dithered, pa_degrees, xc, yc, 
     &     osim_scale
      character star_catalogue*120
c
      parameter(max_stars=10000, max_filters = 54)
c
      dimension x_stars(max_stars), y_stars(max_stars), 
     *     mag_stars(max_stars, max_filters)
      dimension array(max_filters)
c
      common /stars/ x_stars, y_stars, mag_stars, nstars
c
c
c      if(subarray(1:4) .eq.'FULL') then
c         ixmin = 5
c         ixmax = naxis1- 4
c         iymin = 5
c         iymax = naxis2 - 4
c      else
c         ixmin = colcornr
c         ixmax = ixmin + naxis1
c         iymin = rowcornr
c         iymax = iymin + naxis2
c      end if
c
      open(1,file=star_catalogue)
      nstars = 0
      do i = 1, max_stars
         read(1,110,end= 1000) indx, rra, ddec, x_sca, y_sca,
     *        (array(j),j=1,nfilters)
 110     format(i5,2(2x, f16.12), 2(1x,f8.3), 58(2x,f8.3))
c
c     verify that star is contained within field
c     (could be done later for objects at the edges)
c     
c         if(old_style.eq.1) then
c            call ra_dec_to_sca(sca_id, 
c     *           ra_dithered, dec_dithered, 
c     *           rra, ddec, pa_degrees, 
c     *           xc, yc,  osim_scale, x_sca, y_sca)
c            if(x_sca.lt.ixmin .or. x_sca.gt.ixmax) go to 100
c            if(y_sca.lt.iymin .or. y_sca.gt.iymax) go to 100
c         end if
         nstars = nstars + 1
         x_stars(nstars) = x_sca
         y_stars(nstars) = y_sca
         do j = 1, nfilters
            mag_stars(nstars,j) = array(j)
         end do
         if(verbose.gt.1) then
             print 110,indx, rra, ddec, x_sca, y_sca,array(1)
         end if
 100     continue
      end do
 1000 close(1)
      return
      end
