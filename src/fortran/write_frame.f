c
      subroutine write_frame(filename, iunit, bitpix, level, 
     *     naxis, naxes)
c      subroutine write_frame(filename, iunit, bitpix, level, 
c     *     naxis, naxes, scratch, n_image_x, n_image_y)
c
c     write/add an image "plane" to a data cube, at level "level";
c     In this way, there is no need to store the entire data cube
c     in memory. The data cube file (accessed through logical unit
c     "iunit" already must be open.
c
c     cnaw 2015-01-27
c     Steward Observatory, University of Arizona
c
      implicit none
      character filename*120
      integer   bitpix,status, fpixels, lpixels, group, level, naxis,
     *     naxes, nnn, n_image_x, n_image_y, iunit, i, j, max_order,
     *     order
      real    image, scratch, well_depth, bias, linearity, accum
      real    value, tmax
      double precision linearity_gain, lincut, well_fact
      integer(kind=4) int_image
      parameter (nnn=2048,max_order=7)
      dimension scratch(nnn,nnn), image(nnn,nnn),accum(nnn,nnn)
      dimension int_image(nnn,nnn), naxes(3), fpixels(3), lpixels(3)
      dimension well_depth(nnn,nnn), linearity(nnn,nnn,max_order), 
     *     bias(nnn,nnn)
c     
      common /well_d/ well_depth, bias, linearity, linearity_gain,
     *     lincut, well_fact, order
c
      common /images/ accum, image, n_image_x, n_image_y
      common /scratch_/ scratch 
c
c     open file (which should have been initialised already)
c     
c      print *,'enter write_frame:', iunit, bitpix,naxes,level
c      print *, naxis, naxes, n_image_x, n_image_y
c
c     starting point 
c
      status  = 0
      group   = 0
      fpixels(1) = 1
      fpixels(2) = 1
      fpixels(3) = level
      lpixels(1) = nnn
      lpixels(2) = nnn
      lpixels(3) = level
      tmax       = 2.**16-1.0
c
c     Unsigned integers (bitpix has to be 20)
c     note that CFITSIO will write the integer (kind=4) array as
c     integer (kind=2). 
c     
c     Compare the number of counts (which are in ADU) with
c     well-depth for individual pixels.
c
      if(bitpix.eq.16.or.bitpix.eq.20) then
         do j = 1, n_image_y
            do i = 1, n_image_x
               value = scratch(i,j)
               if(value.gt.well_depth(i,j) .and.
     *              well_depth(i,j).gt.tmax) value = well_depth(i,j)
c     
c     This gets rid of NaNs - otherwise arrays will trigger the following:
c     FITSIO Error Status =         412 : datatype conversion overflow
c
               if(value.ne.value) then
                  value = 0.0
               end if
c     further cleaning up - values above 65K and < 0
               if(value .ge. tmax)  value = tmax-1.0
               if(value .lt. 0.0)   value = 1.0e-6
               int_image(i,j) = int(value)
c               print *, i, j, value, int_image(i,j)
            end do
         end do
         status = 0
         call ftpssj(iunit, group, naxis, naxes, fpixels, lpixels, 
     *        int_image,status)
      end if
c     
c     integer*4
c
      if(bitpix.eq.32) then
         do j = 1, n_image_y
            do i = 1, n_image_x
               int_image(i,j) = nint(scratch(i,j))
            end do
         end do
         call ftpssj(iunit, group, naxis, naxes, fpixels, lpixels, 
     *        int_image,status)
      end if
c
c     floating point
c
      if(bitpix.eq.-32) then
         call ftpsse(iunit, group, naxis, naxes, fpixels, lpixels, 
     *        scratch,status)
      end if
      if (status .gt. 0) then 
         call printerror(status)
         print *, "pause"
         read(*,'(A)')
      endif
c      print *,' exit read_fits: write_frame'
      return
      end
