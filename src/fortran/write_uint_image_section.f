c
c----------------------------------------------------------------------     
c
      subroutine write_uint_image_section(iunit, image, nx, ny, 
     *     nz, level, naxis, naxes)
      implicit none
      integer(kind=2) image
      integer fpixels, lpixels, nnn, naxis, naxes, iunit, nx, ny, nz,
     *     level
      integer status, group
      parameter (nnn=2048)
c     
      dimension image(nnn,nnn),naxes(3), fpixels(3), lpixels(3)
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
      call ftpssi(iunit, group, naxis, naxes, fpixels, lpixels, image,
     *     status)
      if (status .gt. 0) then 
         call printerror(status)
      endif
      return
      end
