c
c----------------------------------------------------------------------     
c
      subroutine write_float_image_section(iunit, image, nx, ny, 
     *     nz, level, naxis, naxes)
      implicit none
      real image
      integer fpixels, lpixels, iunit, nx, ny, nz, level, naxis, naxes
      integer status, group, nnn
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
      call ftpsse(iunit, group, naxis, naxes, fpixels, lpixels, image,
     *     status)
      if (status .gt. 0) then 
         call printerror(status)
      endif
      return
      end
c
