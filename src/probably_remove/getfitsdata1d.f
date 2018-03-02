c
c-----------------------------------------------------------------------
c
      subroutine getfitsdata1d(unit,x0,dx,x,y,npixels)
      parameter (nnn=20000)
      integer unit, status, fpixel,naxes(3), group
c      logical*1 chabu
      logical chabu
      real x, y, x0, dx
      dimension x(nnn), y(nnn)
      status = 0
C     determine the size of the image
      call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)

C     check that it found both NAXIS1 and NAXIS2 keywords
c 5    if (nfound .ne. 3)then
c      print *,' This is the number of NAXISn keywords.', 
c     *     nfound, (naxes(i),i=1,nfound)
      if(nfound.eq.1) then
         npixels = naxes(1)
         npts = naxes(1)
      else
         npixels=naxes(1)*naxes(2)
         npts = npixels
      end if
      group=1
      fpixel = 1
      call ftgpve(unit, group, fpixel, npixels,0, y, chabu,status)
      do i=1, npts
         x(i) = x0 + (i-1) * dx
      end do
      if(status .ne.0) call printerror(status)
      return
      end
