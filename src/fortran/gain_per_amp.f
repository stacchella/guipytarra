      subroutine gain_per_amp(amp_gain)
      implicit none
      double precision amp_gain, array, xmed, xsig, xmed1, xsig1
      real gain_image
      integer nnn, i, j, n, debug
c
      parameter(nnn=2048)
c
      dimension  gain_image(nnn,nnn), amp_gain(4), array(nnn*512)
c
      common /gain_/ gain_image
      
      debug = 0
      n     = 0
      do j = 5, 2044
         do i = 5, 512
            if(gain_image(i,j)+0.0d0 .eq. gain_image(i,j)) then
               n = n + 1 
               array(n) = gain_image(i,j)
            end if
         end do
      end do
c
      call xbiwt(array, n, xmed, xsig, xmed1, xsig1, debug)
      amp_gain(1) = xmed
c
      n  = 0
      do j = 5, 2044
         do i = 513, 1024
            n = n + 1 
            array(n) = gain_image(i,j)
         end do
      end do
c
      call xbiwt(array, n, xmed, xsig, xmed1, xsig1, debug)
      amp_gain(2) = xmed
c
      n     = 0
      do j = 5, 2044
         do i = 1025, 1536
            n = n + 1 
            array(n) = gain_image(i,j)
         end do
      end do
c
      call xbiwt(array, n, xmed, xsig, xmed1, xsig1, debug)
      amp_gain(3) = xmed
c
      n     = 0
      do j = 5, 2044
         do i = 1537, 2044
            n = n + 1 
            array(n) = gain_image(i,j)
         end do
      end do
c
      call xbiwt(array, n, xmed, xsig, xmed1, xsig1, debug)
      amp_gain(4) = xmed
      do i = 1, 4
         amp_gain(i) = 1.d0
      end do
      return
      end
