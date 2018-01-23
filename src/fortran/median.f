c
c     CNAW 2013-06-13
c
c     Calculate the median
c
      subroutine median(npts, x, xmed)
      implicit none
      double precision x, xmed, x1, x2, temp
      integer npts, i, npt,lll
      parameter(lll=50000)
c      parameter(lll=67108864)
      dimension temp(lll), x(lll)
      do i = 1, npts
         temp(i) = x(i)
      end do
c
      call dsort(temp,x,npts,1)
      if(mod(npts,2) .eq. 1) then
         npt = npts/2 + 1
         xmed = temp(npt)
      else
         npt  = npts/2
         x1   = temp(npt)
         x2   = temp(npt+1)
         xmed = 0.5d0*(x1+x2)
      end if
      return
      end
c-----------------------------------------------------------------------
c
c     front end for real*4
c
      subroutine fmedian(npts, fx, fxmed)
      implicit none
      double precision x, xmed
      real fx, fxmed
      integer npts, i, lll
      parameter(lll=50000)
c      parameter(lll=67108864)
      dimension x(lll), fx(lll)
      do i = 1, npts
         x(i) = dble(x(i))
      end do
      call median(npts, x, xmed)
      fxmed = real(xmed)
      return
      end
