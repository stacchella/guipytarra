      subroutine read_gclf(mag, counts, npts, nnn)
      implicit none
      double precision am1, am2, rms1, rms2, mms1,mms2, bms1, bms2, 
     *     red,err_red, mid,err_mid, blue, err_blue, mag, total,
     *     counts
      integer i, nnn, npts
c
      dimension mag(nnn), counts(nnn)
c     
c      open(1,file='ngc2808_lf.dat')
c      read(1,*)
c      total = 0.d0
      npts = 0 
c      do i = 1, 10
c         read(1,*) am1, am2, rms1, rms2, mms1,mms2, bms1, bms2, 
c     *        red,err_red, mid,err_mid, blue, err_blue
c         npts = npts + 1
c         mag(npts) = (am1+ am2)/2.d0 -16.125d0
c         counts(npts) = red
c         total = total + red
c      end do
c      close(1)
cc
c      do i = 1, npts
c         counts(i) = counts(i) /total
c         print *, mag(i), counts(i)
c      end do
c
      open(2,file='ngc2419_lf.dat')
      read(2,*)
      read(2,*)
      total = 0.d0
      do i = 1, 20
         read(2, *,err=10,end=100) am1, rms1
         go to 20
 10      print *, am1, rms1
         stop
 20      npts = npts + 1
         mag(npts) = am1
         red = 10.d0**rms1
         counts(npts) = red
         total     = total+ red
      end do
 100  close(2)
      do i = 1, npts
         counts(i) = counts(i) /total
c         print *, mag(i), counts(i)
      end do
c
      return
      end
