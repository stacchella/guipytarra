C  By Jason Pinkney  5/94, 5/3/95, 7/23/97 (UNLV), 10/29/97 (UNLV)
c
c     Although this came from J. Pinkney's cluster analysis software,
c     this code was writen by  Beers, Flynn & Gebhardt 1990 AJ as
c     part of the Robust Statistics package:
c     http://www.pa.msu.edu/ftp/pub/beers/posts/rostat/rostat.f
c
c     Modified so that the input data are NOT sorted AND removed
c     Numerical recipes routines.
c     
c     2013-06-11 CNAW
c      
        SUBROUTINE XBIWT (XDATAIN,N,XLBIWT,XSBIWT,XLBIWT1,XSBIWT1,DEBUG)
c-----------------------------------------------------------------------------

c--- The subroutine XBIWT provides an estimator of the location and
c    scale of the data set XDATA.  The scale uses the Biweight function
c    in the general formula of "A-estimators." This formula is given
c    on page of 416 in UREDA (formula 4). The BIWEIGHT scale estimate
c    is returned as the value XSBIWT. The BIWEIGHT function is given
c    by:
c
c                                  u((1-u*u)**2)     abs(u) <= 1
c                         f(u) =
c                                  0                 abs(u) >  1
c
c    where u is defined by
c
c                         u = (XDATA(I) - M) / c*MAD  .
c
c    M, MAD, and c are the median, the median absolute deviation from
c    the median, and the tuning constant respectively. The tuning
c    constant is a parameter which is chosen depending on the sample
c    size and the specific function being used for the scale estimate.
c    (See page 417 in UREDA).  Here we take c = 9.0.
c
c--- The biweght location is found using the formula:
c
c                         T = M + (sums)
c
c                         where M is the sample median and sums are
c                         as given on page 421 in UREDA
c
c                         the tuning constant c is set to 6.0 for calculation
c                         of the location as reccommended by Tukey ()
c
c--- NOTE that the biweight is meant to be an iterated estimator, but one
c    commonly only takes the first step of the iteration.  Here we report
c    both the one-step estimators (XLBIWT1, XSBIWT1) and the preferred
c    fully iterated versions (XLBIWT, XSBIWT).
c
c****************************************************************************
       implicit none
       double precision xdata,u1,u2,c1,c2,xlbiwt1,xm,xmadm
       double precision zero,d6,d1,d5,d9,s1,s2,s3,s4
       double precision  xlb,xsb,xlbiwt,xsbiwt,xsbiwt1,xmm
       double precision temp, xdatain
       integer i,n,j
       integer nnn, debug
       parameter (nnn=1048576)
c        parameter (nnn=67108864)
       dimension xdata(nnn),u1(nnn),u2(nnn),xlb(11),xsb(11)
        dimension xdatain(nnn), temp(nnn)
        data zero,d6,d1,d5,d9/0.0,6.0,1.0,5.0,9.0/
c
c     copy input data 
c
        do i = 1, n
           xdata(i) = xdatain(i)
        end do
c---    sort the data and find the median
c 
c        CALL MDIAN1(XDATA,N,XM)
        do i = 1, n
           temp(i) = xdata(i)
        end do
        if(debug.eq.1) print *,'xbiwt enter dsort 1',n
        call dsort(temp, xdata,n,1)
        if(debug.eq.1) print *,'xbiwt enter median '
        call median(n, xdata, xm)
        if(debug.eq.1) print *,'xbiwt exit median:xm ',xm
c
c---    call xmad to find the median absolute deviation
 
        if(debug.eq.1) print *,'xbiwt enter XMAD'
        CALL XMAD(XDATA,N,XM,XMADM)
        if(debug.eq.1) print *,'xbiwt exit XMAD: xmadm',xmadm
 
c---    must choose value of the tuning constant "c"
c       here c = 6.0 for the location estimator and
c       9.0 for the scale estimator
 
        c1 = d6
        c2 = d9
 
        if (xmadm.le..0001d0) then
        xlbiwt=xm
        xlbiwt1=xm
        xsbiwt=xmadm
        xsbiwt1=xmadm
        goto 20
        endif
 
        do 11 i = 1,n
        u1(i) = (xdata(i) - xm)/(c1*xmadm)
        u2(i) = (xdata(i) - xm)/(c2*xmadm)
11      continue
 
        s1 = zero
        s2 = zero
        s3 = zero
        s4 = zero
 
        do 12 i = 1,n
        if (dabs(u2(i)) .lt. d1) then
            s1 = s1+(((xdata(i)-xm)**2)*(d1-(u2(i)*u2(i)))**4)
            s2 = s2+((d1-u2(i)*u2(i))*(d1-(d5*u2(i)*u2(i))))
        endif  
        if (dabs(u1(i)) .lt. d1) then
            s3 = s3+(xdata(i)-xm)*(d1-u1(i)*u1(i))**2
            s4 = s4+(d1-u1(i)*u1(i))**2
        endif  
12      continue
 
c--- here are the one-step estimators
 
        xlbiwt1 = xm+s3/s4
        xsbiwt1 = dble(n)/(dble(n-1))**0.5*s1**0.5/dabs(s2)
 
c--- now obtain the fully-iterated versions
 
c--- solve for new estimates of u1 and u2

        xlb(1) = xlbiwt1
        xsb(1) = xsbiwt1
 
        do 15 j = 2,11   !assuming 10 iterations is sufficient
 
        xmm = xlb(j-1)
        do 13 i = 1,n
        u1(i) = (xdata(i) - xmm)/(c1*xmadm)
        u2(i) = (xdata(i) - xmm)/(c2*xmadm)
13      continue
 
        s1 = zero
        s2 = zero
        s3 = zero
        s4 = zero
 
        do 14 i = 1,n
        if (dabs(u2(i)) .lt. d1) then
            s1 = s1+(((xdata(i)-xmm)**2)*(d1-(u2(i)*u2(i)))**4)
            s2 = s2+((d1-u2(i)*u2(i))*(d1-(d5*u2(i)*u2(i))))
        endif  
        if (dabs(u1(i)) .lt. d1) then
            s3 = s3+(xdata(i)-xmm)*(d1-u1(i)*u1(i))**2
            s4 = s4+(d1-u1(i)*u1(i))**2
        endif  
14      continue
 
        xlb(j) = xlb(j-1)+s3/s4
        xsb(j) = dble(n)/(dble(n-1))**0.5*s1**0.5/dabs(s2)
 
15      continue
 
        xlbiwt = xlb(11)
        xsbiwt = xsb(11)
 
20      continue
        return
        end
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------------
        SUBROUTINE XMAD (XDATA,N,XMED,XMADM)
c-----------------------------------------------------------------------------

c--- The XMAD subroutine calculates the Median Absolute Deviation from
c    the sample median. The median, M , is subtracted from each
c    ORDERED statistic and then the absolute value is taken. This new
c    set of of statistics is then resorted so that they are ORDERED
c    statistics. The MAD is then defined to be the median of this
c    new set of statistics and is returned as XMADM. The MAD can
c    be defined:   
c
c                   XMADM = median{ abs(x(i) - M) }
c
c    where the x(i) are the values passed in the array XDATA, and
c    the median, M, is passed in the array XLETTER. The set of stats
c    in the brackets is assumed to be resorted. For more information
c    see page 408 in UREDA.
c
c
c****************************************************************************
       implicit none 
       double precision xdata2,xdata, temp
       double precision dhalf,xmadm,xmed
       integer n,i,i1,i2,n1,n2
       integer nnn
       parameter (nnn = 1048576)
       dimension  xdata2(nnn),xdata(nnn), temp(nnn)
       data dhalf,n1,n2/0.5,1,2/

       do i = 1,n
          xdata2(i)   = dabs(xdata(i) - xmed)
          temp(i)    = xdata2(i)
       end do
        
       call dsort(xdata2,temp, n,1)
       
       if (dble(n)/dble(n2) - int(n/n2) .eq. 0) then
          i1 = n/n2
          i2 = n/n2 + n1
          xmadm = dhalf*(xdata2(i1) + xdata2(i2))
       else       
          i1 = int(n/n2) + n1
          xmadm = xdata2(i1)
       endif      
       return
       end

