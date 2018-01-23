c      double precision zs, wlobs, igm, dwl, matrix, wl_rest
c      integer nsed,i, j
c      dimension wl_rest(2000),wlobs(2000), igm(2000), matrix(2000,10)
cc
c      nsed = 1301
c      dwl = (1400.d0-100.d0)/(nsed-1)
c      do i = 1, nsed
c         wl_rest(i) = 100. + (i-1)*dwl
c      end do
c      do j = 1, 7
c         zs = j
c         do i = 1, nsed
c            wlobs(i) = wl_rest(i) * (1.d0+zs)
c         end do
c         call igm_inoue_2014(zs, wlobs, igm, nsed)
c         do i = 1, nsed
c            matrix(i,j) = igm(i)
c         end do
c      end do
c      open(1,file='lixo')
c      do i = 1, nsed
c         write(1,100) wl_rest(i), (matrix(i,j),j=1,7)
c 100     format(10(1x,e12.6))
cc         write(*,*) wlobs(i), igm(i)
c      end do
c      close(1)
c      stop
c      end
c
c-----------------------------------------------------------------------
c
      subroutine igm_convolve_inoue(zs, wlobs, fluxin, fluxout, nwl)
      implicit none
      double precision zs, wlobs, fluxin, fluxout, igm
      integer nwl, i
      dimension wlobs(nwl), fluxin(nwl), fluxout(nwl), igm(nwl)
c
c      print *,'igm_convolve_inoue: zs, nwl ', zs, nwl
      call igm_inoue_2014(zs, wlobs, igm, nwl)
      do i = 1, nwl
         fluxout(i) = fluxin(i) * igm(i)
      end do
      return
      end
c
c-----------------------------------------------------------------------
c
c     Calculate the IGM attenuation using the model of
c     Inoue, Shimizu, Iwata & Tanaka 2014, MNRAS, 442, 1805.
c     The coefficients come from f90 code by A. Inoue, and part
c     of this F77 code is based on it.
c
c     The inputs are the source redshift, an array with the observed
c     wavelengths, returning an array with IGM attenuation.
c
c     cnaw@as.arizona.edu
c     2015-10-15
c
      subroutine igm_inoue_2014(zs, wlobs, igm, nwl)
      implicit none
      double precision wl, alaf1, alaf2, alaf3, dla1, dla2
      double precision lyman_limit, tau_laf, tau_dla, tau_laf_lc,
     &     tau_dla_lc, zplus1, wl_lyl, temp
      double precision zs, wlobs, igm
      integer nwl
      integer nl, i, j, jnd, init
c
      parameter(nl=39)
c
      dimension wlobs(nwl), igm(nwl)
      dimension wl(nl), alaf1(nl), alaf2(nl), alaf3(nl),
     &     dla1(nl),dla2(nl), jnd(nl)
      save init, wl, alaf1, alaf2, alaf3, dla1, dla2
      data init/1/
c
      lyman_limit = 911.8d-4    ! microns
c
c     read parameters
c
      if(init .eq.1) then 
         open(1,file='inoue_2014_mnras.dat')
         read(1,*)
         do i = 1, nl
            read(1, *) jnd(i), wl(i), alaf1(i), alaf2(i), alaf3(i),
     &           dla1(i),dla2(i)
            wl(i) = wl(i)*1.d-04 ! convert into microns
c            print 100, jnd(i), wl(i), alaf1(i), alaf2(i), alaf3(i),
c     &           dla1(i),dla2(i)
 100        format(i3,2x, f10.3,6(1x,e16.4))
         end do
         close(1)
         init = 0
      end if
c
      zplus1 = 1.d0+zs
      do i = 1, nwl
         wl_lyl  = (wlobs(i)/lyman_limit)
         tau_laf    = 0.0d0
         tau_dla    = 0.0d0
         tau_laf_lc = 10000.d0
         tau_dla_lc = 10000.d0
c
c     loop over Lyman series lines
c
         do j = 1, nl
c
c     contribution due to Lyman series absorption in Lyman-alpha forest
c
            if(wlobs(i).gt.wl(j) .and. wlobs(i).lt.wl(j)*zplus1) then
               if(wlobs(i).lt.wl(j)*2.2d0) then
                  tau_laf = tau_laf + alaf1(j) *(wlobs(i)/wl(j))**1.2d0
               end if
c     
               if(wlobs(i).ge.wl(j)*2.2d0.and.
     &              wlobs(i).lt.wl(j)*5.7d0) then
                  tau_laf = tau_laf + alaf2(j) *(wlobs(i)/wl(j))**3.7d0
               end if
c
               if(wlobs(i).ge.wl(j)*5.7d0) then
                  tau_laf = tau_laf + alaf3(j) *(wlobs(i)/wl(j))**5.5d0
               end if
            end if
c
c     contribution due to DLA systems
c
            if(wlobs(i).gt.wl(1) .and. wlobs(i).lt.wl(j)*zplus1) then
               if(wlobs(i).lt.wl(j)*3.0d0) then
                  tau_dla = tau_dla + dla1(j) *(wlobs(i)/wl(j))**2
               end if
c     
               if(wlobs(i).ge.wl(j)*3.0d0) then
                  tau_dla = tau_dla + dla2(j) *(wlobs(i)/wl(j))**3
               end if
            end if
         end do 
c
c     Continuum absorption due to Lyman-alpha forest
c     
         if(wlobs(i) .ge. lyman_limit*zplus1) then
            tau_laf_lc = 0.0d0
         else
c     
c     z < 1.2
c     
            if(zs.lt.1.2d0) then
               temp = wl_lyl**1.2d0 -
     &              (wl_lyl**2.1d0) / (zplus1**0.9)
               tau_laf_lc = 0.3248d0 * temp
            end if
c     
c     1.2 <= z < 4.7
c     
            if(zs.ge.1.2d0 .and. zs.lt.4.7d0) then
               if(wlobs(i) .lt. lyman_limit*2.2d0) then
                  tau_laf_lc = 
     &                 2.545d-02 * zplus1**1.6d0 * wl_lyl**2.1d0
     &                 + 0.3248d0 * wl_lyl**1.2d0 
     &                 - 0.2496d0 * wl_lyl**2.1d0 
               else
                  temp = zplus1**1.6d0 * wl_lyl**2.1d0
     &                 -  wl_lyl**3.7d0
                  tau_laf_lc = 2.545d-02 * temp
               end if
            end if
c     
c     4.7 <= z 
c     
            if(zs.ge. 4.7d0) then
               if(wlobs(i) .lt. lyman_limit*2.2d0) then
                  tau_laf_lc = 5.221d-04*(zplus1**3.4d0)*wl_lyl**2.1d0
     &                 + 0.325d0 * wl_lyl**1.2d0 
     &                 - 3.14d-02 * wl_lyl**2.1d0
               end if
c     
               if(wlobs(i) .ge. 2.2d0*lyman_limit .and. 
     &              wlobs(i).lt. 5.7d0*lyman_limit) then
                  tau_laf_lc = 5.221d-04*(zplus1**3.4d0)*wl_lyl**2.1d0
     &                 + 0.2182d0  * (wl_lyl**2.1d0 )
     &                 - 2.545d-02 * (wl_lyl**3.7d0)
               end if
            end if
c     
            if(wlobs(i) .ge. 5.7d0 * lyman_limit .and. 
     &           wlobs(i).lt. lyman_limit*zplus1) then
               temp = zplus1**3.4d0 * wl_lyl**2.1d0
     &              - wl_lyl**5.5d0
               tau_laf_lc  = temp * 5.221d-04
            end if
c     
         end if
c     
c     DLA continuum component
c     
         if(wlobs(i).gt.lyman_limit*zplus1) then
            tau_dla_lc =0.0d0
         else
            if(zs.lt.2.0d0) then
               tau_dla_lc = 0.2113d0 * zplus1*zplus1 
     &              - 7.661d-2 * (zplus1**2.3d0) / (wl_lyl**0.3d0)
     &              - 0.1347d0*wl_lyl*wl_lyl
            else
               if(wlobs(i) .ge. lyman_limit*3.0d0) then
                  tau_dla_lc  = 4.696d-02*zplus1**3 
     &                 - 1.779d-02 * (zplus1**3.3d0) / (wl_lyl**0.3)
     &                 - 2.916d-02 * wl_lyl**3
               else
                  tau_dla_lc = 0.634d0 + 4.696d-02*zplus1**3 
     &                 - 1.779d-02*(zplus1**3.3d0)/(wl_lyl**0.3)
     &                 - 0.1347d0*(wl_lyl**2)-0.2905d0/(wl_lyl**0.3)
               end if
            end if
         end if
c     print *,i, wlobs(i), tau_laf, tau_dla, tau_laf_lc,tau_dla_lc
      igm(i) = tau_laf + tau_dla + tau_laf_lc + tau_dla_lc
c         igm(i ) = tau_laf
c         igm(i) =  tau_dla
c         igm(i) = tau_laf_lc
c         igm(i) = tau_dla_lc
c     igm(i) = 1
         igm(i) = dexp(-igm(i))
c     print *,zs, wlobs(i)/zplus1, igm(i)
      end do
      return
      end
c
