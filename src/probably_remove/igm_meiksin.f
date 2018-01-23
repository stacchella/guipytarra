c     
c-----------------------------------------------------------------------
c
c     read the transmission curves of the Universe calculated by
c     Avery Meikin (2006) from 1.5 <= z <= 7
c     Note that wavelengths are converted into microns
c
      subroutine read_ztrans_meiksin
      implicit double precision (a-h,o-z)
      parameter (lll=8000,nz=15)
      dimension zz_ztrans(nz),wl_ztrans(lll,nz), ztrans(lll,nz)
      dimension wl(lll), tau(lll), exp_tau(lll)
      common /igm2/ zz_ztrans, wl_ztrans, ztrans, n_ztrans, nwl_ztrans
      character file*80
      n_ztrans   = 12
      nwl_ztrans = lll 
      open(1,file='./igm/meiksin_igm_transmission_table.dat')
      open(2,file='lixo1')
      read(1,*)
      do i = 1, n_ztrans 
         zz_ztrans(i) = 1.5d0 + dble(i-1) * 0.5d0
      end do
      do i = 1, nwl_ztrans
         read(1,*) wll, (ztrans(i,j),j=1,n_ztrans)
         do j = 1, n_ztrans
            zplus1 = zz_ztrans(j) +1.0d0
            wl_ztrans(i,j) = wll/zplus1/1.d04
         end do
         write(2,100) (wl_ztrans(i,j)*1.d04, ztrans(i,j),j=1,n_ztrans)
 100     format(20(f9.3,1x,f9.7))
      end do
      close(1)
      close(2)
c
      return
      end
c
c-----------------------------------------------------------------------
c
c     Follow the Dickinson/Papovich procedure of convolving the IGM
c     with an SED
c
      subroutine igm_convolve_meiksin(zz, wlsed, sed, sed_igm, nsed)
      implicit double precision (a-h,o-z)
      integer debug
      parameter (lll=8000,nz=15,nnn=25000)
      dimension zz_ztrans(nz),wl_ztrans(lll,nz), ztrans(lll,nz)
      dimension wl(lll), tau(lll), exp_tau(lll)
      dimension wlsed(nnn), sed(nnn), sed_igm(nnn)
      dimension spline_a(nnn), spline_b(nnn), spline_c(nnn)
      dimension indices(nnn), iend(nnn)
      dimension x(nnn), y(nnn), deriv(nnn), trans(nnn)
      common /igm2/ zz_ztrans, wl_ztrans, ztrans, n_ztrans, nwl_ztrans
      character file*80
c
      debug = 0
c
c     Find closest redshift to input value. According to Stiavelli,
c     Giavalisco & Carollo (STScI Instrument Science Report WFC3 2000-19
c     for z > 5, shifting the z =5 by (1+z) should work
c
      index = idnint((zz-1.5d0)/0.5d0) + 1
      if(debug .eq.1) then 
         print *, 'igm_meiksin: igm_convolve: zz, index',zz, index
      end if
c
c     if z > 7 make all flux bluer than 1216 equal to zero
c
      if(index .gt. 11)  then
         do i = 1, nsed
c     atten        = 0.00058d0 * (1.d0 + zz)**4.5d0
c     trans(i)     = dexp(-atten)
            wl_rest = (wlsed(i)/(1.d0+zz))*1.d04
            trans(i) = 1.0d0
            if(wl_rest.le.1216.d0) trans(i) = 0.0d0
         end do
      else
c
c     This will be the transmission curve
c
         do i = 1, nnn
            trans(i)  = 1.0d0
         end do
c     
c     The Meiksin curves are sampled at every Angstrom (observed) and
c     although there are discontinuities, each wavelength has a unique
c     transmission, different from the Madau version implemented by
c     Dickinson & Papovich. Thus, resample the observed spectrum.
c
c      k = 0 
         do i = 1, nwl_ztrans
            x(i) = wl_ztrans(i, index)
            y(i) = ztrans(i,index)
         end do
c     
c         call spline(x, y, nwl_ztrans, 0.0d0, 0.0d0, deriv)
c         call spline_interp(nwl_ztrans, x, y, 
c     *        spline_a, spline_b, spline_c)
c     
c     use the wavelengths of the galaxy spectrum to find what is the
c     transmission due to the IGM
c     
         do i = 1, nsed
            wl_rest = wlsed(i)/(1.d0+zz)
            if(wl_rest.lt. x(1).and.zz .le.3.5d0)  trans(i) = y(1)
            if(wl_rest.lt. x(1).and.zz .gt.3.5d0)  trans(i) = 0.0d0
            if(wl_rest.ge. x(1) .and. wl_rest.le.x(nwl_ztrans)) then
               call linear_interpolation(nwl_trans, x, y, wl_rest,
     *              trans(i), nnn)
c               trans(i) =  seval(nwl_trans, wl_rest, x, y,
c     *        spline_a, spline_b, spline_c)
c               call splint(x, y, deriv, nwl_ztrans, wl_rest, trans(i))
               if(trans(i) .gt. 1.d0) then
                  print *, zz, wl_rest, trans(i),x(1), y(1)
                  read(*,'(A)')
               end if
            end if
         end do
c
      end if
c     
c     now convolve with SED
c     
      do i = 1, nsed
         sed_igm(i) = trans(i) * sed(i)
c         if(debug .eq.1) then 
            wl_rest = wlsed(i) * (1.d0+zz)
c            print *, wlsed(i), wl_rest, sed(i), trans(i),sed_igm(i)
c         end if
      end do
c
      return
      end
