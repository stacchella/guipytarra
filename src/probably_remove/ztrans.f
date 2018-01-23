c     
c-----------------------------------------------------------------------
c
c     read the transmission curves of the Universe calculated by
c     Piero Madau ranging from 0 <= z <= 5
c     Note that wavelengths are converted into microns
c
      subroutine read_ztrans
      implicit double precision (a-h,o-z)
      parameter (lll=100)
      dimension zz_ztrans(lll),wl_ztrans(lll,lll), ztrans(lll,lll)
      dimension wl(lll), tau(lll), exp_tau(lll)
      common /igm/ zz_ztrans, wl_ztrans, ztrans, n_ztrans, nwl_ztrans
      character file*80
      n_ztrans   = 51
      nwl_ztrans = 82 
      do k = 1, n_ztrans
         zz = (k-1.d0)/10.d0
         zz_ztrans(k) = zz
         write(file, 40) zz
 40      format('/home/cnaw/igm/trans_z',f3.1)
c         print 80, k, file
 80      format(i3,2x,a50)
         open(1,file=file)
         do i = 1, 3
            read(1,*)
         end do
c
         do i = 1, nwl_ztrans
            read(1,*) wll, t, expt
            wl(i) = wll/1.d04   ! convert into microns
            tau(i) = t
            exp_tau(i) = expt
         end do
         close(1)
c
         do i = 1, nwl_ztrans
            wl_ztrans(k,i) = wl(nwl_ztrans-i+1)
            ztrans(k,i)    = exp_tau(nwl_ztrans-i+1)
         end do
      end do
      open(2,file ='lixo2')
      do i = 1, nwl_ztrans
         write(2, 100) (wl_ztrans(k,i)*1.d04/(1.d0+((k-1.d0)/10.d0)),
     *        ztrans(k,i),k =1, n_ztrans)
 100     format(100(f9.3,1x,f9.7))
      end do
      close(2)
      return
      end
c
c-----------------------------------------------------------------------
c
c     Follow the Dickinson/Papovich procedure of convolving the IGM
c     with an SED
c
      subroutine igm_convolve(zz, wlsed, sed, sed_igm, nsed)
      implicit double precision (a-h,o-z)
      integer debug
      parameter (lll=100,nnn=25000)
      dimension wlsed(nnn), sed(nnn), sed_igm(nnn)
      dimension indices(nnn), iend(nnn)
      dimension x(nnn), y(nnn), deriv(nnn), trans(nnn)
      dimension spline_a(nnn), spline_b(nnn), spline_c(nnn)
       dimension zz_ztrans(lll),wl_ztrans(lll,lll), ztrans(lll,lll)
      common /igm/ zz_ztrans, wl_ztrans, ztrans, n_ztrans, nwl_ztrans
c
      debug = 0
c
c     Find closest redshift to input value. According to Stiavelli,
c     Giavalisco & Carollo (STScI Instrument Science Report WFC3 2000-19
c     for z > 5, shifting the z = 5 by (1+z) should work. However, the
c     results differ from Avery Meiksin's, which cover to z=7.
c
      index = idnint(zz*10.d0) + 1
      if(index .gt. 51) index = 51
      if(debug .eq.1) then 
         print *, 'ztrans: igm_convolve: zz, index',zz, index
      end if
c
c     This will be the transmission curve
c
      do i = 1, nnn
         trans(i)  = 1.0d0
      end do
c
c     The trick is to find the discontinuities in the IGM transmission 
c     curves. This is done following Mark D's SM script
c
      j  = 1
      indices(j) = 1
      do i = 2, nwl_ztrans
         dwl = wl_ztrans(index,i) - wl_ztrans(index,i-1)
         if(dwl.eq.0.0d0) then
            j = j + 1
            iend(j-1)  = i-1
            indices(j) = i
            if(debug .eq.1) then 
               print *, ' break at ',i, wl_ztrans(index,i)
            end if
         end if
      end do
      iend(j) = nwl_ztrans
c     
c     For each each piece, use spline interpolations to convolve
c     the IGM transmission with spectrum
c
      displacement = 1.0d0
      atten        = 1.0d0
      if(zz.gt.5.d0) then
         displacement = (1.d0+zz)/(1.d0+5.d0)
         if(zz.le.5.5)  then
            atten        = 0.00058d0 * (1.d0 + zz)**3.2d0
         else
            atten        = 0.00058d0 * (1.d0 + zz)**4.0d0
         end if
         atten        = dexp(-atten)
      end if
      nparts = j
      do k = 1, nparts
         nn = iend(k) - indices(k) +1
c
c     copy the continuous segment so a spline may be calculated
c
         i = 0
         do l = indices(k), iend(k)
            i = i + 1
            x(i) = wl_ztrans(index,l) * displacement
            y(i) = ztrans(index,l) *atten
         end do
c         call spline(x, y, nn, 0.0d0, 0.0d0, deriv)
         call spline_interp(nn, x, y, 
     *        spline_a, spline_b, spline_c)
c     
c     use the wavelengths of the galaxy spectrum to find what is the
c     transmission due to the IGM
c     
         do i = 1, nsed
            if(wlsed(i).ge. x(1) .and. wlsed(i).le.x(nn)) then
c               call splint(x, y, deriv, nn, wlsed(i), trans(i))
               trans(i) =  seval(nwl_trans, wlsed(i), x, y,
     *        spline_a, spline_b, spline_c)
               if(debug .eq.1) then 
                  print *, wlsed(i), trans(i)
               end if
            end if
         end do
      end do
c
c     set transmission to  1 for wavelengths bluer than  wl_ztrans(index,1)
c     and wavelengths > (1+z) * 1216
c
      wlmin = (1.d0 + zz) *  912.d0/1.d04
      wlmax = (1.d0 + zz) * 1216.d0/1.d04
      do i = 1, nsed
         if(zz.lt.3.d0) then
            if(wlsed(i) .lt. wl_ztrans(index,1)) trans(i) = 1.d0
         else
            if(wlsed(i) .lt.  wlmin)  trans(i) = 0.d0
         end if
         if(wlsed(i) .gt. wlmax) trans(i) = 1.d0
c         print *, i, wlsed(i), wlsed(i)/(1.d0 + zz), trans(i)
      end do
c     
c     now convolve with SED
c
      do i = 1, nsed
         sed_igm(i) = trans(i) * sed(i)
         if(debug .eq.1) then 
            wl_rest = wlsed(i) * (1.d0+zz)
            ratio   = sed(i)/sed_igm(i)
c            print *, wlsed(i), wl_rest, sed(i), sed_igm(i),ratio
         end if
      end do
c
      return
      end
