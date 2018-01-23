      subroutine read_calibrated_filters(nf)
      implicit double precision (a-h,o-z)
      double precision nominal
      character filterid*20, filename*80
      parameter (nnn=25000,mmm=200,lll=1000,npar=30)
      dimension wlsed(nnn), sed(nnn), deriv(nnn)
      dimension filters(mmm,nnn),filtpars(mmm,npar), filterid(mmm)
      common /filter/filters, filtpars, filterid
c
c      filtpars(nf,1)  = wl(1)   - dwl
c      filtpars(nf,2)  = wl(nwl) + dwl
c      filtpars(nf,3)  = vega0   ! zero-point relative to vega
c      filtpars(nf,4)  = ab0     ! AB zero-point 
c      filtpars(nf,5)  = effective ! effective_lam
c      filtpars(nf,6)  = fwhm
c      filtpars(nf,7)  = effective_nu
c      filtpars(nf,8)  = response ! average system response (filter+mirror+qe+op)
c      filtpars(nf,9)  = nominal_reach ! nominal WL according to Reach
c      filtpars(nf,10) = fnu_zp ! Flux in Jy for zero magnitudes
c      filtpars(nf,11) = dwl
c      filtpars(nf,12) = nominal   ! nominal WL according to GHR
c      filtpars(nf,13) = zmag    ! magnitude zero-point including offset (vega0)
c      filtpars(nf,14) = flambda_zp ! flambda corresponding to mag=0.0
c      filtpars(nf,15) = flambda_zp ! effective flux 
c      filtpars(nf,16) = bandwidth ! integral of normalised response
c      filtpars(nf,17) = response_normalised ! normalised response
c      filtpars(nf,18) = pivot        ! pivot wavelength
c      filtpars(nf,19) = effective_wl ! effective wavelength based on Vega
c      filtpars(nf,20) = lambda_iso ! isophotal wavelength for Vega
c      filtpars(nf,21) = half_power_l ! half-power wavelength (blue)
c      filtpars(nf,22) = half_power_r ! half-power wavelength (red)
c      filtpars(nf,23) = effective_response ! Hilbert&Stansberry definition 
c      filtpars(nf,24) = st_mag
c      filtpars(nf,25) = ab_mag
c      filtpars(nf,26) = f_lambda  ! f_lambda corresponding to 1 e-/sec
c      filtpars(nf,27) = wave_resp ! um**2 to convert from f_lambda to electron flux

      open(1,file='calibrated_filters.dat')
      nf = 0
      do 1000 k = 1, mmm
         read(1,10,end=1010,err=20) filename
 10      format(a80)
         go to 30
 20      print 10, filename
 30      continue
         nf = nf + 1
         open(3,file=filename)
         read(3,290) len, filterid(nf)
 290     format(i6,1x,a20)
c         print 290,  len, filterid(nf)
         read(3,300) (filtpars(nf,l),l=1, 27)
 300     format(10(1x,e21.12))
c 300     format(10(1x,e16.10))
c         print 300,  (filtpars(nf,l),l=1, 26)
         do ll =1 , len
c     read wavelength and sensitivity curve
            read(3,300) filters(nf,ll)
         end do
         close(3)
 1000 continue
 1010 close(1)
      return
      end

