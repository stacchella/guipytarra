c
c-----------------------------------------------------------------------
c
      subroutine read_filter_parameters(nf, verbose)
c
      implicit none
c
c     are all of these really needed ? One can use the actual
c     throughput curves:
c
      double precision effective_wl_nircam, width_nircam, 
     *     system_transmission
      double precision filters, filtpars
c
      integer nnn, mmm, lll, npar, nfilters, i, nf,l, ll, len, nf_used,
     &     verbose
c
      character filterid*20, tempfile*80
c
      parameter (nnn=25000,mmm=200,lll=1000,npar=30, nfilters = 54)
c
c     filter parameters, average responses from filters, and qe and mirror
c     averaged over filter.
c
      dimension effective_wl_nircam(nfilters), width_nircam(nfilters), 
     *     system_transmission(nfilters)
c
      dimension filters(mmm,nnn),filtpars(mmm,npar), filterid(mmm)
      common /filter/filters, filtpars, filterid
c     
      common /throughput/effective_wl_nircam, width_nircam,
     *     system_transmission
c
c     read filter parameters
c
c      open(1,file='calibrated_filters.dat')
      open(1,file='nircam_calib.list')
      nf = 0
      do i = 1, nfilters
         read(1,110,end=200, err=120) tempfile
 110     format(a80)
         go to 140
 120     print 130, i, tempfile
 130     format('problem with ',i3,2x,a80)
         stop
 140     nf = nf + 1
         open(3,file= tempfile)
         read(3,290) len, filterid(nf)
 290     format(i6,1x,a20)
c         print 290,  len, filters(nf)
         read(3,300) (filtpars(i,l),l=1, 27)
 300     format(10(1x,e21.10))
c 301     print 300,  (filtpars(nf,l),l=1, 27)
c     These are the parameters saved for filters. These are all included:
c     1. telescope reflectivity
c     2. dichroic reflectance/throughput
c     3. filter throughput (which combines 2 filters for filters on pupil Wheel)
c     4. NIRCam optics throughput
c     5. detector QE
c
c     filtpars(nf,1)  = wl(1)   - dwl
c     filtpars(nf,2)  = wl(nwl) + dwl
c     filtpars(nf,3)  = vega0   ! zero-point relative to vega
c     filtpars(nf,4)  = ab0     ! AB zero-point 
c     filtpars(nf,5)  = wl_mean ! = effective_wl for FIS, Rieke
c     filtpars(nf,6)  = fwhm    ! full-width half maximum
c     filtpars(nf,7)  = effective_nu
c     filtpars(nf,8)  = response ! average system response (filter+mirror+qe+op)
c     filtpars(nf,9)  = nominal_reach ! nominal WL according to Reach
c     filtpars(nf,10) = fnu_zp ! Flux in Jy for zero magnitudes
c     filtpars(nf,11) = dwl ! saved by rebin.f
c     filtpars(nf,12) = nominal   ! nominal WL according to GHR
c     filtpars(nf,13) = zmag    ! magnitude zero-point including offset (vega0)
c     filtpars(nf,14) = flambda_zp ! flambda corresponding to mag=0.0
c     filtpars(nf,15) = flambda_zp ! effective flux 
c     filtpars(nf,16) = bandwidth ! integral of normalised response
c     filtpars(nf,17) = response_normalised ! normalised response
c     filtpars(nf,18) = pivot        ! pivot wavelength
c     filtpars(nf,19) = effective_wl ! effective wavelength based on Vega
c     filtpars(nf,20) = lambda_iso ! isophotal wavelength for Vega
c     filtpars(nf,21) = half_power_l ! half-power wavelength (blue)
c     filtpars(nf,22) = half_power_r ! half-power wavelength (red)
c     filtpars(nf,23) = effective_response ! Hilbert & Stansberry definition 
c     filtpars(nf,24) = st_mag
c     filtpars(nf,25) = ab_mag    ! ab mag   corresponding to 1 e-/sec
c     filtpars(nf,26) = f_lambda  ! f_lambda corresponding to 1 e-/sec
c     filtpars(nf,27) = wave_resp ! um**2 to convert from f_lambda to electron f
c          
         effective_wl_nircam(nf) = filtpars(i, 5)
         width_nircam(nf)        = filtpars(i, 16)
         system_transmission(nf) = filtpars(i,8)
         if(verbose.gt.1) then
            print 310, filterid(nf),filtpars(i, 5),
     &           filtpars(i, 16), filtpars(i,8)
 310     format(a20,3(3x,f10.6),' read_filter_parameters')
         end if
c     
c     read wavelength and sensitivity curve
c
         do ll =1 , len
            read(3,300) filters(nf,ll)
         end do
         close(3)
      end do
 200  close(1)
c      stop
      return
      end
