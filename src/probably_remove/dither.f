c     Create data ramps for an SCA, a module or both modules at
c     a given dither offset. At each dither position loop through
c     the set of filters (obviously this is what will _not_ happen
c     in orbit !)
c
c     cnaw 2015-01-27
c     Steward Observatory, University of Arizona
c     
      subroutine dither(idither, ra_dithered, dec_dithered,
     *     sca, module, brain_dead_test,
     *     x0, y0, pa_v3, osim_scale, bkg,
     *     include_ktc, include_dark, include_readnoise, 
     &     include_reference,
     *     include_1_over_f, include_latents, include_non_linear,
     *     include_cr, cr_mode, include_bg,
     *     include_stars, nstars, include_galaxies, ngal,
     *     bitpix, ngroups, nframe, nskip, tframe, tgroup, object,
     *     subarray, colcornr, rowcornr, naxis1, naxis2,
     *     use_filter, filterid, nf, nf_used,
     *     cat_filter, ncat_f,
     *     psf_file, verbose)
      
      implicit none
c
c     variables are defined as double precision
c     arrays used for images are either real*4 or integer
c
      double precision x0, y0, pa_v3, ra_dithered, dec_dithered, 
     *     osim_scale
      double precision mirror_area,  integration_time,  gain,
     *     time_since_previous, decay_rate, read_noise,
     *     dark_mean, dark_sigma, ktc, voltage_offset

      double precision effective_wl_nircam, width_nircam, response
      double precision wavelength, bandwidth, system_transmission,
     *     photplam, photflam, stmag, abmag
      double precision bkg, background, scale
c
      real tframe, tgroup
c
      integer bitpix, ngroups, nframe, nskip, cr_mode
      integer include_ktc, include_bg, include_cr, include_dark,
     *     include_latents, include_readnoise, include_non_linear,
     *     include_stars, include_galaxies, brain_dead_test,
     *     include_1_over_f, include_reference
      integer i, j, k, loop, i1, i2, nf, idither, verbose, sca_id, sca
      integer nnn, nfilters, mmm, lll, npar, ngal, nstars, use_filter,
     *     filter_index, indx , nf_used, cat_filter, ncat_f, icat_f
c
      character filter_id*5,module*4, psf_file*120, object*20,
     *     temp*20, filterid*20

      logical subarray
      integer  colcornr, rowcornr, naxis1, naxis2
c
      parameter (nnn=2048, nfilters=54)
      parameter (lll=25000,mmm=200,npar=30)
c
      dimension psf_file(nfilters), cat_filter(nfilters)
c
c     SCA parameters
c
      dimension dark_mean(10), dark_sigma(10), gain(10),
     *     read_noise(10)
c      
      dimension effective_wl_nircam(nfilters), width_nircam(nfilters), 
     *     response(nfilters), use_filter(nfilters), bkg(nfilters)
c      
      dimension filterid(mmm)
 
      common /parameters/ mirror_area, integration_time, gain,
     *     decay_rate, time_since_previous, read_noise, 
     *     dark_mean, dark_sigma, ktc, voltage_offset
      common /throughput/effective_wl_nircam, width_nircam, response
c
c     Parameters
c
c     select SCAs to simulate
c
      if(sca .lt. 481) then 
         if(module .eq.'none') then
            i1 =  1
            i2 =  1
         end if

         if(module .eq.'all') then
            i1 =  1
            i2 = 10
         end if
c     
         if(module .eq. 'modA') then
            i1 =  1
            i2 =  5
         end if
c     
         if(module .eq. 'modB') then
            i1 =  6
            i2 = 10
         end if
c     
      else
         i1 = sca - 480
         i2 = i1
      end if
c
c     This is where parallel processing can be done
c
      do i = i1, i2
         sca_id = 480 + i
         scale               =   0.0317d0
         if(sca_id .eq. 485 .or. sca_id .eq.490) then
            scale = 0.0648d0
         end if
c     
c     run through filters, setting up filter-dependent parameters
c     such as system transmission, filter bandwidth, wavelength
c     
         do indx = 1, nf
            print *,'dither: indx ', indx
            if(nf_used.eq.54 ) then
               j= use_filter(indx)*2 -1
            else
               j = use_filter(indx)
            end if
            temp                = filterid(j)
c            filter_id           = temp(8:12)
            filter_id           = temp(1:5)
            wavelength          = effective_wl_nircam(j)
            bandwidth           = width_nircam(j)
            system_transmission = response(j) 
            background          = bkg(j)
            icat_f              = cat_filter(indx)
c            photplam            =
c            photflam            = 
c            stmag               =
c            abmag               =
c     
c     use appropriate filters for SW/LW
c
            if(scale.eq.0.0317d0.and.wavelength.gt.2.4d0) go to 999
            if(scale.eq.0.0648d0.and.wavelength.lt.2.4d0) go to 999
c     
            print 120, j, filter_id, sca_id, scale, wavelength,NF
 120        format(i3,2x,a5,2x,i3,2x,f9.4,2x,f9.4,2x,i3)
            print 130, j, psf_file(j)
 130        format('dither.f: ',i3,1x,a120)
c     
c     psf should be really constructed for SCA + filter
c            
            call sca_image(idither, ra_dithered, dec_dithered,
     *           sca_id, module, brain_dead_test,
     *           x0, y0, pa_v3, osim_scale, scale,
     *           include_ktc, include_dark, include_readnoise, 
     &           include_reference,
     *           include_1_over_f, include_latents, include_non_linear,
     *           include_cr, cr_mode, include_bg,
     *           include_stars, include_galaxies, nstars, ngal,
     *           bitpix, ngroups, nframe, nskip, tframe, tgroup, object,
     *           subarray, colcornr, rowcornr, naxis1, naxis2,
     *           filter_id, wavelength, bandwidth, system_transmission,
     *           photplam, photflam, stmag, abmag,
     *           background, icat_f,j, psf_file(j), verbose)
            go to 1000
 999        continue
            print *,' in dither.f'
            print *,'scale : ',scale,' wavelength: ',wavelength
            print *,'something is inconsistent; posible error:'
            print *,'read_filter_parameters not reading appropriate',
     &           ' list'
 1000       continue
         end do
      end do
c     
      return
      end
