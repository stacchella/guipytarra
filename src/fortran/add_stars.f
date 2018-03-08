      subroutine add_stars(ra_dithered, dec_dithered, pa_degrees,
     *     xc, yc, osim_scale, sca_id, filter_index,
     *     seed, subarray, colcornr, rowcornr, naxis1, naxis2,
     *     wavelength, bandwidth,system_transmission,
     *     mirror_area, integration_time, 
     *     noiseless, psf_add, ipc_add, verbose)
c
c     Add "stars" to an image. The NIRCam footprint is centered at
c     ra_dithered, dec_dithered 
c     
      implicit none
      double precision ra_dithered, dec_dithered, xc, yc, pa_degrees,
     *     osim_scale
      double precision ra_stars, dec_stars, mag_stars
      double precision rra, ddec, rra_pix, ddec_pix
      double precision intensity, total_per_cycle, stellar_photons,
     *     xs, ys, xx, yy, xhit, yhit, x_osim, y_osim
      double precision ab_mag_to_photon_flux, zbqlnor
      double precision mirror_area, wavelength, bandwidth,
     *     system_transmission, integration_time
c
      real gain_image
c
      integer expected, zbqlpoi
      integer colcornr, rowcornr, naxis1, naxis2, filter_index, junk
      integer verbose, sca_id, nnn, max_stars, nfilters, ix, iy,
     *     i, j, nstars, ixmin, ixmax, iymin, iymax, seed, invert
c
      character subarray*8
      logical noiseless, psf_add, ipc_add
c
      parameter (max_stars=10000,nnn=2048,nfilters=54)
c
      dimension gain_image(nnn,nnn)
      dimension ra_stars(max_stars), dec_stars(max_stars), 
     *     mag_stars(max_stars,nfilters)

      common / gain_/ gain_image
      common /stars/ ra_stars, dec_stars, mag_stars, nstars

      if (verbose .gt. 1) then
         print *,'enter add_stars'
      end if
c
      if(subarray(1:4) .eq.'FULL') then
         ixmin = 5
         ixmax = naxis1- 4
         iymin = 5
         iymax = naxis2 - 4
      else
         ixmin = colcornr
         ixmax = ixmin + naxis1
         iymin = rowcornr
         iymax = iymin + naxis2
      end if
c
      xhit     = 0.0d0
      yhit     = 0.0d0
c
      if(sca_id.lt.481) then
         junk = 1
      else
         junk = sca_id
      end if
c
      print *, 'going to add ', nstars, ' stars to image!'
c
      do i = 1, nstars
c
c     find SCA coordinates for this object 
c
        print *, 'star RA, DEC: ', ra_stars(i), dec_stars(i)
c
        call ra_dec_to_sca(junk, 
     *        ra_dithered, dec_dithered, 
     *        ra_stars(i), dec_stars(i), pa_degrees, 
     *        xc, yc,  osim_scale, xs, ys)
        if(verbose.gt.1) then
            print *,'add_stars: ixmin, ixmax, iymin, iymax, xs, ys',
     &           ixmin, ixmax, iymin, iymax, xs, ys
        end if
c
c     calculate the number of photo-electrons per second
c     (system_transmission contains the quantum efficiency term)
c
        stellar_photons = 
     *        ab_mag_to_photon_flux(mag_stars(i, filter_index),
     *        mirror_area,  wavelength, bandwidth, system_transmission)
c         print *,stellar_photons, integration_time, filter_index,
c     *        mag_stars(i,filter_index), mirror_area
c     
c      Find expected number of photo-electrons
c     
        total_per_cycle = stellar_photons * integration_time
c     
c     for noiseless
c     
        if(noiseless .eqv. .true.) then
           expected = total_per_cycle !+ zbqlnor(0.0d0,0.5d0)
        else
           expected = zbqlpoi(total_per_cycle)
        end if
c     write(17, 10) ra_stars(i), dec_stars(i), xs, ys,
c     *           mag_stars(i,filter_index), stellar_photons,
c     *           total_per_cycle, expected
        if(verbose.gt.1) then
            print 10, ra_stars(i), dec_stars(i), xs, ys,
     *          mag_stars(i, filter_index), stellar_photons,
     *          total_per_cycle, expected
 10         format('add_stars ', 4(1x,f12.6), f8.3, 
     *           2(2x,f12.2),2x,i10)
        end if
        if(expected .gt.0 ) then
            do j = 1, expected
               if(psf_add .eqv. .true.) 
     &              call psf_convolve(seed, xhit, yhit)
c     
               ix = idnint(xs - xhit)
               iy = idnint(ys - yhit)
c
c     add this photo-electron
c     
               if(ix.gt.ixmin .and.ix.lt.ixmax .and. 
     *              iy.gt.iymin .and. iy .lt.ixmax) then
                  intensity = 1.d0
                  if(subarray(1:4) .ne.'FULL') then
                     ix = ix - colcornr
                     iy = iy - rowcornr
                  end if
                  call add_ipc(ix, iy, intensity, ipc_add)
c                  PRINT *, 'ADD IPC ', IX, IY
               end if
            end do
        end if
 100    continue
      end do
      return
      end
