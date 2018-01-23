      subroutine get_bc_spectrum(model, z, zform, age_index, smass, 
     *     wlsed, sed, nsed)
      implicit none
      double precision wlsed, sed, alpha, zf, zi, tt, z
      double precision dl, tlookback, tform, dust_attenuation, th_in_gyr
      double precision l_sun, extinction, zform, t0,t0gyr, th_gyr,
     *     t_lookback_0, t_lookback_i, ti, dt, dtlb, diff, fact, parsec, 
     *     pi, four_pi, t_hubble,tiGyr, age_gyr, smass
      double precision h0, omega_m, omega_l
c
      real par, ages, wl, flux
      integer nsteps, nwl
      integer age_index, model, j
      character fits_file*80
c
      real mbol, stellar_mass, gas_mass, galaxy_mass,
     *     sfr, a_v, tau, zmetal
      integer number_par
c
      integer nsed, index, debug, nplane
      integer ibest, nnn, mmm, lll, npar, nx, ny, nz, np, nbc, mbc,
     *     init, i
      character file*80
      parameter (nnn=25000,mmm=500,lll=1000,npar=30)
      parameter (nx=1221, ny=221, nz=1626, np= 8)
      parameter (nbc=12000, mbc=500)
      dimension par(np,ny,nz), wl(nbc), flux(nbc,mbc), ages(mbc)
      dimension wlsed(nnn), sed(nnn)
c      dimension mbol(mmm), stellar_mass(mmm), 
c     *     gas_mass(mmm), galaxy_mass(mmm), sfr(mmm),zmetal(mmm),
c     *     a_v(mmm), tau(mmm)
c
      common /cosmology/ h0, omega_m, omega_l
      common /bc_models/ par, ages, nsteps, wl, flux, nwl
      common /bc_params/mbol,stellar_mass,gas_mass,galaxy_mass,
     *     sfr, a_v, tau, zmetal, number_par
      external dl
      save init, fits_file
c
      data init/1/
C
      pi          = dacos(-1.0d0)
      four_pi     = 4.d0 * pi
      parsec      = 3.08567802d18    ! cm
c
c     BC2003  units are solar luminosities per Angstrom
c
      l_sun       = 3.828d33    ! erg/s (from 2014 LBNL list of constants)
c 
c     This will give the measured flux for a source at 10 pc away
c     in ergs/(steradian * second * cm**2 * A)
c
      fact        = l_sun / (four_pi* (10.d0*parsec)**2)
c
c     find model
c
      fits_file = 'bc2003.fits'
      if(init. eq. 1) then
         call read_bc_cube(fits_file, nplane, 0)
         init = 0
      end if
      call read_bc_cube(fits_file, nplane, model)
c
c     find age closest to redshift
c
      call closest_age(zform, z, age_index)
      if(age_index.eq. 0) return
c
      call get_bc_2003_params(age_index, model)
c
c     Scale the BC models to the desired stellar mass as models are
c     calculated for 1 solar mass
c      print *, 'read_bc:get_spectrum: model, smass, fact ',
c     *     model, smass, fact, stellar_mass, gas_mass, mbol
      nsed = nwl
      do j = 1, nwl
         wlsed(j) = wl(j)/10000.d0 ! convert from A --> microns
         sed(j) = flux(j,age_index) * fact
         sed(j) = sed(j) * smass  * stellar_mass
      end do
      return
      end
