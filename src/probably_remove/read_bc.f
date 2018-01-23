c
c-----------------------------------------------------------------------
c
c     Given a value of internal absorption and burst length, find
c     the closest family of models
c     
      subroutine find_model(av_in, tau_in, zm, model)
      implicit none
      double precision av_in, tau_in, zm, d_a_v, d_tau, dz_m,
     *     tau, a_v, zmetal
      double precision chimin, chi2
      real par, ages, wl, flux
      integer nwl, nsteps
      integer model
      integer nnn, mmm, nx, ny, nz, np, i
      parameter (nnn=12000, mmm=500)
      parameter (nx=1221, ny=221, nz=1626, np= 8)
      dimension par(np,ny,nz), wl(nnn), flux(nnn,mmm), ages(mmm)
      common /bc_models/ par, ages, nsteps, wl, flux, nwl
      
      chimin  = 10000.d0
      do i = 1, nz
         a_v    = dble(par(6,2,i))
         tau    = dble(par(7,2,i))
         zmetal = dble(par(8,2,i))
         chi2   = (a_v-av_in)**2 + (tau-tau_in)**2 +(zmetal-zm)**2
c         print 220, i,(par(k,2,i),k = 1, np)
 220     format(i10,8(1x,f12.8))
c         print 220, i, chi2, chimin, a_v, tau, zmetal
         if(chi2.le.chimin) then
            chimin = chi2
            model = i
         end if
      end do
c      print *,'find_model : av_in, tau_in, zm_in',av_in,tau_in, zm
c      print *,'find_model : a_v, tau, zmetal    ', par(6,2,model),
c     *      par(7,2,model),  par(8,2,model), model
      return
      end
c
c-----------------------------------------------------------------------
c
c     Recover model parameters (mbol, stellar mass etc.) given
c     model and age
c
c     2015-09-14
c     cnaw@as.arizona.edu
c
      subroutine get_bc_2003_params(age_index, model)
      implicit none
c
      real par
      real mbol, stellar_mass, gas_mass, galaxy_mass,
     *     sfr, a_v, tau, zmetal
      real ages, wl, flux
c
      integer age_index, model
      integer nx, ny, nz, np, nsteps, nwl, nnn, mmm
      integer number_par
c
      parameter (nx=1221, ny=221, nz=1626, np= 8)
      parameter (nnn=12000, mmm=500)
      dimension par(np,ny,nz), wl(nnn), flux(nnn,mmm), ages(mmm)
c      dimension mbol(mmm), stellar_mass(mmm), 
c     *     gas_mass(mmm), galaxy_mass(mmm), sfr(mmm),zmetal(mmm),
c     *     a_v(mmm), tau(mmm)
c
      common /bc_models/ par, ages, nsteps, wl, flux, nwl
      common /bc_params/mbol,stellar_mass,gas_mass,galaxy_mass,
     *     sfr, a_v, tau, zmetal, number_par

      mbol            = par(1, age_index, model)
      stellar_mass    = par(2, age_index, model)
      gas_mass        = par(3, age_index, model)
      galaxy_mass     = par(4, age_index, model)
      sfr             = par(5, age_index, model)
      a_v             = par(6, age_index, model)
      tau             = par(7, age_index, model)
      zmetal          = par(8, age_index, model)
      return
      end
c
c-----------------------------------------------------------------------
c
c     read B&C models from a FITS data cube, probably useless
c
      subroutine fits_bc_2003(index)
      implicit none
      integer index, level, nplane
      character fits_file*120
c
      fits_file = 'bc2003.fits'
c 
c     if index = 0 read parameters, ages and wavelengths (i.e., initialise)
      level = index
      call read_bc_cube(fits_file,  nplane, level)
      return
      end
c
c
c-----------------------------------------------------------------------
c
c     Bruzual & Charlot SEDs come in units of Solar luminosities per
c     angstrom. This routine will convert them into erg/(sec A cm**2)
c     by "measuring" the flux at a distance of 10 pc.
c     
      subroutine read_bc_sed(wlsed, sed, nsed, alpha, zf, zi, index,
     *     tt, ibest, debug)
      implicit none
      double precision wlsed, sed, alpha, zf, zi, tt
      double precision dl, tlookback, tform, dust_attenuation, th_in_gyr
      double precision l_sun, ltot, extinction, zform, t0,t0gyr, th_gyr,
     *     t_lookback_0, t_lookback_i, ti, dt, dtlb, diff, fact, parsec, 
     *     pi, four_pi, t_hubble,tiGyr, age_gyr
      double precision h0, omega_m, omega_l
c
      real par, ages, wl, flux
      integer nsteps, nwl
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
c
      data init/1/
c
      pi          = dacos(-1.0d0)
      four_pi     = 4.d0 * pi
      parsec      = 3.08567802d18    ! cm
      l_sun       = 3.842d33    ! erg/s (http://www.pas.rochester.edu/~emamajek/sun.txt
      ltot        = 1.0d0   ! Solar luminosities
c      ltot        = 10.d0**8
c     This will give the measured flux for a source at 10 pc away
      fact        = ltot*l_sun / (four_pi* (10.d0*parsec)**2)
      file        = 'bc2003.fits'
c
      if(init.eq.1) then
         init = 0
         call read_bc_cube(file,  nplane, init)
c         do i = 1, nz/6
c            print 111,i, par(1, 10,i),  par(4,10,i),
c     *           par(6,10,i), par(7,10,i), par(8,10,i)
c 111        format(i5,8(1x,f9.5))
c         end do
c         pause
      end if
c
c     Number of Giga-years corresponding to a Hubble Time
c
      t_hubble = 1.d0/h0
      th_gyr   = th_in_gyr(1.d0)
c
c     This is the age of the Universe at the galaxy formation redshift
c
      zform = zf
      t0    = tform(zform)      ! in Hubble times
      t0Gyr =  th_in_gyr(t0/t_hubble) ! in Gyr
      t_lookback_0 = th_in_gyr(tlookback(zform)/t_hubble)
c     
c     age of the universe at the observed redshift
c
      ti           = tform(zi)
      tiGyr        = th_in_gyr(ti/t_hubble)
      t_lookback_i = th_in_gyr(tlookback(zi)/t_hubble) !* th_gyr
c
c     This is the difference in Gyrs since formation (== population Age)
c
      dt =  t_lookback_0 - t_lookback_i
c
c     Read the family of models corresponding to "index"; these families
c     consist of spectra corresponding to 221 ages of the stellar population
c
      call read_bc_cube(file,  nplane, index)
c
c     Find closest age, given formation redshift (zf) and current redshift
c
      ibest = 0
      diff       = 1.d030
      do i = 1, nsteps
         tt = ages(i)/1.d09   ! convert into Gyrs
         if(dabs(tt-dt).le.diff.and.tt.lt.dt) then
            diff = dabs(tt-dt)
            ibest = i
         end if
      end do
c      print *,'index, ibest ', index, ibest, diff,zi
c
      if(diff.eq.1.d030) then
         print 60, zi, dt, zf
 60      format('Galaxies were not formed yet !',
     *        4(1x,f8.4))
         index = 0
         return
      end if
         
      tt = ages(ibest)/1.d09
c
c      age_gyr = ages(index)/1.d09
c      if(age_gyr.gt.dt) then
c
      if(debug.eq.1) then
         print *, ' '
c      print *,'index, ibest ', index, ibest
c         print *, ' read_bc_sed: Hubble time in Gyr ',th_gyr
         print 40,  zf, t0Gyr,  t_lookback_0
 40      format('read_bc_sed: zform, Age of Universe, lookback time ',
     *       3(1x,f9.5))
         print 50, zi, tiGyr,t_lookback_i
 50      format('read_bc_sed: z    , Age of Universe, lookback time ',
     *       3(1x,f9.5))
         print *, 'read_bc_sed : Time since formation         ', dt
         print *, 'Selected age:                              ', tt
      end if
c
c
c     These are the parameters for this a_v, z_form, and age
c
      mbol               = par(1, ibest, index)
      stellar_mass       = par(2, ibest, index)
      gas_mass           = par(3, ibest, index)
      galaxy_mass        = par(4, ibest, index)
      sfr                = par(5, ibest, index)
      a_v                = par(6, ibest, index)
      tau                = par(7, ibest, index)
      zmetal             = par(8, ibest, index)

      if(debug.eq.1) then
c         print 300, zform, zi, tt, ti, t0, diff, ibest,index
         print 300, zform, zi, tt, mbol, stellar_mass, a_v, tau, zmetal
 300     format('Zform = ', f9.5, ' z = ', f9.5,
     *        '   Pop age =',f9.5 ,' Gyr',2x,'Mbol  = ', f9.5, 
     *        ' Stellar mass = ', f9.5,1x,' a_v = ',f9.5,2x,
     *        'tau   = ', 1Pe9.2,1x,' Zmetal = ', e9.2)
c         print *, 'Zform = ', zf, ' z = ', zi, ' Pop age = ',tt, ' Gyr'
c         print *, 'mbol = ' , mbol      ,' stellar mass =',
c     *        stellar_mass      , ' a_v  = ', a_v       
c         print *, ' tau = ', tau      ,' ZMetal = ',zmetal      
      end if
      nsed = 0
      do i = 1, nwl
         nsed = nsed + 1
         wlsed(nsed) = wl(i)/10000.d0 ! convert from A --> microns
         extinction = 1.0d0
         extinction  = dust_attenuation(wlsed(nsed),1.0d0)
c
c     convert flux from units of solar Luminosity per Angstrom
c     into erg/(second Angstrom cm**2)
c
c        if(mod(i,100).eq.1)  print *, i, ibest, flux(i,ibest)
         sed(nsed)   = flux(i,ibest) * fact !* extinction
      end do
      return
      end
c
c-----------------------------------------------------------------------
c
c     read a single BC 2003 model and store. Note that all 
c     composite spectra (those passed through csp), the first
c     time step has an SED = 0, which probably means that one should be
c     using the original SSP.
c
      subroutine get_bc_2003(index, ages, nsteps, wl, nwl,flux)
      implicit none
      real ages, wl, flux
      integer index, nsteps, nwl
      real mbol, log_age, m_over_lb, stellar_mass, gas_mass, 
     *     galaxy_mass, sfr
      integer mmm, nnn
      integer i, npar
      character files*80, params*80
      parameter (nnn=12000, mmm=500)
      dimension ages(mmm), wl(nnn), flux(nnn,mmm),
     *     files(6), params(6)
      dimension log_age(mmm), mbol(mmm), m_over_lb(mmm), 
     *     stellar_mass(mmm), gas_mass(mmm), galaxy_mass(mmm),
     *     sfr(mmm)
      data log_age, mbol, m_over_lb/mmm*0.0, mmm*0.0,mmm*0.0/
      data stellar_mass,gas_mass,galaxy_mass/mmm*0.0,mmm*0.0,mmm*0.0/
      data sfr/mmm*0.0/
c
c     Solar metallicity is z=0.02
c
      open(10,file='bc_2003_grid.list')
c      open(10,file='bc_2003_lowres.list')
c      open(10,file='bc_2003_highres.list')
      do i = 1, 6
         read(10,10,end=100) files(i)
         read(10,10,end=100) params(i)
 10      format(a80)
      end do
 100  close(10)

      call read_bc_ascii_file(files(index), ages, nsteps, wl, nwl,flux)
      call read_bc_param_file(params(index), log_age, mbol, 
     *     stellar_mass, gas_mass, galaxy_mass, sfr, npar)
      print 10, files(index)
c      do i = 2, nsteps
c         print *, i, alog10(ages(i)), log_age(i)
c      end do
c      print *, nsteps, nwl
      return
      end
c
c-----------------------------------------------------------------------
c
c     Read the BC2003 ASCII SSP files. It also reads the composite
c     populations generated by csp_galaxev.sh.
c
      subroutine read_bc_ascii_file(file, ages, nsteps, wl, nwl,
     *     flux)
      implicit none
c     Array declarations
      logical stelib
      character id*80,id2*80,id3*80,name*256, file*80
c
      real ages, wl, flux
      integer nsteps, nwl
      real  totm, totn, avs, tau, tau1, tau2, tau3, tau4
      real xml, xmu
      integer iseg, ix
      integer jo, iop
      integer imw, imf, mmm, n, i, ik

c     Maximum number of wavelength points
      PARAMETER (imw=12000)
      parameter (imf=10)
      parameter (mmm=500)
      real w(imw),h(imw),tb(0:500),f(100)
      real ml, mu
c     Array declaration for IMF data
      real xx(imf),lm(imf),um(imf),baux(imf),cc(imf),cn(imf),taupar(4)
      dimension ages(mmm), wl(mmm), flux(imw,mmm)

      open (1,file=file)
c
c     These are the commands to read files resulting from  csp_galaxev.sh
c
c     The ages range form 0 (=t_form) to 20 Gyr (t_form+20)
c
      read(1,*) nsteps, (ages(i),i=1, nsteps)
c
c     xml = lower mass cutoff (should be 0.1 M_sol)
c     xmu = upper mass cutoff (should be 100 M_sol)
c
      read(1,*) xml,xmu, iseg,(xx(i),lm(i),um(i),baux(i),cn(i),
     *     cc(i),i=1,iseg)
c
c     Total mass,
c       
      read(1,*)  totm,totn,avs,jo,tau,tau1,tau2,tau3,tau4,iop,stelib
c      print *,  totm,totn,avs,jo,tau,tau1,tau2,tau3,tau4,iop,stelib
c      read(1,'(a)') stelib
      read (1,'(a)') id
      read (1,'(a)') id2
      read (1,'(a)') id3
c      print 10, stelib, id, id2, id3
 10   format(a80)
      read(1,*)  nwl, (wl(i),i=1, nwl)
c     
c     Read seds from input file
c     Include fluxes for fitting functions
c
      do n=1,nsteps
         read(1,*)   ik,(h(i),i=1,ik),ix,(f(i),i=1,ix)
c         print *,'nstep ', n, ages(n), ik, ix
c
c     create array of age x flux . Units are L_sun/A with L_sun=3.826d33 erg/s
c
         do i = 1, ik
            flux(i,n) = h(i)
         end do
      end do
c      print *,'nsteps, ik, ix ', nsteps, ik, ix
c
c     Copy 12 records after the seds.
c
      do n=1,12
         i = 0
         read  (1,*,end=100)   ik,(h(i),i=1,ik)
c         print *, n,ik, h(1), h(100), h(ik)
      end do
 100  close (1)
      return
      end
c
c-----------------------------------------------------------------------
c
c     Read the BC2003 model parameters. These are generated as
c     a function of age
c
      subroutine read_bc_param_file(file, log_age, mbol, stellar_mass,
     *    gas_mass, galaxy_mass, sfr, nn)
c     Array declarations
      implicit none
      real log_age, mbol, m_over_lb, n_ly_alpha, neutron_stars,
     *     stellar_mass, gas_mass, galaxy_mass, sfr, sn_rate,
     *     pn_rate, black_holes, white_dwarfs, remnants
      real  age, ambol, amb, amv, am_over_lb, am_over_lv, 
     *        stmass, gas, gal, sfryr
      integer nn, mmm, npar, i, k
      logical stelib
      character id*80,id2*80,id3*80,name*256, file*80
      
c     Maximum number of wavelength points
      PARAMETER (mmm=500, npar=220)
      dimension log_age(mmm), mbol(mmm), m_over_lb(mmm), 
     *     stellar_mass(mmm), gas_mass(mmm), galaxy_mass(mmm),
     *     sfr(mmm), n_ly_alpha(mmm), sn_rate(mmm), PN_rate(mmm),
     *     black_holes(mmm), neutron_stars(mmm),white_dwarfs(mmm),
     *     remnants(mmm)
c
c     This is the *4color file
c
      open(1,file = file)
      do i = 1, 29
        read(1,*)
      end do
c
      do k =1,npar
         read(1,*,err=10) age, ambol, amb, amv, am_over_lb, am_over_lv, 
     *        stmass, gas, gal, sfryr
         go to 30
 10      print 20, age, ambol, amb, amv, am_over_lb, am_over_lv, 
     *        stmass, gas, gal, sfryr
 20      format(10(1x,e13.2))
 30      continue
         i = k + 1
         log_age(i)       = age
         mbol(i)          = ambol
         m_over_lb(i)     = am_over_lb
         stellar_mass(i)  = stmass
         gas_mass(i)      = gas
         galaxy_mass(i)   = gal
         sfr(i)           = sfryr
      end do
      close(1)
      nn = i
c
c     read the *3color file
c
c      open(1,file = file)
c      do i = 1, 29
c        read(1,*)
c      end do
cc
c      do i =1,npar
c         read(1,*) age, b4000, b4_vn, b4_sdss, b912, anly_alpha,
c     *        snr, pnr, bhn, starn, wd, rem
c         n_ly_alpha(i)    = anly_alpha
c         sn_rate(i)       = snr
c         PN_rate(i)       = pnr
c         black_holes(i)   = bhn
c         neutron_stars(i) = starn
c         white_dwarfs(i)  = wd
c         remnants(i)      = rem
c      end do
c      close(1)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine closest_age(zform, zi, age_index)
      implicit none
      double precision zi, zform
      integer age_index
c
      double precision h0, omega_m, omega_l
      double precision  th_in_gyr,  tform, tlookback
      double precision t_hubble, th_gyr, zf, t0, t0gyr,
     *  t_lookback_0, ti, tigyr, t_lookback_i, dt, diff, tt
c
      real par, wl, flux, ages
c
      integer np, ny, nz, nbc, mbc, nwl, nsteps
      integer ibest, i
c
      parameter (ny=221, nz=1626, np= 8)
      parameter (nbc=12000, mbc=500)
c
      dimension par(np,ny,nz), wl(nbc), flux(nbc,mbc), ages(mbc)

      common /cosmology/ h0, omega_m, omega_l
      common /bc_models/ par, ages, nsteps, wl, flux, nwl
c
c     Number of Giga-years corresponding to a Hubble Time
c
      t_hubble = 1.d0/h0
      th_gyr   = th_in_gyr(1.d0)
c
c     This is the age of the Universe at the galaxy formation redshift
c
      t0    = tform(zform)      ! in Hubble times
      t0Gyr =  th_in_gyr(t0/t_hubble) ! in Gyr
      t_lookback_0 = th_in_gyr(tlookback(zform)/t_hubble)
c     
c     age of the universe at the observed redshift
c
      ti           = tform(zi)
      tiGyr        = th_in_gyr(ti/t_hubble)
      t_lookback_i = th_in_gyr(tlookback(zi)/t_hubble) !* th_gyr
c
c     This is the difference in Gyrs since formation (== population Age)
c
      dt =  t_lookback_0 - t_lookback_i
c
c      print *,'read_bc:closest_age ', zform, t0, zi, t_lookback_i,dt
c
c     Read the family of models corresponding to "index"; these families
c     consist of spectra corresponding to 221 ages of the stellar population
c
c      call read_bc_cube(file,  nplane, index)
c
c     Find closest age, given formation redshift (zf) and current redshift
c
      ibest = 0
      diff       = 1.d030
      do i = 1, nsteps
         tt = ages(i)/1.d09   ! convert into Gyrs
         if(dabs(tt-dt).le.diff.and.tt.lt.dt) then
            diff = dabs(tt-dt)
            ibest = i
         end if
      end do
c      print *,'i, ibest ', i, ibest, diff,zi,tt
c
      if(diff.eq.1.d030) then
         print 60, zi, dt, zf
 60      format('Galaxies were not formed yet !',
     *        4(1x,f8.4))
         age_index = 0
         return
      end if
      age_index = ibest
      return
      end
