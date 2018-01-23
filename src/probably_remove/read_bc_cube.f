c
c-----------------------------------------------------------------------
c
c     Recover one series of models at a time from the FITS image. For
c     a given SSP, A_V, tau, it will recover the SEDs for all 221 ages.
c     Thus, it is not necessary to store the full data cube in memory,
c     though the entire table of parameters IS stored.
c
      subroutine read_bc_cube(filename, nplane, level)
      implicit none
      real par, ages, wl, flux, cube
      integer nsteps, nwl
      integer unit, status, bitpix, naxes, naxis,pcount, gcount,group,
     *     colnum
      integer frow,felem, fpixels, lpixels, incs
      integer hdu_type, debug, nfound, readwrite
      integer mmm, nnn, nx, ny, nz, np
      integer i, j, init, null, nplane, level
      character comment*20, filename*80
      logical simple, extend, anyf, flagvals
c
      parameter (mmm=500, nnn=12000)
      parameter (nx=1221, ny=221, nz=1626, np=8)
c
      dimension naxes(3), fpixels(3), lpixels(3), incs(3)
      dimension cube(nx,ny,1) 
      dimension par(np,ny,nz), ages(mmm),
     *     wl(nnn),flagvals(nnn), flux(nnn,mmm)
c     
      common /bc_models/ par, ages, nsteps, wl, flux, nwl
      data init/1/
c
c     There are 16 values of a_v and 18 of tau:
c     a_v =  [0.000, 0.010, 0.020, 0.030, 0.050, 0.080, 0.100, 0.200, $
c         0.300, 0.500, 0.800, 1.000, 2.000, 3.000, 5.000, 8.000]
c
c     tau = [0.000, 0.010, 0.020, 0.030, 0.050, 0.080, 0.100, 0.200, 0.300, $
c         0.500, 0.800,1.000, 2.000, 3.000, 5.000, 8.000,10.00, 20.00]
c
c     Each of the lists starts with (a_v, tau) = (0.0,0.0) and then
c     follows with (a_v, tau) = (0.01-8, 0.0-20), so there is a total
c     of 15 x 18 + 1 = 271 models for each metallicity. Since these are
c     a total of 6, there are 1626 models
c
c     
      debug  = 0
      status = 0
      group  = 0
      null   = 0
      readwrite = 0
c
c     open fits file
c
      unit = 88
      call openfits(unit, readwrite, filename)
      if(debug.eq.1) print *,'unit ',unit
c
c     on the first round read parameters, wavelengths and ages
c
      if(init.eq.1) then
         init = 0
c
c     read parameters, which are a fits cube in the second HDU
c
         call ftmahd(unit,2, hdu_type, status)
         if(debug.eq.1) print *,'ftmahd ', status
         call ftgkyj(unit,"BITPIX",bitpix,comment, status)
         call printerror(status)
         status = 0 
         call ftgkyj(unit,"NAXIS",naxis,comment, status)
         call printerror(status)
         status = 0 
         call ftgknj(unit,'NAXIS',1,naxis,naxes,nfound,status)
         call printerror(status)
         status = 0 
         if(debug.eq.1) print *, 'par ', bitpix, naxis, naxes
         call ftg3de(unit, group, null, np, ny, naxes(1),
     *        naxes(2), naxes(3), par, anyf, status)
         if(debug.eq.1) then
            print *,'ftg3de ', status
            print *,(par(i,2,1),i = 1, np)
         end if
c
c     read wavelength table
c
         call ftmahd(unit,3, hdu_type, status)
         call ftgkyj(unit,"NAXIS",naxis,comment, status)
         call printerror(status)
         status = 0 
         call ftgknj(unit,'NAXIS',1,naxis,naxes,nfound,status)
         call printerror(status)
         status = 0 
c
         frow   = 1
         colnum = 1
         felem  = 1
c     naxes(1) contains the number of bytes...
         nwl    = naxes(2)
         call ftgcfe(unit, colnum, frow, felem, nwl, 
     *        wl,flagvals,anyf,status)
         if(debug.eq.1) print *,'ftgcl 1', status
c
c     read ages table
c     
         call ftmahd(unit,4, hdu_type, status)
         call ftgkyj(unit,"NAXIS",naxis,comment, status)
         call printerror(status)
         status = 0 
         call ftgknj(unit,'NAXIS',1,naxis,naxes,nfound,status)
         call printerror(status)
         status = 0 
c
         frow   = 1
         colnum = 1
         nsteps = naxes(2)
         call ftgcfe(unit, colnum, frow, felem, nsteps, 
     *        ages,flagvals,anyf,status)
         if(debug.eq.1) print *,'ftgcl 2', status
      end if
      if(level.eq.0) return
c
c     read data cube of SEDs x Ages x (Z + A_V + Tau)
c     which is located in the first HDU. Code may be changed so that
c     only a portion of the image is read.
c
      call ftmahd(unit,1, hdu_type, status)
      call ftgkyj(unit,"BITPIX",bitpix,comment, status)
      call printerror(status)
      status = 0 
      call ftgkyj(unit,"NAXIS",naxis,comment, status)
      call printerror(status)
      status = 0 
      call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
      call printerror(status)
      status = 0 
      if(debug .eq.1) print *, bitpix, naxis, naxes
      nplane      = naxes(3)
c
c     fpixels and lpixels indicate the locations of the first and
c     and last pixels that must be retrieved
c
      fpixels(1) = 1
      lpixels(1) = naxes(1) 
      incs(1)    = 1
c
      fpixels(2) = 1
      lpixels(2) = naxes(2) 
      incs(2)    = 1
c
      fpixels(3) = (level - 1) + 1
      lpixels(3) = level
      incs(3)    = 1
c
c     this works to read entire array
c
c      fpixels(1)  = 1
c      lpixels(1)  = naxes(1)
c      incs(1)     = 1
cc
c      fpixels(2)  = 1
c      lpixels(2)  = naxes(2)
c      incs(2)     = 1
cc
c      fpixels(3)  = 1
c      lpixels(3)  = naxes(3)
c      incs(3)     = 1

c      print 210,fpixels, lpixels
 210  format(10(1x,i10))
c
c     The 3-D image will only be populated in its first
c     lpixels(3) - fpixels(3) levels (the same holds for
c     other dimensions)
c
c      print *,'enter ftgsve'
      call ftgsve(unit, group, naxis, naxes, fpixels,
     *     lpixels, incs, null, cube, anyf, status)
c      print *,'exit ftgsve'
c
      do i = 1, naxes(1)
         do j = 1, naxes(2)
            flux(i,j) = cube(i,j,1)
         end do
      end do
c      print *,'stuffed flux array'
c     
c     Ex. print values of pixel x=200, y=100, z=level in the original 
c     matrix:
c      print 220, status, (cube(200,100,1))
 220  format(i10,20(1x,f10.7))
      call printerror(status)
      call closefits(unit)
      return
      end
