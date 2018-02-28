c
c----------------------------------------------------------------------     
c
      subroutine write_2d(filename, image, nx, ny,
     *     nframe, tframe, groupgap, tgroup, ngroups, nints,
     *     object, targ_ra, targ_dec, equinox,
     *     crpix1, crpix2, crval1, crval2, cdelt1,
     *     cdelt2, cd1_1, cd1_2, cd2_1, cd2_2,
     *     sca_id, module, filter, 
     *     subarray, colcornr, rowcornr)

c
      implicit none
      character telescop*20, instrume*20, filter*20,
     *     module*20, partname*4,comment*40, object*20      
      double precision equinox,crpix1, crpix2, crval1, crval2, cdelt1,
     &     cdelt2, cd1_1, cd1_2, cd2_1, cd2_2
      double precision targ_ra, targ_dec
      real image, tframe, tgroup
      integer status, bitpix, naxes, naxis, pcount, gcount, block,
     *     groupgap, group, sca_id, nx, ny, nframe, ngroups,
     *     nnn, job, iunit, nints
      logical simple,extend
      logical subarray
      integer colcornr, rowcornr, naxis1, naxis2
      character filename*(*), cunit1*8, cunit2*8, string*8
c
      parameter (nnn=2048)
c
      dimension image(nnn,nnn),naxes(2)
c
c     define the required primary array keywords
c
      bitpix   = -32
      simple   = .true.
      extend   = .true.
      naxis    =  2
      naxes(1) = nx
      naxes(2) = ny
c
      status = 0
      iunit = 91
      call ftgiou(iunit, status)
      if (status .gt. 0) then 
         call printerror(status)
         print *,'iunit ',iunit, status
         stop
      end if
c     
c     delete previous version of the file, if it exists
c
      call ftopen(iunit, filename, 1, block, status)
      if (status .eq. 0) then
         call ftdelt(iunit, status)
      else
c     clear the error message stack
         call ftcmsg
      end if
      status = 0
c
c     create the fits file
c
      call ftinit(iunit, filename, 1, status)
      if (status .gt. 0) then 
         call printerror(status)
         print *,'iunit ',iunit
         print *,'pause: enter return to continue'
         read(*,'(A)')
      endif
      status =  0
c     
      call ftphpr(iunit,simple, bitpix, naxis, naxes, 
     & 0,1,extend,status)
      if (status .gt. 0) then
         print *, 'simple,bitpix,naxis'
         call printerror(status)
      end if
      status =  0
c
      comment  = 'Name of coordinate reference frame '
      call ftpkys(iunit,'RADESYS','FK5',comment,status)
      if (status .gt. 0) then
         call printerror(status)
          print *, 'RADESYS'
       end if
       status =  0
c
      string = 'generic'
      comment='Target type (fixed, moving, generic)   '
      call ftpkys(iunit,'TARGTYPE',string,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print 10, string
 10       format(a8)
       end if
      status =  0
c      
      comment = 'Target RA at mid time of exposure (deg)'
      call ftpkyd(iunit,'TARG_RA', targ_ra,6,
     *     comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'TARG_RA'
      end if
      status =  0
c
      comment = 'Target Dec at mid time of exposure (deg)'
      call ftpkyd(iunit,'TARG_DEC', targ_dec,6,
     *     comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'TARG_RA'
      end if
      status =  0
c     
      comment = 'Number of integrations in exposure   '
      call ftpkyj(iunit,'NINTS',nints,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'NINTS'
         status = 0
      end if
      status =  0
c     
      comment='Number groups in an integration'
      call ftpkyj(iunit,'NGROUPS',ngroups,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'NGROUPS'
      end if
      status =  0
c     
      comment='Number of frames in group'
      call ftpkyj(iunit,'NFRAMES',nframe,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'NFRAMES'
      end if
      status =  0
c     
      comment='Number of frames skipped'
      call ftpkyj(iunit,'GROUPGAP',groupgap,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'GROUPGAP'
      end if
      status =  0
c
      comment='Time in seconds between frames'
      call ftpkyf(iunit,'TFRAME',tframe,7,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'TFRAME'
      end if
      status =  0
c     
      comment='Delta time between groups'
      call ftpkyf(iunit,'TGROUP',tgroup,7,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'TGROUP'
      end if
      status =  0
c     
      comment='Axis 1 coordinate of reference pixel   '
      call ftpkyd(iunit,'CRPIX1',crpix1,7,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CRPIX1'
      end if
      status =  0
c     
      comment='Axis 2 coordinate of reference pixel   '
      call ftpkyd(iunit,'CRPIX2',crpix2,7,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CRPIX2'
      end if
      status =  0
c      
      comment='RA at reference pixel (degrees)        '
      call ftpkyd(iunit,'CRVAL1',crval1,7,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CRVAL1'
      end if
      status =  0
c     
      comment='DEC at reference pixel (degrees)       '
      call ftpkyd(iunit,'CRVAL2',crval2,7,comment,status)
      if (status .gt. 0) then
         call printerror(status)
          print *, 'CRVAL2'
       end if
       status =  0
c
      comment='First axis increment per pixel          '      
      call ftpkyd(iunit,'CDELT1',cdelt1,7,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CDELT1'
      end if
      status =  0
c     
      comment='Second axis increment per pixel         '      
      call ftpkyd(iunit,'CDELT2',cdelt2,7,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CDELT2'
      end if
      status =  0
c     
      comment='Projection type                        '
      call ftpkys(iunit,'CTYPE1','RA---TAN',comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CTYPE1'
      end if
      status =  0
c     
      call ftpkys(iunit,'CTYPE2','DEC--TAN',comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CTYPE2'
      end if
      status =  0
c     
      cunit1 = 'deg'
      comment='First axis units                       '
      call ftpkys(iunit,'CUNIT1',cunit1,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CUNIT1'
      end if
      status =  0
c     
      cunit2 = 'deg'
      comment='Second axis units                      '
      call ftpkys(iunit,'CUNIT2',cunit2,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CUNIT2'
      end if
      status =  0
c     
      comment='                                       '
      call ftpkyf(iunit,'EQUINOX',real(equinox),1,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'EQUINOX'
      end if
      status =  0
c     
c     These are non-STScI standard
c     
      comment='                                       '
      call ftpkyd(iunit,'CD1_1',cd1_1,7,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CD1_1'
      end if
      status =  0
c     
      comment='                                       '
      call ftpkyd(iunit,'CD1_2',cd1_2,7,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CD1_2'
      end if
      status =  0
c     
      comment='                                       '
      call ftpkyd(iunit,'CD2_1',cd2_1,7,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CD2_1'
      end if
      status =  0
c     
      comment='                                       '
      call ftpkyd(iunit,'CD2_2',cd2_2,7,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'CD2_2'
      end if
      status =  0
c     
c
c     Write out image
c
      print *,'write float matrix'
      group=0
      call ftp2de(iunit,group,nnn, naxes(1),naxes(2),
     *     image,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'writing images'
      end if
c     
      call closefits(iunit)
      return
      end
