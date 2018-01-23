c
c----------------------------------------------------------------------
c
      subroutine write_nircam_keywords(iunit,nframe, tframe, 
     *     groupgap, tgroup, ngroup, object, partname, sca_id, 
     *     module, filter, dhas_subarray, colcornr, rowcornr, job)
      implicit none
      logical dhas_subarray
      integer colcornr, rowcornr
      integer ibrefrow, itrefrow, llrefcol, rrefcol,
     *     nframe, ngroup, iunit, job
      integer status, groupgap, drop_frame_1,sca_id
      character telescop*20, instrume*20, filter*5,
     *     module*4, partname*5,comment*40, object*20,
     *     cunit1*4, cunit2*4
      real tframe, tgroup
      double precision 
     *     equinox, crpix1, crpix2, crval1, crval2,
     *      cdelt1, cdelt2, cd1_1, cd1_2, cd2_1, cd2_2, cd3_3
      common /wcs/ equinox, crpix1, crpix2, crval1, crval2,
     *     cdelt1,cdelt2, cd1_1, cd1_2, cd2_1, cd2_2


c     Put information into keywords and write them
c     
c     The following are fixed...
c
      instrume = 'NIRCam'
      telescop = 'Az_Lab'
      if(dhas_subarray .eqv. .True.) then
         ibrefrow = 0
         itrefrow = 0
      else
         ibrefrow = 4
         itrefrow = 4
      end if
      drop_frame_1 = 0
      cd3_3        = 1.00000
c
      status =  0
      comment='Number of frames in group'
      call ftpkyj(iunit,'NFRAME',nframe,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'NFRAME'
       end if
      status =  0
c
      comment='Number groups in an integration'
      call ftpkyj(iunit,'NGROUP',ngroup,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'NGROUP'
       end if
      status =  0

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
      comment='Number of frames skipped'
      call ftpkyj(iunit,'GROUPGAP',groupgap,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'GROUPGAP'
       end if
      status =  0
c
c
c
c     job number so it is compatible with NIRCam tests
c
      comment='                                       '
      call ftpkyj(iunit,'GS_JOBID',job,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'GS_JOBID'
       end if
      status =  0
c
      comment='                                       '
      call ftpkys(iunit,'OBJECT ',object,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'object'
       end if
      status =  0
c
      comment='                                       '
      call ftpkys(iunit,'TELESCOP',telescop,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'telescop'
       end if
      status =  0
c
      comment='                                       '
      call ftpkys(iunit,'INSTRUME',instrume,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'instrume'
       end if
      status =  0
c
      comment='                                       '
      call ftpkys(iunit,'FILTER',filter,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'FILTER'
          print 40, filter
 40       format(a20)
       end if
      status =  0
c
      comment='                                       '
      call ftpkys(iunit,'PARTNAME',partname,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'PARTNAME'
          print 40, partname
       end if
      status =  0
c
      comment='                                       '
      call ftpkyj(iunit,'SCA_ID',sca_id,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'SCA_ID'
       end if
      status =  0
c
      comment='                                       '
      call ftpkys(iunit,'MODULE',module,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'MODULE'
          print 40, module
       end if
      status =  0
c
c
c     Sub-array related keywords
c
      comment='                                       '
      call ftpkys(iunit,'SUBARRAY',dhas_subarray,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'SUBARRAY'
       end if
      status =  0
c
c      comment='                                       '
c      call ftpkyj(iunit,'SUBARSZX',subarszx,comment,status)
c       if (status .gt. 0) then
c          call printerror(status)
c          print *, 'SUBARSZX'
c       end if
c      status =  0
cc
c      comment='                                       '
c      call ftpkyj(iunit,'SUBARSZY',subarszy,comment,status)
c       if (status .gt. 0) then
c          call printerror(status)
c          print *, 'SUBARSZY'
c       end if
      status =  0
c
      comment='                                       '
c      call ftpkyj(iunit,'SUBAR_X1',subar_x1,comment,status)
      call ftpkyj(iunit,'COLCORNR',colcornr,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'COLCORNR'
       end if
      status =  0
c
      comment='                                       '
c      call ftpkyj(iunit,'SUBAR_Y1',subar_y1,comment,status)
      call ftpkyj(iunit,'ROWCORNR',rowcornr,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'ROWCORNR'
       end if
      status =  0
c
      comment='bottom reference pixel rows'
      call ftpkyj(iunit,'BREFROW',ibrefrow,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'BREFROW'
       end if
      status =  0
c
      comment='top reference pixel rows'
      call ftpkyj(iunit,'TREFROW',itrefrow,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'TREFROW'
       end if
c              1234567890123456789012345678901234567890
      comment='Number of frame skipped prior to first i'
      call ftpkyj(iunit,'DRPFRMS1',drop_frame_1,comment,status)
      if (status .gt. 0) then
         call printerror(status)
         print *, 'drop_frame_1'
       end if
c       print *,'end keywords'
c
c     fake WCS keywords
c
      comment='                                       '
      call ftpkys(iunit,'CTYPE1','RA---TAN',comment,status)
      if (status .gt. 0) then
         call printerror(status)
          print *, 'CTYPE1'
       end if
       status =  0
c     
       cunit1 = 'deg'
       call ftpkys(iunit,'CUNIT1',cunit1,comment,status)
       if (status .gt. 0) then
          call printerror(status)
          print *, 'CUNIT1'
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
       cunit2 = 'deg'
       call ftpkys(iunit,'CUNIT2',cunit2,comment,status)
      if (status .gt. 0) then
         call printerror(status)
          print *, 'CUNIT2'
       end if
       status =  0
c
      call ftpkys(iunit,'RADESYS','FK5',comment,status)
      if (status .gt. 0) then
         call printerror(status)
          print *, 'RADESYS'
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
c      comment='                                       '
c      call ftpkye(iunit,'CRVAL1',real(crval1),7,comment,status)
c      if (status .gt. 0) then
c         call printerror(status)
c          print *, 'CRVAL1'
c       end if
c       status =  0
cc
c      comment='                                       '
c      call ftpkye(iunit,'CRPIX1',real(crpix1),7,comment,status)
c      if (status .gt. 0) then
c         call printerror(status)
c          print *, 'CRPIX1'
c       end if
c       status =  0
cc
c      comment='                                       '
c      call ftpkye(iunit,'CRVAL2',real(crval2),7,comment,status)
c      if (status .gt. 0) then
c         call printerror(status)
c          print *, 'CRVAL2'
c       end if
c       status =  0
cc
c      comment='                                       '
c      call ftpkye(iunit,'CRPIX2',real(crpix2),7,comment,status)
c      if (status .gt. 0) then
c         call printerror(status)
c          print *, 'CRPIX2'
c       end if
c       status =  0
cc
c      comment='                                       '
c      call ftpkye(iunit,'CD1_1',real(cd1_1),7,comment,status)
c      if (status .gt. 0) then
c         call printerror(status)
c          print *, 'CD1_1'
c       end if
c       status =  0
cc
c      comment='                                       '
c      call ftpkye(iunit,'CD1_2',real(cd1_2),7,comment,status)
c      if (status .gt. 0) then
c         call printerror(status)
c          print *, 'CD1_2'
c       end if
c       status =  0
cc
c      comment='                                       '
c      call ftpkye(iunit,'CD2_1',real(cd2_1),7,comment,status)
c      if (status .gt. 0) then
c         call printerror(status)
c          print *, 'CD2_1'
c       end if
c       status =  0
cc
c      comment='                                       '
c      call ftpkye(iunit,'CD2_2',real(cd2_2),7,comment,status)
c      if (status .gt. 0) then
c         call printerror(status)
c          print *, 'CD2_2'
c       end if
c       status =  0
c
      comment='                                       '
      call ftpkyd(iunit,'CRVAL1',crval1,7,comment,status)
      if (status .gt. 0) then
         call printerror(status)
          print *, 'CRVAL1'
       end if
       status =  0
c
      comment='                                       '
      call ftpkyd(iunit,'CRPIX1',crpix1,7,comment,status)
      if (status .gt. 0) then
         call printerror(status)
          print *, 'CRPIX1'
       end if
       status =  0
c
      comment='                                       '
      call ftpkyd(iunit,'CRVAL2',crval2,7,comment,status)
      if (status .gt. 0) then
         call printerror(status)
          print *, 'CRVAL2'
       end if
       status =  0
c
      comment='                                       '
      call ftpkyd(iunit,'CRPIX2',crpix2,7,comment,status)
      if (status .gt. 0) then
         call printerror(status)
          print *, 'CRPIX2'
       end if
       status =  0
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
      comment='                                       '
      call ftpkye(iunit,'CD3_3',real(cd3_3),4,comment,status)
      if (status .gt. 0) then
         call printerror(status)
          print *, 'CD2_2'
       end if
       status =  0
c     
c     Fake keywords to enable multi-drizzle
c     IDCTAB, ADCGAIN, EXPEND,
c     SAMP_SEQ, NSAMP, CENTERA1, CENTERA2, SIZAXIS1, SIZAXIS2,
c     These already exist:
c     CRPIX1, CRPIX2, CD1_1, CD1_2, CD2_1, CD2_2, NAXIS1, NAXIS2
c
c     IDCTAB  = 'iref$u1r16228i_idc.fits' / image distortion correction table
c     SAMP_SEQ= 'SPARS100'           / MultiAccum exposure time sequence name
c     NSAMP   =                   16 / number of MULTIACCUM samples
c     / READOUT DEFINITION PARAMETERS
c     CENTERA1=   513 / subarray axis1 center pt in unbinned dect. pix
c     CENTERA2=   513 / subarray axis2 center pt in unbinned dect. pix
c     SIZAXIS1=  1024 / subarray axis1 size in unbinned detector pixels
c     SIZAXIS2=  1024 / subarray axis2 size in unbinned detector pixels       
       return
       end

cPHOTFLAM is the flux of a source with constant flux per unit wavelength (in erg s-1 cm-2 Å-1) which produces a count rate of 1 DN per second. This keyword is generated by the synthetic photometry package synphot, which you may also find useful for a wide range of photometric and spectroscopic analyses. Using PHOTFLAM, it is easy to convert instrumental magnitude to flux density, and thus determine a magnitude in a flux-based system such as AB or STMAG (see previous Section); the actual steps required are detailed below. 
