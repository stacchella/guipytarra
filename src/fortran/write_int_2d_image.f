      subroutine write_int_2d_image(filename, image, nx, ny, bitpix,
     *     nframe, tframe, groupgap, tgroup, ngroup, object, partname, 
     *     sca_id, module, filter, 
     *     subarray, colcornr, rowcornr, naxis1, naxis2, job)
c
      implicit none
      character telescop*20, instrume*20, filter*20,
     *     module*20, partname*4,comment*40, object*20      
      real tframe, tgroup
      integer image, nnn, nx, ny, job, iunit
      integer status, bitpix, naxes, naxis, pcount, gcount, block,
     *     groupgap,sca_id, nframe, ngroup, group
      logical simple,extend
      logical subarray
      integer colcornr, rowcornr, naxis1, naxis2
       character filename*120
      double precision  bscale, bzero
c
      parameter (nnn=2048)
c
      dimension image(nnn,nnn),naxes(2)
c
c     define the required primary array keywords
c
      bzero    = 32678.d0
      bscale   =     1.d0
      simple   = .true.
      extend   = .true.
      naxis    =  2
      naxes(1) = nx
      naxes(2) = ny
c
      status = 0
      call ftgiou(iunit, status)
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
c     write more header keywords
c
      call write_nircam_keywords(iunit,nframe, tframe, groupgap, 
     *     tgroup, ngroup, object, partname, sca_id,module, filter, 
     *     subarray, colcornr, rowcornr, naxis1, naxis2, job)
c
c     Write out image
c
      print *,'write 2d int'
      group=0
      if(bitpix.eq.16) then
         status =  0
c         call ftpkyd(iunit,'BSCALE',bscale,' ',status)
         if (status .gt. 0) then
            call printerror(status)
            print *, 'BSCALE'
         end if
         status =  0
c         call ftpkyd(iunit,'BZERO',bzero,' ',status)
         if (status .gt. 0) then
            call printerror(status)
            print *, 'BSCALE'
         end if
c         call ftpscl(iunit,bscale,bzero, status)
         call ftp2dj(iunit,group,nnn, naxes(1),naxes(2),
     *        image,status)
      end if
      if(bitpix.eq.32) then
         call ftp2dj(iunit,group,nnn, naxes(1),naxes(2),
     *        image,status)
      end if

      if (status .gt. 0) then
         call printerror(status)
         print *, 'writing images'
      end if
c     
      call closefits(iunit)
      return
      end
c
