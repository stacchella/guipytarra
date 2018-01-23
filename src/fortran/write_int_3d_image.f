c
c----------------------------------------------------------------------     
c
      subroutine write_int_3d_image(filename, cube, nx, ny, bitpix,
     *     nframe, tframe, groupgap, tgroup, ngroup, object, partname, 
     *     sca_id, module, filter, 
     *     subarray, colcornr, rowcornr, naxis1, naxis2,job)
c
      implicit none
      character telescop*20, instrume*20, filter*20,
     *     module*20, partname*4,comment*40, object*20      
      integer unit, status, bitpix, naxes, naxis, pcount, gcount,nnn,
     *     block, groupgap,  group, nx, ny, nframe, ngroup,job,
     *     iunit
      integer cube, sca_id
      real tframe, tgroup
      logical simple,extend
      logical subarray
      integer colcornr, rowcornr, naxis1, naxis2
      character filename*120
c
      parameter (nnn=2048)
c
      dimension cube(nnn,nnn,30), naxes(3)
c
c     define the required primary array keywords
c
      simple = .true.
      extend = .true.
      naxis    =  3
      naxes(1) = nx
      naxes(2) = ny
      naxes(3) = ngroup

      status = 0
      call ftgiou(iunit, status)
c
c     delete previous version of the file, if it exists
c
      call ftopen(iunit, filename, 1, block, status)
      if (status .eq. 0) then
         call ftdelt(iunit, status)
      else
c
c     clear the error message stack
c     
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
c     write basic header
c
      call ftphpr(iunit,simple, bitpix, naxis, naxes, 
     & 0,1,extend,status)
      if(status .ne.0) then
         print *, 'simple,bitpix,naxis'
         call printerror(status)
      end if
      status =  0
c
c     write additional keywords
c
      call write_nircam_keywords(iunit,nframe, tframe, groupgap, tgroup, 
     *     ngroup, object, partname, sca_id, module, filter, 
     *     subarray,colcornr, rowcornr, naxis1, naxis2, job)
c
c     write image
c
      group=0

      print *,'write int 3D matrix'
      if(bitpix.eq.16) then
         call ftp3dj(iunit,group,nnn, nnn, naxes(1), naxes(2),
     *        naxes(3), cube, status)
      end if
      
      if(bitpix.eq.32) then
         call ftp3dj(iunit,group,nnn, nnn, naxes(1), naxes(2),
     *        naxes(3), cube, status)
      end if
      if (status .gt. 0) then
         call printerror(status)
         print *, 'writing images'
      end if
c     
      call closefits(iunit)
      close (iunit)
      return
      end
