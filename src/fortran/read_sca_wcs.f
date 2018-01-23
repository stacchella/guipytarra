c
c-----------------------------------------------------------------------
c
      subroutine read_sca_wcs(sca_id,x_sci_ref, y_sci_ref,
     &     x_sci_scale, y_sci_scale,
     &     v2_ref, v3_ref, v3_parity, v3idl_angle)
      implicit none
      double precision x_sci_ref, y_sci_ref, x_sci_scale, y_sci_scale,
     &     v2_ref, v3_ref, v3_parity, v3idl_angle
      integer i, sca_id
      character aperture*10, filename*120
c
      dimension aperture(10)
      data aperture/'NRCA1_FULL','NRCA2_FULL','NRCA3_FULL',
     &     'NRCA4_FULL','NRCA5_FULL','NRCB1_FULL','NRCB2_FULL',
     &     'NRCB3_FULL','NRCB4_FULL','NRCB5_FULL'/
    
      write(filename,10) aperture(sca_id-480)
 10   format(a10,'_wcs_parameters.dat')
      print 20,filename
 20   format(a120)
      open(1,file=filename)
      do i = 1, 10
         read(1,*)
      end do
      read(1,*) x_sci_ref
      read(1,*) y_sci_ref
      read(1,*) x_sci_scale
      read(1,*) y_sci_scale
      read(1,*) v2_ref
      read(1,*) v3_ref
      read(1,*) v3idl_angle
      read(1,*) v3_parity
      close(1)
c      print *,'read_sca_wcs',sca_id,x_sci_ref, y_sci_ref,
c     &     x_sci_scale, y_sci_scale,
c     &     v2_ref, v3_ref, v3_parity, v3idl_angle
      return
      end
      
