      subroutine wcs_keywords_new(sca_id, ra_sca, dec_sca, 
     &     x_sci_ref, y_sci_ref, x_sci_scale, y_sci_scale,
     &     v3idl_angle, pa_degrees, verbose)
c     &     v2_ref, v3_ref, v3_parity, v3idl_angle,
c     *     ra_dithered, dec_dithered, pa_degrees, verbose,
c     *     ,x_sci_ref, y_sci_ref,
c     &     x_sci_scale, y_sci_scale,
c     &     v2_ref, v3_ref, v3_parity, v3idl_angle)
c
c     (C) Copyright 2015 Christopher N. A. Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-10
c     Corrected the rotation matrix
c     2016-10-12
c     changed for SIAF distortions
c     2017-03-31 
c     
c------------------------------------------------------------------------------
c     This routine is free software; you can redistribute it and/or modify it
c     under the terms of the GNU General Public License as published by the Free
c     Software Foundation; either version 3 of the License, or (at your option)
c     any later version.
c     This program is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c     GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      double precision osim_scale, xc, yc, q, sw_scale, lw_scale,
     *     cospa, sinpa, cosx, sinx, cosy, siny, ra_sca, dec_sca,
     *     x_sca, y_sca, ra_dithered, dec_dithered, pa_degrees,
     *     xshift, yshift, xmag, ymag, xrot, yrot,
     *     oxshift, oyshift, oxmag, oymag, oxrot, oyrot,
     *     cosdec,
     *     x_sci_ref, y_sci_ref,
     &     x_sci_scale, y_sci_scale,
     &     v2_ref, v3_ref, v3_parity, v3idl_angle
c
      double precision a11, a12, a21, a22
      double precision
     *     cd1_1, cd1_2, cd2_1, cd2_2, crval1, crval2, crpix1, crpix2,
     *     cdelt1, cdelt2, equinox,
     *     pc1_1, pc1_2, pc2_1, pc2_2
c      real cd1_1, cd1_2, cd2_1, cd2_2, crval1, crval2, crpix1, crpix2,
c     *     cdelt1, cdelt2, equinox
      integer sca_id, i, verbose
c     
c      dimension xshift(10), yshift(10), xmag(10), ymag(10), xrot(10),
c     *     yrot(10)
c      dimension oxshift(10), oyshift(10), oxmag(10), oymag(10),
c     *     oxrot(10), oyrot(10)
c     
      common /wcs/ equinox, crpix1, crpix2, crval1, crval2,
     *     cdelt1,cdelt2, cd1_1, cd1_2, cd2_1, cd2_2,
     *     pc1_1, pc1_2, pc2_1, pc2_2
c      common /transform/ xshift, yshift, xmag, ymag, xrot, yrot
c      common /otransform/ oxshift, oyshift, oxmag, oymag, oxrot, oyrot
c     
c     data xc, yc /0.0d0, -315.6d0/
c     data osim_scale/1.59879d0/ 
c     
      q  = dacos(-1.d0)/180.0d0
c     
c-----------------------------------------------------------------------
c     
      cospa = dcos((pa_degrees)*q)
      sinpa = dsin((pa_degrees)*q)
c     
      cosdec = dcos(dec_sca*q)
c     
      cosx  = dcos(v3idl_angle*q)
      sinx  = dsin(v3idl_angle*q)
c     
      cosy  = dcos(v3idl_angle*q)
      siny  = dsin(v3idl_angle*q)
c
c     These transformations give very small residuals when
c     matching fields at different rotations. They
c     are systematic and of the order of 0.2" in RA and Dec.,
c     but infinitely better than the previous rotation matrix.
c
      a11 =  x_sci_scale * cosx * cospa - x_sci_scale * sinx * sinpa
      a12 =  y_sci_scale * siny * cospa + y_sci_scale * cosy * sinpa
      a21 = -x_sci_scale * cosx * sinpa - x_sci_scale * sinx * cospa
      a22 = -y_sci_scale * siny * sinpa + y_sci_scale * cosy * cospa
c
c      a11 =  x_sci_scale * cosx * cospa - x_sci_scale * sinx * sinpa
c      a12 =  y_sci_scale * siny * cospa + y_sci_scale * cosy * sinpa
c      a21 = -x_sci_scale * cosx * sinpa - x_sci_scale * sinx * cospa
c      a22 = -y_sci_scale * siny * sinpa + y_sci_scale * cosy * cospa
c
      a11 = -a11
      a21 = -a21
c
      cd1_1  =  a11 / 3600.d0
      cd1_2  =  a12 / 3600.d0
      cd2_1  =  a21 / 3600.d0
      cd2_2  =  a22 / 3600.d0
c
      equinox   = 2000.00000000
      crval1    = ra_sca        ! CRVALn  refer to the centre of NIRCam
      crval2    = dec_sca       ! at ra_dithered, dec_dithered
      cdelt1    = x_sci_scale/3600.d0 ! scale in degrees per pixel
      cdelt2    = y_sci_scale/3600.d0 
      crpix1    = x_sci_ref     ! CRPIXn are the SCA coordinates for
      crpix2    = y_sci_ref     ! XC, YC
c     
      if(verbose.gt.1) then
         print *,'wcs_keywords_new ', ra_sca, dec_sca, 
     &        x_sci_ref, y_sci_ref
         print *,'wcs_keywords_new ', x_sci_scale, y_sci_scale, 
     &       v3idl_angle
         print *,'wcs_keywords_new ',crval1, crval2, crpix1, crpix2,
     &        cdelt1, cdelt2
      end if
c     
      return
      end
