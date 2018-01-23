      subroutine globular_cluster (nstars, ra0, dec0, dr,
     *     ra, dec, mag, sbright, rcore_arcsec, seed, nnn, nbands)
c
c     create a random globular cluster following a King profile
c     with a LF based on NGC2419
c
c
c     (C) Copyright 2015 Christopher N. A. Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-02
c     added core radius as an input parameter
c     2015-05-28
c
c------------------------------------------------------------------------------
c This routine is free software; you can redistribute it and/or modify it
c under the terms of the GNU General Public License as published by the Free
c Software Foundation; either version 3 of the License, or (at your option) any
c later version.
c This program is distributed in the hope that it will be useful,
c but WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c GNU General Public License for more details.
c-----------------------------------------------------------------------------
c
      implicit none
      integer nstars, i, k, nnn, nbands, nbins, nobj, seed, j,n0,n1
      double precision ra, dec, chi, r, ra0, dec0, xx, yy, mag,
     *     pi, dr, total, scale_lf, ran, int_lf, magbins, n_central, 
     *     background, rcore_arcmin, lf, zbqlu01, colour, sbright,
     8     rcore_arcsec
c     
      dimension ra(nnn), dec(nnn), mag(nnn,nbands), int_lf(nnn),
     *     magbins(nnn), lf(nnn)
c     
c     create a star cluster with a King profile
c     units degrees
c     
      do i = 1, nnn
         int_lf(i) = 0.0d0
      end do
      n_central     = 10000.
      background    = 300.
      nobj          = nstars
      rcore_arcmin  = rcore_arcsec/60.0d0 ! core radius of 10 arc sec
c     
      call king_profile(rcore_arcmin, n_central, background, nobj, 
     *     seed, ra0, dec0, ra, dec, nnn)
c     
c     Create a magnitude distribution for stars:
c     use ngc2419 luminosity function for magnitude generation
c     
      call read_gclf(magbins, lf, nbins, nnn)
c     
c     calculate the integral luminosity function
c     
      scale_lf = 0.d0
      do k = 1, nbins
         scale_lf  = scale_lf + lf(k)
         int_lf(k) = scale_lf
      end do
c     
      do i = 1, nstars
         ran        = zbqlu01(seed) * scale_lf
         n0       = 0
         n1       = 0
         do j = 1, nnn
            if(int_lf(j).ne.0) n0 = n0 + 1
            if(magbins(j).ne.0) n1 = n1 + 1
         end do
         print *, i, n0, n1, nnn
         call linear_interpolation(nbins, int_lf, magbins, ran,
     *        mag(i,1), nnn)
         mag(i,1) = mag(i,1) + sbright
         do k = 2, nbands
            colour = (k-1)* 0.20d0
            mag(i,k) = mag(i,1) - colour
         end do
c     print *, ra(i), dec(i), mag(i,1),mag(i,2), mag(i,3)
      end do
c     
      return
      end
