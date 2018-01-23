      subroutine star_grid(nstars, ra1, dec1, ra2, dec2, sfaint, dmag,
     *     ra, dec, mag, max_stars, nbands)
c
c     create an exponential spiral distribution
c
c
c     (C) Copyright 2015 Christopher N. A. Willmer
c     Steward Observatory, University of Arizona
c     cnaw@as.arizona.edu 
c     2015-02-02
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

      implicit none
      integer nstars, k, nnn, nbands, i, j, max_stars
      double precision ra1, dec1, ra2, dec2, x0, y0, xx, yy, mag, pi,
     *     sfaint, dmag, colour, ra, dec, dra, ddec, dx
c
      dimension ra(max_stars), dec(max_stars), mag(max_stars, 54)
c
      dx =   dsqrt(dble(nstars))
      x0   = (ra1+ra2)/2.d0
      y0   = (dec1+dec2)/2.d0
      dra  = (ra2 - ra1) /dx
      ddec = (dec2-dec1)/dx
c
      print *,'star_grid: nbands', nbands
      nstars = 0
      do k = 1, idint(dx)
         do j = 1, idint(dx)
            nstars = nstars + 1
            dec(nstars) = dec1 + (k-1)*ddec
            ra(nstars)  = ra1  + (j-1)*dra
            do i = 1, nbands
               colour    = (i-1)*0.2d0
               mag(nstars,i)  = sfaint  !- (nstars-i+1) * dmag ! - colour
            end do
         end do
      end do
      return
      end
