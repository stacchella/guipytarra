      subroutine exponential_spiral(nstars, x0, y0, sfaint, dmag, dr,
     *     phase, ra, dec, mag, nnn, nbands)
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
      integer nstars, k, nnn, nbands, i
      double precision ra, dec, chi, r, x0, y0, xx, yy, mag, pi,
     *     sfaint, phase, dmag, dr, colour,q, cosdec
c
      dimension ra(nnn), dec(nnn), mag(nnn, nbands)
c
      pi = dacos(-1.d0)
      q  = pi/180.d0
      chi = 0.0d0
      do k = 1, nstars
         do i = 1, nbands
            colour    = (i-1)*0.2d0
            mag(k,i)  = sfaint - k * dmag! - colour
c            mag(k,i)  = sfaint - k * dmag - colour
c            mag(k,i)  = sfaint -6.50d0
         end do
         chi       =  chi + (pi/7.0d0)
         r         =  dr + chi/2000.d0
c     
         if(k.eq.1) then
            xx = x0
            yy = y0
         else
            yy        =  y0 + r * dsin(chi+phase)
            cosdec    =  dcos(yy*q)
            xx        =  x0 + r * dcos(chi+phase)/cosdec
         end if
         ra(k) = xx
         dec(k) = yy
      end do
      return
      end
