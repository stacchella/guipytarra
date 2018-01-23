      subroutine king_profile(rcore_arcmin, n_central, background, nobj, 
     *     profile_seed, x0, y0, x_r, y_r, nnn)
c
c     create a random distribution following a radial King profile
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
      double precision int_profile, n_central, zbqlu01, rcore, 
     *     background, x0, y0, x_r, y_r, two_pi, sigma0,q,pi,
     *     sum, radius, profile, fmax, fz, xx, yy, rad, angle, cosdec,
     *     rcore_arcmin, rmax, dr
      integer profile_seed, nobj, i, k, nnn
      dimension x_r(nnn), y_r(nnn)
      dimension radius(nnn), profile(nnn), int_profile(nnn)
c
      pi     = dacos(-1.0d0)
      q      = pi/180.d0
      two_pi = 2.d0 * pi
c
c     Create the King profile
c     create samples out to 10 core radii
c
      rcore    = rcore_arcmin/60.d0
      rmax     = rcore*10.d0 
      dr       = rmax/nnn
      sigma0   = 2.D0* n_central * rcore
      sum = 0.0d0
      do i = 1, nnn
         radius(i) = (i-1)* dr
         profile(i) = sigma0/(1.d0+(radius(i)/rcore)**2) + background
c     
c     Cumulative distribution
c
          sum = sum + profile(i)
          int_profile(i)  = sum
       end do
c
c     Create a random sample with random angles and radii 
c     described by the King profile
c
       fmax = int_profile(nnn)
       k = 0 
  10   continue
       if(k.ge.nobj) return
       fz    = zbqlu01(profile_seed) * fmax
       call linear_interpolation(nnn,int_profile, radius, fz, rad, nnn)
       angle  =  zbqlu01(profile_seed) * two_pi
       yy = y0 + rad * dsin(angle) 
       cosdec = dcos(yy*q)
       xx = x0 + rad * dcos(angle)/cosdec
       k = k + 1
       x_r(k) = xx
       y_r(k) = yy
       go to 10
       end
