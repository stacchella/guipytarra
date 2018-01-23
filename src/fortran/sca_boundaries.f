      subroutine sca_boundaries(option, v2_min, v2_max, v3_min, v3_max)
c
c     load approximate boundaries of individual SCAs using OSIM coordinates
c     derived from CV2
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
c
      implicit none
      integer option
      double precision v2_min, v2_max, v3_min, v3_max
c
c     SCA limits
c
c     rough limits for SCAs (OSIM coordinates)
c     481 : 57 < V2 < 96, -356 < V3 < -318
c     482 : 58 < V2 < 96, -314 < V3 < -276
c     483 : 14 < V2 < 53, -356 < V3 < -318
c     484 : 14 < V2 < 53, -314 < V3 < -276
c     485 : 14 < V2 < 96, -355 < V3 < -276
c
c     486 : -97 < V2 < -58, -313 < V3 < -274
c     487 : -97 < V2 < -57, -355 < V3 < -317
c     488 : -54 < V2 < -16, -313 < V3 < -274
c     489 : -54 < V2 < -16, -355 < V3 < -317
c     490 : -97 < V2 < -16, -355 < V3 < -274
c      
      if(option.eq.0) then
         v2_min = -110.d0
         v2_max =  110.d0
         v3_min = -360.d0
         v3_max = -270.d0
      end if
c
      if(option.eq.481) then
         v2_min =   57.d0
         v2_max =   96.d0
         v3_min = -356.d0
         v3_max = -318.d0
      end if
c
      if(option.eq.482) then
         v2_min =   58.d0
         v2_max =   98.d0
         v3_min = -315.d0
         v3_max = -274.d0
      end if
c
      if(option.eq.483) then
         v2_min =   14.d0
         v2_max =   53.d0
         v3_min = -356.d0
         v3_max = -318.d0
      end if
c
      if(option.eq.484) then
         v2_min =   14.d0
         v2_max =   53.d0
         v3_min = -314.d0
         v3_max = -276.d0
      end if
c
      if(option.eq.485 .or.option.eq.1) then
         v2_min =   14.d0
         v2_max =   96.d0
         v3_min = -356.d0
         v3_max = -276.d0
      end if
c
      if(option.eq.486) then
         v2_min =  -97.d0
         v2_max =  -58.d0
         v3_min = -313.d0
         v3_max = -274.d0
      end if
c
      if(option.eq.487) then
         v2_min =  -97.d0
         v2_max =  -58.d0
         v3_min = -355.d0
         v3_max = -317.d0
      end if
c
      if(option.eq.488) then
         v2_min =  -54.d0
         v2_max =  -16.d0
         v3_min = -313.d0
         v3_max = -274.d0
      end if
c
      if(option.eq.489) then
         v2_min =  -54.d0
         v2_max =  -16.d0
         v3_min = -355.d0
         v3_max = -317.d0
      end if
c
      if(option.eq.490 .or. option.eq.2) then
         v2_min =  -97.d0
         v2_max =  -16.d0
         v3_min = -355.d0
         v3_max = -274.d0
      end if
c
      return
      end
