c
c-----------------------------------------------------------------------
c
c
      subroutine find_sca_coords(x_osim, y_osim, module, sca_sw, x_sw, 
     *     y_sw, sca_lw, x_lw, y_lw)
      implicit none
      double precision x_osim, y_osim, x_sw, y_sw, x_lw, y_lw,
     *     d1, d2, d3, d4
      integer module, sca_sw, sca_lw, index_lw, in_lw, in_sw, ll, indx, 
     *     id1, id2, id3, id4
c
c     For a given OSIM position, find on what module and SCAs it falls 
c     - if LW and SW, SW or LW only, or none at all. 
c     
      module = 0
      if(x_osim.lt.0d0) then 
         module = 2
         index_lw =  10
      else
         module = 1
         index_lw = 5
      end if
c     
c     find LW coordinates for this position
c     
      in_lw  = 0
      sca_lw = 0
      call sca_coords_from_osim(index_lw, x_osim, y_osim, x_lw, y_lw) 
      if(idnint(x_lw) .ge. 5 
     *     .and.  idnint(x_lw) .le. 2044
     *     .and.  idnint(y_lw) .ge. 5 
     *     .and.  idnint(y_lw) .le. 2044) then 
         in_lw = 1
         if(index_lw.eq.5) sca_lw = 485
         if(index_lw.eq.10) sca_lw = 490
      end if
c      print *, 'LW',sca_lw, x_osim, y_osim, x_lw, y_lw,in_lw
c     
c     find SW coordinates for this position and check that it falls
c     within the boundaries of a SW SCA
c     
      in_sw = 0
      do ll= 1, 9
         indx = (module-1)*5 + ll
         if(ll.ne.5) then
            call sca_coords_from_osim(ll, x_osim, y_osim, x_sw, y_sw) 
c            print *, ' '
c            if(idnint(x_sw)  .ge. 5 .and. idnint(x_sw) .le. 2044) then
c               d1 = 0.d0
c               id1 = ll + 480
c               d2 = 0.d0
c               id2 = ll + 480
c            else
c               if(idnint(x_sw) .lt. 5)  then
c                  d1  = dabs(x_sw - 5)
c                  id1 = ll+480
c               end if
c               if(idnint(x_sw) .gt. 2044)  then
c                  d2  = dabs(x_sw-2044)
c                  id2 = ll+480
c               end if
c            end if
c            
c            if(idnint(y_sw) .ge. 5 .and. idnint(y_sw) .le. 2044) then
c               d3 = 0.d0
c               id3 = ll + 480
c               d4 = 0.d0
c               id4 = ll + 480
c            else
c               if(idnint(y_sw) .lt. 5)  then
c                  d3  = dabs(y_sw - 5)
c                  id3 = ll+480
c               end if
c               if(idnint(y_sw) .gt. 2044)  then
c                  d4  = dabs(y_sw-2044)
c                  id4 = ll+480
c               end if
c            end if
c            print 390, id1, d1, id2, d2, id3, d3, id4, d4
c 390        format(4(1x,i4,1x, f10.2))
c
            if(  idnint(x_sw) .ge. 5    .and. 
     *           idnint(x_sw) .le. 2044 .and.
     *           idnint(y_sw) .ge. 5    .and. 
     *           idnint(y_sw) .le. 2044) then 
               in_sw = 1 
c               print *, 'SW',480+ll, x_osim, y_osim, x_sw, y_sw,in_sw
               if(in_lw.eq.1 .and.in_sw.eq.1) then
                  sca_sw = 480 + ll
                  go to 405
               end if
            end if
         end if
      end do
 405  continue
      if(x_lw .eq. -1.0d0 .or. y_lw .eq.-1.0d0) sca_lw = 0
      if(in_lw.eq.0) then
         x_lw = -1.d0
         y_lw = -1.d0
         sca_lw  =  0
      end if
      if(in_sw.eq.0) then
         x_sw = -1.d0
         y_sw = -1.d0
         sca_sw  =  0
      end if
      if(in_lw.eq.0 .and. in_sw.eq.0) module = 0
      return
      end
