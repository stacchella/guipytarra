c-----------------------------------------------------------------------
      subroutine getkeywords(unit,name,x0,dx,status)
      integer unit, status
      character name*20, comment*80
      status=0
      call ftgkys(unit,"IRAFNAME", name, comment,status)
      if(status .ne.0) call printerror(status)
      if(status.eq.202) then
         status= 0
         call ftgkys(unit,"OBJECT", name, comment,status)
c         if(status .ne.0) call printerror(status)
         status= 0
      end if
      call ftgkye(unit,"CRVAL1",x0, comment,status)
      call ftgkye(unit,"CD1_1",dx, comment,status)
      if(status.ne.0) then
         status= 0
         call ftgkye(unit,"CDELT1",dx, comment,status)
         if(status .ne.0) call printerror(status)
      end if
c      print 10, name, x0, dx
 10   format(1x,'OBJECT   ', a20,' Wl1 ', f10.2, ' Dwl/Dpix ', f8.2 )

      if(status .ne.0) call printerror(status)
      return
      end
