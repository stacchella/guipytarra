      double precision function dust_attenuation(wl,A_V)
      implicit double precision (a-h,o-z)
      double precision k_lambda
c
      rv = 4.05d0
      if(wl.lt.0.12d0) k_lambda = 100.d0
c    
      if(wl.le.0.63d0) then
c      if(wl.ge.0.12d0 .and. wl.le.0.63d0) then
         term = -2.156d0 + 1.509d0/wl-0.198/wl**2+0.011/wl**3
         k_lambda = 2.659*term + rv
      end if
c
c      if(wl.ge.0.63d0 .and. wl.le.2.20d0) then
      if(wl.ge.0.63d0) then
         term = -1.857d0 + 1.040d0/wl 
         k_lambda = 2.659*term +rv
      end if
c
      es_BV = A_V/rv
c     
c
c     from Fitzpatrick 2004, ASP Conf. Ser 309, 33
c
c      if(wl.gt.2.2d0) then
c         k_lambda = 0.63d0*rv -0.84
c         es_BV = k_lambda / wl**1.84d0
c      end if
c
      expo = 0.4d0*es_BV * k_lambda
      dust_attenuation = dexp(-expo)
      return
      end
