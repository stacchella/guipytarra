c
c----------------------------------------------------------------------
c
      double precision function int_schechter(z)
      implicit double precision (a-h,o-z)
      double precision mstar, mmin, mmax, eps, final_acc
      integer m1, mfinal
      parameter (m1=14, eps = 1.d-16)
      common /schechter_params/ mstar, alpha, phistar, p_evol, q_evol
      common /sample_params/   bright, faint, apm_bright, apm_faint, 
     *     zmin, zmax,solid_angle
      external schechter
      dm     =  distance_modulus(z)
      
      mmin   = apm_bright - dm 
      mmax   = apm_faint  - dm
      
      if(mmin .lt. bright) mmin = bright
      if(mmax .gt. faint ) mmax = faint
      
      if(mmin .ge. bright .and. mmax .le. faint) then
         int_schechter =  hrvint(schechter,mmin, mmax,
     *        m1, eps, final_acc, mfinal)
c         call qromb(schechter, mmin, mmax, int_schechter)
      else 
         int_schechter = 0.d0
      end if
      return
      end
