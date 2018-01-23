c
c----------------------------------------------------------------------
c
      double precision function schechter (mag)
      implicit double precision (a-h,o-z)
      double precision mag, mstar
      common /schechter_params/ mstar, alpha, phistar, p_evol, q_evol
      xm         = 0.4D0 * ( mstar - mag)
      xm         = 10.d0**xm
      term1      = 0.4D0 * dlog10(10.d0) * phistar
      term2      = xm **(alpha+1.d0)
      term3      = exp(-xm)
      schechter  = term1 * term2 * term3
      return
      end
c
