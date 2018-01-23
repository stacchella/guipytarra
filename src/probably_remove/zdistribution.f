c
c----------------------------------------------------------------------
c
c
c Calculate the redshift distribution for MC sampling
c
      subroutine zdistribution(zdist, expected, int_zdist, nz)
      implicit double precision (a-h,o-z)
      double precision int_zdist, eps, final_acc
      integer m1, mfinal
      parameter(nnn=2048, eps = 1.d-16, m1=10)
      dimension zdist(nnn), int_zdist(nnn), expected(nnn)
      common /sample_params/   bright, faint, apm_bright, apm_faint, 
     *     zmin, zmax,solid_angle
      external expected_z

      print *,'zdistribution: zmin, zmax, apm_bright, apm_faint, nz',
     *     zmin, zmax, apm_bright, apm_faint, nz

      dz   = (zmax-zmin)/(nz-1.d0)
      ytot = 0.d0
      zdist(1)     = getz(5.d0)
      int_zdist(1) = 0.0d0
      expected(1)  = 0.0d0
      do i = 2, nz
         zdist(i) = dble(i-1)* dz
         expected(i) = hrvint(expected_z, zdist(i-1), zdist(i),
     *        m1, eps, final_acc, mfinal)
         if(expected(i) .lt. 0.d0)  expected(i) = 0.d0
         ytot         = ytot + expected(i)
         int_zdist(i) = ytot
         print *,'zdistribution', i, zdist(i), expected(i)
      end do
      print *,'zdistribution: total number of galaxies ', ytot
      return
      end

