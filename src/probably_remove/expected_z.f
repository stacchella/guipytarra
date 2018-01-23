c----------------------------------------------------------------------
c
      double precision function expected_z (z)
      implicit double precision (a-h,o-z)
      double precision int_schechter
      common /sample_params/   bright, faint, apm_bright, apm_faint, 
     *     zmin, zmax,solid_angle
      external int_schechter
      dvdz       = vol_elm(z, solid_angle)
c      print *,'expected_z ', z, dvdz, solid_angle
      expected_z = int_schechter(z) * dvdz
      return
      end
c
