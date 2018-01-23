c-----------------------------------------------------------------------
c
      subroutine make_fake_cat(
     &     npts, id, ra, dec, hmag, hmag_err, 
     &     semi_major, semi_major_err, semi_minor, semi_minor_err,
     &     n_sersic, n_sersic_err, pa, pa_err, z, nnn)
      
      implicit none
      integer npts, id, nnn
      double precision ra, dec, hmag, hmag_err, 
     &     semi_major, semi_major_err,semi_minor, semi_minor_err,
     &     n_sersic, n_sersic_err,  pa, pa_err, z
       double precision tra, tdec, tz, tmagnitude, tnsersic,
     *     ttheta, tlambda, tsemi_major, tsemi_minor,
     *     tra1, tdec1, tra2,tdec2, q, cosdec, angle, dra, ddec,
     *     fit, fit_a
      double precision ref_mag, abmag,magnitude
c
      integer ref_filter, nfilters, label
      integer max_objects, maxfilters, nsub
      integer ngal, ncomponents,i, j, nc, l
      integer debug, in_range, in_range_a
c
      character line*100
c     
      parameter (max_objects=50000, maxfilters=54, nsub=4)
c     
c      dimension ra(max_objects), dec(max_objects), z(max_objects),
c     *     magnitude(max_objects), ncomponents(max_objects), 
c     *     nsersic(max_objects),
c     *     re(max_objects, nsub), theta(max_objects,nsub),
c     *     flux_ratio(max_objects, nsub), label(maxfilters)
      dimension id(nnn), ra(nnn), dec(nnn), hmag(nnn), hmag_err(nnn), 
     &     semi_major(nnn), semi_major_err(nnn), 
     &     semi_minor(nnn), semi_minor_err(nnn), 
     &     n_sersic(nnn), n_sersic_err(nnn), pa(nnn), pa_err(nnn),
     &     z(nnn)
      dimension abmag(maxfilters), magnitude(maxfilters), 
     &     label(maxfilters)
c     
c      common /galaxy/ra, dec, z, magnitude, nsersic, ellipticity, re,
c     *     theta, flux_ratio, ncomponents
c
      q = dacos(-1.0d0)/180.d0
      nfilters = 20
c
      do j = 1, nfilters
         label(j) = j
      end do
c
      debug = 0
      open(42,file='candels_with_fake_mag.cat')
      write(42,20) (label(j), j = 1, nfilters)
 20   format('# id   ra    dec   ref_mag  z  semi_major semi_minor',
     &     '  pa  nsersic ',54(3x,'mag_',i2.2))
      ngal = 0
      ref_filter = 3
      open(2,file='fake_mag_cat.reg')

      do i = 1, npts
         tra         = ra(i)
         tdec        = dec(i)
         tsemi_major = semi_major(i)
         tsemi_minor = semi_minor(i)
         ttheta      = pa(i)
         tmagnitude  = hmag(i)
         tnsersic    = n_sersic(i)
c
         cosdec = dcos(tdec*q)
         angle  = ttheta
         dra    = 0.5d0*tsemi_major * dcos(ttheta*q)/cosdec
         dra    = dra/3600.d0
         ddec   = 0.5d0*tsemi_major * dsin(ttheta*q)
         ddec   = ddec/3600.d0
         tra1 = tra - dra
         tra2 = tra + dra

         tdec1 = tdec - ddec
         tdec2 = tdec + ddec
c
         write(line,55) tra, tdec
 55      format('fk5;point(',f12.6,',',f12.6,
     *        ') # point=circle color=yellow')
         write(2,70) line
         write(line,60) tra1, tdec1, tra2, tdec2 !,l
 60      format('fk5;line(',f12.6,',',f12.6,',',f12.6,',',f12.6,
     *        ' # line = 0 0 color=magenta')! text={',i5.5,'}')
         write(2,70) line
 70      format(a100)
c     call fake_spectrum_mags(tz, ref_mag, ref_filter,
c     *     abmag, nfilters)
         do j = 1, nfilters
           magnitude(j)   = tmagnitude 
c     magnitude(ngal,j)   = abmag(j) !
         end do

         write(42, 80) id(i), tra, tdec, tmagnitude,
     *        tz, tsemi_major, tsemi_minor, ttheta, tnsersic,
     *        (magnitude(j), j = 1, nfilters)
 80      format(i11,1x, 1x, f11.6 ,1x, f11.6, 2x, f6.3, 2x,f8.1,
     *        4(1x,f7.3),54(1x,f8.4))
 90      continue
      end do
 110  continue
      close(42)
      print 120, npts
 120  format('read ',i6)
      close(2)
      return
      end
