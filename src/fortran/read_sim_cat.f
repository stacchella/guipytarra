c      character filename*80
c      filename = 'candels_nircat.cat'
c      call read_sim_cat(filename)
c      stop
c      end
c
      subroutine read_sim_cat(filename,ngal)
      implicit none
      double precision ra, dec, z, magnitude, nsersic,
     *     ellipticity, re, theta, flux_ratio
      double precision tra, tdec, tz, tmagnitude, tnsersic,
     *     ttheta, tlambda, semi_major, semi_minor,
     *     tra1, tdec1, tra2,tdec2, q, cosdec, angle, dra, ddec
c
      integer max_objects, nfilters, nsub
      integer ngal, ncomponents,i, j, nc
c
      character filename*80,line*100
c     
      parameter (max_objects=50000, nfilters=54, nsub=4)
c     
      dimension ra(max_objects), dec(max_objects), z(max_objects),
     *     magnitude(max_objects,nfilters), ncomponents(max_objects), 
     *     nsersic(max_objects,nsub),ellipticity(max_objects,nsub), 
     *     re(max_objects, nsub), theta(max_objects,nsub),
     *     flux_ratio(max_objects, nsub)
c     
      common /galaxy/ra, dec, z, magnitude, nsersic, ellipticity, re,
     *     theta, flux_ratio, ncomponents
c
      q = dacos(-1.0d0)/180.d0
c
      open(1,file=filename)
      read(1,*)
      ngal = 0
      open(2,file='cat.reg')
      do i = 1, max_objects
         read(1, *, err= 30, end=110) j, tra, tdec, tmagnitude,
     *        tz, semi_major, semi_minor, ttheta, tnsersic
         go to 50
 30      continue
         print *,'skipping '
         print 40, j, tra, tdec, tmagnitude,
     *        tz, semi_major, semi_minor, ttheta, tnsersic
 40      format(i8,9(1x,f15.6))
         go to 90
 50      continue
c         write(line,60) tra, tdec, tmagnitude
c 60      format('fk5;point(',f12.6,',',f12.6,') # point=box',
c     *        ' color=cyan text={',f5.2,'}')
         cosdec = dcos(tdec*q)
         angle  = ttheta
         dra    = 0.5d0*semi_major * dcos(ttheta*q)/cosdec
         dra    = dra/3600.d0
         ddec   = 0.5d0*semi_major * dsin(ttheta*q)
         ddec   = ddec/3600.d0
         tra1 = tra - dra
         tra2 = tra + dra

         tdec1 = tdec - ddec
         tdec2 = tdec + ddec

         write(line,60) tra1, tdec1, tra2, tdec2,tmagnitude
 60      format('fk5;line(',f12.6,',',f12.6,',',f12.6,',',f12.6,
     *        ' # line = 0 0 color=black text={',f5.2,'}')
         write(2,70) line
 70      format(a100)
c         if(tmagnitude.gt.25.d0) go to 90
         ngal = ngal + 1
         ra(ngal)           = tra
         dec(ngal)          = tdec
         z(ngal)            = tz
c     axial_ratio = 1.d0-ellipticity
         do j = 1, nfilters
            magnitude(ngal,j)   = tmagnitude !
         end do
         ncomponents(ngal)   = 1
         ellipticity(ngal,1) = semi_minor/semi_major - 1.d0
         nsersic(ngal,1)     = tnsersic
         re(ngal,1)          = dsqrt(semi_major*semi_minor)
         flux_ratio(ngal,1)  = 1.d0
         theta(ngal,1)       = ttheta + 90.d0
 90      continue
      end do
 110  close(1)
      print 120, ngal, filename
 120  format('read ',i6,' objects from',/,a80)
      close(2)
      return
      end
