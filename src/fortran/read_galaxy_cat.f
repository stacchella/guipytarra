c
c     Read galaxy catalogue with ID, RA, DEC, mag, redshift,
c     semi major, semi minor, theta, list of mags
c
      subroutine read_galaxy_cat(filename, filters_in_cat)
      implicit none
      double precision ra_galaxies, dec_galaxies, z, magnitude, nsersic,
     *     ellipticity, re, theta, flux_ratio
      double precision tra, tdec, tz, tmagnitude, tnsersic,
     *     ttheta, tlambda, semi_major, semi_minor,
     *     tra1, tdec1, tra2,tdec2, q, cosdec, angle, dra, ddec,
     *     abmag
c
      integer max_objects, nfilters, nsub, indx
      integer ngal, ncomponents, i, j, nc, l, filters_in_cat, id
c
      character filename*180,line*100
c     
      parameter (max_objects=50000, nfilters=54, nsub=4)
c     
      dimension ra_galaxies(max_objects), dec_galaxies(max_objects), 
     *     z(max_objects), magnitude(max_objects,nfilters),
     *     ncomponents(max_objects), nsersic(max_objects,nsub),
     *     ellipticity(max_objects,nsub), 
     *     re(max_objects, nsub), theta(max_objects,nsub),
     *     flux_ratio(max_objects, nsub), abmag(nfilters),
     *     id(max_objects)
c     
      common /galaxy/ra_galaxies, dec_galaxies, z, magnitude,
     *     nsersic, ellipticity, re, theta, flux_ratio, ncomponents,
     *     id, ngal
c
      q = dacos(-1.0d0)/180.d0
c
      print *, filters_in_cat
      print *, filename
      open(1,file=filename)
      read(1,*)
      ngal = 0
      open(2,file='cat.reg')
      do i = 1, max_objects
         read(1, *, err= 30, end=110) l, tra, tdec, tmagnitude,
     *        tz, semi_major, semi_minor, ttheta, tnsersic,
     *        (abmag(j), j = 1, filters_in_cat)
         go to 50
 30      continue
         print *,'skipping '
         print 40, l, tra, tdec, tmagnitude,
     *        tz, semi_major, semi_minor, ttheta, tnsersic,
     *        (abmag(j), j=1, filters_in_cat)
 40      format(i8,20(1x,f15.6))
         go to 90
 50      continue
c
         cosdec = dcos(tdec*q)
         angle  = ttheta
         dra    = 0.5d0*semi_major * dcos(ttheta*q)/cosdec
         dra    = dra/3600.d0
         ddec   = 0.5d0*semi_major * dsin(ttheta*q)
         ddec   = ddec/3600.d0
         tra1 = tra - dra
         tra2 = tra + dra
c
         tdec1 = tdec - ddec
         tdec2 = tdec + ddec
c
c         if(tmagnitude.gt.25.d0) go to 90
c
         ngal = ngal + 1
         id(ngal)           = l
         ra_galaxies(ngal)  = tra
         dec_galaxies(ngal) = tdec
         z(ngal)            = tz
c     axial_ratio = 1.d0-ellipticity
c
c     As the source catalogue may use a different set of filters,
c     this allows using the correct column for the filter
c
         do j = 1, filters_in_cat
            magnitude(ngal,j)   = abmag(j) !
         end do
c
         ncomponents(ngal)   = 1
         ellipticity(ngal,1) = 1.d0 - semi_minor/semi_major 
         nsersic(ngal,1)     = tnsersic
         re(ngal,1)          = dsqrt(semi_major*semi_minor)
         flux_ratio(ngal,1)  = 1.d0
         theta(ngal,1)       = ttheta !+ 90.d0
c
         write(line,60) tra1, tdec1, tra2, tdec2
     *        ,magnitude(ngal,filters_in_cat)
 60      format('fk5;line(',f12.6,',',f12.6,',',f12.6,',',f12.6,
     *        ' # line = 0 0 color=black text={',f5.2,'}')
         write(2,70) line
 70      format(a100)
 90      continue
c
      end do
 110  close(1)
      print 120, ngal, filename
 120  format('read ',i6,' objects from',/,a80)
      close(2)
      return
      end
