c
c-----------------------------------------------------------------------
c
      subroutine read_psf_list(psf_file)
c
      implicit none
c
      integer i, nf, nfilters
      character psf_file*180, temp_psf*80
      character path_guitarra*100
c
      parameter (nfilters = 54)
c
      dimension psf_file(nfilters)
c
c      open(1,file='psf.list')
      call getenv('GUITARRA_HOME',path_guitarra)
      open(1,file=path_guitarra(1:len_trim(path_guitarra))
     +     //'data/psf.list')
      nf = 0
      do i = 1, nfilters
         read(1,10,end=100,err=20) temp_psf
 10      format(a160)
         go to 40
 20      print 30, i, temp_psf
 30      format('problem with ',i3,2x,a160)
         stop
 40      nf = nf + 1
         psf_file(nf) = path_guitarra(1:len_trim(path_guitarra))
     +                  //temp_psf
c
c     for now just repeat the file name
c
c         nf = nf + 1
c         psf_file(nf) = temp ! module B
      end do
 100  close(1)
      print *,'number of PSFs read', nf
c
      return
      end
