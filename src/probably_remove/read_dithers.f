c
c     Read dither tables calculated by Jay Anderson 2011 JWST-STScI-001738.
c     The sub-pixel dithers are converted into arc seconds
c
      subroutine read_dithers(smallx, smally, sca_x, sca_y, full3x,
     *     full3y, tight3x, tight3y, full6x, full6y,full9x, full9y,
     *     full15x, full15y, full21x,full21y, full45x,full45y,
     *     intrasca_small_x,intrasca_small_y,
     *     intrasca_medium_x,intrasca_medium_y,
     *     intrasca_large_x,intrasca_large_y,
     *     gen_x, gen_y)
      implicit none
      integer i, j,k
      double precision smallx, smally, sca_x, sca_y, full3x,
     *     full3y, tight3x, tight3y, full6x, full6y,full9x, full9y,
     *     full15x, full15y, full21x,full21y, full45x,full45y,
     *     intrasca_small_x,intrasca_small_y,
     *     intrasca_medium_x,intrasca_medium_y,
     *     intrasca_large_x,intrasca_large_y,
     *     gen_x, gen_y
      character dither_file*80
      dimension smallx(9,9), smally(9,9), sca_x(16), sca_y(16), 
     *     full3x(3), full3y(3), tight3x(3),tight3y(3), 
     *     full6x(6), full6y(6), full9x(9), full9y(9),
     *     full15x(15), full15y(15), full21x(21),full21y(21),
     *     full45x(45),full45y(45),
     *     intrasca_small_x(25),intrasca_small_y(25),
     *     intrasca_medium_x(25),intrasca_medium_y(25),
     *     intrasca_large_x(25),intrasca_large_y(25),
     *     gen_x(64), gen_y(64)
      dither_file = 'dithers.dat'
      open(1,file=dither_file)
c
c     small dithers (SW pixels)
c
c     All dithers are positive in the S.I. coordinate system. 
c     Convert into moves in the V2, V3 system which has an orientation
c     such that
c     V2 = - sci_x
c     V3 = + sci_y
c
c     In the case of sub-pixel dithers, the values are in the
c     detector reference frame and are measured in SW pixels (0.0317")
c     For 
c     A5 these would need to be sci_x = -smallx, sci_y =  smally
c     b5 these would need to be sci_x =  smallx, sci_y = -smally
c     Moves as implemented mean that for PA = 0, a +NIRCam move is a -sky move
c     and a +NIRCam move is a + sky move
c
      read(1,*)
      do i = 1, 9
         read(1,*) j, (smallx(j,k), smally(j,k),k=1,j)
c         print *, j, (smallx(j,k), smally(j,k),k=1,j)
         do k = 1, j
            smallx(j,k) = -smallx(j,k) * 0.0317d0
            smally(j,k) = -smally(j,k) * 0.0317d0
         end do
c         print *, j, (smallx(j,k), smally(j,k),k=1,j)
      end do
c
c
c     intra-module (arc sec)
      read(1,*)
      read(1,*)
      do i = 1, 16
         read(1,*) j, sca_x(i), sca_y(i)
         sca_x(i) = -sca_x(i)
      end do
c
c     full field tile
c
      read(1,*)
      read(1,*)
      do i = 1, 3
         read(1,*) j, full3x(i), full3y(i)
         full3x(i) = -full3x(i)
      end do
      
c
c     full field tight
c
      read(1,*)
      read(1,*)
      do i = 1, 3
         read(1,*) j, tight3x(i), tight3y(i)
         tight3x(i) = -tight3x(i)
      end do
      
c
c     6 point full field
c
      read(1,*)
      read(1,*)
      do i = 1, 6
         read(1,*) j, full6x(i), full6y(i)
         full6x(i) = -full6x(i)
      end do
c
c     9 point full field
c
      read(1,*)
      read(1,*)
      do i = 1, 9
         read(1,*) j, full9x(i), full9y(i)
         full9x(i) = -full9x(i)
      end do
c
c     read intra-sca dithers
c
      read(1,*)
      read(1,*)
      read(1,*)
      do i = 1, 25
         read(1, *) j, intrasca_small_x(i), intrasca_small_y(i),
     *        intrasca_medium_x(i), intrasca_medium_y(i),
     *        intrasca_large_x(i), intrasca_large_y(i)
         intrasca_small_x(i)  = -intrasca_small_x(i)
         intrasca_medium_x(i) = -intrasca_medium_x(i)
         intrasca_large_x(i)  = -intrasca_large_x(i)
      end do
c
c     read general secondary dither pattern
c     units of SW pixels
c
      read(1,*)
      read(1,*)
      do i = 1, 64
         read(1,*) j, gen_x(i), gen_y(i)
         gen_x(i)  = - gen_x(i)
      end do
c
c     read full field 15-point
c
      read(1,*)
      read(1,*)
      do i = 1, 15
         read(1,*) j, full15x(i), full15y(i)
         full15x(i) = -full15x(i)
      end do
c
c     read full field 21-point
c
      read(1,*)
      read(1,*)
      do i = 1, 21
         read(1,*) j, full21x(i), full21y(i)
         full21x(i) = -full21x(i)
      end do
c
c     read full field 45-point
c
      read(1,*)
      read(1,*)
      do i = 1, 45
         read(1,*) j, full45x(i), full45y(i)
         full45x(i) = -full45x(i)
      end do
      
      close(1)
      return
      end

      
