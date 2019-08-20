cc Copyright (C) 2017: Travis Askham, Leslie Greengard, Zydrunas Gimbutas
cc email: askhamwhat@gmail.com      
cc 
cc This software is being released under a modified FreeBSD license
cc (see licenses folder in home directory). 
      
C**********************************************************************
C
C     This subroutine sets an adaptive tree for testing purposes.
C     Note that a 1 level tree only has colleague interactions
c
c     The 2 level adaptive tree has small to big and big to small
c     interactions
c
c     The 3 level tree adds child to parent, parent to child, and
c     more planewave interactions (compared to the 2 level which
c     has only small to big far interactions).
c
c     The 4 level tree is available for good measure.
C
C     INPUT:
C     
C     IMODE:  type of test tree to build 
c
c             IMODE = 2 gives a 2 level adaptive tree
C             IMODE = 3 gives a 3 level adaptive tree
C             Otherwise gives a 4 level adaptive tree
C
C
C     OUTPUT:
C
C     LEVELBOX is an array defining the level of each box
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     NBOXES is the total number of boxes
C
C     NLEV is the finest level
C
C**********************************************************************
      subroutine lrt2d_set(levelbox,icolbox,irowbox,nboxes,nlev,
     1     ichildbox, iparentbox, nblevel, istartlev, iboxlev,imode)
      implicit none
c-----Global variables
      integer levelbox(1)
      integer icolbox(1), irowbox(1)
      integer ichildbox(4,*), iparentbox(*), nblevel(0:1)
      integer istartlev(0:1), iboxlev(*)
      integer nboxes, nlev, imode

c-----Local variables
      integer ibox, i

      if (imode .eq. 2) go to 1000
      if (imode .eq. 3) go to 2000

C     Mock tree (4 level adaptive):

      nboxes = 17
      nlev = 4

      do i = 1,17
         iboxlev(i) = i
      enddo

      istartlev(0) = 1
      nblevel(0) = 1

      ibox = 1
      levelbox(ibox) = 0
      icolbox(ibox) = 1
      irowbox(ibox) = 1
      iparentbox(ibox) = -1
      ichildbox(1,ibox) = 4
      ichildbox(2,ibox) = 5
      ichildbox(3,ibox) = 3
      ichildbox(4,ibox) = 2


      istartlev(1) = 2
      nblevel(1) = 4

      ibox = 2
      levelbox(ibox) = 1
      icolbox(ibox) = 1
      irowbox(ibox) = 1
      iparentbox(ibox) = 1
      ichildbox(1,ibox) = 8
      ichildbox(2,ibox) = 9
      ichildbox(3,ibox) = 7
      ichildbox(4,ibox) = 6


      ibox = 3
      levelbox(ibox) = 1
      icolbox(ibox) = 2
      irowbox(ibox) = 1
      iparentbox(ibox) = 1
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 4
      levelbox(ibox) = 1
      icolbox(ibox) = 1
      irowbox(ibox) = 2
      iparentbox(ibox) = 1
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 5
      levelbox(ibox) = 1
      icolbox(ibox) = 2
      irowbox(ibox) = 2
      iparentbox(ibox) = 1
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      istartlev(2) = 6
      nblevel(2) = 4

      ibox = 6
      levelbox(ibox) = 2
      icolbox(ibox) = 1
      irowbox(ibox) = 1
      iparentbox(ibox) = 2
      ichildbox(1,ibox) = 12
      ichildbox(2,ibox) = 13
      ichildbox(3,ibox) = 11
      ichildbox(4,ibox) = 10


      ibox = 7
      levelbox(ibox) = 2
      icolbox(ibox) = 2
      irowbox(ibox) = 1
      iparentbox(ibox) = 2
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 8
      levelbox(ibox) = 2
      icolbox(ibox) = 1
      irowbox(ibox) = 2
      iparentbox(ibox) = 2
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 9 
      levelbox(ibox) = 2
      icolbox(ibox) = 2
      irowbox(ibox) = 2
      iparentbox(ibox) = 2
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      istartlev(3) = 10
      nblevel(3) = 4

      ibox = 10
      levelbox(ibox) = 3
      icolbox(ibox) = 1
      irowbox(ibox) = 1
      iparentbox(ibox) = 6
      ichildbox(1,ibox) = 16
      ichildbox(2,ibox) = 17
      ichildbox(3,ibox) = 15
      ichildbox(4,ibox) = 14


      ibox = 11
      levelbox(ibox) = 3
      icolbox(ibox) = 2
      irowbox(ibox) = 1
      iparentbox(ibox) = 6
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 12
      levelbox(ibox) = 3
      icolbox(ibox) = 1
      irowbox(ibox) = 2
      iparentbox(ibox) = 6
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 13
      levelbox(ibox) = 3
      icolbox(ibox) = 2
      irowbox(ibox) = 2
      iparentbox(ibox) = 6
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      istartlev(4) = 14
      nblevel(4) = 4

      ibox = 14
      levelbox(ibox) = 4
      icolbox(ibox) = 1
      irowbox(ibox) = 1
      iparentbox(ibox) = 10
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 15
      levelbox(ibox) = 4
      icolbox(ibox) = 2
      irowbox(ibox) = 1
      iparentbox(ibox) = 10
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 16
      levelbox(ibox) = 4
      icolbox(ibox) = 1
      irowbox(ibox) = 2
      iparentbox(ibox) = 10
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 17
      levelbox(ibox) = 4
      icolbox(ibox) = 2
      irowbox(ibox) = 2
      iparentbox(ibox) = 10
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      return

c     Another mock tree (2 level, adaptive)

 1000 continue

      nboxes = 9
      nlev = 2

      do i = 1,9
         iboxlev(i) = i
      enddo

      istartlev(0) = 1
      nblevel(0) = 1

      ibox = 1
      levelbox(ibox) = 0
      icolbox(ibox) = 1
      irowbox(ibox) = 1
      iparentbox(ibox) = -1
      ichildbox(1,ibox) = 4
      ichildbox(2,ibox) = 5
      ichildbox(3,ibox) = 3
      ichildbox(4,ibox) = 2


      istartlev(1) = 2
      nblevel(1) = 4

      ibox = 2
      levelbox(ibox) = 1
      icolbox(ibox) = 1
      irowbox(ibox) = 1
      iparentbox(ibox) = 1
      ichildbox(1,ibox) = 8
      ichildbox(2,ibox) = 9
      ichildbox(3,ibox) = 7
      ichildbox(4,ibox) = 6


      ibox = 3
      levelbox(ibox) = 1
      icolbox(ibox) = 2
      irowbox(ibox) = 1
      iparentbox(ibox) = 1
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 4
      levelbox(ibox) = 1
      icolbox(ibox) = 1
      irowbox(ibox) = 2
      iparentbox(ibox) = 1
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 5
      levelbox(ibox) = 1
      icolbox(ibox) = 2
      irowbox(ibox) = 2
      iparentbox(ibox) = 1
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1



      istartlev(2) = 6
      nblevel(2) = 4

      ibox = 6
      levelbox(ibox) = 2
      icolbox(ibox) = 1
      irowbox(ibox) = 1
      iparentbox(ibox) = 2
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 7
      levelbox(ibox) = 2
      icolbox(ibox) = 2
      irowbox(ibox) = 1
      iparentbox(ibox) = 2
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 8
      levelbox(ibox) = 2
      icolbox(ibox) = 1
      irowbox(ibox) = 2
      iparentbox(ibox) = 2
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 9 
      levelbox(ibox) = 2
      icolbox(ibox) = 2
      irowbox(ibox) = 2
      iparentbox(ibox) = 2
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1



      return


 2000 continue

C     Mock tree (3 level adaptive):

      nboxes = 25
      nlev = 3

      do i = 1,13
         iboxlev(i) = i
      enddo

      istartlev(0) = 1
      nblevel(0) = 1

      ibox = 1
      levelbox(ibox) = 0
      icolbox(ibox) = 1
      irowbox(ibox) = 1
      iparentbox(ibox) = -1
      ichildbox(1,ibox) = 4
      ichildbox(2,ibox) = 5
      ichildbox(3,ibox) = 3
      ichildbox(4,ibox) = 2


      istartlev(1) = 2
      nblevel(1) = 4

      ibox = 2
      levelbox(ibox) = 1
      icolbox(ibox) = 1
      irowbox(ibox) = 1
      iparentbox(ibox) = 1
      ichildbox(1,ibox) = 8
      ichildbox(2,ibox) = 9
      ichildbox(3,ibox) = 7
      ichildbox(4,ibox) = 6


      ibox = 3
      levelbox(ibox) = 1
      icolbox(ibox) = 2
      irowbox(ibox) = 1
      iparentbox(ibox) = 1
      ichildbox(1,ibox) = 12
      ichildbox(2,ibox) = 13
      ichildbox(3,ibox) = 11
      ichildbox(4,ibox) = 10


      ibox = 4
      levelbox(ibox) = 1
      icolbox(ibox) = 1
      irowbox(ibox) = 2
      iparentbox(ibox) = 1
      ichildbox(1,ibox) = 16
      ichildbox(2,ibox) = 17
      ichildbox(3,ibox) = 15
      ichildbox(4,ibox) = 14


      ibox = 5
      levelbox(ibox) = 1
      icolbox(ibox) = 2
      irowbox(ibox) = 2
      iparentbox(ibox) = 1
      ichildbox(1,ibox) = 20
      ichildbox(2,ibox) = 21
      ichildbox(3,ibox) = 19
      ichildbox(4,ibox) = 18


      istartlev(2) = 6
      nblevel(2) = 16

      ibox = 6
      levelbox(ibox) = 2
      icolbox(ibox) = 1
      irowbox(ibox) = 1
      iparentbox(ibox) = 2
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 7
      levelbox(ibox) = 2
      icolbox(ibox) = 2
      irowbox(ibox) = 1
      iparentbox(ibox) = 2
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 8
      levelbox(ibox) = 2
      icolbox(ibox) = 1
      irowbox(ibox) = 2
      iparentbox(ibox) = 2
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 9 
      levelbox(ibox) = 2
      icolbox(ibox) = 2
      irowbox(ibox) = 2
      iparentbox(ibox) = 2
      ichildbox(1,ibox) = 24
      ichildbox(2,ibox) = 25
      ichildbox(3,ibox) = 23
      ichildbox(4,ibox) = 22



      ibox = 10
      levelbox(ibox) = 2
      icolbox(ibox) = 3
      irowbox(ibox) = 1
      iparentbox(ibox) = 3
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 11
      levelbox(ibox) = 2
      icolbox(ibox) = 4
      irowbox(ibox) = 1
      iparentbox(ibox) = 3
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 12
      levelbox(ibox) = 2
      icolbox(ibox) = 3
      irowbox(ibox) = 2
      iparentbox(ibox) = 3
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 13
      levelbox(ibox) = 2
      icolbox(ibox) = 4
      irowbox(ibox) = 2
      iparentbox(ibox) = 3
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1

      ibox = 14
      levelbox(ibox) = 2
      icolbox(ibox) = 1
      irowbox(ibox) = 3
      iparentbox(ibox) = 4
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 15
      levelbox(ibox) = 2
      icolbox(ibox) = 2
      irowbox(ibox) = 3
      iparentbox(ibox) = 4
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 16
      levelbox(ibox) = 2
      icolbox(ibox) = 1
      irowbox(ibox) = 4
      iparentbox(ibox) = 4
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 17
      levelbox(ibox) = 2
      icolbox(ibox) = 2
      irowbox(ibox) = 4
      iparentbox(ibox) = 4
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1

      ibox = 18
      levelbox(ibox) = 2
      icolbox(ibox) = 3
      irowbox(ibox) = 3
      iparentbox(ibox) = 5
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 19
      levelbox(ibox) = 2
      icolbox(ibox) = 4
      irowbox(ibox) = 3
      iparentbox(ibox) = 5
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 20
      levelbox(ibox) = 2
      icolbox(ibox) = 3
      irowbox(ibox) = 4
      iparentbox(ibox) = 5
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 21
      levelbox(ibox) = 2
      icolbox(ibox) = 4
      irowbox(ibox) = 4
      iparentbox(ibox) = 5
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      istartlev(3) = 22
      nblevel(3) = 4

      ibox = 22
      levelbox(ibox) = 3
      icolbox(ibox) = 3
      irowbox(ibox) = 3
      iparentbox(ibox) = 9
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 23
      levelbox(ibox) = 3
      icolbox(ibox) = 4
      irowbox(ibox) = 3
      iparentbox(ibox) = 9
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 24
      levelbox(ibox) = 3
      icolbox(ibox) = 3
      irowbox(ibox) = 4
      iparentbox(ibox) = 9
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      ibox = 25
      levelbox(ibox) = 3
      icolbox(ibox) = 4
      irowbox(ibox) = 4
      iparentbox(ibox) = 9
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1


      return

      end

c**********************************************************************
C     The following subroutine is used to generate a uniform tree.
C     This is only for testing purposes.
C
C
C     INPUT:
C
C     NLEV is the finest level
C
C     OUTPUT:
C
C     LEVELBOX is an array determining the level of each box
C
C     ICOLBOX denotes the column of each box
C
C     IROWBOX denotes the row of each box
C
C     NBOXES is the total number of boxes
C
C
C**********************************************************************
      subroutine lrt2d_uni(levelbox,icolbox,irowbox,nboxes,nlev, 
     1     ichildbox, iparentbox, nblevel, istartlev, iboxlev)
      implicit none
c-----Global variables
      integer  levelbox(1)
      integer  icolbox(1), irowbox(1)
      integer  nboxes, nlev
      integer  ichildbox(4,1), iparentbox(1)
      integer  nblevel(0:1), istartlev(0:1), iboxlev(*)
c-----Local variables
      integer  i, ibox, j, k, l, nside
      integer  icol, irow
      integer  ic1, ic2, ic3, ic4
      integer  ilength,  ilev, istart, idim
      integer  mybox, istartc

c     create a uniform tree, given NLEV as the input:
C     (In this loop, initialize the parents
C     and children arrays to be -1.)
      ibox = 1
      nside = 1
      do i = 0, nlev
         do j = 1, nside
            do k = 1, nside
               levelbox(ibox) = i
               icolbox(ibox)  = j
               irowbox(ibox)  = k
               iparentbox(ibox) = -1
               do l = 1, 4
                  ichildbox(l,ibox) = -1
               end do
               ibox = ibox + 1
            end do
         end do
         nside = 2 * nside
      end do
      nboxes = ibox - 1
C
C     initialize the nblevel and istartlev arrays:
C
      j = 4
      istartlev(0) = 1
      istartlev(1) = 2
      nblevel(0) = 1
      nblevel(1) = 4
      do i = 2, nlev
        istartlev(i) = istartlev(i-1) +  j
        j = 4*j
        nblevel(i) = j
      end do
C
C     set all of the parents and children:
C
      do i = 0, nlev - 1
         istart = istartlev(nlev-i-1)-1
         ilength = nblevel(nlev-i-1)
         ilev = nlev-i-1
         istartc = istartlev(nlev-i)-1
         idim = 2**ilev
         nside = 2*idim
         do j = 1, ilength
            icol = 1 + mod(j-1,idim )
            irow = 1 + (j-1)/idim
            mybox = istart+j
            ic1 = istartc + (2*icol) + (2*irow-2)*nside
            ic2 = istartc + (2*icol) + (2*irow-1)*nside
            ic3 = istartc + (2*icol-1) + (2*irow-1)*nside
            ic4 = istartc + (2*icol-1) + (2*irow-2)*nside
            ichildbox(1,mybox) = ic1
            ichildbox(2,mybox) = ic2
            ichildbox(3,mybox) = ic3
            ichildbox(4,mybox) = ic4
            iparentbox(ic1) = mybox
            iparentbox(ic2) = mybox
            iparentbox(ic3) = mybox
            iparentbox(ic4) = mybox
         enddo
      enddo

      do i = 1,nboxes
         iboxlev(i) = i
      enddo

      return
      end
C
c**********************************************************************
c     subroutine mkcolls
c**********************************************************************
C     The following subroutine is used to generate the colleagues
C     for all of the boxes in the tree structure.  If a colleague
C     doesn't exist it is set to -1.  Each box has nine colleagues
C     and they are ordered as follows:
C
C                        7     8     9
C               
C                        4     5     6
C
C                        1     2     3
C
C     You are your own colleague number 5.
C     The algorithm used here is recursive and takes advantage of
C     the fact that your colleagues can only be the children of 
C     your parents colleagues.  There is no need to scan all of the
C     boxes.  
C
C     INPUT:
C
C     ICOLBOX denotes the column of each box
C     IROWBOX denotes the row of each box
C     NBOXES is the total number of boxes
C     NLEV is the finest level
C     IPARENTBOX denotes the parent of each box
C     ICHILDBOX denotes the four children of each box
C     NBLEVEL is the total number of boxes per level
C     IBOXLEV is the array in which the boxes are arranged
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C
C     OUTPUT:
C
C     ICOLLEAGBOX denotes the colleagues of a given box
C**********************************************************************
      subroutine lrt2d_mkcolls(icolbox,irowbox,icolleagbox,nboxes,nlev,
     2      iparentbox,ichildbox,nblevel,iboxlev,istartlev)
      implicit none
c-----Global variables
      integer icolleagbox(9,1)
      integer icolbox(1), irowbox(1)
      integer nboxes, nlev, iparentbox(1)
      integer ichildbox(4,1)
      integer nblevel(0:1), iboxlev(1), istartlev(0:1)
c-----Local variables
      integer colleague, partemp
      integer jcntr, ibox, itest
      integer icntr, ilev, j, l, nside
      integer irowtemp, icoltemp
      integer irowtest, icoltest
C
c     Initialize colleague number 5 to self and other colleagues to
C     -1.  -1 is used (in general) to indicate the colleague doesn't
C     exist.
C   
      do ibox = 1, nboxes
         icolleagbox(5,ibox) = ibox
         do j = 1, 4
            icolleagbox(j,ibox) = -1
         enddo
         do j = 6, 9
            icolleagbox(j,ibox) = -1
         enddo
      enddo
C
c     Scan through levels (ignoring coarsest level which has
C     no colleagues in the free space case.
C
      do ilev = 1, nlev
C
C      Scan through all boxes on each level.  For each
C      box, scan children of parent's colleagues and see if they
C      are in contact.
C      Each colleague (if it exists) is placed in the standard order
C
         do l = istartlev(ilev), istartlev(ilev) + nblevel(ilev) - 1
            ibox    = iboxlev(l)
            partemp = iparentbox(ibox)
C
c        IROWTEMP and ICOLTEMP denote the row and column of
C        the box under consideration.
C
            irowtemp = irowbox(ibox)
            icoltemp = icolbox(ibox)
            do 100 jcntr = 1, 9
C
c          COLLEAGUE denotes the colleague of the parent box.
C
               !write(*,*) 'partemp ', partemp, jcntr
           COLLEAGUE = ICOLLEAGBOX(JCNTR,PARTEMP)
C
C          If parent's colleague doesn't exist
C          or is childless, skip it:
C
               if (colleague .lt. 0)goto 100
               if (ichildbox(1,colleague) .lt. 0)goto 100
               do icntr = 1, 4
                  j = ichildbox(icntr,colleague)
C
C            IROWTEST, ICOLTEST used as row, col of test box.
C
                  irowtest = irowbox(j)
                  icoltest = icolbox(j)
                  if(irowtemp .eq. irowtest+1)then
                    if(icoltemp .eq. icoltest+1)then
                      icolleagbox(1,ibox) = j
                    elseif(icoltemp .eq. icoltest)then
                      icolleagbox(2,ibox) = j
                    elseif(icoltemp .eq. icoltest-1)then
                      icolleagbox(3,ibox) = j
                    endif
                  elseif(irowtemp .eq. irowtest)then
                    if(icoltemp .eq. icoltest+1)then
                      icolleagbox(4,ibox) = j
                    elseif(icoltemp .eq. icoltest-1)then
                      icolleagbox(6,ibox) = j
                    endif
                  elseif(irowtemp .eq. irowtest-1)then
                    if(icoltemp .eq. icoltest+1)then
                      icolleagbox(7,ibox) = j
                    elseif(icoltemp .eq. icoltest)then
                      icolleagbox(8,ibox) = j
                    elseif(icoltemp .eq. icoltest-1)then
                      icolleagbox(9,ibox) = j
                    endif
                  endif
               enddo 
100         continue
         enddo
      enddo 
      return
      end

c***********************************************************************
C      subroutine restriction
C***********************************************************************
C
C     determines whether or not a given tree satisfies 
C     the level restriction. 
C     If it doesn't, call FIXTREE to fix the tree.
C
C***********************************************************************
      subroutine lrt2d_restrict(levelbox,iparentbox,ichildbox,icolbox, 
     1             irowbox,icolleagbox,nboxes,nlev,
     2             nblevel,iboxlev,istartlev,ifixflag)
      implicit none
c-----Global variables
      integer  levelbox(1), icolleagbox(9,1)
      integer  iparentbox(1), ichildbox(4,1)
      integer  icolbox(1), irowbox(1)
      integer  nboxes, nlev
      integer  nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer  ifixflag
c-----Local variables
      integer  ichild(4),icoll(9), ibox
      integer  i, ipar, itest, j, nb
      integer  itemp
C
c     Sort boxes by level. (Takes place of "ladder" structure
C     in uniform case).  All boxes are sorted
C     into the array and placed in their proper places.
C
      call lrt2d_sortboxes(levelbox,nboxes,nlev,
     1           nblevel,iboxlev,istartlev)
C
c     First let's call a subroutine that will
C     generate all of the colleagues for each
C     box.  The colleagues are generated in the
C     correct order so there is no need to 'shuffle'
C     them later on.
C
      call lrt2d_mkcolls(icolbox,
     1       irowbox,icolleagbox,nboxes,nlev,
     2       iparentbox,ichildbox,nblevel,
     3       iboxlev, istartlev)
      do i = nlev, 2, -1
         do j = istartlev(i), istartlev(i) + nblevel(i) - 1
            ibox = iboxlev(j)
            ipar  = iparentbox(ibox)
            itest = iparentbox(ipar)
            icoll(1) = icolleagbox(1,itest)
            icoll(2) = icolleagbox(2,itest)
            icoll(3) = icolleagbox(3,itest)
            icoll(4) = icolleagbox(4,itest)
            icoll(5) = icolleagbox(5,itest)
            icoll(6) = icolleagbox(6,itest)
            icoll(7) = icolleagbox(7,itest)
            icoll(8) = icolleagbox(8,itest)
            icoll(9) = icolleagbox(9,itest)
            ichild(1) = ichildbox(1,itest)
            ichild(2) = ichildbox(2,itest)
            ichild(3) = ichildbox(3,itest)
            ichild(4) = ichildbox(4,itest)
            do nb = 1, 9
               itemp = icoll(nb)
               if(itemp .gt. 0 .and. ichildbox(1,itemp) .lt. 0)then
c             The neighboring box is not divided
C             we could have problems.
                 if (nb .eq. 1)then
                   if(ipar .eq. ichild(4))then
                       ifixflag = 1
                   end if
                 elseif (nb .eq. 2)then
                   if(ipar .eq. ichild(3) .or. ipar .eq. ichild(4))then
                       ifixflag = 1
                   end if
                 elseif (nb .eq. 3)then
                   if(ipar .eq. ichild(3))then
                       ifixflag = 1
                   end if
                 elseif (nb .eq. 4)then
                   if(ipar .eq. ichild(4) .or. ipar .eq. ichild(1))then
                       ifixflag = 1
                   end if
                 elseif (nb .eq. 6)then
                   if(ipar .eq. ichild(2) .or. ipar .eq. ichild(3))then
                       ifixflag = 1
                   end if
                 elseif (nb .eq. 7)then 
                   if(ipar .eq. ichild(1))then
                       ifixflag = 1
                   end if
                 elseif (nb .eq. 8)then
                   if(ipar .eq. ichild(1) .or. ipar .eq. ichild(2))then
                       ifixflag = 1
                   end if
                 elseif (nb .eq. 9)then
                   if(ipar .eq. ichild(2))then
                       ifixflag = 1
                   end if
                 endif
               endif
            enddo
         enddo
      enddo
      return
      end


c***********************************************************************
c      subroutine sortboxes
c***********************************************************************

C     sets up a structure that is analogous
C     to the 'LADDER' structure in the nonadaptive case.
C     It is just a way of organizing the boxes by level in one long
C     array and denoting where in the array the levels change.
C
C     INPUT:
C
C     LEVELBOX is an array determining the level of each box
C     NBOXES is the total number of boxes
C     NLEV is the finest level
C
C     OUTPUT:
C
C     NBLEVEL is the total number of boxes per level
C     IBOXLEV is the array in which the boxes are arranged
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C
C***********************************************************************
      subroutine lrt2d_sortboxes(levelbox,nboxes,nlev,
     1           nblevel,iboxlev,istartlev)
      implicit none
C-----Global variables
      integer levelbox(1), nboxes, nlev
      integer nblevel(0:1), iboxlev(1), istartlev(0:1)
C-----Local variables
      integer ncntr, i, j
c
      ncntr = 1
      do i = 0, nlev
         nblevel(i) = 0
         istartlev(i) = ncntr
         do j = 1, nboxes
            if(levelbox(j) .eq. i)then
              iboxlev(ncntr) = j
              ncntr = ncntr + 1
              nblevel(i) = nblevel(i) + 1
            endif
         enddo
      enddo
      return
      end


      subroutine lrt2d_mktpts_query(pts,pts_sort,isort,npts,
     1     maxnodes,zll,blength,nstart,maxiter,nboxes,nlev,ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     figure out storage needs
c
c     npts - the number of points
c     maxnodes - maximum number of nodes per box
c     pts(2,*) - coordinates of points
c     pts_sort(2,*) - storage of same length
c     isort(*) - integer storage of length npts
c     nstart - number of boxes to start with, usually a guess at
c              the number (20*npts isn't a bad idea)
c     maxiter - number of times to try doubling number of boxes
c     
c     OUTPUT
c
c     nboxes - number of boxes in final tree
c     nlev - number of levels in final tree
c     ier = 0 - normal execution
c         = 8 - final tree didn't have enough boxes 
c               with nstart*2**(maxiter-1)
c         = 4 - needed more than 100 levels
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$    use omp_lib      
      implicit none
c     global variables
      real *8 pts(2,*), pts_sort(2,*), zll(2), blength
      integer isort(*), npts, maxnodes, nstart, maxiter
      integer nboxes, nlev, ier
c     local variables
      integer maxboxes, maxlevel, ilevelbox, iicolbox, iirowbox
      integer iiparentbox, iichildbox, iiboxlev, inblevel
      integer iistartlev, iirregfl, liused
      integer i, ifixflag, ifenough, ier2
      integer, allocatable :: iworktemp(:), itemparray(:)
      integer, allocatable :: iflag(:), icolleagbox(:,:), icladder(:,:)

      ier = 0
 
c     carve out work array for tree storage, other integer arrays

      maxboxes = nstart
      maxlevel = 30
      ilevelbox = 21
      iicolbox = ilevelbox+maxboxes
      iirowbox = iicolbox+maxboxes
      iiparentbox = iirowbox+maxboxes
      iichildbox = iiparentbox+maxboxes
      iiboxlev = iichildbox+4*maxboxes
      inblevel = iiboxlev+maxboxes
      iistartlev = inblevel+maxlevel+10
      iirregfl = iistartlev+maxlevel+10
      liused = iirregfl+maxboxes

c     figure out how many boxes are needed in tree

      ifenough = 0

      do i = 1,maxiter

c         write(*,*) 'iter ', i, zll, blength
         
         allocate(iworktemp(liused))
         allocate(itemparray(maxboxes+20))
         allocate(iflag(maxboxes))
         allocate(icolleagbox(9,maxboxes))
         allocate(icladder(2,maxboxes))

         ier2 = 0
         call lrt2d_mktpts(iworktemp(ilevelbox),iworktemp(iicolbox),
     1        iworktemp(iirowbox),nboxes,nlev,iworktemp(iiparentbox),
     2        iworktemp(iichildbox),iworktemp(inblevel),
     3        iworktemp(iiboxlev),iworktemp(iistartlev),maxboxes,
     4        itemparray,maxlevel,pts,pts_sort,isort,icladder,
     5        npts,maxnodes,zll,blength,ier2)

c         write(*,*) 'ier2 ', ier2

c     see if there are more boxes than maxboxes

         if (ier2 .eq. 0 .or. ier2 .eq. 4) then

c     if maxboxes not exceeded, proceed to level-restriction of tree

            if (ier2 .eq. 4) then
               write(*,*) 'IN MKTPTS_QUERY '
               write(*,*) 'WARNING: maxlevels ', maxlevel, 'exceeded '
               write(*,*) 'maybe maxnodes is too small, '
               write(*,*) 'or points are too tightly spaced '
            endif

            ifixflag = 0
            call lrt2d_restrict(iworktemp(ilevelbox),
     1           iworktemp(iiparentbox),iworktemp(iichildbox),
     2           iworktemp(iicolbox),iworktemp(iirowbox),
     3           icolleagbox,nboxes,nlev,
     4           iworktemp(inblevel),iworktemp(iiboxlev),
     5           iworktemp(iistartlev),ifixflag)

            if (ifixflag.eq.1) then
c               write(*,*) 'FIXING TREE'
               call lrt2d_fix(iworktemp(ilevelbox),
     1              iworktemp(iiparentbox),iworktemp(iichildbox),
     2              iworktemp(iicolbox),iworktemp(iirowbox),
     3              icolleagbox,nboxes,nlev,
     4              iworktemp(inblevel),iworktemp(iiboxlev),
     5              iworktemp(iistartlev),
     6              iflag,maxboxes,itemparray)
c               write(*,*) 'DONE FIXING TREE '
            endif

            if (nboxes .le. maxboxes) then
c     maxboxes was big enough
               ifenough = 1
               goto 1000
            else
c     maxboxes exceeded, increase maxboxes
               maxboxes = 2*maxboxes
            endif
            
         else
c     maxboxes exceeded, increase maxboxes
            maxboxes = 2*maxboxes
         endif

         deallocate(iworktemp)
         deallocate(itemparray)
         deallocate(iflag)
         deallocate(icolleagbox)
         deallocate(icladder)

         ilevelbox = 21
         iicolbox = ilevelbox+maxboxes
         iirowbox = iicolbox+maxboxes
         iiparentbox = iirowbox+maxboxes
         iichildbox = iiparentbox+maxboxes
         iiboxlev = iichildbox+4*maxboxes
         inblevel = iiboxlev+maxboxes
         iistartlev = inblevel+maxlevel+1
         iirregfl = iistartlev+maxlevel+1
         liused = iirregfl+maxboxes
         
      enddo

c     after maxiter iterations, arrays still not big enough

      ier = 8
      return
      
 1000 continue

      if (ifenough .eq. 1 .and. ier2 .eq. 0) ier = 0
      if (ifenough .eq. 1 .and. ier2 .eq. 4) ier = 4

      return
      end
      

     
C
C
c***********************************************************************
C     subroutine lrt2d_mktpts
C***********************************************************************
C     generates a tree given a set of boundary pts XS, YS
C
C     The algorithm works by subdividing boxes until there are less
C     than MAXNODES points of the boundary in each box. 
C
C     THE TREE GENERATED BY THIS ALGORITHM MAY NOT SATISFY THE LEVEL 
C     RESTRICTION.  The routine FIXTREE will convert this data structure 
C     to one which is suitably level restricted and compatible
C     with ADAPFMM4.
C
C     INPUT:
C
C     MAXBOXES    denotes the maximum number of boxes allowed
C     ITEMPARRAY  a workspace array of dimension MAXBOXES
C     MAXLEVEL    denotes the deepest level allowed
C     XS, YS      arrays determining the boundary pts
C     NPTS        the number of pts in the XS, YS arrays
C     MAXNODES    the maximum number of boundary pts for each box
C
C     OUTPUT:
C
C     LEVELBOX    is an array determining the level of each box
C     ICOLBOX     denotes the column of each box
C     IROWBOX     denotes the row of each box
C     NBOXES      is the total number of boxes
C     NLEV        is the finest level
C     IPARENTBOX  denotes the parent of each box
C     ICHILDBOX   denotes the four children of each box
C     NBLEVEL     is the total number of boxes per level
C     IBOXLEV     is an array in which the boxes are arranged
C     ISTARTLEV   is the pointer to where each level
C                 begins in the IBOXLEV array
C
C     IER         = 0 for normal execution
C                 = 4 maxlevel exceeded
C                 = 8 maxboxes exceeded
C                 = 16 failure in assigning points
C
C***********************************************************************
      subroutine lrt2d_mktpts(levelbox, icolbox, irowbox, nboxes, nlev,
     1     iparentbox, ichildbox, nblevel, iboxlev, istartlev,
     2     maxboxes, itemparray, maxlevel, src, srcsort, isrcsort,
     3     isrcladder, ns, maxnodes, zll, blength, ier)
c-----Global variables
      implicit none
      real *8 xlength
      integer levelbox(*), maxboxes, maxnodes, ns
      integer nlev, nboxes,  maxlevel
      integer icolbox(*), irowbox(*)
      integer iparentbox(*), ichildbox(4,*)
      integer nblevel(0:1), iboxlev(*), istartlev(0:1)
      integer itemparray(*)
      integer isrcsort(*), isrcladder(2,*)
      integer ier
      real *8 src(2,*), srcsort(2,*), zll(2), blength
c-----Local variables
      integer i, ibox, iflag
      integer j, istart, iend, l, jj
      integer levflag, nnodes
      integer ier2
      real *8 xstart, xincrem
      real *8, allocatable :: srctemp(:,:)
      integer, allocatable :: itemp(:)
      integer, allocatable :: itemp2(:)

c      write(*,*) 'zll ', zll, blength
c      write(*,*) ns

C
      allocate(srctemp(2,ns))
      allocate(itemp(ns))
      allocate(itemp2(ns))

      ier = 0

      do i = 0, maxlevel
         nblevel(i) = 0
         istartlev(i) = 2
      end do
      do i = 1, maxboxes
         iboxlev(i) = 0
         isrcladder(1,i) = 0
         isrcladder(2,i) = -1
      end do

      do i = 1,ns
         isrcsort(i) = i
         srcsort(1,i) = src(1,i)
         srcsort(2,i) = src(2,i)
      enddo
C
c     First set parameters for the single unit cell at level 0.
C
      ibox = 1
      levelbox(ibox) = 0
      icolbox(ibox) = 1
      irowbox(ibox) = 1
      iparentbox(ibox) = -1
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1

      isrcladder(1,ibox) = 1
      isrcladder(2,ibox) = ns

      nboxes = 1
      nlev = 0
C
c     We also need to initialize the adaptive 'Ladder'
C     structures to the correct initial values:
C
      nblevel(0) = 1
      istartlev(0) = 1
      iboxlev(1) = 1
      xlength = blength
      do i = 0, maxlevel - 1

         iflag = 0
         levflag = 0
         istart = istartlev(i)
         iend = nboxes
         istartlev(i+1) = nboxes+1

         do jj = istart, iend
            ibox = iboxlev(jj)
            nnodes = isrcladder(2,ibox)-isrcladder(1,ibox)+1
            if(nnodes .gt. maxnodes)then

c     box has too many nodes

c     check that there's room

               if (nboxes + 4 .gt. maxboxes) then
                  ier = 8
                  return
               endif

c     subdivide box and update ladder structure

               iboxlev(nboxes+1) = nboxes+1
               iboxlev(nboxes+2) = nboxes+2
               iboxlev(nboxes+3) = nboxes+3
               iboxlev(nboxes+4) = nboxes+4

               nblevel(i+1) = nblevel(i+1)+4

               call lrt2d_subdiv1(ibox,iparentbox,ichildbox,
     1              nboxes,irowbox,icolbox,levelbox)



c     assign points of this box to the appropriate child box

               call lrt2d_pts2child(ibox,isrcladder,isrcsort,
     1              itemp,itemp2,iparentbox,ichildbox,irowbox,icolbox,
     2              srcsort,srctemp,ier2,zll,xlength)

c     check if error from assigning points to child

               if (ier2 .ne. 0) then
                  ier = 16
                  return
               endif

c     increase level count if not already

               if(levflag .eq. 0)then
                  nlev = nlev + 1
                  levflag = 1
               endif
               iflag = 1
            endif
         enddo
         if(iflag .eq. 0)then
c     Nothing was divided, so exit loop 
            return
         endif
         xlength = xlength/2.0d0

c         write(*,*) i, xlength

      enddo

c     Check to see if max levels was too small

      istart = istartlev(maxlevel)
      iend = nboxes

c      do i = 0,maxlevel
c         write(*,*) nblevel(i)
c      enddo
         


      do jj = istart,iend
         ibox = iboxlev(jj)
         nnodes = isrcladder(2,ibox)-isrcladder(1,ibox)+1
c         write(*,*) 'ibox ', ibox, nnodes
         if (nnodes .gt. maxnodes) then
            ier = 4
            return
         endif
      enddo


      return
      end

C
C
c***********************************************************************
C     subroutine lrt2d_mktst
C***********************************************************************
C     generates a tree given a set of sources and targets
C
C     The algorithm works by subdividing boxes until there are less
C     than MAXNODES points (either source or target) in each box
C
C     THE TREE GENERATED BY THIS ALGORITHM MAY NOT SATISFY THE LEVEL 
C     RESTRICTION.  The routine FIXTREE will convert this data structure 
C     to one which is suitably level restricted and compatible
C     with ADAPFMM4.
C
C     INPUT:
C
C     MAXBOXES    denotes the maximum number of boxes allowed
C     ITEMPARRAY  a workspace array of dimension MAXBOXES
C     MAXLEVEL    denotes the deepest level allowed
C     XS, YS      arrays determining the boundary pts
C     NPTS        the number of pts in the XS, YS arrays
C     MAXNODES    the maximum number of boundary pts for each box
C
C     OUTPUT:
C
C     LEVELBOX    is an array determining the level of each box
C     ICOLBOX     denotes the column of each box
C     IROWBOX     denotes the row of each box
C     NBOXES      is the total number of boxes
C     NLEV        is the finest level
C     IPARENTBOX  denotes the parent of each box
C     ICHILDBOX   denotes the four children of each box
C     NBLEVEL     is the total number of boxes per level
C     IBOXLEV     is an array in which the boxes are arranged
C     ISTARTLEV   is the pointer to where each level
C                 begins in the IBOXLEV array
C
C     IER         = 0 for normal execution
C                 = 4 maxlevel exceeded
C                 = 8 maxboxes exceeded
C                 = 16 failure in assigning points
C
C***********************************************************************
      subroutine lrt2d_mktst(levelbox, icolbox, irowbox, nboxes, nlev,
     1     iparentbox, ichildbox, nblevel, iboxlev, istartlev,
     2     maxboxes, itemparray, maxlevel, src, srcsort, isrcsort,
     3     isrcladder, ns, targ, targsort, itargsort, itargladder, nt,
     4     maxnodes, zll, blength, ier)
c-----Global variables
      implicit none
      real *8 xlength
      integer levelbox(*), maxboxes, maxnodes, ns, nt
      integer nlev, nboxes,  maxlevel
      integer icolbox(*), irowbox(*)
      integer iparentbox(*), ichildbox(4,*)
      integer nblevel(0:1), iboxlev(*), istartlev(0:1)
      integer itemparray(*)
      integer isrcsort(*), isrcladder(2,*)
      integer itargsort(*), itargladder(2,*)
      integer ier
      real *8 src(2,*), srcsort(2,*), zll(2), blength
      real *8 targ(2,*), targsort(2,*)
c-----Local variables
      integer i, ibox, iflag, ntemp
      integer j, istart, iend, l, jj
      integer levflag, nnodes
      integer ier2
      real *8 xstart, xincrem
      real *8, allocatable :: srctemp(:,:)
      real *8, allocatable :: targtemp(:,:)      
      integer, allocatable :: itemp(:)
      integer, allocatable :: itemp2(:)

c      write(*,*) 'zll ', zll, blength
c      write(*,*) ns

C

      ntemp = ns
      if (nt .gt. ntemp) ntemp = nt
      
      allocate(srctemp(2,ns))
      allocate(targtemp(2,nt))
      allocate(itemp(ntemp))
      allocate(itemp2(ntemp))

      ier = 0

      do i = 0, maxlevel
         nblevel(i) = 0
         istartlev(i) = 2
      end do
      do i = 1, maxboxes
         iboxlev(i) = 0
         isrcladder(1,i) = 0
         isrcladder(2,i) = -1
         itargladder(1,i) = 0
         itargladder(2,i) = -1
      end do

      do i = 1,ns
         isrcsort(i) = i
         srcsort(1,i) = src(1,i)
         srcsort(2,i) = src(2,i)
      enddo
C
c     First set parameters for the single unit cell at level 0.
C
      ibox = 1
      levelbox(ibox) = 0
      icolbox(ibox) = 1
      irowbox(ibox) = 1
      iparentbox(ibox) = -1
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1

      isrcladder(1,ibox) = 1
      isrcladder(2,ibox) = ns
      itargladder(1,ibox) = 1
      itargladder(2,ibox) = nt

      nboxes = 1
      nlev = 0
C
c     We also need to initialize the adaptive 'Ladder'
C     structures to the correct initial values:
C
      nblevel(0) = 1
      istartlev(0) = 1
      iboxlev(1) = 1
      xlength = blength
      do i = 0, maxlevel - 1

         iflag = 0
         levflag = 0
         istart = istartlev(i)
         iend = nboxes
         istartlev(i+1) = nboxes+1

         do jj = istart, iend
            ibox = iboxlev(jj)
            nnodes = isrcladder(2,ibox)-isrcladder(1,ibox)
     1           + itargladder(2,ibox)-itargladder(1,ibox) + 2
            if(nnodes .gt. maxnodes)then

c     box has too many nodes

c     check that there's room

               if (nboxes + 4 .gt. maxboxes) then
                  ier = 8
                  return
               endif

c     subdivide box and update ladder structure

               iboxlev(nboxes+1) = nboxes+1
               iboxlev(nboxes+2) = nboxes+2
               iboxlev(nboxes+3) = nboxes+3
               iboxlev(nboxes+4) = nboxes+4

               nblevel(i+1) = nblevel(i+1)+4

               call lrt2d_subdiv1(ibox,iparentbox,ichildbox,
     1              nboxes,irowbox,icolbox,levelbox)

c     assign points of this box to the appropriate child box

               call lrt2d_pts2child(ibox,isrcladder,isrcsort,
     1              itemp,itemp2,iparentbox,ichildbox,irowbox,icolbox,
     2              srcsort,srctemp,ier2,zll,xlength)

               call lrt2d_pts2child(ibox,itargladder,itargsort,
     1              itemp,itemp2,iparentbox,ichildbox,irowbox,icolbox,
     2              targsort,targtemp,ier2,zll,xlength)

c     check if error from assigning points to child

               if (ier2 .ne. 0) then
                  ier = 16
                  return
               endif

c     increase level count if not already

               if(levflag .eq. 0)then
                  nlev = nlev + 1
                  levflag = 1
               endif
               iflag = 1
            endif
         enddo
         if(iflag .eq. 0)then
c     Nothing was divided, so exit loop 
            return
         endif
         xlength = xlength/2.0d0

c         write(*,*) i, xlength

      enddo

c     Check to see if max levels was too small

      istart = istartlev(maxlevel)
      iend = nboxes

c      do i = 0,maxlevel
c         write(*,*) nblevel(i)
c      enddo
         


      do jj = istart,iend
         ibox = iboxlev(jj)
         nnodes = isrcladder(2,ibox)-isrcladder(1,ibox)
     1        + itargladder(2,ibox)-itargladder(1,ibox) + 2
c         write(*,*) 'ibox ', ibox, nnodes
         if (nnodes .gt. maxnodes) then
            ier = 4
            return
         endif

      enddo


      return
      end

      subroutine lrt2d_pts2child(ibox,isrcladder,isrcsort,itemp,itemp2,
     1     iparentbox,ichildbox,irowbox,icolbox,srcsort,srctemp,ier,
     2     zll,xlength)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine takes the points located inside the box ibox
c     and assigns them to the appropriate children
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      integer ibox, isrcladder(2,*), isrcsort(*)
      integer itemp(*), itemp2(*)
      integer iparentbox(*), ichildbox(4,*), nboxes
      integer irowbox(*), icolbox(*)
      integer ier
      real *8 srcsort(2,*), srctemp(2,*), xlength, zll(2)
c     local variables
      integer i, irow, icol, istart, iend, j, istart1
      integer ic(4), ipt
      real *8 x, y, xc, yc
      logical xlt, ylt, ifany

      ier = 0

c     box has no children
      if (ichildbox(1,ibox) .lt. 0) then
         ier = 4
         return
      endif

      ic(1) = ichildbox(1,ibox)
      ic(2) = ichildbox(2,ibox)
      ic(3) = ichildbox(3,ibox)
      ic(4) = ichildbox(4,ibox)

      irow = irowbox(ibox)
      icol = icolbox(ibox)

c     center coordinates of ibox
      xc = zll(1)+xlength*(icol-0.5d0)
      yc = zll(2)+xlength*(irow-0.5d0)

c     location in isrcsort array of the pts in ibox
      istart = isrcladder(1,ibox)
      iend = isrcladder(2,ibox)

c     set up with default vals
      isrcladder(1,ic(1)) = 0
      isrcladder(2,ic(1)) = -1
      isrcladder(1,ic(2)) = 0
      isrcladder(2,ic(2)) = -1
      isrcladder(1,ic(3)) = 0
      isrcladder(2,ic(3)) = -1
      isrcladder(1,ic(4)) = 0
      isrcladder(2,ic(4)) = -1

c     there are no points in this box
      if (iend-istart .lt. 0) then
         return
      endif

c     save indeces in itemp
      do i = istart,iend
         itemp(i) = isrcsort(i)
      enddo

c     find out which child each point goes to, store in itemp2
      do i = istart,iend
         x = srcsort(1,i)
         y = srcsort(2,i)
         xlt = (x.lt.xc)
         ylt = (y.lt.yc)
         if (xlt .and. ylt) then
            itemp2(i) = 4
         else if (xlt .and. (.not. ylt)) then
            itemp2(i) = 1
         else if ((.not. xlt) .and. ylt) then
            itemp2(i) = 3
         else
            itemp2(i) = 2
         endif
      enddo
         
c     tracks position of pts in isrcsort array
      istart1 = istart

c     for each child, find corresponding pts
      do j = 1,4
         ifany = .false.
         do i = istart,iend
            if (itemp2(i) .eq. j) then
               isrcsort(istart1) = itemp(i)
               srctemp(1,istart1) = srcsort(1,i)
               srctemp(2,istart1) = srcsort(2,i)
               if ( .not. ifany) then
                  ifany = .true.
                  isrcladder(1,ic(j)) = istart1
               endif
               istart1 = istart1+1
            endif
         enddo
         if (ifany) isrcladder(2,ic(j)) = istart1-1
      enddo

c     copy source pts over

      do i = istart,iend
         srcsort(1,i) = srctemp(1,i)
         srcsort(2,i) = srctemp(2,i)
      enddo

c     check that all pts were assigned
      if (istart1-1 .ne. iend) ier = 8

               

      return
      end


C
c***********************************************************************
c     subroutine mklists
c***********************************************************************
C     computes all of the north, south,
C     east, and west interaction lists.  Because this is for the
C     adaptive case, all of these lists have to be computed in one
C     loop.  North and south take precedence over the east and
C     west and recall that processing is done at the parent level.
C    
C     lists for boxes at same level
C
C            |----|----|----|----|----|----|
C            | N  | N  | N  | N  | N  | N  |
C            |----|----|----|----|----|----|
C            | W  |N34*|N34 |N34 |N34*| E  |
C            |----|----|----|----|----|----|
C            | W  |W23 | 1  | 2  |E14 | E  |
C            |----|----|----|----|----|----|
C            | W  |W23 | 4  | 3  |E14 | E  |
C            |----|----|----|----|----|----|
C            | W  |S12*|S12 |S12 |S12*| E  |
C            |----|----|----|----|----|----|
C            | S  | S  | S  | S  | S  | S  |
C            |----|----|----|----|----|----|
C
C     lists for small to big far boxes 
C
C            |---------|---------|---------|
C            |         |         |         |
C            |  NBIG34 | NBIG34  | NBIG34  |
C            |  IWBIG2 |         | IEBIG1  |
C            |---------|----|----|---------|
C            |         | 1  | 2  |         |
C            | WBIG23  |----|----| EBIG14  |
C            |         | 4  | 3  |         |
C            |---------|----|----|---------|
C            |         |         |         |
C            | SBIG12  | SBIG12  | SBIG12  |
C            | IWBIG3  |         | IEBIG4  |
C            |---------|---------|---------|
C
C
C     N marks the NORTHALL list
C     S marks the SOUTHALL list
C     E marks the EASTALL list
C     W marks the WESTALL list
C     N34 marks boxes well separated from 3,4 but not 1,2
C     S12 marks boxes well separated from 1,2 but not 3,4
C     E14 marks boxes well separated from 1,4 but not 2,3
C     W23 marks boxes well separated from 2,3 but not 1,4
C
C     The boxes marked with an asterix are special cases. 
C     The lower left S12* box needs data from 3 -> W3 list.
C     The lower right S12* box needs data from 4 -> E4 list.
C     The upper left N34* box needs data from 2 -> W2 list.
C     The upper right N34* box needs data from 1 -> E1 list.
C
C     INPUT:
C
C     IBOX denotes the box being considered
C     LEVEL is the level of IBOX
C
C     ICOLLEAGBOX, IROWBOX, AND ICOLBOX define the tree
C
C     OUTPUT:
C
C     The naming convention for the lists is is as follows:
C
C     INALL is an array that denotes the boxes in the
C     north all list.  
C     NNALL is the number of boxes in the north all
C     list.  
C     IYNALL represents the corresponding offsets of the boxes
C     in the north all list.
C
C     The same convention is used for the south all, east all, west all,
C     north34, south12, east14, west23, west2, west3, east1,
C     and east4 lists.
C
C***********************************************************************
      subroutine lrt2d_mklists(ibox,inall,nnall,iynall, in34,nn34,iy34,
     1    isall,nsall,iysall,is12,ns12,iy12,ieall,neall,
     2    iyeall,ie14,ne14,iy14,iwall,nwall,iywall,
     3    iw23,nw23,iy23,iww3,iwy3,nww3,iww2,iwy2,nww2,
     4    iee4,iey4,nee4,iee1,iey1,nee1,inbig34,isbig12,iebig14,
     6    iwbig23,iebig4,iwbig3,iebig1,iwbig2,icolleagbox,ichildbox,
     8    icolbox, irowbox, iflageast, iflagwest, 
     9    iflagnorth,iflagsouth, localonoff)
      implicit none
C-----Global variables
      integer  inall(6),nnall,iynall(6)
      integer  isall(6),nsall,iysall(6)
      integer  ieall(4),neall,iyeall(4)
      integer  iwall(4),nwall,iywall(4)
      integer  in34(2),nn34,iy34(2)
      integer  is12(2),ns12,iy12(2)
      integer  iw23(2),nw23,iy23(2)
      integer  ie14(2),ne14,iy14(2)
      integer  iee4(1),nee4,iey4(1)
      integer  iee1(1),nee1,iey1(1)
      integer  iww2(1),nww2,iwy2(1)
      integer  iww3(1),nww3,iwy3(1)
      integer  localonoff(1)
      integer  ibox
      integer  inbig34(3), isbig12(3)
      integer  iebig14(1), iwbig23(1)
      integer  iebig4(1), iwbig3(1)
      integer  iebig1(1), iwbig2(1)
      integer  icolleagbox(9,1), ichildbox(4,1)
      integer  icolbox(1), irowbox(1)
      integer  iflageast(1), iflagnorth(1)
      integer  iflagwest(1), iflagsouth(1)
c-----Local variables
      integer  i, j, iout
      integer  ichild, ncntr, scntr, ecntr, wcntr
      integer  n34cntr, s12cntr
      integer  w23cntr, e14cntr
      integer  w4cntr, w2cntr
      integer  ee4cntr, e3cntr
      integer  icoltest, irowtest
      integer  icol, irow
C
C     Initially, set all list entries and offsets to zero.
C
      do j = 1, 6
        inall(j)  = 0
        iynall(j) = 0
        isall(j)  = 0
        iysall(j) = 0
      end do
      do j = 1, 4
        ieall(j)  = 0
        iyeall(j) = 0
        iwall(j)  = 0
        iywall(j) = 0
      end do
      do j = 1, 4
        isbig12(j) = -1
        inbig34(j) = -1
      end do
      do j = 1, 2
        in34(j) = 0
        iy34(j) = 0
        ie14(j) = 0
        iy14(j) = 0
        iw23(j) = 0
        iy23(j) = 0
        is12(j) = 0
        iy12(j) = 0
      end do
      iee4(1) = 0
      iey4(1) = 0
      iww3(1) = 0
      iwy3(1) = 0
      iee1(1) = 0
      iey1(1) = 0
      iww2(1) = 0
      iwy2(1) = 0
      iebig14(1) = -1
      iwbig23(1) = -1
      iebig4(1)  = -1
      iebig1(1)  = -1
      iwbig3(1)  = -1
      iwbig2(1)  = -1
C
C     All of the offsets are based from the
C     col,row of child 4 (lower left corner).
C
      icol = icolbox(ichildbox(4,ibox))
      irow = irowbox(ichildbox(4,ibox))
c
C     First do the free space case:
C     Set all of the counters to 1
      ncntr   =  1
      scntr   =  1
      ecntr   =  1
      wcntr   =  1
      n34cntr =  1
      s12cntr =  1
      w23cntr =  1
      e14cntr =  1
      w4cntr  =  1
      w2cntr  =  1
      ee4cntr  =  1
      e3cntr  =  1

c     First scan through all nine of the boxes colleagues
      do 100 i = 1, 9
         iout = icolleagbox(i,ibox)
c     Test to see if this colleague doesn't exist or is
C     childless, if so, skip it.
         if(iout .lt. 0)goto 100
         if(ichildbox(1,iout) .gt. 0)then
c       Scan all four of the colleagues children.
            DO J = 1, 4
C       ICOLTEST and IROWTEST represent the row and column
C       of the box being checked.
               ichild = ichildbox(j,iout)
               icoltest = icolbox(ichild)
               irowtest = irowbox(ichild)
               if(irowtest .eq. irow+3)then
                  inall(ncntr) = ichild
                  iynall(ncntr) = icoltest - icol
                  ncntr = ncntr + 1
               elseif(irowtest .eq. irow-2)then
                  isall(scntr) = ichild
                  iysall(scntr) = icoltest - icol
                  scntr = scntr + 1
               elseif(icoltest .eq. icol-2)then
                  iwall(wcntr) = ichild
                  iywall(wcntr) = irowtest - irow
                  wcntr = wcntr + 1
               elseif(icoltest .eq. icol+3)then
                  ieall(ecntr) = ichild
                  iyeall(ecntr) = irowtest - irow
                  ecntr = ecntr + 1
               elseif(irowtest .eq. irow+2)then
                  in34(n34cntr) = ichild
                  iy34(n34cntr) = icoltest - icol
                  n34cntr = n34cntr + 1
                  if(icoltest .eq. icol-1)then
                      iww2(w4cntr) = ichild
                      iwy2(w4cntr) = irowtest - irow
                      w4cntr = w4cntr + 1
                  endif
                  if(icoltest .eq. icol+2)then
                     iee1(e3cntr) = ichild
                     iey1(e3cntr) = irowtest - irow
                     e3cntr = e3cntr + 1
                  endif
               elseif(irowtest .eq. irow-1)then
                  is12(s12cntr) = ichild
                  iy12(s12cntr) = icoltest - icol
                  s12cntr = s12cntr + 1
                  if(icoltest .eq. icol-1)then
                     iww3(w2cntr) = ichild
                     iwy3(w2cntr) = irowtest - irow
                     w2cntr = w2cntr + 1
                  endif
                  if(icoltest .eq. icol+2)then
                     iee4(ee4cntr) = ichild
                     iey4(ee4cntr) = irowtest - irow
                     ee4cntr = ee4cntr + 1
                  endif
               elseif(icoltest .eq. icol-1)then
                  iw23(w23cntr) = ichild
                  iy23(w23cntr) = irowtest - irow
                  w23cntr = w23cntr + 1
               elseif(icoltest .eq. icol+2)then
                  ie14(e14cntr) = ichild
                  iy14(e14cntr) = irowtest - irow
                  e14cntr = e14cntr + 1
               endif
            enddo
          elseif(ichildbox(1,iout) .lt. 0)then
             if(i .eq. 1)then
               isbig12(1) = iout
               iwbig3(1) = iout
               iflagsouth(iout) = 1
               iflagwest(iout) = 1
             elseif(i .eq. 2)then
               isbig12(2) = iout
               iflagsouth(iout) = 1
             elseif(i .eq. 3)then
               isbig12(3) = iout
               iebig4(1) = iout
               iflageast(iout) = 1
               iflagsouth(iout) = 1
             elseif(i .eq. 4)then
               iwbig23(1) = iout
               iflagwest(iout) = 1
             elseif(i .eq. 6)then
               iebig14(1) = iout
               iflageast(iout) = 1
             elseif(i .eq. 7)then
               inbig34(1) = iout
               iwbig2(1) = iout
               iflagnorth(iout) = 1
               iflagwest(iout) = 1
             elseif(i .eq. 8)then
               inbig34(2) = iout
               iflagnorth(iout) = 1
             elseif(i .eq. 9)then
               inbig34(3) = iout
               iebig1(1) = iout
               iflagnorth(iout) = 1
               iflageast(iout) = 1
            endif
          endif
100   continue
      nnall = ncntr   -  1
      nsall = scntr   -  1
      neall = ecntr   -  1
      nwall = wcntr   -  1
      nn34  = n34cntr -  1
      ns12  = s12cntr -  1
      nw23  = w23cntr -  1
      ne14  = e14cntr -  1
      nww2   = w4cntr  -  1
      nww3   = w2cntr  -  1
      nee4   = ee4cntr  -  1
      nee1   = e3cntr  -  1
      return
      end 
C
      
      subroutine lrt2d_ptsort(levelbox,icolbox,irowbox,nboxes,nlev,
     1     iparentbox, ichildbox, nblevel, iboxlev, istartlev,
     2     src, srcsort, isrcsort,
     3     isrcladder, ns, zll, blength, ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-----Global variables
      implicit none
      real *8 xlength
      integer levelbox(*), ns
      integer nlev, nboxes
      integer icolbox(*), irowbox(*)
      integer iparentbox(*), ichildbox(4,*)
      integer nblevel(0:1), iboxlev(*), istartlev(0:1)
      integer isrcsort(*), isrcladder(2,*)
      integer ier
      real *8 src(2,*), srcsort(2,*), zll(2), blength
c-----Local variables
      integer i, ibox, iflag
      integer j, istart, iend, l, jj
      integer levflag, nnodes
      integer ier2
      real *8 xstart, xincrem
      real *8, allocatable :: srctemp(:,:)
      integer, allocatable :: itemp(:)
      integer, allocatable :: itemp2(:)

C
      allocate(srctemp(2,ns))
      allocate(itemp(ns))
      allocate(itemp2(ns))

      ier = 0

      do i = 1, nboxes
         isrcladder(1,i) = 0
         isrcladder(2,i) = -1
      end do

      do i = 1,ns
         isrcsort(i) = i
         srcsort(1,i) = src(1,i)
         srcsort(2,i) = src(2,i)
      enddo

c     initialize source ladder

      istart = istartlev(0)
      ibox = iboxlev(istart)
      isrcladder(1,ibox) = 1
      isrcladder(2,ibox) = ns

      xlength = blength
      do i = 0, nlev
         istart = istartlev(i)
         iend = istart + nblevel(i) - 1
         do jj = istart, iend
            ibox = iboxlev(jj)
            nnodes = isrcladder(2,ibox)-isrcladder(1,ibox)+1
            if (ichildbox(1,ibox) .gt. 0 .and. nnodes .gt. 0) then

c     assign points of this box to the appropriate child box
               ier2 = 0
               call lrt2d_pts2child(ibox,isrcladder,isrcsort,
     1              itemp,itemp2,iparentbox,ichildbox,irowbox,icolbox,
     2              srcsort,srctemp,ier2,zll,xlength)
               if (ier2 .ne. 0) ier = 4
               isrcladder(1,ibox) = 0
               isrcladder(2,ibox) = -1
            endif
         enddo
         xlength = xlength/2.0d0
      enddo

      return
      end

      subroutine lrt2d_ptsort_wc(levelbox,icolbox,irowbox,nboxes,nlev,
     1     iparentbox, ichildbox, nblevel, iboxlev, istartlev,
     2     src, srcsort, isrcsort,
     3     isrcladder, ns, zll, blength, ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     In this version, the ladder indicates which points are in the
c     box or its children, children's children, etc.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-----Global variables
      implicit none
      real *8 xlength
      integer levelbox(*), ns
      integer nlev, nboxes
      integer icolbox(*), irowbox(*)
      integer iparentbox(*), ichildbox(4,*)
      integer nblevel(0:1), iboxlev(*), istartlev(0:1)
      integer isrcsort(*), isrcladder(2,*)
      integer ier
      real *8 src(2,*), srcsort(2,*), zll(2), blength
c-----Local variables
      integer i, ibox, iflag
      integer j, istart, iend, l, jj
      integer levflag, nnodes
      integer ier2
      real *8 xstart, xincrem
      real *8, allocatable :: srctemp(:,:)
      integer, allocatable :: itemp(:)
      integer, allocatable :: itemp2(:)

C
      allocate(srctemp(2,ns))
      allocate(itemp(ns))
      allocate(itemp2(ns))

      ier = 0

      do i = 1, nboxes
         isrcladder(1,i) = 0
         isrcladder(2,i) = -1
      end do

      do i = 1,ns
         isrcsort(i) = i
         srcsort(1,i) = src(1,i)
         srcsort(2,i) = src(2,i)
      enddo

c     initialize source ladder

      istart = istartlev(0)
      ibox = iboxlev(istart)
      isrcladder(1,ibox) = 1
      isrcladder(2,ibox) = ns

      xlength = blength
      do i = 0, nlev
         istart = istartlev(i)
         iend = istart + nblevel(i) - 1
         do jj = istart, iend
            ibox = iboxlev(jj)
            nnodes = isrcladder(2,ibox)-isrcladder(1,ibox)+1
            if (ichildbox(1,ibox) .gt. 0 .and. nnodes .gt. 0) then

c     assign points of this box to the appropriate child box
               ier2 = 0
               call lrt2d_pts2child(ibox,isrcladder,isrcsort,
     1              itemp,itemp2,iparentbox,ichildbox,irowbox,icolbox,
     2              srcsort,srctemp,ier2,zll,xlength)
               if (ier2 .ne. 0) ier = 4

            endif
         enddo
         xlength = xlength/2.0d0
      enddo

      return
      end

      subroutine lrt2d_discsort(levelbox,icolbox,irowbox,nboxes,nlev,
     1     iparentbox, ichildbox, nblevel, iboxlev, istartlev,
     2     pts, ptsrad, npts, iboxdiscs, iboxnumdiscs, maxdiscs, 
     3     zll, blength, ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT
c
c     levelbox(nboxes)      : level of each box
c     icolbox(nboxes)       : column of each box
c     irowbox(nboxes)       : row of each box
c     nboxes                : number of boxes
c     levelbox(nboxes)      : parent of each box
c     ichildbox(4,nboxes)   : children of each box
c     nblevel(0:nlev)       : number of boxes on each level 
c     iboxlev(nboxes)       : boxes sorted by level
c     istartlev(0:nlev)     : boxes for level i are stored in 
c                    iboxlev(istartlev(i):istartlev(i)+nblevel(i)-1)
c     pts(2,npts)           : coordinates of disc centers
c     ptsrad(npts)          : radii of discs
c     npts                  : number of discs
c     iboxnumdiscs(nboxes)  : integer storage, length nboxes
c     maxdiscs              : if maxdiscs > 0, then this is the first
c                             dimension of iboxdiscs array
c     zll(2)                : lower-left corner of root box
c     blength               : side length of root box
c
c     OUTPUT
c
c     if maxdiscs < 0 ---- QUERY MODE
c
c     ier                   : error flag, ier = 0 means normal
c                             ier = 8 or 16, problem with tree
c     maxdiscs              : maximum number of discs intersecting a 
c                             leaf
c
c     if maxdiscs > 0
c
c     ier                        : error flag, ier = 0 means normal
c                                  ier = 8 or 16, problem with tree
c                                  ier = 32, maxdiscs is too small
c     iboxdiscs(maxdiscs,nboxes) : discs which intersect each leaf
c                                  if ier = 32, then the first 
c                                  maxdiscs discs are stored
c     iboxnumdiscs(nboxes)       : number of discs intersecting box
c                                  if ier = 32, these numbers are still
c                                  accurate. maxdiscs should be set 
c                                  to max(iboxnumdiscs(i))
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-----Global variables
      implicit none
      real *8 xlength
      integer levelbox(*), npts
      integer nlev, nboxes
      integer icolbox(*), irowbox(*)
      integer iparentbox(*), ichildbox(4,*)
      integer nblevel(0:1), iboxlev(*), istartlev(0:1)
      integer ier, imode, maxdiscs
      real *8 zll(2), blength, pts(2,*), ptsrad(*)
      integer iboxdiscs(1:maxdiscs,*), iboxnumdiscs(*)
c-----Local variables
      integer i, ibox, iflag, iboxnext, maxdiscstemp
      integer j, istart, iend, l, jj, jpt, iii
      integer levflag, nnodes, ndiscsbox
      integer irow, icol, iparentlev, icurrentlev
      integer ier2
      real *8 xp, yp, xc, yc, radp, xl2
      integer, allocatable :: idiscslev(:,:)
      integer, allocatable :: ndiscslev(:)
      integer, allocatable :: ichildlev(:)

      ier = 0

      if (maxdiscs .ge. 0) then
         do i = 1,nboxes
            iboxnumdiscs(i) = -1
         enddo
      endif

c     allocate temporary storage
      allocate(idiscslev(npts,0:nlev))
      allocate(ndiscslev(0:nlev))
      allocate(ichildlev(0:nlev))
      maxdiscstemp = 0
c     all pts (should) intersect root box
      do i = 1,npts
         idiscslev(i,0) = i
      enddo

      ibox = iboxlev(istartlev(0))
      iboxnext = ichildbox(1,ibox)
      if (iboxnext .lt. 0) then
c     root box is leaf
         maxdiscstemp = npts
         if (maxdiscs .ge. 0) then
            iboxnumdiscs(ibox) = npts
            do iii = 1,min(maxdiscs,npts)
               iboxdiscs(iii,ibox) = iii
            enddo
            if (maxdiscstemp .gt. maxdiscs) ier = 32
         else
            maxdiscs = maxdiscstemp
         endif

         return
      else
         ibox = iboxnext
         xlength = blength/2.0d0
         iparentlev = 0
         icurrentlev = 1
         ndiscslev(0) = npts
         ichildlev(1) = 1
      endif

      do i = 1,nboxes*nlev
         
c     find which discs intersect this box
         irow = irowbox(ibox)
         icol = icolbox(ibox)
         call lrt2d_getcenter(xc,yc,irow,icol,zll,xlength)
         
         xl2 = xlength/2.0d0

         ndiscsbox = 0

         do j = 1,ndiscslev(iparentlev)
            jpt = idiscslev(j,iparentlev)
            xp = pts(1,jpt)
            yp = pts(2,jpt)
            radp = ptsrad(jpt)
            call lrt2d_rectintdisc(iflag,xp,yp,radp,
     1           xc,yc,xl2,xl2)
            idiscslev(ndiscsbox+1,icurrentlev) = jpt
            ndiscsbox = ndiscsbox+iflag
         enddo

         ndiscslev(icurrentlev) = ndiscsbox

         if (ichildbox(1,ibox) .lt. 0) then
c     current box is a leaf
            if (ndiscslev(icurrentlev) .gt. maxdiscstemp) then
               maxdiscstemp = ndiscslev(icurrentlev)
            endif

            if (maxdiscs .ge. 0) then
               iboxnumdiscs(ibox) = ndiscslev(icurrentlev)
               do iii = 1,min(maxdiscs,ndiscslev(icurrentlev))
                  iboxdiscs(iii,ibox) = idiscslev(iii,icurrentlev)
               enddo
            endif
            
c     find next box
            if (ichildlev(icurrentlev) .lt. 4) then
c     next box is from same parent
               ichildlev(icurrentlev) = ichildlev(icurrentlev)+1
               iboxnext = ichildbox(ichildlev(icurrentlev),
     1              iparentbox(ibox))
               if (iboxnext .lt. 0) then
                  ier = 8
                  return
               else
c     found next box
                  ibox = iboxnext
                  goto 5000
               endif
            else
c     currentbox is 4th child, go up until next box is found
               do iii = 1,nlev
                  ibox = iparentbox(ibox)
                  icurrentlev = iparentlev
                  iparentlev = iparentlev-1
                  xlength = xlength*2.0d0
                  if (icurrentlev .eq. 0) goto 6000
                  
                  if (ichildlev(icurrentlev) .lt. 4) then
                     ichildlev(icurrentlev) = 
     1                    ichildlev(icurrentlev)+1
                     iboxnext = ichildbox(ichildlev(icurrentlev),
     1                    iparentbox(ibox))
                     if (iboxnext .lt. 0) then
                        ier = 8
                        return
                     else
c     found next box
                        ibox = iboxnext
                        goto 5000
                     endif
                  endif
               enddo
c     failed to find next box
               ier = 16
               return

            endif
         else
c     box has children, next box is first child
            iboxnext = ichildbox(1,ibox)
            icurrentlev = icurrentlev+1
            iparentlev = iparentlev+1
            ichildlev(icurrentlev) = 1
            xlength = xlength/2.0d0
            
            if(iboxnext .lt. 0) then
               ier = 8
               return
            else
c     found next box
               ibox = iboxnext
               goto 5000
            endif
            
         endif

 5000    continue
      enddo

 6000 continue


      if (maxdiscs .lt. 0) then
         maxdiscs = maxdiscstemp
      else if (maxdiscstemp .gt. maxdiscs) then
         ier = 32         
      endif
      

      return
      end

      
      subroutine lrt2d_findboxes(iboxes, pts, npts, icolbox, irowbox,
     1     ichildbox, nlev, nblevel, iboxlev, istartlev, 
     2     levelbox, zll, blength, ier)
C***********************************************************************
C     This subroutine finds the leaf boxes containing the points in pts
C
C     input:
C
C     pts(2,npts) - the coordinates of the points
C
C     levelbox - is an array determining the level of each box
C
C     nlev - is the finest level
C
c     ichildbox - denotes the four children of each box
C
C     nblevel - is the total number of boxes per level
C
C     iboxlev - is the array in which the boxes are arranged
C
C     istartlev - is the pointer to where each level begins
C               in the iboxlev array
c
C     OUTPUT:
C  
C     iboxes - the leaf box containing pts(1:2,i) is iboxes(i)
C
C     ier - if ier = 0, the algorithm executed correctly.
C
C***********************************************************************
      implicit none
c     global
      real *8 pts(2,*), zll(2), blength
      integer iboxes(*), npts
      integer  nlev, ier
      integer  ichildbox(4,*)
      integer icolbox(*), irowbox(*), levelbox(*)
      integer  nblevel(0:1), iboxlev(*), istartlev(0:1)
c     local
      integer i

      do i = 1,npts
         call lrt2d_findbox(iboxes(i),pts(1,i),pts(2,i),
     1        icolbox, irowbox, ichildbox, nlev, nblevel, iboxlev, 
     2        istartlev, levelbox, zll, blength, ier)
      enddo

      return
      end

      subroutine lrt2d_findboxes_2(iboxes, xs, ys, npts, icolbox, 
     1     irowbox, ichildbox, nlev, nblevel, iboxlev, istartlev, 
     2     levelbox, zll, blength, ier)
C***********************************************************************
C     This subroutine finds the leaf boxes containing the points in pts
C
C     input:
C
C     xs(npts), ys(npts) - the coordinates of the points
C
C     levelbox - is an array determining the level of each box
C
C     nlev - is the finest level
C
c     ichildbox - denotes the four children of each box
C
C     nblevel - is the total number of boxes per level
C
C     iboxlev - is the array in which the boxes are arranged
C
C     istartlev - is the pointer to where each level begins
C               in the iboxlev array
c
C     OUTPUT:
C  
C     iboxes - the leaf box containing pts(1:2,i) is iboxes(i)
C
C     ier - if ier = 0, the algorithm executed correctly.
C
C***********************************************************************
      implicit none
c     global
      real *8 xs(*), ys(*), zll(2), blength
      integer iboxes(*), npts
      integer  nlev, ier
      integer  ichildbox(4,*)
      integer icolbox(*), irowbox(*), levelbox(*)
      integer  nblevel(0:1), iboxlev(*), istartlev(0:1)
c     local
      integer i

      do i = 1,npts
         call lrt2d_findbox(iboxes(i),xs(i),ys(i),
     1        icolbox, irowbox, ichildbox, nlev, nblevel, iboxlev, 
     2        istartlev, levelbox, zll, blength, ier)
      enddo

      return
      end
      
      subroutine lrt2d_findbox(ibox, x, y, icolbox, irowbox,
     1     ichildbox, nlev, nblevel, iboxlev, istartlev, 
     2     levelbox, zll, blength, ier)
C***********************************************************************
C     This subroutine finds the leaf box containing the point x, y
C
C     input:
C
C     x, y - the location whose box is to be found
C
C     levelbox - is an array determining the level of each box
C
C     nlev - is the finest level
C
c     ichildbox - denotes the four children of each box
C
C     nblevel - is the total number of boxes per level
C
C     iboxlev - is the array in which the boxes are arranged
C
C     istartlev - is the pointer to where each level begins
C               in the iboxlev array
c
C     OUTPUT:
C  
C     ibox - the leaf box containing (x,y)
C
C     ier - if ier = 0, the algorithm executed correctly.
C
C***********************************************************************
      implicit none
c-----global variables
      real *8 x, y, zll(2), blength
      integer ibox
      integer  nlev, ier
      integer  ichildbox(4,*)
      integer icolbox(*), irowbox(*), levelbox(*)
      integer  nblevel(0:1), iboxlev(*), istartlev(0:1)
c-----local variables
      real *8 xc, yc, quarter_side, xlength
      integer ichild, i, j
      logical xleft, ydown
      integer ibox_col, ibox_row, ibox_level
      real *8 h, a, b, c, d

      ier = 0

c     the algorithm works by recursively selecting
c     the childbox that (x,y) belongs to until it
c     reaches a childless box.

c     quarter_side is set recursively below.

      quarter_side = .25d0*blength
      xc = zll(1)+0.5d0*blength
      yc = zll(2)+0.5d0*blength
      i = istartlev(0)
      ibox = iboxlev(i)
      ibox_level = 0

      if (ichildbox(1,ibox) .lt. 0) go to 1000

c     find box by recursively checking appropriate child

      do i = 1,nlev
         ibox_level = ibox_level + 1
         xleft = (x-xc .lt. 0.0d0)
         ydown = (y-yc .lt. 0.0d0)
         if (xleft .and. ydown) then
            ichild = 4
            xc = xc - quarter_side
            yc = yc - quarter_side
         else if (xleft .and. (.not. ydown)) then
            ichild = 1
            xc = xc - quarter_side
            yc = yc + quarter_side
         else if ( (.not. xleft) .and. ydown) then
            ichild = 3
            xc = xc + quarter_side
            yc = yc - quarter_side
         else
            ichild = 2
            xc = xc + quarter_side
            yc = yc + quarter_side
         endif
         quarter_side = quarter_side/2.0d0
         ibox = ichildbox(ichild,ibox)
         if (ichildbox(1,ibox) .lt. 0) goto 1000
      enddo

 1000 continue

c     check we got the right box

      ibox_col = icolbox(ibox)
      ibox_row = irowbox(ibox)
      h = blength/(2.0d0**ibox_level)
      a = zll(1)+h*(ibox_col-1)
      b = zll(1)+h*ibox_col
      c = zll(2)+h*(ibox_row-1)
      d = zll(2)+h*ibox_row
      if (x .lt. a .or. x .gt. b .or.
     1    y .lt. c .or. y .gt. d) then
         ier = 4

      endif
      
      return
      end

      subroutine lrt2d_getcorners(xc1,yc1,xc2,yc2,xc3,yc3,xc4,yc4,
     1     irow, icol, zll, xlength)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine returns the x and y points for each corner 
c     of the box in the given row and column, with side length
c     xlength.
c
c     The corners of a box are numbered clockwise from the top left
c
c     1 2
c     4 3
c
c     the same as the convention used to number child boxes.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 xc1,yc1,xc2,yc2,xc3,yc3,xc4,yc4
      integer irow, icol
      real *8 xlength, zll(2)

      xc1 = zll(1) + (icol-1)*xlength
      yc1 = zll(2) + irow*xlength
      xc2 = zll(1) + icol*xlength
      yc2 = zll(2) + irow*xlength
      xc3 = zll(1) + icol*xlength
      yc3 = zll(2) + (irow-1)*xlength
      xc4 = zll(1) + (icol-1)*xlength
      yc4 = zll(2) + (irow-1)*xlength

      return
      end


      subroutine lrt2d_getcenter(xc,yc,irow, icol,zll,xlength)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine returns the x and y points for the center 
c     of the box in the given row and column, with side length
c     xlength.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 xc, yc, zll(2)
      integer irow, icol
      real *8 xlength

      xc = zll(1) + (icol-0.5d0)*xlength
      yc = zll(2) + (irow-0.5d0)*xlength

      return
      end



      subroutine lrt2d_isinbox(a,b,c,d,x,y,ifin)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine returns ifin = true if (x,y) is in the interior
c     of the box [a,b] X [c,d]
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 a,b,c,d,x,y
      logical ifin

      ifin = (a .lt. x) .and. (x .lt. b) .and. (c .lt. y) 
     1                  .and. (y .lt. d)

      return
      end

      subroutine lrt2d_mknbrs(neighbors,nnbrs,nboxes,ichildbox,
     1     iparentbox,icolleagbox,icolbox,irowbox)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      integer neighbors(12,*), nnbrs(*), nboxes, ichildbox(4,*) 
      integer iparentbox(*), icolleagbox(9,*), icolbox(*), irowbox(*)
c     local variables
      integer i

      do i = 1,nboxes
         if (ichildbox(1,i) .lt. 0) then
            call lrt2d_getnbrs(neighbors(1,i),nnbrs(i),i,ichildbox,
     1           iparentbox,icolleagbox,icolbox,irowbox)
         else
            nnbrs(i) = 0
         endif
      enddo

      return
      end

      subroutine lrt2d_getnbrs(neighbors,nnbrs,ibox,ichildbox,
     1     iparentbox,icolleagbox,icolbox,irowbox)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine returns the neighbors of the box ibox (assuming
c     a level restricted tree)
c
c     The indeces of the neighbors are stored in the first nnbrs
c     entries of the neighbors array. 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      integer neighbors(*), nnbrs, ibox, ichildbox(4,*), iparentbox(*)
      integer icolleagbox(9,*), icolbox(*), irowbox(*)
c     local variables
      integer iparent, irow, icol, irowparent, icolparent, icoll
      integer irowcoll, icolcoll
      integer irightcol, itoprow
      integer i

      nnbrs = 0
      iparent = iparentbox(ibox)

      if (iparent .lt. 0) then
         return
      endif

      irow = irowbox(ibox)
      icol = icolbox(ibox)
      if ( mod(irow,2) .eq. 1) then
         itoprow = -1
      else
         itoprow = 1
      endif
      if ( mod(icol,2) .eq. 1) then
         irightcol = -1
      else
         irightcol = 1
      endif
      
      irowparent = irowbox(iparent)
      icolparent = icolbox(iparent)

c     go through colleagues of parent, for a childless colleague
c     see if it contacts the box
      do i = 1,9
         if (i .ne. 5) then
            icoll = icolleagbox(i,iparent)
c     check that colleague exists and is childless
c            write(*,*) 'i ', i
c            write(*,*) 'icoll ', icoll, 'iparent ', iparent
            if (icoll .gt. 0 .and. ichildbox(1,icoll) .lt. 0) then
               irowcoll = irowbox(icoll)
               icolcoll = icolbox(icoll)
c     check that the colleague touches the box
               if ( (irowcoll-irowparent)*itoprow .ge. 0 .and.
     1              (icolcoll-icolparent)*irightcol .ge. 0 ) then
                  nnbrs = nnbrs + 1
                  neighbors(nnbrs) = icoll
               endif
            endif
         endif
      enddo

c     go through colleagues of ibox. add childless colleagues to
c     neighbors. for colleagues with children, add children which 
c     contact ibox
      do i = 1,9
         if (i .ne. 5) then
            icoll = icolleagbox(i,ibox)
c     check for a childless colleague
            if (icoll .gt. 0 .and. ichildbox(1,icoll) .lt. 0) then
               nnbrs = nnbrs+1
               neighbors(nnbrs) = icoll
c     for colleagues with children, add the children which touch
c     ibox
            else if (icoll .gt. 0) then
               if (i .eq. 1) then
                  nnbrs = nnbrs + 1
                  neighbors(nnbrs) = ichildbox(2,icoll)
               else if (i.eq.2) then
                  nnbrs = nnbrs + 1
                  neighbors(nnbrs) = ichildbox(1,icoll)
                  nnbrs = nnbrs + 1
                  neighbors(nnbrs) = ichildbox(2,icoll)
               else if (i.eq.3) then
                  nnbrs = nnbrs + 1
                  neighbors(nnbrs) = ichildbox(1,icoll)
               else if (i.eq.4) then
                  nnbrs = nnbrs + 1
                  neighbors(nnbrs) = ichildbox(2,icoll)
                  nnbrs = nnbrs + 1
                  neighbors(nnbrs) = ichildbox(3,icoll)
               else if (i.eq.6) then
                  nnbrs = nnbrs + 1
                  neighbors(nnbrs) = ichildbox(1,icoll)
                  nnbrs = nnbrs + 1
                  neighbors(nnbrs) = ichildbox(4,icoll)
               else if (i.eq.7) then
                  nnbrs = nnbrs + 1
                  neighbors(nnbrs) = ichildbox(3,icoll)
               else if (i.eq.8) then
                  nnbrs = nnbrs + 1
                  neighbors(nnbrs) = ichildbox(3,icoll)
                  nnbrs = nnbrs + 1
                  neighbors(nnbrs) = ichildbox(4,icoll)
               else if (i.eq.9) then
                  nnbrs = nnbrs + 1
                  neighbors(nnbrs) = ichildbox(4,icoll)
               endif
            endif
         endif
      enddo
               
      return
      end
      



C
c*************************************************************************
C      subroutine fixtree
c*************************************************************************

C     The following subroutine is designed to take a correctly defined
C     tree and alter it so that no two boxes that contact each other
C     are more than one level apart.  This is corrected only by adding
C     boxes.  The procedure involves flagging down bigger boxes and
C     dividing them and their children as necessary.
C     This routine also produces an array of colleagues for each box
C     that is numbered in the correct order.  All of the children are set
C     so that they satisfy our ordering convention.
C     The algorithm in the periodic case works the same way, it is just 
C     that upon subdivision the new colleagues must be put down to 
C     account for the periodicity.
C
C     INPUT:
C
C     LEVELBOX is an array determining the level of each box
C     IPARENTBOX denotes the parent of each box
C     ICHILDBOX denotes the four children of each box
C     ICOLBOX denotes the column of each box
C     IROWBOX denotes the row of each box
C     ICOLLEAGBOX denotes the colleagues of a given box
C     NBOXES is the total number of boxes
C     NLEV is the finest level
C     NBLEVEL is the total number of boxes per level
C     IBOXLEV is the array in which the boxes are arranged
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C     IFLAG is just a dummy array
C     MAXBOXES is the maximum number of boxes we have storage for
C
C     OUTPUT:
C
C     ICOLBOX, IROWBOX, ICOLLEAGBOX, NBOXES, and all other
C     arrays describing the tree may be change on output
C
C*************************************************************************
      subroutine lrt2d_fix(levelbox,iparentbox,ichildbox,icolbox, 
     1             irowbox,icolleagbox,nboxes,nlev,
     2             nblevel,iboxlev,istartlev,
     3             iflag, maxboxes,itemparray)
      implicit none
C-----Global variables
      integer levelbox(1), icolleagbox(9,1)
      integer maxboxes
      integer iparentbox(1), ichildbox(4,1)
      integer icolbox(1), irowbox(1)
      integer nboxes, nlev
      integer nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer itemparray(1)
c-----Local variables
      integer iflag(maxboxes)
      integer ichild(4),icoll(9), ibox
      integer i, ipar, itest, j, nb
      integer itemp, ntemp, jcntr, icntr
      integer istart, istop
      integer iboxladdrnew(1:2,0:1000)
C
c     Let's sort all of the boxes by level.
C     This takes the place of the LADDER structure
C     in the uniform case.  boxes are sorted into IBOXLEV...
C
      call lrt2d_sortboxes(levelbox,nboxes,nlev,
     1           nblevel,iboxlev,istartlev)



C
C     Generate colleagues for each box in standard order.
C
      call lrt2d_mkcolls(icolbox,irowbox,icolleagbox,nboxes,nlev,
     1       iparentbox,ichildbox,nblevel,iboxlev, istartlev)
C
C     Initialize all of the box flags to zero.
C
      do i = 1, nboxes
         iflag(i) = 0
      end do
C
c     Flag boxes that violate level restriction by "1".
C     Violation refers to any box that is directly touching 
C     a box that is more than one level finer.
C     Method:  
C     1) Carry out upward pass. For each box B,
C     look at colleagues of B's grandparent. 
C     2) See if any of those colleagues are childless
C     and in contact with B.
C     Note that we only need to get up to level two, as
C     we will not find a violation at a coarser level
C     than that.
      do i = nlev, 2, -1
         do j = istartlev(i), istartlev(i) + nblevel(i) - 1
            ibox = iboxlev(j)
            ipar  = iparentbox(ibox)
            itest = iparentbox(ipar)
            icoll(1) = icolleagbox(1,itest)
            icoll(2) = icolleagbox(2,itest)
            icoll(3) = icolleagbox(3,itest)
            icoll(4) = icolleagbox(4,itest)
            icoll(5) = icolleagbox(5,itest)
            icoll(6) = icolleagbox(6,itest)
            icoll(7) = icolleagbox(7,itest)
            icoll(8) = icolleagbox(8,itest)
            icoll(9) = icolleagbox(9,itest)
            ichild(1) = ichildbox(1,itest)
            ichild(2) = ichildbox(2,itest)
            ichild(3) = ichildbox(3,itest)
            ichild(4) = ichildbox(4,itest)
            do nb = 1, 9
               itemp = icoll(nb)
               if(itemp .gt. 0 .and. ichildbox(1,itemp) .lt. 0)then
c             The neighboring box is not divided
C             we could have problems.
                 if (nb .eq. 1)then
                   if(ipar .eq. ichild(4))then
                       iflag(itemp) = 1
                   end if
                 elseif (nb .eq. 2)then
                   if(ipar .eq. ichild(3) .or. ipar .eq. ichild(4))then
                       iflag(itemp) = 1
                   end if
                 elseif (nb .eq. 3)then
                   if(ipar .eq. ichild(3))then
                       iflag(itemp) = 1
                   end if
                 elseif (nb .eq. 4)then
                   if(ipar .eq. ichild(4) .or. ipar .eq. ichild(1))then
                       iflag(itemp) = 1
                   end if
                 elseif (nb .eq. 6)then
                   if(ipar .eq. ichild(2) .or. ipar .eq. ichild(3))then
                       iflag(itemp) = 1
                   end if
                 elseif (nb .eq. 7)then 
                   if(ipar .eq. ichild(1))then
                       iflag(itemp) = 1
                   end if
                 elseif (nb .eq. 8)then
                   if(ipar .eq. ichild(1) .or. ipar .eq. ichild(2))then
                       iflag(itemp) = 1
                   end if
                 elseif (nb .eq. 9)then
                   if(ipar .eq. ichild(2))then
                       iflag(itemp) = 1
                   end if
                 endif
               endif
            enddo
         enddo
      enddo
c     Find all of the boxes that need to be
C     given a flag+.  A flag+ box will be denoted by 
C     setting IFLAG(box) = 2.
C     This refers to any box that is not already flagged
C     and is bigger than and is contacting a flagged box
C     or another box that has already been given a flag+.
C     It is found by performing an upward pass
C     and looking at a flagged box's parents colleagues
C     and a flag+ box's parents colleagues and seeing if
C     they are childless and present the case where a 
C     bigger box is contacting a flagged or a flag+ box.
      do i = nlev, 2, -1
         do j = istartlev(i), istartlev(i) + nblevel(i) - 1
            ibox = iboxlev(j)
            if(iflag(ibox) .eq. 1 .or. iflag(ibox) .eq. 2)then
               ipar  = iparentbox(ibox)
               icoll(1) = icolleagbox(1,ipar)
               icoll(2) = icolleagbox(2,ipar)
               icoll(3) = icolleagbox(3,ipar)
               icoll(4) = icolleagbox(4,ipar)
               icoll(5) = icolleagbox(5,ipar)
               icoll(6) = icolleagbox(6,ipar)
               icoll(7) = icolleagbox(7,ipar)
               icoll(8) = icolleagbox(8,ipar)
               icoll(9) = icolleagbox(9,ipar)
               ichild(1) = ichildbox(1,ipar)
               ichild(2) = ichildbox(2,ipar)
               ichild(3) = ichildbox(3,ipar)
               ichild(4) = ichildbox(4,ipar)
               do nb = 1, 9
                  itemp = icoll(nb)
c           Let's check using the same criteria as above, but noting that
C           a flag will take precedence over a flag+.
                  if(itemp .gt. 0 .and. ichildbox(1,itemp) .lt. 0 
     1                .and. iflag(itemp) .ne. 1)then
c             The neighboring box is not divided
C             we could have problems.
                    if (nb .eq. 1)then
                      if(ibox .eq. ichild(4))then
                          iflag(itemp) = 2
                      end if
                    elseif (nb .eq. 2)then
                      if(ibox.eq.ichild(3) .or. ibox.eq.ichild(4))then
                          iflag(itemp) = 2
                      end if
                    elseif (nb .eq. 3)then
                      if(ibox .eq. ichild(3))then
                          iflag(itemp) = 2
                      end if
                    elseif (nb .eq. 4)then
                      if(ibox.eq.ichild(4) .or. ibox.eq.ichild(1))then
                          iflag(itemp) = 2
                      end if
                    elseif (nb .eq. 6)then
                      if(ibox.eq.ichild(2) .or. ibox.eq.ichild(3))then
                          iflag(itemp) = 2
                      end if
                    elseif (nb .eq. 7)then
                      if(ibox .eq. ichild(1))then
                          iflag(itemp) = 2
                      end if
                    elseif (nb .eq. 8)then
                      if(ibox.eq.ichild(1) .or. ibox.eq.ichild(2))then
                          iflag(itemp) = 2
                      end if
                    elseif (nb .eq. 9)then
                      if(ibox .eq. ichild(2))then
                          iflag(itemp) = 2
                      end if
                    endif
                  endif
               enddo
            endif
         enddo
      enddo
c     Now let's divide the boxes that need to be immediately
C     divided up.  All of the flagged and flag+ boxes need to
C     be divided one time.  The distinction lies in the fact
C     that the children of a flag+ box will never need to be
C     divided but the children of a flagged box may need to 
C     be divided further.
C     Below, all flagged and flag+ boxes are divided once.  The
C     children of a flag+ box are left unflagged while those of
C     the flagged boxes are given a flag++ (denoted by setting
C     IFLAG(box) = 3) which will be needed in the downward pass.     
      ntemp = nboxes
      do i = 1, ntemp
c      dIVIDE FLAGGED BOXES:
         if (iflag(i) .eq. 1)then
            if(ichildbox(1,i) .lt. 0)then
               call lrt2d_subdiv(i,iparentbox,ichildbox,icolleagbox,
     1              nboxes,irowbox,icolbox,levelbox,nlev,
     2              istartlev, nblevel, iboxlev,
     3              itemparray)
            endif
c        Give flag++ to children of flagged boxes.
            itemp = ichildbox(1,i)
            iflag(itemp) = 3
            itemp = ichildbox(2,i)
            iflag(itemp) = 3
            itemp = ichildbox(3,i)
            iflag(itemp) = 3
            itemp = ichildbox(4,i)
            iflag(itemp) = 3
c      Divide flag+ boxes.
         elseif (iflag(i) .eq. 2)then
            if(ichildbox(1,i) .lt. 0)then
               call lrt2d_subdiv(i,iparentbox,ichildbox,icolleagbox,
     1              nboxes,irowbox,icolbox,levelbox,nlev,
     2              istartlev, nblevel, iboxlev,
     3              itemparray)
               itemp = ichildbox(1,i)
               iflag(itemp) = 0
               itemp = ichildbox(2,i)
               iflag(itemp) = 0
               itemp = ichildbox(3,i)
               iflag(itemp) = 0
               itemp = ichildbox(4,i)
               iflag(itemp) = 0

            endif
         endif
      enddo 

c     fix ladder structure
      call lrt2d_sortboxes(levelbox,nboxes,nlev,
     1           nblevel,iboxlev,istartlev)

      iboxladdrnew(1,0) = 0
      iboxladdrnew(2,0) = -1

c     Now we need to do a downward pass.
C     We will concern ourselves only with the children of
C     flagged boxes and their children.  At each level,
C     for each flag++ box, test colleagues children and see
C     if they have children that are contacting you.  If so,
C     divide and flag++ all children that are created.     
      do i = 0, nlev
         ntemp = nboxes
         istart = istartlev(i)
         istop  = istartlev(i) + nblevel(i) - 1
         iboxladdrnew(1,i+1) = nboxes+1
         do 500 j = istart, istop
            ibox = iboxlev(j)
c      Only be concerned with boxes on this level and
C      boxes that are given a flag++:
            if(iflag(ibox) .ne. 3)goto 500
            icoll(1) = icolleagbox(1,ibox)
            icoll(2) = icolleagbox(2,ibox)
            icoll(3) = icolleagbox(3,ibox)
            icoll(4) = icolleagbox(4,ibox)
            icoll(5) = icolleagbox(5,ibox)
            icoll(6) = icolleagbox(6,ibox)
            icoll(7) = icolleagbox(7,ibox)
            icoll(8) = icolleagbox(8,ibox)
            icoll(9) = icolleagbox(9,ibox)
c       Scan colleagues.
            do 400 jcntr = 1, 9
               if(icoll(jcntr) .lt. 0)goto 400
               if(ichildbox(1,icoll(jcntr)) .lt. 0)goto 400
               ichild(1) = ichildbox(1,icoll(jcntr))
               ichild(2) = ichildbox(2,icoll(jcntr))
               ichild(3) = ichildbox(3,icoll(jcntr))
               ichild(4) = ichildbox(4,icoll(jcntr))
c          Scan colleague's children.
               do 300 icntr = 1, 4
                  if (ichildbox(1,ichild(icntr)) .lt. 0)goto 300 
                  if(jcntr .eq. 1 .and. icntr .eq. 2)then
c                    Call subdiv
                     if(ichildbox(1,ibox) .lt. 0)then
                        call lrt2d_subdiv(ibox,iparentbox,ichildbox,
     1                  icolleagbox,nboxes,irowbox,icolbox,levelbox,
     2                  nlev,istartlev, nblevel, iboxlev,itemparray)
                     endif
c                    flag++ all children created
                     itemp = ichildbox(1,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(2,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(3,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(4,ibox)
                     iflag(itemp) = 3
                  elseif(jcntr .eq. 2 .and. 
     1             (icntr .eq. 1 .or. icntr .eq. 2))then
c                    Call subdiv
                     if(ichildbox(1,ibox) .lt. 0)then
                        call lrt2d_subdiv(ibox,iparentbox,ichildbox,
     1                  icolleagbox,nboxes,irowbox,icolbox,levelbox,
     2                  nlev,istartlev, nblevel, iboxlev,itemparray)
                     endif
c                    flag++ all children created
                     itemp = ichildbox(1,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(2,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(3,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(4,ibox)
                     iflag(itemp) = 3
                  elseif(jcntr .eq. 3 .and. icntr .eq. 1)then
c                    Call subdiv
                     if(ichildbox(1,ibox) .lt. 0)then
                        call lrt2d_subdiv(ibox,iparentbox,ichildbox,
     1                  icolleagbox,nboxes,irowbox,icolbox,levelbox,
     2                  nlev,istartlev, nblevel, iboxlev,itemparray)
                     endif
c                    flag++ all children created
                     itemp = ichildbox(1,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(2,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(3,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(4,ibox)
                     iflag(itemp) = 3
                  elseif(jcntr .eq. 4 .and. 
     1             (icntr .eq. 2 .or. icntr .eq. 3))then
c                    Call subdiv
                     if(ichildbox(1,ibox) .lt. 0)then
                        call lrt2d_subdiv(ibox,iparentbox,ichildbox,
     1                  icolleagbox,nboxes,irowbox,icolbox,levelbox,
     2                  nlev,istartlev, nblevel, iboxlev,itemparray)
                     endif
c                    flag++ all children created
                     itemp = ichildbox(1,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(2,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(3,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(4,ibox)
                     iflag(itemp) = 3
                  elseif(jcntr .eq. 6 .and. 
     1             (icntr .eq. 1 .or. icntr .eq. 4))then
c                    Call subdiv
                     if(ichildbox(1,ibox) .lt. 0)then
                        call lrt2d_subdiv(ibox,iparentbox,ichildbox,
     1                  icolleagbox,nboxes,irowbox,icolbox,levelbox,
     2                  nlev,istartlev, nblevel, iboxlev,itemparray)
                     endif
c                    flag++ all children created
                     itemp = ichildbox(1,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(2,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(3,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(4,ibox)
                     iflag(itemp) = 3
                  elseif(jcntr .eq. 7 .and. icntr .eq. 3)then
c                    Call subdiv
                     if(ichildbox(1,ibox) .lt. 0)then
                        call lrt2d_subdiv(ibox,iparentbox,ichildbox,
     1                  icolleagbox,nboxes,irowbox,icolbox,levelbox,
     2                  nlev,istartlev, nblevel, iboxlev,itemparray)
                     endif
c                    flag++ all children created
                     itemp = ichildbox(1,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(2,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(3,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(4,ibox)
                     iflag(itemp) = 3
                  elseif(jcntr .eq. 8 .and. 
     1             (icntr .eq. 3 .or. icntr .eq. 4))then
c                    Call subdiv
                     if(ichildbox(1,ibox) .lt. 0)then
                        call lrt2d_subdiv(ibox,iparentbox,ichildbox,
     1                  icolleagbox,nboxes,irowbox,icolbox,levelbox,
     2                  nlev,istartlev, nblevel, iboxlev,itemparray)
                     endif
c                    flag++ all children created
                     itemp = ichildbox(1,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(2,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(3,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(4,ibox)
                     iflag(itemp) = 3
                  elseif(jcntr .eq. 9 .and. icntr .eq. 4)then
c                    Call subdiv
                     if(ichildbox(1,ibox) .lt. 0)then
                        call lrt2d_subdiv(ibox,iparentbox,ichildbox,
     1                  icolleagbox,nboxes,irowbox,icolbox,levelbox,
     2                  nlev,istartlev, nblevel, iboxlev,itemparray)
                     endif
c                    flag++ all children created
                     itemp = ichildbox(1,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(2,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(3,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(4,ibox)
                     iflag(itemp) = 3
                  endif
300            continue   
400         continue
500      continue
         

         
         istart = iboxladdrnew(1,i)
         istop = iboxladdrnew(2,i)
         do 800 j = istart, istop
            ibox = j
c      Only be concerned with boxes on this level and
C      boxes that are given a flag++:
            if(iflag(ibox) .ne. 3)goto 800
            icoll(1) = icolleagbox(1,ibox)
            icoll(2) = icolleagbox(2,ibox)
            icoll(3) = icolleagbox(3,ibox)
            icoll(4) = icolleagbox(4,ibox)
            icoll(5) = icolleagbox(5,ibox)
            icoll(6) = icolleagbox(6,ibox)
            icoll(7) = icolleagbox(7,ibox)
            icoll(8) = icolleagbox(8,ibox)
            icoll(9) = icolleagbox(9,ibox)
c       Scan colleagues.
            do 700 jcntr = 1, 9
               if(icoll(jcntr) .lt. 0)goto 700
               if(ichildbox(1,icoll(jcntr)) .lt. 0)goto 700
               ichild(1) = ichildbox(1,icoll(jcntr))
               ichild(2) = ichildbox(2,icoll(jcntr))
               ichild(3) = ichildbox(3,icoll(jcntr))
               ichild(4) = ichildbox(4,icoll(jcntr))
c          Scan colleague's children.
               do 600 icntr = 1, 4
                  if (ichildbox(1,ichild(icntr)) .lt. 0)goto 600 
                  if(jcntr .eq. 1 .and. icntr .eq. 2)then
c                    Call subdiv
                     if(ichildbox(1,ibox) .lt. 0)then
                        call lrt2d_subdiv(ibox,iparentbox,ichildbox,
     1                  icolleagbox,nboxes,irowbox,icolbox,levelbox,
     2                  nlev,istartlev, nblevel, iboxlev,itemparray)
                     endif
c                    flag++ all children created
                     itemp = ichildbox(1,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(2,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(3,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(4,ibox)
                     iflag(itemp) = 3
                  elseif(jcntr .eq. 2 .and. 
     1             (icntr .eq. 1 .or. icntr .eq. 2))then
c                    Call subdiv
                     if(ichildbox(1,ibox) .lt. 0)then
                        call lrt2d_subdiv(ibox,iparentbox,ichildbox,
     1                  icolleagbox,nboxes,irowbox,icolbox,levelbox,
     2                  nlev,istartlev, nblevel, iboxlev,itemparray)
                     endif
c                    flag++ all children created
                     itemp = ichildbox(1,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(2,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(3,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(4,ibox)
                     iflag(itemp) = 3
                  elseif(jcntr .eq. 3 .and. icntr .eq. 1)then
c                    Call subdiv
                     if(ichildbox(1,ibox) .lt. 0)then
                        call lrt2d_subdiv(ibox,iparentbox,ichildbox,
     1                  icolleagbox,nboxes,irowbox,icolbox,levelbox,
     2                  nlev,istartlev, nblevel, iboxlev,itemparray)
                     endif
c                    flag++ all children created
                     itemp = ichildbox(1,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(2,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(3,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(4,ibox)
                     iflag(itemp) = 3
                  elseif(jcntr .eq. 4 .and. 
     1             (icntr .eq. 2 .or. icntr .eq. 3))then
c                    Call subdiv
                     if(ichildbox(1,ibox) .lt. 0)then
                        call lrt2d_subdiv(ibox,iparentbox,ichildbox,
     1                  icolleagbox,nboxes,irowbox,icolbox,levelbox,
     2                  nlev,istartlev, nblevel, iboxlev,itemparray)
                     endif
c                    flag++ all children created
                     itemp = ichildbox(1,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(2,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(3,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(4,ibox)
                     iflag(itemp) = 3
                  elseif(jcntr .eq. 6 .and. 
     1             (icntr .eq. 1 .or. icntr .eq. 4))then
c                    Call subdiv
                     if(ichildbox(1,ibox) .lt. 0)then
                        call lrt2d_subdiv(ibox,iparentbox,ichildbox,
     1                  icolleagbox,nboxes,irowbox,icolbox,levelbox,
     2                  nlev,istartlev, nblevel, iboxlev,itemparray)
                     endif
c                    flag++ all children created
                     itemp = ichildbox(1,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(2,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(3,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(4,ibox)
                     iflag(itemp) = 3
                  elseif(jcntr .eq. 7 .and. icntr .eq. 3)then
c                    Call subdiv
                     if(ichildbox(1,ibox) .lt. 0)then
                        call lrt2d_subdiv(ibox,iparentbox,ichildbox,
     1                  icolleagbox,nboxes,irowbox,icolbox,levelbox,
     2                  nlev,istartlev, nblevel, iboxlev,itemparray)
                     endif
c                    flag++ all children created
                     itemp = ichildbox(1,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(2,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(3,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(4,ibox)
                     iflag(itemp) = 3
                  elseif(jcntr .eq. 8 .and. 
     1             (icntr .eq. 3 .or. icntr .eq. 4))then
c                    Call subdiv
                     if(ichildbox(1,ibox) .lt. 0)then
                        call lrt2d_subdiv(ibox,iparentbox,ichildbox,
     1                  icolleagbox,nboxes,irowbox,icolbox,levelbox,
     2                  nlev,istartlev, nblevel, iboxlev,itemparray)
                     endif
c                    flag++ all children created
                     itemp = ichildbox(1,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(2,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(3,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(4,ibox)
                     iflag(itemp) = 3
                  elseif(jcntr .eq. 9 .and. icntr .eq. 4)then
c                    Call subdiv
                     if(ichildbox(1,ibox) .lt. 0)then
                        call lrt2d_subdiv(ibox,iparentbox,ichildbox,
     1                  icolleagbox,nboxes,irowbox,icolbox,levelbox,
     2                  nlev,istartlev, nblevel, iboxlev,itemparray)
                     endif
c                    flag++ all children created
                     itemp = ichildbox(1,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(2,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(3,ibox)
                     iflag(itemp) = 3
                     itemp = ichildbox(4,ibox)
                     iflag(itemp) = 3
                  endif
600            continue   
700         continue
800      continue

c         write(*,*) 'i+1', i+1
         iboxladdrnew(2,i+1) = nboxes

      enddo



c     fix ladder structure
      call lrt2d_sortboxes(levelbox,nboxes,nlev,
     1     nblevel,iboxlev,istartlev)




      return
      end
C
c**********************************************************************
C     subroutine subdiv
c**********************************************************************
C     The following subroutine is designed to divide up a childless
C     box into four children.  
C     The children are placed in correct order (clockwise starting
C     from the upper left corner) so there is no need to 'shuffle' 
C     the child order later on. In the periodic case, the colleagues
C     must be obtained by looking at the potential colleague numbers
C     and their row and column and seeing if they lie outside of 
C     the domain. If they do it must be readjusted to account for the 
C     periodicity.
C
C     INPUT:
C
C     IPARBOX denotes the box being divided
C     IPARENTBOX denotes the parent of each box
C     ICHILDBOX denotes the four children of each box
C     ICOLLEAGBOX denotes the colleagues of a given box
C     NBOXES is the total number of boxes
C     IROWBOX denotes the row of each box
C     ICOLBOX denotes the column of each box
C     LEVELBOX is an array determining the level of each box
C     NLEV is the finest level
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C     NBLEVEL is the total number of boxes per level
C     IBOXLEV is the array in which the boxes are arranged
C     ITEMPARRAY is just a dummy array
C
C     OUTPUT:
C     
C     NBOXES and ICHILDBOX are altered to
C            reflect the addition of new boxes
C
C**********************************************************************
      subroutine lrt2d_subdiv(iparbox,iparentbox,ichildbox,
     1     icolleagbox,nboxes,irowbox,icolbox,levelbox,nlev,
     2     istartlev,nblevel,iboxlev,itemparray)
      IMPLICIT NONE
C-----Global variables
      integer  iparentbox(1), ichildbox(4,1)
      integer  icolbox(1), irowbox(1)
      integer  levelbox(1), nboxes
      integer  iparbox, nlev
      integer  nblevel(0:1), iboxlev(1), istartlev(0:1)
      integer  icolleagbox(9,1)
c-----Local variables
      integer  level, ibox, itemparray(1)
      integer  icolumn, irow, i, j
      integer  icntr, jcntr, isister
      integer  icoltemp, icoltest
      integer  irowtemp, irowtest
      integer  partemp, colleague
      integer  nside, ilev, itest, l

c     Let's initialize the array ITEMPARRAY to zero:
c      do i = 1, nboxes + 4
c         itemparray(i) = 0
c      end do

c     Level, icolumn, and irow refer to the level, column,
C     and row of the parent box, respectively.
      level   = levelbox(iparbox)
      icolumn = icolbox(iparbox)
      irow    = irowbox(iparbox)
c     Here are the new boxes placed in the
C     correct positions.  They are all childless.
C     Their columns and rows are determined from
C     the parents columns and rows.  The level is
C     obviously one level finer than the parent.
      ibox = nboxes + 1
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 1
      irowbox(ibox) = 2*(irow-1) + 1
c
      ibox = nboxes + 2
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 2
      irowbox(ibox) = 2*(irow-1) + 1
c
      ibox = nboxes + 3
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 1
      irowbox(ibox) = 2*(irow-1) + 2
c
      ibox = nboxes + 4
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 2
      irowbox(ibox) = 2*(irow-1) + 2
c
      ichildbox(1,iparbox) = nboxes + 3
      ichildbox(2,iparbox) = nboxes + 4
      ichildbox(3,iparbox) = nboxes + 2
      ichildbox(4,iparbox) = nboxes + 1

      nboxes = nboxes + 4
c      Now let's go through the process of reforming any 
C      necessary colleagues.  For each of the child boxes 
C      that we just formed, all we need to do is scan through
C      the boxes that are children of the above parent boxes 
C      colleagues and test the column and row numbers.  We can 
C      also take advantage of the fact that for every one of 
C      the newly formed boxes colleagues, that box will list 
C      the newly formed box as one of its colleagues.  
C      The colleague numbers can be found easily if we think 
C      of a 'reflection.'  Colleague 1 and 9 are opposites, 
C      3 and 7 etc.
C      First do the free space case:
       do 200 i = 1, 4
          if(ichildbox(1,iparbox) .lt. 0)goto 200
          ibox = ichildbox(i,iparbox)
          icolleagbox(5,ibox) = ibox
          do j = 1, 4
             icolleagbox(j,ibox) = -1
          enddo
          do j = 6, 9
             icolleagbox(j,ibox) = -1
          enddo
          partemp = iparentbox(ibox)
c        IROWTEMP and ICOLTEMP denote the
C        row and column of the test box.
          irowtemp = irowbox(ibox)
          icoltemp = icolbox(ibox)
          do 100 jcntr = 1, 9
C            COLLEAGUE denotes the colleague of the parent box.
             colleague = icolleagbox(jcntr,partemp)
C            If the colleague doesn't exist
C            or is childless, skip it:
             if (colleague .lt. 0)goto 100
             if (ichildbox(1,colleague) .lt. 0)goto 100
c            Otherwise scan the four children:
             do icntr = 1, 4
                j = ichildbox(icntr,colleague)
c               IROWTEST and ICOLTEST denote the row and column of
C               the box being compared to the test box.
                irowtest = irowbox(j)
                icoltest = icolbox(j)
                if(irowtemp .eq. irowtest+1)then
                   if(icoltemp .eq. icoltest+1)then
                      icolleagbox(1,ibox) = j
                      icolleagbox(9,j) = ibox
                   elseif(icoltemp .eq. icoltest)then
                      icolleagbox(2,ibox) = j
                      icolleagbox(8,j) = ibox
                   elseif(icoltemp .eq. icoltest-1)then
                      icolleagbox(3,ibox) = j
                      icolleagbox(7,j) = ibox
                    endif
                elseif(irowtemp .eq. irowtest)then
                   if(icoltemp .eq. icoltest+1)then
                      icolleagbox(4,ibox) = j
                      icolleagbox(6,j) = ibox
                   elseif(icoltemp .eq. icoltest-1)then
                      icolleagbox(6,ibox) = j
                      icolleagbox(4,j) = ibox
                   endif
                elseif(irowtemp .eq. irowtest-1)then
                   if(icoltemp .eq. icoltest+1)then
                      icolleagbox(7,ibox) = j
                      icolleagbox(3,j) = ibox
                   elseif(icoltemp .eq. icoltest)then
                      icolleagbox(8,ibox) = j
                      icolleagbox(2,j) = ibox
                   elseif(icoltemp .eq. icoltest-1)then
                      icolleagbox(9,ibox) = j
                      icolleagbox(1,j) = ibox
                   endif
                endif
             enddo
100      continue
200   continue
      return
      end


C
c**********************************************************************
C     subroutine subdiv
c**********************************************************************
C     The following subroutine is designed to divide up a childless
C     box into four children.  
C     The children are placed in correct order (clockwise starting
C     from the upper left corner) so there is no need to 'shuffle' 
C     the child order later on. In the periodic case, the colleagues
C     must be obtained by looking at the potential colleague numbers
C     and their row and column and seeing if they lie outside of 
C     the domain. If they do it must be readjusted to account for the 
C     periodicity.
C
C     INPUT:
C
C     IPARBOX denotes the box being divided
C     IPARENTBOX denotes the parent of each box
C     ICHILDBOX denotes the four children of each box
C     ICOLLEAGBOX denotes the colleagues of a given box
C     NBOXES is the total number of boxes
C     IROWBOX denotes the row of each box
C     ICOLBOX denotes the column of each box
C     LEVELBOX is an array determining the level of each box
C     NLEV is the finest level
C     ISTARTLEV is the pointer to where each level begins in the
C               IBOXLEV array
C     NBLEVEL is the total number of boxes per level
C     IBOXLEV is the array in which the boxes are arranged
C     ITEMPARRAY is just a dummy array
C
C     OUTPUT:
C     
C     NBOXES and ICHILDBOX are altered to
C            reflect the addition of new boxes
C
C**********************************************************************
      subroutine lrt2d_subdiv1(iparbox,iparentbox,ichildbox,
     1     nboxes,irowbox,icolbox,levelbox)
      IMPLICIT NONE
C-----Global variables
      integer  iparentbox(1), ichildbox(4,1)
      integer  icolbox(1), irowbox(1)
      integer  levelbox(1), nboxes
      integer  iparbox
c-----Local variables
      integer  level, ibox
      integer  icolumn, irow, i, j
      integer  icntr, jcntr, isister
      integer  icoltemp, icoltest
      integer  irowtemp, irowtest
      integer  partemp, colleague
      integer  nside, ilev, itest, l

c     Let's initialize the array ITEMPARRAY to zero:
c      do i = 1, nboxes + 4
c         itemparray(i) = 0
c      end do

c     Level, icolumn, and irow refer to the level, column,
C     and row of the parent box, respectively.
      level   = levelbox(iparbox)
      icolumn = icolbox(iparbox)
      irow    = irowbox(iparbox)
c     Here are the new boxes placed in the
C     correct positions.  They are all childless.
C     Their columns and rows are determined from
C     the parents columns and rows.  The level is
C     obviously one level finer than the parent.
      ibox = nboxes + 1
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 1
      irowbox(ibox) = 2*(irow-1) + 1
c
      ibox = nboxes + 2
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 2
      irowbox(ibox) = 2*(irow-1) + 1
c
      ibox = nboxes + 3
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 1
      irowbox(ibox) = 2*(irow-1) + 2
c
      ibox = nboxes + 4
      levelbox(ibox) = level + 1
      iparentbox(ibox) = iparbox
      ichildbox(1,ibox) = -1
      ichildbox(2,ibox) = -1
      ichildbox(3,ibox) = -1
      ichildbox(4,ibox) = -1
      icolbox(ibox) = 2*(icolumn-1) + 2
      irowbox(ibox) = 2*(irow-1) + 2
c
      ichildbox(1,iparbox) = nboxes + 3
      ichildbox(2,iparbox) = nboxes + 4
      ichildbox(3,iparbox) = nboxes + 2
      ichildbox(4,iparbox) = nboxes + 1

      nboxes = nboxes + 4

      return
      end



      subroutine lrt2d_rectintdisc(iflag,xp,yp,radp,
     1     xc,yc,w2,h2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Determine if rectangle and disc intersect
c
c     INPUT:
c
c     xp         : x coordinate of disc center
c     yp         : y coordinate of disc center
c     radp       : radius of disc
c     xc         : x coordinate of rectangle center
c     yc         : y coordinate of rectangle center
c     w2         : half the width of rectangle
c     h2         : half the height of rectangle
c
c     OUTPUT:
c
c     iflag      : 1 if intersects, 0 otherwise
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 xp, yp, radp, xc, yc, w2, h2
      integer iflag
c     local
      real *8 xdist, ydist, cdist

      iflag = 0

c     by symmetry, check first quadrant

      xdist = dabs(xp-xc)
      ydist = dabs(yp-yc)
      
c     check if outside [0,xc+w/2+r] x [0,yc+h/2+r]

      if (xdist .gt. w2 + radp) return
      if (ydist .gt. h2 + radp) return

c     check if inside [0,xc+w/2] x [0,yc+h/2+r]

      if (xdist .le. w2) then
         iflag = 1
         return
      endif

c     check if inside [0,xc+w/2+r] x [0,yc+h/2]

      if (ydist .le. h2) then
         iflag = 1
         return
      endif

c     check distance to corner

      cdist = (xdist-w2)**2 + (ydist-h2)**2
      if (cdist .le. radp**2) iflag = 1

      return
      end

      subroutine lrt2d_testtree(levelbox,iparentbox,ichildbox,icolbox, 
     1             irowbox,nboxes,nlev,
     2             nblevel,iboxlev,istartlev)
      implicit none
c-----global variables
      integer  levelbox(1)
      integer  iparentbox(1), ichildbox(4,1)
      integer  icolbox(1), irowbox(1)
      integer  nboxes, nlev
      integer  nblevel(0:1), iboxlev(1), istartlev(0:1)
c-----local variables
      integer  i, j
      integer  itemp
      integer  level, levelchild
      integer  irowpar, icolpar
      integer  irowchild, icolchild
      integer  numboxes
      real *8  dnumboxes

c     first test: sweep through all of the boxes,
c     make sure that the level of the parent is
c     one level less than the level of its four children.
      do i = 1, nboxes
        level = levelbox(i)
        if(ichildbox(1,i) .gt. 0)then
          do j = 1, 4
             levelchild = levelbox(ichildbox(j,i))
             if(levelchild .ne. level+1)then
               write(*,*)'levels not properly defined.'
               write(*,*)'abort.'
               write(*,*) i, j
               stop
             endif
          end do
        endif
      end do


c     second test: sweep through all of the boxes,
c     make sure that the parent and children 
c     listings match up.
      do i = 1, nboxes
        if(ichildbox(1,i) .gt. 0)then
          do j = 1, 4
            itemp = ichildbox(j,i)
            if(iparentbox(itemp) .ne. i)then
              write(*,*)'parents and children do not match up.'
              write(*,*)'abort.'
              stop
            endif
          end do
        endif
      end do


c     third test: sweep through all of the boxes,
c     make sure that the column and row numbers of 
c     the parents and children match up correctly.
      do i = 1, nboxes
        icolpar = 2*(icolbox(i) - 1) + 1
        irowpar = 2*(irowbox(i) - 1) + 1
        if(ichildbox(1,i) .gt. 0)then
c          go through all four children.
           do j = 1, 4
             itemp = ichildbox(j,i)
             if(j .eq. 1)then
               irowchild = irowbox(itemp) - 1
               icolchild = icolbox(itemp)
             elseif(j .eq. 2)then
               irowchild = irowbox(itemp) - 1
               icolchild = icolbox(itemp) - 1
             elseif(j .eq. 3)then
               irowchild = irowbox(itemp)
               icolchild = icolbox(itemp) - 1
             elseif(j .eq. 4)then
               irowchild = irowbox(itemp)
               icolchild = icolbox(itemp)
             endif

             if(icolchild .ne. icolpar .or.
     1          irowchild .ne. irowpar)then
                write(*,*)'the columns and rows'
                write(*,*)'of the parents and'
                write(*,*)'children do not correspond'
                write(*,*)'correctly.'
                write(*,*)'abort.'
                stop
             endif
           end do
            
        endif
      end do

c     now call the routine to sort the boxes.
      call lrt2d_sortboxes(levelbox,nboxes,nlev,
     1           nblevel,iboxlev,istartlev)

c     fourth test: do an downward pass through all of the levels,
c     make sure the total number of boxes at each level doesn't 
c     exceed the maximum possible number.
      dnumboxes = 1.0d0
      do i = 0, nlev
        if(nblevel(i) .gt. dnumboxes)then
           write(*,*)'the total number of boxes at'
           write(*,*)'level ',i,'is greater than the'
           write(*,*)'total number allowed.'
           write(*,*) nblevel(i), numboxes
           write(*,*)'abort.'
           stop
        endif
        dnumboxes = 4.0d0 * dnumboxes
      end do
      return
      end


      subroutine lrt2d_pack(itree,litree,levelbox,icolbox,irowbox,
     1     nboxes,nlev,ichildbox,iparentbox,nblevel,istartlev,iboxlev,
     2     ifcolleag,icolleagbox,ifnbr,neighbors,nnbrs,ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c-----Global variables
      integer itree(*), litree, levelbox(*)
      integer icolbox(*), irowbox(*)
      integer ichildbox(4,*), iparentbox(*), nblevel(0:1)
      integer istartlev(0:1), iboxlev(*), icolleagbox(9,*)
      integer neighbors(12,*), nnbrs(*)
      integer ifcolleag, ifnbr
      integer nboxes, nlev, ier
c     local variables
      integer ilevelbox, iicolbox, iirowbox, iichildbox, iiparentbox
      integer inblevel, iistartlev, iiboxlev, iicolleagbox
      integer ineighbors, innbrs, ltot

c     figure out length and indeces

      iicolleagbox = -1
      ineighbors = -1
      innbrs = -1

      ilevelbox = 21
      iicolbox = ilevelbox+nboxes
      iirowbox = iicolbox+nboxes
      iichildbox = iirowbox+nboxes
      iiparentbox = iichildbox+4*nboxes
      iiboxlev = iiparentbox+nboxes
      inblevel = iiboxlev+nboxes
      iistartlev = inblevel+nlev+1
      if (ifcolleag .eq. 1) then
         if (ifnbr .eq. 1) then
            iicolleagbox = iistartlev+nlev+1
            ineighbors = iicolleagbox+9*nboxes
            innbrs = ineighbors+12*nboxes
            ltot = innbrs+nboxes
         else
            iicolleagbox = iistartlev+nlev+1
            ltot = iicolleagbox+9*nboxes
         endif
      else
         if (ifnbr .eq. 1) then
            ineighbors = iistartlev+nlev+1
            innbrs = ineighbors+12*nboxes
            ltot = innbrs+nboxes
         else
            ltot = iistartlev+nlev+1
         endif
      endif

      if (litree .eq. -1) then
         litree = ltot
         return
      endif

      if (litree .lt. ltot) then
         ier = 4
         return
      endif

c     store info

      itree(1) = nboxes
      itree(2) = nlev
      itree(3) = ifcolleag
      itree(4) = ifnbr
      itree(5) = ilevelbox
      itree(6) = iicolbox
      itree(7) = iirowbox
      itree(8) = iichildbox
      itree(9) = iiparentbox
      itree(10) = iiboxlev
      itree(11) = inblevel
      itree(12) = iistartlev
      itree(13) = iicolleagbox
      itree(14) = ineighbors
      itree(15) = innbrs

c     copy over arrays

      call lrt2d_icopy(itree(ilevelbox),levelbox,nboxes)
      call lrt2d_icopy(itree(iicolbox),icolbox,nboxes)
      call lrt2d_icopy(itree(iirowbox),irowbox,nboxes)
      call lrt2d_icopy(itree(iichildbox),ichildbox,4*nboxes)
      call lrt2d_icopy(itree(iiparentbox),iparentbox,nboxes)
      call lrt2d_icopy(itree(iiboxlev),iboxlev,nboxes)
      call lrt2d_icopy(itree(inblevel),nblevel,nlev+1)
      call lrt2d_icopy(itree(iistartlev),istartlev,nlev+1)
      if (ifcolleag .eq. 1) then
         call lrt2d_icopy(itree(iicolleagbox),icolleagbox,9*nboxes)
      endif
      if (ifnbr .eq. 1) then
         call lrt2d_icopy(itree(ineighbors),neighbors,12*nboxes)
         call lrt2d_icopy(itree(innbrs),nnbrs,nboxes)
      endif

      return
      end

      subroutine lrt2d_unpack(itree,levelbox,icolbox,irowbox,
     1     nboxes,nlev,ichildbox,iparentbox,nblevel,istartlev,iboxlev,
     2     ifcolleag,icolleagbox,ifnbr,neighbors,nnbrs,ier,ifquery)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c-----Global variables
      integer itree(*), levelbox(*)
      integer icolbox(*), irowbox(*)
      integer ichildbox(4,*), iparentbox(*), nblevel(0:1)
      integer istartlev(0:1), iboxlev(*), icolleagbox(9,*)
      integer neighbors(12,*), nnbrs(*)
      integer ifcolleag, ifnbr, ier, ifquery
      integer nboxes, nlev
c     local variables
      integer ilevelbox, iicolbox, iirowbox, iichildbox, iiparentbox
      integer inblevel, iistartlev, iiboxlev, iicolleagbox
      integer ineighbors, innbrs, ifcolleag1, ifnbr1

      ier = 0

c     get indeces

      call lrt2d_unpack1(itree,nboxes,nlev,ifcolleag1,ifnbr1,
     1     ilevelbox,iicolbox,iirowbox,iichildbox,iiparentbox,
     2     iiboxlev,inblevel,iistartlev,iicolleagbox,ineighbors,
     3     innbrs)

      if (ifquery .eq. 1) return

      if (ifcolleag .eq. 1 .and. ifcolleag1 .ne. 1) ier = 4
      if (ifnbr .eq. 1 .and. ifnbr1 .ne. 1) ier = ier + 8
      
      if (ier .ne. 0) return


c     copy over arrays

      call lrt2d_icopy(levelbox,itree(ilevelbox),nboxes)
      call lrt2d_icopy(icolbox,itree(iicolbox),nboxes)
      call lrt2d_icopy(irowbox,itree(iirowbox),nboxes)
      call lrt2d_icopy(ichildbox,itree(iichildbox),4*nboxes)
      call lrt2d_icopy(iparentbox,itree(iiparentbox),nboxes)
      call lrt2d_icopy(iboxlev,itree(iiboxlev),nboxes)
      call lrt2d_icopy(nblevel,itree(inblevel),nlev+1)
      call lrt2d_icopy(istartlev,itree(iistartlev),nlev+1)
      if (ifcolleag .eq. 1) then
         call lrt2d_icopy(icolleagbox,itree(iicolleagbox),9*nboxes)
      endif
      if (ifnbr .eq. 1) then
         call lrt2d_icopy(neighbors,itree(ineighbors),12*nboxes)
         call lrt2d_icopy(nnbrs,itree(innbrs),nboxes)
      endif

      return
      end

      subroutine lrt2d_unpack1(itree,nboxes,nlev,ifcolleag,ifnbr,
     1     ilevelbox,iicolbox,iirowbox,iichildbox,iiparentbox,
     2     iiboxlev,inblevel,iistartlev,iicolleagbox,ineighbors,
     3     innbrs)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c-----Global variables
      integer itree(*)
      integer ifcolleag, ifnbr, nboxes, nlev
      integer ilevelbox, iicolbox, iirowbox, iichildbox, iiparentbox
      integer inblevel, iistartlev, iiboxlev, iicolleagbox
      integer ineighbors, innbrs

c     grab info

      nboxes = itree(1)
      nlev = itree(2)
      ifcolleag = itree(3)
      ifnbr = itree(4)
      ilevelbox = itree(5)
      iicolbox = itree(6)
      iirowbox = itree(7)
      iichildbox = itree(8)
      iiparentbox = itree(9)
      iiboxlev = itree(10)
      inblevel = itree(11)
      iistartlev = itree(12)
      iicolleagbox = itree(13)
      ineighbors = itree(14)
      innbrs = itree(15)

      return
      end

      subroutine lrt2d_icopy(ia,ib,n)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     utility subroutine
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      integer ia(*), ib(*), n
c     local
      integer i

c$OMP PARALLEL DO PRIVATE(i) 
c$OMP& IF(n .gt. 10000)
c$OMP& SCHEDULE(static)
      do i = 1,n
         ia(i) = ib(i)
      enddo
C$OMP END PARALLEL DO

      return
      end


      subroutine lrt2d_divall(levelbox,icolbox,irowbox, nboxes, nlev, 
     1     iparentbox,ichildbox,nblevel,iboxlev,istartlev,maxboxes,
     2     maxlevel,itemparray,ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine subdivides all leaf boxes in the tree
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      integer levelbox(*), icolbox(*), irowbox(*), nboxes, nlev
      integer iparentbox(*), ichildbox(4,*), nblevel(0:1)
      integer iboxlev(*), istartlev(0:1), maxboxes, itemparray(*)
      integer maxlevel, ier
c     local variables
      integer nboxes1, i

      if (nlev+1 .gt. maxlevel) then
         ier = 8
         write(*,*) 'lrt2d_divall: CANNOT DIVIDE, MAXLEVELS EXCEEDED'
         return
      endif

      nboxes1 = nboxes

c     only check existing boxes
      do i = 1,nboxes1
c     subdivide if box is a leaf
         if (ichildbox(1,i) .lt. 0) then
            if (nboxes + 4 .gt. maxboxes) then
               ier = 4
               write(*,*) 'lrt2d_divall: INSUFFICIENT MEMORY'
               write(*,*) 'lrt2d_divall: MAXBOXES = ', maxboxes
               return
            endif
            call lrt2d_subdiv1(i,iparentbox,ichildbox,nboxes,irowbox,
     1           icolbox,levelbox)
         endif
      enddo

c     number of levels has increased by 1
      nlev = nlev+1
      
c     sort boxes into iboxlev
      call lrt2d_sortboxes(levelbox,nboxes,nlev,nblevel,iboxlev,
     1     istartlev)

      return
      end
      
