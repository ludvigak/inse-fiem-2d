cc Copyright (C) 2017: Travis Askham, Leslie Greengard, Zydrunas Gimbutas
cc email: askhamwhat@gmail.com      
cc 
cc This software is being released under a modified FreeBSD license
cc (see licenses folder in home directory). 
C***********************************************************************
      subroutine mbhfmm2d_form(beta,ier,iprec,nlev,levelbox,iparentbox, 
     1     ichildbox, icolbox, irowbox, nboxes, nblevel, iboxlev, 
     2     istartlev, ifalltarg, localonoff, zll, blength, ns, srcsort, 
     3     isrcladder, ifcharge, chargesort, ifdipole, dipstrsort, 
     4     dipvecsort, ifquad, quadstrsort, quadvecsort, ifoct, 
     5     octstrsort, octvecsort, isave, lisave, dsave, ldsave, 
     6     csave, lcsave)
c***********************************************************************
c
c     Form the multipole expansions for the modified biharmonic Green's
c     function
c      
c         G(x,y) = (S_beta(r) - S_0(r))/beta^2
c
c     where S_beta(r) = K_0(beta r)/(2*pi) is the fundamental solution
c     of the Yukawa equation and S_0(r) = log(r)/(2*pi) is the
c     fundamental solution of the Laplace equation, with
c     r = sqrt( (x(1)-y(1))**2 + (x(2)-y(2))**2 )
c
c     The expansions represent the sum
c
c     u(x) = sum_j charge(j)*G(x,y(j))
c     + dipstr(j)*(G_{y1}(x,y(j))*dipvec(1,j) + G_{y2}*dipvec(2,j))
c     + quadstr(j)*(G_{y1y1}(x,y(j))*quadvec(1,j)
c             + G_{y1y2}*quadvec(2,j) + G_{y2y2}*quadvec(2,j))
c     + octstr(j)*(G_{y1y1y1}(x,y(j))*octvec(1,j)
c             + G_{y1y1y2}*octvec(2,j) + G_{y1y2y2}*octvec(2,j)
c             + G_{y1y2y2}*octvec(2,j))
c     
c     where the subscripts denote derivatives of the Green's function
c     and the y(j) are points in R^2 contained in the srcsort input
c     vector (the word sort is dropped from the charge, dipstr, etc
c     vectors above for notational convenience).
c
c     On input, the sources should be in a sorted order so that
c     setting ii = isrcladder(1,j) and jj = isrcladder(2,j) means
c     that sources ii through jj are in box i = iboxlev(j) of the tree
c     (jj < ii means that box i has no sources). Note that chargesort,
c     dipstrsort, etc are assumed to also be sorted in this way.
c
c     INPUT:
c     
C 
C     BETA - real *8, modified biharmonic parameter
C     ICOLBOX, IROWBOX, NBOXES, LEVELBOX, IPARENTBOX, IBOXLEV,
c     NBLEVEL, ISTARTLEV, ICHILDBOX - all 
C     define the tree data structure in the standard lrtree format.
c     See lrtree.f for construction routines
c     IFALLTARG - if this flag is set to one, it assumes targets
c     in later evaluations may appear in any leaf box
c     LOCALONOFF - if ifalltarg is not set to one, localonoff is
c     a per-box flag indicating whether local expansions
c     are desired for that box (so that evaluations can be performed
c     there later).
c     ZLL - the bottom left of the highest level box in the tree
c     BLENGTH - the side length of the highest level box in tree
c     NS - the number of sources
c     SRCSORT - real *8 array (2,NS) containing the sources in
c     sorted order as described above
c     ISRCLADDER - as described above
c     IFCHARGE - include charges in multipole expansions
c     CHARGESORT - real *8 array(NS) sorted charge strengths,
c     as in sum above
c     IFDIPOLE - include dipole in multipole expansions
c     DIPSTRSORT - real *8 array(NS) sorted dipole strengths,
c     as in sum above
c     DIPVECSORT - real *8 array(2, NS)sorted dipole direction
c     vector, as in sum above      
c     IFQUAD - include quadrupole in multipole expansions
c     QUADSTRSORT - real *8 array(NS) sorted quad strengths,
c     as in sum above
c     QUADVECSORT - real *8 array(3,NS) sorted quad direction
c     vector, as in sum above      
c     IFOCT - include octopole in multipole expansions
c     OCTSTRSORT - real *8 array(NS) sorted octopole strengths,
c     as in sum above
c     OCTVECSORT - real *8 array(4,NS) sorted octopole direction
c     vector, as in sum above
c     LISAVE, LDSAVE, LCSAVE - the lengths of the ISAVE,
c     DSAVE, and CSAVE output arrays, respectively (which
c     are integer, real, and complex type, respectively)
c     LISAVE, LDSAVE, LCSAVE - if any of these values are
c     negative on input, the subroutine only runs in a query mode
c     (no multipoles formed) and on return the negative values
c     are replaced with the minimum length needed for ISAVE,
c     DSAVE, CSAVE, accordingly
c
c     On output, the vectors isave, dsave, csave contain
c     the full description required to evaluate the sum
c     at any point within the top level box using MBHFMM2D_TARG
c
c     OUTPUT:
c
c     IER - error flag:
c     0 = normal operation (this is not super robust for now...)
c     1 = LISAVE insufficient
c     1 = LDSAVE insufficient
c     3 = LCSAVE insufficient
c     ISAVE, DSAVE, CSAVE - if not in query mode these arrays
c     contain the multipole expansions needed by the evaluation
c     routines MBHFMM2D_TARG, etc.
c     LISAVE, LDSAVE, LCSAVE - if any of these values are
c     negative on input, the subroutine only runs in a query mode
c     (no multipoles formed) and on return the negative values
c     are replaced with the minimum length needed for ISAVE,
c     DSAVE, CSAVE, accordingly
      
      
C     The following subroutine carves up the workspace array before 
C     it is sent to the mbhfmm2d_form1 subroutine, where the actual work 
C     is done.
      
C***********************************************************************
      implicit none
C-----Global variables
      integer nlev,ns,ier,iprec,lisave,ldsave,lcsave,nboxes
      integer levelbox(*), iparentbox(*)
      integer ichildbox(4,*), ifalltarg, localonoff(*)
      integer irowbox(*), icolbox(*)
      integer nblevel(0:nlev), iboxlev(*)
      integer istartlev(0:nlev), isrcladder(2,*)
      integer ifcharge, ifdipole, ifquad, ifoct
      integer isave(*)
      real *8 dsave(*)
      complex *16 csave(*)
      real *8 srcsort(2,*)
      real *8 zll(2), blength, beta
      real *8 chargesort(*), dipstrsort(*), quadstrsort(*)
      real *8 dipvecsort(2,*), quadvecsort(3,*)
      real *8 octstrsort(*), octvecsort(4,*)
c-----Local variables
      integer nterms, nnodes, lmpole, mmbhmpole, mympole
      integer itot, lrscales, mrscales
      integer mcolleag, lcolleag
      integer msrcflag, mflagghost, lflag
      integer mnnbrs, mneighbors
      integer lusedi, lusedd, lusedc, lneighbors, lnnbrs, lbigexp
      integer llocexp, i, mlloc, mmbhloc, mbiglloc, mbigmbhloc
      integer, allocatable :: iwork(:)
      real *8, allocatable :: work(:)
      complex *16, allocatable :: cwork(:)
C
      ier = 0
      
      if(iprec .eq. 0)then
        nterms = 8
        nnodes = 8
      elseif(iprec .eq. 1)then
        nterms = 16
        nnodes = 16
      elseif(iprec .eq. 2)then
        nterms = 24
        nnodes = 24
      elseif(iprec .eq. 3)then
        nterms = 40
        nnodes = 40
      elseif(iprec .eq. 4)then
        nterms = 52
        nnodes = 52
      endif

c     carve up isave

      lflag = nboxes
      lneighbors = 12*nboxes
      lnnbrs = nboxes

      msrcflag = 101
      mflagghost = msrcflag + lflag
      mneighbors = mflagghost + lflag
      mnnbrs = mneighbors + lneighbors
      lusedi = mnnbrs + lnnbrs
      
c     carve up dsave

      lrscales = nlev+5

      mrscales = 1
      lusedd = mrscales + lrscales

c     carve up csave

      llocexp = nboxes*(nterms+1)
      lbigexp = nboxes*(nterms+1)*4

      mlloc = 1
      mmbhloc  = mlloc + llocexp
      mbiglloc = mmbhloc + llocexp
      mbigmbhloc = mbiglloc + lbigexp
      lusedc = mbigmbhloc + lbigexp

      if (lisave .lt. 0 .or. ldsave .lt. 0 .or. lcsave .lt. 0) then

         if (lisave .le. 0) lisave = lusedi
         if (ldsave .le. 0) ldsave = lusedd
         if (lcsave .le. 0) lcsave = lusedc

         return
      endif

      if (lisave .lt. lusedi) then
         ier = 1
         return
      endif
      if (ldsave .lt. lusedd) then
         ier = 2
         return
      endif
      if (lcsave .lt. lusedc) then
         ier = 3
         return
      endif
      

c     save pointers

      isave(1) = nterms
      isave(2) = nnodes
      isave(3) = 0
      isave(11) = msrcflag
      isave(12) = mflagghost
      isave(13) = 0
      isave(14) = 0
      isave(15) = 0
      isave(16) = mneighbors
      isave(17) = mnnbrs

      isave(21) = mrscales

      isave(31) = mlloc
      isave(32) = mmbhloc
      isave(33) = mbiglloc
      isave(34) = mbigmbhloc

c     carve up integer workspace

      lcolleag   =  9*nboxes
      mcolleag   = 1
      itot = mcolleag+lcolleag
      allocate(iwork(itot))

c     carve up complex workspace

      lmpole   = (nterms+1)*nboxes
c
      mympole    = 1
      mmbhmpole  = mympole + lmpole
      itot     = mmbhmpole + lmpole

      allocate(cwork(itot))

C     carve up real workspace
c     without planewaves, there is no real workspace

c     create colleagues 

      call lrt2d_mkcolls(icolbox,irowbox,iwork(mcolleag),nboxes,nlev,
     2      iparentbox,ichildbox,nblevel,iboxlev,istartlev)

c     create neighbors

      call lrt2d_mknbrs(isave(mneighbors),isave(mnnbrs),nboxes,
     1     ichildbox,iparentbox,iwork(mcolleag),icolbox,irowbox)

c     perform fmm forming stage

      call mbhfmm2d_form1(beta,nlev,levelbox,iparentbox,ichildbox,
     1     icolbox,irowbox,iwork(mcolleag),nboxes,nblevel,
     2     iboxlev,istartlev,nterms,cwork(mympole),cwork(mmbhmpole),
     3     csave(mlloc),csave(mmbhloc),csave(mbigmbhloc),
     4     csave(mbiglloc),isave(mflagghost),dsave(mrscales),
     7     ifalltarg,localonoff,isave(msrcflag),ns,srcsort,isrcladder,
     8     ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     9     ifquad,quadstrsort,quadvecsort,ifoct,octstrsort,octvecsort,
     1     zll,blength)


      return
      end
C
c********************************************************************
c      subroutine mbhfmm2d_form1
c********************************************************************
C
C********************************************************************
      subroutine mbhfmm2d_form1(beta,nlev,levelbox,iparentbox,
     1     ichildbox,icolbox,irowbox,icolleagbox,nboxes,nblevel,
     2     iboxlev,istartlev,nterms,ymp,mbhmp,
     3     lloc,mbhloc,bigmbhloc,biglloc,iflagghost,rscales,
     7     ifalltarg,localonoff,isrcflag,ns,srcsort,isrcladder,
     8     ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     9     ifquad,quadstrsort,quadvecsort,ifoct,octstrsort,octvecsort,
     1     zll,blength)

c$    use omp_lib

      implicit none
c-----Global variables
      integer  nlev,nterms,ns,ifalltarg
      integer nnodes, nboxes
      integer icolbox(*), irowbox(*)
      integer nblevel(0:1), iboxlev(*)
      integer istartlev(0:1)
      integer levelbox(*), iparentbox(*)
      integer ichildbox(4,*), icolleagbox(9,*)
      integer iflagghost(*)
      integer localonoff(*), isrcflag(*), isrcladder(2,*)
      integer ifcharge, ifdipole, ifquad, ifoct
      real *8 zll(2), blength, srcsort(2,*), beta
      real *8 rscales(0:nlev+1)
      complex *16 ymp(0:nterms,nboxes),mbhmp(0:nterms,nboxes)
      complex *16 lloc(0:nterms,nboxes),mbhloc(0:nterms,nboxes)
      complex *16 biglloc(0:nterms,4,nboxes)
      complex *16 bigmbhloc(0:nterms,4,nboxes)
      real *8 chargesort(*), dipstrsort(*), quadstrsort(*)
      real *8 dipvecsort(2,*), quadvecsort(3,*)
      real *8 octstrsort(*), octvecsort(4,*)
c-----Local variables
      integer ifar1, ifar2, ifar3
      integer iclose1, iclose2
      integer i, j, ip
      integer iout, ii, jj, ier
      integer ic1,ic2,ic3,ic4
      integer inall(6),iynall(6),in34(4),iy34(4)
      integer isall(6),iysall(6),is12(4),iy12(4)
      integer nnall,nn34,nsall,ns12
      integer ieall(4),iyeall(4),ie14(2),iy14(2)
      integer iee4(1),iey4(1),iee1(1),iey1(1)
      integer iwall(4),iywall(4),iw23(2),iy23(2)
      integer iww3(1),iwy3(1),iww2(1),iwy2(1)
      integer neall,ne14,nee4,nee1,nwall,nw23,nww3,nww2
      integer nside, nsidemark
      integer inbig34(3), isbig12(3)
      integer iebig14(1), iwbig23(1)
      integer iebig1(1), iwbig2(1)
      integer iebig4(1), iwbig3(1)
      integer nb, istart, iend
      integer idbg, idbg2
      integer nsbox, istartbox, irow, icol, iendbox
      integer ntermsall
      integer iz4(2), isrcc1, isrcc2, isrcc3, isrcc4
      integer, allocatable :: iflageast(:), iflagwest(:),
     1     iflagnorth(:), iflagsouth(:)
      real *8 xlength, xlengthc
      real *8 t(7)
      real *8 xsum, pi2, zero
      real *8 scaletemp
      real *8 x
      real *8 a(10,16)
      real *8 time0, time1, second, zbox(2)
      integer ntermsmax, mcarray
      parameter (ntermsmax = 60, mcarray = 120)
c     precomputes for childpar and parchild
      real *8 dfac(0:ntermsmax+ntermsmax), dfac2(0:ntermsmax+ntermsmax)
      real *8 carray(0:mcarray,0:mcarray), ival(0:ntermsmax+ntermsmax)
      real *8 diffs(0:ntermsmax+ntermsmax), pow(0:ntermsmax+ntermsmax)
      real *8 pow2(0:ntermsmax+ntermsmax)
c     scalings
      real *8, allocatable :: scale(:)
c     precomputes for mploc
      real *8, allocatable :: dfac2all(:,:,:,:), kvecall(:,:,:,:)
      real *8, allocatable :: diffsall(:,:,:,:), powall(:,:,:,:)
      complex *16 zshift, zpot
      complex *16 ftarget1, ftarget2, ftarget3, imag
      complex *16 b(0:ntermsmax)
      complex *16 mbhmpc1(0:ntermsmax), ympc1(0:ntermsmax)
      complex *16 mbhmpc2(0:ntermsmax), ympc2(0:ntermsmax)
      complex *16 mbhmpc3(0:ntermsmax), ympc3(0:ntermsmax)
      complex *16 mbhmpc4(0:ntermsmax), ympc4(0:ntermsmax)
      complex *16 spin
      data zero/0.0d0/
      data imag/(0.0d0,1.0d0)/
      logical ifmadebtos
c
      integer iiii, jjjj, jcol, jrow
      integer ifp
      real *8 rscale, sum1, times(20)
c
c
#ifdef VERBOSE
      ifp = 1
#else
      ifp = 0      
#endif

      allocate(iflagnorth(nboxes),iflageast(nboxes),iflagsouth(nboxes),
     1     iflagwest(nboxes))
      
      allocate(scale(0:nlev+1))

      ntermsall = nterms+nterms+20

      allocate(dfac2all(0:ntermsall,-3:3,-3:3,0:nlev))
      allocate(kvecall(0:ntermsall,-3:3,-3:3,0:nlev))
      allocate(diffsall(0:ntermsall,-3:3,-3:3,0:nlev))
      allocate(powall(0:ntermsall,-3:3,-3:3,0:nlev))

      time0 = second()
C$    time0 = omp_get_wtime()

      idbg = 1
      idbg2 = 11
      pi2 = 8.0d0*datan(1.0d0)
      scaletemp = 1.0d0
C
C***********************************************************************
C     STEP 1: Initialization and precomputation of various
C     tables
C***********************************************************************

C$OMP PARALLEL DO SCHEDULE(static,1000)
      do ii = 1, nboxes
         iflagghost(ii) = 0
         iflageast(ii) = 0
         iflagwest(ii) = 0
         iflagnorth(ii) = 0
         iflagsouth(ii) = 0
         isrcflag(ii) = 0
      end do
C$OMP END PARALLEL DO

C
C     initialize multipole and local expansions to zero.
C

C$OMP PARALLEL DO PRIVATE(i,j) IF(nboxes*nterms .gt. 100)
C$OMP& SCHEDULE(static)
      do i = 1, nboxes
         do j = 0, nterms
            ymp(j,i) = zero 
            mbhmp(j,i) = zero 
            lloc(j,i) = zero 
            mbhloc(j,i) = zero 
         enddo
         do jj = 1,4
            do j = 0,nterms
               biglloc(j,jj,i) = zero 
               bigmbhloc(j,jj,i) = zero 
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO

      time1 = second()

c
c     Initialize LOCALONOFF switch to the correct values
c
      if (ifalltarg .eq. 1) then
         do ii = 1, nboxes
            localonoff(ii) = 1
         end do
      endif

C     Call precomputation routines
C     needed in adapfmm, and mkshifts routines.

C
C     compute arrays ZS, etc needed for shifting plane wave exps
C

C
C     create SCALE array
C
      x = 1.0d0 / blength
      do i = 0,nlev+1
         scale(i) = x
         rscales(i) = min(beta/x,1.0d0)
         x = x*2.0d0
      enddo        

C$OMP PARALLEL DO PRIVATE(xlength,ier)
C$OMP& SCHEDULE(static)
      do i = 0,nlev
         xlength = 1.0d0/scale(i)
         call mbh2dmplocall_pre(beta,rscales(i),nterms,rscales(i),
     1        nterms,xlength,ntermsall,dfac2all(0,-3,-3,i),
     2        kvecall(0,-3,-3,i),diffsall(0,-3,-3,i),powall(0,-3,-3,i),
     3        ier)
      enddo
C$OMP END PARALLEL DO

      time1 = second()
c$    time1 = omp_get_wtime()
      
      if (ifp.eq.1) write(*,*) 'wall time, intialize: ', time1-time0


C**********************************************************************
C     STEP 2:  UPWARD PASS
C**********************************************************************
      
      time0 = second()
c$    time0 = omp_get_wtime()

      do i = nlev, 0, -1
         xlength = 1.0d0/scale(i)
         istart = istartlev(i)
         iend = istart + nblevel(i) - 1
         call mbh2d_childpar_pre(beta,xlength,nterms,nterms,
     1        rscales(i+1),rscales(i),dfac,dfac2,carray,mcarray,
     2        ival,diffs,pow)
C$OMP PARALLEL DO PRIVATE(ii,j,ic1,ic2,ic3,ic4)
C$OMP& PRIVATE(istartbox,iendbox,nsbox,zbox)
C$OMP&  IF(iend-istart .gt. 10)
C$OMP&  SCHEDULE(static,1)
         do ii = istart,iend
            j = iboxlev(ii)
            if(ichildbox(1,j) .lt. 0)then
               istartbox = isrcladder(1,j)
               iendbox = isrcladder(2,j)
               nsbox = iendbox-istartbox+1
               if (nsbox .gt. 0) then
                  irow = irowbox(j)
                  icol = icolbox(j)
                  zbox(1) = zll(1)+xlength*(icol-0.5d0)
                  zbox(2) = zll(2)+xlength*(irow-0.5d0)
                  isrcflag(j) = 1
                  call mbh2dformmp_all(ier,beta,rscales(i),
     1                 srcsort(1,istartbox),ifcharge,
     2                 chargesort(istartbox),ifdipole,
     3                 dipstrsort(istartbox),dipvecsort(1,istartbox),
     4                 ifquad,quadstrsort(istartbox),
     5                 quadvecsort(1,istartbox),ifoct,
     6                 octstrsort(istartbox),octvecsort(1,istartbox),
     7                 nsbox,zbox,nterms,mbhmp(0,j),ymp(0,j))
               endif
            elseif(ichildbox(1,j) .gt. 0)then
               ic1 = ichildbox(1,j)
               ic2 = ichildbox(2,j)
               ic3 = ichildbox(3,j)
               ic4 = ichildbox(4,j)
               if (isrcflag(ic1) .eq. 1 .or. isrcflag(ic2) .eq. 1 .or.
     1              isrcflag(ic3) .eq. 1 .or. isrcflag(ic4) .eq. 1) then
                  isrcflag(j) = 1
                  call mbh2d_childpar(mbhmp(0,j),ymp(0,j),mbhmp(0,ic1),
     1                 ymp(0,ic1),mbhmp(0,ic2),ymp(0,ic2),mbhmp(0,ic3),
     2                 ymp(0,ic3),mbhmp(0,ic4),ymp(0,ic4),nterms,nterms,
     3                 rscales(i+1),rscales(i),beta,xlength,dfac,dfac2,
     4                 carray,mcarray,ival,diffs,pow,isrcflag(j))
               endif
            endif
         enddo
C$OMP END PARALLEL DO
      enddo

      time1 = second()
c$    time1 = omp_get_wtime()

      if(ifp.eq.1) write(*,*) 'wall time, upward pass: ', time1-time0     


C
C**********************************************************************
C     STEP 3: DOWNWARD PASS (MULTIPOLES)
C**********************************************************************



c     this goes over boxes so that none of the interactions
c     within a loop affects the same piece of data

      time0 = second()
c$    time0 = omp_get_wtime()

      do iiii = 0,2
      do jjjj = 0,2
C$OMP PARALLEL DO PRIVATE(ii,j,istartbox,iendbox,nsbox,i)
C$OMP& PRIVATE(irow,icol,zbox)
C$OMP& PRIVATE(jrow,jcol,xlength,xlengthc)
C$OMP& PRIVATE(nb,iout,ic1,ic2,ic3,ic4,ifar1,ifar2,ifar3)
C$OMP& PRIVATE(spin,ftarget1,ftarget2,ftarget3,ifmadebtos)
C$OMP& PRIVATE(inall,nnall,iynall,in34)
C$OMP& PRIVATE(nn34,iy34,isall,nsall,iysall,is12,ns12,iy12)
C$OMP& PRIVATE(ieall,neall,iyeall,ie14,ne14,iy14)
C$OMP& PRIVATE(iwall,nwall,iywall,iw23,nw23,iy23)
C$OMP& PRIVATE(iww3,iwy3,nww3,iww2,iwy2,nww2)
C$OMP& PRIVATE(iee4,iey4,nee4,iee1,iey1,nee1)
C$OMP& PRIVATE(inbig34,isbig12,iebig14,iwbig23)
C$OMP& PRIVATE(iebig1, iwbig2, iebig4, iwbig3)
C$OMP& PRIVATE(iz4,mbhmpc1,ympc1,mbhmpc2,ympc2)
C$OMP& PRIVATE(mbhmpc3,ympc3,mbhmpc4,ympc4,isrcc1,isrcc2)
C$OMP& PRIVATE(isrcc3,isrcc4)        
C$OMP&  IF(nboxes .gt. 1000)
C$OMP&  SCHEDULE(dynamic,100)
         do ii = 1,nboxes
            j = iboxlev(ii)

            jrow = irowbox(j)
            jcol = icolbox(j)

            xlength = 1.0d0/scale(levelbox(j))
            xlengthc = xlength/2.0d0

            i = levelbox(j)

            if ( mod(jrow,3) .eq. iiii .and. mod(jcol,3) .eq. jjjj) then

               if(ichildbox(1,j) .gt. 0)then
c     
c     Box has children: 

c     b) create plane wave expansions for processing
c     c) create adaptive lists for plane wave processing
c     d) process all plane wave directions
c     both for same level and for small-to-big far -> big lists

                  ic1 = ichildbox(1,j)
                  ic2 = ichildbox(2,j)
                  ic3 = ichildbox(3,j)
                  ic4 = ichildbox(4,j)

c     if none of the children contain any sources, don't send!!!

                  if (isrcflag(ic1) .eq. 1 .or. isrcflag(ic2) .eq. 1 
     1                 .or. isrcflag(ic3) .eq. 1 
     2                 .or. isrcflag(ic4) .eq. 1) then

                     call lrt2d_mklists(j,inall,nnall,iynall,in34,
     1                    nn34,iy34,isall,nsall,iysall,is12,ns12,iy12,
     3                    ieall,neall,iyeall,ie14,ne14,iy14,
     4                    iwall,nwall,iywall,iw23,nw23,iy23,
     5                    iww3,iwy3,nww3,iww2,iwy2,nww2,
     6                    iee4,iey4,nee4,iee1,iey1,nee1,
     7                    inbig34,isbig12,iebig14,iwbig23,
     8                    iebig4, iwbig3, iebig1, iwbig2,
     9                    icolleagbox,ichildbox,icolbox, irowbox, 
     1                    iflageast, iflagwest, iflagnorth,
     2                    iflagsouth, localonoff)

                     call mbh2dmplocall(beta,rscales(i+1),xlengthc,
     1                    nterms,mbhmp,ymp,mbhloc,lloc,bigmbhloc,
     2                    biglloc,ntermsall,dfac2all(0,-3,-3,i+1),
     3                    kvecall(0,-3,-3,i+1),diffsall(0,-3,-3,i+1),
     3                    powall(0,-3,-3,i+1),j,inall,nnall,
     5                    iynall,in34,nn34,iy34,isall,nsall,iysall,
     6                    is12,ns12,iy12,ieall,neall,iyeall,ie14,
     7                    ne14,iy14,iwall,nwall,iywall,iw23,nw23,
     8                    iy23,iww3,iwy3,nww3,iww2,iwy2,nww2,iee4,
     9                    iey4,nee4,iee1,iey1,nee1,inbig34,isbig12,
     1                    iebig14,iwbig23,iebig4,iwbig3,iebig1,
     2                    iwbig2,ichildbox,localonoff,isrcflag)

                  endif
                  
               elseif (ichildbox(1,j) .lt. 0 
     1                 .and. isrcflag(j) .eq. 1) then

                  istartbox = isrcladder(1,j)
                  iendbox = isrcladder(2,j)
                  nsbox = iendbox-istartbox+1
                  irow = irowbox(j)
                  icol = icolbox(j)
                  zbox(1) = zll(1)+xlength*(icol-0.5d0)
                  zbox(2) = zll(2)+xlength*(irow-0.5d0)

c     
c     Box is childless: check colleagues.
c     1) if colleague is childless, it's a neighbor
c     2) if colleague has children,
c     a) compute big to small far contributions TO colleagues

                  ifmadebtos = .false.

                  do nb = 1, 9
                     iout = icolleagbox(nb,j)
                     if(iout .lt. 0) goto 250
                     if(ichildbox(1,iout).gt.0 
     1                    .and. localonoff(iout) .eq. 1) then
c     
c     2) colleague has children
c     
                        ic1 = ichildbox(1,iout)
                        ic2 = ichildbox(2,iout)
                        ic3 = ichildbox(3,iout)
                        ic4 = ichildbox(4,iout)
c     
C     Form the four expansions needed for BTOS
C     
                        if (.not. ifmadebtos) then

                           call mbh2d_formbtos(beta,xlengthc,
     1                          srcsort(1,istartbox),ifcharge,
     2                          chargesort(istartbox),ifdipole,
     3                          dipstrsort(istartbox),
     4                          dipvecsort(1,istartbox),ifquad,
     5                          quadstrsort(istartbox),
     6                          quadvecsort(1,istartbox),ifoct,
     7                          octstrsort(istartbox),
     8                          octvecsort(1,istartbox),nsbox,zbox,
     7                          mbhmpc1,ympc1,mbhmpc2,ympc2,mbhmpc3,
     8                          ympc3,mbhmpc4,ympc4,isrcc1,isrcc2,
     9                          isrcc3,isrcc4,rscales(i+1),nterms)

                           ifmadebtos = .true.
                        endif
c     
c     Perform big to small far
c     
                        if(nb .eq. 1)then
c     lower left corner: one box not well sep, 3 are.
                           ifar1 = ic4
                           ifar2 = ic3
                           ifar3 = ic1
c     Finally do the far work, big to small
                           spin = -imag
                           ftarget1 = (-2.0d0,-2.0d0)
                           iz4(1) = -2
                           iz4(2) = -2
                           call mbh2d_processbtosfar(beta,xlengthc,
     1                          rscales(i+1),nterms,mbhloc(0,ifar1),
     2                          lloc(0,ifar1),iz4,mbhmpc1,ympc1,mbhmpc2,
     3                          ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4                          isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3                          dfac2all(0,-3,-3,i+1),
     1                          kvecall(0,-3,-3,i+1),
     1                          diffsall(0,-3,-3,i+1),
     1                          powall(0,-3,-3,i+1),localonoff(ifar1))

                           ftarget2 = (-1.0d0,-2.0d0)
                           iz4(1) = -1
                           iz4(2) = -2
                           call mbh2d_processbtosfar(beta,xlengthc,
     1                          rscales(i+1),nterms,mbhloc(0,ifar2),
     2                          lloc(0,ifar2),iz4,mbhmpc1,ympc1,mbhmpc2,
     3                          ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4                          isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3                          dfac2all(0,-3,-3,i+1),
     1                          kvecall(0,-3,-3,i+1),
     1                          diffsall(0,-3,-3,i+1),
     1                          powall(0,-3,-3,i+1),localonoff(ifar2))

                           spin = 1.0d0
                           ftarget3 = (-2.0d0,-1.0d0)
                           iz4(1) = -2
                           iz4(2) = -1
                           call mbh2d_processbtosfar(beta,xlengthc,
     1                          rscales(i+1),nterms,mbhloc(0,ifar3),
     2                          lloc(0,ifar3),iz4,mbhmpc1,ympc1,mbhmpc2,
     3                          ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4                          isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3                          dfac2all(0,-3,-3,i+1),
     1                          kvecall(0,-3,-3,i+1),
     1                          diffsall(0,-3,-3,i+1),
     1                          powall(0,-3,-3,i+1),localonoff(ifar3))

                        elseif(nb .eq. 2) then
C     immediately below, two boxes well separated
                           ifar1 = ic4
                           ifar2 = ic3
c     Finally do the far work, big to small
                           spin = -imag
                           ftarget1 = (0.0d0,-2.0d0)
                           iz4(1) = 0
                           iz4(2) = -2
                           call mbh2d_processbtosfar(beta,xlengthc,
     1                          rscales(i+1),nterms,mbhloc(0,ifar1),
     2                          lloc(0,ifar1),iz4,mbhmpc1,ympc1,mbhmpc2,
     3                          ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4                          isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3                          dfac2all(0,-3,-3,i+1),
     1                          kvecall(0,-3,-3,i+1),
     1                          diffsall(0,-3,-3,i+1),
     1                          powall(0,-3,-3,i+1),localonoff(ifar1))

                           ftarget2 = (1.0d0,-2.0d0)
                           iz4(1) = 1
                           iz4(2) = -2
                           call mbh2d_processbtosfar(beta,xlengthc,
     1                          rscales(i+1),nterms,mbhloc(0,ifar2),
     2                          lloc(0,ifar2),iz4,mbhmpc1,ympc1,mbhmpc2,
     3                          ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4                          isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3                          dfac2all(0,-3,-3,i+1),
     1                          kvecall(0,-3,-3,i+1),
     1                          diffsall(0,-3,-3,i+1),
     1                          powall(0,-3,-3,i+1),localonoff(ifar2))

                        elseif(nb .eq. 3)then
c     lower right corner, one box not well sep., 3 are.
                           ifar1 = ic4  
                           ifar2 = ic3 
                           ifar3 = ic2
c     Finally do the far work, big to small
                           spin = -imag
                           ftarget1 = (2.0d0,-2.0d0)
                           iz4(1) = 2
                           iz4(2) = -2
                           call mbh2d_processbtosfar(beta,xlengthc,
     1                          rscales(i+1),nterms,mbhloc(0,ifar1),
     2                          lloc(0,ifar1),iz4,mbhmpc1,ympc1,mbhmpc2,
     3                          ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4                          isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3                          dfac2all(0,-3,-3,i+1),
     1                          kvecall(0,-3,-3,i+1),
     1                          diffsall(0,-3,-3,i+1),
     1                          powall(0,-3,-3,i+1),localonoff(ifar1))

                           ftarget2 = (3.0d0,-2.0d0)
                           iz4(1) = 3
                           iz4(2) = -2
                           call mbh2d_processbtosfar(beta,xlengthc,
     1                          rscales(i+1),nterms,mbhloc(0,ifar2),
     2                          lloc(0,ifar2),iz4,mbhmpc1,ympc1,mbhmpc2,
     3                          ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4                          isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3                          dfac2all(0,-3,-3,i+1),
     1                          kvecall(0,-3,-3,i+1),
     1                          diffsall(0,-3,-3,i+1),
     1                          powall(0,-3,-3,i+1),localonoff(ifar2))

                           spin = -1.0d0
                           ftarget3 = (3.0d0,-1.0d0)
                           iz4(1) = 3
                           iz4(2) = -1
                           call mbh2d_processbtosfar(beta,xlengthc,
     1                          rscales(i+1),nterms,mbhloc(0,ifar3),
     2                          lloc(0,ifar3),iz4,mbhmpc1,ympc1,mbhmpc2,
     3                          ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4                          isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3                          dfac2all(0,-3,-3,i+1),
     1                          kvecall(0,-3,-3,i+1),
     1                          diffsall(0,-3,-3,i+1),
     1                          powall(0,-3,-3,i+1),localonoff(ifar3))

                        elseif(nb .eq. 4)then
C     immediate left, two boxes well separated
                           ifar1 = ic4 
                           ifar2 = ic1
c     Finally do the far work, big to small
                           spin = 1.0d0
                           ftarget1 = (-2.0d0,0.0d0)
                           iz4(1) = -2
                           iz4(2) = 0
                           call mbh2d_processbtosfar(beta,xlengthc,
     1                          rscales(i+1),nterms,mbhloc(0,ifar1),
     2                          lloc(0,ifar1),iz4,mbhmpc1,ympc1,mbhmpc2,
     3                          ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4                          isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3                          dfac2all(0,-3,-3,i+1),
     1                          kvecall(0,-3,-3,i+1),
     1                          diffsall(0,-3,-3,i+1),
     1                          powall(0,-3,-3,i+1),localonoff(ifar1))

                           ftarget2 = (-2.0d0,1.0d0)
                           iz4(1) = -2
                           iz4(2) = 1
                           call mbh2d_processbtosfar(beta,xlengthc,
     1                          rscales(i+1),nterms,mbhloc(0,ifar2),
     2                          lloc(0,ifar2),iz4,mbhmpc1,ympc1,mbhmpc2,
     3                          ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4                          isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3                          dfac2all(0,-3,-3,i+1),
     1                          kvecall(0,-3,-3,i+1),
     1                          diffsall(0,-3,-3,i+1),
     1                          powall(0,-3,-3,i+1),localonoff(ifar2))

                        elseif(nb .eq. 6)then
C     immediate right, two boxes well separated
                           ifar1 = ic2
                           ifar2 = ic3
C     Finally do the far work, big to small
                           spin = -1.0d0
                           ftarget1 = (3.0d0,1.0d0)
                           iz4(1) = 3
                           iz4(2) = 1
                           call mbh2d_processbtosfar(beta,xlengthc,
     1                          rscales(i+1),nterms,mbhloc(0,ifar1),
     2                          lloc(0,ifar1),iz4,mbhmpc1,ympc1,mbhmpc2,
     3                          ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4                          isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3                          dfac2all(0,-3,-3,i+1),
     1                          kvecall(0,-3,-3,i+1),
     1                          diffsall(0,-3,-3,i+1),
     1                          powall(0,-3,-3,i+1),localonoff(ifar1))

                           ftarget2 = (3.0d0,0.0d0)
                           iz4(1) = 3
                           iz4(2) = 0
                           call mbh2d_processbtosfar(beta,xlengthc,
     1                          rscales(i+1),nterms,mbhloc(0,ifar2),
     2                          lloc(0,ifar2),iz4,mbhmpc1,ympc1,mbhmpc2,
     3                          ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4                          isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3                          dfac2all(0,-3,-3,i+1),
     1                          kvecall(0,-3,-3,i+1),
     1                          diffsall(0,-3,-3,i+1),
     1                          powall(0,-3,-3,i+1),localonoff(ifar2))

                        elseif(nb .eq. 7)then
c     upper left corner, one box not well sep., 3 are.
                           ifar1 = ic4
                           ifar2 = ic1
                           ifar3 = ic2
C     Finally do the far work, big to small
                           spin = 1.0d0
                           ftarget1 = (-2.0d0,2.0d0)
                           iz4(1) = -2
                           iz4(2) = 2
                           call mbh2d_processbtosfar(beta,xlengthc,
     1                          rscales(i+1),nterms,mbhloc(0,ifar1),
     2                          lloc(0,ifar1),iz4,mbhmpc1,ympc1,mbhmpc2,
     3                          ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4                          isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3                          dfac2all(0,-3,-3,i+1),
     1                          kvecall(0,-3,-3,i+1),
     1                          diffsall(0,-3,-3,i+1),
     1                          powall(0,-3,-3,i+1),localonoff(ifar1))

                           spin = imag
                           ftarget2 = (-2.0d0,3.0d0)
                           iz4(1) = -2
                           iz4(2) = 3
                           call mbh2d_processbtosfar(beta,xlengthc,
     1                          rscales(i+1),nterms,mbhloc(0,ifar2),
     2                          lloc(0,ifar2),iz4,mbhmpc1,ympc1,mbhmpc2,
     3                          ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4                          isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3                          dfac2all(0,-3,-3,i+1),
     1                          kvecall(0,-3,-3,i+1),
     1                          diffsall(0,-3,-3,i+1),
     1                          powall(0,-3,-3,i+1),localonoff(ifar2))

                           ftarget3 = (-1.0d0,3.0d0)
                           iz4(1) = -1
                           iz4(2) = 3
                           call mbh2d_processbtosfar(beta,xlengthc,
     1                          rscales(i+1),nterms,mbhloc(0,ifar3),
     2                          lloc(0,ifar3),iz4,mbhmpc1,ympc1,mbhmpc2,
     3                          ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4                          isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3                          dfac2all(0,-3,-3,i+1),
     1                          kvecall(0,-3,-3,i+1),
     1                          diffsall(0,-3,-3,i+1),
     1                          powall(0,-3,-3,i+1),localonoff(ifar3))

                        elseif(nb .eq. 8)then
C     immediately above, two boxes well separated 
                           ifar1 = ic1
                           ifar2 = ic2
c     Finally do the far work, big to small
                           spin = imag
                           ftarget1 = (0.0d0,3.0d0)
                           iz4(1) = 0
                           iz4(2) = 3
                           call mbh2d_processbtosfar(beta,xlengthc,
     1                          rscales(i+1),nterms,mbhloc(0,ifar1),
     2                          lloc(0,ifar1),iz4,mbhmpc1,ympc1,mbhmpc2,
     3                          ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4                          isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3                          dfac2all(0,-3,-3,i+1),
     1                          kvecall(0,-3,-3,i+1),
     1                          diffsall(0,-3,-3,i+1),
     1                          powall(0,-3,-3,i+1),localonoff(ifar1))

                           ftarget2 = (1.0d0,3.0d0)
                           iz4(1) = 1
                           iz4(2) = 3
                           call mbh2d_processbtosfar(beta,xlengthc,
     1                          rscales(i+1),nterms,mbhloc(0,ifar2),
     2                          lloc(0,ifar2),iz4,mbhmpc1,ympc1,mbhmpc2,
     3                          ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4                          isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3                          dfac2all(0,-3,-3,i+1),
     1                          kvecall(0,-3,-3,i+1),
     1                          diffsall(0,-3,-3,i+1),
     1                          powall(0,-3,-3,i+1),localonoff(ifar2))

                        elseif(nb .eq. 9)then
c     upper right corner, one box not well sep., 3 are
                           ifar1 = ic1  
                           ifar2 = ic2
                           ifar3 = ic3
c     Finally do the far work, big to small
                           spin = imag
                           ftarget1 = (2.0d0,3.0d0)
                           iz4(1) = 2
                           iz4(2) = 3
                           call mbh2d_processbtosfar(beta,xlengthc,
     1                          rscales(i+1),nterms,mbhloc(0,ifar1),
     2                          lloc(0,ifar1),iz4,mbhmpc1,ympc1,mbhmpc2,
     3                          ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4                          isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3                          dfac2all(0,-3,-3,i+1),
     1                          kvecall(0,-3,-3,i+1),
     1                          diffsall(0,-3,-3,i+1),
     1                          powall(0,-3,-3,i+1),localonoff(ifar1))

                           ftarget2 = (3.0d0,3.0d0)
                           iz4(1) = 3
                           iz4(2) = 3
                           call mbh2d_processbtosfar(beta,xlengthc,
     1                          rscales(i+1),nterms,mbhloc(0,ifar2),
     2                          lloc(0,ifar2),iz4,mbhmpc1,ympc1,mbhmpc2,
     3                          ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4                          isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3                          dfac2all(0,-3,-3,i+1),
     1                          kvecall(0,-3,-3,i+1),
     1                          diffsall(0,-3,-3,i+1),
     1                          powall(0,-3,-3,i+1),localonoff(ifar2))

                           spin = -1.0d0
                           ftarget3 = (3.0d0,2.0d0)
                           iz4(1) = 3
                           iz4(2) = 2
                           call mbh2d_processbtosfar(beta,xlengthc,
     1                          rscales(i+1),nterms,mbhloc(0,ifar3),
     2                          lloc(0,ifar3),iz4,mbhmpc1,ympc1,mbhmpc2,
     3                          ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4                          isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3                          dfac2all(0,-3,-3,i+1),
     1                          kvecall(0,-3,-3,i+1),
     1                          diffsall(0,-3,-3,i+1),
     1                          powall(0,-3,-3,i+1),localonoff(ifar3))

                        endif
                     endif

 250                 continue
                  enddo
               endif
            endif
         enddo
C$OMP END PARALLEL DO
      enddo
      enddo

      time1 = second()
c$    time1 = omp_get_wtime()
      if (ifp.eq.1) write(*,*) 'wall time, mploc: ', time1-time0

         
c     
      xlength = blength         

c
c     a) Push local expansions to children.
c     For boxes with children, convert the exponential expansions
C     to a single local expansion:
c     

      time0 = second()
c$    time0 = omp_get_wtime()

      do i = 0,nlev
c     
         istart = istartlev(i)
         iend = istart + nblevel(i) - 1

         call mbh2d_parchild_pre(beta,xlength,nterms,nterms,
     1        rscales(i+1),rscales(i),dfac,dfac2,carray,mcarray,
     2        ival,diffs,pow,pow2)
         

C$OMP PARALLEL DO PRIVATE(jj,b,ii,ic1,ic2,ic3,ic4,j) 
C$OMP& IF(iend-istart .gt. 10) 
C$OMP& SCHEDULE(static,1)
         do jj = istart,iend
            j = iboxlev(jj)
            if(ichildbox(1,j).gt.0 .and. localonoff(j).eq.1)then
               ic1 = ichildbox(1,j)
               ic2 = ichildbox(2,j)
               ic3 = ichildbox(3,j)
               ic4 = ichildbox(4,j)
               
               call mbh2d_parchild_add(mbhloc(0,j),lloc(0,j),
     1              mbhloc(0,ic1),lloc(0,ic1),
     2              mbhloc(0,ic2),lloc(0,ic2),
     3              mbhloc(0,ic3),lloc(0,ic3),
     4              mbhloc(0,ic4),lloc(0,ic4),
     5              nterms,nterms,rscales(i+1),rscales(i),
     6              beta,xlength,dfac,dfac2,carray,mcarray,
     7              ival,diffs,pow,pow2,localonoff(j))

            endif
         enddo
C$OMP END PARALLEL DO

         xlength = xlength / 2.0d0
      enddo

      time1 = second()
c$    time1 = omp_get_wtime()
      if(ifp.eq.1) write(*,*) 'wall time, downward pass: ', time1-time0

C$OMP PARALLEL DO SCHEDULE(static) IF(nboxes .gt. 1000)
      do i = 1,nboxes
         if (iflageast(i) .eq. 1 .or. iflagwest(i) .eq. 1 
     1        .or. iflagnorth(i) .eq. 1 .or. iflagsouth(i) .eq. 1) then
            iflagghost(i) = 1
         endif
      enddo
C$OMP END PARALLEL DO

      return
      end

      subroutine mbhfmm2d_targ(beta, ier, nlev, levelbox, iparentbox, 
     1     ichildbox, icolbox, irowbox, nboxes, nblevel, iboxlev, 
     2     istartlev, zll, blength, ns, srcsort, isrcladder, 
     3     ifcharge, chargesort, ifdipole, dipstrsort, dipvecsort,
     4     ifquad, quadstrsort, quadvecsort, isave, dsave, csave, 
     5     nt, targ, ifpot, pot, ifgrad, grad, ifhess, hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Given a multipole description as created by MBHFMM2D_FORM in
c     isave, dsave, csave (along with the tree structure and source
c     info as described in MBHFMM2D_FORM), this routine evaluates
c
c     u(x) = sum_j charge(j)*G(x,y(j))
c     + dipstr(j)*(G_{y1}(x,y(j))*dipvec(1,j) + G_{y2}*dipvec(2,j))
c     + quadstr(j)*(G_{y1y1}(x,y(j))*quadvec(1,j)
c             + G_{y1y2}*quadvec(2,j) + G_{y2y2}*quadvec(2,j))
c     + octstr(j)*(G_{y1y1y1}(x,y(j))*octvec(1,j)
c             + G_{y1y1y2}*octvec(2,j) + G_{y1y2y2}*octvec(2,j)
c             + G_{y1y2y2}*octvec(2,j))
c
c     at each x in the (2,NT) array targ.
c
c     INPUT:
c
c     FOR MOST VARIABLES: see MBHFMM2D_FORM
c     NT: number of targets
c     TARG: real *8 (2,NT) array of target locations (should be
c     inside top level box)
c     IFPOT: evaluate potential (u(x))
c     IFGRAD: evaluate gradient (\nabla u(x))
c     IFHESS: evaluate Hessian (u_{x1x1}(x),u_{x1x2}(x),u_{x2x2}(x))
c
c     OUTPUT:
c
c     POT: real *8 (NT) array, POT(i) = u(targ(1:2,i))
c     GRAD: real *8 (2,NT) array,
c     GRAD(1:2,i) = (u_{x1}(targ(1:2,i)),u_{x2}(targ(1:2,i)))
c     HESS: real *8 (3,NT) array, 
c     HESS(1:3,i) = (u_{x1x1}(targ(1:2,i)),u_{x1x2}(targ(1:2,i)),
c                         u_{x1x2}(targ(1:2,i)))      
c
c     TODO:
c
c     Implement error flag, currently IER does nothing...
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      integer nlev, levelbox(*), iparentbox(*), ichildbox(4,*)
      integer icolbox(*), irowbox(*), nboxes, nblevel(0:1), iboxlev(*)
      integer istartlev(0:1)
      integer ier, iprec, ns, isrcladder(2,*), ifcharge,ifdipole,ifquad
      integer isave(*), nt, ifpot, ifgrad, ifhess
      real *8 zll(2), blength, srcsort(2,*), dsave(*), targ(2,*), beta
      real *8 pot(*), grad(2,*), hess(3,*)
      real *8 chargesort(*), dipstrsort(*), dipvecsort(2,*)
      real *8 quadstrsort(*), quadvecsort(3,*)
      complex *16 csave(*)
c     local
      integer nterms, nnodes , msrcflag, mflagghost
      integer mneighbors, mnnbrs, mrscales
      integer mlloc, mmbhloc, mbiglloc, mbigmbhloc

c     grab pointers

      nterms = isave(1)
      nnodes = isave(2)

      msrcflag = isave(11)
      mflagghost = isave(12)
      mneighbors = isave(16)
      mnnbrs = isave(17)

      mrscales = isave(21)

      mlloc = isave(31)
      mmbhloc = isave(32)
      mbiglloc = isave(33)
      mbigmbhloc = isave(34)


c     call main routine

      call mbhfmm2d_targ1(beta,nlev,levelbox,iparentbox,
     1     ichildbox,icolbox,irowbox,nboxes,nblevel,
     2     iboxlev,istartlev,isave(mneighbors),isave(mnnbrs),nterms,
     3     csave(mlloc),csave(mmbhloc),
     3     csave(mbigmbhloc),csave(mbiglloc),isave(mflagghost),
     4     dsave(mrscales),isave(msrcflag),ns,srcsort,isrcladder,
     8     ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     9     ifquad,quadstrsort,quadvecsort,zll,blength,nt,targ,
     7     ifpot,pot,ifgrad,grad,ifhess,hess)

      return
      end



      subroutine mbhfmm2d_targ1(beta,nlev,levelbox,iparentbox,
     1     ichildbox,icolbox,irowbox,nboxes,nblevel,
     2     iboxlev,istartlev,neighbors,nnbrs,nterms,lloc,mbhloc,
     3     bigmbhloc,biglloc,iflagghost,rscales,isrcflag,ns,srcsort,
     4     isrcladder,ifcharge,chargesort,ifdipole,dipstrsort,
     5     dipvecsort,ifquad,quadstrsort,quadvecsort,zll,blength,
     6     nt,targ,ifpot,pot,ifgrad,grad,ifhess,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global 
      integer nlev, levelbox(*), iparentbox(*), ichildbox(4,*)
      integer icolbox(*), irowbox(*), nboxes, nblevel(0:1), iboxlev(*)
      integer istartlev(0:1), isrcladder(2,*), nterms, nnodes
      integer ier, isrcflag(*), iflagghost(*)
      integer neighbors(12,*), nnbrs(*)
      integer ns, nt, ifpot, ifgrad, ifhess, ifcharge, ifdipole, ifquad
      real *8 srcsort(2,*), targ(2,*), zll(2), blength, rscales(0:nlev)
      real *8 beta
      complex *16 mbhloc(0:nterms,*), lloc(0:nterms,*)
      complex *16 bigmbhloc(0:nterms,4,*), biglloc(0:nterms,4,*)
      real *8 chargesort(*), dipstrsort(*), dipvecsort(2,*)
      real *8 quadstrsort(*), quadvecsort(3,*)
      real *8 pot(*), grad(2,*), hess(3,*)
c     local
      integer i, ibox, ii, icol, irow, jsrc, nsbox, ntbox, istartbox, ic
      integer iendbox, ier2, ilev, OMP_MIN_VAL
      parameter (OMP_MIN_VAL = 1000)
      real *8 xlength, zbox(2), zboxc(2), rscale, xt, yt
      real *8 dzero
      data dzero / 0.0d0 /

C$OMP PARALLEL DO PRIVATE(i,ibox,xlength,ii,jsrc,istartbox,iendbox,
C$OMP& nsbox, icol, irow, ntbox, rscale, zbox, zboxc,ier2,xt,yt,ic,
C$OMP& ilev) 
C$OMP& IF(nt .gt. OMP_MIN_VAL)
C$OMP& SCHEDULE(dynamic,100)
      do i = 1,nt

         if (ifpot .eq. 1) pot(i) = dzero
         if (ifgrad .eq. 1) then
            grad(1,i) = dzero
            grad(2,i) = dzero
         endif
         if (ifhess .eq. 1) then
            hess(1,i) = dzero
            hess(2,i) = dzero
            hess(3,i) = dzero
         endif

c     find containing box

         call lrt2d_findbox(ibox,targ(1,i),targ(2,i),icolbox,irowbox,
     1        ichildbox,nlev,nblevel,iboxlev,istartlev,levelbox, 
     2        zll,blength, ier2)

c     process target

         ilev = levelbox(ibox)

         xlength = blength/(2.0d0**ilev)

c     self direct 

         if (isrcflag(ibox) .eq. 1) then
            istartbox = isrcladder(1,ibox)
            iendbox = isrcladder(2,ibox)
            nsbox = iendbox-istartbox+1

            call mbhpotgrad2dall_cdq_add(beta,srcsort(1,istartbox),
     1           nsbox,ifcharge,chargesort(istartbox),ifdipole,
     2           dipstrsort(istartbox),dipvecsort(1,istartbox),
     3           ifquad,quadstrsort(istartbox),quadvecsort(1,istartbox),
     4           targ(1,i),ifpot,pot(i),ifgrad,grad(1,i),
     5           ifhess,hess(1,i))

         endif

c     neighbors: direct

         do ii = 1,nnbrs(ibox)
            jsrc = neighbors(ii,ibox)
            if (isrcflag(jsrc) .eq. 1) then

               istartbox = isrcladder(1,jsrc)
               iendbox = isrcladder(2,jsrc)
               nsbox = iendbox-istartbox+1

               call mbhpotgrad2dall_cdq_add(beta,srcsort(1,istartbox),
     1              nsbox,ifcharge,chargesort(istartbox),ifdipole,
     2              dipstrsort(istartbox),dipvecsort(1,istartbox),
     3              ifquad,quadstrsort(istartbox),
     4              quadvecsort(1,istartbox),targ(1,i),ifpot,pot(i),
     5              ifgrad,grad(1,i),ifhess,hess(1,i))

            endif
         enddo
            
c     contributions from small to big far

         if (iflagghost(ibox) .eq. 1) then

            xt = targ(1,i)
            yt = targ(2,i)

            ntbox = 1
            irow = irowbox(ibox)
            icol = icolbox(ibox)
            zbox(1) = zll(1)+xlength*(icol-0.5d0)
            zbox(2) = zll(2)+xlength*(irow-0.5d0)

c     figure out which ghost child is relevant

            if (xt .lt. zbox(1)) then
               if (yt .gt. zbox(2)) then
                  ic = 1
                  zboxc(1) = zbox(1)-xlength/4.0d0
                  zboxc(2) = zbox(2)+xlength/4.0d0
               else
                  ic = 4
                  zboxc(1) = zbox(1)-xlength/4.0d0
                  zboxc(2) = zbox(2)-xlength/4.0d0
               endif
            else
               if (yt .gt. zbox(2)) then
                  ic = 2
                  zboxc(1) = zbox(1)+xlength/4.0d0
                  zboxc(2) = zbox(2)+xlength/4.0d0
               else
                  ic = 3
                  zboxc(1) = zbox(1)+xlength/4.0d0
                  zboxc(2) = zbox(2)-xlength/4.0d0
               endif
            endif

            call mbh2dtaevalall(beta,rscales(ilev+1),zboxc,
     1           bigmbhloc(0,ic,ibox),biglloc(0,ic,ibox),
     1           nterms,targ(1,i),ntbox,ifpot,pot(i),ifgrad,grad(1,i),
     3           ifhess,hess(1,i))


            
         endif

c     evaluate local expansions

         ntbox = 1
         irow = irowbox(ibox)
         icol = icolbox(ibox)
         zbox(1) = zll(1)+xlength*(icol-0.5d0)
         zbox(2) = zll(2)+xlength*(irow-0.5d0)

         call mbh2dtaevalall(beta,rscales(ilev),zbox,mbhloc(0,ibox),
     1        lloc(0,ibox),nterms,targ(1,i),ntbox,ifpot,pot(i),
     2        ifgrad,grad(1,i),ifhess,hess(1,i))

      enddo
C$OMP END PARALLEL DO

      return
      end

      subroutine mbhfmm2d3_targ(beta, ier, nlev, levelbox, iparentbox, 
     1     ichildbox, icolbox, irowbox, nboxes, nblevel, iboxlev, 
     2     istartlev, zll, blength, ns, srcsort, isrcladder, 
     3     ifcharge, chargesort, ifdipole, dipstrsort, dipvecsort,
     4     ifquad, quadstrsort, quadvecsort, ifoct, octstrsort,
     5     octvecsort, isave, dsave, csave, 
     5     nt, targ, ifpot, pot, ifgrad, grad, ifhess, hess,
     6     ifder3, der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     SAME AS MBHFMM2D_TARG but also includes the ability to evaluate
c     3rd order derivatives
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      integer nlev, levelbox(*), iparentbox(*), ichildbox(4,*)
      integer icolbox(*), irowbox(*), nboxes, nblevel(0:1), iboxlev(*)
      integer istartlev(0:1), ifoct
      integer ier, iprec, ns, isrcladder(2,*), ifcharge,ifdipole,ifquad
      integer isave(*), nt, ifpot, ifgrad, ifhess, ifder3
      real *8 zll(2), blength, srcsort(2,*), dsave(*), targ(2,*), beta
      real *8 pot(*), grad(2,*), hess(3,*), der3(4,*)
      real *8 chargesort(*), dipstrsort(*), dipvecsort(2,*)
      real *8 quadstrsort(*), quadvecsort(3,*)
      real *8 octstrsort(*), octvecsort(4,*)
      complex *16 csave(*)
c     local
      integer nterms, nnodes , msrcflag, mflagghost
      integer mneighbors, mnnbrs, mrscales
      integer mlloc, mmbhloc, mbiglloc, mbigmbhloc

c     grab pointers

      nterms = isave(1)
      nnodes = isave(2)

      msrcflag = isave(11)
      mflagghost = isave(12)
      mneighbors = isave(16)
      mnnbrs = isave(17)

      mrscales = isave(21)

      mlloc = isave(31)
      mmbhloc = isave(32)
      mbiglloc = isave(33)
      mbigmbhloc = isave(34)


c     call main routine

      call mbhfmm2d3_targ1(beta,nlev,levelbox,iparentbox,
     1     ichildbox,icolbox,irowbox,nboxes,nblevel,
     2     iboxlev,istartlev,isave(mneighbors),isave(mnnbrs),nterms,
     3     csave(mlloc),csave(mmbhloc),
     3     csave(mbigmbhloc),csave(mbiglloc),isave(mflagghost),
     4     dsave(mrscales),isave(msrcflag),ns,srcsort,isrcladder,
     8     ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     9     ifquad,quadstrsort,quadvecsort,ifoct,octstrsort,
     1     octvecsort,zll,blength,nt,targ,
     7     ifpot,pot,ifgrad,grad,ifhess,hess,ifder3,der3)

      return
      end



      subroutine mbhfmm2d3_targ1(beta,nlev,levelbox,iparentbox,
     1     ichildbox,icolbox,irowbox,nboxes,nblevel,
     2     iboxlev,istartlev,neighbors,nnbrs,nterms,lloc,mbhloc,
     3     bigmbhloc,biglloc,iflagghost,rscales,isrcflag,ns,srcsort,
     4     isrcladder,ifcharge,chargesort,ifdipole,dipstrsort,
     5     dipvecsort,ifquad,quadstrsort,quadvecsort,ifoct,octstrsort,
     5     octvecsort,zll,blength,
     6     nt,targ,ifpot,pot,ifgrad,grad,ifhess,hess,ifder3,der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global 
      integer nlev, levelbox(*), iparentbox(*), ichildbox(4,*)
      integer icolbox(*), irowbox(*), nboxes, nblevel(0:1), iboxlev(*)
      integer istartlev(0:1), isrcladder(2,*), nterms, nnodes
      integer ier, isrcflag(*), iflagghost(*)
      integer neighbors(12,*), nnbrs(*)
      integer ns, nt, ifpot, ifgrad, ifhess, ifcharge, ifdipole, ifquad
      integer ifder3, ifoct
      real *8 srcsort(2,*), targ(2,*), zll(2), blength, rscales(0:nlev)
      real *8 beta
      complex *16 mbhloc(0:nterms,*), lloc(0:nterms,*)
      complex *16 bigmbhloc(0:nterms,4,*), biglloc(0:nterms,4,*)
      real *8 chargesort(*), dipstrsort(*), dipvecsort(2,*)
      real *8 quadstrsort(*), quadvecsort(3,*)
      real *8 octstrsort(*), octvecsort(4,*)
      real *8 pot(*), grad(2,*), hess(3,*), der3(4,*)
c     local
      integer i, ibox, ii, icol, irow, jsrc, nsbox, ntbox, istartbox, ic
      integer iendbox, ier2, ilev, OMP_MIN_VAL
      parameter (OMP_MIN_VAL = 1000)
      real *8 xlength, zbox(2), zboxc(2), rscale, xt, yt
      real *8 dzero
      data dzero / 0.0d0 /

C$OMP PARALLEL DO PRIVATE(ibox,xlength,ii,jsrc,istartbox,iendbox,nsbox,
C$OMP& icol, irow, ntbox, rscale, zbox, zboxc, ier2, ic, xt, yt, ilev) 
C$OMP& IF(nt .gt. OMP_MIN_VAL)
C$OMP& SCHEDULE(dynamic,100)
      do i = 1,nt

         if (ifpot .eq. 1) pot(i) = dzero
         if (ifgrad .eq. 1) then
            grad(1,i) = dzero
            grad(2,i) = dzero
         endif
         if (ifhess .eq. 1) then
            hess(1,i) = dzero
            hess(2,i) = dzero
            hess(3,i) = dzero
         endif
         if (ifder3 .eq. 1) then
            der3(1,i) = dzero
            der3(2,i) = dzero
            der3(3,i) = dzero
            der3(4,i) = dzero
         endif

c     find containing box

         call lrt2d_findbox(ibox,targ(1,i),targ(2,i),icolbox,irowbox,
     1        ichildbox,nlev,nblevel,iboxlev,istartlev,levelbox, 
     2        zll,blength, ier2)

c     process target

         ilev = levelbox(ibox)

         xlength = blength/(2.0d0**ilev)

c     self direct 

         if (isrcflag(ibox) .eq. 1) then
            istartbox = isrcladder(1,ibox)
            iendbox = isrcladder(2,ibox)
            nsbox = iendbox-istartbox+1

            call mbhpotgrad2dall_cdqo3_add(beta,srcsort(1,istartbox),
     1           nsbox,ifcharge,chargesort(istartbox),ifdipole,
     2           dipstrsort(istartbox),dipvecsort(1,istartbox),
     3           ifquad,quadstrsort(istartbox),quadvecsort(1,istartbox),
     3           ifoct,octstrsort(istartbox),octvecsort(1,istartbox),
     4           targ(1,i),ifpot,pot(i),ifgrad,grad(1,i),
     5           ifhess,hess(1,i),ifder3,der3(1,i))

         endif

c     neighbors: direct

         do ii = 1,nnbrs(ibox)
            jsrc = neighbors(ii,ibox)
            if (isrcflag(jsrc) .eq. 1) then

               istartbox = isrcladder(1,jsrc)
               iendbox = isrcladder(2,jsrc)
               nsbox = iendbox-istartbox+1

               call mbhpotgrad2dall_cdqo3_add(beta,srcsort(1,istartbox),
     1              nsbox,ifcharge,chargesort(istartbox),ifdipole,
     2              dipstrsort(istartbox),dipvecsort(1,istartbox),
     3              ifquad,quadstrsort(istartbox),
     4              quadvecsort(1,istartbox),ifoct,
     4              octstrsort(istartbox),octvecsort(1,istartbox),
     6              targ(1,i),ifpot,pot(i),
     5              ifgrad,grad(1,i),ifhess,hess(1,i),ifder3,der3(1,i))

            endif
         enddo
            
c     contributions from small to big far

         if (iflagghost(ibox) .eq. 1) then

            xt = targ(1,i)
            yt = targ(2,i)

            ntbox = 1
            irow = irowbox(ibox)
            icol = icolbox(ibox)
            zbox(1) = zll(1)+xlength*(icol-0.5d0)
            zbox(2) = zll(2)+xlength*(irow-0.5d0)

c     figure out which ghost child is relevant

            if (xt .lt. zbox(1)) then
               if (yt .gt. zbox(2)) then
                  ic = 1
                  zboxc(1) = zbox(1)-xlength/4.0d0
                  zboxc(2) = zbox(2)+xlength/4.0d0
               else
                  ic = 4
                  zboxc(1) = zbox(1)-xlength/4.0d0
                  zboxc(2) = zbox(2)-xlength/4.0d0
               endif
            else
               if (yt .gt. zbox(2)) then
                  ic = 2
                  zboxc(1) = zbox(1)+xlength/4.0d0
                  zboxc(2) = zbox(2)+xlength/4.0d0
               else
                  ic = 3
                  zboxc(1) = zbox(1)+xlength/4.0d0
                  zboxc(2) = zbox(2)-xlength/4.0d0
               endif
            endif

            call mbh2dtaeval3all(beta,rscales(ilev+1),zboxc,
     1           bigmbhloc(0,ic,ibox),biglloc(0,ic,ibox),
     1           nterms,targ(1,i),ntbox,ifpot,pot(i),ifgrad,grad(1,i),
     3           ifhess,hess(1,i),ifder3,der3(1,i))


            
         endif

c     evaluate local expansions

         ntbox = 1
         irow = irowbox(ibox)
         icol = icolbox(ibox)
         zbox(1) = zll(1)+xlength*(icol-0.5d0)
         zbox(2) = zll(2)+xlength*(irow-0.5d0)

         call mbh2dtaeval3all(beta,rscales(ilev),zbox,mbhloc(0,ibox),
     1        lloc(0,ibox),nterms,targ(1,i),ntbox,ifpot,pot(i),
     2        ifgrad,grad(1,i),ifhess,hess(1,i),ifder3,der3(1,i))

      enddo
C$OMP END PARALLEL DO

      return
      end

      subroutine mbhfmm2d3_targ_exrad(beta,ier,nlev,levelbox,iparentbox, 
     1     ichildbox, icolbox, irowbox, nboxes, nblevel, iboxlev, 
     2     istartlev, zll, blength, ns, srcsort, isrcladder, exrad,
     3     ifcharge, chargesort, ifdipole, dipstrsort, dipvecsort,
     4     ifquad, quadstrsort, quadvecsort, ifoct, octstrsort,
     5     octvecsort, isave, dsave, csave, 
     5     nt, targ, ifpot, pot, ifgrad, grad, ifhess, hess,
     6     ifder3, der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     SAME AS MBHFMM2D3_TARG but adds the ability to exclude a
c     radius about each source point in *DIRECT* calculations
c
c     EXRAD - radius to exclude per source (in src sorted order)
c
c     Note that if the tree has boxes so small that the influence of
c     a source at a target is computed via multipoles/local expansions
c     then this routine will fail.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      integer nlev, levelbox(*), iparentbox(*), ichildbox(4,*)
      integer icolbox(*), irowbox(*), nboxes, nblevel(0:1), iboxlev(*)
      integer istartlev(0:1), ifoct
      integer ier, iprec, ns, isrcladder(2,*), ifcharge,ifdipole,ifquad
      integer isave(*), nt, ifpot, ifgrad, ifhess, ifder3
      real *8 zll(2), blength, srcsort(2,*), dsave(*), targ(2,*), beta
      real *8 exrad(*)
      real *8 pot(*), grad(2,*), hess(3,*), der3(4,*)
      real *8 chargesort(*), dipstrsort(*), dipvecsort(2,*)
      real *8 quadstrsort(*), quadvecsort(3,*)
      real *8 octstrsort(*), octvecsort(4,*)
      complex *16 csave(*)
c     local
      integer nterms, nnodes , msrcflag, mflagghost
      integer mneighbors, mnnbrs, mrscales
      integer mlloc, mmbhloc, mbiglloc, mbigmbhloc

c     grab pointers

      nterms = isave(1)
      nnodes = isave(2)

      msrcflag = isave(11)
      mflagghost = isave(12)
      mneighbors = isave(16)
      mnnbrs = isave(17)

      mrscales = isave(21)

      mlloc = isave(31)
      mmbhloc = isave(32)
      mbiglloc = isave(33)
      mbigmbhloc = isave(34)


c     call main routine

      call mbhfmm2d3_targ_exrad1(beta,nlev,levelbox,iparentbox,
     1     ichildbox,icolbox,irowbox,nboxes,nblevel,
     2     iboxlev,istartlev,isave(mneighbors),isave(mnnbrs),nterms,
     3     csave(mlloc),csave(mmbhloc),
     3     csave(mbigmbhloc),csave(mbiglloc),isave(mflagghost),
     4     dsave(mrscales),isave(msrcflag),ns,srcsort,isrcladder,
     5     exrad,
     8     ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     9     ifquad,quadstrsort,quadvecsort,ifoct,octstrsort,
     1     octvecsort,zll,blength,nt,targ,
     7     ifpot,pot,ifgrad,grad,ifhess,hess,ifder3,der3)

      return
      end



      subroutine mbhfmm2d3_targ_exrad1(beta,nlev,levelbox,iparentbox,
     1     ichildbox,icolbox,irowbox,nboxes,nblevel,
     2     iboxlev,istartlev,neighbors,nnbrs,nterms,lloc,mbhloc,
     3     bigmbhloc,biglloc,iflagghost,rscales,isrcflag,ns,srcsort,
     4     isrcladder,exrad,ifcharge,chargesort,ifdipole,dipstrsort,
     5     dipvecsort,ifquad,quadstrsort,quadvecsort,ifoct,octstrsort,
     5     octvecsort,zll,blength,
     6     nt,targ,ifpot,pot,ifgrad,grad,ifhess,hess,ifder3,der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global 
      integer nlev, levelbox(*), iparentbox(*), ichildbox(4,*)
      integer icolbox(*), irowbox(*), nboxes, nblevel(0:1), iboxlev(*)
      integer istartlev(0:1), isrcladder(2,*), nterms, nnodes
      integer ier, isrcflag(*), iflagghost(*)
      integer neighbors(12,*), nnbrs(*)
      integer ns, nt, ifpot, ifgrad, ifhess, ifcharge, ifdipole, ifquad
      integer ifder3, ifoct
      real *8 srcsort(2,*), targ(2,*), zll(2), blength, rscales(0:nlev)
      real *8 beta, exrad(*)
      complex *16 mbhloc(0:nterms,*), lloc(0:nterms,*)
      complex *16 bigmbhloc(0:nterms,4,*), biglloc(0:nterms,4,*)
      real *8 chargesort(*), dipstrsort(*), dipvecsort(2,*)
      real *8 quadstrsort(*), quadvecsort(3,*)
      real *8 octstrsort(*), octvecsort(4,*)
      real *8 pot(*), grad(2,*), hess(3,*), der3(4,*)
c     local
      integer i, ibox, ii, icol, irow, jsrc, nsbox, ntbox, istartbox, ic
      integer iendbox, ier2, ilev, OMP_MIN_VAL
      parameter (OMP_MIN_VAL = 1000)
      real *8 xlength, zbox(2), zboxc(2), rscale, xt, yt
      real *8 dzero
      data dzero / 0.0d0 /

C$OMP PARALLEL DO PRIVATE(ibox,xlength,ii,jsrc,istartbox,iendbox,nsbox,
C$OMP& icol, irow, ntbox, rscale, zbox, zboxc, ier2, ic, xt, yt, ilev) 
C$OMP& IF(nt .gt. OMP_MIN_VAL)
C$OMP& SCHEDULE(dynamic,100)
      do i = 1,nt

         if (ifpot .eq. 1) pot(i) = dzero
         if (ifgrad .eq. 1) then
            grad(1,i) = dzero
            grad(2,i) = dzero
         endif
         if (ifhess .eq. 1) then
            hess(1,i) = dzero
            hess(2,i) = dzero
            hess(3,i) = dzero
         endif
         if (ifder3 .eq. 1) then
            der3(1,i) = dzero
            der3(2,i) = dzero
            der3(3,i) = dzero
            der3(4,i) = dzero
         endif

c     find containing box

         call lrt2d_findbox(ibox,targ(1,i),targ(2,i),icolbox,irowbox,
     1        ichildbox,nlev,nblevel,iboxlev,istartlev,levelbox, 
     2        zll,blength, ier2)

c     process target

         ilev = levelbox(ibox)

         xlength = blength/(2.0d0**ilev)

c     self direct 

         if (isrcflag(ibox) .eq. 1) then
            istartbox = isrcladder(1,ibox)
            iendbox = isrcladder(2,ibox)
            nsbox = iendbox-istartbox+1

            call mbhpotgrad2dall_cdqo3_add_exrad(beta,
     1           srcsort(1,istartbox),exrad(istartbox),
     2           nsbox,ifcharge,chargesort(istartbox),ifdipole,
     2           dipstrsort(istartbox),dipvecsort(1,istartbox),
     3           ifquad,quadstrsort(istartbox),quadvecsort(1,istartbox),
     3           ifoct,octstrsort(istartbox),octvecsort(1,istartbox),
     4           targ(1,i),ifpot,pot(i),ifgrad,grad(1,i),
     5           ifhess,hess(1,i),ifder3,der3(1,i))

         endif

c     neighbors: direct

         do ii = 1,nnbrs(ibox)
            jsrc = neighbors(ii,ibox)
            if (isrcflag(jsrc) .eq. 1) then

               istartbox = isrcladder(1,jsrc)
               iendbox = isrcladder(2,jsrc)
               nsbox = iendbox-istartbox+1

               call mbhpotgrad2dall_cdqo3_add_exrad(beta,
     1              srcsort(1,istartbox),exrad(istartbox),
     1              nsbox,ifcharge,chargesort(istartbox),ifdipole,
     2              dipstrsort(istartbox),dipvecsort(1,istartbox),
     3              ifquad,quadstrsort(istartbox),
     4              quadvecsort(1,istartbox),ifoct,
     4              octstrsort(istartbox),octvecsort(1,istartbox),
     6              targ(1,i),ifpot,pot(i),
     5              ifgrad,grad(1,i),ifhess,hess(1,i),ifder3,der3(1,i))

            endif
         enddo
            
c     contributions from small to big far

         if (iflagghost(ibox) .eq. 1) then

            xt = targ(1,i)
            yt = targ(2,i)

            ntbox = 1
            irow = irowbox(ibox)
            icol = icolbox(ibox)
            zbox(1) = zll(1)+xlength*(icol-0.5d0)
            zbox(2) = zll(2)+xlength*(irow-0.5d0)

c     figure out which ghost child is relevant

            if (xt .lt. zbox(1)) then
               if (yt .gt. zbox(2)) then
                  ic = 1
                  zboxc(1) = zbox(1)-xlength/4.0d0
                  zboxc(2) = zbox(2)+xlength/4.0d0
               else
                  ic = 4
                  zboxc(1) = zbox(1)-xlength/4.0d0
                  zboxc(2) = zbox(2)-xlength/4.0d0
               endif
            else
               if (yt .gt. zbox(2)) then
                  ic = 2
                  zboxc(1) = zbox(1)+xlength/4.0d0
                  zboxc(2) = zbox(2)+xlength/4.0d0
               else
                  ic = 3
                  zboxc(1) = zbox(1)+xlength/4.0d0
                  zboxc(2) = zbox(2)-xlength/4.0d0
               endif
            endif

            call mbh2dtaeval3all(beta,rscales(ilev+1),zboxc,
     1           bigmbhloc(0,ic,ibox),biglloc(0,ic,ibox),
     1           nterms,targ(1,i),ntbox,ifpot,pot(i),ifgrad,grad(1,i),
     3           ifhess,hess(1,i),ifder3,der3(1,i))


            
         endif

c     evaluate local expansions

         ntbox = 1
         irow = irowbox(ibox)
         icol = icolbox(ibox)
         zbox(1) = zll(1)+xlength*(icol-0.5d0)
         zbox(2) = zll(2)+xlength*(irow-0.5d0)

         call mbh2dtaeval3all(beta,rscales(ilev),zbox,mbhloc(0,ibox),
     1        lloc(0,ibox),nterms,targ(1,i),ntbox,ifpot,pot(i),
     2        ifgrad,grad(1,i),ifhess,hess(1,i),ifder3,der3(1,i))

      enddo
C$OMP END PARALLEL DO

      return
      end

      subroutine mbhfmm2d3_targrc(beta, ier, nlev, levelbox, iparentbox, 
     1     ichildbox, icolbox, irowbox, nboxes, nblevel, iboxlev, 
     2     istartlev, zll, blength, ns, srcsort, isrcladder, 
     3     ifcharge, chargesort, ifdipole, dipstrsort, dipvecsort,
     4     ifquad, quadstrsort, quadvecsort, ifoct, octstrsort, 
     5     octvecsort, isave, dsave, csave, 
     5     nt, targ, centers, ifpot, pot, ifgrad, grad, ifhess, hess,
     6     ifder3, der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     SAME AS MBHFMM2D3_TARG, targets are described as the sum of
c     targ(1:2,i) and centers(1:2,i). The idea is to avoid cancellation
c     when evaluating targets that lie on a small disc around points
c     in the domain (this was part of a quick and dirty QBX
c     implementation
c     
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      integer nlev, levelbox(*), iparentbox(*), ichildbox(4,*)
      integer icolbox(*), irowbox(*), nboxes, nblevel(0:1), iboxlev(*)
      integer istartlev(0:1)
      integer ier, iprec, ns, isrcladder(2,*), ifcharge,ifdipole,ifquad
      integer isave(*), nt, ifpot, ifgrad, ifhess, ifder3, ifoct
      real *8 zll(2), blength, srcsort(2,*), dsave(*), targ(2,*), beta
      real *8 pot(*), grad(2,*), hess(3,*), der3(4,*), centers(2,*)
      real *8 chargesort(*), dipstrsort(*), dipvecsort(2,*)
      real *8 quadstrsort(*), quadvecsort(3,*), octstrsort(*)
      real *8 octvecsort(4,*)
      complex *16 csave(*)
c     local
      integer nterms, nnodes , msrcflag, mflagghost
      integer mneighbors, mnnbrs, mrscales
      integer mlloc, mmbhloc, mbiglloc, mbigmbhloc

c     grab pointers

      nterms = isave(1)
      nnodes = isave(2)

      msrcflag = isave(11)
      mflagghost = isave(12)
      mneighbors = isave(16)
      mnnbrs = isave(17)

      mrscales = isave(21)

      mlloc = isave(31)
      mmbhloc = isave(32)
      mbiglloc = isave(33)
      mbigmbhloc = isave(34)


c     call main routine

      call mbhfmm2d3_targrc1(beta,nlev,levelbox,iparentbox,
     1     ichildbox,icolbox,irowbox,nboxes,nblevel,
     2     iboxlev,istartlev,isave(mneighbors),isave(mnnbrs),nterms,
     3     csave(mlloc),csave(mmbhloc),
     3     csave(mbigmbhloc),csave(mbiglloc),isave(mflagghost),
     4     dsave(mrscales),isave(msrcflag),ns,srcsort,isrcladder,
     8     ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     9     ifquad,quadstrsort,quadvecsort,ifoct,octstrsort,
     1     octvecsort,zll,blength,nt,targ,centers,
     7     ifpot,pot,ifgrad,grad,ifhess,hess,ifder3,der3)

      return
      end



      subroutine mbhfmm2d3_targrc1(beta,nlev,levelbox,iparentbox,
     1     ichildbox,icolbox,irowbox,nboxes,nblevel,
     2     iboxlev,istartlev,neighbors,nnbrs,nterms,lloc,mbhloc,
     3     bigmbhloc,biglloc,iflagghost,rscales,isrcflag,ns,srcsort,
     4     isrcladder,ifcharge,chargesort,ifdipole,dipstrsort,
     5     dipvecsort,ifquad,quadstrsort,quadvecsort,ifoct,
     6     octstrsort,octvecsort,zll,blength,
     6     nt,targ,centers,ifpot,pot,ifgrad,grad,ifhess,hess,
     7     ifder3,der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global 
      integer nlev, levelbox(*), iparentbox(*), ichildbox(4,*)
      integer icolbox(*), irowbox(*), nboxes, nblevel(0:1), iboxlev(*)
      integer istartlev(0:1), isrcladder(2,*), nterms, nnodes
      integer ier, isrcflag(*), iflagghost(*)
      integer neighbors(12,*), nnbrs(*)
      integer ns, nt, ifpot, ifgrad, ifhess, ifcharge, ifdipole, ifquad
      integer ifder3, ifoct
      real *8 srcsort(2,*), targ(2,*), zll(2), blength, rscales(0:nlev)
      real *8 beta, centers(2,*)
      complex *16 mbhloc(0:nterms,*), lloc(0:nterms,*)
      complex *16 bigmbhloc(0:nterms,4,*), biglloc(0:nterms,4,*)
      real *8 chargesort(*), dipstrsort(*), dipvecsort(2,*)
      real *8 quadstrsort(*), quadvecsort(3,*)
      real *8 octstrsort(*), octvecsort(4,*)
      real *8 pot(*), grad(2,*), hess(3,*), der3(4,*)
c     local
      integer i, ibox, ii, icol, irow, jsrc, nsbox, ntbox, istartbox, ic
      integer j, nsbox1, itemp
      integer iendbox, ier2, ilev, OMP_MIN_VAL
      parameter (OMP_MIN_VAL = 1000)
      real *8 xlength, zbox(2), zboxc(2), rscale, xt, yt, zx(2), zy(2)
      real *8 dzero
      data dzero / 0.0d0 /

C$OMP PARALLEL DO PRIVATE(ibox,xlength,ii,jsrc,istartbox,iendbox,nsbox,
C$OMP& icol,irow,ntbox,rscale,zbox,ier2,zx,zy,j,nsbox1,itemp,zboxc,
C$OMP& xt,yt,ic,ilev) 
C$OMP& IF(nt .gt. OMP_MIN_VAL)
C$OMP& SCHEDULE(dynamic,100)
      do i = 1,nt

         if (ifpot .eq. 1) pot(i) = dzero
         if (ifgrad .eq. 1) then
            grad(1,i) = dzero
            grad(2,i) = dzero
         endif
         if (ifhess .eq. 1) then
            hess(1,i) = dzero
            hess(2,i) = dzero
            hess(3,i) = dzero
         endif
         if (ifder3 .eq. 1) then
            der3(1,i) = dzero
            der3(2,i) = dzero
            der3(3,i) = dzero
            der3(4,i) = dzero
         endif

c     find containing box

         zx(1) = targ(1,i) + centers(1,i)
         zx(2) = targ(2,i) + centers(2,i)

         call lrt2d_findbox(ibox,zx(1),zx(2),icolbox,irowbox,
     1        ichildbox,nlev,nblevel,iboxlev,istartlev,levelbox, 
     2        zll,blength, ier2)

c     process target

         ilev = levelbox(ibox)

         xlength = blength/(2.0d0**ilev)

c     self direct 

         if (isrcflag(ibox) .eq. 1) then
            istartbox = isrcladder(1,ibox)
            iendbox = isrcladder(2,ibox)
            nsbox = iendbox-istartbox+1

            nsbox1 = 1

            do j = 1,nsbox

               itemp = istartbox+j-1

               zy(1) = srcsort(1,itemp)-centers(1,i)
               zy(2) = srcsort(2,itemp)-centers(2,i)

               call mbhpotgrad2dall_cdqo3_add(beta,zy,
     1              nsbox1,ifcharge,chargesort(itemp),ifdipole,
     2              dipstrsort(itemp),dipvecsort(1,itemp),
     3              ifquad,quadstrsort(itemp),quadvecsort(1,itemp),
     4              ifoct,octstrsort(itemp),octvecsort(1,itemp),
     4              targ(1,i),ifpot,pot(i),ifgrad,grad(1,i),
     5              ifhess,hess(1,i),ifder3,der3(1,i))

            enddo
           
         endif

c     neighbors: direct

         do ii = 1,nnbrs(ibox)
            jsrc = neighbors(ii,ibox)
            if (isrcflag(jsrc) .eq. 1) then

               istartbox = isrcladder(1,jsrc)
               iendbox = isrcladder(2,jsrc)
               nsbox = iendbox-istartbox+1

               nsbox1 = 1

               do j = 1,nsbox

                  itemp = istartbox+j-1

                  zy(1) = srcsort(1,itemp)-centers(1,i)
                  zy(2) = srcsort(2,itemp)-centers(2,i)

                  call mbhpotgrad2dall_cdqo3_add(beta,zy,
     1                 nsbox1,ifcharge,chargesort(itemp),ifdipole,
     2                 dipstrsort(itemp),dipvecsort(1,itemp),
     3                 ifquad,quadstrsort(itemp),quadvecsort(1,itemp),
     4                 ifoct,octstrsort(itemp),octvecsort(1,itemp),
     4                 targ(1,i),ifpot,pot(i),ifgrad,grad(1,i),
     5                 ifhess,hess(1,i),ifder3,der3(1,i))

               enddo

            endif
         enddo
            
c     contributions from small to big far

         if (iflagghost(ibox) .eq. 1) then

            xt = targ(1,i)
            yt = targ(2,i)

            ntbox = 1
            irow = irowbox(ibox)
            icol = icolbox(ibox)
            zbox(1) = zll(1)+xlength*(icol-0.5d0)-centers(1,i)
            zbox(2) = zll(2)+xlength*(irow-0.5d0)-centers(2,i)

c     figure out which ghost child is relevant

            if (xt .lt. zbox(1)) then
               if (yt .gt. zbox(2)) then
                  ic = 1
                  zboxc(1) = zbox(1)-xlength/4.0d0
                  zboxc(2) = zbox(2)+xlength/4.0d0
               else
                  ic = 4
                  zboxc(1) = zbox(1)-xlength/4.0d0
                  zboxc(2) = zbox(2)-xlength/4.0d0
               endif
            else
               if (yt .gt. zbox(2)) then
                  ic = 2
                  zboxc(1) = zbox(1)+xlength/4.0d0
                  zboxc(2) = zbox(2)+xlength/4.0d0
               else
                  ic = 3
                  zboxc(1) = zbox(1)+xlength/4.0d0
                  zboxc(2) = zbox(2)-xlength/4.0d0
               endif
            endif

            call mbh2dtaeval3all(beta,rscales(ilev+1),zboxc,
     1           bigmbhloc(0,ic,ibox),biglloc(0,ic,ibox),
     1           nterms,targ(1,i),ntbox,ifpot,pot(i),ifgrad,grad(1,i),
     3           ifhess,hess(1,i),ifder3,der3(1,i))

         endif

c     evaluate local expansions

         ntbox = 1
         irow = irowbox(ibox)
         icol = icolbox(ibox)
         zbox(1) = zll(1)+xlength*(icol-0.5d0)-centers(1,i)
         zbox(2) = zll(2)+xlength*(irow-0.5d0)-centers(2,i)

         call mbh2dtaeval3all(beta,rscales(ilev),zbox,mbhloc(0,ibox),
     1        lloc(0,ibox),nterms,targ(1,i),ntbox,ifpot,pot(i),
     2        ifgrad,grad(1,i),ifhess,hess(1,i),ifder3,der3(1,i))

      enddo
C$OMP END PARALLEL DO

      return
      end


      subroutine mbhfmm2d3_srcsrc(beta, ier, nlev, levelbox, iparentbox, 
     1     ichildbox, icolbox, irowbox, nboxes, nblevel, iboxlev, 
     2     istartlev, zll, blength, ns, srcsort, isrcladder, 
     3     ifcharge, chargesort, ifdipole, dipstrsort, dipvecsort,
     4     ifquad, quadstrsort, quadvecsort, ifoct, octstrsort,
     5     octvecsort, isave, dsave, csave, 
     5     ifpot, pot, ifgrad, grad, ifhess, hess,
     6     ifder3, der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     Computes the sum of all charges at a source location (excluding
c     the self)
c
c     THESE ARE RETURNED IN THE SORTED ORDER
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      integer nlev, levelbox(*), iparentbox(*), ichildbox(4,*)
      integer icolbox(*), irowbox(*), nboxes, nblevel(0:1), iboxlev(*)
      integer istartlev(0:1), ifoct
      integer ier, iprec, ns, isrcladder(2,*), ifcharge,ifdipole,ifquad
      integer isave(*), ifpot, ifgrad, ifhess, ifder3
      real *8 zll(2), blength, srcsort(2,*), dsave(*), beta
      real *8 pot(*), grad(2,*), hess(3,*), der3(4,*)
      real *8 chargesort(*), dipstrsort(*), dipvecsort(2,*)
      real *8 quadstrsort(*), quadvecsort(3,*)
      real *8 octstrsort(*), octvecsort(4,*)
      complex *16 csave(*)
c     local
      integer nterms, nnodes , msrcflag, mflagghost
      integer mneighbors, mnnbrs, mrscales
      integer mlloc, mmbhloc, mbiglloc, mbigmbhloc

c     grab pointers

      nterms = isave(1)
      nnodes = isave(2)

      msrcflag = isave(11)
      mflagghost = isave(12)
      mneighbors = isave(16)
      mnnbrs = isave(17)

      mrscales = isave(21)

      mlloc = isave(31)
      mmbhloc = isave(32)
      mbiglloc = isave(33)
      mbigmbhloc = isave(34)


c     call main routine

      call mbhfmm2d3_srcsrc1(beta,nlev,levelbox,iparentbox,
     1     ichildbox,icolbox,irowbox,nboxes,nblevel,
     2     iboxlev,istartlev,isave(mneighbors),isave(mnnbrs),nterms,
     3     csave(mlloc),csave(mmbhloc),
     3     csave(mbigmbhloc),csave(mbiglloc),isave(mflagghost),
     4     dsave(mrscales),isave(msrcflag),ns,srcsort,isrcladder,
     8     ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     9     ifquad,quadstrsort,quadvecsort,ifoct,octstrsort,
     1     octvecsort,zll,blength,
     7     ifpot,pot,ifgrad,grad,ifhess,hess,ifder3,der3)

      return
      end



      subroutine mbhfmm2d3_srcsrc1(beta,nlev,levelbox,iparentbox,
     1     ichildbox,icolbox,irowbox,nboxes,nblevel,
     2     iboxlev,istartlev,neighbors,nnbrs,nterms,lloc,mbhloc,
     3     bigmbhloc,biglloc,iflagghost,rscales,isrcflag,ns,srcsort,
     4     isrcladder,ifcharge,chargesort,ifdipole,dipstrsort,
     5     dipvecsort,ifquad,quadstrsort,quadvecsort,ifoct,octstrsort,
     5     octvecsort,zll,blength,
     6     ifpot,pot,ifgrad,grad,ifhess,hess,ifder3,der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global 
      integer nlev, levelbox(*), iparentbox(*), ichildbox(4,*)
      integer icolbox(*), irowbox(*), nboxes, nblevel(0:1), iboxlev(*)
      integer istartlev(0:1), isrcladder(2,*), nterms, nnodes
      integer ier, isrcflag(*), iflagghost(*)
      integer neighbors(12,*), nnbrs(*)
      integer ns, ifpot, ifgrad, ifhess, ifcharge, ifdipole, ifquad
      integer ifder3, ifoct
      real *8 srcsort(2,*), zll(2), blength, rscales(0:nlev)
      real *8 beta
      complex *16 mbhloc(0:nterms,*), lloc(0:nterms,*)
      complex *16 bigmbhloc(0:nterms,4,*), biglloc(0:nterms,4,*)
      real *8 chargesort(*), dipstrsort(*), dipvecsort(2,*)
      real *8 quadstrsort(*), quadvecsort(3,*)
      real *8 octstrsort(*), octvecsort(4,*)
      real *8 pot(*), grad(2,*), hess(3,*), der3(4,*)
c     local
      integer i, ibox, ii, icol, irow, jsrc, nsbox, ntbox, istartbox, ic
      integer iendbox, ier2, ilev, OMP_MIN_VAL, jstartbox, jendbox
      integer nsboxj, iii
      parameter (OMP_MIN_VAL = 50)
      real *8 xlength, zbox(2), zboxc(2), rscale, xt, yt
      real *8 dzero
      data dzero / 0.0d0 /
      integer nsrcboxes
      integer, allocatable :: isrcboxes(:)


c     find all boxes with sources

      allocate(isrcboxes(nboxes))

      nsrcboxes = 0
      do i = 1,nboxes
         if (isrcflag(i) .eq. 1) then
            nsrcboxes = nsrcboxes+1
            isrcboxes(nsrcboxes) = i
         endif
      enddo

C$OMP PARALLEL DO PRIVATE(ibox,xlength,ii,jsrc,istartbox,iendbox,nsbox,
C$OMP& icol, irow, ntbox, rscale, zbox, zboxc, ier2, ic, xt, yt, ilev,
C$OMP& jstartbox,jendbox,nsboxj,iii    ) 
C$OMP& IF(nsrcboxes .gt. OMP_MIN_VAL)
C$OMP& SCHEDULE(dynamic,10)
      do i = 1,nsrcboxes

         ibox = isrcboxes(i)
         
         istartbox = isrcladder(1,ibox)
         iendbox = isrcladder(2,ibox)
         nsbox = iendbox-istartbox+1

c     initialize

         if (ifpot .eq. 1) then
            do ii = istartbox,iendbox
               pot(ii) = dzero
            enddo
         endif
         if (ifgrad .eq. 1) then
            do ii = istartbox,iendbox
               grad(1,ii) = dzero
               grad(2,ii) = dzero
            enddo
         endif
         if (ifhess .eq. 1) then
            do ii = istartbox,iendbox
               hess(1,ii) = dzero
               hess(2,ii) = dzero
               hess(3,ii) = dzero
            enddo
         endif
         if (ifder3 .eq. 1) then
            do ii = istartbox,iendbox
               der3(1,ii) = dzero
               der3(2,ii) = dzero
               der3(3,ii) = dzero
               der3(4,ii) = dzero
            enddo
         endif

c     process source box

         ilev = levelbox(ibox)
         xlength = blength/(2.0d0**ilev)

c     self direct 

         call mbhpotgrad2dall_cdqo3_self_add(beta,srcsort(1,istartbox),
     1        nsbox,ifcharge,chargesort(istartbox),ifdipole,
     2        dipstrsort(istartbox),dipvecsort(1,istartbox),
     3        ifquad,quadstrsort(istartbox),quadvecsort(1,istartbox),
     3        ifoct,octstrsort(istartbox),octvecsort(1,istartbox),
     4        ifpot,pot(istartbox),ifgrad,grad(1,istartbox),
     5        ifhess,hess(1,istartbox),ifder3,der3(1,istartbox))

c     neighbors: direct

         do ii = 1,nnbrs(ibox)
            jsrc = neighbors(ii,ibox)
            if (isrcflag(jsrc) .eq. 1) then

               jstartbox = isrcladder(1,jsrc)
               jendbox = isrcladder(2,jsrc)
               nsboxj = jendbox-jstartbox+1

               do iii = istartbox,iendbox
                  
                  call mbhpotgrad2dall_cdqo3_add(beta,
     1                 srcsort(1,jstartbox),
     1                 nsboxj,ifcharge,chargesort(jstartbox),ifdipole,
     2                 dipstrsort(jstartbox),dipvecsort(1,jstartbox),
     3                 ifquad,quadstrsort(jstartbox),
     4                 quadvecsort(1,jstartbox),ifoct,
     4                 octstrsort(jstartbox),octvecsort(1,jstartbox),
     6                 srcsort(1,iii),ifpot,pot(iii),
     5                 ifgrad,grad(1,iii),ifhess,hess(1,iii),ifder3,
     8                 der3(1,iii))
               enddo
               
            endif
         enddo
         
c     contributions from small to big far

         if (iflagghost(ibox) .eq. 1) then

            irow = irowbox(ibox)
            icol = icolbox(ibox)
            zbox(1) = zll(1)+xlength*(icol-0.5d0)
            zbox(2) = zll(2)+xlength*(irow-0.5d0)
            
            do ii = istartbox,iendbox
            
               xt = srcsort(1,ii)
               yt = srcsort(2,ii)

               ntbox = 1

c     figure out which ghost child is relevant
               
               if (xt .lt. zbox(1)) then
                  if (yt .gt. zbox(2)) then
                     ic = 1
                     zboxc(1) = zbox(1)-xlength/4.0d0
                     zboxc(2) = zbox(2)+xlength/4.0d0
                  else
                     ic = 4
                     zboxc(1) = zbox(1)-xlength/4.0d0
                     zboxc(2) = zbox(2)-xlength/4.0d0
                  endif
               else
                  if (yt .gt. zbox(2)) then
                     ic = 2
                     zboxc(1) = zbox(1)+xlength/4.0d0
                     zboxc(2) = zbox(2)+xlength/4.0d0
                  else
                     ic = 3
                     zboxc(1) = zbox(1)+xlength/4.0d0
                     zboxc(2) = zbox(2)-xlength/4.0d0
                  endif
               endif

               call mbh2dtaeval3all(beta,rscales(ilev+1),zboxc,
     1              bigmbhloc(0,ic,ibox),biglloc(0,ic,ibox),
     1              nterms,srcsort(1,ii),ntbox,ifpot,pot(ii),
     2              ifgrad,grad(1,ii),
     3              ifhess,hess(1,ii),ifder3,der3(1,ii))
            enddo


            
         endif

c     evaluate local expansions

         irow = irowbox(ibox)
         icol = icolbox(ibox)
         zbox(1) = zll(1)+xlength*(icol-0.5d0)
         zbox(2) = zll(2)+xlength*(irow-0.5d0)

         call mbh2dtaeval3all(beta,rscales(ilev),zbox,mbhloc(0,ibox),
     1        lloc(0,ibox),nterms,srcsort(1,istartbox),nsbox,
     2        ifpot,pot(istartbox),
     2        ifgrad,grad(1,istartbox),ifhess,hess(1,istartbox),
     3        ifder3,der3(1,istartbox))

      enddo
C$OMP END PARALLEL DO

      return
      end
