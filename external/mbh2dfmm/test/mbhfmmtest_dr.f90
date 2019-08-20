
program mbhfmmtest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine tests the fmm against direct calculation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit real *8 (a-h,o-z)
  
  integer maxboxes, maxlevel, nboxes, nlev, maxsrc
  parameter (maxboxes = 30000, maxlevel = 30)
  integer, dimension(maxboxes) :: levelbox,icolbox,irowbox,iparentbox
  integer, dimension(maxboxes) :: itemparray, iflag
  integer ichildbox(4,maxboxes), istartlev(0:maxlevel)
  integer iboxlev(maxboxes),icolleagbox(9,maxboxes),nblevel(0:maxlevel)
  integer neighbors(12,maxboxes), nnbrs(maxboxes)
  integer iboxes(maxboxes), localonoff(maxboxes)
  integer, allocatable :: isave(:)
  real *8, allocatable :: dsave(:)
  complex *16, allocatable :: csave(:)
  parameter (maxsrc = 1000000)
  integer :: isrcsort(maxsrc), isrcladder(2,maxboxes)
  integer :: itargsort(maxsrc), itargladder(2,maxboxes)  
  real *8 :: src(2,maxsrc), srcsort(2,maxsrc), zll(2), blength
  real *8 :: targ(2,maxsrc), targsort(2,maxsrc)
  real *8 :: dsdt(maxsrc), rnx(maxsrc), rny(maxsrc)
  real *8 :: charge(maxsrc), dipstr(maxsrc), chargesort(maxsrc)
  real *8 :: dipstrsort(maxsrc), dipvec(2,maxsrc), dipvecsort(2,maxsrc)
  real *8 :: quadstr(maxsrc), quadstrsort(maxsrc), quadvec(3,maxsrc)
  real *8 :: quadvecsort(3,maxsrc)
  real *8 :: octstr(maxsrc), octstrsort(maxsrc), octvec(4,maxsrc)
  real *8 :: octvecsort(4,maxsrc)
  real *8 :: pot(maxsrc), grad(2,maxsrc), hess(3,maxsrc)
  real *8 :: pottarg(maxsrc), gradtarg(2,maxsrc), hesstarg(3,maxsrc)
  real *8 gradtemp(2), gradtargtemp(2), err2(2)
  real *8 hesstemp(3), hesstargtemp(3), err3(3)
  parameter (iseed = 281+3308004)

  beta = 1.0d-15

  temp = hkrand(iseed)
  call prini(6,13)

  ! Test accuracy against direct calculation

  ifdirect = 1

  ! Print out geometry, targets, and leaf boxes to various files

  ifprint = 0

  ! Type of calculation 

  ifcharge = 0
  ifdipole = 0
  ifquad = 1
  ifoct = 0
  ifpot = 1
  ifgrad = 1
  ifhess = 1

  ! Geometry  

  pi = 4.0d0*atan(1.0d0)


  ns = 10000

  if (1 .eq. 0) then

     ! ellipse axes
     a = 1.123d0
     b = 0.7d0

     ! ellipse center
     xc = 6.1d0
     yc = -1.3d0

     ! box length and lower left corner of computational
     ! domain (a square) should contain all src and targ

     blength = 3.0d0*max(a,b)

     zll(1) = xc-blength/2.0d0
     zll(2) = yc-blength/2.0d0
     
     ! sources
     h = 2.0d0*pi/ns

     do i = 1,ns

        t = (i-1)*h
        src(1,i) = xc+a*cos(t)
        src(2,i) = yc+b*sin(t)-.01d0
        dx = -a*sin(t)
        dy = b*cos(t)
        dsdt(i) = sqrt(dx**2+dy**2)
        rnx(i) = dy/dsdt(i)
        rny(i) = -dx/dsdt(i)
        if (ifprint .eq. 1) write(21,*) src(1,i), src(2,i), &
             rnx(i), rny(i)
     enddo

     !targets
     nt = 2*(1000)
     ht = 2.0d0*pi/(nt/2)
     curvenrm = h*dsqrt(a**2+b**2)
     do i = 1,nt
        t = (i-1)*ht
        targ(1,i) = xc+ (a+10.0d0*curvenrm)*cos(t)
        targ(2,i) = yc+ (b+10.0d0*curvenrm)*sin(t) - .01d0
        targ(1,i+nt/2) = xc+ (a-10.0d0*curvenrm)*cos(t)
        targ(2,i+nt/2) = yc+ (b-10.0d0*curvenrm)*sin(t) - .01d0
     enddo
     
  else

     ! box length and lower left corner of computational
     ! domain (a square) should contain all src and targ
     blength = 1.0d0
     zll(1:2) = (/ -0.5d0, -0.5d0 /)

     ! sources
     do i = 1,ns
        src(1,i) = zll(1) + 0*blength/4.0d0 + blength*hkrand(0)/4.0d0
        src(2,i) = zll(2) + 0*blength/4.0d0 + blength*hkrand(0)/4.0d0
        src(1,i) = zll(1) + blength*hkrand(0)
        src(2,i) = zll(2) + blength*hkrand(0)
        dx = -0.5d0+hkrand(0)
        dy = -0.5d0+hkrand(0)
        dsdt(i) = sqrt(dx**2+dy**2)
        rnx(i) = dy/dsdt(i)
        rny(i) = -dx/dsdt(i)
     enddo

     ! targets
     nt = 1000
     do i = 1,nt
        targ(1,i) = zll(1) + 0*blength/4.0d0 + blength*hkrand(0)/4.0d0
        targ(2,i) = zll(2) + 2*blength/4.0d0 + blength*hkrand(0)/4.0d0
        targ(1,i) = zll(1) + blength*hkrand(0)
        targ(2,i) = zll(2) + blength*hkrand(0)
     enddo
     
  endif

  do i = 1,nt
     if (ifprint .eq. 1) write(22,*) targ(1,i), targ(2,i)
  enddo

  ! set up charges and dipoles 

  do i = 1,ns
     t = (i-1)*h
     charge(i) = beta**2*dsdt(i)*h*cos(3*t)
     dipstr(i) = beta*dsdt(i)*h*sin(4*t)
     quadstr(i) = dsdt(i)*h*cos(4*t)
     dipvec(1:2,i) = (/ rnx(i), rny(i) /)
     ! d_n d_tau
     quadvec(1:3,i) = (/ -rnx(i)*rny(i), rnx(i)*rnx(i)-rny(i)*rny(i), &
          rny(i)*rnx(i) /)
  enddo

  !call prin2('quadstr *',quadstr,3*ns)
  !call prin2('quadvec *',quadvec,3*ns)

  maxnodes = 10
  time1 = second()
!  call lrt2d_mktpts(levelbox, icolbox, irowbox, nboxes, nlev, &
!       iparentbox, ichildbox, nblevel, iboxlev, istartlev, &
!       maxboxes, itemparray, maxlevel, src, srcsort, isrcsort, &
  !       isrcladder, ns, maxnodes, zll, blength, ier)
  call lrt2d_mktst(levelbox, icolbox, irowbox, nboxes, nlev, &
       iparentbox, ichildbox, nblevel, iboxlev, istartlev, &
       maxboxes, itemparray, maxlevel, src, srcsort, isrcsort, &
       isrcladder, ns, targ, targsort, itargsort, itargladder, nt, &
       maxnodes, zll, blength, ier, localonoff)
  
  time2 = second()

  call prin2('TIME TO MAKE TREE *',time2-time1,1)

  ifixflag = 0
  call lrt2d_restrict(levelbox,iparentbox,ichildbox,icolbox, &
       irowbox,icolleagbox,nboxes,nlev, &
       nblevel,iboxlev,istartlev,ifixflag)

  write(*,*) 'nboxes   ', nboxes
  write(*,*) 'nlev     ', nlev
  
  if (ifixflag .eq. 1) then
     write(*,*) 'fixing ... '
     time1 = second()
     call lrt2d_fix(levelbox,iparentbox,ichildbox,icolbox, &
          irowbox,icolleagbox,nboxes,nlev, &
          nblevel,iboxlev,istartlev, &
          iflag, maxboxes,itemparray)
     time2 = second()
     call prin2('TIME FOR RESTRICTION *',time2-time1,1)
     write(*,*) 'nboxes ', nboxes
  endif
  
  call lrt2d_testtree(levelbox,iparentbox,ichildbox, &
       icolbox,irowbox,nboxes,nlev,nblevel,iboxlev, &
       istartlev)

  call lrt2d_mkcolls(icolbox,irowbox,icolleagbox,nboxes,nlev, &
       iparentbox,ichildbox,nblevel,iboxlev,istartlev)

  time1 = second()
  call lrt2d_mknbrs(neighbors,nnbrs,nboxes,ichildbox, &
       iparentbox,icolleagbox,icolbox,irowbox)
  time2 = second()

  call prin2('TIME TO MAKE NEIGHBORS *',time2-time1,1)

!  nlev = 2
!  call lrt2d_uni(levelbox,icolbox,irowbox,nboxes,nlev, &
!       ichildbox, iparentbox, nblevel, istartlev, iboxlev)

!  imode = 4
!  call lrt2d_set(levelbox,icolbox,irowbox,nboxes,nlev, &
!      ichildbox, iparentbox, nblevel, istartlev, iboxlev,imode)

  call lrt2d_testtree(levelbox,iparentbox,ichildbox, &
       icolbox,irowbox,nboxes,nlev,nblevel,iboxlev, &
       istartlev)

  write(*,*) 'nboxes ', nboxes
  write(*,*) 'nlev ', nlev

  do i = 1,nboxes
     if (ichildbox(1,i) .lt. 0) then
        xlength = blength/(2.0d0**levelbox(i))
        irow = irowbox(i)
        icol = icolbox(i)
        x = zll(1)+xlength*(icol-1)
        y = zll(2)+xlength*(irow-1)
        if (ifprint .eq. 1) write(23,*) x, y, xlength, xlength
     endif
  enddo
  
  call lrt2d_ptsort(levelbox,icolbox,irowbox,nboxes,nlev, &
       iparentbox, ichildbox, nblevel, iboxlev, istartlev, &
       src, srcsort, isrcsort, &
       isrcladder, ns, zll, blength, ier)

  call lrt2d_ptsort_wc(levelbox,icolbox,irowbox,nboxes,nlev, &
       iparentbox, ichildbox, nblevel, iboxlev, istartlev, &
       targ, targsort, itargsort, &
       itargladder, ns, zll, blength, ier)

  do i = 1,nboxes
     localonoff(i) = 0
     if (itargladder(2,i)-itargladder(1,i)+1 .gt. 0) localonoff(i) = 1
  enddo

  do i = 1,ns
     chargesort(i) = charge(isrcsort(i))
     dipstrsort(i) = dipstr(isrcsort(i))
     dipvecsort(1:2,i) = dipvec(1:2,isrcsort(i))
     quadstrsort(i) = quadstr(isrcsort(i))
     quadvecsort(1:3,i) = quadvec(1:3,isrcsort(i))
  enddo

  ! find field at targets directly 

  call prinf('START DIRECT ......*',ifpot,0)


  time1 = second()
  if (ifdirect .eq. 1) then
     do i = 1,nt
        call mbhpotgrad2dall_cdq(beta,src,ns,ifcharge, &
             charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec, &
             targ(1,i),ifpot,pot(i),ifgrad,grad(1,i),ifhess,hess(1,i))
     enddo
  endif
  time2 = second()

  call prinf('END DIRECT ........*',ifpot,0)

  call prin2('TIME FOR DIRECT*',time2-time1,1)  

  iprec = 4

  call prinf('START FMM .........*',ifpot,0)
  
  time1 = second()
  lisave = -1
  ldsave = -1
  lcsave = -1
  
  write(*,*) 'query storage '

  ifalltarg = 0
  call mbhfmm2d_form(beta,ier,iprec,nlev,levelbox,iparentbox, &
       ichildbox, icolbox, irowbox, nboxes, nblevel, iboxlev, &
       istartlev, ifalltarg, localonoff, zll, blength, ns, srcsort, &
       isrcladder, ifcharge, chargesort, ifdipole, dipstrsort, &
       dipvecsort, ifquad, quadstrsort, quadvecsort, ifoct, &
       octstrsort, octvecsort, isave, lisave, dsave, ldsave, &
       csave, lcsave)

  write(*,*) lisave, ldsave, lcsave

  allocate(isave(lisave),dsave(ldsave),csave(lcsave))

  time1 = second()

  !$ time1 = omp_get_wtime()

  call mbhfmm2d_form(beta,ier,iprec,nlev,levelbox,iparentbox, &
       ichildbox, icolbox, irowbox, nboxes, nblevel, iboxlev, &
       istartlev, ifalltarg, localonoff, zll, blength, ns, srcsort, &
       isrcladder, ifcharge, chargesort, ifdipole, dipstrsort, &
       dipvecsort, ifquad, quadstrsort, quadvecsort, ifoct, &
       octstrsort, octvecsort, isave, lisave, dsave, ldsave, &
       csave, lcsave)

  time2 = second()

  !$ time2 = omp_get_wtime()

  call prin2('TIME FOR FORMING EVERYTHING *',time2-time1,1)

  time1 = second()

  !$ time1 = omp_get_wtime()

  write(*,*) 'beta ', beta

  call mbhfmm2d_targ(beta, ier, nlev, levelbox, iparentbox, &
       ichildbox, icolbox, irowbox, nboxes, nblevel, iboxlev, &
       istartlev, zll, blength, ns, srcsort, isrcladder, &
       ifcharge, chargesort, ifdipole, dipstrsort, dipvecsort, &
       ifquad, quadstrsort, quadvecsort, isave, dsave, csave, &
       nt, targ, ifpot, pottarg, ifgrad, gradtarg, ifhess, hesstarg)

  time2 = second()

  !$ time2 = omp_get_wtime()

  write(*,*) 'beta ', beta

  call prin2('TIME FOR EVALS *',time2-time1,1)
  call prin2('EVALS PER SECOND *',nt/(time2-time1),1)

  call prinf('END FMM ...........*',ifpot,0)



  if (ifdirect .eq. 1) then
     
     if (ifpot .eq. 1) then
        
        absnormmax = 0.0d0
        abserrmax = 0.0d0

        do i = 1,nt
           temp1 = pot(i)
           temp2 = pottarg(i)
           !write(*,*) i, temp1/temp2
           absnorm = dabs(temp1)
           abserr = dabs(temp1-temp2)
           if (absnorm .gt. absnormmax) then
              absnormmax = absnorm
              !write(*,*) temp1, temp2
           endif
           if (abserr .gt. abserrmax) then
              abserrmax = abserr
              write(*,*) i, temp1, temp2
           endif
        enddo

        call prin2('ABS ERR (POT) *',abserrmax,1)
        call prin2('REL ERR (POT) *',abserrmax/absnormmax,1)
     endif

     if (ifgrad .eq. 1) then
        
        absnormmax = 0.0d0
        abserrmax = 0.0d0

        do i = 1,nt
           gradtemp = grad(1:2,i)
           gradtargtemp = gradtarg(1:2,i)
           err2 = gradtemp-gradtargtemp
           absnorm = dnorm_dr(gradtemp,2)
           abserr = dnorm_dr(err2,2)
           if (absnorm .gt. absnormmax) absnormmax = absnorm
           if (abserr .gt. abserrmax) abserrmax = abserr
        enddo

        call prin2('ABS ERR (GRAD) *',abserrmax,1)
        call prin2('REL ERR (GRAD) *',abserrmax/absnormmax,1)
     endif

     if (ifhess .eq. 1) then
        absnormmax = 0.0d0
        abserrmax = 0.0d0

        do i = 1,nt
           hesstemp = hess(1:3,i)
           hesstargtemp = hesstarg(1:3,i)
           err3 = hesstemp-hesstargtemp
           absnorm = dnorm_dr(hesstemp,3)
           abserr = dnorm_dr(err3,3)
           if (absnorm .gt. absnormmax) absnormmax = absnorm
           if (abserr .gt. abserrmax) abserrmax = abserr
        enddo

        call prin2('ABS ERR (HESS)*',abserrmax,1)
        call prin2('REL ERR (HESS)*',abserrmax/absnormmax,1)
     endif
  endif

  stop
end program mbhfmmtest

real *8 function dnorm_dr(v,lv)
  implicit none
  ! global
  real *8 v(lv)
  integer lv
  ! local
  integer i

  dnorm_dr = 0.0d0

  do i = 1,lv
     dnorm_dr = dnorm_dr + v(i)**2
  enddo

  dnorm_dr = sqrt(dnorm_dr)

  return
end function dnorm_dr
