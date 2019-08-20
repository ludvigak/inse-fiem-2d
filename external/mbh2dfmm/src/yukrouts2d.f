cc Copyright (C) 2017: Travis Askham, Leslie Greengard, Zydrunas Gimbutas
cc email: askhamwhat@gmail.com      
cc 
cc This software is being released under a modified FreeBSD license
cc (see licenses folder in home directory). 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date: 2011-04-25 22:45:52 -0400 (Mon, 25 Apr 2011) $
c    $Revision: 1889 $
c
c
c     Changed to Yukawa potential subroutines Sep. 2014, T. Askham
c
c      This file contains the basic subroutines for 
c      forming and evaluating multipole (partial wave) expansions
c      in two dimensions.
c
c      Remarks on scaling conventions.
c
c      1)  Hankel and Bessel functions are consistently scaled as
c       	kvec(n)= (2/pi) K_n(z)*rscale^(n)
c       	ivec(n)= I_n(z)/rscale^(n)
c
c          rscale should be of the order of |z| if |z| < 1. Otherwise,
c          rscale should be set to 1.
c
c      2) The potential scaling is that required of the delta function
c      response:  pot = H_0(k*r)*(eye/4). (This is only seen when 
c      evaluating)
c
c-----------------------------------------------------------------------
c
c      Y2DMPEVAL: computes potential and grad(potential)
c                 due to a multipole expansion.
c
c      Y2DFORMMP: creates multipole expansion (outgoing) due to 
c                 a collection of charge sources.
c
c      Y2DTAEVAL: computes potential and grad(potential) 
c                  due to local expansion.
c
c      Y2DFORMTA: creates local expansion due to 
c                 a collection of charge sources.
c
c      Y2DMPEVALALL: computes potential and grad(potential)
c                 due to a multipole expansion for a collection of targets
c
c      Y2DTAEVALALL: computes potential and grad(potential) 
c                  due to local expansion for a collection of targets
c
c      Y2DMPMP:     Converts multipole expansion to a multipole expansion.
c      Y2DMPLOC:     Converts multipole expansion to a local expansion.
c      Y2DLOCLOC:     Converts local expansion to a local expansion.
c
c      YPOTGRAD2DALL:  direct calculation for a collection of charge sources
c      YPOTGRAD2D : direct calculation for a single charge source
c
c
c
c      Y2DFORMMP_DP: creates multipole expansion (outgoing) due to 
c                 a collection of dipole sources.
c
c      Y2DFORMTA_DP: creates local expansion due to 
c                 a collection of dipole sources.
c
c      YPOTGRAD2DALL_DP:  direct calculation for a collection of dipoles
c      YPOTGRAD2D_DP : direct calculation for a single dipole
c
c
c      YPOTGRAD2DALL_SDP:  direct calculation for 
c                 a collection of charges and dipoles
c      YPOTGRAD2D_SDP : direct calculation for a single charge and a dipole
c
c-----------------------------------------------------------------------
c
c
c
c
c**********************************************************************
      subroutine y2dmpeval(dk,rscale,center,mpole,nterms,ztarg,
     1      pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an outgoing partial wave expansion.
c
c     pot  = (1/(2pi)) K_0(kr) mpole_0 + 
c                  
c                  nterms
c       (1/(2pi))   sum   K_n(k r) rscale^n [real(mpole_n) cos(n theta) 
c                   n=1                    + imag(mpole_n) sin(n theta)]
c
c                                   
c     grad  = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     dk     :    Yukawa parameter (real)
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 zk,mpole(0:nterms),zmul,ztemp1
        real *8 pot, grad(2), hess(3)
        real *8 center(2),ztarg(2),zdiff(2)
        real *8 dk, pi
c
        complex *16, allocatable :: hval(:)
        complex *16, allocatable :: hder(:)

        real *8, allocatable :: kval(:)
c
        complex *16, allocatable :: mpolex(:)
        complex *16, allocatable :: mpoley(:)
c
        complex *16, allocatable :: mpolexx(:)
        complex *16, allocatable :: mpolexy(:)
        complex *16, allocatable :: mpoleyy(:)
c
        complex *16, allocatable :: mptemp(:)
c
        complex *16 ima,ima4,z,pieye2
        data ima/(0.0d0,1.0d0)/
c
c
        ima4=-4*ima
        pi = 4.0d0*datan(1.0d0)
        pieye2=pi*ima/2.0d0
c
        allocate(hval(0:nterms+10), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(kval(0:nterms+10), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(hder(0:nterms+10), stat=ier)
        if (ier.eq.1) then
        return
        endif
c
        zdiff(1)=ztarg(1)-center(1)
        zdiff(2)=ztarg(2)-center(2)
        call y2cart2polar(zdiff,r,theta)
c
        zk = ima*dk
        z=zk*r
        ifder=0
        call y2dall(nterms+3,z,rscale,hval,ifder,hder)
c
        allocate(mptemp(0:nterms+2), stat=ier)
        if (ier.eq.1) then
          return
        endif
c

c     convert Hankel H_n to modified Bessel K_n
c     (up to scale of pi/2)
        do n = 0,nterms+2,4
           kval(n) = -dimag(hval(n))
           kval(n+1) = -dreal(hval(n+1))
           kval(n+2) = dimag(hval(n+2))
           kval(n+3) = dreal(hval(n+3))
        enddo

c     at this point kval(n) = 2/pi * K_n(dk*r)
           
        mptemp(0)=kval(0)
        zmul=exp(ima*theta)
        ztemp1=zmul
        do n=1,nterms+2
           mptemp(n)=dcmplx(kval(n)*dreal(ztemp1),kval(n)*dimag(ztemp1))
           ztemp1 = ztemp1*zmul
        enddo

c
c
        if( ifgrad .eq. 1 .or. ifhess .eq. 1 ) then
c
           allocate(mpolex(0:nterms+1), stat=ier)
           if (ier.eq.1) then
           return
           endif
           allocate(mpoley(0:nterms+1), stat=ier)
           if (ier.eq.1) then
           return
           endif
c
           do i=0,nterms+1
              mpolex(i)=0
              mpoley(i)=0
           enddo

           mpolex(1) = -dk/rscale*mpole(0)
           mpoley(1) = -dk/rscale*mpole(0)*ima

           do i=1,nterms
              mpolex(i-1) = mpolex(i-1) -dk/2.0d0*mpole(i)*rscale
              mpolex(i+1) = mpolex(i+1) -dk/2.0d0*mpole(i)/rscale
              mpoley(i-1) = mpoley(i-1) +dk/2.0d0*mpole(i)*rscale*ima
              mpoley(i+1) = mpoley(i+1) -dk/2.0d0*mpole(i)/rscale*ima
           enddo
        endif
c
        if( ifhess .eq. 1 ) then
           allocate(mpolexx(0:nterms+2), stat=ier)
           if (ier.eq.1) then
           return
           endif
           allocate(mpolexy(0:nterms+2), stat=ier)
           if (ier.eq.1) then
           return
           endif
           allocate(mpoleyy(0:nterms+2), stat=ier)
           if (ier.eq.1) then
           return
           endif
c
           do i=0,nterms+2
              mpolexx(i)=0
              mpolexy(i)=0
              mpoleyy(i)=0
           enddo

           mpolexx(1) = -dk/1.0d0/rscale*dreal(mpolex(0))
           mpolexy(1) = -dk/1.0d0/rscale*dreal(mpolex(0))*ima
           mpoleyy(1) = -dk/1.0d0/rscale*dreal(mpoley(0))*ima

           do i=1,nterms+1
              mpolexx(i-1) = mpolexx(i-1) -dk/2.0d0*mpolex(i)*rscale
              mpolexx(i+1) = mpolexx(i+1) -dk/2.0d0*mpolex(i)/rscale
              mpolexy(i-1) = mpolexy(i-1) +dk/2.0d0*mpolex(i)*rscale*ima
              mpolexy(i+1) = mpolexy(i+1) -dk/2.0d0*mpolex(i)/rscale*ima
              mpoleyy(i-1) = mpoleyy(i-1) +dk/2.0d0*mpoley(i)*rscale*ima
              mpoleyy(i+1) = mpoleyy(i+1) -dk/2.0d0*mpoley(i)/rscale*ima
           enddo


        endif
c
c
        pot=kval(0)*dreal(mpole(0))
        do n=1,nterms
           pot=pot+dreal(mpole(n))*dreal(mptemp(n))
     1          +dimag(mpole(n))*dimag(mptemp(n))
        enddo
        pot=pot*0.25d0
c

        if( ifgrad .eq. 1 ) then
           grad(1)=0
           grad(2)=0
c
           grad(1)=kval(0)*dreal(mpolex(0))
           grad(2)=kval(0)*dreal(mpoley(0))

           do n=1,nterms+1
              grad(1)=grad(1)+dreal(mpolex(n))*dreal(mptemp(n))+
     1             dimag(mpolex(n))*dimag(mptemp(n))
              grad(2)=grad(2)+dreal(mpoley(n))*dreal(mptemp(n))+
     1             dimag(mpoley(n))*dimag(mptemp(n))
           enddo
           grad(1)=grad(1)*.25d0
           grad(2)=grad(2)*.25d0

        endif
c
c
        if( ifhess .eq. 1 ) then
           hess(1)=0
           hess(2)=0
           hess(3)=0
c
           hess(1)=kval(0) *dreal(mpolexx(0))
           hess(2)=kval(0) *dreal(mpolexy(0))
           hess(3)=kval(0) *dreal(mpoleyy(0))
           do n=1,nterms+2
              hess(1)=hess(1)+
     1             dreal(mpolexx(n))*dreal(mptemp(n))
     2             +dimag(mpolexx(n))*dimag(mptemp(n))
              hess(2)=hess(2)+
     1             dreal(mpolexy(n))*dreal(mptemp(n))
     2             +dimag(mpolexy(n))*dimag(mptemp(n))
              hess(3)=hess(3)+
     1             dreal(mpoleyy(n))*dreal(mptemp(n))
     2             +dimag(mpoleyy(n))*dimag(mptemp(n))
           enddo
           hess(1)=hess(1)*.25d0
           hess(2)=hess(2)*.25d0
           hess(3)=hess(3)*.25d0
        endif
c
        return
        end
c
c
c
c
c
c**********************************************************************
      subroutine y2dtaeval(dk,rscale,center,mpole,nterms,ztarg,
     1      pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an incoming partial wave expansion.
c               +nterms
c     pot  =      sum   I_n(k r)/rscale^n [real(mpole_n)*cos(n*theta)+
c                 n=0                       imag(mpole_n)*sin(n*theta)]
c     grad = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c     (r,theta) are the polar coordinates of ztarg-center
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk     :    Helmholtz parameter
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     ifgrad  :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        real *8 pot,grad(2),hess(3),dk
        complex *16 zk, mpole(0:nterms)
        real *8 center(2),ztarg(2),zdiff(2)
        integer itt,iisc,lival
        parameter (itt=30000)
        parameter (iisc=1000)
        parameter (lival=1000)
        complex *16, allocatable :: cw1(:)
        integer, allocatable :: iscale1(:)
        real *8, allocatable :: ival1(:)
        complex *16 ima,ima4,ima4inv,z,zmull,zmullinv,zfac
        complex *16 cw(0:itt)
        real *8 ival(0:lival)
        integer iscale(0:iisc)
        data ima/(0.0d0,1.0d0)/
c
c
        zk = ima*dk

        ima4=-4*ima
        ima4inv=ima/4
c
        lwfjs = nterms+5 + 4*nterms + 100
        ijval = 0
        ijder = ijval + lwfjs+4
        imptemp = ijder + lwfjs+4
        impolex = imptemp + nterms+8
        impoley = impolex + nterms+8
        impolexx = impoley + nterms+8
        impolexy = impolexx + nterms+8
        impoleyy = impolexy + nterms+8
        itot = impoleyy + nterms+8

        lival1 = nterms+3
c
        ialloc = 0
        if ((itot.gt.itt).or.(lwfjs+10 .gt. iisc)
     1       .or.(lival1.gt.lival)) then
           allocate(cw1(0:itot))
           allocate(iscale1(0:lwfjs+10))
           allocate(ival1(0:lival1))
           ialloc = 1
        endif

c        write(*,*) ' ialloc is ',ialloc
c
        zdiff(1)=ztarg(1)-center(1)
        zdiff(2)=ztarg(2)-center(2)
        call y2cart2polar(zdiff,r,theta)
c
        z=zk*r
        ifder=0
        if (ialloc.eq.0) then 
           call jfuns2d(ier,nterms+3,z,rscale,cw(ijval),ifder,cw(ijder),
     1        lwfjs,iscale,ntop)
           call mkmptempyk(cw(imptemp),theta,cw(ijval),ival,nterms)
           if( ifgrad .eq. 1 .or. ifhess .eq. 1 ) then
              call mkmpole12yk(cw(impolex),cw(impoley),
     1             ifhess,cw(impolexx),cw(impolexy),
     2             cw(impoleyy),mpole,dk,rscale,nterms)
           endif
           call taevalsyk(mpole,cw(impolex),cw(impoley),cw(impolexx),
     1             cw(impolexy),cw(impoleyy),rscale,nterms,
     2             cw(ijval),cw(imptemp),pot,ifgrad,grad,
     3             ifhess,hess,ima4inv)
        else
           call jfuns2d(ier,nterms+3,z,rscale,cw1(ijval),ifder,
     1        cw1(ijder),lwfjs,iscale,ntop)
           call mkmptempyk(cw1(imptemp),theta,cw1(ijval),ival1,nterms)
           if( ifgrad .eq. 1 .or. ifhess .eq. 1 ) then
              call mkmpole12yk(cw1(impolex),cw1(impoley),
     1             ifhess,cw1(impolexx),cw1(impolexy),
     2             cw1(impoleyy),mpole,dk,rscale,nterms)
           endif
           call taevalsyk(mpole,cw1(impolex),cw1(impoley),cw1(impolexx),
     1             cw1(impolexy),cw1(impoleyy),rscale,nterms,
     2             cw1(ijval),cw1(imptemp),pot,ifgrad,grad,
     3             ifhess,hess,ima4inv)
        endif
c
        return
        end
c
c
c
c
c
c**********************************************************************
      subroutine y2dmpevalall
     $     (dk,rscale,center,mpole,nterms,ztarg,ntarg,
     1      ifpot,pot1,ifgrad,grad1,ifhess,hess1)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an outgoing partial wave expansion.
c               +nterms
c     pot  =      sum   mpole_n H_n(k r) exp(i n theta)
c              n=-nterms  
c     grad = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk     :    Helmholtz parameter
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the multipole expansion
c     ztarg  :    target locations
c     ntarg  :    number of targets
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot    :    potential at ztarg
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        real *8 pot1(*),grad1(2,*),hess1(3,*)
        real *8 center(2),ztarg(2,*),zdiff(2)
        complex *16 zk,mpole(0:nterms),zmul,ztemp1
        real *8 pot, grad(2), hess(3)
        real *8 dk, pi
c
        complex *16, allocatable :: hval(:)
        complex *16, allocatable :: hder(:)

        real *8, allocatable :: kval(:)
c
        complex *16, allocatable :: mpolex(:)
        complex *16, allocatable :: mpoley(:)
c
        complex *16, allocatable :: mpolexx(:)
        complex *16, allocatable :: mpolexy(:)
        complex *16, allocatable :: mpoleyy(:)
c
        complex *16, allocatable :: mptemp(:)
c
        complex *16 ima,ima4,z,pieye2
        data ima/(0.0d0,1.0d0)/
c
c
        ima4=-4*ima
        pi = 4.0d0*datan(1.0d0)
        pieye2=pi*ima/2.0d0
c
        allocate(hval(0:nterms+10), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(kval(0:nterms+10), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(hder(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(mptemp(0:nterms+2), stat=ier)
        if (ier.eq.1) then
          return
        endif
c
c
c
        if( ifgrad .eq. 1 .or. ifhess .eq. 1 ) then
c
           allocate(mpolex(0:nterms+1), stat=ier)
           if (ier.eq.1) then
           return
           endif
           allocate(mpoley(0:nterms+1), stat=ier)
           if (ier.eq.1) then
           return
           endif
c
           do i=0,nterms+1
              mpolex(i)=0
              mpoley(i)=0
           enddo

           mpolex(1) = -dk/rscale*mpole(0)
           mpoley(1) = -dk/rscale*mpole(0)*ima

           do i=1,nterms
              mpolex(i-1) = mpolex(i-1) -dk/2.0d0*mpole(i)*rscale
              mpolex(i+1) = mpolex(i+1) -dk/2.0d0*mpole(i)/rscale
              mpoley(i-1) = mpoley(i-1) +dk/2.0d0*mpole(i)*rscale*ima
              mpoley(i+1) = mpoley(i+1) -dk/2.0d0*mpole(i)/rscale*ima
           enddo
        endif
c
        if( ifhess .eq. 1 ) then
           allocate(mpolexx(0:nterms+2), stat=ier)
           if (ier.eq.1) then
           return
           endif
           allocate(mpolexy(0:nterms+2), stat=ier)
           if (ier.eq.1) then
           return
           endif
           allocate(mpoleyy(0:nterms+2), stat=ier)
           if (ier.eq.1) then
           return
           endif
c
           do i=0,nterms+2
              mpolexx(i)=0
              mpolexy(i)=0
              mpoleyy(i)=0
           enddo

           mpolexx(1) = -dk/1.0d0/rscale*dreal(mpolex(0))
           mpolexy(1) = -dk/1.0d0/rscale*dreal(mpolex(0))*ima
           mpoleyy(1) = -dk/1.0d0/rscale*dreal(mpoley(0))*ima

           do i=1,nterms+1
              mpolexx(i-1) = mpolexx(i-1) -dk/2.0d0*mpolex(i)*rscale
              mpolexx(i+1) = mpolexx(i+1) -dk/2.0d0*mpolex(i)/rscale
              mpolexy(i-1) = mpolexy(i-1) +dk/2.0d0*mpolex(i)*rscale*ima
              mpolexy(i+1) = mpolexy(i+1) -dk/2.0d0*mpolex(i)/rscale*ima
              mpoleyy(i-1) = mpoleyy(i-1) +dk/2.0d0*mpoley(i)*rscale*ima
              mpoleyy(i+1) = mpoleyy(i+1) -dk/2.0d0*mpoley(i)/rscale*ima
           enddo


        endif
c
c

        do k = 1,ntarg

c
           zdiff(1)=ztarg(1,k)-center(1)
           zdiff(2)=ztarg(2,k)-center(2)
           call y2cart2polar(zdiff,r,theta)
c     
           zk = ima*dk
           z=zk*r
           ifder=0
           call y2dall(nterms+3,z,rscale,hval,ifder,hder)
c     


c     convert Hankel H_n to modified Bessel K_n
c     (up to scale of pi/2)
           do n = 0,nterms+2,4
              kval(n) = -dimag(hval(n))
              kval(n+1) = -dreal(hval(n+1))
              kval(n+2) = dimag(hval(n+2))
              kval(n+3) = dreal(hval(n+3))
           enddo

           mptemp(0)=kval(0)
           zmul=exp(ima*theta)
           ztemp1=zmul
           do n=1,nterms+2
              mptemp(n)=dcmplx(kval(n)*dreal(ztemp1),
     1                         kval(n)*dimag(ztemp1))
              ztemp1 = ztemp1*zmul
           enddo

           if (ifpot .eq. 1) then

              pot=kval(0)*dreal(mpole(0))
              do n=1,nterms
                 pot=pot+dreal(mpole(n))*dreal(mptemp(n))
     1                +dimag(mpole(n))*dimag(mptemp(n))
              enddo
              pot=pot*0.25d0

              pot1(k) = pot1(k) + pot
           endif
c     

           if( ifgrad .eq. 1 ) then
              grad(1)=0
              grad(2)=0
c     
              grad(1)=kval(0)*dreal(mpolex(0))
              grad(2)=kval(0)*dreal(mpoley(0))

              do n=1,nterms+1
                 grad(1)=grad(1)+dreal(mpolex(n))*dreal(mptemp(n))+
     1                dimag(mpolex(n))*dimag(mptemp(n))
                 grad(2)=grad(2)+dreal(mpoley(n))*dreal(mptemp(n))+
     1                dimag(mpoley(n))*dimag(mptemp(n))
              enddo
              grad(1)=grad(1)*.25d0
              grad(2)=grad(2)*.25d0

              grad1(1,k) = grad1(1,k)+grad(1)
              grad1(2,k) = grad1(2,k)+grad(2)

           endif
c     
c     
           if( ifhess .eq. 1 ) then
              hess(1)=0
              hess(2)=0
              hess(3)=0
c     
              hess(1)=kval(0) *dreal(mpolexx(0))
              hess(2)=kval(0) *dreal(mpolexy(0))
              hess(3)=kval(0) *dreal(mpoleyy(0))
              do n=1,nterms+2
                 hess(1)=hess(1)+
     1                dreal(mpolexx(n))*dreal(mptemp(n))
     2                +dimag(mpolexx(n))*dimag(mptemp(n))
                 hess(2)=hess(2)+
     1                dreal(mpolexy(n))*dreal(mptemp(n))
     2                +dimag(mpolexy(n))*dimag(mptemp(n))
                 hess(3)=hess(3)+
     1                dreal(mpoleyy(n))*dreal(mptemp(n))
     2                +dimag(mpoleyy(n))*dimag(mptemp(n))
              enddo
              hess(1)=hess(1)*.25d0
              hess(2)=hess(2)*.25d0
              hess(3)=hess(3)*.25d0

              hess1(1,k) = hess1(1,k)+hess(1)
              hess1(2,k) = hess1(2,k)+hess(2)
              hess1(3,k) = hess1(3,k)+hess(3)

           endif

        enddo

        return
        end
c
c
c
c
c
c**********************************************************************
      subroutine y2dtaevalall
     $     (dk,rscale,center,mpole,nterms,ztarg,ntarg,
     1      ifpot,pot1,ifgrad,grad1,ifhess,hess1)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an incoming partial wave expansion.
c               +nterms
c     pot  =      sum   I_n(k r) [real(mpole_n)*cos(n theta)
c                 n=0            +imag(mpole_n)*sin(n theta)]
c
c     grad = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk     :    Helmholtz parameter
c     rscale :    scaling parameter 
c     center :    expansion center
c     mpole  :    multipole expansion 
c     nterms :    order of the multipole expansion
c     ztarg  :    target location
c     ntarg  :    number of targets
c     ifgrad :   flag controlling evaluation of gradient:
c                   ifgrad = 0, do not compute gradient.
c                   ifgrad = 1, compute gradient.
c     ifhess :   flag for computing Hessian:
c	            ifhess = 0 -> don't compute 
c		    ifhess = 1 -> do compute 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     ier    :    error return code
c                     ier=0  successful execution
c                     ier=8  insufficient work space (not implemented yet)
c                     ier=16 ztarg too close to center
c     pot    :    potential at ztarg
c     grad   :    gradient at ztarg (if requested)
c     hess   :    hessian at ztarg (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 zk
        real *8 dk
        real *8 pot,grad(2),hess(3)
        complex *16 mpole(0:nterms)
        real *8 pot1(*),grad1(2,*),hess1(3,*)
        real *8 center(2),ztarg(2,1),zdiff(2)
c
        complex *16, allocatable :: jval(:)
        real *8, allocatable :: ival(:)
        complex *16, allocatable :: jder(:)
        integer, allocatable :: iscale(:)
c
        complex *16, allocatable :: mpolex(:)
        complex *16, allocatable :: mpoley(:)
c
        complex *16, allocatable :: mpolexx(:)
        complex *16, allocatable :: mpolexy(:)
        complex *16, allocatable :: mpoleyy(:)
c
        complex *16, allocatable :: mptemp(:)
c
        complex *16 ima,ima4,ima4inv,z,z1scale,z2scale,z3scale,z4scale
        complex *16 zmul,zinv,ztemp1,ztemp2
        data ima/(0.0d0,1.0d0)/
c
c
        ima4=-4*ima
        ima4inv=ima/4

        zk = ima*dk

c
        lwfjs = nterms+5 + 4*nterms + 100
        allocate(jval(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(ival(0:nterms+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jder(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(iscale(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
c
        if( ifgrad .eq. 1 .or. ifhess .eq. 1 ) then
           allocate(mpolex(0:nterms+1), stat=ier)
           if (ier.eq.1) then
           return
           endif
           allocate(mpoley(0:nterms+1), stat=ier)
           if (ier.eq.1) then
           return
           endif

           if( ifhess .eq. 1 ) then
              allocate(mpolexx(0:nterms+2), stat=ier)
              if (ier.eq.1) then
                 return
              endif
              allocate(mpolexy(0:nterms+2), stat=ier)
              if (ier.eq.1) then
                 return
              endif
              allocate(mpoleyy(0:nterms+2), stat=ier)
              if (ier.eq.1) then
                 return
              endif
           endif
c
           call mkmpole12yk(mpolex,mpoley,
     1          ifhess,mpolexx,mpolexy,
     2          mpoleyy,mpole,dk,rscale,nterms)
           

        endif
c
c
        allocate(mptemp(0:nterms+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
c
        do k=1,ntarg
           zdiff(1)=ztarg(1,k)-center(1)
           zdiff(2)=ztarg(2,k)-center(2)
           call y2cart2polar(zdiff,r,theta)
c
           z=zk*r
           ifder=0
           call jfuns2d(ier,nterms+3,z,rscale,jval,ifder,jder,
     1           lwfjs,iscale,ntop)
c

           call mkmptempyk(mptemp,theta,jval,ival,nterms)

           pot = 0.0d0
           grad(1) = 0.0d0
           grad(2) = 0.0d0
           hess(1) = 0.0d0
           hess(2) = 0.0d0
           hess(3) = 0.0d0

           call taevalsyk(mpole,mpolex,mpoley,mpolexx,
     1             mpolexy,mpoleyy,rscale,nterms,
     2             jval,mptemp,pot,ifgrad,grad,
     3             ifhess,hess,ima4inv)


           if( ifpot .eq. 1 ) then
              pot1(k)=pot1(k)+pot
           endif
           if( ifgrad .eq. 1 ) then
              grad1(1,k)=grad1(1,k)+grad(1)
              grad1(2,k)=grad1(2,k)+grad(2)
           endif
           if( ifhess .eq. 1 ) then
              hess1(1,k)=hess1(1,k)+hess(1)
              hess1(2,k)=hess1(2,k)+hess(2)
              hess1(3,k)=hess1(3,k)+hess(3)
           endif
        enddo
        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine y2dformmp(ier,dk,rscale,source,charge,ns,center,
     1                       nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a multipole (k) expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c               
c     mpole_n  =  2 * sum charge_j  I_n(k r) exp(i n theta_j) /rscale^n
c                      j 
c
c     mpole_0 does not include the factor of 2
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     dk              : Yukawa parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     charge(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the k-expansion
c
        complex *16 zk,mpole(0:nterms)
        real *8 charge(*)
        real *8 dk
        real *8 center(2),source(2,1),zdiff(2)

        complex *16, allocatable :: jval(:)
        real *8, allocatable :: ival(:)
        complex *16, allocatable :: jder(:)
        integer, allocatable :: iscale(:)
        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,ima4,z
        data ima/(0.0d0,1.0d0)/
c
c
ccc        lwfjs = min(10000,nterms+5 + 4*nterms + 100)
        lwfjs = nterms+5 + 4*nterms + 100
        allocate(jval(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(ival(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jder(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(iscale(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
c
        do n=0,nterms
           mpole(n)=0
        enddo

        zk = ima*dk

        do j=1,ns
           zdiff(1)=source(1,j)-center(1)
           zdiff(2)=source(2,j)-center(2)
           call y2cart2polar(zdiff,r,theta)
           z=zk*r
           ifder=0
           write(*,*) 'z = ', z
           write(*,*) 'rscale = ', rscale
           call jfuns2d(ier,nterms+1,z,rscale,jval,ifder,jder,
     1           lwfjs,iscale,ntop)

           write(*,*) 'ier = ', ier

c     convert Bessel J_n to modified Bessel I_n
           do n = 0,nterms,4
              ival(n) = dreal(jval(n))
              ival(n+1) = dimag(jval(n+1))
              ival(n+2) = -dreal(jval(n+2))
              ival(n+3) = -dimag(jval(n+3))
           enddo
           mpole(0)=mpole(0)+charge(j)*jval(0)
           zmul=exp(ima*theta)
           ztemp1= zmul*2.0d0*charge(j)
           do n=1,nterms
              mpole(n)=mpole(n)+ival(n)*ztemp1
              ztemp1= ztemp1*zmul
           enddo
        enddo
        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine y2dformmp_add(ier,dk,rscale,source,charge,ns,center,
     1                       nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a multipole (h) expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c               
c     mpole_n  =  sum charge_j  J_n(k r) exp(-i n theta_j) /rscale_j^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     charge(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the h-expansion
c
        complex *16 zk,mpole(0:nterms)
        real *8 charge(*)
        real *8 dk
        real *8 center(2),source(2,1),zdiff(2)

        complex *16, allocatable :: jval(:)
        real *8, allocatable :: ival(:)
        complex *16, allocatable :: jder(:)
        integer, allocatable :: iscale(:)
        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,ima4,z
        data ima/(0.0d0,1.0d0)/
c
c
ccc        lwfjs = min(10000,nterms+5 + 4*nterms + 100)
        lwfjs = nterms+5 + 4*nterms + 100
        allocate(jval(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(ival(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jder(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(iscale(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
c
        zk = ima*dk

        do j=1,ns
           zdiff(1)=source(1,j)-center(1)
           zdiff(2)=source(2,j)-center(2)
           call y2cart2polar(zdiff,r,theta)
           z=zk*r
           ifder=0
           call jfuns2d(ier,nterms+1,z,rscale,jval,ifder,jder,
     1           lwfjs,iscale,ntop)

c     convert Bessel J_n to modified Bessel I_n
           do n = 0,nterms,4
              ival(n) = dreal(jval(n))
              ival(n+1) = dimag(jval(n+1))
              ival(n+2) = -dreal(jval(n+2))
              ival(n+3) = -dimag(jval(n+3))
           enddo
           mpole(0)=mpole(0)+charge(j)*jval(0)
           zmul=exp(ima*theta)
           ztemp1= zmul*2.0d0*charge(j)
           do n=1,nterms
              mpole(n)=mpole(n)+ival(n)*ztemp1
              ztemp1= ztemp1*zmul
           enddo
        enddo
        return
        end
c
c
c
c
C***********************************************************************
        subroutine y2dformta(ier,dk,rscale,source,charge,ns,center,
     1                       nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a local (j) expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c               
c     mpole_n  =  2*sum charge_j  K_n(k r) exp(i n theta_j) *rscale^n
c                    j  
c
c     for n = 0 the above is divided by 2
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     charge(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : epxansion center
c     nterms          : order of local expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the j-expansion
c
        complex *16 zk,mpole(0:nterms)
        real *8 charge(*), dk
        real *8 center(2),source(2,1),zdiff(2)

        complex *16, allocatable :: hval(:)
        complex *16, allocatable :: hder(:)
        real *8, allocatable :: kval(:)
        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,ima4,z
        data ima/(0.0d0,1.0d0)/
c
c
        allocate(hval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(hder(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(kval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif

        zk = ima*dk
c
        do n=0,nterms
           mpole(n)=0
        enddo

        do j=1,ns
           zdiff(1)=source(1,j)-center(1)
           zdiff(2)=source(2,j)-center(2)
           call y2cart2polar(zdiff,r,theta)
           z=zk*r
           ifder=0
           call y2dall(nterms+1,z,rscale,hval,ifder,hder)
           do n = 0,nterms,4
              kval(n) = -dimag(hval(n))
              kval(n+1) = -dreal(hval(n+1))
              kval(n+2) = dimag(hval(n+2))
              kval(n+3) = dreal(hval(n+3))
           enddo
           mpole(0)=mpole(0)+charge(j)*kval(0)
           zmul=exp(ima*theta)
           ztemp1= zmul*charge(j)*2.0d0
           do n=1,nterms
              mpole(n)=mpole(n)+kval(n)*ztemp1
              ztemp1= ztemp1*zmul
           enddo
        enddo
        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine y2dformta_add(ier,dk,rscale,source,charge,ns,center,
     1                       nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a local (j) expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c               
c     mpole_n  =  2*sum charge_j  K_n(k r) exp(i n theta_j) *rscale^n
c                    j  
c
c     for n = 0 the above is divided by 2
c
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     charge(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : epxansion center
c     nterms          : order of local expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the j-expansion
c
        complex *16 zk,mpole(0:nterms)
        real *8 charge(*), dk
        real *8 center(2),source(2,1),zdiff(2)

        complex *16, allocatable :: hval(:)
        complex *16, allocatable :: hder(:)
        real *8, allocatable :: kval(:)
        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,ima4,z
        data ima/(0.0d0,1.0d0)/
c
c
        allocate(hval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(hder(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(kval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif

        zk = ima*dk
c
        do j=1,ns
           zdiff(1)=source(1,j)-center(1)
           zdiff(2)=source(2,j)-center(2)
           call y2cart2polar(zdiff,r,theta)
           z=zk*r
           ifder=0
           call y2dall(nterms+1,z,rscale,hval,ifder,hder)
           do n = 0,nterms,4
              kval(n) = -dimag(hval(n))
              kval(n+1) = -dreal(hval(n+1))
              kval(n+2) = dimag(hval(n+2))
              kval(n+3) = dreal(hval(n+3))
           enddo
           mpole(0)=mpole(0)+charge(j)*kval(0)
           zmul=exp(ima*theta)
           ztemp1= zmul*charge(j)*2.0d0
           do n=1,nterms
              mpole(n)=mpole(n)+kval(n)*ztemp1
              ztemp1= ztemp1*zmul
           enddo
        enddo
        return
        end
c
c
c
        subroutine y2dmpmp(dk,
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts multipole expansion to a multipole expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           zk      = Helmholtz parameter
C           rscale1 = scaling parameter for original multipole expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original multipole expansion
C           rscale2 = scaling parameter for shifted multipole expansion
C           center2 = center of shifted multipole expansion
C           nterms2 = order of shifted multipole expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted multipole expansion
c
        real *8 dk
        complex *16 zk,hexp(0:nterms1),jexp(0:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,ima, zmul,zinv,ztemp1,ztemp2
        complex *16, allocatable :: jval(:), jder(:), jtemp(:)
        real *8, allocatable :: ival(:)
        integer, allocatable :: iscale(:)
c
        data ima/(0.0d0,1.0d0)/
c

        zk = ima*dk

        done=1
        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        lwfjs = nterms+5 + 4*nterms + 100
        allocate(jval(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(ival(0:nterms+3), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jder(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(iscale(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jtemp(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
        call y2cart2polar(zdiff,r,theta)
        theta=theta-pi
        z=zk*r
        ifder=0
        call jfuns2d(ier,nterms+3,z,rscale1,jval,ifder,jder,
     1        lwfjs,iscale,ntop)
c

c     convert Bessel J_n to modified Bessel I_n
        do n = 0,nterms,4
           ival(n) = dreal(jval(n))
           ival(n+1) = dimag(jval(n+1))
           ival(n+2) = -dreal(jval(n+2))
           ival(n+3) = -dimag(jval(n+3))
        enddo


        do i = 0,nterms2
           jexp(i) = 0
        enddo
c
        jtemp(0) = ival(0)
        zmul=exp(ima*theta)
        ztemp1= zmul
        do j = 1,nterms
           jtemp( j) = ztemp1*ival(j)
           ztemp1= ztemp1*zmul
        enddo
c
        jexp(0) = jexp(0) + hexp(0)*jtemp(0)
        rsj=rscale1
        do j = 1,nterms1
           jexp(0) = jexp(0)+(dreal(hexp(j))*dreal(jtemp(j))
     1          +dimag(hexp(j))*dimag(jtemp(j)))*rsj**2
           rsj=rsj*rscale1
        enddo
c
        rsi=rscale1
        rsi7=rscale2
        rsi5=rscale1/rscale2
        do i = 1,nterms2
           jexp(i) = jexp(i) + 2.0d0*(dreal(hexp(0))*dreal(jtemp(i))
     1          +ima*dreal(hexp(0))*dimag(jtemp(i)))*rsi5
           rsj=rscale1
        do j = 1,min(nterms1,i)
           jexp(i) = jexp(i)+(dreal(hexp(j))*dreal(jtemp(i-j))
     1          -dimag(hexp(j))*dimag(jtemp(i-j))
     2          +ima*(dreal(hexp(j))*dimag(jtemp(i-j))
     3          +dimag(hexp(j))*dreal(jtemp(i-j))))*rsi5
           jexp(i) = jexp(i)+(dreal(hexp(j))*dreal(jtemp(i+j))
     1          +dimag(hexp(j))*dimag(jtemp(i+j))
     2          +ima*(dreal(hexp(j))*dimag(jtemp(i+j))
     3          -dimag(hexp(j))*dreal(jtemp(i+j))))*rsj**2*rsi5
           rsj=rsj*rscale1
        enddo
        rsj=rscale1**(i+1)
        fs2=rsi5*rscale1**2
        do j = i+1,nterms1
           jexp(i) = jexp(i)+(dreal(hexp(j))*dreal(jtemp(j-i))
     1          +dimag(hexp(j))*dimag(jtemp(j-i))
     2          -ima*(dreal(hexp(j))*dimag(jtemp(j-i))
     3          -dimag(hexp(j))*dreal(jtemp(j-i))))*fs2
           jexp(i) = jexp(i)+(dreal(hexp(j))*dreal(jtemp(i+j))
     1          +dimag(hexp(j))*dimag(jtemp(i+j))
     2          +ima*(dreal(hexp(j))*dimag(jtemp(i+j))
     3          -dimag(hexp(j))*dreal(jtemp(i+j))))*rsj**2*rsi5
           rsj=rsj*rscale1
           fs2=fs2*rscale1**2
        enddo
        rsi=rsi*rscale1
        rsi7=rsi7*rscale2
        rsi5=rsi5*rscale1/rscale2
        enddo
        return
        end
c
c
c
c
c
        subroutine y2dlocloc(dk,
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts local expansion to a local expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           zk      = Helmholtz parameter
C           rscale1 = scaling parameter for original multipole expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original multipole expansion
C           rscale2 = scaling parameter for shifted multipole expansion
C           center2 = center of shifted multipole expansion
C           nterms2 = order of shifted multipole expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted multipole expansion
c
        real *8 dk
        complex *16 zk,hexp(0:nterms1),jexp(0:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,ima, zmul,zinv,ztemp1,ztemp2
        complex *16, allocatable :: jval(:), jder(:), jtemp(:)
        real *8, allocatable :: ival(:)
        integer, allocatable :: iscale(:)
c
        data ima/(0.0d0,1.0d0)/
c

        zk = ima*dk

        done=1
        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        lwfjs = nterms+5 + 4*nterms + 100
        allocate(jval(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(ival(0:nterms+3), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jder(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(iscale(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jtemp(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
        call y2cart2polar(zdiff,r,theta)
        theta=theta-pi
        z=zk*r
        ifder=0
        call jfuns2d(ier,nterms+3,z,rscale1,jval,ifder,jder,
     1        lwfjs,iscale,ntop)
c

c     convert Bessel J_n to modified Bessel I_n
        do n = 0,nterms,4
           ival(n) = dreal(jval(n))
           ival(n+1) = dimag(jval(n+1))
           ival(n+2) = -dreal(jval(n+2))
           ival(n+3) = -dimag(jval(n+3))
        enddo


        do i = 0,nterms2
           jexp(i) = 0
        enddo
c
        jtemp(0) = ival(0)
        zmul=-exp(ima*theta)
        ztemp1= zmul
        do j = 1,nterms
           jtemp( j) = ztemp1*ival(j)
           ztemp1= ztemp1*zmul
        enddo
c
        jexp(0) = jexp(0) + hexp(0)*jtemp(0)
        rsj=rscale1
        do j = 1,nterms1
           jexp(0) = jexp(0)+(dreal(hexp(j))*dreal(jtemp(j))
     1          +dimag(hexp(j))*dimag(jtemp(j)))
        enddo
c
        rsi=rscale1
        rsi7=rscale2
        rsi5=rscale2/rscale1
        do i = 1,nterms2
           jexp(i) = jexp(i) + 2.0d0*(dreal(hexp(0))*dreal(jtemp(i))
     1          +ima*dreal(hexp(0))*dimag(jtemp(i)))*rsi7*rsi
           fs2=rsi5
           if( nterms1 .le. i ) fs2=fs2*rscale1**(2*(i-nterms1))
           do j = min(nterms1,i),1,-1
              jexp(i) = jexp(i)+(dreal(hexp(j))*dreal(jtemp(i-j))
     1             -dimag(hexp(j))*dimag(jtemp(i-j))
     2             +ima*(dreal(hexp(j))*dimag(jtemp(i-j))
     3             +dimag(hexp(j))*dreal(jtemp(i-j))))*fs2
              jexp(i) = jexp(i)+(dreal(hexp(j))*dreal(jtemp(i+j))
     1             +dimag(hexp(j))*dimag(jtemp(i+j))
     2             +ima*(dreal(hexp(j))*dimag(jtemp(i+j))
     3             -dimag(hexp(j))*dreal(jtemp(i+j))))*rsi*rsi7
              fs2=fs2*rscale1**2
           enddo
           do j = i+1,nterms1
              jexp(i) = jexp(i)+(dreal(hexp(j))*dreal(jtemp(j-i))
     1             +dimag(hexp(j))*dimag(jtemp(j-i))
     2             -ima*(dreal(hexp(j))*dimag(jtemp(j-i))
     3             -dimag(hexp(j))*dreal(jtemp(j-i))))*rsi5
              jexp(i) = jexp(i)+(dreal(hexp(j))*dreal(jtemp(i+j))
     1             +dimag(hexp(j))*dimag(jtemp(i+j))
     2             +ima*(dreal(hexp(j))*dimag(jtemp(i+j))
     3             -dimag(hexp(j))*dreal(jtemp(i+j))))*rsi*rsi7
              
           enddo
           rsi=rsi*rscale1
           rsi7=rsi7*rscale2
           rsi5=rsi5*rscale2/rscale1
        enddo
        return
        end
c
c
c
c
c
        subroutine y2dmploc(dk,
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           zk      = Helmholtz parameter
C           rscale1 = scaling parameter for original multipole expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original multipole expansion
C           rscale2 = scaling parameter for shifted local expansion
C           center2 = center of shifted local expansion
C           nterms2 = order of shifted local expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted local expansion
c
        real *8 dk
        complex *16 zk,hexp(0:nterms1),jexp(0:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,ima, zmul,zinv,ztemp1,ztemp2
        complex *16, allocatable :: hval(:), hder(:), htemp(:)
        real *8, allocatable :: kval(:)
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        allocate(hval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(hder(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(htemp(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(kval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
        call y2cart2polar(zdiff,r,theta)
c
        theta=theta-pi
c
        zk = ima*dk
        z=zk*r
        ifder=0
        call y2dall(nterms+3,z,rscale1,hval,ifder,hder)


c     convert Hankel H_n to modified Bessel K_n
c     (up to scale of pi/2)
        do n = 0,nterms,4
           kval(n) = -dimag(hval(n))
           kval(n+1) = -dreal(hval(n+1))
           kval(n+2) = dimag(hval(n+2))
           kval(n+3) = dreal(hval(n+3))
        enddo

c        
        do i = 0,nterms2
           jexp(i) = 0
        enddo
c
        htemp(0) = kval(0)
        zmul=-exp(ima*theta)
        ztemp1= zmul
        do j = 1,nterms
           htemp( j) = ztemp1*kval(j)
           ztemp1= ztemp1*zmul
        enddo
c
        jexp(0) = jexp(0) + hexp(0)*htemp(0)
        do j = 1,nterms1
           jexp(0) = jexp(0)+(dreal(hexp(j))*dreal(htemp(j))
     1          +dimag(hexp(j))*dimag(htemp(j)))
        enddo
c
        rsi=rscale1
        rsi2=rscale1**2
        rsi5=rscale2/rscale1
        rsi52=rsi5*rscale2/rscale1
        isign = -1
        do i = 1,nterms2
           jexp(i) = jexp(i) + 2.0d0*(dreal(hexp(0))*dreal(htemp(i))
     1          +ima*dreal(hexp(0))*dimag(htemp(i)))
           rsj=rscale1
           rsj2=rscale1**2
           do j = 1,min(nterms1,i)
              jexp(i) = jexp(i)+(dreal(hexp(j))*dreal(htemp(i-j))
     1             -dimag(hexp(j))*dimag(htemp(i-j))
     2             +ima*(dreal(hexp(j))*dimag(htemp(i-j))
     3             +dimag(hexp(j))*dreal(htemp(i-j))))*rsj2
              jexp(i) = jexp(i)+(dreal(hexp(j))*dreal(htemp(i+j))
     1             +dimag(hexp(j))*dimag(htemp(i+j))
     2             +ima*(dreal(hexp(j))*dimag(htemp(i+j))
     3             -dimag(hexp(j))*dreal(htemp(i+j))))

              rsj=rsj*rscale1
              rsj2=rsj2*rscale1**2
           enddo
           do j = i+1,nterms1
              jexp(i) = jexp(i)+(dreal(hexp(j))*dreal(htemp(j-i))
     1             +dimag(hexp(j))*dimag(htemp(j-i))
     2             -ima*(dreal(hexp(j))*dimag(htemp(j-i))
     3             -dimag(hexp(j))*dreal(htemp(j-i))))*rsi2
              jexp(i) = jexp(i)+(dreal(hexp(j))*dreal(htemp(i+j))
     1             +dimag(hexp(j))*dimag(htemp(i+j))
     2             +ima*(dreal(hexp(j))*dimag(htemp(i+j))
     3             -dimag(hexp(j))*dreal(htemp(i+j))))

           enddo
           jexp(i)=jexp(i)*isign
           isign = -isign
           rsi=rsi*rscale1
           rsi2=rsi2*rscale1**2
        enddo
        rsi=rscale2/rscale1
        do i = 1,nterms2
           jexp(+i) = jexp(+i)*rsi
           rsi=rsi*rscale2/rscale1
        enddo
        return
        end
c
c
c
c
c
        subroutine y2dmploc_add(dk,
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts multipole expansion to a local expansion.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           zk      = Helmholtz parameter
C           rscale1 = scaling parameter for original multipole expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original multipole expansion
C           rscale2 = scaling parameter for shifted local expansion
C           center2 = center of shifted local expansion
C           nterms2 = order of shifted local expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted local expansion
c
        real *8 dk
        complex *16 zk,hexp(0:nterms1),jexp(0:nterms2)
        real *8 center1(2),center2(2),zdiff(2)
        complex *16 z,ima, zmul,zinv,ztemp1,ztemp2
        complex *16, allocatable :: hval(:), hder(:), htemp(:)
        real *8, allocatable :: kval(:)
        data ima/(0.0d0,1.0d0)/
c
        done=1
        pi=4*atan(done)
c
        nterms = nterms1+nterms2
c
        allocate(hval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(hder(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(htemp(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(kval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
c
        zdiff(1)=center2(1)-center1(1)
        zdiff(2)=center2(2)-center1(2)
        call y2cart2polar(zdiff,r,theta)
c
        theta=theta-pi
c
        zk = ima*dk
        z=zk*r
        ifder=0
        call y2dall(nterms+3,z,rscale1,hval,ifder,hder)


c     convert Hankel H_n to modified Bessel K_n
c     (up to scale of pi/2)
        do n = 0,nterms,4
           kval(n) = -dimag(hval(n))
           kval(n+1) = -dreal(hval(n+1))
           kval(n+2) = dimag(hval(n+2))
           kval(n+3) = dreal(hval(n+3))
        enddo

c        
c
        htemp(0) = kval(0)
        zmul=-exp(ima*theta)
        ztemp1= zmul
        do j = 1,nterms
           htemp( j) = ztemp1*kval(j)
           ztemp1= ztemp1*zmul
        enddo
c
        jexp(0) = jexp(0) + hexp(0)*htemp(0)
        do j = 1,nterms1
           jexp(0) = jexp(0)+(dreal(hexp(j))*dreal(htemp(j))
     1          +dimag(hexp(j))*dimag(htemp(j)))
        enddo
c
        rsi=rscale1
        rsi2=rscale1**2
        rsi5=rscale2/rscale1
        rsi52=rsi5*rscale2/rscale1
        isign = -1
        do i = 1,nterms2
           jexp(i) = jexp(i)+isign*2.0d0*(dreal(hexp(0))*dreal(htemp(i))
     1          +ima*dreal(hexp(0))*dimag(htemp(i)))*rsi5
           rsj=rscale1
           rsj2=rscale1**2
           do j = 1,min(nterms1,i)
              jexp(i) = jexp(i)+isign*(dreal(hexp(j))*dreal(htemp(i-j))
     1             -dimag(hexp(j))*dimag(htemp(i-j))
     2             +ima*(dreal(hexp(j))*dimag(htemp(i-j))
     3             +dimag(hexp(j))*dreal(htemp(i-j))))*rsj2*rsi5
              jexp(i) = jexp(i)+isign*(dreal(hexp(j))*dreal(htemp(i+j))
     1             +dimag(hexp(j))*dimag(htemp(i+j))
     2             +ima*(dreal(hexp(j))*dimag(htemp(i+j))
     3             -dimag(hexp(j))*dreal(htemp(i+j))))*rsi5

              rsj=rsj*rscale1
              rsj2=rsj2*rscale1**2
           enddo
           do j = i+1,nterms1
              jexp(i) = jexp(i)+isign*(dreal(hexp(j))*dreal(htemp(j-i))
     1             +dimag(hexp(j))*dimag(htemp(j-i))
     2             -ima*(dreal(hexp(j))*dimag(htemp(j-i))
     3             -dimag(hexp(j))*dreal(htemp(j-i))))*rsi2*rsi5
              jexp(i) = jexp(i)+isign*(dreal(hexp(j))*dreal(htemp(i+j))
     1             +dimag(hexp(j))*dimag(htemp(i+j))
     2             +ima*(dreal(hexp(j))*dimag(htemp(i+j))
     3             -dimag(hexp(j))*dreal(htemp(i+j))))*rsi5

           enddo

           isign = -isign
           rsi=rsi*rscale1
           rsi2=rsi2*rscale1**2
           rsi5 = rsi5*rscale2/rscale1
        enddo

        return
        end
c
c
c
c
        subroutine y2dmploc_add_trunc(dk,
     $     rscale1,center1,hexp,nterms1,nterms1_trunc,
     $     rscale2,center2,jexp,nterms2,nterms2_trunc)
        implicit real *8 (a-h,o-z)
C
C     Usage:
C
C           Converts multipole expansion to a local expansion
C           with truncation, and add to jexp.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           zk      = Helmholtz parameter
C           rscale1 = scaling parameter for original multipole expansion
C           center1 = center of original multiple expansion
C           hexp    = coefficients of original multiple expansion
C           nterms1 = order of original multipole expansion
C           nterms1_trunc = number of terms used in the multipole expansion
C           rscale2 = scaling parameter for shifted local expansion
C           center2 = center of shifted local expansion
C           nterms2 = order of shifted local expansion
C           nterms2_trunc = number of terms used in the shifted local expansion
C
C     OUTPUT:
C
C           jexp    = coefficients of shifted local expansion
c
c
        real *8 dk
        complex *16 zk,hexp(0:nterms1),jexp(0:nterms2)
        real *8 center1(2),center2(2)
        complex *16, allocatable :: hexp1(:), jexp1(:)

        allocate(hexp1(0:nterms1_trunc), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(jexp1(0:nterms2_trunc), stat=ier)
        if (ier.eq.1) then
        return
        endif
c
        do i = 0,nterms1_trunc
           hexp1(i) = hexp(i)
        enddo
c
        call y2dmploc(dk,
     $     rscale1,center1,hexp1,nterms1_trunc,
     $     rscale2,center2,jexp1,nterms2_trunc)
c
        do i = 0,nterms2_trunc
           jexp(i) = jexp(i) + jexp1(i)
        enddo
        return
        end
c
c
c
c
c
c
c**********************************************************************
      subroutine ypotgrad2dall(ifgrad,ifhess,sources,charge,ns,
     1                   target,dk,pot,grad,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD, and
c     Hessian at the target point TARGET, due to a collection of 
c     charges at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot = H_0(k*r)*(eye/4)
c		grad = gradient
c		hess = Hessian
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifgrad         : flag for computing gradient
c	                 	   ifgrad = 0 -> don't compute 
c		                   ifgrad = 1 -> do compute 
c     ifhess        : flag for computing Hessian
c	                 	   ifhess = 0 -> don't compute 
c		                   ifhess = 1 -> do compute 
c     sources(2,*)  : location of the sources
c     charge        : charge strengths
c     ns            : number of sources
c     target        : location of the target
c     dk            : Yukawa parameter
c     wavek         : Helmholtz parameter
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (complex *16)      : calculated potential
c     grad  (complex *16)      : calculated gradient
c     hess  (complex *16)      : calculated Hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 dk
      real *8 sources(2,ns),target(2)
      complex *16 wavek
      real *8 pot,grad(2),hess(3),potloc,gradloc(2),hessloc(3)
      complex *16 h0,h1,cd,eye,z
      real *8 charge(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      wavek = eye*dk

      pot = 0.0d0
      if (ifgrad.eq.1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0
      endif
      if (ifhess.eq.1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0
      endif
c
      do i = 1,ns
         call ypotgrad2d(ifgrad,ifhess,sources(1,i),charge(i),target,
     1        dk,potloc,gradloc,hessloc)
         pot = pot + potloc
         if (ifgrad.eq.1) then
         grad(1) = grad(1) + gradloc(1)
         grad(2) = grad(2) + gradloc(2)
         endif
         if (ifhess.eq.1) then
         hess(1) = hess(1) + hessloc(1)
         hess(2) = hess(2) + hessloc(2)
         hess(3) = hess(3) + hessloc(3)
         endif
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine ypotgrad2d(ifgrad,ifhess,source,charge,target,dk,
     1                pot,grad,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD
c     and Hessian HESS at the target point TARGET, due to a charge at 
c     SOURCE. The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot = H_0(k*r)*(eye/4)
c		grad = gradient(pot)
c		hess = Hessian
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     target    : location of the target
c     wavek     : Helmholtz parameter
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      real *8 source(2),target(2)
      complex *16 wavek
      real *8 pot,grad(2),hess(3)
      real *8 charge
      complex *16 z, h0, h1, h2, cd, h2z, zk, ima, ima4, ima4inv
c
      data ima/(0.0d0,1.0d0)/

      wavek = ima*dk

c
c ... Calculate offsets and distance
c
      ima4=-4*ima
      ima4inv=ima/4
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      rr=xdiff*xdiff+ydiff*ydiff
      r=sqrt(rr)
      z=wavek*r
c
      ifexpon = 1
      call hank103(z, h0, h1, ifexpon)
      pot = dreal(h0*charge*ima4inv)
c
      if (ifgrad.eq.1) then
         cd = -h1*(wavek*charge*ima4inv/r)
         grad(1) = dreal(cd*xdiff)
         grad(2) = dreal(cd*ydiff)
c         ctheta=xdiff/r
c         stheta=ydiff/r
c         cd = -h1*(wavek*charge*ima4inv)
c         grad(1) = cd*ctheta
c         grad(2) = cd*stheta
      endif
c
      if (ifhess.eq.1) then
         cd = (wavek*charge*ima4inv/r)/rr
         h2z=(-z*h0+2*h1)
         hess(1) = dreal(cd*(h2z*xdiff*xdiff-rr*h1))
         hess(2) = dreal(cd*(h2z*xdiff*ydiff      ))
         hess(3) = dreal(cd*(h2z*ydiff*ydiff-rr*h1))
c         ctheta=xdiff/r
c         stheta=ydiff/r
c         cd = (wavek*charge*ima4inv/r)
c         h2z=(-z*h0+2*h1)
c         hess(1) = cd*(h2z*ctheta*ctheta-h1)
c         hess(2) = cd*(h2z*ctheta*stheta   )
c         hess(3) = cd*(h2z*stheta*stheta-h1)
      endif
c
      return
      end
c
c
c
c**********************************************************************
      subroutine y2cart2polar(zat,r,theta)
c**********************************************************************
c
c     Convert form Cartesian to polar coordinates.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c       zat   :  Cartesian vector
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c       r      :  |zat|
c       theta  : angle of (zat(1),zat(2)) subtended with 
c                 respect to x-axis
c
c-----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 zat(2)
      complex *16 eye
      data eye/(0.0d0,1.0d0)/
c
c 
      r= sqrt(zat(1)**2+zat(2)**2)
      if( abs(zat(1)) .eq. 0 .and. abs(zat(2)) .eq. 0 ) then
      theta = 0
      else
      theta = datan2(zat(2),zat(1))
      endif
      return
      end
c
c
c
c**********************************************************************
      subroutine y2dall(nterms,z,rscale,hvec,ifder,hder)
c**********************************************************************
c
c     This subroutine computes scaled versions of the Hankel 
c     functions H_n of orders 0 to nterms.
c
c       	hvec(n)= H_n(z)*rscale^(n)
c
c     The parameter SCALE is useful when |z| < 1, in which case
c     it damps out the rapid growth of H_n as n increases. In such 
c     cases, we recommend setting 
c                                 
c               rscale approx |z|
c
c     or something close. If |z| > 1, set scale = 1.
c
c     If the flag IFDER is set to one, it also computes the 
c     derivatives of h_n.
c
c		hder(n)= H_n'(z)*rscale^(n)
c
c     NOTE: If |z| < 1.0d-200, the subroutine returns zero.
c     
c-----------------------------------------------------------------------
c     INPUT:
c
c     nterms  : highest order of the Hankel functions to be computed.
c     z       : argument of the Hankel functions.
c     rscale   : scaling parameter discussed above
c     ifder   : flag indcating whether derivatives should be computed.
c		ifder = 1   ==> compute 
c		ifder = 0   ==> do not compute
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     hvec    : the vector of Hankel functions 
c     hder    : the derivatives of the Hankel functions 
c
c-----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      complex *16 hvec(0:*),hder(0:*)
      complex *16 zk2,z,zinv,ztmp,fhextra,h0,h1
c
      data thresh/1.0d-200/,done/1.0d0/
c
c     If |z| < thresh, return zeros.
c
      if (abs(z).lt.thresh) then
         do i=0,nterms
            hvec(i)=0
            hder(i)=0
         enddo
         return
      endif
c
c     Otherwise, get H_0 and H_1 via hank103 and the rest via
c     recursion.
c
      ifexpon=1
      call hank103(z,h0,h1,ifexpon)
      hvec(0)=h0
      hvec(1)=h1*rscale
c
c
c     From Abramowitz and Stegun (9.1.27)
c
c       H_{n-1}(z) + H_{n+1}(z) = 2*n/z * H_n(z)
c
c     With scaling:
c
c       H_{n-1}(z) *rscale + H_{n+1}(z)/rscale = 2*n/z * H_n(z)
c       H_{n+1}(z) = rscale*(2*n/z*H_n(z) - H_{n-1}(z)*rscale)
c
      scal2=rscale*rscale
      zinv=rscale/z
      do i=1,nterms-1
         dtmp=2*i
         ztmp=zinv*dtmp
         hvec(i+1)=ztmp*hvec(i)-scal2*hvec(i-1)
      enddo
c
c
c     From Abramowitz and Stegun (9.1.27)
c
c     H_{n}'(z)= H_{n-1}(z) - (n)/z * H_n(z)
c
c     With scaling:
c
c     hder(n)=scale* hvec(n-1) - n/z * hvec(n)
c
      if (ifder.eq.1) then
c
         hder(0)=-hvec(1)/rscale
         zinv=1.0d0/z
         do i=1,nterms
            dtmp=(i)
            ztmp=zinv*dtmp
            hder(i)=rscale*hvec(i-1)-ztmp*hvec(i)
         enddo
c
      endif
c
      return
      end
c
c
c
c
c
C***********************************************************************
        subroutine y2dformmp_dp(ier,dk,rscale,source,dipstr,dipvec,
     $     ns,center,nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a multipole (h) expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c               
c     mpole_n  =  sum dipstr_j  J_n(k r) exp(-i n theta_j) /rscale_j^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     dipstr(ns)      : dipole strengths
c     dipvec(2,ns)    : dipole orientation vectors
c     ns              : number of sources
c     center(2)       : epxansion center
c     nterms          : order of multipole expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the h-expansion
c
        real *8 dk
        real *8 dipstr(*)
        complex *16 zk,mpole(0:nterms)
        real *8 center(2),source(2,1),zdiff(2),dipvec(2,*)

        real *8, allocatable :: ival(:)
        complex *16, allocatable :: jtemp(:)
        complex *16, allocatable :: jval(:)
        complex *16, allocatable :: jder(:)
        integer, allocatable :: iscale(:)
        complex *16 zmul,zinv,ztemp1,ztemp2,ztemp3,ztemp4,ztemp5,ztemp6
c
        complex *16 ima,ima4,z,hexp1
        data ima/(0.0d0,1.0d0)/

        pi = 4.0d0*datan(1.0d0)
c
c
ccc        lwfjs = min(10000,nterms+5 + 4*nterms + 100)
        lwfjs = nterms+5 + 4*nterms + 100
        allocate(jval(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(ival(0:nterms+4), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jtemp(0:nterms+4), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jder(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(iscale(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
c
        do n=0,nterms
           mpole(n)=0
        enddo

        zk = ima*dk

        do j=1,ns
c
           zdiff(1)=source(1,j)-center(1)
           zdiff(2)=source(2,j)-center(2)
           call y2cart2polar(zdiff,r,theta)
           theta = theta-pi
           z=zk*r
           ifder=0
           call jfuns2d(ier,nterms+1,z,rscale,jval,ifder,jder,
     1           lwfjs,iscale,ntop)
c     convert Bessel J_n to modified Bessel I_n
           do n = 0,nterms+1,4
              ival(n) = dreal(jval(n))
              ival(n+1) = dimag(jval(n+1))
              ival(n+2) = -dreal(jval(n+2))
              ival(n+3) = -dimag(jval(n+3))
           enddo

           jtemp(0) = -ival(0)
           zmul=-exp(ima*theta)
           ztemp1= -zmul
           do i = 1,nterms+1
              jtemp( i) = ztemp1*ival(i)
              ztemp1= ztemp1*zmul
           enddo

           hexp1 = -dk*dipstr(j)*dipvec(1,j)/rscale 
     1          - ima*dk*dipstr(j)*dipvec(2,j)/rscale 

           rsj = rscale

           mpole(0)=mpole(0)+(dreal(hexp1)*dreal(jtemp(1))
     1          +dimag(hexp1)*dimag(jtemp(1)))*rsj**2
c
           hexp1 = hexp1

           do i=1,nterms
              mpole(i) = mpole(i)+(dreal(hexp1)*dreal(jtemp(i-1))
     1             -dimag(hexp1)*dimag(jtemp(i-1))
     2             +ima*(dreal(hexp1)*dimag(jtemp(i-1))
     3             +dimag(hexp1)*dreal(jtemp(i-1))))
              mpole(i) = mpole(i)+(dreal(hexp1)*dreal(jtemp(i+1))
     1             +dimag(hexp1)*dimag(jtemp(i+1))
     2             +ima*(dreal(hexp1)*dimag(jtemp(i+1))
     3             -dimag(hexp1)*dreal(jtemp(i+1))))*rsj**2
           enddo
        enddo
        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine y2dformmp_dp_add(ier,dk,rscale,source,dipstr,dipvec,
     $     ns,center,nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a multipole (h) expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c               
c     mpole_n  =  sum dipstr_j  J_n(k r) exp(-i n theta_j) /rscale_j^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     dipstr(ns)      : dipole strengths
c     dipvec(2,ns)    : dipole orientation vectors
c     ns              : number of sources
c     center(2)       : epxansion center
c     nterms          : order of multipole expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the h-expansion
c
       real *8 dk
        real *8 dipstr(*)
        complex *16 zk,mpole(0:nterms)
        real *8 center(2),source(2,1),zdiff(2),dipvec(2,*)

        real *8, allocatable :: ival(:)
        complex *16, allocatable :: jtemp(:)
        complex *16, allocatable :: jval(:)
        complex *16, allocatable :: jder(:)
        integer, allocatable :: iscale(:)
        complex *16 zmul,zinv,ztemp1,ztemp2,ztemp3,ztemp4,ztemp5,ztemp6
c
        complex *16 ima,ima4,z,hexp1
        data ima/(0.0d0,1.0d0)/

        pi = 4.0d0*datan(1.0d0)
c
c
ccc        lwfjs = min(10000,nterms+5 + 4*nterms + 100)
        lwfjs = nterms+5 + 4*nterms + 100
        allocate(jval(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(ival(0:nterms+4), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jtemp(0:nterms+4), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(jder(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
        allocate(iscale(0:lwfjs+10), stat=ier)
        if (ier.eq.1) then
          return
        endif
c
        zk = ima*dk

        do j=1,ns
c
           zdiff(1)=source(1,j)-center(1)
           zdiff(2)=source(2,j)-center(2)
           call y2cart2polar(zdiff,r,theta)
           theta = theta-pi
           z=zk*r
           ifder=0
           call jfuns2d(ier,nterms+1,z,rscale,jval,ifder,jder,
     1           lwfjs,iscale,ntop)
c     convert Bessel J_n to modified Bessel I_n
           do n = 0,nterms+1,4
              ival(n) = dreal(jval(n))
              ival(n+1) = dimag(jval(n+1))
              ival(n+2) = -dreal(jval(n+2))
              ival(n+3) = -dimag(jval(n+3))
           enddo

           jtemp(0) = -ival(0)
           zmul=-exp(ima*theta)
           ztemp1= -zmul
           do i = 1,nterms+1
              jtemp( i) = ztemp1*ival(i)
              ztemp1= ztemp1*zmul
           enddo

           hexp1 = -dk*dipstr(j)*dipvec(1,j)/rscale 
     1          - ima*dk*dipstr(j)*dipvec(2,j)/rscale 

           rsj = rscale

           mpole(0)=mpole(0)+(dreal(hexp1)*dreal(jtemp(1))
     1          +dimag(hexp1)*dimag(jtemp(1)))*rsj**2
c
           hexp1 = hexp1

           do i=1,nterms
              mpole(i) = mpole(i)+(dreal(hexp1)*dreal(jtemp(i-1))
     1             -dimag(hexp1)*dimag(jtemp(i-1))
     2             +ima*(dreal(hexp1)*dimag(jtemp(i-1))
     3             +dimag(hexp1)*dreal(jtemp(i-1))))
              mpole(i) = mpole(i)+(dreal(hexp1)*dreal(jtemp(i+1))
     1             +dimag(hexp1)*dimag(jtemp(i+1))
     2             +ima*(dreal(hexp1)*dimag(jtemp(i+1))
     3             -dimag(hexp1)*dreal(jtemp(i+1))))*rsj**2
           enddo
        enddo
        return
        end
c
c
c
C***********************************************************************
        subroutine y2dformta_dp(ier,dk,rscale,source,dipstr,dipvec,
     $     ns,center,nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a local (j) expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c               
c     mpole_n  =  sum dipstr_j  H_n(k r) exp(-i n theta_j) *rscale_j^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     dipstr(ns)      : dipole strengths
c     dipvec(2,ns)    : dipole orientation vectors
c     ns              : number of sources
c     center(2)       : epxansion center
c     nterms          : order of local expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the j-expansion
c
        real *8 dk
        complex *16 zk,mpole(0:nterms)
        real *8 dipstr(*)
        real *8 center(2),source(2,1),zdiff(2),dipvec(2,*)

        complex *16, allocatable :: hval(:)
        complex *16, allocatable :: hder(:)
        complex *16, allocatable :: htemp(:)
        real *8, allocatable :: kval(:)
        complex *16 zmul,zinv,ztemp1,ztemp2,ztemp3,ztemp4,ztemp5,ztemp6
        complex *16 hexp1
c
        complex *16 ima,ima4,z
        data ima/(0.0d0,1.0d0)/

        pi = 4.0d0*datan(1.0d0)
c
c
        allocate(hval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(hder(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(htemp(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(kval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif

c
        do n=0,nterms
           mpole(n)=0
        enddo

        zk = ima*dk

        do j=1,ns
           zdiff(1)=source(1,j)-center(1)
           zdiff(2)=source(2,j)-center(2)
           call y2cart2polar(zdiff,r,theta)
           theta = theta-pi
           z=zk*r
           ifder=0
           call y2dall(nterms+2,z,rscale,hval,ifder,hder)

c     convert Hankel H_n to modified Bessel K_n
c     (up to scale of pi/2)
           do n = 0,nterms+1,4
              kval(n) = -dimag(hval(n))
              kval(n+1) = -dreal(hval(n+1))
              kval(n+2) = dimag(hval(n+2))
              kval(n+3) = dreal(hval(n+3))
           enddo
           
           htemp(0) = kval(0)
           zmul=-exp(ima*theta)
           ztemp1= zmul
           do n = 1,nterms+1
              htemp( n) = ztemp1*kval(n)
              ztemp1= ztemp1*zmul
           enddo

           hexp1 = -dk*dipstr(j)*dipvec(1,j)/rscale 
     1          - ima*dk*dipstr(j)*dipvec(2,j)/rscale 

           mpole(0) = mpole(0)+(dreal(hexp1)*dreal(htemp(1))
     1          +dimag(hexp1)*dimag(htemp(1)))

           rsj2 = rscale**2

           do i = 1,nterms
              mpole(i) = mpole(i)+
     1             (dreal(hexp1)*dreal(htemp(i-1))
     1             -dimag(hexp1)*dimag(htemp(i-1))
     2             +ima*(dreal(hexp1)*dimag(htemp(i-1))
     3             +dimag(hexp1)*dreal(htemp(i-1))))*rsj2
              mpole(i) = mpole(i)+
     1             (dreal(hexp1)*dreal(htemp(i+1))
     1             +dimag(hexp1)*dimag(htemp(i+1))
     2             +ima*(dreal(hexp1)*dimag(htemp(i+1))
     3             -dimag(hexp1)*dreal(htemp(i+1))))
           enddo

        enddo
        return
        end
c
c
c
c
c
C***********************************************************************
        subroutine y2dformta_dp_add(ier,dk,rscale,source,dipstr,dipvec,
     $     ns,center,nterms,mpole)
        implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a local (j) expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c               
c     mpole_n  =  sum dipstr_j  H_n(k r) exp(-i n theta_j) *rscale_j^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     dipstr(ns)      : dipole strengths
c     dipvec(2,ns)    : dipole orientation vectors
c     ns              : number of sources
c     center(2)       : epxansion center
c     nterms          : order of local expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the j-expansion
c
        real *8 dk
        complex *16 zk,mpole(0:nterms)
        real *8 dipstr(*)
        real *8 center(2),source(2,1),zdiff(2),dipvec(2,*)

        complex *16, allocatable :: hval(:)
        complex *16, allocatable :: hder(:)
        complex *16, allocatable :: htemp(:)
        real *8, allocatable :: kval(:)
        complex *16 zmul,zinv,ztemp1,ztemp2,ztemp3,ztemp4,ztemp5,ztemp6
        complex *16 hexp1
c
        complex *16 ima,ima4,z
        data ima/(0.0d0,1.0d0)/

        pi = 4.0d0*datan(1.0d0)
c
c
        allocate(hval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(hder(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(htemp(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif
        allocate(kval(0:nterms+5), stat=ier)
        if (ier.eq.1) then
        return
        endif

c
        zk = ima*dk

        do j=1,ns
           zdiff(1)=source(1,j)-center(1)
           zdiff(2)=source(2,j)-center(2)
           call y2cart2polar(zdiff,r,theta)
           theta = theta-pi
           z=zk*r
           ifder=0
           call y2dall(nterms+2,z,rscale,hval,ifder,hder)

c     convert Hankel H_n to modified Bessel K_n
c     (up to scale of pi/2)
           do n = 0,nterms+1,4
              kval(n) = -dimag(hval(n))
              kval(n+1) = -dreal(hval(n+1))
              kval(n+2) = dimag(hval(n+2))
              kval(n+3) = dreal(hval(n+3))
           enddo
           
           htemp(0) = kval(0)
           zmul=-exp(ima*theta)
           ztemp1= zmul
           do n = 1,nterms
              htemp( n) = ztemp1*kval(n)
              ztemp1= ztemp1*zmul
           enddo

           hexp1 = -dk*dipstr(j)*dipvec(1,j)/rscale 
     1          - ima*dk*dipstr(j)*dipvec(2,j)/rscale 

           mpole(0) = mpole(0)+(dreal(hexp1)*dreal(htemp(1))
     1          +dimag(hexp1)*dimag(htemp(1)))

           rsj2 = rscale**2

           do i = 1,nterms
              mpole(i) = mpole(i)+
     1             (dreal(hexp1)*dreal(htemp(i-1))
     1             -dimag(hexp1)*dimag(htemp(i-1))
     2             +ima*(dreal(hexp1)*dimag(htemp(i-1))
     3             +dimag(hexp1)*dreal(htemp(i-1))))*rsj2
              mpole(i) = mpole(i)+
     1             (dreal(hexp1)*dreal(htemp(i+1))
     1             +dimag(hexp1)*dimag(htemp(i+1))
     2             +ima*(dreal(hexp1)*dimag(htemp(i+1))
     3             -dimag(hexp1)*dreal(htemp(i+1))))
           enddo

        enddo
        return
        end
c
c
c
c
c**********************************************************************
      subroutine ypotgrad2dall_dp(ifgrad,ifhess,sources,dipstr,dipvec,
     1                   ns,target,dk,pot,grad,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD, and
c     Hessian at the target point TARGET, due to a collection of 
c     dipoles at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot  = H_0(k*r)*(eye/4)
c		grad  = gradient(pot)
c		hess = Hessian
c---------------------------------------------------------------------
c     INPUT:
c
c     ifgrad         : flag for computing gradient
c	                 	   ifgrad = 0 -> don't compute 
c		                   ifgrad = 1 -> do compute 
c     ifhess        : flag for computing Hessian
c	                 	   ifhess = 0 -> don't compute 
c		                   ifhess = 1 -> do compute 
c     sources(2,*)  : location of the sources
c     dipstr        : dipole strengths
c     dipvec        : dipole vectors
c     ns            : number of sources
c     target        : location of the target
c     wavek         : Helmholtz parameter
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (complex *16)      : calculated potential
c     grad  (complex *16)      : calculated gradient
c     hess  (complex *16)      : calculated Hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 dk
      real *8 sources(2,ns),target(2),dipvec(2,ns)
      complex *16 wavek
      real *8 pot,grad(2),hess(3),potloc,gradloc(2),hessloc(3)
      complex *16 h0,h1,cd,eye,z
      real *8 dipstr(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      pot = 0.0d0
      if (ifgrad.eq.1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0
      endif
      if (ifhess.eq.1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0
      endif
c
      do i = 1,ns
         call ypotgrad2d_dp(ifgrad,ifhess,sources(1,i),dipstr(i),
     1        dipvec(1,i),target,dk,potloc,gradloc,hessloc)
         pot = pot + potloc
         if (ifgrad.eq.1) then
         grad(1) = grad(1) + gradloc(1)
         grad(2) = grad(2) + gradloc(2)
         endif
         if (ifhess.eq.1) then
         hess(1) = hess(1) + hessloc(1)
         hess(2) = hess(2) + hessloc(2)
         hess(3) = hess(3) + hessloc(3)
         endif
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine ypotgrad2d_dp(ifgrad,ifhess,source,dipstr,dipvec,
     $     target,dk,pot,grad,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine calculates the potential POT, gradient GRAD
c     and Hessian HESS at the target point TARGET, due to a dipole at 
c     SOURCE. The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot = H_0(k*r)*(eye/4)
c		grad = gradient(pot)
c		hess = Hessian
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     dipstr    : dipole strength
c     dipvec    : dipole orientation
c     target    : location of the target
c     wavek     : Helmholtz parameter
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      real *8 source(2),target(2)
      real *8 dk
      complex *16 wavek
      real *8 pot,grad(2),hess(3)
      real *8  dipstr
      real *8 dipvec(2)
      complex *16 z, h0, h1, h2, h3, h4, cd, h2z, zk, ima, ima4, ima4inv
      complex *16 hx,hy
      complex *16 hxx,hxy,hyy
      complex *16 hxxx,hxxy,hxyy,hyyy
c
      data ima/(0.0d0,1.0d0)/

      wavek = ima*dk
c
c ... Calculate offsets and distance
c
      ima4=-4*ima
      ima4inv=ima/4
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      rr=xdiff*xdiff+ydiff*ydiff
      r=sqrt(rr)
      z=wavek*r
c
      ifexpon = 1
      call hank103(z, h0, h1, ifexpon)
      ctheta=xdiff/r
      stheta=ydiff/r

      h2 = 2*h1/z-h0
      h3 = 4*h2/z-h1
c
      cd = h1/r*wavek*dipstr*ima4inv
      dotprod = xdiff*dipvec(1)+ydiff*dipvec(2)
      pot = dreal(cd*dotprod)
c
      if (ifgrad.eq.1) then
         cd = -wavek**2*dipstr*ima4inv
         hxx = cd*(h2*(ctheta*ctheta-0.5d0)-h0/2)
         hxy = cd*(h2*ctheta*stheta             )
         hyy = cd*(h2*(stheta*stheta-0.5d0)-h0/2)
         grad(1) = dreal(hxx*dipvec(1)+hxy*dipvec(2))
         grad(2) = dreal(hxy*dipvec(1)+hyy*dipvec(2))
      endif
c
      if (ifhess.eq.1) then
c         hess(1) = 0
c         hess(2) = 0
c         hess(3) = 0
         cd = -wavek**3*dipstr*ima4inv
         hxxx = cd*(
     $      -h1/2*(-1.5d0)
     $      -h3/2*(+ctheta*ctheta-0.5d0-stheta*stheta)
     $      )*ctheta
         hxxy = cd*(
     $      -h1/2*(-0.5d0)
     $      -h3/2*(1.5d0*ctheta*ctheta-0.5d0*stheta*stheta)
     $      )*stheta
         hxyy = cd*(
     $      -h1/2*(-0.5d0)
     $      -h3/2*(1.5d0*stheta*stheta-0.5d0*ctheta*ctheta)
     $      )*ctheta
         hyyy = cd*(
     $      -h1/2*(-1.5d0)
     $      -h3/2*(+stheta*stheta-0.5d0-ctheta*ctheta)
     $      )*stheta
         hess(1) = dreal(hxxx*dipvec(1)+hxxy*dipvec(2))
         hess(2) = dreal(hxxy*dipvec(1)+hxyy*dipvec(2))
         hess(3) = dreal(hxyy*dipvec(1)+hyyy*dipvec(2))
      endif
c
      return
      end
c
c
c
c**********************************************************************
c
c     ... direct calculation for a collection of charges and dipoles
c
c**********************************************************************
      subroutine ypotgrad2dall_sdp(dk,sources,ns,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     1     target,ifpot,pot,ifgrad,grad,ifhess,hess)
c**********************************************************************
c
c     This subroutine calculates the potential POT, field GRAD, and
c     Hessian at the target point TARGET, due to a collection of 
c     charges and dipoles at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c     
c              	pot  = H_0(k*r)*(eye/4)
c		grad  = gradient(pot)
c		hess = Hessian
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot         : flag for computing 
c	                 	   ifpot = 0 -> don't compute 
c		                   ifpot = 1 -> do compute 
c     ifgrad        : flag for computing gradient
c	                 	   ifgrad = 0 -> don't compute 
c		                   ifgrad = 1 -> do compute 
c     ifhess        : flag for computing Hessian
c	                 	   ifhess = 0 -> don't compute 
c		                   ifhess = 1 -> do compute 
c     sources(2,*)  : location of the sources
c     charge        : charge strengths
c     dipstr        : dipole strengths
c     dipvec        : dipole vectors
c     ns            : number of sources
c     target        : location of the target
c     wavek         : Helmholtz parameter
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot   (complex *16)      : calculated potential
c     grad  (complex *16)      : calculated gradient
c     hess  (complex *16)      : calculated Hessian
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 dk
      real *8 sources(2,ns),target(2),dipvec(2,ns)
      complex *16 wavek
      real *8 pot,grad(2),hess(3),potloc,gradloc(2),hessloc(3)
      complex *16 h0,h1,cd,eye,z
      real *8 dipstr(ns),charge(ns)
c
      data eye/(0.0d0,1.0d0)/
c
      if (ifpot.eq.1) pot = 0.0d0
      if (ifgrad.eq.1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0
      endif
      if (ifhess.eq.1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0
      endif
c
      do i = 1,ns
         call ypotgrad2d_sdp(dk,sources(1,i),
     1   ifcharge,charge(i),ifdipole,dipstr(i),dipvec(1,i),
     $   target,ifpot,potloc,ifgrad,gradloc,ifhess,hessloc)
         if (ifpot.eq.1) pot = pot + potloc
         if (ifgrad.eq.1) then
            grad(1) = grad(1) + gradloc(1)
            grad(2) = grad(2) + gradloc(2)
         endif
         if (ifhess.eq.1) then
            hess(1) = hess(1) + hessloc(1)
            hess(2) = hess(2) + hessloc(2)
            hess(3) = hess(3) + hessloc(3)
         endif
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine ypotgrad2d_sdp(dk,source,
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     target,ifpot,pot,ifgrad,grad,ifhess,hess)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     This subroutine calculates the potential POT, field GRAD and
c     Hessian HESS at the target point TARGET, due to a charge and a
c     dipole dipole at SOURCE. The scaling is that required of the delta
c     function response: i.e.,
c     
c              	pot = H_0(k*r)*(eye/4)
c		grad = gradient(pot)
c		hess = Hessian
c
c---------------------------------------------------------------------
c     INPUT:
c
c     ifpot      : flag for computing 
c	                	   ifpot = 0 -> don't compute 
c		                   ifpot = 1 -> do compute 
c     ifgrad     : flag for computing gradient
c	                 	ifgrad = 0 -> don't compute 
c		                ifgrad = 1 -> do compute 
c     ifhess     : flag for computing hessian
c	                 	ifhess = 0 -> don't compute 
c		                ifhess = 1 -> do compute 
c     source    : location of the source 
c     charge    : charge strength
c     dipstr    : dipole strength
c     dipvec    : dipole orientation
c     target    : location of the target
c     wavek     : Helmholtz parameter
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot       : calculated potential
c     grad      : calculated gradient
c     hess      : calculated hessian
c
c---------------------------------------------------------------------
      real *8 source(2),target(2)
      real *8 dk
      complex *16 wavek
      real *8 pot,grad(2),hess(3)
      real *8 charge,dipstr
      real *8 dipvec(2)
      complex *16 z, h0, h1, h2, h3, h4, cd, h2z, zk, ima, ima4, ima4inv
      complex *16 hx,hy
      complex *16 hxx,hxy,hyy
      complex *16 hxxx,hxxy,hxyy,hyyy
c
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      wavek = dk*ima
      
      ima4=-4*ima
      ima4inv=ima/4
c
      xdiff=target(1)-source(1)
      ydiff=target(2)-source(2)
      rr=xdiff*xdiff+ydiff*ydiff
      r=sqrt(rr)
      z=wavek*r
c
      ifexpon = 1
      call hank103(z, h0, h1, ifexpon)
c
c
      if (ifpot.eq.1) then
         pot = 0
      endif
c
      if (ifgrad.eq.1) then
         grad(1) = 0
         grad(2) = 0
      endif
c
      if (ifhess.eq.1) then
         hess(1) = 0
         hess(2) = 0
         hess(3) = 0
      endif
c
c
      if( ifcharge .eq. 1 ) then
c
      if( ifpot.eq.1) pot = dreal(h0*charge*ima4inv)
c
      if (ifgrad.eq.1) then
         cd = -h1*(wavek*charge*ima4inv/r)
         grad(1) = dreal(cd*xdiff)
         grad(2) = dreal(cd*ydiff)
      endif
c
      if (ifhess.eq.1) then
         cd = (wavek*charge*ima4inv/r)/rr
         h2z=(-z*h0+2*h1)
         hess(1) = dreal(cd*(h2z*xdiff*xdiff-rr*h1))
         hess(2) = dreal(cd*(h2z*xdiff*ydiff      ))
         hess(3) = dreal(cd*(h2z*ydiff*ydiff-rr*h1))
      endif
c
      endif
c
c
      if( ifdipole .eq. 1 ) then
      
      ctheta=xdiff/r
      stheta=ydiff/r

      h2 = 2*h1/z-h0
      h3 = 4*h2/z-h1
c
      cd = h1/r*wavek*dipstr*ima4inv
      dotprod = xdiff*dipvec(1)+ydiff*dipvec(2)

      if(ifpot.eq.1) pot = pot+dreal(cd*dotprod)
c
      if (ifgrad.eq.1) then
         cd = -wavek**2*dipstr*ima4inv
         hxx = cd*(h2*(ctheta*ctheta-0.5d0)-h0/2)
         hxy = cd*(h2*ctheta*stheta             )
         hyy = cd*(h2*(stheta*stheta-0.5d0)-h0/2)
         grad(1) = grad(1) + dreal(hxx*dipvec(1)+hxy*dipvec(2))
         grad(2) = grad(2) + dreal(hxy*dipvec(1)+hyy*dipvec(2))
      endif
c
      if (ifhess.eq.1) then
         cd = -wavek**3*dipstr*ima4inv
         hxxx = cd*(
     $      -h1/2*(-1.5d0)
     $      -h3/2*(+ctheta*ctheta-0.5d0-stheta*stheta)
     $      )*ctheta
         hxxy = cd*(
     $      -h1/2*(-0.5d0)
     $      -h3/2*(1.5d0*ctheta*ctheta-0.5d0*stheta*stheta)
     $      )*stheta
         hxyy = cd*(
     $      -h1/2*(-0.5d0)
     $      -h3/2*(1.5d0*stheta*stheta-0.5d0*ctheta*ctheta)
     $      )*ctheta
         hyyy = cd*(
     $      -h1/2*(-1.5d0)
     $      -h3/2*(+stheta*stheta-0.5d0-ctheta*ctheta)
     $      )*stheta
         hess(1) = hess(1) + dreal(hxxx*dipvec(1)+hxxy*dipvec(2))
         hess(2) = hess(2) + dreal(hxxy*dipvec(1)+hxyy*dipvec(2))
         hess(3) = hess(3) + dreal(hxyy*dipvec(1)+hyyy*dipvec(2))
      endif
c
      endif
c
      return
      end
c
c
c 
c**********************************************************************
      subroutine mkmptempyk(mptemp,theta,jval,ival,nterms)
c**********************************************************************
c
c     This subroutine computes the complex terms in a partial
c     wave expansion:
c                       I_n(k r) exp(i n theta)
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     theta   :    angular argument
c     jval    :    array of J_n values (up to nterms+3)
c     nterms  :    order of expansion is [-nterms-2:nterms+2]
c     ival    :    real *8 storage of length at least nterms+3
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     mptemp  :    I_n(kr)*exp(i n theta)
c
c
      implicit real *8 (a-h,o-z)
      complex *16 mptemp(0:nterms+2)
      complex *16 jval(0:nterms+3)
      real *8 ival(0:nterms+6)
      real *8 theta
      complex *16 ima,zfac,zfac2,zmull,zmullinv
      data ima/(0.0d0,1.0d0)/
c
c     convert Bessel J_n to modified Bessel I_n
      do n = 0,nterms+2,4
         ival(n) = dreal(jval(n))
         ival(n+1) = dimag(jval(n+1))
         ival(n+2) = -dreal(jval(n+2))
         ival(n+3) = -dimag(jval(n+3))
      enddo

      mptemp(0)=ival(0)
      zfac = exp(ima*theta)
      zmull = zfac

      do n=1,nterms+2
         mptemp(n)=ival(n)*zmull
         zmull = zmull*zfac
      enddo
      return
      end
c
c
c
       subroutine mkmpole12yk(mpolex,mpoley,ifhess,mpolexx,
     1               mpolexy,mpoleyy,mpole,dk,rscale,nterms)
c
c     This subroutine converts multipole expansion for
c     potential into expansions for various derivatives.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     ifhess  :    flag determines whether 2nd derivs are desired.
c     mpole   :    multipole expansion 
c     zk      :    Helmholtz parameter
c     rscale  :    scaling parameter for expansion
c     nterms  :    order of expansion
c                  for mpole [-nterms:nterms] 
c                  for mpolex,mpoley [-nterms-1:nterms+1] 
c                  for 2nd derivs [-nterms-2:nterms+2] 
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     mpolex,mpoley   : expansions for x,y derivatives
c     mpolexx,mpolexy,
c     mpoleyy         : expansions for 2nd derivatives if requested
c
c
c
       implicit real *8 (a-h,o-z)
       complex *16 mpolex(0:nterms+1)
       complex *16 mpoley(0:nterms+1)
       complex *16 mpolexx(0:nterms+2)
       complex *16 mpolexy(0:nterms+2)
       complex *16 mpoleyy(0:nterms+2)
       complex *16 mpole(0:nterms)
       complex *16 ima,zk
       real *8 dk
       data ima/(0.0d0,1.0d0)/
c
       do i=0,nterms+1
          mpolex(i)=0
          mpoley(i)=0
       enddo
       mpolex(1) = dreal(mpole(0))*dk*rscale
       mpoley(1) = dreal(mpole(0))*dk*rscale*ima
       do i=1,nterms
          mpolex(i-1)=mpolex(i-1)+dk/2/rscale*mpole(i)
          mpolex(i+1)=mpolex(i+1)+dk/2*rscale*mpole(i)
          mpoley(i-1)=mpoley(i-1)-dk/2*(ima)/rscale*mpole(i)
          mpoley(i+1)=mpoley(i+1)+dk/2*(ima)*rscale*mpole(i)
       enddo
c
       if (ifhess.eq.1) then
          do i=0,nterms+2
             mpolexx(i)=0
             mpolexy(i)=0
             mpoleyy(i)=0
          enddo
          mpolexx(1) = dreal(mpolex(0))*dk*rscale
          mpolexy(1) = dreal(mpolex(0))*dk*rscale*ima
          mpoleyy(1) = dreal(mpoley(0))*dk*rscale*ima

          do i=1,nterms+1
             mpolexx(i-1)=mpolexx(i-1)+dk/2/rscale*mpolex(i)
             mpolexx(i+1)=mpolexx(i+1)+dk/2*rscale*mpolex(i)
             mpolexy(i-1)=mpolexy(i-1)-dk/2*(ima)/rscale*mpolex(i)
             mpolexy(i+1)=mpolexy(i+1)+dk/2*(ima)*rscale*mpolex(i)
             mpoleyy(i-1)=mpoleyy(i-1)-dk/2*(ima)/rscale*mpoley(i)
             mpoleyy(i+1)=mpoleyy(i+1)+dk/2*(ima)*rscale*mpoley(i)
          enddo
       endif
       return
       end

      subroutine taevalsyk(mpole,mpolex,mpoley,mpolexx,mpolexy,mpoleyy,
     1                   rscale,nterms,jval,mptemp,pot,
     2                   ifgrad,grad,ifhess,hess,ima4inv)
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     mpole           : multipole expansion 
c     mpolex,mpoley   : expansions for x,y derivatives if requested
c     mpolexx,mpolexy,
c     mpoleyy         : expansions for 2nd derivatives if requested
c     rscale          :    scaling parameter for expansion
c     nterms          :    order of expansion
c                          for mpole [-nterms:nterms] 
c                          for mpolex,mpoley [-nterms-1:nterms+1] 
c                          for 2nd derivs [-nterms-2:nterms+2] 
c     jval            : Bessel expansion values
c     mptemp          : J_n exp(in theta) values for target from
c                       prior call to mkmptemp
c     ifgrad          : flag to request gradient computation
c     ifhess          : flag to request Hessian computation
c     ima4inv         : scaling factor (1/4i)
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot             : potential
c     grad            : gradient (if requested)  (pot_x,pot_y)
c     hess            : Hessian (if requested) (pot_xx,pot_xy,pot_yy)
c
c
      implicit real *8 (a-h,o-z)
      complex *16 mpole(0:nterms)
      complex *16 mptemp(0:nterms+2)
      complex *16 mpolex(0:nterms+1)
      complex *16 mpoley(0:nterms+1)
      complex *16 mpolexx(0:nterms+2)
      complex *16 mpolexy(0:nterms+2)
      complex *16 mpoleyy(0:nterms+2)
      complex *16 jval(0:nterms+3)
      real *8 pot,grad(2),hess(3)
      complex *16 ima4inv
c
      pot=dreal(mptemp(0))*dreal(mpole(0))
      do n=1,nterms
         pot=pot+dreal(mpole(n))*dreal(mptemp(n))
     1          +dimag(mpole(n))*dimag(mptemp(n))
      enddo
      pot=pot*0.25d0
c
      if (ifgrad .eq. 1) then
         grad(1)=0
         grad(2)=0
         grad(1)=dreal(mptemp(0))*dreal(mpolex(0))
         grad(2)=dreal(mptemp(0))*dreal(mpoley(0))
         do n=1,nterms+1
            grad(1)=grad(1)+dreal(mpolex(n))*dreal(mptemp(n))
     1           +dimag(mpolex(n))*dimag(mptemp(n))
            grad(2)=grad(2)+dreal(mpoley(n))*dreal(mptemp(n))
     1           +dimag(mpoley(n))*dimag(mptemp(n))
         enddo
         grad(1)=grad(1)*0.25d0
         grad(2)=grad(2)*0.25d0
      endif
c
      if (ifhess .eq. 1) then
         hess(1)=0
         hess(2)=0
         hess(3)=0
         hess(1)=dreal(mptemp(0))*dreal(mpolexx(0))
         hess(2)=dreal(mptemp(0))*dreal(mpolexy(0))
         hess(3)=dreal(mptemp(0))*dreal(mpoleyy(0))
         do n=1,nterms+2
            hess(1)=hess(1)+dreal(mpolexx(n))*dreal(mptemp(n))
     1           +dimag(mpolexx(n))*dimag(mptemp(n))
            hess(2)=hess(2)+dreal(mpolexy(n))*dreal(mptemp(n))
     1           +dimag(mpolexy(n))*dimag(mptemp(n))
            hess(3)=hess(3)+dreal(mpoleyy(n))*dreal(mptemp(n))
     1           +dimag(mpoleyy(n))*dimag(mptemp(n))
         enddo
         hess(1)=hess(1)*0.25d0
         hess(2)=hess(2)*0.25d0
         hess(3)=hess(3)*0.25d0
      endif
      return
      end
c
c
c

c**********************************************************************
      subroutine y2dmpeval_mult(dk,rscale,center,mpole,nterms,ztarg,
     1      pot,ifgrad,grad,ifhess,hess,nmpole)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     WARNING: for speed, the maximum NTERMS allowed is hardcoded
c              (NTERMS must be less than 500)
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to several outgoing partial wave expansions
c     (from the same source point)
c
c          
c     pot_i=   (1/(2pi)) K_0(kr) mpole_0,i + 
c                  
c                  nterms
c       (1/(2pi))   sum  K_n(k r)rscale^n [real(mpole_n,i)cos(n theta) 
c                   n=1                 + imag(mpole_n,i)sin(n theta)]
c                                                                      
c                                   
c     grad  = gradient(pot) if ifgrad = 1.
c     hess = hessian if ifhess = 1.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     dk          :    Yukawa parameter (real)
c     rscale      :    scaling parameter 
c     center      :    expansion center
c     mpole(:,i)  :    i-th multipole expansion 
c     nterms      :    order of the multipole expansion
c     ztarg       :    target location
c     ifgrad      :    flag controlling evaluation of gradient:
c                        ifgrad = 0, do not compute gradient.
c                        ifgrad = 1, compute gradient.
c     ifhess      :    flag for computing Hessian:
c	                 ifhess = 0 -> don't compute 
c		         ifhess = 1 -> do compute 
c     nmpole      :    number of multipole expansions to evaluate
c                      at target point
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot(i)        :    potential at ztarg due to i-th expansion
c     grad(1:2,i)   :    gradient due to i-th expansion (if requested)
c     hess(1:3,i)   :    hessian due to i-th expansion (if requested)
c
c-----------------------------------------------------------------------
c
        complex *16 zk,mpole(0:nterms,nmpole),zmul,ztemp1
        real *8 pot(nmpole), grad(2,nmpole), hess(3,nmpole)
        real *8 center(2),ztarg(2),zdiff(2)
        real *8 dk, pi, dkor, dktr
c
        complex *16 hval(0:500)
        complex *16 hder(0:500)

        real *8 kval(0:500)
c
        complex *16 mpolex(0:500)
        complex *16 mpoley(0:500)
c
        complex *16 mpolexx(0:500)
        complex *16 mpolexy(0:500)
        complex *16 mpoleyy(0:500)
c
        complex *16 mptemp(0:500)
c
        complex *16 ima,ima4,z,pieye2
        data ima/(0.0d0,1.0d0)/
c
c
        ima4=-4*ima
c        pi = 4.0d0*datan(1.0d0)
c        pieye2=pi*ima/2.0d0
c
c
        zdiff(1)=ztarg(1)-center(1)
        zdiff(2)=ztarg(2)-center(2)
        call y2cart2polar(zdiff,r,theta)
c
        zk = ima*dk
        z=zk*r
        ifder=0
        call y2dall(nterms+3,z,rscale,hval,ifder,hder)
c
c

c     convert Hankel H_n to modified Bessel K_n
c     (up to scale of pi/2)
        do n = 0,nterms+2,4
           kval(n) = -dimag(hval(n))
           kval(n+1) = -dreal(hval(n+1))
           kval(n+2) = dimag(hval(n+2))
           kval(n+3) = dreal(hval(n+3))
        enddo

c     at this point kval(n) = 2/pi * K_n(dk*r)

        mptemp(0)=kval(0)
        zmul=exp(ima*theta)
        ztemp1=zmul
        do n=1,nterms+2
           mptemp(n)=dcmplx(kval(n)*dreal(ztemp1),
     1          kval(n)*dimag(ztemp1))
           ztemp1 = ztemp1*zmul
        enddo
        
        dkor = dk/rscale
        dktr = dk*rscale
        dkorh = 0.5d0*dkor
        dktrh = 0.5d0*dktr

        do iii = 1,nmpole
           

           if( ifgrad .eq. 1 .or. ifhess .eq. 1 ) then
c     
c     
              do i=0,nterms+1
                 mpolex(i)=0
                 mpoley(i)=0
              enddo

              mpolex(1) = -dkor*mpole(0,iii)
              mpoley(1) = -dkor*mpole(0,iii)*ima

              do i=1,nterms
                 mpolex(i-1) = mpolex(i-1) -dktrh*mpole(i,iii)
                 mpolex(i+1) = mpolex(i+1) -dkorh*mpole(i,iii)
                 mpoley(i-1) = mpoley(i-1) +dktrh*mpole(i,iii)*ima
                 mpoley(i+1) = mpoley(i+1) -dkorh*mpole(i,iii)*ima
              enddo
           endif
c     
           if( ifhess .eq. 1 ) then
c     
              do i=0,nterms+2
                 mpolexx(i)=0
                 mpolexy(i)=0
                 mpoleyy(i)=0
              enddo

              mpolexx(1) = -dk/1.0d0/rscale*dreal(mpolex(0))
              mpolexy(1) = -dk/1.0d0/rscale*dreal(mpolex(0))*ima
              mpoleyy(1) = -dk/1.0d0/rscale*dreal(mpoley(0))*ima

              do i=1,nterms+1
                 mpolexx(i-1) = mpolexx(i-1) -dk/2.0d0*mpolex(i)*rscale
                 mpolexx(i+1) = mpolexx(i+1) -dk/2.0d0*mpolex(i)/rscale
                 mpolexy(i-1) = mpolexy(i-1) +dk/2.0d0*mpolex(i)*rscale
     1                *ima
                 mpolexy(i+1) = mpolexy(i+1) -dk/2.0d0*mpolex(i)/rscale
     1                *ima
                 mpoleyy(i-1) = mpoleyy(i-1) +dk/2.0d0*mpoley(i)*rscale
     1                *ima
                 mpoleyy(i+1) = mpoleyy(i+1) -dk/2.0d0*mpoley(i)/rscale
     1                *ima
              enddo


           endif
c     
c     
           pot(iii)=kval(0)*dreal(mpole(0,iii))
           do n=1,nterms
              pot(iii)=pot(iii)+dreal(mpole(n,iii))*dreal(mptemp(n))
     1             +dimag(mpole(n,iii))*dimag(mptemp(n))
           enddo
           pot(iii)=pot(iii)*0.25d0
c     

           if( ifgrad .eq. 1 ) then
              grad(1,iii)=0
              grad(2,iii)=0
c     
              grad(1,iii)=kval(0)*dreal(mpolex(0))
              grad(2,iii)=kval(0)*dreal(mpoley(0))

              do n=1,nterms+1
                 grad(1,iii)=grad(1,iii)+dreal(mpolex(n))*
     1                dreal(mptemp(n))+dimag(mpolex(n))*dimag(mptemp(n))
                 grad(2,iii)=grad(2,iii)+dreal(mpoley(n))*
     1                dreal(mptemp(n))+dimag(mpoley(n))*dimag(mptemp(n))
              enddo
              grad(1,iii)=grad(1,iii)*.25d0
              grad(2,iii)=grad(2,iii)*.25d0

           endif
c     
c     
           if( ifhess .eq. 1 ) then
              hess(1,iii)=0
              hess(2,iii)=0
              hess(3,iii)=0
c     
              hess(1,iii)=kval(0) *dreal(mpolexx(0))
              hess(2,iii)=kval(0) *dreal(mpolexy(0))
              hess(3,iii)=kval(0) *dreal(mpoleyy(0))
              do n=1,nterms+2
                 hess(1,iii)=hess(1,iii)+
     1                dreal(mpolexx(n))*dreal(mptemp(n))
     2                +dimag(mpolexx(n))*dimag(mptemp(n))
                 hess(2,iii)=hess(2,iii)+
     1                dreal(mpolexy(n))*dreal(mptemp(n))
     2                +dimag(mpolexy(n))*dimag(mptemp(n))
                 hess(3,iii)=hess(3,iii)+
     1                dreal(mpoleyy(n))*dreal(mptemp(n))
     2                +dimag(mpoleyy(n))*dimag(mptemp(n))
              enddo
              hess(1,iii)=hess(1,iii)*.25d0
              hess(2,iii)=hess(2,iii)*.25d0
              hess(3,iii)=hess(3,iii)*.25d0
           endif

        enddo
c     
        return
        end
c
c

C***********************************************************************
      subroutine y2dformmp_qp(ier,dk,rscale,source,quadstr,quadvec,
     $     ns,center,nterms,mpole)
      implicit none
C***********************************************************************
c     
c     This subroutine constructs a multipole (h) expansion about 
c     CENTER due to NS quadrupole sources located at SOURCES(2,*).
c     
c     When evaluated at the point p (assumed to be well separated), 
c     the multipole expansion approximates the function
c
c     sum quadstr(j)*( quadvec(1,j)* (K_0)_xx (k ||p-source_j||)
c                    + quadvec(2,j)* (K_0)_xy (k ||p-source_j||)
c                    + quadvec(3,j)* (K_0)_yy (k ||p-source_j||))
c     
c-----------------------------------------------------------------------
c     INPUT:
c     
c     k              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     quadstr(ns)     : quadrupole strengths
c     quadvec(3,ns)   : quadrupole ``orientation''
c                     quadvec(1,i) - (H_0)_xx term
c                     quadvec(2,i) - (H_0)_xy term
c                     quadvec(3,i) - (H_0)_yy term
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c     
c     OUTPUT:
c     
c     ier       : error return code
c     ier=0 returned successfully;
c     
c     mpole     : coeffs for the h-expansion
c     
      
c     global variables
      complex *16 quadstr(*),mpole(0:nterms)
      real *8 dk
      integer ns, nterms, ier
      real *8 center(2),source(2,*),quadvec(3,*), rscale
c     local variables
      complex *16 z,ima, zmul,zinv,ztemp1,ztemp2,hexp(-2:2)
      complex *16 zk
      real *8 zdiff(2), theta, r, rs2, rs4, pi, done, dk2
      integer iii, i, j, ntop, lwfjs, ifder, n
      complex *16, allocatable :: jval(:), jder(:), jtemp(:)
      integer, allocatable :: iscale(:)
      real *8, allocatable :: ival(:)
c     
      data ima/(0.0d0,1.0d0)/
c     

      dk2 = dk*dk
      zk = ima*dk

      ier = 0

      done=1
      pi=4*atan(done)
c     
      lwfjs = nterms+5 + 4*nterms + 100
      allocate(jval(0:lwfjs+10), stat=ier)
      if (ier.eq.1) then
         return
      endif
      allocate(ival(0:nterms+10), stat=ier)
      if (ier.eq.1) then
         return
      endif
      allocate(jder(0:lwfjs+10), stat=ier)
      if (ier.eq.1) then
         return
      endif
      allocate(iscale(0:lwfjs+10), stat=ier)
      if (ier.eq.1) then
         return
      endif
      allocate(jtemp(0:nterms+5), stat=ier)
      if (ier.eq.1) then
         return
      endif
c     
      do i = 0,nterms
         mpole(i) = 0
      enddo

      do iii = 1,ns

c     get j vals for this source

         zdiff(1)=center(1)-source(1,iii)
         zdiff(2)=center(2)-source(2,iii)
         call h2cart2polar(zdiff,r,theta)
         theta=theta-pi
         z=zk*r
         ifder=0
         call jfuns2d(ier,nterms+2,z,rscale,jval,ifder,jder,
     1        lwfjs,iscale,ntop)

c     convert Bessel J_n to modified Bessel I_n
         do n = 0,nterms+2,4
            ival(n) = dreal(jval(n))
            ival(n+1) = dimag(jval(n+1))
            ival(n+2) = -dreal(jval(n+2))
            ival(n+3) = -dimag(jval(n+3))
         enddo
c     
         jtemp(0) = ival(0)
         zmul=exp(ima*theta)
         ztemp1= zmul
         do j = 1,nterms+2
            jtemp( j) = ztemp1*ival(j)
            ztemp1= ztemp1*zmul
         enddo

c     set-up quadrupole as if centered at the source

         rs2 = rscale*rscale 

         hexp(0) = quadstr(iii)*(dk2*quadvec(1,iii)
     1                           +dk2*quadvec(3,iii))/2.0d0
         hexp(1) = 0.0d0
         hexp(2) = quadstr(iii)*(dk2*quadvec(1,iii)
     1                          +dk2*quadvec(2,iii)*ima
     2                          -dk2*quadvec(3,iii))/(rs2*2.0d0)


c     shift the quadrupole

         mpole(0) = mpole(0) + hexp(0)*jtemp(0)

         rs4 = rs2*rs2

         mpole(0) = mpole(0)+(dreal(hexp(2))*dreal(jtemp(2)) 
     1        +dimag(hexp(2))*dimag(jtemp(2)))*rs4


         do i = 1,nterms

c     due to zero-th term

            mpole(i) = mpole(i) + 2.0d0*(dreal(hexp(0))*dreal(jtemp(i))
     1           +ima*dreal(hexp(0))*dimag(jtemp(i)))

c     due to 2nd term

            if (i .ge. 2) then
               mpole(i) = mpole(i)+(dreal(hexp(2))*dreal(jtemp(i-2))
     1              -dimag(hexp(2))*dimag(jtemp(i-2))
     2              +ima*(dreal(hexp(2))*dimag(jtemp(i-2))
     3              +dimag(hexp(2))*dreal(jtemp(i-2))))
               mpole(i) = mpole(i)+(dreal(hexp(2))*dreal(jtemp(i+2))
     1              +dimag(hexp(2))*dimag(jtemp(i+2))
     2              +ima*(dreal(hexp(2))*dimag(jtemp(i+2))
     3              -dimag(hexp(2))*dreal(jtemp(i+2))))*rs4
            else
c     i = 1 here
               mpole(i) = mpole(i)+(dreal(hexp(2))*dreal(jtemp(2-i))
     1              +dimag(hexp(2))*dimag(jtemp(2-i))
     2              -ima*(dreal(hexp(2))*dimag(jtemp(2-i))
     3              -dimag(hexp(2))*dreal(jtemp(2-i))))*rs2
               mpole(i) = mpole(i)+(dreal(hexp(2))*dreal(jtemp(i+2))
     1              +dimag(hexp(2))*dimag(jtemp(i+2))
     2              +ima*(dreal(hexp(2))*dimag(jtemp(i+2))
     3              -dimag(hexp(2))*dreal(jtemp(i+2))))*rs4
            endif
            
         enddo
         
      enddo

      return
      end
c
c
C***********************************************************************
      subroutine y2dformmp_qp_add(ier,dk,rscale,source,quadstr,quadvec,
     $     ns,center,nterms,mpole)
      implicit none
C***********************************************************************
c     
c     This subroutine constructs a multipole (h) expansion about 
c     CENTER due to NS quadrupole sources located at SOURCES(2,*).
c     
c     When evaluated at the point p (assumed to be well separated), 
c     the multipole expansion approximates the function
c
c     sum quadstr(j)*( quadvec(1,j)* (K_0)_xx (k ||p-source_j||)
c                    + quadvec(2,j)* (K_0)_xy (k ||p-source_j||)
c                    + quadvec(3,j)* (K_0)_yy (k ||p-source_j||))
c     
c-----------------------------------------------------------------------
c     INPUT:
c     
c     k              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     quadstr(ns)     : quadrupole strengths
c     quadvec(3,ns)   : quadrupole ``orientation''
c                     quadvec(1,i) - (H_0)_xx term
c                     quadvec(2,i) - (H_0)_xy term
c                     quadvec(3,i) - (H_0)_yy term
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c     
c     OUTPUT:
c     
c     ier       : error return code
c     ier=0 returned successfully;
c     
c     mpole     : coeffs for the h-expansion
c     
      
c     global variables
      complex *16 quadstr(*),mpole(0:nterms)
      real *8 dk
      integer ns, nterms, ier
      real *8 center(2),source(2,*),quadvec(3,*), rscale
c     local variables
      complex *16 z,ima, zmul,zinv,ztemp1,ztemp2,hexp(-2:2)
      complex *16 zk
      real *8 zdiff(2), theta, r, rs2, rs4, pi, done, dk2
      integer iii, i, j, ntop, lwfjs, ifder, n
      complex *16, allocatable :: jval(:), jder(:), jtemp(:)
      integer, allocatable :: iscale(:)
      real *8, allocatable :: ival(:)
c     
      data ima/(0.0d0,1.0d0)/
c     

      zk = ima*dk

      dk2 = dk*dk

      ier = 0

      done=1
      pi=4*atan(done)
c     
      lwfjs = nterms+5 + 4*nterms + 100
      allocate(jval(0:lwfjs+10), stat=ier)
      if (ier.eq.1) then
         return
      endif
      allocate(ival(0:nterms+10), stat=ier)
      if (ier.eq.1) then
         return
      endif
      allocate(jder(0:lwfjs+10), stat=ier)
      if (ier.eq.1) then
         return
      endif
      allocate(iscale(0:lwfjs+10), stat=ier)
      if (ier.eq.1) then
         return
      endif
      allocate(jtemp(0:nterms+5), stat=ier)
      if (ier.eq.1) then
         return
      endif
c     


      do iii = 1,ns

c     get j vals for this source

         zdiff(1)=center(1)-source(1,iii)
         zdiff(2)=center(2)-source(2,iii)
         call h2cart2polar(zdiff,r,theta)
         theta=theta-pi
         z=zk*r
         ifder=0
         call jfuns2d(ier,nterms+2,z,rscale,jval,ifder,jder,
     1        lwfjs,iscale,ntop)

c     convert Bessel J_n to modified Bessel I_n
         do n = 0,nterms+2,4
            ival(n) = dreal(jval(n))
            ival(n+1) = dimag(jval(n+1))
            ival(n+2) = -dreal(jval(n+2))
            ival(n+3) = -dimag(jval(n+3))
         enddo
c     
         jtemp(0) = ival(0)
         zmul=exp(ima*theta)
         ztemp1= zmul
         do j = 1,nterms+2
            jtemp( j) = ztemp1*ival(j)
            ztemp1= ztemp1*zmul
         enddo

c     set-up quadrupole as if centered at the source

         rs2 = rscale*rscale 

         hexp(0) = quadstr(iii)*(dk2*quadvec(1,iii)
     1                           +dk2*quadvec(3,iii))/2.0d0
         hexp(1) = 0.0d0
         hexp(2) = quadstr(iii)*(dk2*quadvec(1,iii)
     1                          +dk2*quadvec(2,iii)*ima
     2                          -dk2*quadvec(3,iii))/(rs2*2.0d0)


c     shift the quadrupole

         mpole(0) = mpole(0) + hexp(0)*jtemp(0)

         rs4 = rs2*rs2

         mpole(0) = mpole(0)+(dreal(hexp(2))*dreal(jtemp(2)) 
     1        +dimag(hexp(2))*dimag(jtemp(2)))*rs4


         do i = 1,nterms

c     due to zero-th term

            mpole(i) = mpole(i) + 2.0d0*(dreal(hexp(0))*dreal(jtemp(i))
     1           +ima*dreal(hexp(0))*dimag(jtemp(i)))

c     due to 2nd term

            if (i .ge. 2) then
               mpole(i) = mpole(i)+(dreal(hexp(2))*dreal(jtemp(i-2))
     1              -dimag(hexp(2))*dimag(jtemp(i-2))
     2              +ima*(dreal(hexp(2))*dimag(jtemp(i-2))
     3              +dimag(hexp(2))*dreal(jtemp(i-2))))
               mpole(i) = mpole(i)+(dreal(hexp(2))*dreal(jtemp(i+2))
     1              +dimag(hexp(2))*dimag(jtemp(i+2))
     2              +ima*(dreal(hexp(2))*dimag(jtemp(i+2))
     3              -dimag(hexp(2))*dreal(jtemp(i+2))))*rs4
            else
c     i = 1 here
               mpole(i) = mpole(i)+(dreal(hexp(2))*dreal(jtemp(2-i))
     1              +dimag(hexp(2))*dimag(jtemp(2-i))
     2              -ima*(dreal(hexp(2))*dimag(jtemp(2-i))
     3              -dimag(hexp(2))*dreal(jtemp(2-i))))*rs2
               mpole(i) = mpole(i)+(dreal(hexp(2))*dreal(jtemp(i+2))
     1              +dimag(hexp(2))*dimag(jtemp(i+2))
     2              +ima*(dreal(hexp(2))*dimag(jtemp(i+2))
     3              -dimag(hexp(2))*dreal(jtemp(i+2))))*rs4
            endif
            
         enddo
         
      enddo

      return
      end
c
c
c
C***********************************************************************
      subroutine y2dformta_qp(ier,dk,rscale,source,quadstr,quadvec,
     $     ns,center,nterms,mpole)
      implicit none
C***********************************************************************
c
c     This subroutine constructs a local (j) expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c               
c     mpole_n  =  sum quadstr_j  H_n(k r) exp(-i n theta_j) *rscale_j^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     dk              : Yukawa parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     quadstr(ns)      : dipole strengths
c     quadvec(2,ns)    : dipole orientation vectors
c     ns              : number of sources
c     center(2)       : epxansion center
c     nterms          : order of local expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the j-expansion
c
      
c     global variables
      complex *16 mpole(0:nterms),quadstr(*)
      real *8 center(2),source(2,*),quadvec(3,*)
      real *8 rscale, dk
      integer ns, nterms, ier
c     local variables
      complex *16 hexp(0:2)
      complex *16, allocatable :: hval(:), htemp(:)
      complex *16, allocatable :: hder(:)
      real *8, allocatable :: kval(:)
      complex *16 zmul,zinv,ztemp1,ztemp2,zk
      real *8 theta, r, rs2, rs4, zdiff(2), dk2, pi
      integer i, iii, ifder, n
c     
      complex *16 ima,ima4,z
      data ima/(0.0d0,1.0d0)/

      pi = 4.0d0*datan(1.0d0)

      dk2 = dk*dk

      zk = ima*dk
c     
c     
      allocate(hval(0:nterms+5), stat=ier)
      if (ier.eq.1) then
         return
      endif
      allocate(kval(0:nterms+10), stat=ier)
      if (ier.eq.1) then
         return
      endif
      allocate(hder(0:nterms+5), stat=ier)
      if (ier.eq.1) then
         return
      endif
      allocate(htemp(0:nterms+5), stat=ier)
      if (ier.eq.1) then
         return
      endif

c     
      do i=0,nterms
         mpole(i)=0
      enddo

      do iii=1,ns

c     get h vals for this source

         zdiff(1)=source(1,iii)-center(1)
         zdiff(2)=source(2,iii)-center(2)
         call h2cart2polar(zdiff,r,theta)
         z=zk*r
         ifder=0
         call h2dall(nterms+2,z,rscale,hval,ifder,hder)

c     convert Hankel H_n to modified Bessel K_n
c     (up to scale of pi/2)
         do n = 0,nterms+2,4
            kval(n) = -dimag(hval(n))
            kval(n+1) = -dreal(hval(n+1))
            kval(n+2) = dimag(hval(n+2))
            kval(n+3) = dreal(hval(n+3))
         enddo


         htemp(0) = kval(0)
         zmul=exp(ima*theta)
         ztemp1= zmul
         do i = 1,nterms+2
            htemp( i) = ztemp1*kval(i)
            ztemp1= ztemp1*zmul
         enddo

c     set-up quadrupole as if centered at the source

         rs2 = rscale*rscale 

         hexp(0) = quadstr(iii)*(dk2*quadvec(1,iii)
     1                           +dk2*quadvec(3,iii))/2.0d0
         hexp(1) = 0.0d0
         hexp(2) = quadstr(iii)*(dk2*quadvec(1,iii)
     1                          +dk2*quadvec(2,iii)*ima
     2                          -dk2*quadvec(3,iii))/(rs2*2.0d0)

c     shift the quadrupole

         mpole(0) = mpole(0) + hexp(0)*htemp(0)

         mpole(0) = mpole(0)+(dreal(hexp(2))*dreal(htemp(2))
     1        +dimag(hexp(2))*dimag(htemp(2)))
c     
         rs2 = rscale*rscale
         rs4 = rs2*rs2
         do i = 1,nterms
            mpole(i) = mpole(i) + 2.0d0*(dreal(hexp(0))*dreal(htemp(i))
     1           +ima*dreal(hexp(0))*dimag(htemp(i)))
            if (i .ge. 2) then
               mpole(i) = mpole(i)+(dreal(hexp(2))*dreal(htemp(i-2))
     1              -dimag(hexp(2))*dimag(htemp(i-2)) 
     2              +ima*(dreal(hexp(2))*dimag(htemp(i-2)) 
     3              +dimag(hexp(2))*dreal(htemp(i-2))))*rs4
               mpole(i) = mpole(i)+(dreal(hexp(2))*dreal(htemp(i+2))
     1              +dimag(hexp(2))*dimag(htemp(i+2))
     2              +ima*(dreal(hexp(2))*dimag(htemp(i+2))
     3              -dimag(hexp(2))*dreal(htemp(i+2))))
            else
c     here i = 1
               mpole(i) = mpole(i)+(dreal(hexp(2))*dreal(htemp(2-i))
     1              +dimag(hexp(2))*dimag(htemp(2-i))
     2              -ima*(dreal(hexp(2))*dimag(htemp(2-i))
     3              -dimag(hexp(2))*dreal(htemp(2-i))))*rs2
               mpole(i) = mpole(i)+(dreal(hexp(2))*dreal(htemp(i+2))
     1              +dimag(hexp(2))*dimag(htemp(i+2))
     2              +ima*(dreal(hexp(2))*dimag(htemp(i+2))
     3              -dimag(hexp(2))*dreal(htemp(i+2))))

            endif
         enddo
      enddo

      return
      end
c     
c     
c     
c
c
C***********************************************************************
      subroutine y2dformta_qp_add(ier,dk,rscale,source,quadstr,quadvec,
     $     ns,center,nterms,mpole)
      implicit none
C***********************************************************************
c
c     This subroutine constructs a local (j) expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c               
c     mpole_n  =  sum quadstr_j  H_n(k r) exp(-i n theta_j) *rscale_j^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     dk              : Yukawa parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     quadstr(ns)      : dipole strengths
c     quadvec(2,ns)    : dipole orientation vectors
c     ns              : number of sources
c     center(2)       : epxansion center
c     nterms          : order of local expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the j-expansion
c
      
c     global variables
      complex *16 mpole(0:nterms),quadstr(*)
      real *8 center(2),source(2,*),quadvec(3,*)
      real *8 rscale, dk
      integer ns, nterms, ier
c     local variables
      complex *16 hexp(0:2)
      complex *16, allocatable :: hval(:), htemp(:)
      complex *16, allocatable :: hder(:)
      real *8, allocatable :: kval(:)
      complex *16 zmul,zinv,ztemp1,ztemp2,zk
      real *8 theta, r, rs2, rs4, zdiff(2), dk2, pi
      integer i, iii, ifder, n
c     
      complex *16 ima,ima4,z
      data ima/(0.0d0,1.0d0)/

      pi = 4.0d0*datan(1.0d0)

      dk2 = dk*dk

      zk = ima*dk
c     
c     
      allocate(hval(0:nterms+5), stat=ier)
      if (ier.eq.1) then
         return
      endif
      allocate(kval(0:nterms+10), stat=ier)
      if (ier.eq.1) then
         return
      endif
      allocate(hder(0:nterms+5), stat=ier)
      if (ier.eq.1) then
         return
      endif
      allocate(htemp(0:nterms+5), stat=ier)
      if (ier.eq.1) then
         return
      endif

c     
      do iii=1,ns

c     get h vals for this source

         zdiff(1)=source(1,iii)-center(1)
         zdiff(2)=source(2,iii)-center(2)
         call h2cart2polar(zdiff,r,theta)
         z=zk*r
         ifder=0
         call h2dall(nterms+2,z,rscale,hval,ifder,hder)

c     convert Hankel H_n to modified Bessel K_n
c     (up to scale of pi/2)
         do n = 0,nterms+2,4
            kval(n) = -dimag(hval(n))
            kval(n+1) = -dreal(hval(n+1))
            kval(n+2) = dimag(hval(n+2))
            kval(n+3) = dreal(hval(n+3))
         enddo


         htemp(0) = kval(0)
         zmul=exp(ima*theta)
         ztemp1= zmul
         do i = 1,nterms+2
            htemp( i) = ztemp1*kval(i)
            ztemp1= ztemp1*zmul
         enddo

c     set-up quadrupole as if centered at the source

         rs2 = rscale*rscale 

         hexp(0) = quadstr(iii)*(dk2*quadvec(1,iii)
     1                           +dk2*quadvec(3,iii))/2.0d0
         hexp(1) = 0.0d0
         hexp(2) = quadstr(iii)*(dk2*quadvec(1,iii)
     1                          +dk2*quadvec(2,iii)*ima
     2                          -dk2*quadvec(3,iii))/(rs2*2.0d0)

c     shift the quadrupole

         mpole(0) = mpole(0) + hexp(0)*htemp(0)

         mpole(0) = mpole(0)+(dreal(hexp(2))*dreal(htemp(2))
     1        +dimag(hexp(2))*dimag(htemp(2)))
c     
         rs2 = rscale*rscale
         rs4 = rs2*rs2
         do i = 1,nterms
            mpole(i) = mpole(i) + 2.0d0*(dreal(hexp(0))*dreal(htemp(i))
     1           +ima*dreal(hexp(0))*dimag(htemp(i)))
            if (i .ge. 2) then
               mpole(i) = mpole(i)+(dreal(hexp(2))*dreal(htemp(i-2))
     1              -dimag(hexp(2))*dimag(htemp(i-2)) 
     2              +ima*(dreal(hexp(2))*dimag(htemp(i-2)) 
     3              +dimag(hexp(2))*dreal(htemp(i-2))))*rs4
               mpole(i) = mpole(i)+(dreal(hexp(2))*dreal(htemp(i+2))
     1              +dimag(hexp(2))*dimag(htemp(i+2))
     2              +ima*(dreal(hexp(2))*dimag(htemp(i+2))
     3              -dimag(hexp(2))*dreal(htemp(i+2))))
            else
c     here i = 1
               mpole(i) = mpole(i)+(dreal(hexp(2))*dreal(htemp(2-i))
     1              +dimag(hexp(2))*dimag(htemp(2-i))
     2              -ima*(dreal(hexp(2))*dimag(htemp(2-i))
     3              -dimag(hexp(2))*dreal(htemp(2-i))))*rs2
               mpole(i) = mpole(i)+(dreal(hexp(2))*dreal(htemp(i+2))
     1              +dimag(hexp(2))*dimag(htemp(i+2))
     2              +ima*(dreal(hexp(2))*dimag(htemp(i+2))
     3              -dimag(hexp(2))*dreal(htemp(i+2))))

            endif
         enddo
      enddo

      return
      end
c     
c     
