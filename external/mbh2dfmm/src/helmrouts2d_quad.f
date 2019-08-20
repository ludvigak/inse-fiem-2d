cc Copyright (C) 2017: Travis Askham, Leslie Greengard, Zydrunas Gimbutas
cc email: askhamwhat@gmail.com      
cc 
cc This software is being released under a modified FreeBSD license
cc (see licenses folder in home directory). 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      Remarks on scaling conventions.
c
c      1)  Hankel and Bessel functions are consistently scaled as
c       	hvec(n)= H_n(z)*rscale^(n)
c       	jvec(n)= J_n(z)/rscale^(n)
c
c          rscale should be of the order of |z| if |z| < 1. Otherwise,
c          rscale should be set to 1.
c
c      2) The potential scaling is that required of the delta function
c      response:  pot = H_0(k*r)*(eye/4)
c
c-----------------------------------------------------------------------
c
c
c
C***********************************************************************
      subroutine h2dformmp_qp(ier,zk,rscale,source,quadstr,quadvec,
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
c     sum quadstr(j)*( quadvec(1,j)* (H_0)_xx (zk ||p-source_j||)
c                    + quadvec(2,j)* (H_0)_xy (zk ||p-source_j||)
c                    + quadvec(3,j)* (H_0)_yy (zk ||p-source_j||))
c     
c-----------------------------------------------------------------------
c     INPUT:
c     
c     zk              : Helmholtz parameter 
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
      complex *16 zk,quadstr(*),mpole(-nterms:nterms)
      integer ns, nterms, ier
      real *8 center(2),source(2,*),quadvec(3,*), rscale
c     local variables
      complex *16 z,ima, zmul,zinv,ztemp1,ztemp2,hexp(-2:2),zk2
      real *8 zdiff(2), theta, r, rs2, rs4, pi, done
      integer iii, i, j, ntop, lwfjs, ifder
      complex *16, allocatable :: jval(:), jder(:), jtemp(:)
      integer, allocatable :: iscale(:)
c     
      data ima/(0.0d0,1.0d0)/
c     
      ier = 0

      done=1
      pi=4*atan(done)
c     
      lwfjs = nterms+5 + 4*nterms + 100
      allocate(jval(0:lwfjs+10), stat=ier)
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
      allocate(jtemp(-nterms-5:nterms+5), stat=ier)
      if (ier.eq.1) then
         return
      endif
c     
      zk2 = zk*zk

      do i = -nterms,nterms
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
c     
         jtemp(0) = jval(0)
         zmul=exp(-ima*theta)
         zinv=conjg(zmul)
         ztemp1= zmul
         ztemp2=-zinv
         do j = 1,nterms+2
            jtemp( j) = ztemp1*jval(j)
            jtemp(-j) = ztemp2*jval(j)
            ztemp1= ztemp1*zmul
            ztemp2=-ztemp2*zinv
         enddo

c     set-up quadrupole as if centered at the source

         rs2 = rscale*rscale 

         hexp(-2) = quadstr(iii)*(zk2*quadvec(1,iii)
     1                           +zk2*quadvec(2,iii)*ima
     2                           -zk2*quadvec(3,iii))/(rs2*4.0d0)

         hexp(-1) = 0.0d0
         hexp(0) = quadstr(iii)*(-zk2*quadvec(1,iii)
     1                           -zk2*quadvec(3,iii))/2.0d0
         hexp(1) = 0.0d0
         hexp(2) = quadstr(iii)*(zk2*quadvec(1,iii)
     1                          -zk2*quadvec(2,iii)*ima
     2                          -zk2*quadvec(3,iii))/(rs2*4.0d0)


c     shift the quadrupole

         mpole(0) = mpole(0) + hexp(0)*jtemp(0)

         rs4 = rs2*rs2

         mpole(0) = mpole(0)+(hexp(2)*jtemp(-2))*rs4
         mpole(0) = mpole(0)+(hexp(-2)*jtemp(2))*rs4

         do i = 1,nterms

c     due to zero-th term

            mpole(i) = mpole(i) + hexp(0)*jtemp(i)
            mpole(-i) = mpole(-i) + hexp(0)*jtemp(-i)

c     due to 2nd term

            if (i .ge. 2) then
               mpole(i) = mpole(i)+(hexp(2)*jtemp(i-2))
               mpole(i) = mpole(i)+(hexp(-2)*jtemp(i+2))*rs4
               mpole(-i) = mpole(-i)+(hexp(2)*jtemp(-i-2))*rs4
               mpole(-i) = mpole(-i)+(hexp(-2)*jtemp(-i+2))
            else
c     i = 1 here
               mpole(i) = mpole(i)+(hexp(2)*jtemp(i-2))*rs2
               mpole(i) = mpole(i)+(hexp(-2)*jtemp(i+2))*rs4
               mpole(-i) = mpole(-i)+(hexp(2)*jtemp(-i-2))*rs4
               mpole(-i) = mpole(-i)+(hexp(-2)*jtemp(-i+2))*rs2
            endif
            
         enddo
         
      enddo

      return
      end
c
c
c
c
C***********************************************************************
      subroutine h2dformmp_qp_add(ier,zk,rscale,source,quadstr,
     1     quadvec,ns,center,nterms,mpole)
      implicit none
C***********************************************************************
c     
c     This subroutine constructs a multipole (h) expansion about 
c     CENTER due to NS quadrupole sources located at SOURCES(2,*).
c     
c     When evaluated at the point p (assumed to be well separated), 
c     the multipole expansion approximates the function
c
c     sum quadstr(j)*( quadvec(1,j)* (H_0)_xx (zk ||p-source_j||)
c                    + quadvec(2,j)* (H_0)_xy (zk ||p-source_j||)
c                    + quadvec(3,j)* (H_0)_yy (zk ||p-source_j||))
c     
c-----------------------------------------------------------------------
c     INPUT:
c     
c     zk              : Helmholtz parameter 
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
      complex *16 zk,quadstr(*),mpole(-nterms:nterms)
      integer ns, nterms, ier
      real *8 center(2),source(2,*),quadvec(3,*), rscale
c     local variables
      complex *16 z,ima, zmul,zinv,ztemp1,ztemp2,hexp(-2:2),zk2
      real *8 zdiff(2), theta, r, rs2, rs4, pi, done
      integer iii, i, j, ntop, lwfjs, ifder
      complex *16, allocatable :: jval(:), jder(:), jtemp(:)
      integer, allocatable :: iscale(:)
c     
      data ima/(0.0d0,1.0d0)/
c     
      ier = 0

      done=1
      pi=4*atan(done)
c     
      lwfjs = nterms+5 + 4*nterms + 100
      allocate(jval(0:lwfjs+10), stat=ier)
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
      allocate(jtemp(-nterms-5:nterms+5), stat=ier)
      if (ier.eq.1) then
         return
      endif
c     
      zk2 = zk*zk

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
c     
         jtemp(0) = jval(0)
         zmul=exp(-ima*theta)
         zinv=conjg(zmul)
         ztemp1= zmul
         ztemp2=-zinv
         do j = 1,nterms+2
            jtemp( j) = ztemp1*jval(j)
            jtemp(-j) = ztemp2*jval(j)
            ztemp1= ztemp1*zmul
            ztemp2=-ztemp2*zinv
         enddo

c     set-up quadrupole as if centered at the source

         rs2 = rscale*rscale 

         hexp(-2) = quadstr(iii)*(zk2*quadvec(1,iii)
     1                           +zk2*quadvec(2,iii)*ima
     2                           -zk2*quadvec(3,iii))/(rs2*4.0d0)
         hexp(-1) = 0.0d0
         hexp(0) = quadstr(iii)*(-zk2*quadvec(1,iii)
     1                           -zk2*quadvec(3,iii))/2.0d0
         hexp(1) = 0.0d0
         hexp(2) = quadstr(iii)*(zk2*quadvec(1,iii)
     1                          -zk2*quadvec(2,iii)*ima
     2                          -zk2*quadvec(3,iii))/(rs2*4.0d0)


c     shift the quadrupole

         mpole(0) = mpole(0) + hexp(0)*jtemp(0)

         rs4 = rs2*rs2

         mpole(0) = mpole(0)+(hexp(2)*jtemp(-2))*rs4
         mpole(0) = mpole(0)+(hexp(-2)*jtemp(2))*rs4

         do i = 1,nterms

c     due to zero-th term

            mpole(i) = mpole(i) + hexp(0)*jtemp(i)
            mpole(-i) = mpole(-i) + hexp(0)*jtemp(-i)

c     due to 2nd term

            if (i .ge. 2) then
               mpole(i) = mpole(i)+(hexp(2)*jtemp(i-2))
               mpole(i) = mpole(i)+(hexp(-2)*jtemp(i+2))*rs4
               mpole(-i) = mpole(-i)+(hexp(2)*jtemp(-i-2))*rs4
               mpole(-i) = mpole(-i)+(hexp(-2)*jtemp(-i+2))
            else
c     i = 1 here
               mpole(i) = mpole(i)+(hexp(2)*jtemp(i-2))*rs2
               mpole(i) = mpole(i)+(hexp(-2)*jtemp(i+2))*rs4
               mpole(-i) = mpole(-i)+(hexp(2)*jtemp(-i-2))*rs4
               mpole(-i) = mpole(-i)+(hexp(-2)*jtemp(-i+2))*rs2
            endif
            
         enddo
         
      enddo

      return
      end
c
c
c
C***********************************************************************
      subroutine h2dformta_qp(ier,zk,rscale,source,quadstr,quadvec,
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
c     zk              : Helmholtz parameter 
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
      complex *16 zk,mpole(-nterms:nterms),quadstr(*)
      real *8 center(2),source(2,*),quadvec(3,*)
      real *8 rscale
      integer ns, nterms, ier
c     local variables
      complex *16 hexp(-2:2)
      complex *16, allocatable :: hval(:), htemp(:)
      complex *16, allocatable :: hder(:)
      complex *16 zmul,zinv,ztemp1,ztemp2,zk2
      real *8 theta, r, rs2, rs4, zdiff(2)
      integer i, iii, ifder
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
      allocate(htemp(-nterms-5:nterms+5), stat=ier)
      if (ier.eq.1) then
         return
      endif

c     
      do i=-nterms,nterms
         mpole(i)=0
      enddo

      zk2 = zk*zk
      
      do iii=1,ns

c     get h vals for this source

         zdiff(1)=source(1,iii)-center(1)
         zdiff(2)=source(2,iii)-center(2)
         call h2cart2polar(zdiff,r,theta)
         z=zk*r
         ifder=0
         call h2dall(nterms+2,z,rscale,hval,ifder,hder)
         htemp(0) = hval(0)
         zmul=exp(-ima*theta)
         zinv=conjg(zmul)
         ztemp1= zmul
         ztemp2=-zinv
         do i = 1,nterms+2
            htemp( i) = ztemp1*hval(i)
            htemp(-i) = ztemp2*hval(i)
            ztemp1= ztemp1*zmul
            ztemp2=-ztemp2*zinv
         enddo

c     set-up quadrupole as if centered at the source

         rs2 = rscale*rscale 

         hexp(-2) = quadstr(iii)*(zk2*quadvec(1,iii)
     1                           +zk2*quadvec(2,iii)*ima
     2                           -zk2*quadvec(3,iii))/(rs2*4.0d0)
         hexp(-1) = 0.0d0
         hexp(0) = quadstr(iii)*(-zk2*quadvec(1,iii)
     1                           -zk2*quadvec(3,iii))/2.0d0
         hexp(1) = 0.0d0
         hexp(2) = quadstr(iii)*(zk2*quadvec(1,iii)
     1                          -zk2*quadvec(2,iii)*ima
     2                          -zk2*quadvec(3,iii))/(rs2*4.0d0)

c     shift the quadrupole

         mpole(0) = mpole(0) + hexp(0)*htemp(0)

         mpole(0) = mpole(0)+(hexp(2)*htemp(-2))
         mpole(0) = mpole(0)+(hexp(-2)*htemp(2))
c     
         rs2 = rscale*rscale
         rs4 = rs2*rs2
         do i = 1,nterms
            mpole(i) = mpole(i) + hexp(0)*htemp(i)
            mpole(-i) = mpole(-i) + hexp(0)*htemp(-i)
            if (i .ge. 2) then
               mpole(i) = mpole(i)+(hexp(2)*htemp(i-2))*rs4
               mpole(i) = mpole(i)+(hexp(-2)*htemp(i+2))
               mpole(-i) = mpole(-i)+(hexp(2)*htemp(-i-2))
               mpole(-i) = mpole(-i)+(hexp(-2)*htemp(-i+2))*rs4
            else
c     here i = 1
               mpole(i) = mpole(i)+(hexp(2)*htemp(i-2))*rs2
               mpole(i) = mpole(i)+(hexp(-2)*htemp(i+2))
               mpole(-i) = mpole(-i)+(hexp(2)*htemp(-i-2))
               mpole(-i) = mpole(-i)+(hexp(-2)*htemp(-i+2))*rs2
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
      subroutine h2dformta_qp_add(ier,zk,rscale,source,quadstr,
     1     quadvec,ns,center,nterms,mpole)
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
c     zk              : Helmholtz parameter 
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
      complex *16 zk,mpole(-nterms:nterms),quadstr(*)
      real *8 center(2),source(2,*),quadvec(3,*)
      real *8 rscale
      integer ns, nterms, ier
c     local variables
      complex *16 hexp(-2:2)
      complex *16, allocatable :: hval(:), htemp(:)
      complex *16, allocatable :: hder(:)
      complex *16 zmul,zinv,ztemp1,ztemp2,zk2
      real *8 theta, r, rs2, rs4, zdiff(2)
      integer i, iii, ifder
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
      allocate(htemp(-nterms-5:nterms+5), stat=ier)
      if (ier.eq.1) then
         return
      endif

      zk2 = zk*zk
      
      do iii=1,ns

c     get h vals for this source

         zdiff(1)=source(1,iii)-center(1)
         zdiff(2)=source(2,iii)-center(2)
         call h2cart2polar(zdiff,r,theta)
         z=zk*r
         ifder=0
         call h2dall(nterms+2,z,rscale,hval,ifder,hder)
         htemp(0) = hval(0)
         zmul=exp(-ima*theta)
         zinv=conjg(zmul)
         ztemp1= zmul
         ztemp2=-zinv
         do i = 1,nterms+2
            htemp( i) = ztemp1*hval(i)
            htemp(-i) = ztemp2*hval(i)
            ztemp1= ztemp1*zmul
            ztemp2=-ztemp2*zinv
         enddo

c     set-up quadrupole as if centered at the source

         rs2 = rscale*rscale 

         hexp(-2) = quadstr(iii)*(zk2*quadvec(1,iii)
     1                           +zk2*quadvec(2,iii)*ima
     2                           -zk2*quadvec(3,iii))/(rs2*4.0d0)
         hexp(-1) = 0.0d0
         hexp(0) = quadstr(iii)*(-zk2*quadvec(1,iii)
     1                           -zk2*quadvec(3,iii))/2.0d0
         hexp(1) = 0.0d0
         hexp(2) = quadstr(iii)*(zk2*quadvec(1,iii)
     1                          -zk2*quadvec(2,iii)*ima
     2                          -zk2*quadvec(3,iii))/(rs2*4.0d0)

c     shift the quadrupole

         mpole(0) = mpole(0) + hexp(0)*htemp(0)

         mpole(0) = mpole(0)+(hexp(2)*htemp(-2))
         mpole(0) = mpole(0)+(hexp(-2)*htemp(2))
c     
         rs2 = rscale*rscale
         rs4 = rs2*rs2
         do i = 1,nterms
            mpole(i) = mpole(i) + hexp(0)*htemp(i)
            mpole(-i) = mpole(-i) + hexp(0)*htemp(-i)
            if (i .ge. 2) then
               mpole(i) = mpole(i)+(hexp(2)*htemp(i-2))*rs4
               mpole(i) = mpole(i)+(hexp(-2)*htemp(i+2))
               mpole(-i) = mpole(-i)+(hexp(2)*htemp(-i-2))
               mpole(-i) = mpole(-i)+(hexp(-2)*htemp(-i+2))*rs4
            else
c     here i = 1
               mpole(i) = mpole(i)+(hexp(2)*htemp(i-2))*rs2
               mpole(i) = mpole(i)+(hexp(-2)*htemp(i+2))
               mpole(-i) = mpole(-i)+(hexp(2)*htemp(-i-2))
               mpole(-i) = mpole(-i)+(hexp(-2)*htemp(-i+2))*rs2
            endif
         enddo
      enddo

      return
      end
c     
c     
c     
