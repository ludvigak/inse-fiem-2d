cc Copyright (C) 2009-2012: Leslie Greengard, Travis Askham, and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This software is being released under a modified FreeBSD license
cc (see licenses folder in home directory). 
c
c
c      This file contains some basic subroutines for 
c      forming and evaluating multipole (partial wave) expansions
c      in two dimensions.
c
c-----------------------------------------------------------------------
c
C***********************************************************************
      subroutine l2dformmp_add(ier,rscale,source,charge,ns,center,
     1     nterms,mpole)
      implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a multipole expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c     mpole_0  =  sum charge_j 
c                  j  
c
c     mpole_n  = -sum charge_j 1/n (z)^n /rscale^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
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
c     mpole     : coeffs for the multipole-expansion
c
        complex *16 mpole(0:nterms),charge(ns)
        real *8 center(2),source(2,ns),zdiff(2)

        complex *16 zmul,zinv,ztemp1,ztemp2
c
        complex *16 ima,z,z0
        data ima/(0.0d0,1.0d0)/
c
c

        do j=1,ns
c
        zdiff(1)=source(1,j)-center(1)
        zdiff(2)=source(2,j)-center(2)
ccc        call h2cart2polar(zdiff,r,theta)
        z0=dcmplx(zdiff(1),zdiff(2))
c
        mpole(0)=mpole(0)+charge(j)
c
        ztemp1=z0/rscale
        ztemp2=ztemp1
        do n=1,nterms
        mpole(n)=mpole(n)-charge(j)*ztemp1/n
        ztemp1=ztemp1*ztemp2
        enddo

        enddo

        return
        end
c
c      
c     
      subroutine l2dformmp_dp_add(ier,rscale,source,dipstr,ns,center,
     1     nterms,mpole)
      implicit real *8 (a-h,o-z)
C***********************************************************************
c
c     This subroutine constructs a multipole expansion about CENTER due
c     to NS sources located at SOURCES(2,*).
c
c     mpole_0  =  0
c                 
c
c     mpole_n  =  sum dipstr_j (z_0)^(n-1)/z^n /rscale^n
c                  j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     dipstr(ns)      : source strengths
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c
c     OUTPUT:
c
c     ier       : error return code
c                 ier=0 returned successfully;
c
c     mpole     : coeffs for the multipole-expansion
c     
      complex *16 mpole(0:nterms),dipstr(ns)
      real *8 center(2),source(2,ns),zdiff(2)

      complex *16 zmul,zinv,ztemp1,ztemp2
c     
      complex *16 ima,z,z0
      data ima/(0.0d0,1.0d0)/
c     
c     
      do j=1,ns
c     
         zdiff(1)=source(1,j)-center(1)
         zdiff(2)=source(2,j)-center(2)
ccc   call h2cart2polar(zdiff,r,theta)
         z0=dcmplx(zdiff(1),zdiff(2))
c     
c     
         zmul=z0/rscale
         ztemp1=1/rscale
         do n=1,nterms
            mpole(n)=mpole(n)+dipstr(j)*ztemp1
            ztemp1=ztemp1*zmul
         enddo

      enddo

      return
      end
c     
c     
C***********************************************************************
      subroutine l2dformmp_qp(ier,rscale,source,quadstr,quadvec,
     $     ns,center,nterms,mpole)
      implicit none
C***********************************************************************
c     
c     This subroutine constructs a multipole expansion about 
c     CENTER due to NS quadrupole sources located at SOURCES(2,*).
c     
c     When evaluated at the point p (assumed to be well separated), 
c     the multipole expansion approximates the function
c
c     sum quadstr(j)*( quadvec(1,j)* (S_0)_xx (||p-source_j||)
c                    + quadvec(2,j)* (S_0)_xy (||p-source_j||)
c                    + quadvec(3,j)* (S_0)_yy (||p-source_j||))
c     
c-----------------------------------------------------------------------
c     INPUT:
c     
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     quadstr(ns)     : quadrupole strengths
c     quadvec(3,ns)   : quadrupole ``orientation''
c                     quadvec(1,i) - (S_0)_xx term
c                     quadvec(2,i) - (S_0)_xy term
c                     quadvec(3,i) - (S_0)_yy term
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c     
c     OUTPUT:
c     
c     ier       : error return code
c     ier=0 returned successfully;
c     
c     mpole     : coeffs for the multipole expansion
c     
      
c     global variables
      complex *16 quadstr(*),mpole(0:nterms)
      integer ns, nterms, ier
      real *8 center(2),source(2,*),quadvec(3,*), rscale
c     local variables
      complex *16 z0,ima, zmul,zinv,ztemp1,ztemp2,hexp(0:2)
      real *8 zdiff(2), rinv2
      integer iii, i, j, nmax, m, l
      complex *16, allocatable :: z0pow(:)
      real *8, allocatable :: carray(:,:)
      integer, allocatable :: iscale(:)
c     
      data ima/(0.0d0,1.0d0)/
c     
      ier = 0

      nmax = max(nterms,2)
      allocate( carray(0:nmax,0:2) )
      allocate( z0pow(0:nmax) )
      
      do l = 0,nmax
         carray(l,0) = 1.0d0
      enddo
      do m=1,2
         carray(m,m) = 1.0d0
         do l=m+1,nmax
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
         enddo
      enddo

      do i = 0,nterms
         mpole(i) = 0
      enddo
         
      rinv2 = 1.0d0/(rscale**2)

      do iii = 1,ns

c     set-up quadrupole as if centered at source
         
         zdiff(1) = center(1)-source(1,iii)
         zdiff(2) = center(2)-source(2,iii)
         z0 = dcmplx(-zdiff(1),-zdiff(2))

         hexp(0) = 0.0d0
         hexp(1) = 0.0d0
         hexp(2) = quadstr(iii)*
     1        (quadvec(1,iii)+ima*quadvec(2,iii)-quadvec(3,iii))*rinv2
         
c     shift expansion to center 

         ztemp1=z0/rscale
         ztemp2=ztemp1
         z0pow(0)=1
         do i=1,nmax
            z0pow(i)=ztemp1
            ztemp1=ztemp1*ztemp2
         enddo
c     
         do i = 1,nterms
            do j = 2,min(i,2)
               mpole(i) = mpole(i) +
     1              hexp(j)*z0pow(i)/z0pow(j)*carray(i-1,j-1)
            enddo
         enddo
         
      enddo

      return
      end
c
c
c
c
C***********************************************************************
      subroutine l2dformmp_qp_add(ier,rscale,source,quadstr,
     1     quadvec,ns,center,nterms,mpole)
      implicit none
C***********************************************************************
c     
c     This subroutine constructs a multipole expansion about 
c     CENTER due to NS quadrupole sources located at SOURCES(2,*).
c     
c     When evaluated at the point p (assumed to be well separated), 
c     the multipole expansion approximates the function
c
c     sum quadstr(j)*( quadvec(1,j)* (S_0)_xx (||p-source_j||)
c                    + quadvec(2,j)* (S_0)_xy (||p-source_j||)
c                    + quadvec(3,j)* (S_0)_yy (||p-source_j||))
c     
c-----------------------------------------------------------------------
c     INPUT:
c     
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     quadstr(ns)     : quadrupole strengths
c     quadvec(3,ns)   : quadrupole ``orientation''
c                     quadvec(1,i) - (S_0)_xx term
c                     quadvec(2,i) - (S_0)_xy term
c                     quadvec(3,i) - (S_0)_yy term
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c     
c     OUTPUT:
c     
c     ier       : error return code
c     ier=0 returned successfully;
c     
c     mpole     : coeffs for the multipole expansion
c     
      
c     global variables
      complex *16 quadstr(*),mpole(0:nterms)
      integer ns, nterms, ier
      real *8 center(2),source(2,*),quadvec(3,*), rscale
c     local variables
      complex *16 z0,ima, zmul,zinv,ztemp1,ztemp2,hexp(0:2)
      real *8 zdiff(2), rinv2
      integer iii, i, j, nmax, m, l
      complex *16, allocatable :: z0pow(:)
      real *8, allocatable :: carray(:,:)
      integer, allocatable :: iscale(:)
c     
      data ima/(0.0d0,1.0d0)/
c     
      ier = 0

      nmax = max(nterms,2)
      allocate( carray(0:nmax,0:2) )
      allocate( z0pow(0:nmax) )
      
      do l = 0,nmax
         carray(l,0) = 1.0d0
      enddo
      do m=1,2
         carray(m,m) = 1.0d0
         do l=m+1,nmax
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
         enddo
      enddo

      rinv2 = 1.0d0/(rscale**2)

      do iii = 1,ns

c     set-up quadrupole as if centered at source
         
         zdiff(1) = center(1)-source(1,iii)
         zdiff(2) = center(2)-source(2,iii)
         z0 = dcmplx(-zdiff(1),-zdiff(2))

         hexp(0) = 0.0d0
         hexp(1) = 0.0d0
         hexp(2) = quadstr(iii)*
     1        (quadvec(1,iii)+ima*quadvec(2,iii)-quadvec(3,iii))*rinv2
         
c     shift expansion to center 

         ztemp1=z0/rscale
         ztemp2=ztemp1
         z0pow(0)=1
         do i=1,nmax
            z0pow(i)=ztemp1
            ztemp1=ztemp1*ztemp2
         enddo
c     
         do i = 1,nterms
            do j = 2,min(i,2)
               mpole(i) = mpole(i) +
     1              hexp(j)*z0pow(i)/z0pow(j)*carray(i-1,j-1)
            enddo
         enddo
         
      enddo

      return
      end
c
c
C***********************************************************************
      subroutine l2dformta_qp(ier,rscale,source,quadstr,quadvec,
     $     ns,center,nterms,loc)
      implicit none
C***********************************************************************
c     
c     This subroutine constructs a multipole expansion about 
c     CENTER due to NS quadrupole sources located at SOURCES(2,*).
c     
c     When evaluated at the point p (assumed to be well separated), 
c     the multipole expansion approximates the function
c
c     sum quadstr(j)*( quadvec(1,j)* (S_0)_xx (||p-source_j||)
c                    + quadvec(2,j)* (S_0)_xy (||p-source_j||)
c                    + quadvec(3,j)* (S_0)_yy (||p-source_j||))
c     
c-----------------------------------------------------------------------
c     INPUT:
c     
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     quadstr(ns)     : quadrupole strengths
c     quadvec(3,ns)   : quadrupole ``orientation''
c                     quadvec(1,i) - (S_0)_xx term
c                     quadvec(2,i) - (S_0)_xy term
c                     quadvec(3,i) - (S_0)_yy term
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c     
c     OUTPUT:
c     
c     ier       : error return code
c     ier=0 returned successfully;
c     
c     mpole     : coeffs for the multipole expansion
c     
      
c     global variables
      complex *16 quadstr(*),loc(0:nterms)
      integer ns, nterms, ier
      real *8 center(2),source(2,*),quadvec(3,*), rscale
c     local variables
      complex *16 z0,ima, zmul,zinv,ztemp1,ztemp2,hexp(0:2)
      real *8 zdiff(2), rinv2
      integer iii, i, j, nmax, m, l
      complex *16, allocatable :: z0pow(:), z0powm(:)
      real *8, allocatable :: carray(:,:)
      integer, allocatable :: iscale(:)
c     
      data ima/(0.0d0,1.0d0)/
c     
      ier = 0

      nmax = nterms+2
      allocate( carray(0:nmax,0:2) )
      allocate( z0pow(0:nmax) )
      allocate( z0powm(0:nmax) )
      
      do l = 0,nmax
         carray(l,0) = 1.0d0
      enddo
      do m=1,2
         carray(m,m) = 1.0d0
         do l=m+1,nmax
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
         enddo
      enddo

      do i = 0,nterms
         loc(i) = 0
      enddo
         
      rinv2 = 1.0d0/(rscale**2)

      do iii = 1,ns

c     set-up quadrupole as if centered at source
         
         zdiff(1) = center(1)-source(1,iii)
         zdiff(2) = center(2)-source(2,iii)
         z0 = dcmplx(-zdiff(1),-zdiff(2))

         hexp(0) = 0.0d0
         hexp(1) = 0.0d0
         hexp(2) = quadstr(iii)*
     1        (quadvec(1,iii)+ima*quadvec(2,iii)-quadvec(3,iii))*rinv2
         
c     shift expansion to center 

         ztemp1=rscale/z0
         ztemp2=ztemp1
         z0pow(0)=1
         do i=1,nmax
            z0pow(i)=ztemp1
            ztemp1=ztemp1*ztemp2
         enddo

         ztemp1=-rscale/z0
         ztemp2=ztemp1
         z0powm(0)=1
         do i=1,2
            z0powm(i)=ztemp1
            ztemp1=ztemp1*ztemp2
         enddo
c     
         loc(0) = loc(0) + hexp(2)*z0powm(2)

         do i = 1,nterms
            loc(i) = loc(i) + hexp(2)*z0powm(2)*z0pow(i)*
     1           carray(i+2-1,2-1)
         enddo
         
      enddo

      return
      end
c
c     
c     
C***********************************************************************
      subroutine l2dformta_qp_add(ier,rscale,source,quadstr,
     1     quadvec,ns,center,nterms,loc)
      implicit none
C***********************************************************************
c     
c     This subroutine constructs a multipole expansion about 
c     CENTER due to NS quadrupole sources located at SOURCES(2,*).
c     
c     When evaluated at the point p (assumed to be well separated), 
c     the multipole expansion approximates the function
c
c     sum quadstr(j)*( quadvec(1,j)* (S_0)_xx (||p-source_j||)
c                    + quadvec(2,j)* (S_0)_xy (||p-source_j||)
c                    + quadvec(3,j)* (S_0)_yy (||p-source_j||))
c     
c-----------------------------------------------------------------------
c     INPUT:
c     
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     quadstr(ns)     : quadrupole strengths
c     quadvec(3,ns)   : quadrupole ``orientation''
c                     quadvec(1,i) - (S_0)_xx term
c                     quadvec(2,i) - (S_0)_xy term
c                     quadvec(3,i) - (S_0)_yy term
c     ns              : number of sources
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c     
c     OUTPUT:
c     
c     ier       : error return code
c     ier=0 returned successfully;
c     
c     mpole     : coeffs for the multipole expansion
c     
      
c     global variables
      complex *16 quadstr(*),loc(0:nterms)
      integer ns, nterms, ier
      real *8 center(2),source(2,*),quadvec(3,*), rscale
c     local variables
      complex *16 z0,ima, zmul,zinv,ztemp1,ztemp2,hexp(0:2)
      real *8 zdiff(2), rinv2
      integer iii, i, j, nmax, m, l
      complex *16, allocatable :: z0pow(:), z0powm(:)
      real *8, allocatable :: carray(:,:)
      integer, allocatable :: iscale(:)
c     
      data ima/(0.0d0,1.0d0)/
c     
      ier = 0

      nmax = nterms+2
      allocate( carray(0:nmax,0:2) )
      allocate( z0pow(0:nmax) )
      allocate( z0powm(0:nmax) )
      
      do l = 0,nmax
         carray(l,0) = 1.0d0
      enddo
      do m=1,2
         carray(m,m) = 1.0d0
         do l=m+1,nmax
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
         enddo
      enddo

      rinv2 = 1.0d0/(rscale**2)

      do iii = 1,ns

c     set-up quadrupole as if centered at source
         
         zdiff(1) = center(1)-source(1,iii)
         zdiff(2) = center(2)-source(2,iii)
         z0 = dcmplx(-zdiff(1),-zdiff(2))

         hexp(0) = 0.0d0
         hexp(1) = 0.0d0
         hexp(2) = quadstr(iii)*
     1        (quadvec(1,iii)+ima*quadvec(2,iii)-quadvec(3,iii))*rinv2
         
c     shift expansion to center 

         ztemp1=rscale/z0
         ztemp2=ztemp1
         z0pow(0)=1
         do i=1,nmax
            z0pow(i)=ztemp1
            ztemp1=ztemp1*ztemp2
         enddo

         ztemp1=-rscale/z0
         ztemp2=ztemp1
         z0powm(0)=1
         do i=1,2
            z0powm(i)=ztemp1
            ztemp1=ztemp1*ztemp2
         enddo
c     
         loc(0) = loc(0) + hexp(2)*z0powm(2)

         do i = 1,nterms
            loc(i) = loc(i) + hexp(2)*z0powm(2)*z0pow(i)*
     1           carray(i+2-1,2-1)
         enddo
         
      enddo

      return
      end
c
