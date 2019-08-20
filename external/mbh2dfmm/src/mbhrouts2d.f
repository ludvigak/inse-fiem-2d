cc Copyright (C) 2017: Travis Askham, Leslie Greengard, Zydrunas Gimbutas
cc email: askhamwhat@gmail.com      
cc 
cc This software is being released under a modified FreeBSD license
cc (see licenses folder in home directory). 


      subroutine mbhpotgrad2dall_cdq(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     target,ifpot,pot,ifgrad,grad,ifhess,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*), charge(*), dipstr(*)
      real *8 dipvec(2,*), quadstr(*), quadvec(3,*)
      real *8 target(2)
      real *8 pot, grad(2), hess(3)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      ifpotloc = 1
      ifgradloc = 1
      ifhessloc = 1
      ifder3 = 1
      ifder4 = 1
      ifder5 = 0

      if (ifpot .eq. 1) pot = 0.0d0
      if (ifgrad .eq. 1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0
      endif
      if (ifhess .eq. 1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0
      endif

      do i = 1,ns
         call modbhgreen_all(beta,target,source(1,i),ifpot,potloc,
     1        ifgrad,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,der4,
     2        ifder5,der5)

         if (ifcharge .eq. 1) then
            if (ifpot .eq. 1) then
               pot = pot + potloc*charge(i)
            endif
            if (ifgrad .eq. 1) then
               grad(1) = grad(1) + gradloc(1)*charge(i)
               grad(2) = grad(2) + gradloc(2)*charge(i)
            endif
            if (ifhess .eq. 1) then
               hess(1) = hess(1) + hessloc(1)*charge(i)
               hess(2) = hess(2) + hessloc(2)*charge(i)
               hess(3) = hess(3) + hessloc(3)*charge(i)
            endif
         endif

         if (ifdipole .eq. 1) then
            if (ifpot .eq. 1) then
               pot = pot - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1              + gradloc(2)*dipvec(2,i))
            endif
            if (ifgrad .eq. 1) then
               grad(1) = grad(1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1              + hessloc(2)*dipvec(2,i))
               grad(2) = grad(2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1              + hessloc(3)*dipvec(2,i))
            endif
            if (ifhess .eq. 1) then
               hess(1) = hess(1) - dipstr(i)*(der3(1)*dipvec(1,i)
     1              + der3(2)*dipvec(2,i))
               hess(2) = hess(2) - dipstr(i)*(der3(2)*dipvec(1,i)
     1              + der3(3)*dipvec(2,i))
               hess(3) = hess(3) - dipstr(i)*(der3(3)*dipvec(1,i)
     1              + der3(4)*dipvec(2,i))
            endif
         endif

         if (ifquad .eq. 1) then
            if (ifpot .eq. 1) then
               pot = pot + quadstr(i)*(hessloc(1)*quadvec(1,i)
     1              + hessloc(2)*quadvec(2,i) + hessloc(3)*quadvec(3,i))
            endif
            if (ifgrad .eq. 1) then
               grad(1) = grad(1) + quadstr(i)*(der3(1)*quadvec(1,i)
     1              + der3(2)*quadvec(2,i) + der3(3)*quadvec(3,i))
               grad(2) = grad(2) + quadstr(i)*(der3(2)*quadvec(1,i)
     1              + der3(3)*quadvec(2,i) + der3(4)*quadvec(3,i))
            endif
            if (ifhess .eq. 1) then
               hess(1) = hess(1) + quadstr(i)*(der4(1)*quadvec(1,i)
     1              + der4(2)*quadvec(2,i) + der4(3)*quadvec(3,i))
               hess(2) = hess(2) + quadstr(i)*(der4(2)*quadvec(1,i)
     1              + der4(3)*quadvec(2,i) + der4(4)*quadvec(3,i))
               hess(3) = hess(3) + quadstr(i)*(der4(3)*quadvec(1,i)
     1              + der4(4)*quadvec(2,i) + der4(5)*quadvec(3,i))
            endif
         endif
      enddo

      return
      end

      subroutine mbhpotgrad2dall_cdq_add(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     target,ifpot,pot,ifgrad,grad,ifhess,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*), charge(*), dipstr(*)
      real *8 dipvec(2,*), quadstr(*), quadvec(3,*)
      real *8 target(2)
      real *8 pot, grad(2), hess(3)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
c     local variables
      real *8 pottemp, gradtemp(2), hesstemp(3)

      call mbhpotgrad2dall_cdq(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     target,ifpot,pottemp,ifgrad,gradtemp,ifhess,hesstemp)

      if (ifpot .eq. 1) pot = pot + pottemp
      if (ifgrad .eq. 1) then
         grad(1) = grad(1) + gradtemp(1)
         grad(2) = grad(2) + gradtemp(2)
      endif
      if (ifhess .eq. 1) then
         hess(1) = hess(1) + hesstemp(1)
         hess(2) = hess(2) + hesstemp(2)
         hess(3) = hess(3) + hesstemp(3)
      endif

      return
      end
  
      subroutine mbhpotgrad2dall_cdq3(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     target,ifpot,pot,ifgrad,grad,ifhess,hess,ifder3,der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*), charge(*), dipstr(*)
      real *8 dipvec(2,*), quadstr(*), quadvec(3,*)
      real *8 target(2)
      real *8 pot, grad(2), hess(3), der3(4)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
      integer ifder3
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3loc(4),der4(5),der5(6)
      integer ifpotloc, ifgradloc, ifhessloc, ifder3loc, ifder4, ifder5
      integer i

      ifpotloc = 1
      ifgradloc = 1
      ifhessloc = 1
      ifder3loc = 1
      ifder4 = 1
      ifder5 = 1

      if (ifpot .eq. 1) pot = 0.0d0
      if (ifgrad .eq. 1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0
      endif
      if (ifhess .eq. 1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0
      endif
      if (ifder3 .eq. 1) then
         der3(1) = 0.0d0
         der3(2) = 0.0d0
         der3(3) = 0.0d0
         der3(4) = 0.0d0
      endif

      do i = 1,ns
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3loc,der3loc,
     2        ifder4,der4,ifder5,der5)

         if (ifcharge .eq. 1) then
            if (ifpot .eq. 1) then
               pot = pot + potloc*charge(i)
            endif
            if (ifgrad .eq. 1) then
               grad(1) = grad(1) + gradloc(1)*charge(i)
               grad(2) = grad(2) + gradloc(2)*charge(i)
            endif
            if (ifhess .eq. 1) then
               hess(1) = hess(1) + hessloc(1)*charge(i)
               hess(2) = hess(2) + hessloc(2)*charge(i)
               hess(3) = hess(3) + hessloc(3)*charge(i)
            endif
            if (ifder3 .eq. 1) then
               der3(1) = der3(1) + der3loc(1)*charge(i)
               der3(2) = der3(2) + der3loc(2)*charge(i)
               der3(3) = der3(3) + der3loc(3)*charge(i)
               der3(4) = der3(4) + der3loc(4)*charge(i)
            endif
         endif

         if (ifdipole .eq. 1) then
            if (ifpot .eq. 1) then
               pot = pot - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1              + gradloc(2)*dipvec(2,i))
            endif
            if (ifgrad .eq. 1) then
               grad(1) = grad(1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1              + hessloc(2)*dipvec(2,i))
               grad(2) = grad(2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1              + hessloc(3)*dipvec(2,i))
            endif
            if (ifhess .eq. 1) then
               hess(1) = hess(1) - dipstr(i)*(der3loc(1)*dipvec(1,i)
     1              + der3loc(2)*dipvec(2,i))
               hess(2) = hess(2) - dipstr(i)*(der3loc(2)*dipvec(1,i)
     1              + der3loc(3)*dipvec(2,i))
               hess(3) = hess(3) - dipstr(i)*(der3loc(3)*dipvec(1,i)
     1              + der3loc(4)*dipvec(2,i))
            endif
            if (ifder3 .eq. 1) then
               der3(1) = der3(1) - dipstr(i)*(der4(1)*dipvec(1,i)
     1              + der4(2)*dipvec(2,i))
               der3(2) = der3(2) - dipstr(i)*(der4(2)*dipvec(1,i)
     1              + der4(3)*dipvec(2,i))
               der3(3) = der3(3) - dipstr(i)*(der4(3)*dipvec(1,i)
     1              + der4(4)*dipvec(2,i))
               der3(4) = der3(4) - dipstr(i)*(der4(4)*dipvec(1,i)
     1              + der4(5)*dipvec(2,i))
            endif
         endif

         if (ifquad .eq. 1) then
            if (ifpot .eq. 1) then
               pot = pot + quadstr(i)*(hessloc(1)*quadvec(1,i)
     1              + hessloc(2)*quadvec(2,i) + hessloc(3)*quadvec(3,i))
            endif
            if (ifgrad .eq. 1) then
               grad(1) = grad(1) + quadstr(i)*(der3loc(1)*quadvec(1,i)
     1              + der3loc(2)*quadvec(2,i) + der3loc(3)*quadvec(3,i))
               grad(2) = grad(2) + quadstr(i)*(der3loc(2)*quadvec(1,i)
     1              + der3loc(3)*quadvec(2,i) + der3loc(4)*quadvec(3,i))
            endif
            if (ifhess .eq. 1) then
               hess(1) = hess(1) + quadstr(i)*(der4(1)*quadvec(1,i)
     1              + der4(2)*quadvec(2,i) + der4(3)*quadvec(3,i))
               hess(2) = hess(2) + quadstr(i)*(der4(2)*quadvec(1,i)
     1              + der4(3)*quadvec(2,i) + der4(4)*quadvec(3,i))
               hess(3) = hess(3) + quadstr(i)*(der4(3)*quadvec(1,i)
     1              + der4(4)*quadvec(2,i) + der4(5)*quadvec(3,i))
            endif
            if (ifder3 .eq. 1) then
               der3(1) = der3(1) + quadstr(i)*(der5(1)*quadvec(1,i)
     1              + der5(2)*quadvec(2,i) + der5(3)*quadvec(3,i))
               der3(2) = der3(2) + quadstr(i)*(der5(2)*quadvec(1,i)
     1              + der5(3)*quadvec(2,i) + der5(4)*quadvec(3,i))
               der3(3) = der3(3) + quadstr(i)*(der5(3)*quadvec(1,i)
     1              + der5(4)*quadvec(2,i) + der5(5)*quadvec(3,i))
               der3(4) = der3(4) + quadstr(i)*(der5(4)*quadvec(1,i)
     1              + der5(5)*quadvec(2,i) + der5(6)*quadvec(3,i))
            endif
         endif
      enddo
      
      return
      end

      subroutine mbhpotgrad2dall_cdq3_add(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     target,ifpot,pot,ifgrad,grad,ifhess,hess,ifder3,der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*), charge(*), dipstr(*)
      real *8 dipvec(2,*), quadstr(*), quadvec(3,*)
      real *8 target(2)
      real *8 pot, grad(2), hess(3), der3(4)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
      integer ifder3
c     local variables
      real *8 pottemp, gradtemp(2), hesstemp(3), der3temp(4)

      call mbhpotgrad2dall_cdq3(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     target,ifpot,pottemp,ifgrad,gradtemp,ifhess,hesstemp,
     3     ifder3,der3temp)

      if (ifpot .eq. 1) pot = pot + pottemp
      if (ifgrad .eq. 1) then
         grad(1) = grad(1) + gradtemp(1)
         grad(2) = grad(2) + gradtemp(2)
      endif
      if (ifhess .eq. 1) then
         hess(1) = hess(1) + hesstemp(1)
         hess(2) = hess(2) + hesstemp(2)
         hess(3) = hess(3) + hesstemp(3)
      endif
      if (ifder3 .eq. 1) then
         der3(1) = der3(1) + der3temp(1)
         der3(2) = der3(2) + der3temp(2)
         der3(3) = der3(3) + der3temp(3)
         der3(4) = der3(4) + der3temp(4)
      endif

      return
      end
  
      subroutine mbhpotgrad2dall_cdqo(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     ifoct,octstr,octvec,target,ifpot,pot,ifgrad,grad,
     3     ifhess,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*), charge(*), dipstr(*)
      real *8 dipvec(2,*), quadstr(*), quadvec(3,*)
      real *8 octstr(*), octvec(4,*)
      real *8 target(2)
      real *8 pot, grad(2), hess(3)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
      integer ifoct
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      ifpotloc = 1
      ifgradloc = 1
      ifhessloc = 1
      ifder3 = 1
      ifder4 = 1
      ifder5 = 1

      if (ifpot .eq. 1) pot = 0.0d0
      if (ifgrad .eq. 1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0
      endif
      if (ifhess .eq. 1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0
      endif

      do i = 1,ns
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         if (ifcharge .eq. 1) then
            if (ifpot .eq. 1) then
               pot = pot + potloc*charge(i)
            endif
            if (ifgrad .eq. 1) then
               grad(1) = grad(1) + gradloc(1)*charge(i)
               grad(2) = grad(2) + gradloc(2)*charge(i)
            endif
            if (ifhess .eq. 1) then
               hess(1) = hess(1) + hessloc(1)*charge(i)
               hess(2) = hess(2) + hessloc(2)*charge(i)
               hess(3) = hess(3) + hessloc(3)*charge(i)
            endif
         endif

         if (ifdipole .eq. 1) then
            if (ifpot .eq. 1) then
               pot = pot - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1              + gradloc(2)*dipvec(2,i))
            endif
            if (ifgrad .eq. 1) then
               grad(1) = grad(1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1              + hessloc(2)*dipvec(2,i))
               grad(2) = grad(2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1              + hessloc(3)*dipvec(2,i))
            endif
            if (ifhess .eq. 1) then
               hess(1) = hess(1) - dipstr(i)*(der3(1)*dipvec(1,i)
     1              + der3(2)*dipvec(2,i))
               hess(2) = hess(2) - dipstr(i)*(der3(2)*dipvec(1,i)
     1              + der3(3)*dipvec(2,i))
               hess(3) = hess(3) - dipstr(i)*(der3(3)*dipvec(1,i)
     1              + der3(4)*dipvec(2,i))
            endif
         endif

         if (ifquad .eq. 1) then
            if (ifpot .eq. 1) then
               pot = pot + quadstr(i)*(hessloc(1)*quadvec(1,i)
     1              + hessloc(2)*quadvec(2,i) + hessloc(3)*quadvec(3,i))
            endif
            if (ifgrad .eq. 1) then
               grad(1) = grad(1) + quadstr(i)*(der3(1)*quadvec(1,i)
     1              + der3(2)*quadvec(2,i) + der3(3)*quadvec(3,i))
               grad(2) = grad(2) + quadstr(i)*(der3(2)*quadvec(1,i)
     1              + der3(3)*quadvec(2,i) + der3(4)*quadvec(3,i))
            endif
            if (ifhess .eq. 1) then
               hess(1) = hess(1) + quadstr(i)*(der4(1)*quadvec(1,i)
     1              + der4(2)*quadvec(2,i) + der4(3)*quadvec(3,i))
               hess(2) = hess(2) + quadstr(i)*(der4(2)*quadvec(1,i)
     1              + der4(3)*quadvec(2,i) + der4(4)*quadvec(3,i))
               hess(3) = hess(3) + quadstr(i)*(der4(3)*quadvec(1,i)
     1              + der4(4)*quadvec(2,i) + der4(5)*quadvec(3,i))
            endif
         endif

         if (ifoct .eq. 1) then
            if (ifpot .eq. 1) then
               pot = pot - octstr(i)*(der3(1)*octvec(1,i)
     1              + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1              + der3(4)*octvec(4,i))
            endif
            if (ifgrad .eq. 1) then
               grad(1) = grad(1) - octstr(i)*(der4(1)*octvec(1,i)
     1              + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1              + der4(4)*octvec(4,i))
               grad(2) = grad(2) - octstr(i)*(der4(2)*octvec(1,i)
     1              + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1              + der4(5)*octvec(4,i))
            endif
            if (ifhess .eq. 1) then
               hess(1) = hess(1) - octstr(i)*(der5(1)*octvec(1,i)
     1              + der5(2)*octvec(2,i) + der5(3)*octvec(3,i)
     1              + der5(4)*octvec(4,i))
               hess(2) = hess(2) - octstr(i)*(der5(2)*octvec(1,i)
     1              + der5(3)*octvec(2,i) + der5(4)*octvec(3,i)
     1              + der5(5)*octvec(4,i))
               hess(3) = hess(3) - octstr(i)*(der5(3)*octvec(1,i)
     1              + der5(4)*octvec(2,i) + der5(5)*octvec(3,i)
     1              + der5(6)*octvec(4,i))
            endif
         endif

      enddo
      
      return
      end

      subroutine mbhpotgrad2dall_cdqo_add(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     ifoct,octstr,octvec,target,ifpot,pot,ifgrad,grad,
     3     ifhess,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*), charge(*), dipstr(*)
      real *8 dipvec(2,*), quadstr(*), quadvec(3,*), octstr(*)
      real *8 octvec(4,*)
      real *8 target(2)
      real *8 pot, grad(2), hess(3)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
      integer ifoct
c     local variables
      real *8 pottemp, gradtemp(2), hesstemp(3)

      call mbhpotgrad2dall_cdqo(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     ifoct,octstr,octvec,target,ifpot,pottemp,ifgrad,gradtemp,
     3     ifhess,hesstemp)

      if (ifpot .eq. 1) pot = pot + pottemp
      if (ifgrad .eq. 1) then
         grad(1) = grad(1) + gradtemp(1)
         grad(2) = grad(2) + gradtemp(2)
      endif
      if (ifhess .eq. 1) then
         hess(1) = hess(1) + hesstemp(1)
         hess(2) = hess(2) + hesstemp(2)
         hess(3) = hess(3) + hesstemp(3)
      endif

      return
      end
  
      subroutine mbhpotgrad2dall_cdqo3(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     ifoct,octstr,octvec,target,ifpot,pot,ifgrad,grad,
     3     ifhess,hess,ifder3,der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*), charge(*), dipstr(*)
      real *8 dipvec(2,*), quadstr(*), quadvec(3,*), octstr(*)
      real *8 octvec(4,*)
      real *8 target(2)
      real *8 pot, grad(2), hess(3), der3(4)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
      integer ifder3, ifoct
c     local variables
      integer i, ntermstemp
      real *8 pottemp, gradtemp(2), hesstemp(3), der3temp(4)
      real *8 rscale, rs2, pi, xtemp, ytemp, rtemp, tiny
      complex *16 ima, hexpmbh(0:3), hexpy(0:3)
      data ima /(0.0d0,1.0d0)/
      data tiny / 1.0d-200 /

      pi = 4.0d0*datan(1.0d0)

      if (ifpot .eq. 1) pot = 0.0d0
      if (ifgrad .eq. 1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0
      endif
      if (ifhess .eq. 1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0
      endif
      if (ifder3 .eq. 1) then
         der3(1) = 0.0d0
         der3(2) = 0.0d0
         der3(3) = 0.0d0
         der3(4) = 0.0d0
      endif

      ntermstemp = 0
      if (ifdipole .eq. 1) ntermstemp = 1
      if (ifquad .eq. 1) ntermstemp = 2
      if (ifoct .eq. 1) ntermstemp = 3

      do i = 1,ns
         
         xtemp = target(1)-source(1,i)
         ytemp = target(2)-source(2,i)

         rtemp = dsqrt(xtemp**2+ytemp**2)

         rscale = min(rtemp*beta,1.0d0)
         if (rscale .lt. tiny) return

         rs2 = rscale*rscale

         hexpmbh(0) = 0.0d0
         hexpmbh(1) = 0.0d0
         hexpmbh(2) = 0.0d0
         hexpmbh(3) = 0.0d0

         hexpy(0) = 0.0d0
         hexpy(1) = 0.0d0
         hexpy(2) = 0.0d0
         hexpy(3) = 0.0d0

         if (ifcharge .eq. 1) then
            hexpmbh(0) = hexpmbh(0) + charge(i)/(2.0d0*pi*beta**2)
         endif

         if (ifdipole .eq. 1) then
            hexpmbh(1) = hexpmbh(1) + dipstr(i)*(dipvec(1,i)+
     1           ima*dipvec(2,i))/(rscale*2.0d0*pi*beta)
         endif

         if (ifquad .eq. 1) then
            hexpmbh(2) = hexpmbh(2) + quadstr(i)*(quadvec(1,i)
     1           +quadvec(2,i)*ima
     2           -quadvec(3,i))/(rs2*4.0d0*pi)

            hexpy(0) = quadstr(i)*(quadvec(1,i)
     1           +quadvec(3,i))/(4.0d0*pi)
         endif

         if (ifoct .eq. 1) then
            hexpmbh(3) = hexpmbh(3) + beta*octstr(i)*(octvec(1,i) 
     1           + ima*octvec(2,i) - octvec(3,i) 
     2           -ima*octvec(4,i))/(rs2*rscale*8.0d0*pi)
            hexpy(1) = hexpy(1) + beta*octstr(i)*(3.0d0*octvec(1,i)
     1           + ima*octvec(2,i) + octvec(3,i) 
     2           + 3.0d0*ima*octvec(4,i))/(rscale*8.0d0*pi)
         endif

         call mbh2dmpeval3(beta,rscale,source(1,i),hexpmbh,hexpy,
     1        ntermstemp,target,pottemp,ifgrad,gradtemp,ifhess,hesstemp,
     2        ifder3,der3temp)

         if (ifpot .eq. 1) pot = pot+pottemp
         if (ifgrad .eq. 1) then
            grad(1) = grad(1)+gradtemp(1)
            grad(2) = grad(2)+gradtemp(2)
         endif
         if (ifhess .eq. 1) then
            hess(1) = hess(1)+hesstemp(1)
            hess(2) = hess(2)+hesstemp(2)
            hess(3) = hess(3)+hesstemp(3)
         endif
         if (ifder3 .eq. 1) then
            der3(1) = der3(1)+der3temp(1)
            der3(2) = der3(2)+der3temp(2)
            der3(3) = der3(3)+der3temp(3)
            der3(4) = der3(4)+der3temp(4)
         endif
         
      enddo
      
      return
      end

      subroutine mbhpotgrad2dall_cdqo3_self(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     ifoct,octstr,octvec,ifpot,pot,ifgrad,grad,
     3     ifhess,hess,ifder3,der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This performs the sum over all charges at the source locations
c     (excluding the infinite self interaction)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*), charge(*), dipstr(*)
      real *8 dipvec(2,*), quadstr(*), quadvec(3,*), octstr(*)
      real *8 octvec(4,*)
      real *8 pot(*), grad(2,*), hess(3,*), der3(4,*)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
      integer ifder3, ifoct
c     local variables
      integer i, ntermstemp, j
      real *8 pottemp, gradtemp(2), hesstemp(3), der3temp(4)
      real *8 pottemp2, gradtemp2(2), hesstemp2(3), der3temp2(4)      
      real *8 rscale, rs2, pi, xtemp, ytemp, rtemp, tiny, dzero
      complex *16 ima, hexpmbh(0:3), hexpy(0:3), zero
      complex *16 hexpmbh2(0:3), hexpy2(0:3)
      data ima /(0.0d0,1.0d0)/
      data tiny / 1.0d-200 /
      data dzero /0.0d0/
      data zero /(0.0d0,0.0d0)/

      pi = 4.0d0*datan(1.0d0)

      if (ifpot .eq. 1) then
         do i = 1,ns
            pot(i) = dzero
         enddo
      endif
      if (ifgrad .eq. 1) then
         do i = 1,ns
            grad(1,i) = dzero
            grad(2,i) = dzero
         enddo
      endif
      if (ifhess .eq. 1) then
         do i = 1,ns
            hess(1,i) = dzero
            hess(2,i) = dzero
            hess(3,i) = dzero
         enddo
      endif
      if (ifder3 .eq. 1) then
         do i = 1,ns
            der3(1,i) = dzero
            der3(2,i) = dzero
            der3(3,i) = dzero
            der3(4,i) = dzero
         enddo
      endif

      ntermstemp = 0
      if (ifdipole .eq. 1) ntermstemp = 1
      if (ifquad .eq. 1) ntermstemp = 2
      if (ifoct .eq. 1) ntermstemp = 3

      do i = 1,ns
         do j = i+1,ns
            xtemp = source(1,j)-source(1,i)
            ytemp = source(2,j)-source(2,i)
            
            rtemp = dsqrt(xtemp**2+ytemp**2)
            
            rscale = min(rtemp*beta,1.0d0)
            if (rscale .lt. tiny) return
            
            rs2 = rscale*rscale
            
            hexpmbh(0) = zero
            hexpmbh(1) = zero
            hexpmbh(2) = zero
            hexpmbh(3) = zero

            hexpy(0) = zero
            hexpy(1) = zero
            hexpy(2) = zero
            hexpy(3) = zero

            hexpmbh2(0) = zero
            hexpmbh2(1) = zero
            hexpmbh2(2) = zero
            hexpmbh2(3) = zero

            hexpy2(0) = zero
            hexpy2(1) = zero
            hexpy2(2) = zero
            hexpy2(3) = zero

            if (ifcharge .eq. 1) then
               hexpmbh(0) = hexpmbh(0) + charge(i)/(2.0d0*pi*beta**2)
               hexpmbh2(0) = hexpmbh2(0) + charge(j)/(2.0d0*pi*beta**2)
            endif

            if (ifdipole .eq. 1) then
               hexpmbh(1) = hexpmbh(1) + dipstr(i)*(dipvec(1,i)+
     1              ima*dipvec(2,i))/(rscale*2.0d0*pi*beta)
               hexpmbh2(1) = hexpmbh2(1) + dipstr(j)*(dipvec(1,j)+
     1              ima*dipvec(2,j))/(rscale*2.0d0*pi*beta)
            endif

            if (ifquad .eq. 1) then
               hexpmbh(2) = hexpmbh(2) + quadstr(i)*(quadvec(1,i)
     1              +quadvec(2,i)*ima
     2              -quadvec(3,i))/(rs2*4.0d0*pi)
               hexpmbh2(2) = hexpmbh2(2) + quadstr(j)*(quadvec(1,j)
     1              +quadvec(2,j)*ima
     2              -quadvec(3,j))/(rs2*4.0d0*pi)

               hexpy(0) = quadstr(i)*(quadvec(1,i)
     1              +quadvec(3,i))/(4.0d0*pi)
               hexpy2(0) = quadstr(j)*(quadvec(1,j)
     1              +quadvec(3,j))/(4.0d0*pi)
            endif

            if (ifoct .eq. 1) then
               hexpmbh(3) = hexpmbh(3) + beta*octstr(i)*(octvec(1,i) 
     1              + ima*octvec(2,i) - octvec(3,i) 
     2              -ima*octvec(4,i))/(rs2*rscale*8.0d0*pi)
               hexpy(1) = hexpy(1) + beta*octstr(i)*(3.0d0*octvec(1,i)
     1              + ima*octvec(2,i) + octvec(3,i) 
     2              + 3.0d0*ima*octvec(4,i))/(rscale*8.0d0*pi)
               hexpmbh2(3) = hexpmbh2(3) + beta*octstr(j)*(octvec(1,j) 
     1              + ima*octvec(2,j) - octvec(3,j) 
     2              -ima*octvec(4,j))/(rs2*rscale*8.0d0*pi)
               hexpy2(1) = hexpy2(1) + beta*octstr(j)*(3.0d0*octvec(1,j)
     1              + ima*octvec(2,j) + octvec(3,j) 
     2              + 3.0d0*ima*octvec(4,j))/(rscale*8.0d0*pi)
            endif

            call mbh2dmpeval3(beta,rscale,source(1,j),hexpmbh2,hexpy2,
     1           ntermstemp,source(1,i),pottemp,
     2           ifgrad,gradtemp,ifhess,hesstemp,
     3           ifder3,der3temp)

            call mbh2dmpeval3(beta,rscale,source(1,i),hexpmbh,hexpy,
     1           ntermstemp,source(1,j),pottemp2,
     2           ifgrad,gradtemp2,ifhess,
     3           hesstemp2,ifder3,der3temp2)

            if (ifpot .eq. 1) then
               pot(i) = pot(i)+pottemp
               pot(j) = pot(j)+pottemp2
            endif
            if (ifgrad .eq. 1) then
               grad(1,i) = grad(1,i)+gradtemp(1)
               grad(2,i) = grad(2,i)+gradtemp(2)
               grad(1,j) = grad(1,j)+gradtemp2(1)
               grad(2,j) = grad(2,j)+gradtemp2(2)
            endif
            if (ifhess .eq. 1) then
               hess(1,i) = hess(1,i)+hesstemp(1)
               hess(2,i) = hess(2,i)+hesstemp(2)
               hess(3,i) = hess(3,i)+hesstemp(3)
               hess(1,j) = hess(1,j)+hesstemp2(1)
               hess(2,j) = hess(2,j)+hesstemp2(2)
               hess(3,j) = hess(3,j)+hesstemp2(3)
            endif
            if (ifder3 .eq. 1) then
               der3(1,i) = der3(1,i)+der3temp(1)
               der3(2,i) = der3(2,i)+der3temp(2)
               der3(3,i) = der3(3,i)+der3temp(3)
               der3(4,i) = der3(4,i)+der3temp(4)
               der3(1,j) = der3(1,j)+der3temp2(1)
               der3(2,j) = der3(2,j)+der3temp2(2)
               der3(3,j) = der3(3,j)+der3temp2(3)
               der3(4,j) = der3(4,j)+der3temp2(4)
            endif
         enddo         
      enddo
      
      return
      end

      subroutine mbhpotgrad2dall_cdqo3_self_add(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     ifoct,octstr,octvec,ifpot,pot,ifgrad,grad,
     3     ifhess,hess,ifder3,der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This performs the sum over all charges at the source locations
c     (excluding the infinite self interaction)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*), charge(*), dipstr(*)
      real *8 dipvec(2,*), quadstr(*), quadvec(3,*), octstr(*)
      real *8 octvec(4,*)
      real *8 pot(*), grad(2,*), hess(3,*), der3(4,*)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
      integer ifder3, ifoct
c     local variables
      integer i, ntermstemp, j
      real *8 pottemp, gradtemp(2), hesstemp(3), der3temp(4)
      real *8 pottemp2, gradtemp2(2), hesstemp2(3), der3temp2(4)      
      real *8 rscale, rs2, pi, xtemp, ytemp, rtemp, tiny, dzero
      complex *16 ima, hexpmbh(0:3), hexpy(0:3), zero
      complex *16 hexpmbh2(0:3), hexpy2(0:3)
      data ima /(0.0d0,1.0d0)/
      data tiny / 1.0d-200 /
      data dzero /0.0d0/
      data zero /(0.0d0,0.0d0)/

      pi = 4.0d0*datan(1.0d0)

      ntermstemp = 0
      if (ifdipole .eq. 1) ntermstemp = 1
      if (ifquad .eq. 1) ntermstemp = 2
      if (ifoct .eq. 1) ntermstemp = 3

      do i = 1,ns
         do j = i+1,ns
            xtemp = source(1,j)-source(1,i)
            ytemp = source(2,j)-source(2,i)
            
            rtemp = dsqrt(xtemp**2+ytemp**2)
            
            rscale = min(rtemp*beta,1.0d0)
            if (rscale .lt. tiny) return
            
            rs2 = rscale*rscale
            
            hexpmbh(0) = zero
            hexpmbh(1) = zero
            hexpmbh(2) = zero
            hexpmbh(3) = zero

            hexpy(0) = zero
            hexpy(1) = zero
            hexpy(2) = zero
            hexpy(3) = zero

            hexpmbh2(0) = zero
            hexpmbh2(1) = zero
            hexpmbh2(2) = zero
            hexpmbh2(3) = zero

            hexpy2(0) = zero
            hexpy2(1) = zero
            hexpy2(2) = zero
            hexpy2(3) = zero

            if (ifcharge .eq. 1) then
               hexpmbh(0) = hexpmbh(0) + charge(i)/(2.0d0*pi*beta**2)
               hexpmbh2(0) = hexpmbh2(0) + charge(j)/(2.0d0*pi*beta**2)
            endif

            if (ifdipole .eq. 1) then
               hexpmbh(1) = hexpmbh(1) + dipstr(i)*(dipvec(1,i)+
     1              ima*dipvec(2,i))/(rscale*2.0d0*pi*beta)
               hexpmbh2(1) = hexpmbh2(1) + dipstr(j)*(dipvec(1,j)+
     1              ima*dipvec(2,j))/(rscale*2.0d0*pi*beta)
            endif

            if (ifquad .eq. 1) then
               hexpmbh(2) = hexpmbh(2) + quadstr(i)*(quadvec(1,i)
     1              +quadvec(2,i)*ima
     2              -quadvec(3,i))/(rs2*4.0d0*pi)
               hexpmbh2(2) = hexpmbh2(2) + quadstr(j)*(quadvec(1,j)
     1              +quadvec(2,j)*ima
     2              -quadvec(3,j))/(rs2*4.0d0*pi)

               hexpy(0) = quadstr(i)*(quadvec(1,i)
     1              +quadvec(3,i))/(4.0d0*pi)
               hexpy2(0) = quadstr(j)*(quadvec(1,j)
     1              +quadvec(3,j))/(4.0d0*pi)
            endif

            if (ifoct .eq. 1) then
               hexpmbh(3) = hexpmbh(3) + beta*octstr(i)*(octvec(1,i) 
     1              + ima*octvec(2,i) - octvec(3,i) 
     2              -ima*octvec(4,i))/(rs2*rscale*8.0d0*pi)
               hexpy(1) = hexpy(1) + beta*octstr(i)*(3.0d0*octvec(1,i)
     1              + ima*octvec(2,i) + octvec(3,i) 
     2              + 3.0d0*ima*octvec(4,i))/(rscale*8.0d0*pi)
               hexpmbh2(3) = hexpmbh2(3) + beta*octstr(j)*(octvec(1,j) 
     1              + ima*octvec(2,j) - octvec(3,j) 
     2              -ima*octvec(4,j))/(rs2*rscale*8.0d0*pi)
               hexpy2(1) = hexpy2(1) + beta*octstr(j)*(3.0d0*octvec(1,j)
     1              + ima*octvec(2,j) + octvec(3,j) 
     2              + 3.0d0*ima*octvec(4,j))/(rscale*8.0d0*pi)
            endif

            call mbh2dmpeval3(beta,rscale,source(1,j),hexpmbh2,hexpy2,
     1           ntermstemp,source(1,i),pottemp,
     2           ifgrad,gradtemp,ifhess,hesstemp,
     3           ifder3,der3temp)

            call mbh2dmpeval3(beta,rscale,source(1,i),hexpmbh,hexpy,
     1           ntermstemp,source(1,j),pottemp2,
     2           ifgrad,gradtemp2,ifhess,
     3           hesstemp2,ifder3,der3temp2)

            if (ifpot .eq. 1) then
               pot(i) = pot(i)+pottemp
               pot(j) = pot(j)+pottemp2
            endif
            if (ifgrad .eq. 1) then
               grad(1,i) = grad(1,i)+gradtemp(1)
               grad(2,i) = grad(2,i)+gradtemp(2)
               grad(1,j) = grad(1,j)+gradtemp2(1)
               grad(2,j) = grad(2,j)+gradtemp2(2)
            endif
            if (ifhess .eq. 1) then
               hess(1,i) = hess(1,i)+hesstemp(1)
               hess(2,i) = hess(2,i)+hesstemp(2)
               hess(3,i) = hess(3,i)+hesstemp(3)
               hess(1,j) = hess(1,j)+hesstemp2(1)
               hess(2,j) = hess(2,j)+hesstemp2(2)
               hess(3,j) = hess(3,j)+hesstemp2(3)
            endif
            if (ifder3 .eq. 1) then
               der3(1,i) = der3(1,i)+der3temp(1)
               der3(2,i) = der3(2,i)+der3temp(2)
               der3(3,i) = der3(3,i)+der3temp(3)
               der3(4,i) = der3(4,i)+der3temp(4)
               der3(1,j) = der3(1,j)+der3temp2(1)
               der3(2,j) = der3(2,j)+der3temp2(2)
               der3(3,j) = der3(3,j)+der3temp2(3)
               der3(4,j) = der3(4,j)+der3temp2(4)
            endif
         enddo         
      enddo
      
      return
      end

      subroutine mbhpotgrad2dall_cdqo3_add(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     ifoct,octstr,octvec,target,ifpot,pot,ifgrad,grad,
     3     ifhess,hess,ifder3,der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*), charge(*), dipstr(*)
      real *8 dipvec(2,*), quadstr(*), quadvec(3,*), octstr(*)
      real *8 target(2), octvec(4,*)
      real *8 pot, grad(2), hess(3), der3(4)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
      integer ifder3, ifoct
c     local variables
      real *8 pottemp, gradtemp(2), hesstemp(3), der3temp(4)

      call mbhpotgrad2dall_cdqo3(beta,source,ns,ifcharge,
     1     charge,ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     ifoct,octstr,octvec,target,ifpot,pottemp,ifgrad,gradtemp,
     3     ifhess,hesstemp,ifder3,der3temp)

      if (ifpot .eq. 1) pot = pot + pottemp
      if (ifgrad .eq. 1) then
         grad(1) = grad(1) + gradtemp(1)
         grad(2) = grad(2) + gradtemp(2)
      endif
      if (ifhess .eq. 1) then
         hess(1) = hess(1) + hesstemp(1)
         hess(2) = hess(2) + hesstemp(2)
         hess(3) = hess(3) + hesstemp(3)
      endif
      if (ifder3 .eq. 1) then
         der3(1) = der3(1) + der3temp(1)
         der3(2) = der3(2) + der3temp(2)
         der3(3) = der3(3) + der3temp(3)
         der3(4) = der3(4) + der3temp(4)
      endif

      return
      end
  
      subroutine mbhpotgrad2dall_cdqo3_add_exrad(beta,src,exrad,ns,
     1     ifcharge,charge,ifdipole,dipstr,dipvec,ifquad,quadstr,
     2     quadvec,ifoct,octstr,octvec,target,ifpot,pot,ifgrad,grad,
     3     ifhess,hess,ifder3,der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, src(2,*), charge(*), dipstr(*), exrad(*)
      real *8 dipvec(2,*), quadstr(*), quadvec(3,*), octstr(*)
      real *8 target(2), octvec(4,*)
      real *8 pot, grad(2), hess(3), der3(4)
      integer ns, ifcharge, ifdipole, ifquad, ifpot, ifgrad, ifhess
      integer ifder3, ifoct
c     local variables
      real *8 pottemp, gradtemp(2), hesstemp(3), der3temp(4), rs

      integer nstemp, i

      nstemp = 1
      
      do i = 1,ns
         rs = dsqrt((src(1,i)-target(1))**2 + (src(2,i)-target(2))**2)
         if (rs .gt. exrad(i)) then
            call mbhpotgrad2dall_cdqo3(beta,src(1,i),nstemp,ifcharge,
     1           charge(i),ifdipole,dipstr(i),dipvec(1,i),ifquad,
     2           quadstr(i),quadvec(1,i),
     2           ifoct,octstr(i),octvec(1,i),target,ifpot,pottemp,
     2           ifgrad,gradtemp,ifhess,hesstemp,ifder3,der3temp)

            if (ifpot .eq. 1) pot = pot + pottemp
            if (ifgrad .eq. 1) then
               grad(1) = grad(1) + gradtemp(1)
               grad(2) = grad(2) + gradtemp(2)
            endif
            if (ifhess .eq. 1) then
               hess(1) = hess(1) + hesstemp(1)
               hess(2) = hess(2) + hesstemp(2)
               hess(3) = hess(3) + hesstemp(3)
            endif
            if (ifder3 .eq. 1) then
               der3(1) = der3(1) + der3temp(1)
               der3(2) = der3(2) + der3temp(2)
               der3(3) = der3(3) + der3temp(3)
               der3(4) = der3(4) + der3temp(4)
            endif
         endif
      enddo

      return
      end
  
      subroutine mbh2dmpeval(beta,rscale,center,mbhmpole,ympole,
     1     nterms,ztarg,pot,ifgrad,grad,ifhess,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     beta               : Yukawa parameter
c     rscale             : the scaling factor
c     center(2)          : expansion center
c     mbhmpole(0:nterms) : coefficients for difference-type functions
c     ympole(nterms)     : coefficients for modified Bessel functions
c     nterms             : number of terms in expansions
c     ztarg(2)           : coordinates of target
c     ifgrad             : flag, IFGRAD = 1 means compute gradient
c     ifhess             : flag, IFHESS = 1 means compute Hessian
c
c     OUTPUT:
c
c     pot             : value of potential at target
c     grad(2)         : value of gradient at target (if requested)
c     hess(3)         : value of Hessian at target (if requested)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      complex *16 mbhmpole(0:nterms), ympole(0:nterms)
      real *8 rscale, center(2), ztarg(2), beta
      real *8 pot, grad(2), hess(3)
      integer nterms, ifgrad, ifhess
c     local
      real *8 zdiff(2), r, theta, kvec(0:200)
      real *8 diffs(0:200), ders(0:200), pih, pow(0:200), dpow(0:200)
      complex *16 hvec(0:200), hder(0:200), z, eye, ztemp1, ztemp2
      complex *16 mptemp1(0:200),mptemp2(0:200)
      complex *16 ympolex(0:200), ympoley(0:200)
      complex *16 ympolexx(0:200), ympolexy(0:200), ympoleyy(0:200)
      complex *16 mbhmpolex(0:200), mbhmpoley(0:200)
      complex *16 mbhmpolexx(0:200),mbhmpolexy(0:200),mbhmpoleyy(0:200)
      real *8 dc, ds, dc2, ds2, dcx, dsx, dcy, dsy
      integer ifder, ifders, l, j, i
      data eye /(0.0d0,1.0d0)/

      if( ifgrad .eq. 1 .or. ifhess .eq. 1 ) then
c     
         do i=0,nterms+1
            ympolex(i)=0
            ympoley(i)=0
            mbhmpolex(i)=0
            mbhmpoley(i)=0
         enddo

         ympolex(1) = -beta/rscale*ympole(0)
         ympoley(1) = -beta/rscale*ympole(0)*eye

         mbhmpolex(1) = -beta/rscale*mbhmpole(0)
         mbhmpoley(1) = -beta/rscale*mbhmpole(0)*eye

         do i=1,nterms
            ympolex(i-1) = ympolex(i-1) -beta/2.0d0*ympole(i)*rscale
            ympolex(i+1) = ympolex(i+1) -beta/2.0d0*ympole(i)/rscale
            ympoley(i-1) = ympoley(i-1) +beta/2.0d0*ympole(i)*rscale*eye
            ympoley(i+1) = ympoley(i+1) -beta/2.0d0*ympole(i)/rscale*eye

            ympolex(i-1) = ympolex(i-1) -beta/2.0d0*mbhmpole(i)*rscale
            mbhmpolex(i+1) = mbhmpolex(i+1) 
     1           -beta/2.0d0*mbhmpole(i)/rscale
            ympoley(i-1) = ympoley(i-1) 
     1           +beta/2.0d0*mbhmpole(i)*rscale*eye
            mbhmpoley(i+1) = mbhmpoley(i+1) 
     1           -beta/2.0d0*mbhmpole(i)/rscale*eye
         enddo
      endif
c     
      if( ifhess .eq. 1 ) then
c     
         do i=0,nterms+2
            ympolexx(i)=0
            ympolexy(i)=0
            ympoleyy(i)=0
            mbhmpolexx(i)=0
            mbhmpolexy(i)=0
            mbhmpoleyy(i)=0
         enddo

         ympolexx(1) = -beta/1.0d0/rscale*dreal(ympolex(0))
         ympolexy(1) = -beta/1.0d0/rscale*dreal(ympolex(0))*eye
         ympoleyy(1) = -beta/1.0d0/rscale*dreal(ympoley(0))*eye

         mbhmpolexx(1) = -beta/1.0d0/rscale*dreal(mbhmpolex(0))
         mbhmpolexy(1) = -beta/1.0d0/rscale*dreal(mbhmpolex(0))*eye
         mbhmpoleyy(1) = -beta/1.0d0/rscale*dreal(mbhmpoley(0))*eye

         do i=1,nterms+1
            ympolexx(i-1) = ympolexx(i-1) -beta/2.0d0*ympolex(i)*rscale
            ympolexx(i+1) = ympolexx(i+1) -beta/2.0d0*ympolex(i)/rscale
            ympolexy(i-1) = ympolexy(i-1) 
     1           +beta/2.0d0*ympolex(i)*rscale*eye
            ympolexy(i+1) = ympolexy(i+1) 
     1           -beta/2.0d0*ympolex(i)/rscale*eye
            ympoleyy(i-1) = ympoleyy(i-1) 
     1           +beta/2.0d0*ympoley(i)*rscale*eye
            ympoleyy(i+1) = ympoleyy(i+1) 
     1           -beta/2.0d0*ympoley(i)/rscale*eye

            ympolexx(i-1) = ympolexx(i-1) 
     1           -beta/2.0d0*mbhmpolex(i)*rscale
            mbhmpolexx(i+1) = mbhmpolexx(i+1) 
     1           -beta/2.0d0*mbhmpolex(i)/rscale
            ympolexy(i-1) = ympolexy(i-1) 
     1           +beta/2.0d0*mbhmpolex(i)*rscale*eye
            mbhmpolexy(i+1) = mbhmpolexy(i+1) 
     1           -beta/2.0d0*mbhmpolex(i)/rscale*eye
            ympoleyy(i-1) = ympoleyy(i-1) 
     1           +beta/2.0d0*mbhmpoley(i)*rscale*eye
            mbhmpoleyy(i+1) = mbhmpoleyy(i+1) 
     1           -beta/2.0d0*mbhmpoley(i)/rscale*eye
         enddo


      endif


      zdiff(1) = ztarg(1)-center(1)
      zdiff(2) = ztarg(2)-center(2)
      call h2cart2polar(zdiff,r,theta)

c     get values of difference functions
      ifders = 0
      call diffslogbk_fast(r,beta,rscale,diffs,ifders,ders,kvec,
     1     nterms+2)       

      mptemp1(0)=diffs(0)
      mptemp2(0)=kvec(0)
      ztemp2=exp(eye*theta)
      ztemp1=ztemp2
      do j=1,nterms+2
         mptemp1(j)=dcmplx(diffs(j)*dreal(ztemp1),
     1        diffs(j)*dimag(ztemp1))
         mptemp2(j)=dcmplx(kvec(j)*dreal(ztemp1),kvec(j)*dimag(ztemp1))
         ztemp1 = ztemp1*ztemp2
      enddo

c     evaluate potential and derivatives

      pot = 0.0d0

      pot = pot + dreal(mptemp1(0))*dreal(mbhmpole(0))
     1     + dreal(mptemp2(0))*dreal(ympole(0))

      do j = 1,nterms
         pot = pot + dreal(mptemp1(j))*dreal(mbhmpole(j)) +
     1        dimag(mptemp1(j))*dimag(mbhmpole(j))
         pot = pot + dreal(mptemp2(j))*dreal(ympole(j)) +
     1        dimag(mptemp2(j))*dimag(ympole(j))
      enddo

      if (ifgrad .eq. 1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0

         grad(1) = grad(1) + dreal(mptemp1(0))*dreal(mbhmpolex(0))
     1        + dreal(mptemp2(0))*dreal(ympolex(0))
         grad(2) = grad(2) + dreal(mptemp1(0))*dreal(mbhmpoley(0))
     1        + dreal(mptemp2(0))*dreal(ympoley(0))

         do j = 1,nterms+1
            grad(1) = grad(1) + dreal(mptemp1(j))*dreal(mbhmpolex(j)) +
     1           dimag(mptemp1(j))*dimag(mbhmpolex(j))
            grad(1) = grad(1) + dreal(mptemp2(j))*dreal(ympolex(j)) +
     1           dimag(mptemp2(j))*dimag(ympolex(j))
            grad(2) = grad(2) + dreal(mptemp1(j))*dreal(mbhmpoley(j)) +
     1           dimag(mptemp1(j))*dimag(mbhmpoley(j))
            grad(2) = grad(2) + dreal(mptemp2(j))*dreal(ympoley(j)) +
     1           dimag(mptemp2(j))*dimag(ympoley(j))
         enddo
      endif

      if (ifhess .eq. 1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0

         hess(1) = hess(1) + dreal(mptemp1(0))*dreal(mbhmpolexx(0))
     1        + dreal(mptemp2(0))*dreal(ympolexx(0))
         hess(2) = hess(2) + dreal(mptemp1(0))*dreal(mbhmpolexy(0))
     1        + dreal(mptemp2(0))*dreal(ympolexy(0))
         hess(3) = hess(3) + dreal(mptemp1(0))*dreal(mbhmpoleyy(0))
     1        + dreal(mptemp2(0))*dreal(ympoleyy(0))

         do j = 1,nterms+1
            hess(1) = hess(1) + dreal(mptemp1(j))*dreal(mbhmpolexx(j)) +
     1           dimag(mptemp1(j))*dimag(mbhmpolexx(j))
            hess(1) = hess(1) + dreal(mptemp2(j))*dreal(ympolexx(j)) +
     1           dimag(mptemp2(j))*dimag(ympolexx(j))
            hess(2) = hess(2) + dreal(mptemp1(j))*dreal(mbhmpolexy(j)) +
     1           dimag(mptemp1(j))*dimag(mbhmpolexy(j))
            hess(2) = hess(2) + dreal(mptemp2(j))*dreal(ympolexy(j)) +
     1           dimag(mptemp2(j))*dimag(ympolexy(j))
            hess(3) = hess(3) + dreal(mptemp1(j))*dreal(mbhmpoleyy(j)) +
     1           dimag(mptemp1(j))*dimag(mbhmpoleyy(j))
            hess(3) = hess(3) + dreal(mptemp2(j))*dreal(ympoleyy(j)) +
     1           dimag(mptemp2(j))*dimag(ympoleyy(j))
         enddo
      endif

      return 
      end

      subroutine mbh2dmpeval3(beta,rscale,center,mbhmpole,ympole,
     1     nterms,ztarg,pot,ifgrad,grad,ifhess,hess,ifder3,der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     beta               : Yukawa parameter
c     rscale             : the scaling factor
c     center(2)          : expansion center
c     mbhmpole(0:nterms) : coefficients for difference-type functions
c     ympole(nterms)     : coefficients for modified Bessel functions
c     nterms             : number of terms in expansions
c     ztarg(2)           : coordinates of target
c     ifgrad             : flag, IFGRAD = 1 means compute gradient
c     ifhess             : flag, IFHESS = 1 means compute Hessian
c
c     OUTPUT:
c
c     pot             : value of potential at target
c     grad(2)         : value of gradient at target (if requested)
c     hess(3)         : value of Hessian at target (if requested)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      complex *16 mbhmpole(0:nterms), ympole(0:nterms)
      real *8 rscale, center(2), ztarg(2), beta
      real *8 pot, grad(2), hess(3), der3(4)
      integer nterms, ifgrad, ifhess, ifder3
c     local
      real *8 zdiff(2), r, theta, kvec(0:200)
      real *8 diffs(0:200), ders(0:200), pih, pow(0:200), dpow(0:200)
      complex *16 hvec(0:200), hder(0:200), z, eye, ztemp1, ztemp2
      complex *16 mptemp1(0:200),mptemp2(0:200)
      complex *16 ympolex(0:200), ympoley(0:200)
      complex *16 ympolexx(0:200), ympolexy(0:200), ympoleyy(0:200)
      complex *16 mbhmpolex(0:200), mbhmpoley(0:200)
      complex *16 mbhmpolexx(0:200),mbhmpolexy(0:200),mbhmpoleyy(0:200)
      complex *16 mbhmpolexxx(0:200),mbhmpolexxy(0:200)
      complex *16 mbhmpolexyy(0:200), mbhmpoleyyy(0:200)
      complex *16 ympolexxx(0:200),ympolexxy(0:200)
      complex *16 ympolexyy(0:200), ympoleyyy(0:200)
      real *8 dc, ds, dc2, ds2, dcx, dsx, dcy, dsy
      integer ifder, ifders, l, j, i
      data eye /(0.0d0,1.0d0)/

      if( ifgrad .eq. 1 .or. ifhess .eq. 1 .or. ifder3 .eq. 1 ) then
c     
         do i=0,nterms+1
            ympolex(i)=0
            ympoley(i)=0
            mbhmpolex(i)=0
            mbhmpoley(i)=0
         enddo

         ympolex(1) = -beta/rscale*ympole(0)
         ympoley(1) = -beta/rscale*ympole(0)*eye

         mbhmpolex(1) = -beta/rscale*mbhmpole(0)
         mbhmpoley(1) = -beta/rscale*mbhmpole(0)*eye

         do i=1,nterms
            ympolex(i-1) = ympolex(i-1) -beta/2.0d0*ympole(i)*rscale
            ympolex(i+1) = ympolex(i+1) -beta/2.0d0*ympole(i)/rscale
            ympoley(i-1) = ympoley(i-1) +beta/2.0d0*ympole(i)*rscale*eye
            ympoley(i+1) = ympoley(i+1) -beta/2.0d0*ympole(i)/rscale*eye

            ympolex(i-1) = ympolex(i-1) -beta/2.0d0*mbhmpole(i)*rscale
            mbhmpolex(i+1) = mbhmpolex(i+1) 
     1           -beta/2.0d0*mbhmpole(i)/rscale
            ympoley(i-1) = ympoley(i-1) 
     1           +beta/2.0d0*mbhmpole(i)*rscale*eye
            mbhmpoley(i+1) = mbhmpoley(i+1) 
     1           -beta/2.0d0*mbhmpole(i)/rscale*eye
         enddo
      endif
c     
      if( ifhess .eq. 1 .or. ifder3 .eq. 1 ) then
c     
         do i=0,nterms+2
            ympolexx(i)=0
            ympolexy(i)=0
            ympoleyy(i)=0
            mbhmpolexx(i)=0
            mbhmpolexy(i)=0
            mbhmpoleyy(i)=0
         enddo

         ympolexx(1) = -beta/1.0d0/rscale*dreal(ympolex(0))
         ympolexy(1) = -beta/1.0d0/rscale*dreal(ympolex(0))*eye
         ympoleyy(1) = -beta/1.0d0/rscale*dreal(ympoley(0))*eye

         mbhmpolexx(1) = -beta/1.0d0/rscale*dreal(mbhmpolex(0))
         mbhmpolexy(1) = -beta/1.0d0/rscale*dreal(mbhmpolex(0))*eye
         mbhmpoleyy(1) = -beta/1.0d0/rscale*dreal(mbhmpoley(0))*eye

         do i=1,nterms+1
            ympolexx(i-1) = ympolexx(i-1) -beta/2.0d0*ympolex(i)*rscale
            ympolexx(i+1) = ympolexx(i+1) -beta/2.0d0*ympolex(i)/rscale
            ympolexy(i-1) = ympolexy(i-1) 
     1           +beta/2.0d0*ympolex(i)*rscale*eye
            ympolexy(i+1) = ympolexy(i+1) 
     1           -beta/2.0d0*ympolex(i)/rscale*eye
            ympoleyy(i-1) = ympoleyy(i-1) 
     1           +beta/2.0d0*ympoley(i)*rscale*eye
            ympoleyy(i+1) = ympoleyy(i+1) 
     1           -beta/2.0d0*ympoley(i)/rscale*eye

            ympolexx(i-1) = ympolexx(i-1) 
     1           -beta/2.0d0*mbhmpolex(i)*rscale
            mbhmpolexx(i+1) = mbhmpolexx(i+1) 
     1           -beta/2.0d0*mbhmpolex(i)/rscale
            ympolexy(i-1) = ympolexy(i-1) 
     1           +beta/2.0d0*mbhmpolex(i)*rscale*eye
            mbhmpolexy(i+1) = mbhmpolexy(i+1) 
     1           -beta/2.0d0*mbhmpolex(i)/rscale*eye
            ympoleyy(i-1) = ympoleyy(i-1) 
     1           +beta/2.0d0*mbhmpoley(i)*rscale*eye
            mbhmpoleyy(i+1) = mbhmpoleyy(i+1) 
     1           -beta/2.0d0*mbhmpoley(i)/rscale*eye
         enddo


      endif

      if( ifder3 .eq. 1 ) then
c     
         do i=0,nterms+3
            ympolexxx(i)=0
            ympolexxy(i)=0
            ympolexyy(i)=0
            ympoleyyy(i)=0
            mbhmpolexxx(i)=0
            mbhmpolexxy(i)=0
            mbhmpolexyy(i)=0
            mbhmpoleyyy(i)=0
         enddo

         ympolexxx(1) = -beta/1.0d0/rscale*dreal(ympolexx(0))
         ympolexxy(1) = -beta/1.0d0/rscale*dreal(ympolexy(0))
         ympolexyy(1) = -beta/1.0d0/rscale*dreal(ympolexy(0))*eye
         ympoleyyy(1) = -beta/1.0d0/rscale*dreal(ympoleyy(0))*eye

         mbhmpolexxx(1) = -beta/1.0d0/rscale*dreal(mbhmpolexx(0))
         mbhmpolexxy(1) = -beta/1.0d0/rscale*dreal(mbhmpolexy(0))
         mbhmpolexyy(1) = -beta/1.0d0/rscale*dreal(mbhmpolexy(0))*eye
         mbhmpoleyyy(1) = -beta/1.0d0/rscale*dreal(mbhmpoleyy(0))*eye

         do i=1,nterms+2
            ympolexxx(i-1) = ympolexxx(i-1) 
     1           -beta/2.0d0*ympolexx(i)*rscale
            ympolexxx(i+1) = ympolexxx(i+1) 
     1           -beta/2.0d0*ympolexx(i)/rscale
            ympolexxy(i-1) = ympolexxy(i-1) 
     1           -beta/2.0d0*ympolexy(i)*rscale
            ympolexxy(i+1) = ympolexxy(i+1) 
     1           -beta/2.0d0*ympolexy(i)/rscale
            ympolexyy(i-1) = ympolexyy(i-1) 
     1           +beta/2.0d0*ympolexy(i)*rscale*eye
            ympolexyy(i+1) = ympolexyy(i+1) 
     1           -beta/2.0d0*ympolexy(i)/rscale*eye
            ympoleyyy(i-1) = ympoleyyy(i-1) 
     1           +beta/2.0d0*ympoleyy(i)*rscale*eye
            ympoleyyy(i+1) = ympoleyyy(i+1) 
     1           -beta/2.0d0*ympoleyy(i)/rscale*eye

            ympolexxx(i-1) = ympolexxx(i-1) 
     1           -beta/2.0d0*mbhmpolexx(i)*rscale
            mbhmpolexxx(i+1) = mbhmpolexxx(i+1) 
     1           -beta/2.0d0*mbhmpolexx(i)/rscale
            ympolexxy(i-1) = ympolexxy(i-1) 
     1           -beta/2.0d0*mbhmpolexy(i)*rscale
            mbhmpolexxy(i+1) = mbhmpolexxy(i+1) 
     1           -beta/2.0d0*mbhmpolexy(i)/rscale
            ympolexyy(i-1) = ympolexyy(i-1) 
     1           +beta/2.0d0*mbhmpolexy(i)*rscale*eye
            mbhmpolexyy(i+1) = mbhmpolexyy(i+1) 
     1           -beta/2.0d0*mbhmpolexy(i)/rscale*eye
            ympoleyyy(i-1) = ympoleyyy(i-1) 
     1           +beta/2.0d0*mbhmpoleyy(i)*rscale*eye
            mbhmpoleyyy(i+1) = mbhmpoleyyy(i+1) 
     1           -beta/2.0d0*mbhmpoleyy(i)/rscale*eye
         enddo


      endif

      zdiff(1) = ztarg(1)-center(1)
      zdiff(2) = ztarg(2)-center(2)
      call h2cart2polar(zdiff,r,theta)

c     get values of difference functions
      ifders = 0
      call diffslogbk_fast(r,beta,rscale,diffs,ifders,ders,kvec,
     1     nterms+3)       

      mptemp1(0)=diffs(0)
      mptemp2(0)=kvec(0)
      ztemp2=exp(eye*theta)
      ztemp1=ztemp2
      do j=1,nterms+3
         mptemp1(j)=dcmplx(diffs(j)*dreal(ztemp1),
     1        diffs(j)*dimag(ztemp1))
         mptemp2(j)=dcmplx(kvec(j)*dreal(ztemp1),kvec(j)*dimag(ztemp1))
         ztemp1 = ztemp1*ztemp2
      enddo

c     evaluate potential and derivatives

      pot = 0.0d0

      pot = pot + dreal(mptemp1(0))*dreal(mbhmpole(0))
     1     + dreal(mptemp2(0))*dreal(ympole(0))

      do j = 1,nterms
         pot = pot + dreal(mptemp1(j))*dreal(mbhmpole(j)) +
     1        dimag(mptemp1(j))*dimag(mbhmpole(j))
         pot = pot + dreal(mptemp2(j))*dreal(ympole(j)) +
     1        dimag(mptemp2(j))*dimag(ympole(j))
      enddo

      if (ifgrad .eq. 1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0

         grad(1) = grad(1) + dreal(mptemp1(0))*dreal(mbhmpolex(0))
     1        + dreal(mptemp2(0))*dreal(ympolex(0))
         grad(2) = grad(2) + dreal(mptemp1(0))*dreal(mbhmpoley(0))
     1        + dreal(mptemp2(0))*dreal(ympoley(0))

         do j = 1,nterms+1
            grad(1) = grad(1) + dreal(mptemp1(j))*dreal(mbhmpolex(j)) +
     1           dimag(mptemp1(j))*dimag(mbhmpolex(j))
            grad(1) = grad(1) + dreal(mptemp2(j))*dreal(ympolex(j)) +
     1           dimag(mptemp2(j))*dimag(ympolex(j))
            grad(2) = grad(2) + dreal(mptemp1(j))*dreal(mbhmpoley(j)) +
     1           dimag(mptemp1(j))*dimag(mbhmpoley(j))
            grad(2) = grad(2) + dreal(mptemp2(j))*dreal(ympoley(j)) +
     1           dimag(mptemp2(j))*dimag(ympoley(j))
         enddo
      endif

      if (ifhess .eq. 1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0

         hess(1) = hess(1) + dreal(mptemp1(0))*dreal(mbhmpolexx(0))
     1        + dreal(mptemp2(0))*dreal(ympolexx(0))
         hess(2) = hess(2) + dreal(mptemp1(0))*dreal(mbhmpolexy(0))
     1        + dreal(mptemp2(0))*dreal(ympolexy(0))
         hess(3) = hess(3) + dreal(mptemp1(0))*dreal(mbhmpoleyy(0))
     1        + dreal(mptemp2(0))*dreal(ympoleyy(0))

         do j = 1,nterms+2
            hess(1) = hess(1) + dreal(mptemp1(j))*dreal(mbhmpolexx(j)) +
     1           dimag(mptemp1(j))*dimag(mbhmpolexx(j))
            hess(1) = hess(1) + dreal(mptemp2(j))*dreal(ympolexx(j)) +
     1           dimag(mptemp2(j))*dimag(ympolexx(j))
            hess(2) = hess(2) + dreal(mptemp1(j))*dreal(mbhmpolexy(j)) +
     1           dimag(mptemp1(j))*dimag(mbhmpolexy(j))
            hess(2) = hess(2) + dreal(mptemp2(j))*dreal(ympolexy(j)) +
     1           dimag(mptemp2(j))*dimag(ympolexy(j))
            hess(3) = hess(3) + dreal(mptemp1(j))*dreal(mbhmpoleyy(j)) +
     1           dimag(mptemp1(j))*dimag(mbhmpoleyy(j))
            hess(3) = hess(3) + dreal(mptemp2(j))*dreal(ympoleyy(j)) +
     1           dimag(mptemp2(j))*dimag(ympoleyy(j))
         enddo
      endif

      if (ifder3 .eq. 1) then
         der3(1) = 0.0d0
         der3(2) = 0.0d0
         der3(3) = 0.0d0
         der3(4) = 0.0d0

         der3(1) = der3(1) + dreal(mptemp1(0))*dreal(mbhmpolexxx(0))
     1        + dreal(mptemp2(0))*dreal(ympolexxx(0))
         der3(2) = der3(2) + dreal(mptemp1(0))*dreal(mbhmpolexxy(0))
     1        + dreal(mptemp2(0))*dreal(ympolexxy(0))
         der3(3) = der3(3) + dreal(mptemp1(0))*dreal(mbhmpolexyy(0))
     1        + dreal(mptemp2(0))*dreal(ympolexyy(0))
         der3(4) = der3(4) + dreal(mptemp1(0))*dreal(mbhmpoleyyy(0))
     1        + dreal(mptemp2(0))*dreal(ympoleyyy(0))

         do j = 1,nterms+3
            der3(1) = der3(1) + dreal(mptemp1(j))*dreal(mbhmpolexxx(j))+
     1           dimag(mptemp1(j))*dimag(mbhmpolexxx(j))
            der3(1) = der3(1) + dreal(mptemp2(j))*dreal(ympolexxx(j)) +
     1           dimag(mptemp2(j))*dimag(ympolexxx(j))
            der3(2) = der3(2) + dreal(mptemp1(j))*dreal(mbhmpolexxy(j))+
     1           dimag(mptemp1(j))*dimag(mbhmpolexxy(j))
            der3(2) = der3(2) + dreal(mptemp2(j))*dreal(ympolexxy(j)) +
     1           dimag(mptemp2(j))*dimag(ympolexxy(j))
            der3(3) = der3(3) + dreal(mptemp1(j))*dreal(mbhmpolexyy(j))+
     1           dimag(mptemp1(j))*dimag(mbhmpolexyy(j))
            der3(3) = der3(3) + dreal(mptemp2(j))*dreal(ympolexyy(j)) +
     1           dimag(mptemp2(j))*dimag(ympolexyy(j))
            der3(4) = der3(4) + dreal(mptemp1(j))*dreal(mbhmpoleyyy(j))+
     1           dimag(mptemp1(j))*dimag(mbhmpoleyyy(j))
            der3(4) = der3(4) + dreal(mptemp2(j))*dreal(ympoleyyy(j)) +
     1           dimag(mptemp2(j))*dimag(ympoleyyy(j))
         enddo
      endif

      return 
      end

      subroutine mbh2dtaeval(beta,rscale,center,mbhloc,lloc,
     1     nterms,ztarg,pot,ifgrad,grad,ifhess,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     beta               : Yukawa parameter
c     rscale             : the scaling factor
c     center(2)          : expansion center
c     mbhloc(0:nterms)   : coefficients for difference-type functions
c     lloc(0:nterms)     : coefficients for Taylor series
c     nterms             : number of terms in expansions
c     ztarg(2)           : coordinates of target
c     ifgrad             : flag, IFGRAD = 1 means compute gradient
c     ifhess             : flag, IFHESS = 1 means compute Hessian
c
c     OUTPUT:
c
c     pot             : value of potential at target
c     grad(2)         : value of gradient at target (if requested)
c     hess(3)         : value of Hessian at target (if requested)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      complex *16 mbhloc(0:nterms), lloc(0:nterms)
      real *8 rscale, center(2), ztarg(2), beta
      real *8 pot, grad(2), hess(3)
      integer nterms, ifgrad, ifhess
c     local
      real *8 zdiff(2), r, theta, ivec(0:200)
      real *8 diffs(0:200), ders(0:200), pih, pow(0:200), dpow(0:200)
      complex *16 hvec(0:200), hder(0:200), z, eye, ztemp1, ztemp2
      complex *16 mptemp1(0:200),mptemp2(0:200)
      complex *16 llocx(0:200), llocy(0:200)
      complex *16 llocxx(0:200), llocxy(0:200), llocyy(0:200)
      complex *16 mbhlocx(0:200), mbhlocy(0:200)
      complex *16 mbhlocxx(0:200),mbhlocxy(0:200),mbhlocyy(0:200)
      real *8 dtemp
      integer ifder, ifders, l, j, i
      data eye /(0.0d0,1.0d0)/

      if (ifgrad .eq. 1 .or. ifhess .eq. 1) then
         do i=0,nterms+1
            mbhlocx(i)=0
            mbhlocy(i)=0
            llocx(i)=0
            llocy(i)=0
         enddo

         dtemp = 1.0d0/2.0d0

         mbhlocx(1) = dreal(mbhloc(0))*beta*rscale
         mbhlocy(1) = dreal(mbhloc(0))*beta*rscale*eye
         llocx(1) = dreal(mbhloc(0))*beta*dtemp*rscale
         llocy(1) = dreal(mbhloc(0))*beta*dtemp*rscale*eye

         do i=1,nterms
            dtemp = dtemp/(2*(i+1))
            mbhlocx(i-1)=mbhlocx(i-1)+beta/2/rscale*mbhloc(i)
            mbhlocx(i+1)=mbhlocx(i+1)+beta/2*rscale*mbhloc(i)
            mbhlocy(i-1)=mbhlocy(i-1)-beta/2*(eye)/rscale*mbhloc(i)
            mbhlocy(i+1)=mbhlocy(i+1)+beta/2*(eye)*rscale*mbhloc(i)
            
            llocx(i-1)=llocx(i-1)+beta*i/rscale*lloc(i)
            llocy(i-1)=llocy(i-1)-beta*i*(eye)/rscale*lloc(i)

            llocx(i+1)=llocx(i+1)+beta/2*rscale*mbhloc(i)*dtemp
            llocy(i+1)=llocy(i+1)+beta/2*(eye)*rscale*mbhloc(i)*dtemp
         enddo
      endif
c
      if (ifhess.eq.1) then
         do i=0,nterms+2
            mbhlocxx(i)=0
            mbhlocxy(i)=0
            mbhlocyy(i)=0
            llocxx(i)=0
            llocxy(i)=0
            llocyy(i)=0
         enddo

         dtemp = 1.0d0/2.0d0

         mbhlocxx(1) = dreal(mbhlocx(0))*beta*rscale
         mbhlocxy(1) = dreal(mbhlocx(0))*beta*rscale*eye
         mbhlocyy(1) = dreal(mbhlocy(0))*beta*rscale*eye
         llocxx(1) = dreal(mbhlocx(0))*beta*dtemp*rscale
         llocxy(1) = dreal(mbhlocx(0))*beta*dtemp*rscale*eye
         llocyy(1) = dreal(mbhlocy(0))*beta*dtemp*rscale*eye

         do i=1,nterms+1
            dtemp = dtemp/(2*(i+1))
            mbhlocxx(i-1)=mbhlocxx(i-1)+beta/2/rscale*mbhlocx(i)
            mbhlocxx(i+1)=mbhlocxx(i+1)+beta/2*rscale*mbhlocx(i)
            mbhlocxy(i-1)=mbhlocxy(i-1)-beta/2*(eye)/rscale*mbhlocx(i)
            mbhlocxy(i+1)=mbhlocxy(i+1)+beta/2*(eye)*rscale*mbhlocx(i)
            mbhlocyy(i-1)=mbhlocyy(i-1)-beta/2*(eye)/rscale*mbhlocy(i)
            mbhlocyy(i+1)=mbhlocyy(i+1)+beta/2*(eye)*rscale*mbhlocy(i)
            
            llocxx(i-1)=llocxx(i-1)+beta*i/rscale*llocx(i)
            llocxy(i-1)=llocxy(i-1)-beta*i*(eye)/rscale*llocx(i)
            llocyy(i-1)=llocyy(i-1)-beta*i*(eye)/rscale*llocy(i)

            llocxx(i+1)=llocxx(i+1)+beta/2*rscale*mbhlocx(i)*dtemp
            llocxy(i+1)=llocxy(i+1)+beta/2*(eye)*rscale*mbhlocx(i)*dtemp
            llocyy(i+1)=llocyy(i+1)+beta/2*(eye)*rscale*mbhlocy(i)*dtemp
         enddo
      endif


      zdiff(1) = ztarg(1)-center(1)
      zdiff(2) = ztarg(2)-center(2)
      call h2cart2polar(zdiff,r,theta)

c     get values of (beta*r)^k
      call mbh2d_rk(pow,dpow,r,beta,rscale,nterms+2)

c     get values of difference functions
      ifders = 0
      call diffszkik_fast(r,beta,rscale,diffs,ifders,ders,ivec,
     1     nterms+2)

      mptemp1(0)=diffs(0)
      mptemp2(0)=pow(0)
      ztemp2=exp(eye*theta)
      ztemp1=ztemp2
      do j=1,nterms+2
         mptemp1(j)=dcmplx(diffs(j)*dreal(ztemp1),
     1        diffs(j)*dimag(ztemp1))
         mptemp2(j)=dcmplx(pow(j)*dreal(ztemp1),pow(j)*dimag(ztemp1))
         ztemp1 = ztemp1*ztemp2
      enddo

c     evaluate potential and derivatives

      pot = 0.0d0

      pot = pot + dreal(mptemp1(0))*dreal(mbhloc(0))
     1     + dreal(mptemp2(0))*dreal(lloc(0))

      do j = 1,nterms
         pot = pot + dreal(mptemp1(j))*dreal(mbhloc(j)) +
     1        dimag(mptemp1(j))*dimag(mbhloc(j))
         pot = pot + dreal(mptemp2(j))*dreal(lloc(j)) +
     1        dimag(mptemp2(j))*dimag(lloc(j))
      enddo

      if (ifgrad .eq. 1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0

         grad(1) = grad(1) + dreal(mptemp1(0))*dreal(mbhlocx(0))
     1        + dreal(mptemp2(0))*dreal(llocx(0))
         grad(2) = grad(2) + dreal(mptemp1(0))*dreal(mbhlocy(0))
     1        + dreal(mptemp2(0))*dreal(llocy(0))

         do j = 1,nterms+1
            grad(1) = grad(1) + dreal(mptemp1(j))*dreal(mbhlocx(j)) +
     1           dimag(mptemp1(j))*dimag(mbhlocx(j))
            grad(1) = grad(1) + dreal(mptemp2(j))*dreal(llocx(j)) +
     1           dimag(mptemp2(j))*dimag(llocx(j))
            grad(2) = grad(2) + dreal(mptemp1(j))*dreal(mbhlocy(j)) +
     1           dimag(mptemp1(j))*dimag(mbhlocy(j))
            grad(2) = grad(2) + dreal(mptemp2(j))*dreal(llocy(j)) +
     1           dimag(mptemp2(j))*dimag(llocy(j))
         enddo
      endif

      if (ifhess .eq. 1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0

         hess(1) = hess(1) + dreal(mptemp1(0))*dreal(mbhlocxx(0))
     1        + dreal(mptemp2(0))*dreal(llocxx(0))
         hess(2) = hess(2) + dreal(mptemp1(0))*dreal(mbhlocxy(0))
     1        + dreal(mptemp2(0))*dreal(llocxy(0))
         hess(3) = hess(3) + dreal(mptemp1(0))*dreal(mbhlocyy(0))
     1        + dreal(mptemp2(0))*dreal(llocyy(0))

         do j = 1,nterms+1
            hess(1) = hess(1) + dreal(mptemp1(j))*dreal(mbhlocxx(j)) +
     1           dimag(mptemp1(j))*dimag(mbhlocxx(j))
            hess(1) = hess(1) + dreal(mptemp2(j))*dreal(llocxx(j)) +
     1           dimag(mptemp2(j))*dimag(llocxx(j))
            hess(2) = hess(2) + dreal(mptemp1(j))*dreal(mbhlocxy(j)) +
     1           dimag(mptemp1(j))*dimag(mbhlocxy(j))
            hess(2) = hess(2) + dreal(mptemp2(j))*dreal(llocxy(j)) +
     1           dimag(mptemp2(j))*dimag(llocxy(j))
            hess(3) = hess(3) + dreal(mptemp1(j))*dreal(mbhlocyy(j)) +
     1           dimag(mptemp1(j))*dimag(mbhlocyy(j))
            hess(3) = hess(3) + dreal(mptemp2(j))*dreal(llocyy(j)) +
     1           dimag(mptemp2(j))*dimag(llocyy(j))
         enddo
      endif

      return 
      end

      subroutine mbh2dtaeval3(beta,rscale,center,mbhloc,lloc,
     1     nterms,ztarg,pot,ifgrad,grad,ifhess,hess,ifder3,der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     beta               : Yukawa parameter
c     rscale             : the scaling factor
c     center(2)          : expansion center
c     mbhloc(0:nterms)   : coefficients for difference-type functions
c     lloc(0:nterms)     : coefficients for Taylor series
c     nterms             : number of terms in expansions
c     ztarg(2)           : coordinates of target
c     ifgrad             : flag, IFGRAD = 1 means compute gradient
c     ifhess             : flag, IFHESS = 1 means compute Hessian
c
c     OUTPUT:
c
c     pot             : value of potential at target
c     grad(2)         : value of gradient at target (if requested)
c     hess(3)         : value of Hessian at target (if requested)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      complex *16 mbhloc(0:nterms), lloc(0:nterms)
      real *8 rscale, center(2), ztarg(2), beta
      real *8 pot, grad(2), hess(3), der3(4)
      integer nterms, ifgrad, ifhess, ifder3
c     local
      real *8 zdiff(2), r, theta, ivec(0:200)
      real *8 diffs(0:200), ders(0:200), pih, pow(0:200), dpow(0:200)
      complex *16 hvec(0:200), hder(0:200), z, eye, ztemp1, ztemp2
      complex *16 mptemp1(0:200),mptemp2(0:200)
      complex *16 llocx(0:200), llocy(0:200)
      complex *16 llocxx(0:200), llocxy(0:200), llocyy(0:200)
      complex *16 mbhlocx(0:200), mbhlocy(0:200)
      complex *16 mbhlocxx(0:200),mbhlocxy(0:200),mbhlocyy(0:200)
      complex *16 mbhlocxxx(0:200),mbhlocxxy(0:200),mbhlocxyy(0:200)
      complex *16 mbhlocyyy(0:200)
      complex *16 llocxxx(0:200),llocxxy(0:200),llocxyy(0:200)
      complex *16 llocyyy(0:200)
      real *8 dtemp
      integer ifder, ifders, l, j, i
      data eye /(0.0d0,1.0d0)/

      if (ifgrad .eq. 1 .or. ifhess .eq. 1 .or. ifder3 .eq. 1) then
         do i=0,nterms+1
            mbhlocx(i)=0
            mbhlocy(i)=0
            llocx(i)=0
            llocy(i)=0
         enddo

         dtemp = 1.0d0/2.0d0

         mbhlocx(1) = dreal(mbhloc(0))*beta*rscale
         mbhlocy(1) = dreal(mbhloc(0))*beta*rscale*eye
         llocx(1) = dreal(mbhloc(0))*beta*dtemp*rscale
         llocy(1) = dreal(mbhloc(0))*beta*dtemp*rscale*eye

         do i=1,nterms
            dtemp = dtemp/(2*(i+1))
            mbhlocx(i-1)=mbhlocx(i-1)+beta/2/rscale*mbhloc(i)
            mbhlocx(i+1)=mbhlocx(i+1)+beta/2*rscale*mbhloc(i)
            mbhlocy(i-1)=mbhlocy(i-1)-beta/2*(eye)/rscale*mbhloc(i)
            mbhlocy(i+1)=mbhlocy(i+1)+beta/2*(eye)*rscale*mbhloc(i)
            
            llocx(i-1)=llocx(i-1)+beta*i/rscale*lloc(i)
            llocy(i-1)=llocy(i-1)-beta*i*(eye)/rscale*lloc(i)

            llocx(i+1)=llocx(i+1)+beta/2*rscale*mbhloc(i)*dtemp
            llocy(i+1)=llocy(i+1)+beta/2*(eye)*rscale*mbhloc(i)*dtemp
         enddo
      endif
c
      if (ifhess.eq.1 .or. ifder3 .eq. 1) then
         do i=0,nterms+2
            mbhlocxx(i)=0
            mbhlocxy(i)=0
            mbhlocyy(i)=0
            llocxx(i)=0
            llocxy(i)=0
            llocyy(i)=0
         enddo

         dtemp = 1.0d0/2.0d0

         mbhlocxx(1) = dreal(mbhlocx(0))*beta*rscale
         mbhlocxy(1) = dreal(mbhlocx(0))*beta*rscale*eye
         mbhlocyy(1) = dreal(mbhlocy(0))*beta*rscale*eye
         llocxx(1) = dreal(mbhlocx(0))*beta*dtemp*rscale
         llocxy(1) = dreal(mbhlocx(0))*beta*dtemp*rscale*eye
         llocyy(1) = dreal(mbhlocy(0))*beta*dtemp*rscale*eye

         do i=1,nterms+1
            dtemp = dtemp/(2*(i+1))
            mbhlocxx(i-1)=mbhlocxx(i-1)+beta/2/rscale*mbhlocx(i)
            mbhlocxx(i+1)=mbhlocxx(i+1)+beta/2*rscale*mbhlocx(i)
            mbhlocxy(i-1)=mbhlocxy(i-1)-beta/2*(eye)/rscale*mbhlocx(i)
            mbhlocxy(i+1)=mbhlocxy(i+1)+beta/2*(eye)*rscale*mbhlocx(i)
            mbhlocyy(i-1)=mbhlocyy(i-1)-beta/2*(eye)/rscale*mbhlocy(i)
            mbhlocyy(i+1)=mbhlocyy(i+1)+beta/2*(eye)*rscale*mbhlocy(i)
            
            llocxx(i-1)=llocxx(i-1)+beta*i/rscale*llocx(i)
            llocxy(i-1)=llocxy(i-1)-beta*i*(eye)/rscale*llocx(i)
            llocyy(i-1)=llocyy(i-1)-beta*i*(eye)/rscale*llocy(i)

            llocxx(i+1)=llocxx(i+1)+beta/2*rscale*mbhlocx(i)*dtemp
            llocxy(i+1)=llocxy(i+1)+beta/2*(eye)*rscale*mbhlocx(i)*dtemp
            llocyy(i+1)=llocyy(i+1)+beta/2*(eye)*rscale*mbhlocy(i)*dtemp
         enddo
      endif

      if (ifder3 .eq. 1) then
         do i=0,nterms+3
            mbhlocxxx(i)=0
            mbhlocxxy(i)=0
            mbhlocxyy(i)=0
            mbhlocyyy(i)=0
            llocxxx(i)=0
            llocxxy(i)=0
            llocxyy(i)=0
            llocyyy(i)=0
         enddo

         dtemp = 1.0d0/2.0d0

         mbhlocxxx(1) = dreal(mbhlocxx(0))*beta*rscale
         mbhlocxxy(1) = dreal(mbhlocxy(0))*beta*rscale
         mbhlocxyy(1) = dreal(mbhlocxy(0))*beta*rscale*eye
         mbhlocyyy(1) = dreal(mbhlocyy(0))*beta*rscale*eye
         llocxxx(1) = dreal(mbhlocxx(0))*beta*dtemp*rscale
         llocxxy(1) = dreal(mbhlocxy(0))*beta*dtemp*rscale
         llocxyy(1) = dreal(mbhlocxy(0))*beta*dtemp*rscale*eye
         llocyyy(1) = dreal(mbhlocyy(0))*beta*dtemp*rscale*eye

         do i=1,nterms+2
            dtemp = dtemp/(2*(i+1))
            mbhlocxxx(i-1)=mbhlocxxx(i-1)+beta/2/rscale*mbhlocxx(i)
            mbhlocxxx(i+1)=mbhlocxxx(i+1)+beta/2*rscale*mbhlocxx(i)
            mbhlocxxy(i-1)=mbhlocxxy(i-1)+beta/2/rscale*mbhlocxy(i)
            mbhlocxxy(i+1)=mbhlocxxy(i+1)+beta/2*rscale*mbhlocxy(i)
            mbhlocxyy(i-1)=mbhlocxyy(i-1)
     1           -beta/2*(eye)/rscale*mbhlocxy(i)
            mbhlocxyy(i+1)=mbhlocxyy(i+1)
     1           +beta/2*(eye)*rscale*mbhlocxy(i)
            mbhlocyyy(i-1)=mbhlocyyy(i-1)
     1           -beta/2*(eye)/rscale*mbhlocyy(i)
            mbhlocyyy(i+1)=mbhlocyyy(i+1)
     1           +beta/2*(eye)*rscale*mbhlocyy(i)
            
            llocxxx(i-1)=llocxxx(i-1)+beta*i/rscale*llocxx(i)
            llocxxy(i-1)=llocxxy(i-1)+beta*i/rscale*llocxy(i)
            llocxyy(i-1)=llocxyy(i-1)-beta*i*(eye)/rscale*llocxy(i)
            llocyyy(i-1)=llocyyy(i-1)-beta*i*(eye)/rscale*llocyy(i)

            llocxxx(i+1)=llocxxx(i+1)+beta/2*rscale*mbhlocxx(i)*dtemp
            llocxxy(i+1)=llocxxy(i+1)+beta/2*rscale*mbhlocxy(i)*dtemp
            llocxyy(i+1)=llocxyy(i+1)
     1           +beta/2*(eye)*rscale*mbhlocxy(i)*dtemp
            llocyyy(i+1)=llocyyy(i+1)
     1           +beta/2*(eye)*rscale*mbhlocyy(i)*dtemp
         enddo
      endif


      zdiff(1) = ztarg(1)-center(1)
      zdiff(2) = ztarg(2)-center(2)
      call h2cart2polar(zdiff,r,theta)

c     get values of (beta*r)^k
      call mbh2d_rk(pow,dpow,r,beta,rscale,nterms+3)

c     get values of difference functions
      ifders = 0
      call diffszkik_fast(r,beta,rscale,diffs,ifders,ders,ivec,
     1     nterms+3)

      mptemp1(0)=diffs(0)
      mptemp2(0)=pow(0)
      ztemp2=exp(eye*theta)
      ztemp1=ztemp2
      do j=1,nterms+3
         mptemp1(j)=dcmplx(diffs(j)*dreal(ztemp1),
     1        diffs(j)*dimag(ztemp1))
         mptemp2(j)=dcmplx(pow(j)*dreal(ztemp1),pow(j)*dimag(ztemp1))
         ztemp1 = ztemp1*ztemp2
      enddo

c     evaluate potential and derivatives

      pot = 0.0d0

      pot = pot + dreal(mptemp1(0))*dreal(mbhloc(0))
     1     + dreal(mptemp2(0))*dreal(lloc(0))

      do j = 1,nterms
         pot = pot + dreal(mptemp1(j))*dreal(mbhloc(j)) +
     1        dimag(mptemp1(j))*dimag(mbhloc(j))
         pot = pot + dreal(mptemp2(j))*dreal(lloc(j)) +
     1        dimag(mptemp2(j))*dimag(lloc(j))
      enddo

      if (ifgrad .eq. 1) then
         grad(1) = 0.0d0
         grad(2) = 0.0d0

         grad(1) = grad(1) + dreal(mptemp1(0))*dreal(mbhlocx(0))
     1        + dreal(mptemp2(0))*dreal(llocx(0))
         grad(2) = grad(2) + dreal(mptemp1(0))*dreal(mbhlocy(0))
     1        + dreal(mptemp2(0))*dreal(llocy(0))

         do j = 1,nterms+1
            grad(1) = grad(1) + dreal(mptemp1(j))*dreal(mbhlocx(j)) +
     1           dimag(mptemp1(j))*dimag(mbhlocx(j))
            grad(1) = grad(1) + dreal(mptemp2(j))*dreal(llocx(j)) +
     1           dimag(mptemp2(j))*dimag(llocx(j))
            grad(2) = grad(2) + dreal(mptemp1(j))*dreal(mbhlocy(j)) +
     1           dimag(mptemp1(j))*dimag(mbhlocy(j))
            grad(2) = grad(2) + dreal(mptemp2(j))*dreal(llocy(j)) +
     1           dimag(mptemp2(j))*dimag(llocy(j))
         enddo
      endif

      if (ifhess .eq. 1) then
         hess(1) = 0.0d0
         hess(2) = 0.0d0
         hess(3) = 0.0d0

         hess(1) = hess(1) + dreal(mptemp1(0))*dreal(mbhlocxx(0))
     1        + dreal(mptemp2(0))*dreal(llocxx(0))
         hess(2) = hess(2) + dreal(mptemp1(0))*dreal(mbhlocxy(0))
     1        + dreal(mptemp2(0))*dreal(llocxy(0))
         hess(3) = hess(3) + dreal(mptemp1(0))*dreal(mbhlocyy(0))
     1        + dreal(mptemp2(0))*dreal(llocyy(0))

         do j = 1,nterms+2
            hess(1) = hess(1) + dreal(mptemp1(j))*dreal(mbhlocxx(j)) +
     1           dimag(mptemp1(j))*dimag(mbhlocxx(j))
            hess(1) = hess(1) + dreal(mptemp2(j))*dreal(llocxx(j)) +
     1           dimag(mptemp2(j))*dimag(llocxx(j))
            hess(2) = hess(2) + dreal(mptemp1(j))*dreal(mbhlocxy(j)) +
     1           dimag(mptemp1(j))*dimag(mbhlocxy(j))
            hess(2) = hess(2) + dreal(mptemp2(j))*dreal(llocxy(j)) +
     1           dimag(mptemp2(j))*dimag(llocxy(j))
            hess(3) = hess(3) + dreal(mptemp1(j))*dreal(mbhlocyy(j)) +
     1           dimag(mptemp1(j))*dimag(mbhlocyy(j))
            hess(3) = hess(3) + dreal(mptemp2(j))*dreal(llocyy(j)) +
     1           dimag(mptemp2(j))*dimag(llocyy(j))
         enddo
      endif

      if (ifder3 .eq. 1) then
         der3(1) = 0.0d0
         der3(2) = 0.0d0
         der3(3) = 0.0d0
         der3(4) = 0.0d0

         der3(1) = der3(1) + dreal(mptemp1(0))*dreal(mbhlocxxx(0))
     1        + dreal(mptemp2(0))*dreal(llocxxx(0))
         der3(2) = der3(2) + dreal(mptemp1(0))*dreal(mbhlocxxy(0))
     1        + dreal(mptemp2(0))*dreal(llocxxy(0))
         der3(3) = der3(3) + dreal(mptemp1(0))*dreal(mbhlocxyy(0))
     1        + dreal(mptemp2(0))*dreal(llocxyy(0))
         der3(4) = der3(4) + dreal(mptemp1(0))*dreal(mbhlocyyy(0))
     1        + dreal(mptemp2(0))*dreal(llocyyy(0))

         do j = 1,nterms+3
            der3(1) = der3(1) + dreal(mptemp1(j))*dreal(mbhlocxxx(j)) +
     1           dimag(mptemp1(j))*dimag(mbhlocxxx(j))
            der3(1) = der3(1) + dreal(mptemp2(j))*dreal(llocxxx(j)) +
     1           dimag(mptemp2(j))*dimag(llocxxx(j))
            der3(2) = der3(2) + dreal(mptemp1(j))*dreal(mbhlocxxy(j)) +
     1           dimag(mptemp1(j))*dimag(mbhlocxxy(j))
            der3(2) = der3(2) + dreal(mptemp2(j))*dreal(llocxxy(j)) +
     1           dimag(mptemp2(j))*dimag(llocxxy(j))
            der3(3) = der3(3) + dreal(mptemp1(j))*dreal(mbhlocxyy(j)) +
     1           dimag(mptemp1(j))*dimag(mbhlocxyy(j))
            der3(3) = der3(3) + dreal(mptemp2(j))*dreal(llocxyy(j)) +
     1           dimag(mptemp2(j))*dimag(llocxyy(j))
            der3(4) = der3(4) + dreal(mptemp1(j))*dreal(mbhlocyyy(j)) +
     1           dimag(mptemp1(j))*dimag(mbhlocyyy(j))
            der3(4) = der3(4) + dreal(mptemp2(j))*dreal(llocyyy(j)) +
     1           dimag(mptemp2(j))*dimag(llocyyy(j))
         enddo
      endif

      return 
      end

      subroutine mbh2dmpevalall(beta,rscale,center,mbhmpole,ympole,
     1     nterms,ztarg,ntarg,ifpot,pot,ifgrad,grad,ifhess,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     beta               : Yukawa parameter
c     rscale             : the scaling factor
c     center(2)          : expansion center
c     mbhmpole(0:nterms) : coefficients for difference-type functions
c     ympole(nterms)     : coefficients for modified Bessel functions
c     nterms             : number of terms in expansions
c     ztarg(2)           : coordinates of target
c     ifgrad             : flag, IFGRAD = 1 means compute gradient
c     ifhess             : flag, IFHESS = 1 means compute Hessian
c
c     OUTPUT:
c
c     pot             : value of potential at target
c     grad(2)         : value of gradient at target 
c                       (if requested)
c     hess(3)         : value of Hessian at target 
c                       (if requested)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      complex *16 mbhmpole(0:nterms), ympole(0:nterms)
      real *8 rscale, center(2), ztarg(2,*), beta
      real *8 pot(*), grad(2,*), hess(3,*)
      integer nterms, ifpot, ifgrad, ifhess, ntarg
c     local
      integer iii
      real *8 potloc, gradloc(2), hessloc(3)

      do iii = 1,ntarg
         call mbh2dmpeval(beta,rscale,center,mbhmpole,ympole,
     1     nterms,ztarg(1,iii),potloc,ifgrad,gradloc,ifhess,hessloc)
         if (ifpot .eq. 1) pot(iii) = potloc
         if (ifgrad .eq. 1) then
            grad(1,iii) = gradloc(1)
            grad(2,iii) = gradloc(2)
         endif
         if (ifhess .eq. 1) then
            hess(1,iii) = hessloc(1)
            hess(2,iii) = hessloc(2)
            hess(3,iii) = hessloc(3)
         endif
      enddo

      return
      end
     
      subroutine mbh2dmpeval3all(beta,rscale,center,mbhmpole,ympole,
     1     nterms,ztarg,ntarg,ifpot,pot,ifgrad,grad,ifhess,hess,
     2     ifder3,der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     beta               : Yukawa parameter
c     rscale             : the scaling factor
c     center(2)          : expansion center
c     mbhmpole(0:nterms) : coefficients for difference-type functions
c     ympole(nterms)     : coefficients for modified Bessel functions
c     nterms             : number of terms in expansions
c     ztarg(2)           : coordinates of target
c     ifgrad             : flag, IFGRAD = 1 means compute gradient
c     ifhess             : flag, IFHESS = 1 means compute Hessian
c
c     OUTPUT:
c
c     pot             : value of potential at target
c     grad(2)         : value of gradient at target 
c                       (if requested)
c     hess(3)         : value of Hessian at target 
c                       (if requested)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      complex *16 mbhmpole(0:nterms), ympole(0:nterms)
      real *8 rscale, center(2), ztarg(2,*), beta
      real *8 pot(*), grad(2,*), hess(3,*), der3(4,*)
      integer nterms, ifpot, ifgrad, ifhess, ntarg, ifder3
c     local
      integer iii
      real *8 potloc, gradloc(2), hessloc(3), der3loc(4)

      do iii = 1,ntarg
         call mbh2dmpeval3(beta,rscale,center,mbhmpole,ympole,
     1        nterms,ztarg(1,iii),potloc,ifgrad,gradloc,ifhess,hessloc,
     2        ifder3,der3loc)
         if (ifpot .eq. 1) pot(iii) = potloc
         if (ifgrad .eq. 1) then
            grad(1,iii) = gradloc(1)
            grad(2,iii) = gradloc(2)
         endif
         if (ifhess .eq. 1) then
            hess(1,iii) = hessloc(1)
            hess(2,iii) = hessloc(2)
            hess(3,iii) = hessloc(3)
         endif
         if (ifder3 .eq. 1) then
            der3(1,iii) = der3loc(1)
            der3(2,iii) = der3loc(2)
            der3(3,iii) = der3loc(3)
            der3(4,iii) = der3loc(4)
         endif
      enddo

      return
      end
     
      subroutine mbh2dtaevalall(beta,rscale,center,mbhloc,lloc,
     1     nterms,ztarg,ntarg,ifpot,pot,ifgrad,grad,ifhess,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     beta               : Yukawa parameter
c     rscale             : the scaling factor
c     center(2)          : expansion center
c     mbhmpole(0:nterms) : coefficients for difference-type functions
c     ympole(nterms)     : coefficients for modified Bessel functions
c     nterms             : number of terms in expansions
c     ztarg(2)           : coordinates of target
c     ifgrad             : flag, IFGRAD = 1 means compute gradient
c     ifhess             : flag, IFHESS = 1 means compute Hessian
c
c     OUTPUT:
c
c     pot             : incremented value of potential at target
c     grad(2)         : incremented value of gradient at target 
c                       (if requested)
c     hess(3)         : incremented value of Hessian at target 
c                       (if requested)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      complex *16 mbhloc(0:nterms), lloc(0:nterms)
      real *8 rscale, center(2), ztarg(2,*), beta
      real *8 pot(*), grad(2,*), hess(3,*)
      integer nterms, ifpot, ifgrad, ifhess, ntarg
c     local
      integer iii
      real *8 potloc, gradloc(2), hessloc(3)

      do iii = 1,ntarg
         call mbh2dtaeval(beta,rscale,center,mbhloc,lloc,
     1     nterms,ztarg(1,iii),potloc,ifgrad,gradloc,ifhess,hessloc)
         if (ifpot .eq. 1) pot(iii) = pot(iii)+potloc
         if (ifgrad .eq. 1) then
            grad(1,iii) = grad(1,iii)+gradloc(1)
            grad(2,iii) = grad(2,iii)+gradloc(2)
         endif
         if (ifhess .eq. 1) then
            hess(1,iii) = hess(1,iii)+hessloc(1)
            hess(2,iii) = hess(2,iii)+hessloc(2)
            hess(3,iii) = hess(3,iii)+hessloc(3)
         endif
      enddo

      return
      end
     
      subroutine mbh2dtaeval3all(beta,rscale,center,mbhloc,lloc,
     1     nterms,ztarg,ntarg,ifpot,pot,ifgrad,grad,ifhess,hess,
     2     ifder3,der3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     beta               : Yukawa parameter
c     rscale             : the scaling factor
c     center(2)          : expansion center
c     mbhmpole(0:nterms) : coefficients for difference-type functions
c     ympole(nterms)     : coefficients for modified Bessel functions
c     nterms             : number of terms in expansions
c     ztarg(2)           : coordinates of target
c     ifgrad             : flag, IFGRAD = 1 means compute gradient
c     ifhess             : flag, IFHESS = 1 means compute Hessian
c
c     OUTPUT:
c
c     pot             : incremented value of potential at target
c     grad(2)         : incremented value of gradient at target 
c                       (if requested)
c     hess(3)         : incremented value of Hessian at target 
c                       (if requested)
c     der3(4)         : incremented value of 3rd order ders at target
c                       (if requested)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      complex *16 mbhloc(0:nterms), lloc(0:nterms)
      real *8 rscale, center(2), ztarg(2,*), beta
      real *8 pot(*), grad(2,*), hess(3,*), der3(4,*)
      integer nterms, ifpot, ifgrad, ifhess, ifder3, ntarg
c     local
      integer iii
      real *8 potloc, gradloc(2), hessloc(3), der3loc(4)

      do iii = 1,ntarg
         call mbh2dtaeval3(beta,rscale,center,mbhloc,lloc,
     1        nterms,ztarg(1,iii),potloc,ifgrad,gradloc,ifhess,hessloc, 
     2        ifder3,der3loc)
         if (ifpot .eq. 1) pot(iii) = pot(iii)+potloc
         if (ifgrad .eq. 1) then
            grad(1,iii) = grad(1,iii)+gradloc(1)
            grad(2,iii) = grad(2,iii)+gradloc(2)
         endif
         if (ifhess .eq. 1) then
            hess(1,iii) = hess(1,iii)+hessloc(1)
            hess(2,iii) = hess(2,iii)+hessloc(2)
            hess(3,iii) = hess(3,iii)+hessloc(3)
         endif
         if (ifder3 .eq. 1) then
            der3(1,iii) = der3(1,iii)+der3loc(1)
            der3(2,iii) = der3(2,iii)+der3loc(2)
            der3(3,iii) = der3(3,iii)+der3loc(3)
            der3(4,iii) = der3(4,iii)+der3loc(4)
         endif
      enddo

      return
      end
     

      subroutine mbh2dformmp_qp(ier,beta,rscale,source,quadstr,quadvec,
     1     ns,center,nterms,mbhmpole,ympole)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     rscale          : the scaling factor
c     source(2,ns)    : coordinates of sources
c     quadstr(ns)     : quadrupole strengths
c     quadvec(3,ns)   : quadrupole vector
c     ns              : number of sources
c     center(2)       : expansion center
c     rad             : radius of circle where field values are matched
c     zcirc(2,nterms) : precomputed evenly-spaced points on unit circle
c     beta            : the modified biharmonic parameter
c     nterms          : order of multipole expansion

c     work(*)         : work array. recommended length 16*nterms + 500
c
c     OUTPUT:
c
c     ier             : error flag. for no error, ier = 0
c     coeffs1         : coefficients for difference-type functions
c     coeffs2         : coefficients for modified Bessel functions
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      integer ier, ns, nterms
      real *8 rscale, source(2,*), quadstr(*), quadvec(3,*), center(2)
      real *8 beta
      complex *16 mbhmpole(0:nterms), ympole(0:nterms)
c     local

c     local variables
      integer iii, nterms1, i
      complex *16 ima,hexpmbh(0:2),hexpy(0:2)
      complex *16 mbhtemp(0:100), ytemp(0:100)
      real *8 pi, rs2
c     
      data ima/(0.0d0,1.0d0)/
c     

      pi = 4.0d0*datan(1.0d0)

      do i = 0,nterms
         mbhmpole(i) = 0
         ympole(i) = 0
      enddo


      do iii = 1,ns

         rs2 = rscale*rscale 

         hexpmbh(0) = 0.0d0
         hexpmbh(1) = 0.0d0
         hexpmbh(2) = quadstr(iii)*(quadvec(1,iii)
     1                          +quadvec(2,iii)*ima
     2                          -quadvec(3,iii))/(rs2*4.0d0*pi)

         hexpy(0) = quadstr(iii)*(quadvec(1,iii)
     1        +quadvec(3,iii))/(4.0d0*pi)
         hexpy(1) = 0.0d0
         hexpy(2) = 0.0d0

         nterms1 = 2

         call mbh2dmpmp(beta,rscale,source(1,iii),hexpmbh,hexpy,
     1        nterms1,rscale,center,mbhtemp,ytemp,nterms)


         do i = 0,nterms
            mbhmpole(i) = mbhmpole(i) + mbhtemp(i)
            ympole(i) = ympole(i) + ytemp(i)
         enddo

        
      enddo
      

      return
      end
      
      subroutine mbh2dformmp_all(ier,beta,rscale,src,ifcharge,charge,
     1     ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,ifoct,octstr,
     2     octvec,ns,center,nterms,mbhmpole,ympole)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      integer ier, ns, nterms, ifcharge, ifdipole, ifquad, ifoct
      real *8 beta, rscale, src(2,*), charge(*), dipstr(*), dipvec(2,*)
      real *8 quadstr(*), quadvec(3,*), octstr(*), octvec(4,*),center(2)
      complex *16 mbhmpole(0:nterms), ympole(0:nterms)
c     local variables
      integer iii, nterms1, i
      complex *16 ima,hexpmbh(0:3),hexpy(0:3)
      complex *16 mbhtemp(0:100), ytemp(0:100)
      real *8 pi, rs2
c     
      data ima/(0.0d0,1.0d0)/
c     

      pi = 4.0d0*datan(1.0d0)

      do i = 0,nterms
         mbhmpole(i) = 0
         ympole(i) = 0
      enddo

      nterms1 = 0
      if (ifdipole .eq. 1) nterms1 = 1
      if (ifquad .eq. 1) nterms1 = 2
      if (ifoct .eq. 1) nterms1 = 3

      rs2 = rscale*rscale 

      do iii = 1,ns

         hexpmbh(0) = 0.0d0
         hexpmbh(1) = 0.0d0
         hexpmbh(2) = 0.0d0
         hexpmbh(3) = 0.0d0

         hexpy(0) = 0.0d0
         hexpy(1) = 0.0d0
         hexpy(2) = 0.0d0
         hexpy(3) = 0.0d0

         if (ifcharge .eq. 1) then
            hexpmbh(0) = hexpmbh(0) + charge(iii)/(2.0d0*pi*beta**2)
         endif

         if (ifdipole .eq. 1) then
            hexpmbh(1) = hexpmbh(1) + dipstr(iii)*(dipvec(1,iii)+
     1           ima*dipvec(2,iii))/(rscale*2.0d0*pi*beta)
         endif

         if (ifquad .eq. 1) then
            hexpmbh(2) = hexpmbh(2) + quadstr(iii)*(quadvec(1,iii)
     1           +quadvec(2,iii)*ima
     2           -quadvec(3,iii))/(rs2*4.0d0*pi)

            hexpy(0) = quadstr(iii)*(quadvec(1,iii)
     1           +quadvec(3,iii))/(4.0d0*pi)
         endif

         if (ifoct .eq. 1) then
            hexpmbh(3) = hexpmbh(3) + beta*octstr(iii)*(octvec(1,iii) 
     1           + ima*octvec(2,iii) - octvec(3,iii) 
     2           -ima*octvec(4,iii))/(rs2*rscale*8.0d0*pi)
            hexpy(1) = hexpy(1) + beta*octstr(iii)*(3.0d0*octvec(1,iii)
     1           + ima*octvec(2,iii) + octvec(3,iii) 
     2           + 3.0d0*ima*octvec(4,iii))/(rscale*8.0d0*pi)
         endif

         call mbh2dmpmp(beta,rscale,src(1,iii),hexpmbh,hexpy,
     1        nterms1,rscale,center,mbhtemp,ytemp,nterms)


         do i = 0,nterms
            mbhmpole(i) = mbhmpole(i) + mbhtemp(i)
            ympole(i) = ympole(i) + ytemp(i)
         enddo

        
      enddo
      

      return
      end

      subroutine mbh2dformta_all(ier,beta,rscale,src,ifcharge,charge,
     1     ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,ifoct,octstr,
     2     octvec,ns,center,nterms,mbhmpole,ympole)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      integer ier, ns, nterms, ifcharge, ifdipole, ifquad, ifoct
      real *8 beta, rscale, src(2,*), charge(*), dipstr(*), dipvec(2,*)
      real *8 quadstr(*), quadvec(3,*), octstr(*), octvec(4,*),center(2)
      complex *16 mbhmpole(0:nterms), ympole(0:nterms)
c     local variables
      integer iii, nterms1, i
      complex *16 ima,hexpmbh(0:3),hexpy(0:3)
      complex *16 mbhtemp(0:100), ytemp(0:100)
      real *8 pi, rs2
c     
      data ima/(0.0d0,1.0d0)/
c     

      pi = 4.0d0*datan(1.0d0)

      do i = 0,nterms
         mbhmpole(i) = 0
         ympole(i) = 0
      enddo

      nterms1 = 0
      if (ifdipole .eq. 1) nterms1 = 1
      if (ifquad .eq. 1) nterms1 = 2
      if (ifoct .eq. 1) nterms1 = 3

      rs2 = rscale*rscale 

      do iii = 1,ns

         hexpmbh(0) = 0.0d0
         hexpmbh(1) = 0.0d0
         hexpmbh(2) = 0.0d0
         hexpmbh(3) = 0.0d0

         hexpy(0) = 0.0d0
         hexpy(1) = 0.0d0
         hexpy(2) = 0.0d0
         hexpy(3) = 0.0d0

         if (ifcharge .eq. 1) then
            hexpmbh(0) = hexpmbh(0) + charge(iii)/(2.0d0*pi*beta**2)
         endif

         if (ifdipole .eq. 1) then
            hexpmbh(1) = hexpmbh(1) + dipstr(iii)*(dipvec(1,iii)+
     1           ima*dipvec(2,iii))/(rscale*2.0d0*pi*beta)
         endif

         if (ifquad .eq. 1) then
            hexpmbh(2) = hexpmbh(2) + quadstr(iii)*(quadvec(1,iii)
     1           +quadvec(2,iii)*ima
     2           -quadvec(3,iii))/(rs2*4.0d0*pi)

            hexpy(0) = quadstr(iii)*(quadvec(1,iii)
     1           +quadvec(3,iii))/(4.0d0*pi)
         endif

         if (ifoct .eq. 1) then
            hexpmbh(3) = hexpmbh(3) + beta*octstr(iii)*(octvec(1,iii) 
     1           + ima*octvec(2,iii) - octvec(3,iii) 
     2           -ima*octvec(4,iii))/(rs2*rscale*8.0d0*pi)
            hexpy(1) = hexpy(1) + beta*octstr(iii)*(3.0d0*octvec(1,iii)
     1           + ima*octvec(2,iii) + octvec(3,iii) 
     2           + 3.0d0*ima*octvec(4,iii))/(rscale*8.0d0*pi)
         endif

         call mbh2dmploc(beta,rscale,src(1,iii),hexpmbh,hexpy,
     1        nterms1,rscale,center,mbhtemp,ytemp,nterms)


         do i = 0,nterms
            mbhmpole(i) = mbhmpole(i) + mbhtemp(i)
            ympole(i) = ympole(i) + ytemp(i)
         enddo

        
      enddo
      

      return
      end


      subroutine mbh2dformta_qp(ier,beta,rscale,source,quadstr,quadvec,
     1     ns,center,nterms,mbhloc,lloc)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     rscale          : the scaling factor
c     source(2,ns)    : coordinates of sources
c     quadstr(ns)     : quadrupole strengths
c     quadvec(3,ns)   : quadrupole vector
c     ns              : number of sources
c     center(2)       : expansion center
c     rad             : radius of circle where field values are matched
c     zcirc(2,nterms) : precomputed evenly-spaced points on unit circle
c     beta            : the modified biharmonic parameter
c     nterms          : order of multipole expansion

c     work(*)         : work array. recommended length 16*nterms + 500
c
c     OUTPUT:
c
c     ier             : error flag. for no error, ier = 0
c     coeffs1         : coefficients for difference-type functions
c     coeffs2         : coefficients for modified Bessel functions

c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      integer ier, ns, nterms
      real *8 rscale, source(2,*), quadstr(*), quadvec(3,*), center(2)
      real *8 beta
      complex *16 mbhloc(0:nterms), lloc(0:nterms)
c     local

c     local variables
      integer iii, nterms1, i
      complex *16 ima,hexpmbh(0:2),hexpy(0:2)
      complex *16 mbhtemp(0:100), ytemp(0:100)
      real *8 pi, rs2
c     
      data ima/(0.0d0,1.0d0)/
c     

      pi = 4.0d0*datan(1.0d0)

      do i = 0,nterms
         mbhloc(i) = 0
         lloc(i) = 0
      enddo


      do iii = 1,ns

         rs2 = rscale*rscale 

         hexpmbh(0) = 0.0d0
         hexpmbh(1) = 0.0d0
         hexpmbh(2) = quadstr(iii)*(quadvec(1,iii)
     1                          +quadvec(2,iii)*ima
     2                          -quadvec(3,iii))/(rs2*4.0d0*pi)

         hexpy(0) = quadstr(iii)*(quadvec(1,iii)
     1        +quadvec(3,iii))/(4.0d0*pi)
         hexpy(1) = 0.0d0
         hexpy(2) = 0.0d0

         nterms1 = 2

         call mbh2dmploc(beta,rscale,source(1,iii),hexpmbh,hexpy,
     1        nterms1,rscale,center,mbhtemp,ytemp,nterms)


         do i = 0,nterms
            mbhloc(i) = mbhloc(i) + mbhtemp(i)
            lloc(i) = lloc(i) + ytemp(i)
         enddo
        
      enddo
      

      return
      end

      


      subroutine mbh2dmpmp(beta,rscale1,center1,mbhmpole1,ympole1,
     1     nterms1,rscale2,center2,mbhmpole2,ympole2,nterms2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     rscale1           : scaling for original expansion
c     center1           : center for original expansion
c     coeffs1_1(nterms1): original coeffs, difference-type functions
c     coeffs2_1(nterms1): original coeffs, modified Bessel functions
c     nterms1           : number of terms in original expansion
c     rscale2           : scaling for new expansion
c     center2           : center for new expansion
c     rad2              : radius of circle to sample on 
c                         for new expansion
c     nterms2           : number of terms in original expansion
c     zcirc2(2,nterms2) : precomp'd evenly-spaced points on unit circle
c     beta              : the modified biharmonic parameter

c     work(*)           : work array. recommended length 16*nterms+500
c
c     OUTPUT:
c
c     coeffs1_2(nterms2): new coeffs, difference-type functions
c     coeffs2_2(nterms2): new coeffs, modified Bessel functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 rscale1, center1(2), beta
      real *8 rscale2, center2(2)
      complex *16 mbhmpole1(0:nterms1), ympole1(0:nterms1)
      complex *16 mbhmpole2(0:nterms2), ympole2(0:nterms2)
      integer nterms1, nterms2
c     local
      integer i, j, ifder, nterms, m, l
      real *8 zdiff(2), r, theta, pi
      real *8 diffs(0:200), ders(0:200), ival(0:200), dfac(0:200)
      real *8 pow(0:200), dpow(0:200), dfac2(0:200)
      real *8 carray(0:200,0:200), twojm1
      real *8 rsj, rsi5, rsi7, rsi, fs2, rtemp
      complex *16 jtemp(0:200), difftemp(0:200), powtemp(0:200)
      complex *16 zk, z, ima, ztemp1, zmul, ztemp
      data ima /(0.0d0,1.0d0)/

      pi = 4.0d0*datan(1.0d0)
      zk = ima*beta

      nterms = nterms1+nterms2

c     
      do i = 0,nterms2
         mbhmpole2(i) = 0
         ympole2(i) = 0
      enddo


      dfac(0) = 1.0d0
      dfac2(0) = 1.0d0
      do i = 1,nterms
         dfac(i) = dfac(i-1)*i
         dfac2(i) = dfac2(i-1)*2.0d0/i
      enddo

      do l = 0,nterms
         carray(l,0) = 1.0d0
      enddo
      do m=1,nterms
         carray(m,m) = 1.0d0
         do l=m+1,nterms
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
         enddo
      enddo
c     

c     get local difference and BesselI vals
      
      zdiff(1)=center2(1)-center1(1)
      zdiff(2)=center2(2)-center1(2)
      call h2cart2polar(zdiff,r,theta)
      theta=theta-pi
      z=zk*r
      ifder=0

      call diffszkik_fast(r,beta,rscale1,diffs,ifder,ders,ival,
     1     nterms)

      rtemp = r
      call mbh2d_rksc(pow,dpow,rtemp,beta,rscale1,nterms)

c     
      jtemp(0) = ival(0)
      difftemp(0) = diffs(0)
      powtemp(0) = pow(0)
      zmul=exp(ima*theta)
      ztemp1= zmul

      do j = 1,nterms
         jtemp( j) = ztemp1*ival(j)
         difftemp(j) = ztemp1*diffs(j)
         powtemp(j) = ztemp1*pow(j)
         ztemp1= ztemp1*zmul
      enddo

c     shift the expansion

      rsj = rscale1

      ympole2(0) = ympole2(0) + dreal(ympole1(0))*jtemp(0)
      ympole2(0) = ympole2(0) + dreal(mbhmpole1(0))*difftemp(0)
      mbhmpole2(0) = mbhmpole2(0) + dreal(mbhmpole1(0))

      do j = 1,nterms1
         ympole2(0) = ympole2(0)
     1        +(dreal(ympole1(j))*dreal(jtemp(j))
     1        +dimag(ympole1(j))*dimag(jtemp(j)))*rsj**2
         ympole2(0) = ympole2(0)
     1        +(dreal(mbhmpole1(j))*dreal(jtemp(j))
     1        +dimag(mbhmpole1(j))*dimag(jtemp(j)))*rsj**2
         rsj=rsj*rscale1
      enddo

c     
      rsi=rscale1
      rsi7=rscale2
      rsi5=rscale1/rscale2
      do i = 1,nterms2
         ympole2(i) = ympole2(i) 
     1        + 2.0d0*(dreal(ympole1(0))*dreal(jtemp(i))
     2        +ima*dreal(ympole1(0))*dimag(jtemp(i)))*rsi5
         ympole2(i) = ympole2(i) 
     1        + 2.0d0*(dreal(mbhmpole1(0))*dreal(difftemp(i))
     2        +ima*dreal(mbhmpole1(0))*dimag(difftemp(i)))*rsi5
         mbhmpole2(i) = mbhmpole2(i) 
     1        + 2.0d0*(dreal(mbhmpole1(0))*dreal(powtemp(i))
     2        +ima*dreal(mbhmpole1(0))*dimag(powtemp(i)))*rsi5
         rsj=rscale1
         twojm1 = 1.0d0
         do j = 1,min(nterms1,i)
            ympole2(i) = ympole2(i)
     1           +(dreal(ympole1(j))*dreal(jtemp(i-j))
     1           -dimag(ympole1(j))*dimag(jtemp(i-j))
     2           +ima*(dreal(ympole1(j))*dimag(jtemp(i-j))
     3           +dimag(ympole1(j))*dreal(jtemp(i-j))))*rsi5
            ympole2(i) = ympole2(i)
     1           +(dreal(ympole1(j))*dreal(jtemp(i+j))
     1           +dimag(ympole1(j))*dimag(jtemp(i+j))
     2           +ima*(dreal(ympole1(j))*dimag(jtemp(i+j))
     3           -dimag(ympole1(j))*dreal(jtemp(i+j))))*rsj**2*rsi5
            mbhmpole2(i) = mbhmpole2(i)
     1           +(dreal(mbhmpole1(j))*dreal(powtemp(i-j))
     1           -dimag(mbhmpole1(j))*dimag(powtemp(i-j))
     2           +ima*(dreal(mbhmpole1(j))*dimag(powtemp(i-j))
     3           +dimag(mbhmpole1(j))*dreal(powtemp(i-j))))
     5           *rsi5
            ztemp = difftemp(i-j)
            ympole2(i) = ympole2(i)
     1           +(dreal(mbhmpole1(j))*dreal(ztemp)
     1           -dimag(mbhmpole1(j))*dimag(ztemp)
     2           +ima*(dreal(mbhmpole1(j))*dimag(ztemp)
     3           +dimag(mbhmpole1(j))*dreal(ztemp)))*rsi5
            ympole2(i) = ympole2(i)
     1           +(dreal(mbhmpole1(j))*dreal(jtemp(i+j))
     1           +dimag(mbhmpole1(j))*dimag(jtemp(i+j))
     2           +ima*(dreal(mbhmpole1(j))*dimag(jtemp(i+j))
     3           -dimag(mbhmpole1(j))*dreal(jtemp(i+j))))
     5           *rsj**2*rsi5
            rsj=rsj*rscale1
            twojm1 = twojm1*2.0d0
         enddo
         rsj=rscale1**(i+1)
         fs2=rsi5*rscale1**2
         do j = i+1,nterms1
            ympole2(i) = ympole2(i)
     1           +(dreal(mbhmpole1(j))*dreal(jtemp(j-i))
     1           +dimag(mbhmpole1(j))*dimag(jtemp(j-i))
     2           -ima*(dreal(mbhmpole1(j))*dimag(jtemp(j-i))
     3           -dimag(mbhmpole1(j))*dreal(jtemp(j-i))))*fs2
            ympole2(i) = ympole2(i)
     1           +(dreal(mbhmpole1(j))*dreal(jtemp(i+j))
     1           +dimag(mbhmpole1(j))*dimag(jtemp(i+j))
     2           +ima*(dreal(mbhmpole1(j))*dimag(jtemp(i+j))
     3           -dimag(mbhmpole1(j))*dreal(jtemp(i+j))))
     5           *rsj**2*rsi5
            ympole2(i) = ympole2(i)
     1           +(dreal(ympole1(j))*dreal(jtemp(j-i))
     1           +dimag(ympole1(j))*dimag(jtemp(j-i))
     2           -ima*(dreal(ympole1(j))*dimag(jtemp(j-i))
     3           -dimag(ympole1(j))*dreal(jtemp(j-i))))*fs2
            ympole2(i) = ympole2(i)
     1           +(dreal(ympole1(j))*dreal(jtemp(i+j))
     1           +dimag(ympole1(j))*dimag(jtemp(i+j))
     2           +ima*(dreal(ympole1(j))*dimag(jtemp(i+j))
     3           -dimag(ympole1(j))*dreal(jtemp(i+j))))
     5           *rsj**2*rsi5
            rsj=rsj*rscale1
            fs2=fs2*rscale1**2
         enddo
         rsi=rsi*rscale1
         rsi7=rsi7*rscale2
         rsi5=rsi5*rscale1/rscale2
      enddo

      return
      end

      subroutine mbh2dmpmp_pre(beta,rscale1,center1,nterms1,
     1     rscale2,center2,nterms2,dfac,dfac2,carray,mcarray,
     2     ival,diffs,pow)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     rscale1           : scaling for original expansion
c     center1           : center for original expansion
c     coeffs1_1(nterms1): original coeffs, difference-type functions
c     coeffs2_1(nterms1): original coeffs, modified Bessel functions
c     nterms1           : number of terms in original expansion
c     rscale2           : scaling for new expansion
c     center2           : center for new expansion
c     rad2              : radius of circle to sample on 
c                         for new expansion
c     nterms2           : number of terms in original expansion
c     zcirc2(2,nterms2) : precomp'd evenly-spaced points on unit circle
c     beta              : the modified biharmonic parameter

c     work(*)           : work array. recommended length 16*nterms+500
c
c     OUTPUT:
c
c     coeffs1_2(nterms2): new coeffs, difference-type functions
c     coeffs2_2(nterms2): new coeffs, modified Bessel functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 rscale1, center1(2), beta
      real *8 rscale2, center2(2)
      real *8 pow(0:1), diffs(0:1), ival(0:1)
      real *8 carray(0:mcarray,0:mcarray)
      real *8 dfac(0:1), dfac2(0:1)
      integer nterms1, nterms2, mcarray
c     local
      integer i, j, ifder, nterms, m, l
      real *8 zdiff(2), r, theta, pi
      real *8 ders(0:200)
      real *8 dpow(0:200)
      real *8 twojm1
      real *8 rsj, rsi5, rsi7, rsi, fs2, rtemp
      complex *16 zk, z, ima, ztemp1, zmul, ztemp
      data ima /(0.0d0,1.0d0)/

      pi = 4.0d0*datan(1.0d0)
      zk = ima*beta

      nterms = nterms1+nterms2

      dfac(0) = 1.0d0
      dfac2(0) = 1.0d0
      do i = 1,nterms
         dfac(i) = dfac(i-1)*i
         dfac2(i) = dfac2(i-1)*2.0d0/i
      enddo

      do l = 0,nterms
         carray(l,0) = 1.0d0
      enddo
      do m=1,nterms
         carray(m,m) = 1.0d0
         do l=m+1,nterms
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
         enddo
      enddo
c     

c     get local difference and BesselI vals
      
      zdiff(1)=center2(1)-center1(1)
      zdiff(2)=center2(2)-center1(2)
      call h2cart2polar(zdiff,r,theta)
      theta=theta-pi
      z=zk*r
      ifder=0

      call diffszkik_fast(r,beta,rscale1,diffs,ifder,ders,ival,
     1     nterms)

      rtemp = r
      call mbh2d_rksc(pow,dpow,rtemp,beta,rscale1,nterms)

      return
      end

      subroutine mbh2dmpmp_wpre(beta,rscale1,center1,mbhmpole1,ympole1,
     1     nterms1,rscale2,center2,mbhmpole2,ympole2,nterms2,dfac,dfac2,
     2     carray,mcarray,ival,diffs,pow)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     rscale1           : scaling for original expansion
c     center1           : center for original expansion
c     coeffs1_1(nterms1): original coeffs, difference-type functions
c     coeffs2_1(nterms1): original coeffs, modified Bessel functions
c     nterms1           : number of terms in original expansion
c     rscale2           : scaling for new expansion
c     center2           : center for new expansion
c     rad2              : radius of circle to sample on 
c                         for new expansion
c     nterms2           : number of terms in original expansion
c     zcirc2(2,nterms2) : precomp'd evenly-spaced points on unit circle
c     beta              : the modified biharmonic parameter

c     work(*)           : work array. recommended length 16*nterms+500
c
c     OUTPUT:
c
c     coeffs1_2(nterms2): new coeffs, difference-type functions
c     coeffs2_2(nterms2): new coeffs, modified Bessel functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 rscale1, center1(2), beta
      real *8 rscale2, center2(2)
      real *8 pow(0:1), diffs(0:1), ival(0:1)
      real *8 carray(0:mcarray,0:mcarray)
      real *8 dfac(0:1), dfac2(0:1)
      complex *16 mbhmpole1(0:nterms1), ympole1(0:nterms1)
      complex *16 mbhmpole2(0:nterms2), ympole2(0:nterms2)
      integer nterms1, nterms2, mcarray
c     local
      integer i, j, ifder, nterms, m, l
      real *8 zdiff(2), r, theta, pi
      real *8 rsj, rsi5, rsi7, rsi, fs2, rtemp, twojm1
      complex *16 jtemp(0:200), difftemp(0:200), powtemp(0:200)
      complex *16 zk, z, ima, ztemp1, zmul, ztemp
      data ima /(0.0d0,1.0d0)/

      pi = 4.0d0*datan(1.0d0)
      zk = ima*beta

      nterms = nterms1+nterms2

c     zero out entries
      
      do i = 0,nterms2
         mbhmpole2(i) = 0.0d0
         ympole2(i) = 0.0d0
      enddo

c     get local difference and BesselI vals
      
      zdiff(1)=center2(1)-center1(1)
      zdiff(2)=center2(2)-center1(2)
      call h2cart2polar(zdiff,r,theta)
      theta=theta-pi
      z=zk*r
      ifder=0
c     
      jtemp(0) = ival(0)
      difftemp(0) = diffs(0)
      powtemp(0) = pow(0)
      zmul=exp(ima*theta)
      ztemp1= zmul

      do j = 1,nterms
         jtemp( j) = ztemp1*ival(j)
         difftemp(j) = ztemp1*diffs(j)
         powtemp(j) = ztemp1*pow(j)
         ztemp1= ztemp1*zmul
      enddo

c     shift the expansion

      rsj = rscale1

      ympole2(0) = ympole2(0) + dreal(ympole1(0))*jtemp(0)
      ympole2(0) = ympole2(0) + dreal(mbhmpole1(0))*difftemp(0)
      mbhmpole2(0) = mbhmpole2(0) + dreal(mbhmpole1(0))

      do j = 1,nterms1
         ympole2(0) = ympole2(0)
     1        +(dreal(ympole1(j))*dreal(jtemp(j))
     1        +dimag(ympole1(j))*dimag(jtemp(j)))*rsj**2
         ympole2(0) = ympole2(0)
     1        +(dreal(mbhmpole1(j))*dreal(jtemp(j))
     1        +dimag(mbhmpole1(j))*dimag(jtemp(j)))*rsj**2
         rsj=rsj*rscale1
      enddo

c     

      rsi=rscale1
      rsi7=rscale2
      rsi5=rscale1/rscale2
      do i = 1,nterms2
         ympole2(i) = ympole2(i) 
     1        + 2.0d0*(dreal(ympole1(0))*dreal(jtemp(i))
     2        +ima*dreal(ympole1(0))*dimag(jtemp(i)))*rsi5
         ympole2(i) = ympole2(i) 
     1        + 2.0d0*(dreal(mbhmpole1(0))*dreal(difftemp(i))
     2        +ima*dreal(mbhmpole1(0))*dimag(difftemp(i)))*rsi5
         mbhmpole2(i) = mbhmpole2(i) 
     1        + 2.0d0*(dreal(mbhmpole1(0))*dreal(powtemp(i))
     2        +ima*dreal(mbhmpole1(0))*dimag(powtemp(i)))*rsi5
         rsj=rscale1
         twojm1 = 1.0d0
         do j = 1,min(nterms1,i)
            ympole2(i) = ympole2(i)
     1           +(dreal(ympole1(j))*dreal(jtemp(i-j))
     1           -dimag(ympole1(j))*dimag(jtemp(i-j))
     2           +ima*(dreal(ympole1(j))*dimag(jtemp(i-j))
     3           +dimag(ympole1(j))*dreal(jtemp(i-j))))*rsi5
            ympole2(i) = ympole2(i)
     1           +(dreal(ympole1(j))*dreal(jtemp(i+j))
     1           +dimag(ympole1(j))*dimag(jtemp(i+j))
     2           +ima*(dreal(ympole1(j))*dimag(jtemp(i+j))
     3           -dimag(ympole1(j))*dreal(jtemp(i+j))))*rsj**2*rsi5
            mbhmpole2(i) = mbhmpole2(i)
     1           +(dreal(mbhmpole1(j))*dreal(powtemp(i-j))
     1           -dimag(mbhmpole1(j))*dimag(powtemp(i-j))
     2           +ima*(dreal(mbhmpole1(j))*dimag(powtemp(i-j))
     3           +dimag(mbhmpole1(j))*dreal(powtemp(i-j))))
     5           *rsi5
            ztemp = difftemp(i-j)
            ympole2(i) = ympole2(i)
     1           +(dreal(mbhmpole1(j))*dreal(ztemp)
     1           -dimag(mbhmpole1(j))*dimag(ztemp)
     2           +ima*(dreal(mbhmpole1(j))*dimag(ztemp)
     3           +dimag(mbhmpole1(j))*dreal(ztemp)))*rsi5
            ympole2(i) = ympole2(i)
     1           +(dreal(mbhmpole1(j))*dreal(jtemp(i+j))
     1           +dimag(mbhmpole1(j))*dimag(jtemp(i+j))
     2           +ima*(dreal(mbhmpole1(j))*dimag(jtemp(i+j))
     3           -dimag(mbhmpole1(j))*dreal(jtemp(i+j))))
     5           *rsj**2*rsi5
            rsj=rsj*rscale1
            twojm1 = twojm1*2.0d0
         enddo
         rsj=rscale1**(i+1)
         fs2=rsi5*rscale1**2
         do j = i+1,nterms1
            ympole2(i) = ympole2(i)
     1           +(dreal(mbhmpole1(j))*dreal(jtemp(j-i))
     1           +dimag(mbhmpole1(j))*dimag(jtemp(j-i))
     2           -ima*(dreal(mbhmpole1(j))*dimag(jtemp(j-i))
     3           -dimag(mbhmpole1(j))*dreal(jtemp(j-i))))*fs2
            ympole2(i) = ympole2(i)
     1           +(dreal(mbhmpole1(j))*dreal(jtemp(i+j))
     1           +dimag(mbhmpole1(j))*dimag(jtemp(i+j))
     2           +ima*(dreal(mbhmpole1(j))*dimag(jtemp(i+j))
     3           -dimag(mbhmpole1(j))*dreal(jtemp(i+j))))
     5           *rsj**2*rsi5
            ympole2(i) = ympole2(i)
     1           +(dreal(ympole1(j))*dreal(jtemp(j-i))
     1           +dimag(ympole1(j))*dimag(jtemp(j-i))
     2           -ima*(dreal(ympole1(j))*dimag(jtemp(j-i))
     3           -dimag(ympole1(j))*dreal(jtemp(j-i))))*fs2
            ympole2(i) = ympole2(i)
     1           +(dreal(ympole1(j))*dreal(jtemp(i+j))
     1           +dimag(ympole1(j))*dimag(jtemp(i+j))
     2           +ima*(dreal(ympole1(j))*dimag(jtemp(i+j))
     3           -dimag(ympole1(j))*dreal(jtemp(i+j))))
     5           *rsj**2*rsi5
            rsj=rsj*rscale1
            fs2=fs2*rscale1**2
         enddo
         rsi=rsi*rscale1
         rsi7=rsi7*rscale2
         rsi5=rsi5*rscale1/rscale2
      enddo

      return
      end

      subroutine mbh2dmpmp_vec(beta,rscale1,center1,mbhmpole1,ympole1,
     1     nterms1,mstep1,rscale2,center2,mbhmpole2,ympole2,nterms2,
     2     mstep2,nexp)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     rscale1           : scaling for original expansion
c     center1           : center for original expansion
c     coeffs1_1(nterms1): original coeffs, difference-type functions
c     coeffs2_1(nterms1): original coeffs, modified Bessel functions
c     nterms1           : number of terms in original expansion
c     rscale2           : scaling for new expansion
c     center2           : center for new expansion
c     rad2              : radius of circle to sample on 
c                         for new expansion
c     nterms2           : number of terms in original expansion
c     zcirc2(2,nterms2) : precomp'd evenly-spaced points on unit circle
c     beta              : the modified biharmonic parameter

c     work(*)           : work array. recommended length 16*nterms+500
c
c     OUTPUT:
c
c     coeffs1_2(nterms2): new coeffs, difference-type functions
c     coeffs2_2(nterms2): new coeffs, modified Bessel functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 rscale1, center1(2), beta
      real *8 rscale2, center2(2)
      complex *16 mbhmpole1(0:mstep1,nexp), ympole1(0:mstep1,nexp)
      complex *16 mbhmpole2(0:mstep2,nexp), ympole2(0:mstep2,nexp)
      integer nterms1, nterms2, nexp, mstep1, mstep2
c     local
      integer i, j, ifder, nterms, m, l, iii
      real *8 zdiff(2), r, theta, pi
      real *8 diffs(0:200), ders(0:200), ival(0:200), dfac(0:200)
      real *8 pow(0:200), dpow(0:200), dfac2(0:200)
      real *8 carray(0:200,0:200), twojm1
      real *8 rsj, rsi5, rsi7, rsi, fs2, rtemp
      complex *16 jtemp(0:200), difftemp(0:200), powtemp(0:200)
      complex *16 zk, z, ima, ztemp1, zmul, ztemp
      data ima /(0.0d0,1.0d0)/

      pi = 4.0d0*datan(1.0d0)
      zk = ima*beta

      nterms = nterms1+nterms2

c     

      dfac(0) = 1.0d0
      dfac2(0) = 1.0d0
      do i = 1,nterms
         dfac(i) = dfac(i-1)*i
         dfac2(i) = dfac2(i-1)*2.0d0/i
      enddo

      do l = 0,nterms
         carray(l,0) = 1.0d0
      enddo
      do m=1,nterms
         carray(m,m) = 1.0d0
         do l=m+1,nterms
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
         enddo
      enddo
c     

c     get local difference and BesselI vals
      
      zdiff(1)=center2(1)-center1(1)
      zdiff(2)=center2(2)-center1(2)
      call h2cart2polar(zdiff,r,theta)
      theta=theta-pi
      z=zk*r
      ifder=0

      call diffszkik_fast(r,beta,rscale1,diffs,ifder,ders,ival,
     1     nterms)

      rtemp = r
      call mbh2d_rksc(pow,dpow,rtemp,beta,rscale1,nterms)

c     
      jtemp(0) = ival(0)
      difftemp(0) = diffs(0)
      powtemp(0) = pow(0)
      zmul=exp(ima*theta)
      ztemp1= zmul

      do j = 1,nterms
         jtemp( j) = ztemp1*ival(j)
         difftemp(j) = ztemp1*diffs(j)
         powtemp(j) = ztemp1*pow(j)
         ztemp1= ztemp1*zmul
      enddo

C$OMP PARALLEL DO PRIVATE(iii,rsj,rsi,rsi7,rsi5,fs2,twojm1)
C$OMP& IF(nexp .gt. 10)
C$OMP& SCHEDULE(static)
      do iii = 1,nexp

         do i = 0,nterms2
            mbhmpole2(i,iii) = 0
            ympole2(i,iii) = 0
         enddo

c     shift the expansion

         rsj = rscale1

         ympole2(0,iii) = ympole2(0,iii) 
     1        + dreal(ympole1(0,iii))*jtemp(0)
         ympole2(0,iii) = ympole2(0,iii) 
     1        + dreal(mbhmpole1(0,iii))*difftemp(0)
         mbhmpole2(0,iii) = mbhmpole2(0,iii) 
     1        + dreal(mbhmpole1(0,iii))

         do j = 1,nterms1
            ympole2(0,iii) = ympole2(0,iii)
     1           +(dreal(ympole1(j,iii))*dreal(jtemp(j))
     1           +dimag(ympole1(j,iii))*dimag(jtemp(j)))*rsj**2
            ympole2(0,iii) = ympole2(0,iii)
     1           +(dreal(mbhmpole1(j,iii))*dreal(jtemp(j))
     1           +dimag(mbhmpole1(j,iii))*dimag(jtemp(j)))*rsj**2
            rsj=rsj*rscale1
         enddo

c     
         rsi=rscale1
         rsi7=rscale2
         rsi5=rscale1/rscale2
         do i = 1,nterms2
            ympole2(i,iii) = ympole2(i,iii) 
     1           + 2.0d0*(dreal(ympole1(0,iii))*dreal(jtemp(i))
     2           +ima*dreal(ympole1(0,iii))*dimag(jtemp(i)))*rsi5
            ympole2(i,iii) = ympole2(i,iii) 
     1           + 2.0d0*(dreal(mbhmpole1(0,iii))*dreal(difftemp(i))
     2           +ima*dreal(mbhmpole1(0,iii))*dimag(difftemp(i)))*rsi5
            mbhmpole2(i,iii) = mbhmpole2(i,iii) 
     1           + 2.0d0*(dreal(mbhmpole1(0,iii))*dreal(powtemp(i))
     2           +ima*dreal(mbhmpole1(0,iii))*dimag(powtemp(i)))*rsi5
            rsj=rscale1
            twojm1 = 1.0d0
            do j = 1,min(nterms1,i)
               ympole2(i,iii) = ympole2(i,iii)
     1              +(dreal(ympole1(j,iii))*dreal(jtemp(i-j))
     1              -dimag(ympole1(j,iii))*dimag(jtemp(i-j))
     2              +ima*(dreal(ympole1(j,iii))*dimag(jtemp(i-j))
     3              +dimag(ympole1(j,iii))*dreal(jtemp(i-j))))*rsi5
               ympole2(i,iii) = ympole2(i,iii)
     1              +(dreal(ympole1(j,iii))*dreal(jtemp(i+j))
     1              +dimag(ympole1(j,iii))*dimag(jtemp(i+j))
     2              +ima*(dreal(ympole1(j,iii))*dimag(jtemp(i+j))
     3              -dimag(ympole1(j,iii))*dreal(jtemp(i+j))))
     5              *rsj**2*rsi5
               mbhmpole2(i,iii) = mbhmpole2(i,iii)
     1              +(dreal(mbhmpole1(j,iii))*dreal(powtemp(i-j))
     1              -dimag(mbhmpole1(j,iii))*dimag(powtemp(i-j))
     2              +ima*(dreal(mbhmpole1(j,iii))*dimag(powtemp(i-j))
     3              +dimag(mbhmpole1(j,iii))*dreal(powtemp(i-j))))
     5              *rsi5
               ztemp = difftemp(i-j)
               ympole2(i,iii) = ympole2(i,iii)
     1              +(dreal(mbhmpole1(j,iii))*dreal(ztemp)
     1              -dimag(mbhmpole1(j,iii))*dimag(ztemp)
     2              +ima*(dreal(mbhmpole1(j,iii))*dimag(ztemp)
     3              +dimag(mbhmpole1(j,iii))*dreal(ztemp)))*rsi5
               ympole2(i,iii) = ympole2(i,iii)
     1              +(dreal(mbhmpole1(j,iii))*dreal(jtemp(i+j))
     1              +dimag(mbhmpole1(j,iii))*dimag(jtemp(i+j))
     2              +ima*(dreal(mbhmpole1(j,iii))*dimag(jtemp(i+j))
     3              -dimag(mbhmpole1(j,iii))*dreal(jtemp(i+j))))
     5              *rsj**2*rsi5
               rsj=rsj*rscale1
               twojm1 = twojm1*2.0d0
            enddo
            rsj=rscale1**(i+1)
            fs2=rsi5*rscale1**2
            do j = i+1,nterms1
               ympole2(i,iii) = ympole2(i,iii)
     1              +(dreal(mbhmpole1(j,iii))*dreal(jtemp(j-i))
     1              +dimag(mbhmpole1(j,iii))*dimag(jtemp(j-i))
     2              -ima*(dreal(mbhmpole1(j,iii))*dimag(jtemp(j-i))
     3              -dimag(mbhmpole1(j,iii))*dreal(jtemp(j-i))))*fs2
               ympole2(i,iii) = ympole2(i,iii)
     1              +(dreal(mbhmpole1(j,iii))*dreal(jtemp(i+j))
     1              +dimag(mbhmpole1(j,iii))*dimag(jtemp(i+j))
     2              +ima*(dreal(mbhmpole1(j,iii))*dimag(jtemp(i+j))
     3              -dimag(mbhmpole1(j,iii))*dreal(jtemp(i+j))))
     5              *rsj**2*rsi5
               ympole2(i,iii) = ympole2(i,iii)
     1              +(dreal(ympole1(j,iii))*dreal(jtemp(j-i))
     1              +dimag(ympole1(j,iii))*dimag(jtemp(j-i))
     2              -ima*(dreal(ympole1(j,iii))*dimag(jtemp(j-i))
     3              -dimag(ympole1(j,iii))*dreal(jtemp(j-i))))*fs2
               ympole2(i,iii) = ympole2(i,iii)
     1              +(dreal(ympole1(j,iii))*dreal(jtemp(i+j))
     1              +dimag(ympole1(j,iii))*dimag(jtemp(i+j))
     2              +ima*(dreal(ympole1(j,iii))*dimag(jtemp(i+j))
     3              -dimag(ympole1(j,iii))*dreal(jtemp(i+j))))
     5              *rsj**2*rsi5
               rsj=rsj*rscale1
               fs2=fs2*rscale1**2
            enddo
            rsi=rsi*rscale1
            rsi7=rsi7*rscale2
            rsi5=rsi5*rscale1/rscale2
         enddo
      enddo
C$OMP END PARALLEL DO

      return
      end

      subroutine mbh2dmploc(beta,rscale1,center1,mbhmpole,ympole,
     1     nterms1,rscale2,center2,mbhloc,lloc,nterms2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT: 
c
c     rscale1           : scaling for original expansion
c     center1           : center for original expansion
c     coeffs1_1(nterms1): original coeffs, difference-type functions
c     coeffs2_1(nterms1): original coeffs, modified Bessel functions
c     nterms1           : number of terms in original expansion
c     rscale2           : scaling for new expansion
c     center2           : center for new expansion
c     rad2              : radius of circle to sample on 
c                         for new expansion
c     nterms2           : number of terms in original expansion
c     zcirc2(2,nterms2) : precomp'd evenly-spaced points on unit circle
c     beta              : the modified biharmonic parameter

c     work(*)           : work array. recommended length 16*nterms+500
c
c     OUTPUT:
c
c     coeffs1_2(nterms2): new coeffs, difference-type functions
c     coeffs2_2(nterms2): new coeffs, modified Bessel functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      complex *16 mbhmpole(0:nterms1), ympole(0:nterms1)
      complex *16 mbhloc(0:nterms2), lloc(0:nterms2)
      real *8 rscale1, rscale2, center1(2), center2(2), beta
      integer nterms1, nterms2
c     local
      real *8 zdiff(2), r, theta, pi, rsi, rsj, rsi5, rsi7, rsj2
      real *8 rsi52, rsi2
      real *8 diffs(0:200), ders(0:200), kvec(0:200)
      real *8 pow(0:200), dpow(0:200), dfac(0:200), dfac2(0:200)
      complex *16 ktemp(0:200), difftemp(0:200), powtemp(0:200)
      complex *16 ztemp1, zmul, ima
      integer nterms, i, j, isign, ifders
      data ima /(0.0d0,1.0d0)/
      
      pi = 4.0d0*datan(1.0d0)

      nterms = nterms1+nterms2

      do i = 0,nterms2
         mbhloc(i) = 0.0d0
         lloc(i) = 0.0d0
      enddo

      dfac2(0) = 1.0d0
      do i = 1,nterms
         dfac2(i) = dfac2(i-1)/(2.0d0*i)
      enddo

      zdiff(1) = center2(1)-center1(1)
      zdiff(2) = center2(2)-center1(2)
      call h2cart2polar(zdiff,r,theta)

      theta = theta-pi

c     get values of (beta*r)^-k
      call mbh2d_rmk(pow,dpow,r,beta,rscale1,nterms)

c     get values of difference functions
      ifders = 0
      call diffslogbk_fast(r,beta,rscale1,diffs,ifders,ders,kvec,
     1     nterms)

c     form terms for shifting 

      ktemp(0) = kvec(0)
      difftemp(0) = diffs(0)
      powtemp(0) = pow(0)
      zmul=-exp(ima*theta)
      ztemp1= zmul
      do j = 1,nterms
         ktemp( j) = ztemp1*kvec(j)
         difftemp(j) = ztemp1*diffs(j)
         powtemp(j) = ztemp1*pow(j)
         ztemp1= ztemp1*zmul
      enddo

c     shift

      mbhloc(0) = mbhloc(0) + mbhmpole(0)*ktemp(0)
      mbhloc(0) = mbhloc(0) + ympole(0)*ktemp(0)
c      lloc(0) = lloc(0) + mbhmpole(0)*difftemp(0)
      lloc(0) = lloc(0) + mbhmpole(0)*(ktemp(0)+pow(0))
      lloc(0) = lloc(0) + ympole(0)*ktemp(0)
      do j = 1,nterms1
         mbhloc(0) = mbhloc(0)+(dreal(mbhmpole(j))*dreal(ktemp(j))
     1        +dimag(mbhmpole(j))*dimag(ktemp(j)))
         mbhloc(0) = mbhloc(0)+(dreal(ympole(j))*dreal(ktemp(j))
     1        +dimag(ympole(j))*dimag(ktemp(j)))
         lloc(0) = lloc(0)+(dreal(mbhmpole(j))*dreal(difftemp(j))
     1        +dimag(mbhmpole(j))*dimag(difftemp(j)))
         lloc(0) = lloc(0)+(dreal(ympole(j))*dreal(ktemp(j))
     1        +dimag(ympole(j))*dimag(ktemp(j)))
      enddo
c     
      rsi=rscale1
      rsi2=rscale1**2
      rsi5=rscale2/rscale1
      rsi52=rsi5*rscale2/rscale1
      isign = -1
      do i = 1,nterms2
         mbhloc(i) = mbhloc(i) 
     1        + 2.0d0*(dreal(mbhmpole(0))*dreal(ktemp(i))
     2        + ima*dreal(mbhmpole(0))*dimag(ktemp(i)))
         mbhloc(i) = mbhloc(i) 
     1        + 2.0d0*(dreal(ympole(0))*dreal(ktemp(i))
     2        +ima*dreal(ympole(0))*dimag(ktemp(i)))
         lloc(i) = lloc(i) 
     1        + 2.0d0*(dreal(mbhmpole(0))*dreal(difftemp(i))
     2        + ima*dreal(mbhmpole(0))*dimag(difftemp(i)))*dfac2(i)
         lloc(i) = lloc(i) 
     1        + 2.0d0*(dreal(ympole(0))*dreal(ktemp(i))
     2        +ima*dreal(ympole(0))*dimag(ktemp(i)))*dfac2(i)
c         lloc(i) = lloc(i)
c     1        - (dreal(mbhmpole(0))*dreal(powtemp(i))
c     2        +ima*dreal(mbhmpole(0))*dimag(powtemp(i)))/i

         rsj=rscale1
         rsj2=rscale1**2
         do j = 1,min(nterms1,i)
            mbhloc(i) = mbhloc(i)+(dreal(mbhmpole(j))*dreal(ktemp(i-j))
     1           -dimag(mbhmpole(j))*dimag(ktemp(i-j))
     2           +ima*(dreal(mbhmpole(j))*dimag(ktemp(i-j))
     3           +dimag(mbhmpole(j))*dreal(ktemp(i-j))))*rsj2
            mbhloc(i) = mbhloc(i)+(dreal(mbhmpole(j))*dreal(ktemp(i+j))
     1           +dimag(mbhmpole(j))*dimag(ktemp(i+j))
     2           +ima*(dreal(mbhmpole(j))*dimag(ktemp(i+j))
     3           -dimag(mbhmpole(j))*dreal(ktemp(i+j))))
            mbhloc(i) = mbhloc(i)+(dreal(ympole(j))*dreal(ktemp(i-j))
     1           -dimag(ympole(j))*dimag(ktemp(i-j))
     2           +ima*(dreal(ympole(j))*dimag(ktemp(i-j))
     3           +dimag(ympole(j))*dreal(ktemp(i-j))))*rsj2
            mbhloc(i) = mbhloc(i)+(dreal(ympole(j))*dreal(ktemp(i+j))
     1           +dimag(ympole(j))*dimag(ktemp(i+j))
     2           +ima*(dreal(ympole(j))*dimag(ktemp(i+j))
     3           -dimag(ympole(j))*dreal(ktemp(i+j))))
            lloc(i) = lloc(i)+(dreal(mbhmpole(j))*dreal(ktemp(i-j))
     1           -dimag(mbhmpole(j))*dimag(ktemp(i-j))
     2           +ima*(dreal(mbhmpole(j))*dimag(ktemp(i-j))
     3           +dimag(mbhmpole(j))*dreal(ktemp(i-j))))*rsj2*dfac2(i)
            lloc(i) = lloc(i)+(dreal(mbhmpole(j))*dreal(difftemp(i+j))
     1           +dimag(mbhmpole(j))*dimag(difftemp(i+j))
     2           +ima*(dreal(mbhmpole(j))*dimag(difftemp(i+j))
     3           -dimag(mbhmpole(j))*dreal(difftemp(i+j))))*dfac2(i)
            lloc(i) = lloc(i)+(dreal(ympole(j))*dreal(ktemp(i-j))
     1           -dimag(ympole(j))*dimag(ktemp(i-j))
     2           +ima*(dreal(ympole(j))*dimag(ktemp(i-j))
     3           +dimag(ympole(j))*dreal(ktemp(i-j))))*rsj2*dfac2(i)
            lloc(i) = lloc(i)+(dreal(ympole(j))*dreal(ktemp(i+j))
     1           +dimag(ympole(j))*dimag(ktemp(i+j))
     2           +ima*(dreal(ympole(j))*dimag(ktemp(i+j))
     3           -dimag(ympole(j))*dreal(ktemp(i+j))))*dfac2(i)

            rsj=rsj*rscale1
            rsj2=rsj2*rscale1**2
         enddo
         do j = i+1,nterms1
            mbhloc(i) = mbhloc(i)+(dreal(mbhmpole(j))*dreal(ktemp(j-i))
     1           +dimag(mbhmpole(j))*dimag(ktemp(j-i))
     2           -ima*(dreal(mbhmpole(j))*dimag(ktemp(j-i))
     3           -dimag(mbhmpole(j))*dreal(ktemp(j-i))))*rsi2
            mbhloc(i) = mbhloc(i)+(dreal(mbhmpole(j))*dreal(ktemp(i+j))
     1           +dimag(mbhmpole(j))*dimag(ktemp(i+j))
     2           +ima*(dreal(mbhmpole(j))*dimag(ktemp(i+j))
     3           -dimag(mbhmpole(j))*dreal(ktemp(i+j))))
            mbhloc(i) = mbhloc(i)+(dreal(ympole(j))*dreal(ktemp(j-i))
     1           +dimag(ympole(j))*dimag(ktemp(j-i))
     2           -ima*(dreal(ympole(j))*dimag(ktemp(j-i))
     3           -dimag(ympole(j))*dreal(ktemp(j-i))))*rsi2
            mbhloc(i) = mbhloc(i)+(dreal(ympole(j))*dreal(ktemp(i+j))
     1           +dimag(ympole(j))*dimag(ktemp(i+j))
     2           +ima*(dreal(ympole(j))*dimag(ktemp(i+j))
     3           -dimag(ympole(j))*dreal(ktemp(i+j))))
            lloc(i) = lloc(i)+(dreal(mbhmpole(j))*dreal(ktemp(j-i))
     1           +dimag(mbhmpole(j))*dimag(ktemp(j-i))
     2           -ima*(dreal(mbhmpole(j))*dimag(ktemp(j-i))
     3           -dimag(mbhmpole(j))*dreal(ktemp(j-i))))*rsi2*dfac2(i)
            lloc(i) = lloc(i)+(dreal(mbhmpole(j))*dreal(difftemp(i+j))
     1           +dimag(mbhmpole(j))*dimag(difftemp(i+j))
     2           +ima*(dreal(mbhmpole(j))*dimag(difftemp(i+j))
     3           -dimag(mbhmpole(j))*dreal(difftemp(i+j))))*dfac2(i)
            lloc(i) = lloc(i)+(dreal(ympole(j))*dreal(ktemp(j-i))
     1           +dimag(ympole(j))*dimag(ktemp(j-i))
     2           -ima*(dreal(ympole(j))*dimag(ktemp(j-i))
     3           -dimag(ympole(j))*dreal(ktemp(j-i))))*rsi2*dfac2(i)
            lloc(i) = lloc(i)+(dreal(ympole(j))*dreal(ktemp(i+j))
     1           +dimag(ympole(j))*dimag(ktemp(i+j))
     2           +ima*(dreal(ympole(j))*dimag(ktemp(i+j))
     3           -dimag(ympole(j))*dreal(ktemp(i+j))))*dfac2(i)

         enddo
         mbhloc(i)=mbhloc(i)*isign
         lloc(i)=lloc(i)*isign
         isign = -isign
         rsi=rsi*rscale1
         rsi2=rsi2*rscale1**2
      enddo
      rsi=rscale2/rscale1
      do i = 1,nterms2
         mbhloc(+i) = mbhloc(+i)*rsi
         lloc(+i) = lloc(+i)*rsi
         rsi=rsi*rscale2/rscale1
      enddo


      return
      end

      subroutine mbh2dmploc_pre(beta,rscale1,center1,nterms1,
     1     rscale2,center2,nterms2,dfac2,kvec,diffs,pow)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT: 
c
c     rscale1           : scaling for original expansion
c     center1           : center for original expansion
c     coeffs1_1(nterms1): original coeffs, difference-type functions
c     coeffs2_1(nterms1): original coeffs, modified Bessel functions
c     nterms1           : number of terms in original expansion
c     rscale2           : scaling for new expansion
c     center2           : center for new expansion
c     rad2              : radius of circle to sample on 
c                         for new expansion
c     nterms2           : number of terms in original expansion
c     zcirc2(2,nterms2) : precomp'd evenly-spaced points on unit circle
c     beta              : the modified biharmonic parameter

c     work(*)           : work array. recommended length 16*nterms+500
c
c     OUTPUT:
c
c     coeffs1_2(nterms2): new coeffs, difference-type functions
c     coeffs2_2(nterms2): new coeffs, modified Bessel functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 rscale1, rscale2, center1(2), center2(2), beta
      real *8 dfac2(0:1), kvec(0:1), diffs(0:1), pow(0:1)
      integer nterms1, nterms2
c     local
      real *8 zdiff(2), r, theta, pi, rsi, rsj, rsi5, rsi7, rsj2
      real *8 rsi52, rsi2
      real *8 ders(0:200)
      real *8 dpow(0:200), dfac(0:200)
      complex *16 ztemp1, zmul, ima
      integer nterms, i, j, isign, ifders
      data ima /(0.0d0,1.0d0)/
      
      pi = 4.0d0*datan(1.0d0)

      nterms = nterms1+nterms2

      dfac2(0) = 1.0d0
      do i = 1,nterms
         dfac2(i) = dfac2(i-1)/(2.0d0*i)
      enddo

      zdiff(1) = center2(1)-center1(1)
      zdiff(2) = center2(2)-center1(2)
      call h2cart2polar(zdiff,r,theta)

      theta = theta-pi

c     get values of (beta*r)^-k
      call mbh2d_rmk(pow,dpow,r,beta,rscale1,nterms)

c     get values of difference functions
      ifders = 0
      call diffslogbk_fast(r,beta,rscale1,diffs,ifders,ders,kvec,
     1     nterms)

      return
      end

      subroutine mbh2dmploc_wpre(beta,rscale1,center1,mbhmpole,ympole,
     1     nterms1,rscale2,center2,mbhloc,lloc,nterms2,dfac2,kvec,
     2     diffs,pow)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT: 
c
c     rscale1           : scaling for original expansion
c     center1           : center for original expansion
c     coeffs1_1(nterms1): original coeffs, difference-type functions
c     coeffs2_1(nterms1): original coeffs, modified Bessel functions
c     nterms1           : number of terms in original expansion
c     rscale2           : scaling for new expansion
c     center2           : center for new expansion
c     rad2              : radius of circle to sample on 
c                         for new expansion
c     nterms2           : number of terms in original expansion
c     zcirc2(2,nterms2) : precomp'd evenly-spaced points on unit circle
c     beta              : the modified biharmonic parameter

c     work(*)           : work array. recommended length 16*nterms+500
c
c     OUTPUT:
c
c     coeffs1_2(nterms2): new coeffs, difference-type functions
c     coeffs2_2(nterms2): new coeffs, modified Bessel functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      complex *16 mbhmpole(0:nterms1), ympole(0:nterms1)
      complex *16 mbhloc(0:nterms2), lloc(0:nterms2)
      real *8 rscale1, rscale2, center1(2), center2(2), beta
      real *8 dfac2(0:1), kvec(0:1), diffs(0:1), pow(0:1)
      integer nterms1, nterms2
c     local
      real *8 zdiff(2), r, theta, pi, rsi, rsj, rsi5, rsi7, rsj2
      real *8 rsi52, rsi2
      complex *16 ktemp(0:200), difftemp(0:200), powtemp(0:200)
      complex *16 ztemp1, zmul, ima
      integer nterms, i, j, isign, ifders
      data ima /(0.0d0,1.0d0)/
      
      pi = 4.0d0*datan(1.0d0)

      nterms = nterms1+nterms2

      do i = 0,nterms2
         mbhloc(i) = 0.0d0
         lloc(i) = 0.0d0
      enddo

      zdiff(1) = center2(1)-center1(1)
      zdiff(2) = center2(2)-center1(2)
      call h2cart2polar(zdiff,r,theta)

      theta = theta-pi

c     form terms for shifting 

      ktemp(0) = kvec(0)
      difftemp(0) = diffs(0)
      powtemp(0) = pow(0)
      zmul=-exp(ima*theta)
      ztemp1= zmul
      do j = 1,nterms
         ktemp( j) = ztemp1*kvec(j)
         difftemp(j) = ztemp1*diffs(j)
         powtemp(j) = ztemp1*pow(j)
         ztemp1= ztemp1*zmul
      enddo

c     shift

      mbhloc(0) = mbhloc(0) + mbhmpole(0)*ktemp(0)
      mbhloc(0) = mbhloc(0) + ympole(0)*ktemp(0)
      lloc(0) = lloc(0) + mbhmpole(0)*difftemp(0)
c      lloc(0) = lloc(0) + mbhmpole(0)*(ktemp(0)+pow(0))
      lloc(0) = lloc(0) + ympole(0)*ktemp(0)

      do j = 1,nterms1
         mbhloc(0) = mbhloc(0)+(dreal(mbhmpole(j))*dreal(ktemp(j))
     1        +dimag(mbhmpole(j))*dimag(ktemp(j)))
         mbhloc(0) = mbhloc(0)+(dreal(ympole(j))*dreal(ktemp(j))
     1        +dimag(ympole(j))*dimag(ktemp(j)))
         lloc(0) = lloc(0)+(dreal(mbhmpole(j))*dreal(difftemp(j))
     1        +dimag(mbhmpole(j))*dimag(difftemp(j)))
         lloc(0) = lloc(0)+(dreal(ympole(j))*dreal(ktemp(j))
     1        +dimag(ympole(j))*dimag(ktemp(j)))
      enddo
c    
      rsi=rscale1
      rsi2=rscale1**2
      rsi5=rscale2/rscale1
      rsi52=rsi5*rscale2/rscale1
      isign = -1
      do i = 1,nterms2
         mbhloc(i) = mbhloc(i) 
     1        + 2.0d0*(dreal(mbhmpole(0))*dreal(ktemp(i))
     2        + ima*dreal(mbhmpole(0))*dimag(ktemp(i)))
         mbhloc(i) = mbhloc(i) 
     1        + 2.0d0*(dreal(ympole(0))*dreal(ktemp(i))
     2        +ima*dreal(ympole(0))*dimag(ktemp(i)))
         lloc(i) = lloc(i) 
     1        + 2.0d0*(dreal(mbhmpole(0))*dreal(difftemp(i))
     2        + ima*dreal(mbhmpole(0))*dimag(difftemp(i)))*dfac2(i)
         lloc(i) = lloc(i) 
     1        + 2.0d0*(dreal(ympole(0))*dreal(ktemp(i))
     2        +ima*dreal(ympole(0))*dimag(ktemp(i)))*dfac2(i)
c         lloc(i) = lloc(i)
c     1        - (dreal(mbhmpole(0))*dreal(powtemp(i))
c     2        +ima*dreal(mbhmpole(0))*dimag(powtemp(i)))/i
         rsj=rscale1
         rsj2=rscale1**2
         do j = 1,min(nterms1,i)
            mbhloc(i) = mbhloc(i)+(dreal(mbhmpole(j))*dreal(ktemp(i-j))
     1           -dimag(mbhmpole(j))*dimag(ktemp(i-j))
     2           +ima*(dreal(mbhmpole(j))*dimag(ktemp(i-j))
     3           +dimag(mbhmpole(j))*dreal(ktemp(i-j))))*rsj2
            mbhloc(i) = mbhloc(i)+(dreal(mbhmpole(j))*dreal(ktemp(i+j))
     1           +dimag(mbhmpole(j))*dimag(ktemp(i+j))
     2           +ima*(dreal(mbhmpole(j))*dimag(ktemp(i+j))
     3           -dimag(mbhmpole(j))*dreal(ktemp(i+j))))
            mbhloc(i) = mbhloc(i)+(dreal(ympole(j))*dreal(ktemp(i-j))
     1           -dimag(ympole(j))*dimag(ktemp(i-j))
     2           +ima*(dreal(ympole(j))*dimag(ktemp(i-j))
     3           +dimag(ympole(j))*dreal(ktemp(i-j))))*rsj2
            mbhloc(i) = mbhloc(i)+(dreal(ympole(j))*dreal(ktemp(i+j))
     1           +dimag(ympole(j))*dimag(ktemp(i+j))
     2           +ima*(dreal(ympole(j))*dimag(ktemp(i+j))
     3           -dimag(ympole(j))*dreal(ktemp(i+j))))
            lloc(i) = lloc(i)+(dreal(mbhmpole(j))*dreal(ktemp(i-j))
     1           -dimag(mbhmpole(j))*dimag(ktemp(i-j))
     2           +ima*(dreal(mbhmpole(j))*dimag(ktemp(i-j))
     3           +dimag(mbhmpole(j))*dreal(ktemp(i-j))))*rsj2*dfac2(i)
            lloc(i) = lloc(i)+(dreal(mbhmpole(j))*dreal(difftemp(i+j))
     1           +dimag(mbhmpole(j))*dimag(difftemp(i+j))
     2           +ima*(dreal(mbhmpole(j))*dimag(difftemp(i+j))
     3           -dimag(mbhmpole(j))*dreal(difftemp(i+j))))*dfac2(i)
            lloc(i) = lloc(i)+(dreal(ympole(j))*dreal(ktemp(i-j))
     1           -dimag(ympole(j))*dimag(ktemp(i-j))
     2           +ima*(dreal(ympole(j))*dimag(ktemp(i-j))
     3           +dimag(ympole(j))*dreal(ktemp(i-j))))*rsj2*dfac2(i)
            lloc(i) = lloc(i)+(dreal(ympole(j))*dreal(ktemp(i+j))
     1           +dimag(ympole(j))*dimag(ktemp(i+j))
     2           +ima*(dreal(ympole(j))*dimag(ktemp(i+j))
     3           -dimag(ympole(j))*dreal(ktemp(i+j))))*dfac2(i)

            rsj=rsj*rscale1
            rsj2=rsj2*rscale1**2
         enddo
         do j = i+1,nterms1
            mbhloc(i) = mbhloc(i)+(dreal(mbhmpole(j))*dreal(ktemp(j-i))
     1           +dimag(mbhmpole(j))*dimag(ktemp(j-i))
     2           -ima*(dreal(mbhmpole(j))*dimag(ktemp(j-i))
     3           -dimag(mbhmpole(j))*dreal(ktemp(j-i))))*rsi2
            mbhloc(i) = mbhloc(i)+(dreal(mbhmpole(j))*dreal(ktemp(i+j))
     1           +dimag(mbhmpole(j))*dimag(ktemp(i+j))
     2           +ima*(dreal(mbhmpole(j))*dimag(ktemp(i+j))
     3           -dimag(mbhmpole(j))*dreal(ktemp(i+j))))
            mbhloc(i) = mbhloc(i)+(dreal(ympole(j))*dreal(ktemp(j-i))
     1           +dimag(ympole(j))*dimag(ktemp(j-i))
     2           -ima*(dreal(ympole(j))*dimag(ktemp(j-i))
     3           -dimag(ympole(j))*dreal(ktemp(j-i))))*rsi2
            mbhloc(i) = mbhloc(i)+(dreal(ympole(j))*dreal(ktemp(i+j))
     1           +dimag(ympole(j))*dimag(ktemp(i+j))
     2           +ima*(dreal(ympole(j))*dimag(ktemp(i+j))
     3           -dimag(ympole(j))*dreal(ktemp(i+j))))
            lloc(i) = lloc(i)+(dreal(mbhmpole(j))*dreal(ktemp(j-i))
     1           +dimag(mbhmpole(j))*dimag(ktemp(j-i))
     2           -ima*(dreal(mbhmpole(j))*dimag(ktemp(j-i))
     3           -dimag(mbhmpole(j))*dreal(ktemp(j-i))))*rsi2*dfac2(i)
            lloc(i) = lloc(i)+(dreal(mbhmpole(j))*dreal(difftemp(i+j))
     1           +dimag(mbhmpole(j))*dimag(difftemp(i+j))
     2           +ima*(dreal(mbhmpole(j))*dimag(difftemp(i+j))
     3           -dimag(mbhmpole(j))*dreal(difftemp(i+j))))*dfac2(i)
            lloc(i) = lloc(i)+(dreal(ympole(j))*dreal(ktemp(j-i))
     1           +dimag(ympole(j))*dimag(ktemp(j-i))
     2           -ima*(dreal(ympole(j))*dimag(ktemp(j-i))
     3           -dimag(ympole(j))*dreal(ktemp(j-i))))*rsi2*dfac2(i)
            lloc(i) = lloc(i)+(dreal(ympole(j))*dreal(ktemp(i+j))
     1           +dimag(ympole(j))*dimag(ktemp(i+j))
     2           +ima*(dreal(ympole(j))*dimag(ktemp(i+j))
     3           -dimag(ympole(j))*dreal(ktemp(i+j))))*dfac2(i)

         enddo
         mbhloc(i)=mbhloc(i)*isign
         lloc(i)=lloc(i)*isign
         isign = -isign
         rsi=rsi*rscale1
         rsi2=rsi2*rscale1**2
      enddo
      rsi=rscale2/rscale1
      do i = 1,nterms2
         mbhloc(+i) = mbhloc(+i)*rsi
         lloc(+i) = lloc(+i)*rsi
         rsi=rsi*rscale2/rscale1
      enddo

      return
      end

      subroutine mbh2dmploc_vec(beta,rscale1,center1,mbhmpole,ympole,
     1     nterms1,rscale2,center2,mbhloc,lloc,nterms2,nexp)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT: 
c
c     rscale1           : scaling for original expansion
c     center1           : center for original expansion
c     coeffs1_1(nterms1): original coeffs, difference-type functions
c     coeffs2_1(nterms1): original coeffs, modified Bessel functions
c     nterms1           : number of terms in original expansion
c     rscale2           : scaling for new expansion
c     center2           : center for new expansion
c     rad2              : radius of circle to sample on 
c                         for new expansion
c     nterms2           : number of terms in original expansion
c     zcirc2(2,nterms2) : precomp'd evenly-spaced points on unit circle
c     beta              : the modified biharmonic parameter

c     work(*)           : work array. recommended length 16*nterms+500
c
c     OUTPUT:
c
c     coeffs1_2(nterms2): new coeffs, difference-type functions
c     coeffs2_2(nterms2): new coeffs, modified Bessel functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      complex *16 mbhmpole(0:nterms1,nexp), ympole(0:nterms1,nexp)
      complex *16 mbhloc(0:nterms2,nexp), lloc(0:nterms2,nexp)
      real *8 rscale1, rscale2, center1(2), center2(2), beta
      integer nterms1, nterms2, nexp
c     local
      real *8 zdiff(2), r, theta, pi, rsi, rsj, rsi5, rsi7, rsj2
      real *8 rsi52, rsi2
      real *8 diffs(0:200), ders(0:200), kvec(0:200)
      real *8 pow(0:200), dpow(0:200), dfac(0:200), dfac2(0:200)
      complex *16 ktemp(0:200), difftemp(0:200), powtemp(0:200)
      complex *16 ztemp1, zmul, ima
      integer nterms, i, j, isign, ifders, iii
      data ima /(0.0d0,1.0d0)/
      
      pi = 4.0d0*datan(1.0d0)

      nterms = nterms1+nterms2

      dfac2(0) = 1.0d0
      do i = 1,nterms
         dfac2(i) = dfac2(i-1)/(2.0d0*i)
      enddo

      zdiff(1) = center2(1)-center1(1)
      zdiff(2) = center2(2)-center1(2)
      call h2cart2polar(zdiff,r,theta)

      theta = theta-pi

c     get values of (beta*r)^-k
      call mbh2d_rmk(pow,dpow,r,beta,rscale1,nterms)

c     get values of difference functions
      ifders = 0
      call diffslogbk_fast(r,beta,rscale1,diffs,ifders,ders,kvec,
     1     nterms)

c     form terms for shifting 

      ktemp(0) = kvec(0)
      difftemp(0) = diffs(0)
      powtemp(0) = pow(0)
      zmul=-exp(ima*theta)
      ztemp1= zmul
      do j = 1,nterms
         ktemp( j) = ztemp1*kvec(j)
         difftemp(j) = ztemp1*diffs(j)
         powtemp(j) = ztemp1*pow(j)
         ztemp1= ztemp1*zmul
      enddo

      do iii = 1,nexp

         do i = 0,nterms2
            mbhloc(i,iii) = 0.0d0
            lloc(i,iii) = 0.0d0
         enddo
         
c     shift

         mbhloc(0,iii) = mbhloc(0,iii) + dreal(mbhmpole(0,iii))*ktemp(0)
         mbhloc(0,iii) = mbhloc(0,iii) + dreal(ympole(0,iii))*ktemp(0)
         lloc(0,iii) = lloc(0,iii) + dreal(mbhmpole(0,iii))*difftemp(0)
         lloc(0,iii) = lloc(0,iii) + dreal(ympole(0,iii))*ktemp(0)
         do j = 1,nterms1
            mbhloc(0,iii) = mbhloc(0,iii)
     1           +(dreal(mbhmpole(j,iii))*dreal(ktemp(j))
     1           +dimag(mbhmpole(j,iii))*dimag(ktemp(j)))
            mbhloc(0,iii) = mbhloc(0,iii)
     1           +(dreal(ympole(j,iii))*dreal(ktemp(j))
     1           +dimag(ympole(j,iii))*dimag(ktemp(j)))
            lloc(0,iii) = lloc(0,iii)
     1           +(dreal(mbhmpole(j,iii))*dreal(difftemp(j))
     1           +dimag(mbhmpole(j,iii))*dimag(difftemp(j)))
            lloc(0,iii) = lloc(0,iii)
     1           +(dreal(ympole(j,iii))*dreal(ktemp(j))
     1           +dimag(ympole(j,iii))*dimag(ktemp(j)))
         enddo
c     
         rsi=rscale1
         rsi2=rscale1**2
         rsi5=rscale2/rscale1
         rsi52=rsi5*rscale2/rscale1
         isign = -1
         do i = 1,nterms2
            mbhloc(i,iii) = mbhloc(i,iii) 
     1           + 2.0d0*(dreal(mbhmpole(0,iii))*dreal(ktemp(i))
     2           + ima*dreal(mbhmpole(0,iii))*dimag(ktemp(i)))
            mbhloc(i,iii) = mbhloc(i,iii) 
     1           + 2.0d0*(dreal(ympole(0,iii))*dreal(ktemp(i))
     2           +ima*dreal(ympole(0,iii))*dimag(ktemp(i)))
            lloc(i,iii) = lloc(i,iii) 
     1           + 2.0d0*(dreal(mbhmpole(0,iii))*dreal(difftemp(i))
     2           + ima*dreal(mbhmpole(0,iii))*dimag(difftemp(i)))
     3           *dfac2(i)
            lloc(i,iii) = lloc(i,iii) 
     1           + 2.0d0*(dreal(ympole(0,iii))*dreal(ktemp(i))
     2           +ima*dreal(ympole(0,iii))*dimag(ktemp(i)))*dfac2(i)

            rsj=rscale1
            rsj2=rscale1**2
            do j = 1,min(nterms1,i)
               mbhloc(i,iii) = mbhloc(i,iii)
     1              +(dreal(mbhmpole(j,iii))*dreal(ktemp(i-j))
     1              -dimag(mbhmpole(j,iii))*dimag(ktemp(i-j))
     2              +ima*(dreal(mbhmpole(j,iii))*dimag(ktemp(i-j))
     3              +dimag(mbhmpole(j,iii))*dreal(ktemp(i-j))))*rsj2
               mbhloc(i,iii) = mbhloc(i,iii)
     1              +(dreal(mbhmpole(j,iii))*dreal(ktemp(i+j))
     1              +dimag(mbhmpole(j,iii))*dimag(ktemp(i+j))
     2              +ima*(dreal(mbhmpole(j,iii))*dimag(ktemp(i+j))
     3              -dimag(mbhmpole(j,iii))*dreal(ktemp(i+j))))
               mbhloc(i,iii) = mbhloc(i,iii)
     1              +(dreal(ympole(j,iii))*dreal(ktemp(i-j))
     1              -dimag(ympole(j,iii))*dimag(ktemp(i-j))
     2              +ima*(dreal(ympole(j,iii))*dimag(ktemp(i-j))
     3              +dimag(ympole(j,iii))*dreal(ktemp(i-j))))*rsj2
               mbhloc(i,iii) = mbhloc(i,iii)
     1              +(dreal(ympole(j,iii))*dreal(ktemp(i+j))
     1              +dimag(ympole(j,iii))*dimag(ktemp(i+j))
     2              +ima*(dreal(ympole(j,iii))*dimag(ktemp(i+j))
     3              -dimag(ympole(j,iii))*dreal(ktemp(i+j))))
               lloc(i,iii) = lloc(i,iii)
     1              +(dreal(mbhmpole(j,iii))*dreal(ktemp(i-j))
     1              -dimag(mbhmpole(j,iii))*dimag(ktemp(i-j))
     2              +ima*(dreal(mbhmpole(j,iii))*dimag(ktemp(i-j))
     3              +dimag(mbhmpole(j,iii))*dreal(ktemp(i-j))))
     5              *rsj2*dfac2(i)
               lloc(i,iii) = lloc(i,iii)
     1              +(dreal(mbhmpole(j,iii))*dreal(difftemp(i+j))
     1              +dimag(mbhmpole(j,iii))*dimag(difftemp(i+j))
     2              +ima*(dreal(mbhmpole(j,iii))*dimag(difftemp(i+j))
     3              -dimag(mbhmpole(j,iii))*dreal(difftemp(i+j))))
     5              *dfac2(i)
               lloc(i,iii) = lloc(i,iii)
     1              +(dreal(ympole(j,iii))*dreal(ktemp(i-j))
     1              -dimag(ympole(j,iii))*dimag(ktemp(i-j))
     2              +ima*(dreal(ympole(j,iii))*dimag(ktemp(i-j))
     3              +dimag(ympole(j,iii))*dreal(ktemp(i-j))))
     5              *rsj2*dfac2(i)
               lloc(i,iii) = lloc(i,iii)
     1              +(dreal(ympole(j,iii))*dreal(ktemp(i+j))
     1              +dimag(ympole(j,iii))*dimag(ktemp(i+j))
     2              +ima*(dreal(ympole(j,iii))*dimag(ktemp(i+j))
     3              -dimag(ympole(j,iii))*dreal(ktemp(i+j))))*dfac2(i)

               rsj=rsj*rscale1
               rsj2=rsj2*rscale1**2
            enddo
            do j = i+1,nterms1
               mbhloc(i,iii) = mbhloc(i,iii)
     1              +(dreal(mbhmpole(j,iii))*dreal(ktemp(j-i))
     1              +dimag(mbhmpole(j,iii))*dimag(ktemp(j-i))
     2              -ima*(dreal(mbhmpole(j,iii))*dimag(ktemp(j-i))
     3              -dimag(mbhmpole(j,iii))*dreal(ktemp(j-i))))*rsi2
               mbhloc(i,iii) = mbhloc(i,iii)
     1              +(dreal(mbhmpole(j,iii))*dreal(ktemp(i+j))
     1              +dimag(mbhmpole(j,iii))*dimag(ktemp(i+j))
     2              +ima*(dreal(mbhmpole(j,iii))*dimag(ktemp(i+j))
     3              -dimag(mbhmpole(j,iii))*dreal(ktemp(i+j))))
               mbhloc(i,iii) = mbhloc(i,iii)
     1              +(dreal(ympole(j,iii))*dreal(ktemp(j-i))
     1              +dimag(ympole(j,iii))*dimag(ktemp(j-i))
     2              -ima*(dreal(ympole(j,iii))*dimag(ktemp(j-i))
     3              -dimag(ympole(j,iii))*dreal(ktemp(j-i))))*rsi2
               mbhloc(i,iii) = mbhloc(i,iii)
     1              +(dreal(ympole(j,iii))*dreal(ktemp(i+j))
     1              +dimag(ympole(j,iii))*dimag(ktemp(i+j))
     2              +ima*(dreal(ympole(j,iii))*dimag(ktemp(i+j))
     3              -dimag(ympole(j,iii))*dreal(ktemp(i+j))))
               lloc(i,iii) = lloc(i,iii)
     1              +(dreal(mbhmpole(j,iii))*dreal(ktemp(j-i))
     1              +dimag(mbhmpole(j,iii))*dimag(ktemp(j-i))
     2              -ima*(dreal(mbhmpole(j,iii))*dimag(ktemp(j-i))
     3              -dimag(mbhmpole(j,iii))*dreal(ktemp(j-i))))
     5              *rsi2*dfac2(i)
               lloc(i,iii) = lloc(i,iii)
     1              +(dreal(mbhmpole(j,iii))*dreal(difftemp(i+j))
     1              +dimag(mbhmpole(j,iii))*dimag(difftemp(i+j))
     2              +ima*(dreal(mbhmpole(j,iii))*dimag(difftemp(i+j))
     3              -dimag(mbhmpole(j,iii))*dreal(difftemp(i+j))))
     5              *dfac2(i)
               lloc(i,iii) = lloc(i,iii)
     1              +(dreal(ympole(j,iii))*dreal(ktemp(j-i))
     1              +dimag(ympole(j,iii))*dimag(ktemp(j-i))
     2              -ima*(dreal(ympole(j,iii))*dimag(ktemp(j-i))
     3              -dimag(ympole(j,iii))*dreal(ktemp(j-i))))
     5              *rsi2*dfac2(i)
               lloc(i,iii) = lloc(i,iii)
     1              +(dreal(ympole(j,iii))*dreal(ktemp(i+j))
     1              +dimag(ympole(j,iii))*dimag(ktemp(i+j))
     2              +ima*(dreal(ympole(j,iii))*dimag(ktemp(i+j))
     3              -dimag(ympole(j,iii))*dreal(ktemp(i+j))))*dfac2(i)

            enddo
            mbhloc(i,iii)=mbhloc(i,iii)*isign
            lloc(i,iii)=lloc(i,iii)*isign
            isign = -isign
            rsi=rsi*rscale1
            rsi2=rsi2*rscale1**2
         enddo
         rsi=rscale2/rscale1
         do i = 1,nterms2
            mbhloc(i,iii) = mbhloc(i,iii)*rsi
            lloc(i,iii) = lloc(i,iii)*rsi
            rsi=rsi*rscale2/rscale1
         enddo

      enddo

      return
      end

      subroutine mbh2dmplocall_pre(beta,rscale1,nterms1,
     1     rscale2,nterms2,xlength,ntermsall,dfac2all,kvecall,diffsall,
     2     powall,ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT: 
c
c     rscale1           : scaling for original expansion
c     center1           : center for original expansion
c     coeffs1_1(nterms1): original coeffs, difference-type functions
c     coeffs2_1(nterms1): original coeffs, modified Bessel functions
c     nterms1           : number of terms in original expansion
c     rscale2           : scaling for new expansion
c     center2           : center for new expansion
c     rad2              : radius of circle to sample on 
c                         for new expansion
c     nterms2           : number of terms in original expansion
c     zcirc2(2,nterms2) : precomp'd evenly-spaced points on unit circle
c     beta              : the modified biharmonic parameter

c     work(*)           : work array. recommended length 16*nterms+500
c
c     OUTPUT:
c
c     coeffs1_2(nterms2): new coeffs, difference-type functions
c     coeffs2_2(nterms2): new coeffs, modified Bessel functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 rscale1, rscale2, beta, xlength
      real *8 dfac2all(0:ntermsall,-3:3,-3:3)
      real *8 kvecall(0:ntermsall,-3:3,-3:3)
      real *8 diffsall(0:ntermsall,-3:3,-3:3)
      real *8 powall(0:ntermsall,-3:3,-3:3)
      integer nterms1, nterms2, ntermsall, ier
c     local
      integer i,j,k
      real *8 center1(2), center2(2)
      
      ier = 0

      if (nterms1+nterms2 .gt. ntermsall) then
         ier = 1
         return
      endif

      do k = -3,3
      do j = -3,3
      do i = 0,ntermsall
         dfac2all(i,j,k) = 0.0d0
         kvecall(i,j,k) = 0.0d0
         diffsall(i,j,k) = 0.0d0
         powall(i,j,k) = 0.0d0
      enddo
      enddo
      enddo
      
      center1(1) = 0.0d0
      center1(2) = 0.0d0

      do k = 0,3
         do j = 0,3
            if (j .gt. 1 .or. k .gt. 1) then
               center2(1) = xlength*j
               center2(2) = xlength*k
               call mbh2dmploc_pre(beta,rscale1,center1,nterms1,
     1              rscale2,center2,nterms2,dfac2all(0,j,k),
     1              kvecall(0,j,k),diffsall(0,j,k),powall(0,j,k))
            endif
         enddo
      enddo

      do k = -3,3
         do j = -3,3
            if (k .lt. 0 .and. j .lt. 0) then
               do i = 0,ntermsall
                  dfac2all(i,j,k) = dfac2all(i,-j,-k)
                  kvecall(i,j,k) = kvecall(i,-j,-k)
                  diffsall(i,j,k) = diffsall(i,-j,-k)
                  powall(i,j,k) = powall(i,-j,-k)
               enddo
            else if (k .lt. 0 .and. j .ge. 0) then
               do i = 0,ntermsall
                  dfac2all(i,j,k) = dfac2all(i,j,-k)
                  kvecall(i,j,k) = kvecall(i,j,-k)
                  diffsall(i,j,k) = diffsall(i,j,-k)
                  powall(i,j,k) = powall(i,j,-k)
               enddo
            else if (j .lt. 0 .and. k .ge. 0) then
               do i = 0,ntermsall
                  dfac2all(i,j,k) = dfac2all(i,-j,k)
                  kvecall(i,j,k) = kvecall(i,-j,k)
                  diffsall(i,j,k) = diffsall(i,-j,k)
                  powall(i,j,k) = powall(i,-j,k)
               enddo
            endif
         enddo
      enddo
               

      return
      end

      subroutine mbh2dmplocall(beta,rscale,xlength,nterms,mbhmp,ymp,
     1     mbhloc,lloc,bigmbhloc,biglloc,
     3     ntermsall,dfac2all,kvecall,diffsall,powall,
     4     jbox,inall,nnall,iynall,in34,nn34,iy34,isall,nsall,
     5     iysall,is12,ns12,iy12,ieall,neall,iyeall,ie14,ne14,iy14,
     6     iwall,nwall,iywall,iw23,nw23,iy23,iww3,iwy3,nww3,iww2,
     7     iwy2,nww2,iee4,iey4,nee4,iee1,iey1,nee1,inbig34,isbig12,
     8     iebig14,iwbig23,iebig4,iwbig3,iebig1,iwbig2,ichildbox,
     9     localonoff,isrcflag)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     processing for the interaction list...
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, rscale, xlength
      integer nterms, ntermsall, localonoff(*), jbox, ichildbox(4,*)
      integer isrcflag(*)
c     lists info
      integer inall(6),nnall,iynall(6)
      integer isall(6),nsall,iysall(6)
      integer ieall(4),neall,iyeall(4)
      integer iwall(4),nwall,iywall(4)
      integer in34(2),nn34,iy34(2)
      integer is12(2),ns12,iy12(2)
      integer iw23(2),nw23,iy23(2)
      integer ie14(2),ne14,iy14(2)
      integer iee4(1),nee4,iey4(1)
      integer iee1(1),nee1,iey1(1)
      integer iww2(1),nww2,iwy2(1)
      integer iww3(1),nww3,iwy3(1)
      integer inbig34(3), isbig12(3)
      integer iebig14(1), iwbig23(1)
      integer iebig4(1), iwbig3(1)
      integer iebig1(1), iwbig2(1)
c     precomputes
      real *8 dfac2all(0:ntermsall,-3:3,-3:3)
      real *8 kvecall(0:ntermsall,-3:3,-3:3)
      real *8 diffsall(0:ntermsall,-3:3,-3:3)
      real *8 powall(0:ntermsall,-3:3,-3:3)
c     multipoles and locals
      complex *16 ymp(0:nterms,*),mbhmp(0:nterms,*)
      complex *16 lloc(0:nterms,*),mbhloc(0:nterms,*)
      complex *16 biglloc(0:nterms,4,*)
      complex *16 bigmbhloc(0:nterms,4,*)
c     local variables
      integer ic1,ic2,ic3,ic4,ibox,i,j,ibigbox,ioff4big,joff4big
      integer isrc1,isrc2,isrc3,isrc4
      integer ioff1, joff1, ioff2, joff2, ioff3, joff3, ioff4, joff4
      integer ioffghosts(4), joffghosts(4)
      real *8 center1(2), center2(2)
      complex *16 mbhloctemp(0:200), lloctemp(0:200)
c     position of ghost children relative to 4th ghost child
      data ioffghosts / 0, 1, 1, 0/
      data joffghosts / 1, 1, 0, 0/
      

      ic1 = ichildbox(1,jbox)
      ic2 = ichildbox(2,jbox)
      ic3 = ichildbox(3,jbox)
      ic4 = ichildbox(4,jbox)

      isrc1 = isrcflag(ic1)
      isrc2 = isrcflag(ic2)
      isrc3 = isrcflag(ic3)
      isrc4 = isrcflag(ic4)

      center1(1) = 0.0d0
      center1(2) = 0.0d0

c     process each child for same size boxes
c     ---> for each listed box find offsets 
c          for each child (based on 4th)
c     ---> send each child to that box

      do i = 1,nnall
         ibox = inall(i)
         ioff4 = iynall(i)
         joff4 = 3
         ioff1 = ioff4
         joff1 = joff4-1
         ioff2 = ioff4-1
         joff2 = joff4-1
         ioff3 = ioff4-1
         joff3 = joff4

         if (localonoff(ibox) .eq. 1) then

            center2(1) = ioff1*xlength
            center2(2) = joff1*xlength
            
            if (isrc1 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic1),
     1              ymp(0,ic1),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff1,joff1),
     3              kvecall(0,ioff1,joff1),diffsall(0,ioff1,joff1),
     4              powall(0,ioff1,joff1))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif
            
            center2(1) = ioff2*xlength
            center2(2) = joff2*xlength

            if (isrc2 .eq. 1) then

               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic2),
     1              ymp(0,ic2),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff2,joff2),
     3              kvecall(0,ioff2,joff2),diffsall(0,ioff2,joff2),
     4              powall(0,ioff2,joff2))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif

            center2(1) = ioff3*xlength
            center2(2) = joff3*xlength
            
            if (isrc3 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic3),
     1              ymp(0,ic3),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff3,joff3),
     3              kvecall(0,ioff3,joff3),diffsall(0,ioff3,joff3),
     4              powall(0,ioff3,joff3))
               
               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif

            center2(1) = ioff4*xlength
            center2(2) = joff4*xlength

            if (isrc4 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic4),
     1              ymp(0,ic4),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff4,joff4),
     3              kvecall(0,ioff4,joff4),diffsall(0,ioff4,joff4),
     4              powall(0,ioff4,joff4))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif
         endif
      enddo

      do i = 1,nsall
         ibox = isall(i)
         ioff4 = iysall(i)
         joff4 = -2
         ioff1 = ioff4
         joff1 = joff4-1
         ioff2 = ioff4-1
         joff2 = joff4-1
         ioff3 = ioff4-1
         joff3 = joff4

         if (localonoff(ibox) .eq. 1) then

            center2(1) = ioff1*xlength
            center2(2) = joff1*xlength

            if (isrc1 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic1),
     1              ymp(0,ic1),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff1,joff1),
     3              kvecall(0,ioff1,joff1),diffsall(0,ioff1,joff1),
     4              powall(0,ioff1,joff1))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif
            
            center2(1) = ioff2*xlength
            center2(2) = joff2*xlength

            if (isrc2 .eq. 1) then

               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic2),
     1              ymp(0,ic2),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff2,joff2),
     3              kvecall(0,ioff2,joff2),diffsall(0,ioff2,joff2),
     4              powall(0,ioff2,joff2))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif

            center2(1) = ioff3*xlength
            center2(2) = joff3*xlength
            
            if (isrc3 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic3),
     1              ymp(0,ic3),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff3,joff3),
     3              kvecall(0,ioff3,joff3),diffsall(0,ioff3,joff3),
     4              powall(0,ioff3,joff3))
               
               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif

            center2(1) = ioff4*xlength
            center2(2) = joff4*xlength

            if (isrc4 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic4),
     1              ymp(0,ic4),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff4,joff4),
     3              kvecall(0,ioff4,joff4),diffsall(0,ioff4,joff4),
     4              powall(0,ioff4,joff4))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif
         endif
      enddo

      do i = 1,neall
         ibox = ieall(i)
         ioff4 = 3
         joff4 = iyeall(i)
         ioff1 = ioff4
         joff1 = joff4-1
         ioff2 = ioff4-1
         joff2 = joff4-1
         ioff3 = ioff4-1
         joff3 = joff4

         if (localonoff(ibox) .eq. 1) then

            center2(1) = ioff1*xlength
            center2(2) = joff1*xlength

            if (isrc1 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic1),
     1              ymp(0,ic1),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff1,joff1),
     3              kvecall(0,ioff1,joff1),diffsall(0,ioff1,joff1),
     4              powall(0,ioff1,joff1))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif
            
            center2(1) = ioff2*xlength
            center2(2) = joff2*xlength

            if (isrc2 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic2),
     1              ymp(0,ic2),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff2,joff2),
     3              kvecall(0,ioff2,joff2),diffsall(0,ioff2,joff2),
     4              powall(0,ioff2,joff2))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif

            center2(1) = ioff3*xlength
            center2(2) = joff3*xlength
            
            if (isrc3 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic3),
     1              ymp(0,ic3),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff3,joff3),
     3              kvecall(0,ioff3,joff3),diffsall(0,ioff3,joff3),
     4              powall(0,ioff3,joff3))
               
               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif

            center2(1) = ioff4*xlength
            center2(2) = joff4*xlength

            if (isrc4 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic4),
     1              ymp(0,ic4),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff4,joff4),
     3              kvecall(0,ioff4,joff4),diffsall(0,ioff4,joff4),
     4              powall(0,ioff4,joff4))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif
         endif
      enddo

      do i = 1,nwall
         ibox = iwall(i)
         ioff4 = -2
         joff4 = iywall(i)
         ioff1 = ioff4
         joff1 = joff4-1
         ioff2 = ioff4-1
         joff2 = joff4-1
         ioff3 = ioff4-1
         joff3 = joff4

         if (localonoff(ibox) .eq. 1) then

            center2(1) = ioff1*xlength
            center2(2) = joff1*xlength

            if (isrc1 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic1),
     1              ymp(0,ic1),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff1,joff1),
     3              kvecall(0,ioff1,joff1),diffsall(0,ioff1,joff1),
     4              powall(0,ioff1,joff1))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif
            
            center2(1) = ioff2*xlength
            center2(2) = joff2*xlength

            if (isrc2 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic2),
     1              ymp(0,ic2),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff2,joff2),
     3              kvecall(0,ioff2,joff2),diffsall(0,ioff2,joff2),
     4              powall(0,ioff2,joff2))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif

            center2(1) = ioff3*xlength
            center2(2) = joff3*xlength
            
            if (isrc3 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic3),
     1              ymp(0,ic3),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff3,joff3),
     3              kvecall(0,ioff3,joff3),diffsall(0,ioff3,joff3),
     4              powall(0,ioff3,joff3))
               
               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif

            center2(1) = ioff4*xlength
            center2(2) = joff4*xlength

            if (isrc4 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic4),
     1              ymp(0,ic4),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff4,joff4),
     3              kvecall(0,ioff4,joff4),diffsall(0,ioff4,joff4),
     4              powall(0,ioff4,joff4))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif
         endif
      enddo

      do i = 1,nn34
         ibox = in34(i)
         ioff4 = iy34(i)
         joff4 = 2
         ioff3 = ioff4-1
         joff3 = joff4

         if (localonoff(ibox) .eq. 1) then

            center2(1) = ioff3*xlength
            center2(2) = joff3*xlength
            
            if (isrc3 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic3),
     1              ymp(0,ic3),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff3,joff3),
     3              kvecall(0,ioff3,joff3),diffsall(0,ioff3,joff3),
     4              powall(0,ioff3,joff3))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif
            
            center2(1) = ioff4*xlength
            center2(2) = joff4*xlength

            if (isrc4 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic4),
     1              ymp(0,ic4),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff4,joff4),
     3              kvecall(0,ioff4,joff4),diffsall(0,ioff4,joff4),
     4              powall(0,ioff4,joff4))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif
         endif
      enddo

      do i = 1,ns12
         ibox = is12(i)
         ioff4 = iy12(i)
         joff4 = -1
         ioff1 = ioff4
         joff1 = joff4-1
         ioff2 = ioff4-1
         joff2 = joff4-1

         if (localonoff(ibox) .eq. 1) then

            center2(1) = ioff1*xlength
            center2(2) = joff1*xlength
            
            if (isrc1 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic1),
     1              ymp(0,ic1),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff1,joff1),
     3              kvecall(0,ioff1,joff1),diffsall(0,ioff1,joff1),
     4              powall(0,ioff1,joff1))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif

            center2(1) = ioff2*xlength
            center2(2) = joff2*xlength
            
            if (isrc2 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic2),
     1              ymp(0,ic2),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff2,joff2),
     3              kvecall(0,ioff2,joff2),diffsall(0,ioff2,joff2),
     4              powall(0,ioff2,joff2))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif
         endif         
      enddo

      do i = 1,nw23
         ibox = iw23(i)
         ioff4 = -1
         joff4 = iy23(i)
         ioff2 = ioff4-1
         joff2 = joff4-1
         ioff3 = ioff4-1
         joff3 = joff4

         if (localonoff(ibox) .eq. 1) then

            center2(1) = ioff2*xlength
            center2(2) = joff2*xlength
            
            if (isrc2 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic2),
     1              ymp(0,ic2),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff2,joff2),
     3              kvecall(0,ioff2,joff2),diffsall(0,ioff2,joff2),
     4              powall(0,ioff2,joff2))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif
            
            center2(1) = ioff3*xlength
            center2(2) = joff3*xlength

            if (isrc3 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic3),
     1              ymp(0,ic3),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff3,joff3),
     3              kvecall(0,ioff3,joff3),diffsall(0,ioff3,joff3),
     4              powall(0,ioff3,joff3))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif
         endif         
      enddo

      do i = 1,ne14
         ibox = ie14(i)
         ioff4 = 2
         joff4 = iy14(i)
         ioff1 = ioff4
         joff1 = joff4-1

         if (localonoff(ibox) .eq. 1) then

            center2(1) = ioff1*xlength
            center2(2) = joff1*xlength
            
            if (isrc1 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic1),
     1              ymp(0,ic1),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff1,joff1),
     3              kvecall(0,ioff1,joff1),diffsall(0,ioff1,joff1),
     4              powall(0,ioff1,joff1))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif
            
            center2(1) = ioff4*xlength
            center2(2) = joff4*xlength
            
            if (isrc4 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic4),
     1              ymp(0,ic4),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff4,joff4),
     3              kvecall(0,ioff4,joff4),diffsall(0,ioff4,joff4),
     4              powall(0,ioff4,joff4))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif
         endif
      enddo

      do i = 1,nee4
         ibox = iee4(i)
         ioff4 = 2
         joff4 = iey4(i)

         if (localonoff(ibox) .eq. 1) then
            
            center2(1) = ioff4*xlength
            center2(2) = joff4*xlength
            
            if (isrc4 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic4),
     1              ymp(0,ic4),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff4,joff4),
     3              kvecall(0,ioff4,joff4),diffsall(0,ioff4,joff4),
     4              powall(0,ioff4,joff4))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif
         endif
      enddo

      do i = 1,nee1
         ibox = iee1(i)
         ioff4 = 2
         joff4 = iey1(i)
         ioff1 = ioff4
         joff1 = joff4-1

         if (localonoff(ibox) .eq. 1) then         

            center2(1) = ioff1*xlength
            center2(2) = joff1*xlength
            
            if (isrc1 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic1),
     1              ymp(0,ic1),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff1,joff1),
     3              kvecall(0,ioff1,joff1),diffsall(0,ioff1,joff1),
     4              powall(0,ioff1,joff1))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif
         endif
      enddo

      do i = 1,nww2
         ibox = iww2(i)
         ioff4 = -1
         joff4 = iwy2(i)
         ioff2 = ioff4-1
         joff2 = joff4-1
         
         if (localonoff(ibox) .eq. 1) then         

            center2(1) = ioff2*xlength
            center2(2) = joff2*xlength
            
            if (isrc2 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic2),
     1              ymp(0,ic2),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff2,joff2),
     3              kvecall(0,ioff2,joff2),diffsall(0,ioff2,joff2),
     4              powall(0,ioff2,joff2))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif
         endif

      enddo

      do i = 1,nww3
         ibox = iww3(i)
         ioff4 = -1
         joff4 = iwy3(i)
         ioff3 = ioff4-1
         joff3 = joff4

         if (localonoff(ibox) .eq. 1) then

            center2(1) = ioff3*xlength
            center2(2) = joff3*xlength
            
            if (isrc3 .eq. 1) then
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic3),
     1              ymp(0,ic3),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff3,joff3),
     3              kvecall(0,ioff3,joff3),diffsall(0,ioff3,joff3),
     4              powall(0,ioff3,joff3))

               call mbh2d_addexp(mbhloctemp,mbhloc(0,ibox),nterms)
               call mbh2d_addexp(lloctemp,lloc(0,ibox),nterms)
            endif
         endif
      enddo

c     process each child for bigger sized boxes
c     ---> create 4 ghost children for each
c          listed bigger box
c     ---> for each listed bigger box find offsets 
c          for each child to 4th ghost (based on 4th)
c     ---> use offsets to 4th ghost to figure out offsets
c          to other ghosts
c     ---> send each child to each of the 4 ghost children

      do i = 1,3
         ibigbox = inbig34(i)
         if (ibigbox .gt. 0 .and. localonoff(ibigbox) .eq. 1) then
            ioff4big = -2 + 2*(i-1)
            joff4big = 2
c     loop over ghost children
            do j = 1,4
               ioff4 = ioff4big + ioffghosts(j)
               joff4 = joff4big + joffghosts(j)
               ioff1 = ioff4
               joff1 = joff4-1
               ioff2 = ioff4-1
               joff2 = joff4-1
               ioff3 = ioff4-1
               joff3 = joff4
               
               if (i .eq. 3 .and. isrc1 .eq. 1) then

                  center2(1) = ioff1*xlength
                  center2(2) = joff1*xlength
                  
                  call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic1),
     1                 ymp(0,ic1),nterms,rscale,center2,mbhloctemp,
     2                 lloctemp,nterms,dfac2all(0,ioff1,joff1),
     3                 kvecall(0,ioff1,joff1),diffsall(0,ioff1,joff1),
     4                 powall(0,ioff1,joff1))

                  call mbh2d_addexp(mbhloctemp,bigmbhloc(0,j,ibigbox),
     1                 nterms)
                  call mbh2d_addexp(lloctemp,biglloc(0,j,ibigbox),
     1                 nterms)
               endif

               if (i .eq. 1 .and. isrc2 .eq. 1) then

                  center2(1) = ioff2*xlength
                  center2(2) = joff2*xlength
                  
                  call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic2),
     1                 ymp(0,ic2),nterms,rscale,center2,mbhloctemp,
     2                 lloctemp,nterms,dfac2all(0,ioff2,joff2),
     3                 kvecall(0,ioff2,joff2),diffsall(0,ioff2,joff2),
     4                 powall(0,ioff2,joff2))

                  call mbh2d_addexp(mbhloctemp,bigmbhloc(0,j,ibigbox),
     1                 nterms)
                  call mbh2d_addexp(lloctemp,biglloc(0,j,ibigbox),
     1                 nterms)

               endif
                  
               center2(1) = ioff3*xlength
               center2(2) = joff3*xlength
               
               if (isrc3 .eq. 1) then

                  call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic3),
     1                 ymp(0,ic3),nterms,rscale,center2,mbhloctemp,
     2                 lloctemp,nterms,dfac2all(0,ioff3,joff3),
     3                 kvecall(0,ioff3,joff3),diffsall(0,ioff3,joff3),
     4                 powall(0,ioff3,joff3))

                  call mbh2d_addexp(mbhloctemp,bigmbhloc(0,j,ibigbox),
     1                 nterms)
                  call mbh2d_addexp(lloctemp,biglloc(0,j,ibigbox),
     1                 nterms)
               endif
               
               center2(1) = ioff4*xlength
               center2(2) = joff4*xlength
               
               if (isrc4 .eq. 1) then
                  call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic4),
     1                 ymp(0,ic4),nterms,rscale,center2,mbhloctemp,
     2                 lloctemp,nterms,dfac2all(0,ioff4,joff4),
     3                 kvecall(0,ioff4,joff4),diffsall(0,ioff4,joff4),
     4                 powall(0,ioff4,joff4))

                  call mbh2d_addexp(mbhloctemp,bigmbhloc(0,j,ibigbox),
     1                 nterms)
                  call mbh2d_addexp(lloctemp,biglloc(0,j,ibigbox),
     1                 nterms)
               endif

            enddo
         endif
      enddo
               
      do i = 1,3
         ibigbox = isbig12(i)
         if (ibigbox .gt. 0 .and. localonoff(ibigbox) .eq. 1) then
            ioff4big = -2 + 2*(i-1)
            joff4big = -2
c     loop over ghost children
            do j = 1,4
               ioff4 = ioff4big + ioffghosts(j)
               joff4 = joff4big + joffghosts(j)
               ioff1 = ioff4
               joff1 = joff4-1
               ioff2 = ioff4-1
               joff2 = joff4-1
               ioff3 = ioff4-1
               joff3 = joff4
               
               if (isrc1 .eq. 1) then
                  center2(1) = ioff1*xlength
                  center2(2) = joff1*xlength
                  
                  call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic1),
     1                 ymp(0,ic1),nterms,rscale,center2,mbhloctemp,
     2                 lloctemp,nterms,dfac2all(0,ioff1,joff1),
     3                 kvecall(0,ioff1,joff1),diffsall(0,ioff1,joff1),
     4                 powall(0,ioff1,joff1))

                  call mbh2d_addexp(mbhloctemp,bigmbhloc(0,j,ibigbox),
     1                 nterms)
                  call mbh2d_addexp(lloctemp,biglloc(0,j,ibigbox),
     1                 nterms)
               endif

               if (isrc2 .eq. 1) then
                  center2(1) = ioff2*xlength
                  center2(2) = joff2*xlength
                  
                  call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic2),
     1                 ymp(0,ic2),nterms,rscale,center2,mbhloctemp,
     2                 lloctemp,nterms,dfac2all(0,ioff2,joff2),
     3                 kvecall(0,ioff2,joff2),diffsall(0,ioff2,joff2),
     4                 powall(0,ioff2,joff2))

                  call mbh2d_addexp(mbhloctemp,bigmbhloc(0,j,ibigbox),
     1                 nterms)
                  call mbh2d_addexp(lloctemp,biglloc(0,j,ibigbox),
     1                 nterms)
               endif

               if (i .eq. 1 .and. isrc3 .eq. 1) then
                  
                  center2(1) = ioff3*xlength
                  center2(2) = joff3*xlength
                  
                  call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic3),
     1                 ymp(0,ic3),nterms,rscale,center2,mbhloctemp,
     2                 lloctemp,nterms,dfac2all(0,ioff3,joff3),
     3                 kvecall(0,ioff3,joff3),diffsall(0,ioff3,joff3),
     4                 powall(0,ioff3,joff3))

                  call mbh2d_addexp(mbhloctemp,bigmbhloc(0,j,ibigbox),
     1                 nterms)
                  call mbh2d_addexp(lloctemp,biglloc(0,j,ibigbox),
     1                 nterms)

               endif

               if (i .eq. 3 .and. isrc4 .eq. 1) then

                  center2(1) = ioff4*xlength
                  center2(2) = joff4*xlength
                  
                  call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic4),
     1                 ymp(0,ic4),nterms,rscale,center2,mbhloctemp,
     2                 lloctemp,nterms,dfac2all(0,ioff4,joff4),
     3                 kvecall(0,ioff4,joff4),diffsall(0,ioff4,joff4),
     4                 powall(0,ioff4,joff4))

                  call mbh2d_addexp(mbhloctemp,bigmbhloc(0,j,ibigbox),
     1                 nterms)
                  call mbh2d_addexp(lloctemp,biglloc(0,j,ibigbox),
     1                 nterms)

               endif
            enddo
         endif
      enddo
               
      ibigbox = iwbig23(1)
      if (ibigbox .gt. 0 .and. localonoff(ibigbox) .eq. 1) then
         ioff4big = -2
         joff4big = 0
c     loop over ghost children
         do j = 1,4
            ioff4 = ioff4big + ioffghosts(j)
            joff4 = joff4big + joffghosts(j)
            ioff1 = ioff4
            joff1 = joff4-1
            ioff2 = ioff4-1
            joff2 = joff4-1
            ioff3 = ioff4-1
            joff3 = joff4
            
            if (isrc2 .eq. 1) then
               center2(1) = ioff2*xlength
               center2(2) = joff2*xlength
               
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic2),
     1              ymp(0,ic2),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff2,joff2),
     3              kvecall(0,ioff2,joff2),diffsall(0,ioff2,joff2),
     4              powall(0,ioff2,joff2))

               call mbh2d_addexp(mbhloctemp,bigmbhloc(0,j,ibigbox),
     1              nterms)
               call mbh2d_addexp(lloctemp,biglloc(0,j,ibigbox),
     1              nterms)
            endif

            if (isrc3 .eq. 1) then
               center2(1) = ioff3*xlength
               center2(2) = joff3*xlength
               
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic3),
     1              ymp(0,ic3),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff3,joff3),
     3              kvecall(0,ioff3,joff3),diffsall(0,ioff3,joff3),
     4              powall(0,ioff3,joff3))

               call mbh2d_addexp(mbhloctemp,bigmbhloc(0,j,ibigbox),
     1              nterms)
               call mbh2d_addexp(lloctemp,biglloc(0,j,ibigbox),
     1              nterms)
            endif

         enddo
      endif
               
      ibigbox = iebig14(1)
      if (ibigbox .gt. 0 .and. localonoff(ibigbox) .eq. 1) then
         ioff4big = 2
         joff4big = 0
c     loop over ghost children
         do j = 1,4
            ioff4 = ioff4big + ioffghosts(j)
            joff4 = joff4big + joffghosts(j)
            ioff1 = ioff4
            joff1 = joff4-1
            ioff2 = ioff4-1
            joff2 = joff4-1
            ioff3 = ioff4-1
            joff3 = joff4
            
            if (isrc1 .eq. 1) then
               center2(1) = ioff1*xlength
               center2(2) = joff1*xlength
               
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic1),
     1              ymp(0,ic1),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff1,joff1),
     3              kvecall(0,ioff1,joff1),diffsall(0,ioff1,joff1),
     4              powall(0,ioff1,joff1))

               call mbh2d_addexp(mbhloctemp,bigmbhloc(0,j,ibigbox),
     1              nterms)
               call mbh2d_addexp(lloctemp,biglloc(0,j,ibigbox),
     1              nterms)
            endif

            if (isrc4 .eq. 1) then
               center2(1) = ioff4*xlength
               center2(2) = joff4*xlength
               
               call mbh2dmploc_wpre(beta,rscale,center1,mbhmp(0,ic4),
     1              ymp(0,ic4),nterms,rscale,center2,mbhloctemp,
     2              lloctemp,nterms,dfac2all(0,ioff4,joff4),
     3              kvecall(0,ioff4,joff4),diffsall(0,ioff4,joff4),
     4              powall(0,ioff4,joff4))

               call mbh2d_addexp(mbhloctemp,bigmbhloc(0,j,ibigbox),
     1              nterms)
               call mbh2d_addexp(lloctemp,biglloc(0,j,ibigbox),
     1              nterms)
            endif

         enddo
      endif
               

      return
      end


      subroutine mbh2dlocloc(beta,rscale1,center1,mbhloc1,lloc1,
     1     nterms1,rscale2,center2,mbhloc2,lloc2,nterms2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT: 
c
c     rscale1           : scaling for original expansion
c     center1           : center for original expansion
c     coeffs1_1(nterms1): original coeffs, difference-type functions
c     coeffs2_1(nterms1): original coeffs, modified Bessel functions
c     nterms1           : number of terms in original expansion
c     rscale2           : scaling for new expansion
c     center2           : center for new expansion
c     rad2              : radius of circle to sample on 
c                         for new expansion
c     nterms2           : number of terms in original expansion
c     zcirc2(2,nterms2) : precomp'd evenly-spaced points on unit circle
c     beta              : the modified biharmonic parameter

c     work(*)           : work array. recommended length 16*nterms+500
c
c     OUTPUT:
c
c     coeffs1_2(nterms2): new coeffs, difference-type functions
c     coeffs2_2(nterms2): new coeffs, modified Bessel functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 rscale1, center1(2), beta, rscale2, center2(2)
      complex *16 mbhloc1(0:nterms1), lloc1(0:nterms1)
      complex *16 mbhloc2(0:nterms2), lloc2(0:nterms2)
      integer nterms1, nterms2
c     local
      complex *16 z,ima, zmul,zinv,ztemp1,ztemp2
      complex *16 jtemp(0:200), difftemp(0:200), powtemp(0:200)
      complex *16 pow2temp(0:200)
      real *8 zdiff(2), r, theta, pi, done
      real *8 rsj, rsi, rsi7, rsi5, fs2
      real *8 ival(0:200), diffs(0:200), ders(0:200)
      real *8 pow(0:200), dpow(0:200), pow2(0:200)
      real *8 dfac(0:200), dfac2(0:200), carray(0:200,0:200)
      integer i, j, m, l, nterms, ifder
c     
      data ima/(0.0d0,1.0d0)/
c     
      done=1
      pi=4*atan(done)
c     
      nterms = nterms1+nterms2
c     
      zdiff(1)=center2(1)-center1(1)
      zdiff(2)=center2(2)-center1(2)
      call y2cart2polar(zdiff,r,theta)
      theta=theta-pi

      ifder = 0
      call diffszkik_fast(r,beta,rscale1,diffs,ifder,ders,ival,
     1     nterms)

      call mbh2d_rk(pow,dpow,r,beta,rscale1,nterms)
      call mbh2d_rksc(pow2,dpow,r,beta,rscale1,nterms)

      dfac2(0) = 1.0d0
      dfac(0) = 1.0d0

      do i = 1,nterms
         dfac2(i) = dfac2(i-1)/(2.0d0*i)
         dfac(i) = dfac(i-1)*i
      enddo

      do l = 0,nterms
         carray(l,0) = 1.0d0
      enddo
      do m=1,nterms
         carray(m,m) = 1.0d0
         do l=m+1,nterms
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
         enddo
      enddo

      do i = 0,nterms2
         mbhloc2(i) = 0
         lloc2(i) = 0
      enddo
c     
      jtemp(0) = ival(0)
      powtemp(0) = pow(0)
      pow2temp(0) = pow2(0)
      difftemp(0) = diffs(0)
      zmul=-exp(ima*theta)
      ztemp1= zmul
      do j = 1,nterms
         jtemp( j) = ztemp1*ival(j)
         difftemp(j) = ztemp1*diffs(j)
         powtemp(j) = ztemp1*pow(j)
         pow2temp(j) = ztemp1*pow2(j)
         ztemp1= ztemp1*zmul
      enddo
c     
      mbhloc2(0) = mbhloc2(0) + mbhloc1(0)*jtemp(0)
      lloc2(0) = lloc2(0) + mbhloc1(0)*difftemp(0)
      lloc2(0) = lloc2(0) + lloc1(0)*powtemp(0)
      rsj=rscale1
      do j = 1,nterms1
         mbhloc2(0) = mbhloc2(0)+(dreal(mbhloc1(j))*dreal(jtemp(j))
     1        +dimag(mbhloc1(j))*dimag(jtemp(j)))
         lloc2(0) = lloc2(0)+(dreal(mbhloc1(j))*dreal(difftemp(j))
     1        +dimag(mbhloc1(j))*dimag(difftemp(j)))
         lloc2(0) = lloc2(0)+(dreal(lloc1(j))*dreal(powtemp(j))
     1        +dimag(lloc1(j))*dimag(powtemp(j)))
      enddo
c     
      rsi=rscale1
      rsi7=rscale2
      rsi5=rscale2/rscale1
      do i = 1,nterms2
         mbhloc2(i) = mbhloc2(i) 
     1        + 2.0d0*(dreal(mbhloc1(0))*dreal(jtemp(i))
     2        +ima*dreal(mbhloc1(0))*dimag(jtemp(i)))*rsi7*rsi
         lloc2(i) = lloc2(i) 
     1        + 2.0d0*(dreal(mbhloc1(0))*dreal(jtemp(i))
     2        +ima*dreal(mbhloc1(0))*dimag(jtemp(i)))*rsi7*rsi*dfac2(i)
         fs2=rsi5*rscale1**2
         if( nterms1 .le. i-1 ) fs2=fs2*rscale1**(2*(i-1-nterms1))
         do j = min(nterms1,i-1),1,-1
            mbhloc2(i) = mbhloc2(i)+(dreal(mbhloc1(j))*dreal(jtemp(i-j))
     1           -dimag(mbhloc1(j))*dimag(jtemp(i-j))
     2           +ima*(dreal(mbhloc1(j))*dimag(jtemp(i-j))
     3           +dimag(mbhloc1(j))*dreal(jtemp(i-j))))*fs2
            mbhloc2(i) = mbhloc2(i)+(dreal(mbhloc1(j))*dreal(jtemp(i+j))
     1           +dimag(mbhloc1(j))*dimag(jtemp(i+j))
     2           +ima*(dreal(mbhloc1(j))*dimag(jtemp(i+j))
     3           -dimag(mbhloc1(j))*dreal(jtemp(i+j))))*rsi*rsi7
            lloc2(i) = lloc2(i)+(dreal(mbhloc1(j))*dreal(jtemp(i-j))
     1           -dimag(mbhloc1(j))*dimag(jtemp(i-j))
     2           +ima*(dreal(mbhloc1(j))*dimag(jtemp(i-j))
     3           +dimag(mbhloc1(j))*dreal(jtemp(i-j))))*fs2*dfac2(i)
            lloc2(i) = lloc2(i)+(dreal(mbhloc1(j))*dreal(jtemp(i+j))
     1           +dimag(mbhloc1(j))*dimag(jtemp(i+j))
     2           +ima*(dreal(mbhloc1(j))*dimag(jtemp(i+j))
     3           -dimag(mbhloc1(j))*dreal(jtemp(i+j))))
     4           *rsi*rsi7*dfac2(i)
            fs2=fs2*rscale1**2
         enddo
         do j = i,nterms1
            mbhloc2(i) = mbhloc2(i)+(dreal(mbhloc1(j))*dreal(jtemp(j-i))
     1           +dimag(mbhloc1(j))*dimag(jtemp(j-i))
     2           -ima*(dreal(mbhloc1(j))*dimag(jtemp(j-i))
     3           -dimag(mbhloc1(j))*dreal(jtemp(j-i))))*rsi5
            lloc2(i) = lloc2(i)+(dreal(mbhloc1(j))*dreal(difftemp(j-i))
     1           +dimag(mbhloc1(j))*dimag(difftemp(j-i))
     2           -ima*(dreal(mbhloc1(j))*dimag(difftemp(j-i))
     3           -dimag(mbhloc1(j))*dreal(difftemp(j-i))))*rsi5*dfac2(i)
            lloc2(i) = lloc2(i)+(dreal(lloc1(j))*dreal(powtemp(j-i))
     1           +dimag(lloc1(j))*dimag(powtemp(j-i))
     2           -ima*(dreal(lloc1(j))*dimag(powtemp(j-i))
     3           -dimag(lloc1(j))*dreal(powtemp(j-i))))*carray(j,i)*rsi5
            mbhloc2(i) = mbhloc2(i)+(dreal(mbhloc1(j))*dreal(jtemp(i+j))
     1           +dimag(mbhloc1(j))*dimag(jtemp(i+j))
     2           +ima*(dreal(mbhloc1(j))*dimag(jtemp(i+j))
     3           -dimag(mbhloc1(j))*dreal(jtemp(i+j))))*rsi*rsi7
            lloc2(i) = lloc2(i)+(dreal(mbhloc1(j))*dreal(jtemp(i+j))
     1           +dimag(mbhloc1(j))*dimag(jtemp(i+j))
     2           +ima*(dreal(mbhloc1(j))*dimag(jtemp(i+j))
     3           -dimag(mbhloc1(j))*dreal(jtemp(i+j))))
     4           *rsi*rsi7*dfac2(i)
            
         enddo
         rsi=rsi*rscale1
         rsi7=rsi7*rscale2
         rsi5=rsi5*rscale2/rscale1
      enddo
      return
      end

      subroutine mbh2dlocloc_pre(beta,rscale1,center1,nterms1,
     1     rscale2,center2,nterms2,dfac,dfac2,carray,mcarray,
     2     ival,diffs,pow,pow2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT: 
c
c     rscale1           : scaling for original expansion
c     center1           : center for original expansion
c     coeffs1_1(nterms1): original coeffs, difference-type functions
c     coeffs2_1(nterms1): original coeffs, modified Bessel functions
c     nterms1           : number of terms in original expansion
c     rscale2           : scaling for new expansion
c     center2           : center for new expansion
c     rad2              : radius of circle to sample on 
c                         for new expansion
c     nterms2           : number of terms in original expansion
c     zcirc2(2,nterms2) : precomp'd evenly-spaced points on unit circle
c     beta              : the modified biharmonic parameter

c     work(*)           : work array. recommended length 16*nterms+500
c
c     OUTPUT:
c
c     coeffs1_2(nterms2): new coeffs, difference-type functions
c     coeffs2_2(nterms2): new coeffs, modified Bessel functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 rscale1, center1(2), beta, rscale2, center2(2)
      integer nterms1, nterms2, mcarray
      real *8 ival(0:1), diffs(0:1), carray(0:mcarray,0:mcarray)
      real *8 dfac(0:1), dfac2(0:1), pow(0:1), pow2(0:1)
c     local
      complex *16 z,ima, zmul,zinv,ztemp1,ztemp2
      real *8 zdiff(2), r, theta, pi, done
      real *8 rsj, rsi, rsi7, rsi5, fs2
      real *8 ders(0:200), dpow(0:200)
      integer i, j, m, l, nterms, ifder
c     
      data ima/(0.0d0,1.0d0)/
c     
      done=1
      pi=4*atan(done)
c     
      nterms = nterms1+nterms2
c     
      zdiff(1)=center2(1)-center1(1)
      zdiff(2)=center2(2)-center1(2)
      call y2cart2polar(zdiff,r,theta)
      theta=theta-pi

      ifder = 0
      call diffszkik_fast(r,beta,rscale1,diffs,ifder,ders,ival,
     1     nterms)

      call mbh2d_rk(pow,dpow,r,beta,rscale1,nterms)
      call mbh2d_rksc(pow2,dpow,r,beta,rscale1,nterms)

      dfac2(0) = 1.0d0
      dfac(0) = 1.0d0

      do i = 1,nterms
         dfac2(i) = dfac2(i-1)/(2.0d0*i)
         dfac(i) = dfac(i-1)*i
      enddo

      do l = 0,nterms
         carray(l,0) = 1.0d0
      enddo
      do m=1,nterms
         carray(m,m) = 1.0d0
         do l=m+1,nterms
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
         enddo
      enddo

      return
      end

      subroutine mbh2dlocloc_wpre(beta,rscale1,center1,mbhloc1,lloc1,
     1     nterms1,rscale2,center2,mbhloc2,lloc2,nterms2,dfac,dfac2,
     2     carray,mcarray,ival,diffs,pow,pow2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT: 
c
c     rscale1           : scaling for original expansion
c     center1           : center for original expansion
c     coeffs1_1(nterms1): original coeffs, difference-type functions
c     coeffs2_1(nterms1): original coeffs, modified Bessel functions
c     nterms1           : number of terms in original expansion
c     rscale2           : scaling for new expansion
c     center2           : center for new expansion
c     rad2              : radius of circle to sample on 
c                         for new expansion
c     nterms2           : number of terms in original expansion
c     zcirc2(2,nterms2) : precomp'd evenly-spaced points on unit circle
c     beta              : the modified biharmonic parameter

c     work(*)           : work array. recommended length 16*nterms+500
c
c     OUTPUT:
c
c     coeffs1_2(nterms2): new coeffs, difference-type functions
c     coeffs2_2(nterms2): new coeffs, modified Bessel functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 rscale1, center1(2), beta, rscale2, center2(2)
      real *8 ival(0:1), diffs(0:1), carray(0:mcarray,0:mcarray)
      real *8 dfac(0:1), dfac2(0:1), pow(0:1), pow2(0:1)
      complex *16 mbhloc1(0:nterms1), lloc1(0:nterms1)
      complex *16 mbhloc2(0:nterms2), lloc2(0:nterms2)
      integer nterms1, nterms2, mcarray
c     local
      complex *16 z,ima, zmul,zinv,ztemp1,ztemp2
      complex *16 jtemp(0:200), difftemp(0:200), powtemp(0:200)
      complex *16 pow2temp(0:200)
      real *8 zdiff(2), r, theta, pi, done
      real *8 rsj, rsi, rsi7, rsi5, fs2
      integer i, j, m, l, nterms, ifder
c     
      data ima/(0.0d0,1.0d0)/
c     
      done=1
      pi=4*atan(done)
c     
      nterms = nterms1+nterms2
c     
      zdiff(1)=center2(1)-center1(1)
      zdiff(2)=center2(2)-center1(2)
      call y2cart2polar(zdiff,r,theta)
      theta=theta-pi

      do i = 0,nterms2
         mbhloc2(i) = 0
         lloc2(i) = 0
      enddo
c     
      jtemp(0) = ival(0)
      powtemp(0) = pow(0)
      pow2temp(0) = pow2(0)
      difftemp(0) = diffs(0)
      zmul=-exp(ima*theta)
      ztemp1= zmul
      do j = 1,nterms
         jtemp( j) = ztemp1*ival(j)
         difftemp(j) = ztemp1*diffs(j)
         powtemp(j) = ztemp1*pow(j)
         pow2temp(j) = ztemp1*pow2(j)
         ztemp1= ztemp1*zmul
      enddo
c     
      mbhloc2(0) = mbhloc2(0) + mbhloc1(0)*jtemp(0)
      lloc2(0) = lloc2(0) + mbhloc1(0)*difftemp(0)
      lloc2(0) = lloc2(0) + lloc1(0)*powtemp(0)
      rsj=rscale1
      do j = 1,nterms1
         mbhloc2(0) = mbhloc2(0)+(dreal(mbhloc1(j))*dreal(jtemp(j))
     1        +dimag(mbhloc1(j))*dimag(jtemp(j)))
         lloc2(0) = lloc2(0)+(dreal(mbhloc1(j))*dreal(difftemp(j))
     1        +dimag(mbhloc1(j))*dimag(difftemp(j)))
         lloc2(0) = lloc2(0)+(dreal(lloc1(j))*dreal(powtemp(j))
     1        +dimag(lloc1(j))*dimag(powtemp(j)))
      enddo
c     
      rsi=rscale1
      rsi7=rscale2
      rsi5=rscale2/rscale1
      do i = 1,nterms2
         mbhloc2(i) = mbhloc2(i) 
     1        + 2.0d0*(dreal(mbhloc1(0))*dreal(jtemp(i))
     2        +ima*dreal(mbhloc1(0))*dimag(jtemp(i)))*rsi7*rsi
         lloc2(i) = lloc2(i) 
     1        + 2.0d0*(dreal(mbhloc1(0))*dreal(jtemp(i))
     2        +ima*dreal(mbhloc1(0))*dimag(jtemp(i)))*rsi7*rsi*dfac2(i)
         fs2=rsi5*rscale1**2
         if( nterms1 .le. i-1 ) fs2=fs2*rscale1**(2*(i-1-nterms1))
         do j = min(nterms1,i-1),1,-1
            mbhloc2(i) = mbhloc2(i)+(dreal(mbhloc1(j))*dreal(jtemp(i-j))
     1           -dimag(mbhloc1(j))*dimag(jtemp(i-j))
     2           +ima*(dreal(mbhloc1(j))*dimag(jtemp(i-j))
     3           +dimag(mbhloc1(j))*dreal(jtemp(i-j))))*fs2
            mbhloc2(i) = mbhloc2(i)+(dreal(mbhloc1(j))*dreal(jtemp(i+j))
     1           +dimag(mbhloc1(j))*dimag(jtemp(i+j))
     2           +ima*(dreal(mbhloc1(j))*dimag(jtemp(i+j))
     3           -dimag(mbhloc1(j))*dreal(jtemp(i+j))))*rsi*rsi7
            lloc2(i) = lloc2(i)+(dreal(mbhloc1(j))*dreal(jtemp(i-j))
     1           -dimag(mbhloc1(j))*dimag(jtemp(i-j))
     2           +ima*(dreal(mbhloc1(j))*dimag(jtemp(i-j))
     3           +dimag(mbhloc1(j))*dreal(jtemp(i-j))))*fs2*dfac2(i)
            lloc2(i) = lloc2(i)+(dreal(mbhloc1(j))*dreal(jtemp(i+j))
     1           +dimag(mbhloc1(j))*dimag(jtemp(i+j))
     2           +ima*(dreal(mbhloc1(j))*dimag(jtemp(i+j))
     3           -dimag(mbhloc1(j))*dreal(jtemp(i+j))))
     4           *rsi*rsi7*dfac2(i)
            fs2=fs2*rscale1**2
         enddo
         do j = i,nterms1
            mbhloc2(i) = mbhloc2(i)+(dreal(mbhloc1(j))*dreal(jtemp(j-i))
     1           +dimag(mbhloc1(j))*dimag(jtemp(j-i))
     2           -ima*(dreal(mbhloc1(j))*dimag(jtemp(j-i))
     3           -dimag(mbhloc1(j))*dreal(jtemp(j-i))))*rsi5
            lloc2(i) = lloc2(i)+(dreal(mbhloc1(j))*dreal(difftemp(j-i))
     1           +dimag(mbhloc1(j))*dimag(difftemp(j-i))
     2           -ima*(dreal(mbhloc1(j))*dimag(difftemp(j-i))
     3           -dimag(mbhloc1(j))*dreal(difftemp(j-i))))*rsi5*dfac2(i)
            lloc2(i) = lloc2(i)+(dreal(lloc1(j))*dreal(powtemp(j-i))
     1           +dimag(lloc1(j))*dimag(powtemp(j-i))
     2           -ima*(dreal(lloc1(j))*dimag(powtemp(j-i))
     3           -dimag(lloc1(j))*dreal(powtemp(j-i))))*carray(j,i)*rsi5
            mbhloc2(i) = mbhloc2(i)+(dreal(mbhloc1(j))*dreal(jtemp(i+j))
     1           +dimag(mbhloc1(j))*dimag(jtemp(i+j))
     2           +ima*(dreal(mbhloc1(j))*dimag(jtemp(i+j))
     3           -dimag(mbhloc1(j))*dreal(jtemp(i+j))))*rsi*rsi7
            lloc2(i) = lloc2(i)+(dreal(mbhloc1(j))*dreal(jtemp(i+j))
     1           +dimag(mbhloc1(j))*dimag(jtemp(i+j))
     2           +ima*(dreal(mbhloc1(j))*dimag(jtemp(i+j))
     3           -dimag(mbhloc1(j))*dreal(jtemp(i+j))))
     4           *rsi*rsi7*dfac2(i)
            
         enddo
         rsi=rsi*rscale1
         rsi7=rsi7*rscale2
         rsi5=rsi5*rscale2/rscale1
      enddo
      return
      end

      subroutine mbh2dlocloc_vec(beta,rscale1,center1,mbhloc1,lloc1,
     1     nterms1,mstep1,rscale2,center2,mbhloc2,lloc2,nterms2,
     1     mstep2,nexp)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT: 
c
c     rscale1           : scaling for original expansion
c     center1           : center for original expansion
c     coeffs1_1(nterms1): original coeffs, difference-type functions
c     coeffs2_1(nterms1): original coeffs, modified Bessel functions
c     nterms1           : number of terms in original expansion
c     rscale2           : scaling for new expansion
c     center2           : center for new expansion
c     rad2              : radius of circle to sample on 
c                         for new expansion
c     nterms2           : number of terms in original expansion
c     zcirc2(2,nterms2) : precomp'd evenly-spaced points on unit circle
c     beta              : the modified biharmonic parameter

c     work(*)           : work array. recommended length 16*nterms+500
c
c     OUTPUT:
c
c     coeffs1_2(nterms2): new coeffs, difference-type functions
c     coeffs2_2(nterms2): new coeffs, modified Bessel functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 rscale1, center1(2), beta, rscale2, center2(2)
      complex *16 mbhloc1(0:mstep1,nexp), lloc1(0:mstep1,nexp)
      complex *16 mbhloc2(0:mstep2,nexp), lloc2(0:mstep2,nexp)
      integer nterms1, nterms2, mstep1, mstep2, nexp
c     local
      complex *16 z,ima, zmul,zinv,ztemp1,ztemp2
      complex *16 jtemp(0:200), difftemp(0:200), powtemp(0:200)
      complex *16 pow2temp(0:200)
      real *8 zdiff(2), r, theta, pi, done
      real *8 rsj, rsi, rsi7, rsi5, fs2
      real *8 ival(0:200), diffs(0:200), ders(0:200)
      real *8 pow(0:200), dpow(0:200), pow2(0:200)
      real *8 dfac(0:200), dfac2(0:200), carray(0:200,0:200)
      integer i, j, m, l, nterms, ifder, iii
c     
      data ima/(0.0d0,1.0d0)/
c     
      done=1
      pi=4*atan(done)
c     
      nterms = nterms1+nterms2
c     
      zdiff(1)=center2(1)-center1(1)
      zdiff(2)=center2(2)-center1(2)
      call y2cart2polar(zdiff,r,theta)
      theta=theta-pi

      call diffszkik_fast(r,beta,rscale1,diffs,ifder,ders,ival,
     1     nterms)

      call mbh2d_rk(pow,dpow,r,beta,rscale1,nterms)
      call mbh2d_rksc(pow2,dpow,r,beta,rscale1,nterms)

      dfac2(0) = 1.0d0
      dfac(0) = 1.0d0

      do i = 1,nterms
         dfac2(i) = dfac2(i-1)/(2.0d0*i)
         dfac(i) = dfac(i-1)*i
      enddo

      do l = 0,nterms
         carray(l,0) = 1.0d0
      enddo
      do m=1,nterms
         carray(m,m) = 1.0d0
         do l=m+1,nterms
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
         enddo
      enddo

c     
      jtemp(0) = ival(0)
      powtemp(0) = pow(0)
      pow2temp(0) = pow2(0)
      difftemp(0) = diffs(0)
      zmul=-exp(ima*theta)
      ztemp1= zmul
      do j = 1,nterms
         jtemp( j) = ztemp1*ival(j)
         difftemp(j) = ztemp1*diffs(j)
         powtemp(j) = ztemp1*pow(j)
         pow2temp(j) = ztemp1*pow2(j)
         ztemp1= ztemp1*zmul
      enddo
c     
      do iii = 1,nexp

         do i = 0,nterms2
            mbhloc2(i,iii) = 0
            lloc2(i,iii) = 0
         enddo


         mbhloc2(0,iii) = mbhloc2(0,iii) 
     1        + dreal(mbhloc1(0,iii))*jtemp(0)
         lloc2(0,iii) = lloc2(0,iii) + dreal(mbhloc1(0,iii))*difftemp(0)
         lloc2(0,iii) = lloc2(0,iii) + dreal(lloc1(0,iii))*powtemp(0)
         rsj=rscale1
         do j = 1,nterms1
            mbhloc2(0,iii) = mbhloc2(0,iii)
     1           +(dreal(mbhloc1(j,iii))*dreal(jtemp(j))
     1           +dimag(mbhloc1(j,iii))*dimag(jtemp(j)))
            lloc2(0,iii) = lloc2(0,iii)
     1           +(dreal(mbhloc1(j,iii))*dreal(difftemp(j))
     1           +dimag(mbhloc1(j,iii))*dimag(difftemp(j)))
            lloc2(0,iii) = lloc2(0,iii)
     1           +(dreal(lloc1(j,iii))*dreal(powtemp(j))
     1           +dimag(lloc1(j,iii))*dimag(powtemp(j)))
         enddo
c     
         rsi=rscale1
         rsi7=rscale2
         rsi5=rscale2/rscale1
         do i = 1,nterms2
            mbhloc2(i,iii) = mbhloc2(i,iii) 
     1           + 2.0d0*(dreal(mbhloc1(0,iii))*dreal(jtemp(i))
     2           +ima*dreal(mbhloc1(0,iii))*dimag(jtemp(i)))*rsi7*rsi
            lloc2(i,iii) = lloc2(i,iii) 
     1           + 2.0d0*(dreal(mbhloc1(0,iii))*dreal(jtemp(i))
     2           +ima*dreal(mbhloc1(0,iii))*dimag(jtemp(i)))
     3           *rsi7*rsi*dfac2(i)
            fs2=rsi5*rscale1**2
            if( nterms1 .le. i-1 ) fs2=fs2*rscale1**(2*(i-1-nterms1))
            do j = min(nterms1,i-1),1,-1
               mbhloc2(i,iii) = mbhloc2(i,iii)
     1              +(dreal(mbhloc1(j,iii))*dreal(jtemp(i-j))
     1              -dimag(mbhloc1(j,iii))*dimag(jtemp(i-j))
     2              +ima*(dreal(mbhloc1(j,iii))*dimag(jtemp(i-j))
     3              +dimag(mbhloc1(j,iii))*dreal(jtemp(i-j))))*fs2
               mbhloc2(i,iii) = mbhloc2(i,iii)
     1              +(dreal(mbhloc1(j,iii))*dreal(jtemp(i+j))
     1              +dimag(mbhloc1(j,iii))*dimag(jtemp(i+j))
     2              +ima*(dreal(mbhloc1(j,iii))*dimag(jtemp(i+j))
     3              -dimag(mbhloc1(j,iii))*dreal(jtemp(i+j))))*rsi*rsi7
               lloc2(i,iii) = lloc2(i,iii)
     1              +(dreal(mbhloc1(j,iii))*dreal(jtemp(i-j))
     1              -dimag(mbhloc1(j,iii))*dimag(jtemp(i-j))
     2              +ima*(dreal(mbhloc1(j,iii))*dimag(jtemp(i-j))
     3              +dimag(mbhloc1(j,iii))*dreal(jtemp(i-j))))
     5              *fs2*dfac2(i)
               lloc2(i,iii) = lloc2(i,iii)
     1              +(dreal(mbhloc1(j,iii))*dreal(jtemp(i+j))
     1              +dimag(mbhloc1(j,iii))*dimag(jtemp(i+j))
     2              +ima*(dreal(mbhloc1(j,iii))*dimag(jtemp(i+j))
     3              -dimag(mbhloc1(j,iii))*dreal(jtemp(i+j))))
     4              *rsi*rsi7*dfac2(i)
               fs2=fs2*rscale1**2
            enddo
            do j = i,nterms1
               mbhloc2(i,iii) = mbhloc2(i,iii)
     1              +(dreal(mbhloc1(j,iii))*dreal(jtemp(j-i))
     1              +dimag(mbhloc1(j,iii))*dimag(jtemp(j-i))
     2              -ima*(dreal(mbhloc1(j,iii))*dimag(jtemp(j-i))
     3              -dimag(mbhloc1(j,iii))*dreal(jtemp(j-i))))*rsi5
               lloc2(i,iii) = lloc2(i,iii)
     1              +(dreal(mbhloc1(j,iii))*dreal(difftemp(j-i))
     1              +dimag(mbhloc1(j,iii))*dimag(difftemp(j-i))
     2              -ima*(dreal(mbhloc1(j,iii))*dimag(difftemp(j-i))
     3              -dimag(mbhloc1(j,iii))*dreal(difftemp(j-i))))
     5              *rsi5*dfac2(i)
               lloc2(i,iii) = lloc2(i,iii)
     1              +(dreal(lloc1(j,iii))*dreal(powtemp(j-i))
     1              +dimag(lloc1(j,iii))*dimag(powtemp(j-i))
     2              -ima*(dreal(lloc1(j,iii))*dimag(powtemp(j-i))
     3              -dimag(lloc1(j,iii))*dreal(powtemp(j-i))))
     5              *carray(j,i)*rsi5
               mbhloc2(i,iii) = mbhloc2(i,iii)
     1              +(dreal(mbhloc1(j,iii))*dreal(jtemp(i+j))
     1              +dimag(mbhloc1(j,iii))*dimag(jtemp(i+j))
     2              +ima*(dreal(mbhloc1(j,iii))*dimag(jtemp(i+j))
     3              -dimag(mbhloc1(j,iii))*dreal(jtemp(i+j))))*rsi*rsi7
               lloc2(i,iii) = lloc2(i,iii)
     1              +(dreal(mbhloc1(j,iii))*dreal(jtemp(i+j))
     1              +dimag(mbhloc1(j,iii))*dimag(jtemp(i+j))
     2              +ima*(dreal(mbhloc1(j,iii))*dimag(jtemp(i+j))
     3              -dimag(mbhloc1(j,iii))*dreal(jtemp(i+j))))
     4              *rsi*rsi7*dfac2(i)
               
            enddo
            rsi=rsi*rscale1
            rsi7=rsi7*rscale2
            rsi5=rsi5*rscale2/rscale1
         enddo
      enddo

      return
      end

      subroutine mbh2d_childpar_pre(beta,xlengthp,ntermsc,ntermsp,
     1     rscalec,rscalep,dfac,dfac2,carray,mcarray,ival,diffs,pow)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta,xlengthp,rscalec,rscalep,dfac(0:1),dfac2(0:1)
      real *8 carray(0:mcarray,0:mcarray),ival(0:1),diffs(0:1),pow(0:1)
      integer ntermsc, ntermsp, mcarray
c     local variables
      real *8 center1(2,4), center2(2)

      center2(1) = 0.0d0
      center2(2) = 0.0d0

      center1(1,1) = -xlengthp*0.25d0
      center1(2,1) = xlengthp*0.25d0
      center1(1,2) = xlengthp*0.25d0
      center1(2,2) = xlengthp*0.25d0
      center1(1,3) = xlengthp*0.25d0
      center1(2,3) = -xlengthp*0.25d0
      center1(1,4) = -xlengthp*0.25d0
      center1(2,4) = -xlengthp*0.25d0

c     get precomputed values for first child (same values as for the
c     other children)

      call mbh2dmpmp_pre(beta,rscalec,center1(1,1),ntermsc,
     1     rscalep,center2,ntermsp,dfac,dfac2,carray,mcarray,
     2     ival,diffs,pow)

      return
      end

      subroutine mbh2d_childpar(mbhmp,ymp,mbhmpc1,ympc1,mbhmpc2,ympc2,
     1     mbhmpc3,ympc3,mbhmpc4,ympc4,ntermsc,ntermsp,rscalec,
     2     rscalep,beta,xlengthp,dfac,dfac2,carray,mcarray,ival,
     3     diffs,pow)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      complex *16 mbhmp(0:1), ymp(0:1)
      complex *16 mbhmpc1(0:1), ympc1(0:1), mbhmpc2(0:1), ympc2(0:1)
      complex *16 mbhmpc3(0:1), ympc3(0:1), mbhmpc4(0:1), ympc4(0:1)
      real *8 beta,xlengthp,rscalec,rscalep,dfac(0:1),dfac2(0:1)
      real *8 carray(0:mcarray,0:mcarray),ival(0:1),diffs(0:1),pow(0:1)
      integer ntermsc, ntermsp, mcarray
c     local variables
      integer i
      real *8 center1(2,4), center2(2)
      complex *16 mbhmptemp(0:200), ymptemp(0:200)

      center2(1) = 0.0d0
      center2(2) = 0.0d0

      center1(1,1) = -xlengthp*0.25d0
      center1(2,1) = xlengthp*0.25d0
      center1(1,2) = xlengthp*0.25d0
      center1(2,2) = xlengthp*0.25d0
      center1(1,3) = xlengthp*0.25d0
      center1(2,3) = -xlengthp*0.25d0
      center1(1,4) = -xlengthp*0.25d0
      center1(2,4) = -xlengthp*0.25d0

      do i = 0,ntermsp
         mbhmp(i) = 0.0d0
         ymp(i) = 0.0d0
      enddo

      call mbh2dmpmp_wpre(beta,rscalec,center1(1,1),mbhmpc1,ympc1,
     1     ntermsc,rscalep,center2,mbhmptemp,ymptemp,ntermsp,dfac,dfac2,
     2     carray,mcarray,ival,diffs,pow)
      
      call mbh2d_addexp(mbhmptemp,mbhmp,ntermsp)
      call mbh2d_addexp(ymptemp,ymp,ntermsp)

      call mbh2dmpmp_wpre(beta,rscalec,center1(1,2),mbhmpc2,ympc2,
     1     ntermsc,rscalep,center2,mbhmptemp,ymptemp,ntermsp,dfac,dfac2,
     2     carray,mcarray,ival,diffs,pow)

      call mbh2d_addexp(mbhmptemp,mbhmp,ntermsp)
      call mbh2d_addexp(ymptemp,ymp,ntermsp)

      call mbh2dmpmp_wpre(beta,rscalec,center1(1,3),mbhmpc3,ympc3,
     1     ntermsc,rscalep,center2,mbhmptemp,ymptemp,ntermsp,dfac,dfac2,
     2     carray,mcarray,ival,diffs,pow)

      call mbh2d_addexp(mbhmptemp,mbhmp,ntermsp)
      call mbh2d_addexp(ymptemp,ymp,ntermsp)
      
      call mbh2dmpmp_wpre(beta,rscalec,center1(1,4),mbhmpc4,ympc4,
     1     ntermsc,rscalep,center2,mbhmptemp,ymptemp,ntermsp,dfac,dfac2,
     2     carray,mcarray,ival,diffs,pow)

      call mbh2d_addexp(mbhmptemp,mbhmp,ntermsp)
      call mbh2d_addexp(ymptemp,ymp,ntermsp)

      return
      end
      
      subroutine mbh2d_parchild_pre(beta,xlengthp,ntermsc,ntermsp,
     1     rscalec,rscalep,dfac,dfac2,carray,mcarray,ival,diffs,
     2     pow,pow2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta,xlengthp,rscalec,rscalep,dfac(0:1),dfac2(0:1)
      real *8 carray(0:mcarray,0:mcarray),ival(0:1),diffs(0:1),pow(0:1)
      real *8 pow2(0:1)
      integer ntermsc, ntermsp, mcarray
c     local variables
      real *8 center2(2,4), center1(2)

      center1(1) = 0.0d0
      center1(2) = 0.0d0

      center2(1,1) = -xlengthp*0.25d0
      center2(2,1) = xlengthp*0.25d0
      center2(1,2) = xlengthp*0.25d0
      center2(2,2) = xlengthp*0.25d0
      center2(1,3) = xlengthp*0.25d0
      center2(2,3) = -xlengthp*0.25d0
      center2(1,4) = -xlengthp*0.25d0
      center2(2,4) = -xlengthp*0.25d0

c     get precomputed values for first child (same values as for the
c     other children)

      call mbh2dlocloc_pre(beta,rscalep,center1,ntermsp,
     1     rscalec,center2(1,1),ntermsc,dfac,dfac2,carray,mcarray,
     2     ival,diffs,pow,pow2)

      return
      end

      subroutine mbh2d_parchild(mbhloc,lloc,mbhlocc1,llocc1,
     1     mbhlocc2,llocc2,mbhlocc3,llocc3,mbhlocc4,llocc4,ntermsc,
     2     ntermsp,rscalec,rscalep,beta,xlengthp,dfac,dfac2,carray,
     3     mcarray,ival,diffs,pow,pow2,iswitch)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      complex *16 mbhloc(0:1), lloc(0:1)
      complex *16 mbhlocc1(0:1), llocc1(0:1), mbhlocc2(0:1), llocc2(0:1)
      complex *16 mbhlocc3(0:1), llocc3(0:1), mbhlocc4(0:1), llocc4(0:1)
      real *8 beta,xlengthp,rscalec,rscalep,dfac(0:1),dfac2(0:1)
      real *8 carray(0:mcarray,0:mcarray),ival(0:1),diffs(0:1),pow(0:1)
      real *8 pow2(0:1)
      integer ntermsc, ntermsp, mcarray, iswitch
c     local variables
      integer i
      real *8 center2(2,4), center1(2)
      complex *16 mbhloctempc(0:200), lloctempc(0:200)

      if (iswitch .eq. 0) return

      center1(1) = 0.0d0
      center1(2) = 0.0d0

      center2(1,1) = -xlengthp*0.25d0
      center2(2,1) = xlengthp*0.25d0
      center2(1,2) = xlengthp*0.25d0
      center2(2,2) = xlengthp*0.25d0
      center2(1,3) = xlengthp*0.25d0
      center2(2,3) = -xlengthp*0.25d0
      center2(1,4) = -xlengthp*0.25d0
      center2(2,4) = -xlengthp*0.25d0

      do i = 0,ntermsc
         mbhlocc1(i) = 0.0d0
         llocc1(i) = 0.0d0
         mbhlocc2(i) = 0.0d0
         llocc2(i) = 0.0d0
         mbhlocc3(i) = 0.0d0
         llocc3(i) = 0.0d0
         mbhlocc4(i) = 0.0d0
         llocc4(i) = 0.0d0
      enddo

      call mbh2dlocloc_wpre(beta,rscalep,center1,mbhloc,lloc,ntermsp,
     1     rscalec,center2(1,1),mbhloctempc,lloctempc,ntermsc,dfac,
     2     dfac2,carray,mcarray,ival,diffs,pow,pow2)
      
      call mbh2d_addexp(mbhloctempc,mbhlocc1,ntermsc)
      call mbh2d_addexp(lloctempc,llocc1,ntermsc)

      call mbh2dlocloc_wpre(beta,rscalep,center1,mbhloc,lloc,ntermsp,
     1     rscalec,center2(1,2),mbhloctempc,lloctempc,ntermsc,dfac,
     2     dfac2,carray,mcarray,ival,diffs,pow,pow2)

      call mbh2d_addexp(mbhloctempc,mbhlocc2,ntermsc)
      call mbh2d_addexp(lloctempc,llocc2,ntermsc)

      call mbh2dlocloc_wpre(beta,rscalep,center1,mbhloc,lloc,ntermsp,
     1     rscalec,center2(1,3),mbhloctempc,lloctempc,ntermsc,dfac,
     2     dfac2,carray,mcarray,ival,diffs,pow,pow2)

      call mbh2d_addexp(mbhloctempc,mbhlocc3,ntermsc)
      call mbh2d_addexp(lloctempc,llocc3,ntermsc)
      
      call mbh2dlocloc_wpre(beta,rscalep,center1,mbhloc,lloc,ntermsp,
     1     rscalec,center2(1,4),mbhloctempc,lloctempc,ntermsc,dfac,
     2     dfac2,carray,mcarray,ival,diffs,pow,pow2)

      call mbh2d_addexp(mbhloctempc,mbhlocc4,ntermsc)
      call mbh2d_addexp(lloctempc,llocc4,ntermsc)

      return
      end
      
      subroutine mbh2d_parchild_add(mbhloc,lloc,mbhlocc1,llocc1,
     1     mbhlocc2,llocc2,mbhlocc3,llocc3,mbhlocc4,llocc4,ntermsc,
     2     ntermsp,rscalec,rscalep,beta,xlengthp,dfac,dfac2,carray,
     3     mcarray,ival,diffs,pow,pow2,iswitch)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      complex *16 mbhloc(0:1), lloc(0:1)
      complex *16 mbhlocc1(0:1), llocc1(0:1), mbhlocc2(0:1), llocc2(0:1)
      complex *16 mbhlocc3(0:1), llocc3(0:1), mbhlocc4(0:1), llocc4(0:1)
      real *8 beta,xlengthp,rscalec,rscalep,dfac(0:1),dfac2(0:1)
      real *8 carray(0:mcarray,0:mcarray),ival(0:1),diffs(0:1),pow(0:1)
      real *8 pow2(0:1)
      integer ntermsc, ntermsp, mcarray, iswitch
c     local variables
      real *8 center2(2,4), center1(2)
      complex *16 mbhloctempc(0:200), lloctempc(0:200)

      if (iswitch .eq. 0) return

      center1(1) = 0.0d0
      center1(2) = 0.0d0

      center2(1,1) = -xlengthp*0.25d0
      center2(2,1) = xlengthp*0.25d0
      center2(1,2) = xlengthp*0.25d0
      center2(2,2) = xlengthp*0.25d0
      center2(1,3) = xlengthp*0.25d0
      center2(2,3) = -xlengthp*0.25d0
      center2(1,4) = -xlengthp*0.25d0
      center2(2,4) = -xlengthp*0.25d0

      call mbh2dlocloc_wpre(beta,rscalep,center1,mbhloc,lloc,ntermsp,
     1     rscalec,center2(1,1),mbhloctempc,lloctempc,ntermsc,dfac,
     2     dfac2,carray,mcarray,ival,diffs,pow,pow2)
      
      call mbh2d_addexp(mbhloctempc,mbhlocc1,ntermsc)
      call mbh2d_addexp(lloctempc,llocc1,ntermsc)

      call mbh2dlocloc_wpre(beta,rscalep,center1,mbhloc,lloc,ntermsp,
     1     rscalec,center2(1,2),mbhloctempc,lloctempc,ntermsc,dfac,
     2     dfac2,carray,mcarray,ival,diffs,pow,pow2)

      call mbh2d_addexp(mbhloctempc,mbhlocc2,ntermsc)
      call mbh2d_addexp(lloctempc,llocc2,ntermsc)

      call mbh2dlocloc_wpre(beta,rscalep,center1,mbhloc,lloc,ntermsp,
     1     rscalec,center2(1,3),mbhloctempc,lloctempc,ntermsc,dfac,
     2     dfac2,carray,mcarray,ival,diffs,pow,pow2)

      call mbh2d_addexp(mbhloctempc,mbhlocc3,ntermsc)
      call mbh2d_addexp(lloctempc,llocc3,ntermsc)
      
      call mbh2dlocloc_wpre(beta,rscalep,center1,mbhloc,lloc,ntermsp,
     1     rscalec,center2(1,4),mbhloctempc,lloctempc,ntermsc,dfac,
     2     dfac2,carray,mcarray,ival,diffs,pow,pow2)

      call mbh2d_addexp(mbhloctempc,mbhlocc4,ntermsc)
      call mbh2d_addexp(lloctempc,llocc4,ntermsc)

      return
      end

      subroutine mbh2d_processbtosfar(beta,xlengthc,rscale,nterms,
     1     mbhloc,lloc,izg4,mbhmpc1,ympc1,mbhmpc2,ympc2,mbhmpc3,ympc3,
     2     mbhmpc4,ympc4,isrcc1,isrcc2,isrcc3,isrcc4,ntermsall,
     3     dfac2all,kvecall,diffsall,powall,iswitch)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      integer nterms, izg4(2), ntermsall
      integer isrcc1, isrcc2, isrcc3, isrcc4, iswitch
      complex *16 mbhloc(0:nterms), lloc(0:nterms)
      complex *16 mbhmpc1(0:nterms), ympc1(0:nterms)
      complex *16 mbhmpc2(0:nterms), ympc2(0:nterms)
      complex *16 mbhmpc3(0:nterms), ympc3(0:nterms)
      complex *16 mbhmpc4(0:nterms), ympc4(0:nterms)
      real *8 beta, xlengthc, rscale
      real *8 dfac2all(0:ntermsall,-3:3,-3:3)
      real *8 kvecall(0:ntermsall,-3:3,-3:3)
      real *8 diffsall(0:ntermsall,-3:3,-3:3)
      real *8 powall(0:ntermsall,-3:3,-3:3)
c     local variables
      integer ioff1, joff1, ioff2, joff2, ioff3, joff3, ioff4, joff4
      real *8 center1(2), center2(2)
      complex *16 mbhloctemp(0:200), lloctemp(0:200)

      if (iswitch .eq. 0) return

      ioff4 = izg4(1)
      joff4 = izg4(2)

c     offsets for other ghosts based on 4th

      ioff1 = ioff4
      joff1 = joff4-1
      ioff2 = ioff4-1
      joff2 = joff4-1
      ioff3 = ioff4-1
      joff3 = joff4

      center1(1) = 0.0d0
      center1(2) = 0.0d0

c     add contribution from each ghost child 
c     if there is any

      center2(1) = ioff1*xlengthc
      center2(2) = joff1*xlengthc
      
      if (isrcc1 .eq. 1) then
         call mbh2dmploc_wpre(beta,rscale,center1,mbhmpc1,
     1        ympc1,nterms,rscale,center2,mbhloctemp,
     2        lloctemp,nterms,dfac2all(0,ioff1,joff1),
     3        kvecall(0,ioff1,joff1),diffsall(0,ioff1,joff1),
     4        powall(0,ioff1,joff1))

         call mbh2d_addexp(mbhloctemp,mbhloc,nterms)
         call mbh2d_addexp(lloctemp,lloc,nterms)
      endif
      
      center2(1) = ioff2*xlengthc
      center2(2) = joff2*xlengthc

      if (isrcc2 .eq. 1) then
         call mbh2dmploc_wpre(beta,rscale,center1,mbhmpc2,
     1        ympc2,nterms,rscale,center2,mbhloctemp,
     2        lloctemp,nterms,dfac2all(0,ioff2,joff2),
     3        kvecall(0,ioff2,joff2),diffsall(0,ioff2,joff2),
     4        powall(0,ioff2,joff2))

         call mbh2d_addexp(mbhloctemp,mbhloc,nterms)
         call mbh2d_addexp(lloctemp,lloc,nterms)
      endif

      center2(1) = ioff3*xlengthc
      center2(2) = joff3*xlengthc
      
      if (isrcc3 .eq. 1) then
         call mbh2dmploc_wpre(beta,rscale,center1,mbhmpc3,
     1        ympc3,nterms,rscale,center2,mbhloctemp,
     2        lloctemp,nterms,dfac2all(0,ioff3,joff3),
     3        kvecall(0,ioff3,joff3),diffsall(0,ioff3,joff3),
     4        powall(0,ioff3,joff3))
         
         call mbh2d_addexp(mbhloctemp,mbhloc,nterms)
         call mbh2d_addexp(lloctemp,lloc,nterms)
      endif

      center2(1) = ioff4*xlengthc
      center2(2) = joff4*xlengthc

      if (isrcc4 .eq. 1) then
         call mbh2dmploc_wpre(beta,rscale,center1,mbhmpc4,
     1        ympc4,nterms,rscale,center2,mbhloctemp,
     2        lloctemp,nterms,dfac2all(0,ioff4,joff4),
     3        kvecall(0,ioff4,joff4),diffsall(0,ioff4,joff4),
     4        powall(0,ioff4,joff4))

         call mbh2d_addexp(mbhloctemp,mbhloc,nterms)
         call mbh2d_addexp(lloctemp,lloc,nterms)
      endif
      

      return
      end

      subroutine mbh2d_formbtos(beta,xlengthc,srcsort,
     1     ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort,
     2     ifquad,quadstrsort,quadvecsort,ifoct,octstrsort,
     3     octvecsort,nsbox,zbox,
     3     mbhmpc1,ympc1,mbhmpc2,ympc2,mbhmpc3,ympc3,mbhmpc4,ympc4,
     4     isrcc1,isrcc2,isrcc3,isrcc4,rscaleg,ntermsg)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     zbox(2) -- coordinates of the center of the given box
c     rscaleg -- scaling for boxes the size of the ghost children
c                of the given box
c     ntermsg -- number of terms used for boxes the size of the ghosts
c     xlengthc -- side length of ghost boxes
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, xlengthc, srcsort(2,*), chargesort(*), dipstrsort(*)
      real *8 dipvecsort(2,*), quadstrsort(*), quadvecsort(3,*)
      real *8 octstrsort(*), octvecsort(4,*)
      real *8 zbox(2), rscaleg
      integer nsbox, ifcharge, ifdipole, ifquad, ifoct, ntermsg
      integer isrcc1, isrcc2, isrcc3, isrcc4
      complex *16 mbhmpc1(0:ntermsg), ympc1(0:ntermsg)
      complex *16 mbhmpc2(0:ntermsg), ympc2(0:ntermsg)
      complex *16 mbhmpc3(0:ntermsg), ympc3(0:ntermsg)
      complex *16 mbhmpc4(0:ntermsg), ympc4(0:ntermsg)
c     local variables
      real *8 zg1(2), zg2(2), zg3(3), zg4(4), x, y
      complex *16 mbhmptemp(0:200), ymptemp(0:200)
      integer nstemp, ier, i

      zg1(1) = zbox(1) - xlengthc/2.0d0
      zg1(2) = zbox(2) + xlengthc/2.0d0
      zg2(1) = zbox(1) + xlengthc/2.0d0
      zg2(2) = zbox(2) + xlengthc/2.0d0
      zg3(1) = zbox(1) + xlengthc/2.0d0
      zg3(2) = zbox(2) - xlengthc/2.0d0
      zg4(1) = zbox(1) - xlengthc/2.0d0
      zg4(2) = zbox(2) - xlengthc/2.0d0

      isrcc1 = 0
      isrcc2 = 0
      isrcc3 = 0
      isrcc4 = 0

      nstemp = 1

      do i = 0,ntermsg
         mbhmpc1(i) = 0.0d0
         ympc1(i) = 0.0d0
         mbhmpc2(i) = 0.0d0
         ympc2(i) = 0.0d0
         mbhmpc3(i) = 0.0d0
         ympc3(i) = 0.0d0
         mbhmpc4(i) = 0.0d0
         ympc4(i) = 0.0d0
      enddo

      do i = 1,nsbox
         x = srcsort(1,i)
         y = srcsort(2,i)
         if ( x .lt. zbox(1) ) then
            if (y .gt. zbox(2)) then
               isrcc1 = 1
               call mbh2dformmp_all(ier,beta,rscaleg,srcsort(1,i),
     1              ifcharge,chargesort(i),ifdipole,dipstrsort(i),
     2              dipvecsort(1,i),ifquad,quadstrsort(i),
     3              quadvecsort(1,i),ifoct,octstrsort(i),
     4              octvecsort(1,i),nstemp,zg1,ntermsg,mbhmptemp,
     4              ymptemp)
               call mbh2d_addexp(mbhmptemp,mbhmpc1,ntermsg)
               call mbh2d_addexp(ymptemp,ympc1,ntermsg)
            else
               isrcc4 = 1
               call mbh2dformmp_all(ier,beta,rscaleg,srcsort(1,i),
     1              ifcharge,chargesort(i),ifdipole,dipstrsort(i),
     2              dipvecsort(1,i),ifquad,quadstrsort(i),
     3              quadvecsort(1,i),ifoct,octstrsort(i),
     4              octvecsort(1,i),nstemp,zg4,ntermsg,mbhmptemp,
     4              ymptemp)
               call mbh2d_addexp(mbhmptemp,mbhmpc4,ntermsg)
               call mbh2d_addexp(ymptemp,ympc4,ntermsg)
            endif
         else
            if (y .gt. zbox(2)) then
               isrcc2 = 1
               call mbh2dformmp_all(ier,beta,rscaleg,srcsort(1,i),
     1              ifcharge,chargesort(i),ifdipole,dipstrsort(i),
     2              dipvecsort(1,i),ifquad,quadstrsort(i),
     3              quadvecsort(1,i),ifoct,octstrsort(i),
     4              octvecsort(1,i),nstemp,zg2,ntermsg,mbhmptemp,
     4              ymptemp)
               call mbh2d_addexp(mbhmptemp,mbhmpc2,ntermsg)
               call mbh2d_addexp(ymptemp,ympc2,ntermsg)
            else
               isrcc3 = 1
               call mbh2dformmp_all(ier,beta,rscaleg,srcsort(1,i),
     1              ifcharge,chargesort(i),ifdipole,dipstrsort(i),
     2              dipvecsort(1,i),ifquad,quadstrsort(i),
     3              quadvecsort(1,i),ifoct,octstrsort(i),
     4              octvecsort(1,i),nstemp,zg3,ntermsg,mbhmptemp,
     4              ymptemp)
               call mbh2d_addexp(mbhmptemp,mbhmpc3,ntermsg)
               call mbh2d_addexp(ymptemp,ympc3,ntermsg)
            endif
         endif
      enddo

      return
      end

      subroutine mbh2d_multar(amat,m,n,x,y)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute y = A*x
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 amat(m,n), x(*), y(*)
      integer m, n
      

      return
      end

      subroutine mbh2d_ge22cp(ier,amat,x,y)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     gaussian elimination with complete pivoting for 2x2 matrices
c     
c     solves:   amat*x = y
c
c     if the matrix is singular, the least squares solution is returned
c
c     INPUT:
c
c     a(2,2)        : system matrix a
c     y(2)          : right-hand side
c
c     OUTPUT:
c
c     x(2)          : solution
c     ier           : flag, IER = 0 means success
c                     IER = 1 means a is very small in norm
c                     IER = 2 means a is nearly not invertible
c                     IER = 3 means a is not invertible
c                      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 amat(2,2), x(2), y(2)
      integer ier
c     local
      integer imax, jmax, i, j
      real *8 a, b, c, d, e, f, g, h, gamma, dtemp, a2, b2, coa
      real *8 epssmall
      parameter (epssmall = 1.0d-200)
      
      imax = 1
      jmax = 1
      a = amat(1,1)
      if (dabs(amat(2,1)) .gt. dabs(a)) then
         imax = 2
         jmax = 1
         a = amat(2,1)
      endif
      if (dabs(amat(1,2)) .gt. dabs(a)) then
         imax = 1
         jmax = 2
         a = amat(1,2)
      endif
      if (dabs(amat(2,2)) .gt. dabs(a)) then
         imax = 2
         jmax = 2
         a = amat(2,2)
      endif

c     special case, norm of amat is tiny

      if (dabs(a) .eq. 0.0d0) then
         ier = 3
         x(1) = 0.0d0
         x(2) = 0.0d0
         return
      endif

      if (dabs(a) .lt. epssmall) ier = 1

      i = mod(imax,2)+1
      j = mod(jmax,2)+1

c     grab other entries in pivoted matrix (a,b;c,d)

      b = amat(imax,j)
      c = amat(i,jmax)
      d = amat(i,j)

c     grab pivoted right hand side

      g = y(imax)
      h = y(i)

c     perform elimination

      dtemp = d-b*c/a

c     special case, nearly not invertible or not invertible
      

c     not invertible, return least squares solution
      if (dabs(dtemp) .eq. 0.0d0) then
         ier = 3
         a2 = a**2
         b2 = b**2
         coa = c/a
         gamma = (g+h*coa)/(a+c*coa)
         e = gamma*a2/(a2+b2)
         f = gamma*b/a
         x(jmax) = e
         x(j) = f
         return
      endif

      if (dabs(dtemp) .lt. epssmall) ier = 2

      h = h-g*c/a

      f = h/dtemp
      e = (g-b*f)/a

c     copy solution to correct entries

      x(jmax) = e
      x(j) = f

      return
      end


      subroutine mbh2d_rk(pow,dpow,r,beta,rscale,nterms)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     r               : radius
c     rscale          : scaling factor
c     nterms          : number of terms in the expansion
c     
c     OUTPUT:
c
c     pow(0:nterms)   : pow(i) = (r/rscale)^i
c     dpow(0:nterms)  : dpow(i) = i*(r/rscale)^(i-1)/rscale
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 pow(0:nterms), dpow(0:nterms), r, rscale, beta
      integer nterms
c     local
      real *8 dtemp1, dtemp2
      integer i

      dtemp2 = (r*beta)/rscale
      dtemp1 = 1.0d0

      pow(0) = 1.0d0
      dpow(0) = 0.0d0

      do i = 1,nterms
         dpow(i) = i*beta*dtemp1/rscale
         dtemp1 = dtemp1*dtemp2
         if (dtemp1 .gt. 1.0d250) dtemp1 = 0.0d0
         pow(i) = dtemp1
      enddo

      return
      end

      subroutine mbh2d_rksc(pow,dpow,r,beta,rscale,nterms)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     r               : radius
c     rscale          : scaling factor
c     nterms          : number of terms in the expansion
c     
c     OUTPUT:
c
c     pow(0:nterms)   : pow(i) = (r*beta/(2*rscale))^i/(i!)
c     dpow(0:nterms)  : dpow(i) = i*(r/rscale)^(i-1)/rscale
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 pow(0:nterms), dpow(0:nterms), r, rscale, beta
      integer nterms
c     local
      real *8 dtemp1, dtemp2
      integer i

      dtemp2 = (r*beta)/(2.0d0*rscale)
      dtemp1 = 1.0d0

      pow(0) = 1.0d0
      dpow(0) = 0.0d0

      do i = 1,nterms
         dpow(i) = beta*dtemp1/(2.0d0*rscale)
         dtemp1 = dtemp1*dtemp2/i
         if (dtemp1 .gt. 1.0d250) dtemp1 = 0.0d0
         pow(i) = dtemp1
      enddo

      return
      end

      subroutine mbh2d_rmk(pow,dpow,r,beta,rscale,nterms)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     r               : radius
c     rscale          : scaling factor
c     nterms          : number of terms in the expansion
c     
c     OUTPUT:
c
c     pow(0:nterms)   : pow(i) = (rscale/r)^(-i)
c     dpow(0:nterms)  : dpow(i) = -i*(rscale/r)^(-i-1)*rscale
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 pow(0:nterms), dpow(0:nterms), r, rscale, beta
      integer nterms
c     local
      real *8 dtemp1, dtemp2
      integer i

      dtemp2 = rscale/(r*beta)
      dtemp1 = dtemp2

      pow(0) = dlog(r)
      dpow(0) = 1.0d0/r

      do i = 1,nterms
         pow(i) = dtemp1
         if (dtemp1 .lt. 1.0d-200) dtemp1 = 0.0d0
         dtemp1 = dtemp1*dtemp2
         dpow(i) = -i*dtemp1*rscale/beta
      enddo

      return
      end


      subroutine l2dmpeval_v2(beta,rscale,center,mpole,nterms,ztarg,
     1     pot,ifgrad,grad,ifhess,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     rscale          : the scaling factor
c     center(2)       : expansion center
c     coeffs1(nterms) : coefficients for difference-type functions
c     coeffs2(nterms) : coefficients for modified Bessel functions
c     nterms          : number of terms in expansions
c     ztarg(2)        : coordinates of target
c     ifgrad          : flag, IFGRAD = 1 means compute gradient
c     ifhess          : flag, IFHESS = 1 means compute Hessian
c
c     OUTPUT:
c
c     pot             : value of potential at target
c     grad(2)         : value of gradient at target (if requested)
c     hess(3)         : value of Hessian at target (if requested)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      complex *16 mpole(0:nterms)
      real *8 rscale, center(2), ztarg(2), beta
      real *8 pot, grad(2), hess(3)
      integer nterms, ifgrad, ifhess
c     local
      real *8 zdiff(2), r, theta, kvec(0:200), kder(0:200)
      real *8 pow(0:200), dpow(0:200)
      complex *16 z, eye, ztemp1, ztemp2
      complex *16 mptemp2(0:200)
      real *8 dc, ds, dc2, ds2, dcx, dsx, dcy, dsy
      integer ifder, ifders, l, j, i
      data eye /(0.0d0,1.0d0)/

      zdiff(1) = ztarg(1)-center(1)
      zdiff(2) = ztarg(2)-center(2)
      call h2cart2polar(zdiff,r,theta)

c     get values of r^-k

      call mbh2d_rmk(pow,dpow,r,beta,rscale,nterms+2)

      mptemp2(0)=pow(0)
      ztemp2=exp(eye*theta)
      ztemp1=ztemp2
      do j=1,nterms+2
         mptemp2(j)=dcmplx(pow(j)*dreal(ztemp1),pow(j)*dimag(ztemp1))
         ztemp1 = ztemp1*ztemp2
      enddo

c     evaluate potential and derivatives

      pot = 0.0d0

      pot = pot + dreal(mptemp2(0))*dreal(mpole(0))

      do j = 1,nterms
         pot = pot + dreal(mptemp2(j))*dreal(mpole(j)) +
     1        dimag(mptemp2(j))*dimag(mpole(j))
      enddo


      return
      end

      

      subroutine l2dmpmp_v2(beta,
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
      implicit real *8 (a-h,o-z)
C     
C     Usage:
C     
C     Converts multipole expansion to a multipole expansion.
C     
C---------------------------------------------------------------------
C     INPUT:
C     
C     zk      = Helmholtz parameter
C     rscale1 = scaling parameter for original multipole expansion
C     center1 = center of original multiple expansion
C     hexp    = coefficients of original multiple expansion
C     nterms1 = order of original multipole expansion
C     rscale2 = scaling parameter for shifted multipole expansion
C     center2 = center of shifted multipole expansion
C     nterms2 = order of shifted multipole expansion
C     
C     OUTPUT:
C     
C     jexp    = coefficients of shifted multipole expansion
c     
      real *8 beta
      complex *16 hexp(0:nterms1),jexp(0:nterms2)
      real *8 center1(2),center2(2),zdiff(2)
      complex *16 z,ima, zmul,zinv,ztemp1,ztemp2
      real *8 pow(0:200), dpow(0:200)
      real *8 carray(0:200,0:200)
      complex *16 powtemp(0:200)
      integer l, m

      data ima/(0.0d0,1.0d0)/
c     

      done=1
      pi=4*atan(done)
c     
      nterms = nterms1+nterms2

      do l = 0,nterms
         carray(l,0) = 1.0d0
      enddo
      do m=1,nterms
         carray(m,m) = 1.0d0
         do l=m+1,nterms
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
         enddo
      enddo


c     
      zdiff(1)=center2(1)-center1(1)
      zdiff(2)=center2(2)-center1(2)
      call y2cart2polar(zdiff,r,theta)
      theta=theta-pi
      z=zk*r
      
      call mbh2d_rk(pow,dpow,r,beta,rscale1,nterms)

      do i = 0,nterms2
         jexp(i) = 0
      enddo
c     
      powtemp(0) = pow(0)
      zmul=exp(ima*theta)
      ztemp1= zmul
      do j = 1,nterms
         powtemp( j) = ztemp1*pow(j)
         ztemp1= ztemp1*zmul
      enddo
c     
      jexp(0) = jexp(0) + hexp(0)*powtemp(0)
c     
      rsi=rscale1
      rsi7=rscale2
      rsi5=rscale1/rscale2
      do i = 1,nterms2
         jexp(i) = jexp(i) - 1.0d0*(dreal(hexp(0))*dreal(powtemp(i))
     1        +ima*dreal(hexp(0))*dimag(powtemp(i)))*rsi5/i
         rsj=rscale1
         do j = 1,min(nterms1,i)
            jexp(i) = jexp(i)+(dreal(hexp(j))*dreal(powtemp(i-j))
     1           -dimag(hexp(j))*dimag(powtemp(i-j))
     2           +ima*(dreal(hexp(j))*dimag(powtemp(i-j))
     3           +dimag(hexp(j))*dreal(powtemp(i-j))))*carray(i-1,j-1)
     4           *rsi5
            rsj=rsj*rscale1
         enddo
         rsi=rsi*rscale1
         rsi7=rsi7*rscale2
         rsi5=rsi5*rscale1/rscale2
      enddo
      return
      end
c     
      subroutine l2dtaeval_v2(beta,rscale,center,loc,nterms,ztarg,
     1     pot,ifgrad,grad,ifhess,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     rscale          : the scaling factor
c     center(2)       : expansion center
c     coeffs1(nterms) : coefficients for difference-type functions
c     coeffs2(nterms) : coefficients for modified Bessel functions
c     nterms          : number of terms in expansions
c     ztarg(2)        : coordinates of target
c     ifgrad          : flag, IFGRAD = 1 means compute gradient
c     ifhess          : flag, IFHESS = 1 means compute Hessian
c
c     OUTPUT:
c
c     pot             : value of potential at target
c     grad(2)         : value of gradient at target (if requested)
c     hess(3)         : value of Hessian at target (if requested)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      complex *16 loc(0:nterms)
      real *8 rscale, center(2), ztarg(2), beta
      real *8 pot, grad(2), hess(3)
      integer nterms, ifgrad, ifhess
c     local
      real *8 zdiff(2), r, theta, kvec(0:200), kder(0:200)
      real *8 pow(0:200), dpow(0:200)
      complex *16 z, eye, ztemp1, ztemp2
      complex *16 mptemp2(0:200)
      real *8 dc, ds, dc2, ds2, dcx, dsx, dcy, dsy
      integer ifder, ifders, l, j, i
      data eye /(0.0d0,1.0d0)/

      zdiff(1) = ztarg(1)-center(1)
      zdiff(2) = ztarg(2)-center(2)
      call h2cart2polar(zdiff,r,theta)

c     get values of r^-k

      call mbh2d_rk(pow,dpow,r,beta,rscale,nterms+2)

      mptemp2(0)=pow(0)
      ztemp2=exp(eye*theta)
      ztemp1=ztemp2
      do j=1,nterms+2
         mptemp2(j)=dcmplx(pow(j)*dreal(ztemp1),pow(j)*dimag(ztemp1))
         ztemp1 = ztemp1*ztemp2
      enddo

c     evaluate potential and derivatives

      pot = 0.0d0

      pot = pot + dreal(mptemp2(0))*dreal(loc(0))

      do j = 1,nterms
         pot = pot + dreal(mptemp2(j))*dreal(loc(j)) +
     1        dimag(mptemp2(j))*dimag(loc(j))
      enddo


      return
      end

      

      subroutine l2dmploc_v2(beta,
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
      implicit real *8 (a-h,o-z)
C     
C     Usage:
C     
C     Converts multipole expansion to a multipole expansion.
C     
C---------------------------------------------------------------------
C     INPUT:
C     
C     zk      = Helmholtz parameter
C     rscale1 = scaling parameter for original multipole expansion
C     center1 = center of original multiple expansion
C     hexp    = coefficients of original multiple expansion
C     nterms1 = order of original multipole expansion
C     rscale2 = scaling parameter for shifted multipole expansion
C     center2 = center of shifted multipole expansion
C     nterms2 = order of shifted multipole expansion
C     
C     OUTPUT:
C     
C     jexp    = coefficients of shifted multipole expansion
c     
      real *8 beta
      complex *16 hexp(0:nterms1),jexp(0:nterms2)
      real *8 center1(2),center2(2),zdiff(2)
      complex *16 z,ima, zmul,zinv,ztemp1
      real *8 pow(0:200), powtempm(0:200), dpow(0:200)
      real *8 carray(0:200,0:200)
      complex *16 powtemp(0:200)
      integer l, m

      data ima/(0.0d0,1.0d0)/
c     
      do i = 0,nterms2
         jexp(i) = 0.0d0
      enddo

      done=1
      pi=4*atan(done)
c     
      nterms = nterms1+nterms2

      do l = 0,nterms
         carray(l,0) = 1.0d0
      enddo
      do m=1,nterms
         carray(m,m) = 1.0d0
         do l=m+1,nterms
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
         enddo
      enddo

      zdiff(1) = center2(1)-center1(1)
      zdiff(2) = center2(2)-center1(2)
      call h2cart2polar(zdiff,r,theta)

      theta = theta-pi

c     get values of (beta*r)^-k
      call mbh2d_rmk(pow,dpow,r,beta,rscale1,nterms)

c     form terms for shifting 

      powtemp(0) = 1.0d0
      zmul=-exp(ima*theta)
      ztemp1= zmul
      do j = 1,nterms
         powtemp(j) = ztemp1*pow(j)
         ztemp1= ztemp1*zmul
      enddo

c     shift

      jexp(0) = jexp(0) + hexp(0)*pow(0)
      do j = 1,nterms1
         jexp(0) = jexp(0)+(dreal(hexp(j))*dreal(powtemp(j))
     1        +dimag(hexp(j))*dimag(powtemp(j)))
      enddo
c     
      rsi=rscale1
      rsi2=rscale1**2
      rsi5=rscale2/rscale1
      rsi52=rsi5*rscale2/rscale1
      isign = -1
      do i = 1,nterms2
         jexp(i) = jexp(i)
     1        - (dreal(hexp(0))*dreal(powtemp(i))
     2        +ima*dreal(hexp(0))*dimag(powtemp(i)))/i
         rsj=rscale1
         rsj2=rscale1**2
         do j = 1,min(nterms1,i)
            jexp(i) = jexp(i)+(dreal(hexp(j))*dreal(powtemp(i+j))
     1           +dimag(hexp(j))*dimag(powtemp(i+j))
     2           +ima*(dreal(hexp(j))*dimag(powtemp(i+j))
     3           -dimag(hexp(j))*dreal(powtemp(i+j))))*carray(i+j-1,j-1)

            rsj=rsj*rscale1
            rsj2=rsj2*rscale1**2
         enddo
         do j = i+1,nterms1
            jexp(i) = jexp(i)+(dreal(hexp(j))*dreal(powtemp(i+j))
     1           +dimag(hexp(j))*dimag(powtemp(i+j))
     2           +ima*(dreal(hexp(j))*dimag(powtemp(i+j))
     3           -dimag(hexp(j))*dreal(powtemp(i+j))))*carray(i+j-1,j-1)

         enddo
         jexp(i)=jexp(i)*isign
         isign = -isign
         rsi=rsi*rscale1
         rsi2=rsi2*rscale1**2
      enddo
      rsi=rscale2/rscale1
      do i = 1,nterms2
         jexp(i) = jexp(i)*rsi
         rsi=rsi*rscale2/rscale1
      enddo


      return
      end
c     
      subroutine l2dlocloc_v2(beta,
     $     rscale1,center1,hexp,nterms1,
     $     rscale2,center2,jexp,nterms2)
      implicit real *8 (a-h,o-z)
C     
C     Usage:
C     
C     Converts multipole expansion to a multipole expansion.
C     
C---------------------------------------------------------------------
C     INPUT:
C     
C     zk      = Helmholtz parameter
C     rscale1 = scaling parameter for original multipole expansion
C     center1 = center of original multiple expansion
C     hexp    = coefficients of original multiple expansion
C     nterms1 = order of original multipole expansion
C     rscale2 = scaling parameter for shifted multipole expansion
C     center2 = center of shifted multipole expansion
C     nterms2 = order of shifted multipole expansion
C     
C     OUTPUT:
C     
C     jexp    = coefficients of shifted multipole expansion
c     
      real *8 beta
      complex *16 hexp(0:nterms1),jexp(0:nterms2)
      real *8 center1(2),center2(2),zdiff(2)
      complex *16 z,ima, zmul,zinv,ztemp1
      real *8 pow(0:200), powtempm(0:200), dpow(0:200)
      real *8 carray(0:200,0:200)
      complex *16 powtemp(0:200)
      integer l, m

      data ima/(0.0d0,1.0d0)/
c     
      do i = 0,nterms2
         jexp(i) = 0.0d0
      enddo

      done=1
      pi=4*atan(done)
c     
      nterms = nterms1+nterms2

      do l = 0,nterms
         carray(l,0) = 1.0d0
      enddo
      do m=1,nterms
         carray(m,m) = 1.0d0
         do l=m+1,nterms
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
         enddo
      enddo

      zdiff(1) = center2(1)-center1(1)
      zdiff(2) = center2(2)-center1(2)
      call h2cart2polar(zdiff,r,theta)

      theta = theta-pi

c     get values of (beta*r)^k
      call mbh2d_rk(pow,dpow,r,beta,rscale1,nterms)

c     form terms for shifting 

      powtemp(0) = pow(0)
      zmul=-exp(ima*theta)
      ztemp1= zmul
      do j = 1,nterms
         powtemp(j) = ztemp1*pow(j)
         ztemp1= ztemp1*zmul
      enddo

c     shift

      jexp(0) = jexp(0) + hexp(0)*powtemp(0)
      do j = 1,nterms1
         jexp(0) = jexp(0)+(dreal(hexp(j))*dreal(powtemp(j))
     1        +dimag(hexp(j))*dimag(powtemp(j)))
      enddo
c     
      rsi5=rscale2/rscale1
      do i = 1,nterms2
         do j = i,nterms1
            jexp(i) = jexp(i)+(dreal(hexp(j))*dreal(powtemp(j-i))
     1           +dimag(hexp(j))*dimag(powtemp(j-i))
     2           -ima*(dreal(hexp(j))*dimag(powtemp(j-i))
     3           -dimag(hexp(j))*dreal(powtemp(j-i))))*carray(j,i)
     4           *rsi5
         enddo
         rsi5=rsi5*rscale2/rscale1
      enddo

      return
      end
c     
c
c


      subroutine mbh2d_addexp(b,a,nterms)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     result: a(0:nterms) = a(0:nterms) + b(0:nterms)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      complex *16 b(0:nterms), a(0:nterms)
      integer nterms
      integer i

      do i = 0,nterms
         a(i) = a(i) + b(i)
      enddo

      return
      end
      
      subroutine mbh2d_circpts(z,npts)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     npts       : the number of points on circle
c
c     OUTPUT:
c
c     z(2,npts)  : equally spaced points on unit circle
c                  z(1,i) = cos((i-1)*2*pi/npts)
c                  z(2,i) = sin((i-1)*2*pi/npts)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global 
      integer npts
      real *8 z(2,*)
c     local
      real *8 pi2
      integer i

      pi2 = 8.0d0*datan(1.0d0)
      
      do i = 1,npts
         z(1,i) = dcos((i-1)*pi2/npts)
         z(2,i) = dsin((i-1)*pi2/npts)
      enddo

      return
      end


