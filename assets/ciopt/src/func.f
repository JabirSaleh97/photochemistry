!!****h* CIOpt/func.f
!  NAME
!    func.f
!*******
!!****f* func.f/func
!  NAME
!    func
!  DESCRIPTION
!    Evaluate function which is to be optimized. The functions func and dfunc MUST
!    be consistent!
!  SYNOPSIS
!    function func(p)
!***
      function func(p,gappart)

      include 'implic'
      include 'efunc.cmn'

      real*8 func, p(MaxDims),gappart
      real*8 grad(MaxDims),grad2(MaxDims),grad3(MaxDims)
      real*8 en,gap,rdiff,gap2,gapij,gapjk,gapik
      real*8 dlambda,violation
      real*8 difftmp
      integer imode,icons,k,i,j
      real*8 EvalCons
      real*8,parameter::Pi=3.14159625
      external EvalCons

      if (zrestart) then
         read(37,*,END=3000) en,gap,func
         return
 3000    continue
         close(37)
         open(16,file='long.out',status='unknown')
         zrestart=.false.
      endif
      
      call writexyz(p)

      imode=0
      call runmp(p,en,gap,grad,grad2,grad3,istate,imode)

      dlambda=dlambdagap
      select case (nefunc)
         case (1)
            func=en
         case (6)
            func=envals(istate)+
     $           dlambda*((gap*gap)/(gap+alpha))
         case (7)
            func=0.5d0*(envals(istate)+envals(jstate))+
     $           dlambda*((gap*gap)/(gap+alpha))
            gappart=((gap*gap)/(gap+alpha))
         case (8)
            func=0.5d0*(envals(istate)+envals(jstate))+
     $           0.5d0*dlambda*gap*gap
            gappart=0.5d0*gap*gap
         case (9)
            func=0.5d0*(envals(istate)+envals(jstate))-
     $           dlambda*gap+0.5d0*gap*spen*gap
            gappart=-dlambda*gap+0.5d0*gap*spen*gap
         case (10)
            gapij=envals(istate)-envals(jstate)
            gapik=envals(istate)-envals(kstate)
            gapjk=envals(jstate)-envals(kstate)
            func=0.3333333333d0*(envals(istate)+envals(jstate)+
     $           envals(kstate))+.3333333333d0*dlambda*
     $           (((gapij*gapij)/(gapij+alpha))+
     $           ((gapik*gapik)/(gapik+alpha))+
     $           ((gapjk*gapjk)/(gapjk+alpha)))
            gappart=.3333333333d0*(((gapij*gapij)/(gapij+alpha))+
     $           ((gapik*gapik)/(gapik+alpha))+
     $           ((gapjk*gapjk)/(gapjk+alpha)))
         case (11)
            rmsd=0.0d0
            do i=1,natoms
               do j=1,3
                  difftmp=p(3*(i-1)+j)-rmsdgeo(3*(i-1)+j)
                  rmsd=rmsd+rmsdweights(i)*difftmp*difftmp
               enddo
            enddo
            rmsd=sqrt(rmsd)
            func=((rmsd*rmsd)/(rmsd+alpharmsd))+
     $           dlambda*((gap*gap)/(gap+alpha))
            gappart=((gap*gap)/(gap+alpha))
         case (12)
            rmsd=0.0d0
            do i=1,natoms
               do j=1,3
                  difftmp=p(3*(i-1)+j)-rmsdgeo(3*(i-1)+j)
                  rmsd=rmsd+rmsdweights(i)*difftmp*difftmp
               enddo
            enddo
            rmsd=sqrt(rmsd)
            func=0.5d0*(envals(istate)+envals(jstate))+
     $           dlambdarmsd*((rmsd*rmsd)/(rmsd+alpharmsd))+
     $           dlambda*((gap*gap)/(gap+alpha))
            gappart=((gap*gap)/(gap+alpha))
         case default
            print *, 'Unknown nefunc in func',nefunc
            stop
      end select

cbgl  Inverse Barrier Method to constrain maximum energy.  This helps
c     with electronic structure problems, especially for nefunc=11

      if (zibf) then
         difftmp=Emaxibf-envals(istate)
         if (difftmp.gt.0.0d0) then
            func=func+ribf/difftmp
            print *,'IBF correction = ',ribf/difftmp
         else
c     yes, this is suppose to be infinity
            func=1.0d0/0.0d0
         endif
      endif

ctjm  Add in geometric constraint terms
      do icons=1,ngeocons
         if (igeocons(icons) .eq. 2 .or. igeocons(icons) .eq. 3) then
            rdiff=EvalCons(p,icons)-EvalCons(pref,icons)
            if (istcallfunc .eq. 1 .and. rdiff .lt. -Pi) then
               violation=EvalCons(p,icons)+2*Pi-dgeocons(icons)
            else if (istcallfunc .eq. 1 .and. rdiff .gt. Pi) then
               violation=EvalCons(p,icons)-2*Pi-dgeocons(icons)
            else
               violation=EvalCons(p,icons)-dgeocons(icons)
            endif
         else
            violation=EvalCons(p,icons)-dgeocons(icons)
         endif 
c         if (zseqpen) then
c            func=func+dlambdag(icons)*violation*violation
c         else
            func=func-dlambdag(icons)*violation+
     $           0.5d0*speng(icons)*violation*violation
c         endif
      enddo

      open(16,file='long.out',status='unknown',position='append')
      write(16,4000) en,gap,func,violation,(p(k),k=1,ndims)
 4000 format(400f15.8)
      close(16)

      end function

!!****f* func.f/dfunc
!  NAME
!    dfunc
!  DESCRIPTION
!    Evaluate derivative of function which is to be optimized. The functions func and dfunc MUST
!    be consistent!
!  SYNOPSIS
!    function dfunc(p,xi,nvar,fval)
!***
      subroutine dfunc(p,xi,nvar,fval,gapnorm,perpnorm)
c----------------------------------------------------------------------------------------
c     Evaluate derivative of function being optimized. 
c     
c     p(1..nvar) - Variables being optimized
c     xi(1..nvar)- Gradient of optimization function
c     nvar       - Number of variables corresponding to atomic coordinates
c     fval       - function value at current point
c     gapnorm    - gapnorm is the norm of the gradient in the gap directions
c     perpnorm   - avenorm is the norm of the gradient in the perpendicular direction
c----------------------------------------------------------------------------------------
      include 'implic'
      include 'efunc.cmn'

      integer nvar,icons,i,imode,idim,j,ixtmp,jxtmp,kxtmp,lxtmp
      real*8 dtmp,dtmp1,dtmp2,dtmp3,dtmp4
      real*8 p(nvar),xi(nvar),fp,fm,fval,violation 
      real*8 gapnorm,perpnorm,xigap(nvar),xidotxigap,xigapnorm,xinorm
      real*8 en,gap,eni,enj,enk,etmp,etmp2,etmp2a,etmp2b,dlambda
      real*8 etmp2bij,etmp2bik,etmp2bjk,etmp2ij,etmp2ik,etmp2jk
      real*8 fpi,fpj,fmi,fmj,x1itest,x1jtest,gap2
      
      real*8 ptmp(MaxDims),x1i(MaxDims),x1k(MaxDims)
      real*8 x1j(MaxDims),stpnd,fval1
      real*8 x1test(MaxDims),xjunk(MaxDims),xjunk2(MaxDims)

      real*8 pcoord(3,4),rij(3),rkj(3),rijmag
      real*8 dsumij,dsumkj,theta,dotijkj,rkjmag
      real*8 dsum1,dsum2,dsum3,dotab,rjk(3),rkl(3)
      real*8 phi,BigA(3),BigB(3),BigC(3),BigS(3),dotcb
      real*8 xgc(MaxDims), dotkjS
      real*8 xii,xjj,xkk,xll
      real*8 gappart,gpp,gpm,gpval
      real*8,parameter::Pi=3.14159625
      real*8 difftmp,xideltadotxibar,xibardotxibar

      real*8 func,funcgc,EvalCons
      external func,funcgc,EvalCons

      do i=1,MaxDims
         xgc(i)=0
      enddo

      dlambda=dlambdagap

ctjm  If we are getting gradients analytically, we need to be careful that we
ctjm  will only get "pieces" of the gradient from the electronic structure code. In 
ctjm  other words, this part of the code better be consistent with what is coded in 
ctjm  the function func.
      if (zangrad) then 
c         if (ngeocons .gt. 0) then
c            print *, 'No analytic derivatives w/geometric constraints'
c            stop
c         endif
         imode=1
         if (zmultigrad) imode=2
         if (zmultigrad.and.nefunc.eq.10) imode=3
         if (nefunc .eq. 1) imode=1
         call runmp(p,en,gap,x1i,x1j,x1k,istate,imode)
         eni=envals(istate)
         if ((imode .eq. 1) .and. (nefunc .ne. 1)) then
            call runmp(p,en,gap,x1j,xjunk,xjunk2,jstate,imode)
            if ((nefunc .eq. 10)) then
               call runmp(p,en,gap,x1k,xjunk,xjunk2,kstate,imode)
            endif
         endif
         if (nefunc .ne. 1) enj=envals(jstate)
         if (nefunc .eq. 10) enk=envals(kstate)
         print *, 'In dfunc; eni,enj=',eni,enj,en
         select case (nefunc)
            case (1) 
               fval=eni
               xinorm=0.0d0
               do idim=1,nvar
                  xi(idim)=x1i(idim)
                  xinorm=xinorm+xi(idim)*xi(idim)
                  gradave(idim)=x1i(idim)
               enddo
               gapnorm=0.0d0
               perpnorm=sqrt(xinorm)
            case (6) 
               etmp=(eni-enj+alpha)*(eni-enj+alpha)
               etmp2=(eni-enj)*(eni-enj)
               etmp2a=2.0d0*(eni-enj)*(eni-enj+alpha)
               etmp2b=(etmp2a-etmp2)/etmp
               fval=eni+dlambda*(etmp2/(eni-enj+alpha))
               do idim=1,nvar
                  xi(idim)=x1i(idim)+
     $                 dlambda*etmp2b*(x1i(idim)-x1j(idim))
                  gradave(idim)=x1i(idim)
               enddo
            case (7)
               etmp=(eni-enj+alpha)*(eni-enj+alpha)
               etmp2=(eni-enj)*(eni-enj)
               etmp2a=2.0d0*(eni-enj)*(eni-enj+alpha)
               etmp2b=(etmp2a-etmp2)/etmp
               fval=0.5d0*(eni+enj)+dlambda*(etmp2/(eni-enj+alpha))
               xideltadotxibar=0.0d0
               xibardotxibar=0.0d0
c               xidotxigap=0.0d0
c               xigapnorm=0.0d0
c               xinorm=0.0d0
               do idim=1,nvar
                  xi(idim)=0.5d0*(x1i(idim)+x1j(idim))+
     $                 dlambda*etmp2b*(x1i(idim)-x1j(idim))
                  xigap(idim)=etmp2b*(x1i(idim)-x1j(idim))
                  xideltadotxibar=xideltadotxibar+.5d0*
     $                 (x1i(idim)-x1j(idim))*(x1i(idim)+x1j(idim))
                  xibardotxibar=xibardotxibar+.25d0*
     $                 (x1i(idim)+x1j(idim))*(x1i(idim)+x1j(idim))
c                  xidotxigap=xidotxigap+xi(idim)*xigap(idim)
c                  xigapnorm=xigapnorm+xigap(idim)*xigap(idim)
c                  xinorm=xinorm+xi(idim)*xi(idim)
                  gradave(idim)=0.5d0*(x1i(idim)+x1j(idim))
               enddo
               dZsmart=xideltadotxibar/xibardotxibar
               print *, 'dZsmart = ',dZsmart
               print *, 'gaps',(eni-enj),alpha*(sqrt(dlambdagap*dZsmart/
     $              (1+dlambdagap*dZsmart))-1.0d0)
c               gapnorm=xidotxigap*xidotxigap/xigapnorm
c               perpnorm=xinorm-gapnorm
c               gapnorm=sqrt(gapnorm)
c               perpnorm=sqrt(perpnorm)
            case (8)
               etmp2=(eni-enj)*(eni-enj)
               fval=0.5d0*(eni+enj)+0.5d0*dlambda*etmp2
c               xidotxigap=0.0d0
c               xigapnorm=0.0d0
c               xinorm=0.0d0
               do idim=1,nvar
                  xi(idim)=0.5d0*(x1i(idim)+x1j(idim))+
     $                 dlambda*(x1i(idim)-x1j(idim))*(eni-enj)
                  xigap(idim)=(x1i(idim)-x1j(idim))*(eni-enj)
c                  xidotxigap=xidotxigap+xi(idim)*xigap(idim)
c                  xigapnorm=xigapnorm+xigap(idim)*xigap(idim)
c                  xinorm=xinorm+xi(idim)*xi(idim)
                  gradave(idim)=0.5d0*(x1i(idim)+x1j(idim))
               enddo
c               gapnorm=xidotxigap*xidotxigap/xigapnorm
c               perpnorm=xinorm-gapnorm
c               gapnorm=sqrt(gapnorm)
c               perpnorm=sqrt(perpnorm)
             case (9)
               etmp2=(eni-enj)*(eni-enj)
               fval=0.5d0*(eni+enj)-dlambda*(eni-enj)+0.5d0*etmp2*spen
               print *, 'In dfunc9; dlambda,spen,etmp2=',dlambda,
     $              spen,etmp2,fval
               do idim=1,nvar
                  xi(idim)=0.5d0*(x1i(idim)+x1j(idim))-
     $                 dlambda*(x1i(idim)-x1j(idim))+
     $                 spen*(eni-enj)*(x1i(idim)-x1j(idim))
                  xigap(idim)=-dlambda*(x1i(idim)-x1j(idim))+
     $                 spen*(eni-enj)*(x1i(idim)-x1j(idim))
                  gradave(idim)=0.5d0*(x1i(idim)+x1j(idim))
               enddo
            case (10)
               etmp=(eni-enj+alpha)*(eni-enj+alpha)
               etmp2ij=(eni-enj)*(eni-enj)
               etmp2a=2.0d0*(eni-enj)*(eni-enj+alpha)
               etmp2bij=(etmp2a-etmp2ij)/etmp

               etmp=(eni-enk+alpha)*(eni-enk+alpha)
               etmp2ik=(eni-enk)*(eni-enk)
               etmp2a=2.0d0*(eni-enk)*(eni-enk+alpha)
               etmp2bik=(etmp2a-etmp2ik)/etmp

               etmp=(enj-enk+alpha)*(enj-enk+alpha)
               etmp2jk=(enj-enk)*(enj-enk)
               etmp2a=2.0d0*(enj-enk)*(enj-enk+alpha)
               etmp2bjk=(etmp2a-etmp2jk)/etmp

               fval=0.3333333333d0*(eni+enj+enk)+dlambda*.3333333333d0*
     $              ((etmp2ij/(eni-enj+alpha))+
     $              (etmp2ik/(eni-enk+alpha))+
     $              (etmp2jk/(enj-enk+alpha)))
c               xidotxigap=0.0d0
c               xigapnorm=0.0d0
c               xinorm=0.0d0
               do idim=1,nvar
                  xi(idim)=0.3333333333d0*(x1i(idim)+x1j(idim)+
     $                 x1k(idim))+dlambda*.3333333333d0*
     $                 (etmp2bij*(x1i(idim)-x1j(idim))+
     $                 etmp2bik*(x1i(idim)-x1k(idim))+
     $                 etmp2bjk*(x1j(idim)-x1k(idim)))
                  xigap(idim)=.3333333333d0*
     $                 (etmp2bij*(x1i(idim)-x1j(idim))+
     $                 etmp2bik*(x1i(idim)-x1k(idim))+
     $                 etmp2bjk*(x1j(idim)-x1k(idim)))
c                  xidotxigap=xidotxigap+xi(idim)*xigap(idim)
c                  xigapnorm=xigapnorm+xigap(idim)*xigap(idim)
c                  xinorm=xinorm+xi(idim)*xi(idim)
                  gradave(idim)=0.3333333333d0*(x1i(idim)+x1j(idim)+
     $                 x1k(idim))
               enddo
c               gapnorm=xidotxigap*xidotxigap/xigapnorm
c               perpnorm=xinorm-gapnorm
c               gapnorm=sqrt(gapnorm)
c               perpnorm=sqrt(perpnorm)
            case (11)
               rmsd=0.0d0
               do i=1,natoms
                  do j=1,3
                     difftmp=p(3*(i-1)+j)-rmsdgeo(3*(i-1)+j)
                     rmsd=rmsd+rmsdweights(i)*difftmp*difftmp
                  enddo
               enddo
               rmsd=sqrt(rmsd)
               etmp=(rmsd+2*alpharmsd)/((rmsd+alpharmsd)*
     $              (rmsd+alpharmsd))
               do i=1,natoms
                  do j=1,3
                     difftmp=p(3*(i-1)+j)-rmsdgeo(3*(i-1)+j)
                     xi(3*(i-1)+j)=etmp*rmsdweights(i)*difftmp
                     print *, 'rmsdweight, alpharmsd ', rmsdweights(i),
     $                    alpharmsd
                  enddo
               enddo

               etmp=(eni-enj+alpha)*(eni-enj+alpha)
               etmp2=(eni-enj)*(eni-enj)
               etmp2a=2.0d0*(eni-enj)*(eni-enj+alpha)
               etmp2b=(etmp2a-etmp2)/etmp
               fval=((rmsd*rmsd)/(rmsd+alpharmsd))+
     $           dlambda*((gap*gap)/(gap+alpha))
               do idim=1,nvar
                  gradave(idim)=xi(idim)
                  xi(idim)=xi(idim)+
     $                 dlambda*etmp2b*(x1i(idim)-x1j(idim))
                  xigap(idim)=etmp2b*(x1i(idim)-x1j(idim))
               enddo
            case (12)
               rmsd=0.0d0
               do i=1,natoms
                  do j=1,3
                     difftmp=p(3*(i-1)+j)-rmsdgeo(3*(i-1)+j)
                     rmsd=rmsd+rmsdweights(i)*difftmp*difftmp
                  enddo
               enddo
               rmsd=sqrt(rmsd)
               etmp=(rmsd+2*alpharmsd)/((rmsd+alpharmsd)*
     $              (rmsd+alpharmsd))
               do i=1,natoms
                  do j=1,3
                     difftmp=p(3*(i-1)+j)-rmsdgeo(3*(i-1)+j)
                     xi(3*(i-1)+j)=dlambdarmsd*
     $                    etmp*rmsdweights(i)*difftmp
                  enddo
               enddo

               etmp=(eni-enj+alpha)*(eni-enj+alpha)
               etmp2=(eni-enj)*(eni-enj)
               etmp2a=2.0d0*(eni-enj)*(eni-enj+alpha)
               etmp2b=(etmp2a-etmp2)/etmp
               fval=((rmsd*rmsd)/(rmsd+alpharmsd))+
     $           dlambda*((gap*gap)/(gap+alpha))
               do idim=1,nvar
                  gradave(idim)=xi(idim)
                  xi(idim)=xi(idim)+.5d0*(x1i(idim)+x1j(idim))+
     $                 dlambda*etmp2b*(x1i(idim)-x1j(idim))
                  xigap(idim)=etmp2b*(x1i(idim)-x1j(idim))
               enddo
            case default
               print *, 'Unknown nefunc in dfunc',nefunc
               stop
         end select
cbgl  Inverse Barrier Method to constrain maximum energy.  This helps
c     with electronic structure problems, especially for nefunc=11

         if (zibf) then
            difftmp=Emaxibf-envals(istate)
            if (difftmp.gt.0.0d0) then
               fval=fval+ribf/difftmp
               do idim=1,nvar
                  xi(idim)=xi(idim)+ribf/(difftmp*difftmp)*x1i(idim)
                  print *,'bgl',idim,gradave(idim),
     $                 dlambda*xigap(idim),
     $                 ribf/(difftmp*difftmp)*x1i(idim),xi(idim)
               enddo
            else
c     yes, this is suppose to be infinity
               fval=1.0d0/0.0d0
               do idim=1,nvar
                  gradave(idim)=0.0d0
                  xi(idim)=0.0d0
                  xigap(idim)=0.0d0
               enddo
            endif
         endif

         if (perpnorm.ne.perpnorm) perpnorm=0.0d0
         write(*,*) 'nefunc = ',nefunc
         write(*,*) 'Gap = ',eni-enj
         write(*,*) 'alpha = ',alpha
         write(*,*) 'dlambdagap = ',dlambdagap
         do idim=1,nvar
            write(*,6060) x1i(idim),x1j(idim),xigap(idim),xi(idim)
         enddo
 6060    format(6f11.7)

c     Add in geometric constraint terms to function value
         do icons=1,ngeocons
            violation=EvalCons(p,icons)-dgeocons(icons)
c            if (zseqpen) then
c               fval=fval+dlambdag(icons)*violation*violation
c            else
               fval=fval-dlambdag(icons)*violation+
     $              0.5d0*speng(icons)*violation*violation
c            endif
            print *, 'c#,v,target,lambda,sigma',
     $           icons,violation,dgeocons(icons),dlambdag(icons),
     $           speng(icons)
         enddo
         print *, 'In dfunc; fval after gcons:',fval

c        Add in geometric constraint terms to derivatives
         do icons=1,ngeocons
            select case (igeocons(icons))

c	    Bond Length
            case (1)
               do j=1,2
                  do i=1,3
                     pcoord(i,j)=p(iacons(j,icons)*3-3+i)
                  enddo
               enddo
               dsumij=0
               do i=1,3
                  rij(i)=pcoord(i,1)-pcoord(i,2)
                  dsumij=dsumij+rij(i)*rij(i)
               enddo
               rijmag=sqrt(dsumij)

c              Calculate the gradient
               do i=1,3
                  dtmp=rij(i)/rijmag
                  ixtmp=(iacons(1,icons)-1)*3+i
                  jxtmp=(iacons(2,icons)-1)*3+i
                  xii=xi(ixtmp)
                  xjj=xi(jxtmp)
                  xi(ixtmp)=xi(ixtmp)-dlambdag(icons)*dtmp
                  xi(jxtmp)=xi(jxtmp)+dlambdag(icons)*dtmp
                  xi(ixtmp)=xi(ixtmp)+speng(icons)*
     $                     (rijmag-dgeocons(icons))*dtmp
                  xi(jxtmp)=xi(jxtmp)-speng(icons)*
     $                      (rijmag-dgeocons(icons))*dtmp
                  xgc(ixtmp)=xgc(ixtmp)+xi(ixtmp)-xii
                  xgc(jxtmp)=xgc(jxtmp)+xi(jxtmp)-xjj
c                  print *, 'GeoCons deriv; ',xi(ixtmp)-xii,
c     $                 xi(jxtmp)-xjj,xii,xjj
c                  print *, 'Check; ',dtmp, dlambdag(icons),
c     $              speng(icons),rijmag-dgeocons(icons)
               enddo

c	    Bond Angle
            case (2)
               do j=1,3
                  do i=1,3
                     pcoord(i,j)=p(iacons(j,icons)*3-3+i)
                  enddo
               enddo
               dsumij=0
               dsumkj=0
               dotijkj=0
               theta=0
               do i=1,3
                  rij(i)=pcoord(i,1)-pcoord(i,2)
                  rkj(i)=pcoord(i,3)-pcoord(i,2)
                  dsumij=dsumij+rij(i)*rij(i)
                  dsumkj=dsumkj+rkj(i)*rkj(i)
                  dotijkj=dotijkj+rij(i)*rkj(i)
               enddo
               rijmag=sqrt(dsumij)
               rkjmag=sqrt(dsumkj)
               theta=acos(dotijkj/(rijmag*rkjmag))

c              Calculate the gradient
               do i=1,3
                  dtmp1=((rij(i)*cos(theta)/rijmag)-(rkj(i)/rkjmag))
     $              /(sin(theta)*rijmag)
                  dtmp2=((rkj(i)*cos(theta)/rkjmag)-(rij(i)/rijmag))
     $              /(sin(theta)*rkjmag)
                  ixtmp=iacons(1,icons)*3-3+i
                  jxtmp=iacons(2,icons)*3-3+i
                  kxtmp=iacons(3,icons)*3-3+i
                  xii=xi(ixtmp)
                  xjj=xi(jxtmp)
                  xkk=xi(kxtmp)
                  xi(ixtmp)=xi(ixtmp)-dlambdag(icons)*dtmp1+
     $              speng(icons)*(theta-dgeocons(icons))*dtmp1
                  xi(jxtmp)=xi(jxtmp)+dlambdag(icons)*(dtmp1+dtmp2)-
     $              speng(icons)*(theta-dgeocons(icons))*(dtmp1+dtmp2)
                  xi(kxtmp)=xi(kxtmp)-dlambdag(icons)*dtmp2+
     $              speng(icons)*(theta-dgeocons(icons))*dtmp2
                  xgc(ixtmp)=xgc(ixtmp)+xi(ixtmp)-xii
                  xgc(jxtmp)=xgc(jxtmp)+xi(jxtmp)-xjj
                  xgc(kxtmp)=xgc(kxtmp)+xi(kxtmp)-xkk
c                  print *, 'GeoCons deriv; ',xi(ixtmp)-xii,
c     $              xi(jxtmp)-xjj,xi(kxtmp)-xkk,xii,xjj,xkk
c                  print *, 'Check; ',dtmp1,dtmp2,dlambdag(icons),
c     $              speng(icons),theta-dgeocons(icons)
               enddo

c           Dihedral Angle
            case (3)
               do j=1,4
                  do i=1,3
                     pcoord(i,j)=p(iacons(j,icons)*3-3+i)
                  enddo
               enddo
               do i=1,3
                  rij(i)=pcoord(i,1)-pcoord(i,2)
                  rjk(i)=pcoord(i,2)-pcoord(i,3)
                  rkl(i)=pcoord(i,3)-pcoord(i,4)
               enddo
               BigA(1)=rij(2)*rjk(3)-rij(3)*rjk(2)
               BigA(2)=rij(3)*rjk(1)-rij(1)*rjk(3)
               BigA(3)=rij(1)*rjk(2)-rij(2)*rjk(1)

               BigB(1)=rjk(2)*rkl(3)-rjk(3)*rkl(2)
               BigB(2)=rjk(3)*rkl(1)-rjk(1)*rkl(3)
               BigB(3)=rjk(1)*rkl(2)-rjk(2)*rkl(1)

               BigC(1)=rjk(2)*BigA(3)-rjk(3)*BigA(2)
               BigC(2)=rjk(3)*BigA(1)-rjk(1)*BigA(3)
               BigC(3)=rjk(1)*BigA(2)-rjk(2)*BigA(1)

               BigS(1)=BigA(2)*BigB(3)-BigA(3)*BigB(2)
               BigS(2)=BigA(3)*BigB(1)-BigA(1)*BigB(3)
               BigS(3)=BigA(1)*BigB(2)-BigA(2)*BigB(1)

               dsum1=0
               dsum2=0
               dsum3=0
               do i=1,3
                  dsum1=dsum1+BigA(i)*BigA(i)
                  dsum2=dsum2+BigB(i)*BigB(i)
                  dsum3=dsum3+BigC(i)*BigC(i)
               enddo
               do i=1,3
                  BigA(i)=BigA(i)/sqrt(dsum1)
                  BigB(i)=BigB(i)/sqrt(dsum2)
                  BigC(i)=BigC(i)/sqrt(dsum3)
               enddo
               dotab=0
	       dotcb=0
               dotkjS=0
               do i=1,3
                  dotab=dotab+BigA(i)*BigB(i)
                  dotcb=dotcb+BigC(i)*BigB(i)
                  dotkjS=dotkjS+(-rjk(i))*BigS(i)
               enddo
               phi=acos(dotab)
               if (dotkjS.lt.0) phi=-phi

c              Calculate the gradient
               do i=1,3
                  ixtmp=iacons(1,icons)*3-3+i
                  jxtmp=iacons(2,icons)*3-3+i
                  kxtmp=iacons(3,icons)*3-3+i
                  lxtmp=iacons(4,icons)*3-3+i
                  xii=xi(ixtmp)
                  xjj=xi(jxtmp)
                  xkk=xi(kxtmp)
                  xll=xi(lxtmp)
                  select case (i)
                     case (1)
                        if(abs(sin(phi))<=0.0001) then
                          print *,'Singularity Condition:'
                          dtmp1=((rjk(2)*rjk(2)+rjk(3)*rjk(3))*(BigB(1)
     $                      -BigC(1)*sin(phi))-rjk(1)*rjk(2)*(BigB(2)
     $                      -BigC(2)*sin(phi))-rjk(1)*rjk(3)*(BigB(3)
     $                      -BigC(3)*sin(phi)))/(cos(phi)*sqrt(dsum3))
                          dtmp2=((-(rjk(2)*rij(2)+rjk(3)*rij(3))*
     $                      (BigB(1)-BigC(1)*sin(phi))+(2*rjk(1)*rij(2)
     $                      -rjk(2)*rij(1))*(BigB(2)-BigC(2)*sin(phi))
     $                      +(2*rjk(1)*rij(3)-rjk(3)*rij(1))*(BigB(3)
     $                      -BigC(3)*sin(phi)))/(cos(phi)*sqrt(dsum3)))
     $                      +((rkl(2)*(BigC(3)-BigB(3)*sin(phi))-
     $                      rkl(3)*(BigC(2)-BigB(2)*sin(phi)))/
     $                      (cos(phi)*sqrt(dsum2)))-dtmp1
                          dtmp3=(rjk(3)*(BigC(2)-BigB(2)*sin(phi))-
     $                      rjk(2)*(BigC(3)-BigB(3)*sin(phi)))/
     $                      (-cos(phi)*sqrt(dsum2))
                          dtmp4=((-(rjk(2)*rij(2)+rjk(3)*rij(3))*
     $                      (BigB(1)-BigC(1)*sin(phi))+(2*rjk(1)*rij(2)
     $                      -rjk(2)*rij(1))*(BigB(2)-BigC(2)*sin(phi))
     $                      +(2*rjk(1)*rij(3)-rjk(3)*rij(1))*(BigB(3)
     $                      -BigC(3)*sin(phi)))/(-cos(phi)*sqrt(dsum3)))
     $                      +((rkl(2)*(BigC(3)-BigB(3)*sin(phi))-
     $                      rkl(3)*(BigC(2)-BigB(2)*sin(phi)))/
     $                      (-cos(phi)*sqrt(dsum2)))-dtmp3

                          xi(ixtmp)=xi(ixtmp)+dlambdag(icons)*dtmp1-
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp1
                          xi(jxtmp)=xi(jxtmp)+dlambdag(icons)*dtmp2-
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp2
                          xi(kxtmp)=xi(kxtmp)+dlambdag(icons)*dtmp4-
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp4
                          xi(lxtmp)=xi(lxtmp)+dlambdag(icons)*dtmp3-
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp3
                        else
                          dtmp1=(rjk(2)*(BigB(3)-BigA(3)*cos(phi))
     $                      -rjk(3)*(BigB(2)-BigA(2)*cos(phi)))
     $                      /(-sin(phi)*sqrt(dsum1))
                          dtmp2=(-(rij(3)*(BigB(2)-BigA(2)*cos(phi))-
     $                      rij(2)*(BigB(3)-BigA(3)*cos(phi)))/(sin(phi)
     $                      *sqrt(dsum1)))+(-(rkl(2)*(BigA(3)-BigB(3)*
     $                      cos(phi))-rkl(3)*(BigA(2)-BigB(2)*cos(phi)))
     $                      /(sin(phi)*sqrt(dsum2)))-dtmp1
                          dtmp3=(-(rjk(2)*(BigA(3)-BigB(3)*cos(phi))-
     $                      rjk(3)*(BigA(2)-BigB(2)*cos(phi)))/(sin(phi)
     $                      *sqrt(dsum2)))
                          dtmp4=((rij(3)*(BigB(2)-BigA(2)*cos(phi))-
     $                      rij(2)*(BigB(3)-BigA(3)*cos(phi)))/(sin(phi)
     $                      *sqrt(dsum1)))+((rkl(2)*(BigA(3)-BigB(3)*
     $                      cos(phi))-rkl(3)*(BigA(2)-BigB(2)*cos(phi)))
     $                      /(sin(phi)*sqrt(dsum2)))-dtmp3

                          xi(ixtmp)=xi(ixtmp)-dlambdag(icons)*dtmp1+
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp1
                          xi(jxtmp)=xi(jxtmp)-dlambdag(icons)*dtmp2+
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp2
                          xi(kxtmp)=xi(kxtmp)-dlambdag(icons)*dtmp4+
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp4
                          xi(lxtmp)=xi(lxtmp)-dlambdag(icons)*dtmp3+
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp3
                        end if

                     case (2)
                        if(abs(sin(phi))<=0.0001) then
                          print *,'Singularity Condition:'
                          dtmp1=((rjk(3)*rjk(3)+rjk(1)*rjk(1))*(BigB(2)
     $                      -BigC(2)*sin(phi))-rjk(2)*rjk(3)*(BigB(3)
     $                      -BigC(3)*sin(phi))-rjk(1)*rjk(2)*(BigB(1)
     $                      -BigC(1)*sin(phi)))/(cos(phi)*sqrt(dsum3))
                          dtmp2=((-(rjk(3)*rij(3)+rjk(1)*rij(1))*
     $                      (BigB(2)-BigC(2)*sin(phi))+(2*rjk(2)*rij(3)
     $                      -rjk(3)*rij(2))*(BigB(3)-BigC(3)*sin(phi))
     $                      +(2*rjk(2)*rij(1)-rjk(1)*rij(2))*(BigB(1)
     $                      -BigC(1)*sin(phi)))/(cos(phi)*sqrt(dsum3)))
     $                      +((rkl(3)*(BigC(1)-BigB(1)*sin(phi))-
     $                      rkl(1)*(BigC(3)-BigB(3)*sin(phi)))/
     $                      (cos(phi)*sqrt(dsum2)))-dtmp1
                          dtmp3=(rjk(1)*(BigC(3)-BigB(3)*sin(phi))-
     $                      rjk(3)*(BigC(1)-BigB(1)*sin(phi)))/
     $                      (-cos(phi)*sqrt(dsum2))
                          dtmp4=((-(rjk(3)*rij(3)+rjk(1)*rij(1))*
     $                      (BigB(2)-BigC(2)*sin(phi))+(2*rjk(2)*rij(3)
     $                      -rjk(3)*rij(2))*(BigB(3)-BigC(3)*sin(phi))
     $                      +(2*rjk(2)*rij(1)-rjk(1)*rij(2))*(BigB(1)
     $                      -BigC(1)*sin(phi)))/(-cos(phi)*sqrt(dsum3)))
     $                      +((rkl(3)*(BigC(1)-BigB(1)*sin(phi))-
     $                      rkl(1)*(BigC(3)-BigB(3)*sin(phi)))/
     $                      (-cos(phi)*sqrt(dsum2)))-dtmp3
                                         
                          xi(ixtmp)=xi(ixtmp)+dlambdag(icons)*dtmp1-
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp1
                          xi(jxtmp)=xi(jxtmp)+dlambdag(icons)*dtmp2-
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp2
                          xi(kxtmp)=xi(kxtmp)+dlambdag(icons)*dtmp4-
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp4
                          xi(lxtmp)=xi(lxtmp)+dlambdag(icons)*dtmp3-
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp3
                        else
                          dtmp1=(rjk(3)*(BigB(1)-BigA(1)*cos(phi))
     $                      -rjk(1)*(BigB(3)-BigA(3)*cos(phi)))
     $                      /(-sin(phi)*sqrt(dsum1))
                          dtmp2=(-(rij(1)*(BigB(3)-BigA(3)*cos(phi))-
     $                      rij(3)*(BigB(1)-BigA(1)*cos(phi)))/(sin(phi)
     $                      *sqrt(dsum1)))+(-(rkl(3)*(BigA(1)-BigB(1)*
     $                      cos(phi))-rkl(1)*(BigA(3)-BigB(3)*cos(phi)))
     $                      /(sin(phi)*sqrt(dsum2)))-dtmp1
                          dtmp3=(-(rjk(3)*(BigA(1)-BigB(1)*cos(phi))-
     $                      rjk(1)*(BigA(3)-BigB(3)*cos(phi)))/(sin(phi)
     $                      *sqrt(dsum2)))
                          dtmp4=((rij(1)*(BigB(3)-BigA(3)*cos(phi))-
     $                      rij(3)*(BigB(1)-BigA(1)*cos(phi)))/(sin(phi)
     $                      *sqrt(dsum1)))+((rkl(3)*(BigA(1)-BigB(1)*
     $                      cos(phi))-rkl(1)*(BigA(3)-BigB(3)*cos(phi)))
     $                      /(sin(phi)*sqrt(dsum2)))-dtmp3

                          xi(ixtmp)=xi(ixtmp)-dlambdag(icons)*dtmp1+
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp1
                          xi(jxtmp)=xi(jxtmp)-dlambdag(icons)*dtmp2+
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp2
                          xi(kxtmp)=xi(kxtmp)-dlambdag(icons)*dtmp4+
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp4
                          xi(lxtmp)=xi(lxtmp)-dlambdag(icons)*dtmp3+
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp3
                        end if

                     case (3)
                        if(abs(sin(phi))<=0.0001) then
                          print *,'Singularity Condition:'
                          dtmp1=((rjk(1)*rjk(1)+rjk(2)*rjk(2))*(BigB(3)
     $                      -BigC(3)*sin(phi))-rjk(3)*rjk(1)*(BigB(1)
     $                      -BigC(1)*sin(phi))-rjk(3)*rjk(2)*(BigB(2)
     $                      -BigC(2)*sin(phi)))/(cos(phi)*sqrt(dsum3))
                          dtmp2=((-(rjk(1)*rij(1)+rjk(2)*rij(2))*
     $                      (BigB(3)-BigC(3)*sin(phi))+(2*rjk(3)*rij(1)
     $                      -rjk(1)*rij(3))*(BigB(1)-BigC(1)*sin(phi))
     $                      +(2*rjk(3)*rij(2)-rjk(2)*rij(3))*(BigB(2)
     $                      -BigC(2)*sin(phi)))/(cos(phi)*sqrt(dsum3)))
     $                      +((rkl(1)*(BigC(2)-BigB(2)*sin(phi))-
     $                      rkl(2)*(BigC(1)-BigB(1)*sin(phi)))/
     $                      (cos(phi)*sqrt(dsum2)))-dtmp1
                          dtmp3=(rjk(2)*(BigC(1)-BigB(1)*sin(phi))-
     $                      rjk(1)*(BigC(2)-BigB(2)*sin(phi)))/
     $                      (-cos(phi)*sqrt(dsum2))
                          dtmp4=((-(rjk(1)*rij(1)+rjk(2)*rij(2))*
     $                      (BigB(3)-BigC(3)*sin(phi))+(2*rjk(3)*rij(1)
     $                      -rjk(1)*rij(3))*(BigB(1)-BigC(1)*sin(phi))
     $                      +(2*rjk(3)*rij(2)-rjk(2)*rij(3))*(BigB(2)
     $                      -BigC(2)*sin(phi)))/(-cos(phi)*sqrt(dsum3)))
     $                      +((rkl(1)*(BigC(2)-BigB(2)*sin(phi))-
     $                      rkl(2)*(BigC(1)-BigB(1)*sin(phi)))/
     $                      (-cos(phi)*sqrt(dsum2)))-dtmp3
        
                          xi(ixtmp)=xi(ixtmp)+dlambdag(icons)*dtmp1-
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp1
                          xi(jxtmp)=xi(jxtmp)+dlambdag(icons)*dtmp2-
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp2
                          xi(kxtmp)=xi(kxtmp)+dlambdag(icons)*dtmp4-
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp4
                          xi(lxtmp)=xi(lxtmp)+dlambdag(icons)*dtmp3-
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp3
                        else
                          dtmp1=(rjk(1)*(BigB(2)-BigA(2)*cos(phi))
     $                      -rjk(2)*(BigB(1)-BigA(1)*cos(phi)))
     $                      /(-sin(phi)*sqrt(dsum1))
                          dtmp2=(-(rij(2)*(BigB(1)-BigA(1)*cos(phi))-
     $                      rij(1)*(BigB(2)-BigA(2)*cos(phi)))/(sin(phi)
     $                      *sqrt(dsum1)))+(-(rkl(1)*(BigA(2)-BigB(2)*
     $                      cos(phi))-rkl(2)*(BigA(1)-BigB(1)*cos(phi)))
     $                      /(sin(phi)*sqrt(dsum2)))-dtmp1
                          dtmp3=(-(rjk(1)*(BigA(2)-BigB(2)*cos(phi))-
     $                      rjk(2)*(BigA(1)-BigB(1)*cos(phi)))/(sin(phi)
     $                      *sqrt(dsum2)))
                          dtmp4=((rij(2)*(BigB(1)-BigA(1)*cos(phi))-
     $                      rij(1)*(BigB(2)-BigA(2)*cos(phi)))/(sin(phi)
     $                      *sqrt(dsum1)))+((rkl(1)*(BigA(2)-BigB(2)*
     $                      cos(phi))-rkl(2)*(BigA(1)-BigB(1)*cos(phi)))
     $                      /(sin(phi)*sqrt(dsum2)))-dtmp3

                          xi(ixtmp)=xi(ixtmp)-dlambdag(icons)*dtmp1+
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp1
                          xi(jxtmp)=xi(jxtmp)-dlambdag(icons)*dtmp2+
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp2
                          xi(kxtmp)=xi(kxtmp)-dlambdag(icons)*dtmp4+
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp4
                          xi(lxtmp)=xi(lxtmp)-dlambdag(icons)*dtmp3+
     $                      speng(icons)*(phi-dgeocons(icons))*dtmp3
                        end if
                     case default
                        print *,'Unknown case in dihedral gradient'
                     end select
                     xgc(ixtmp)=xgc(ixtmp)+xi(ixtmp)-xii
                     xgc(jxtmp)=xgc(jxtmp)+xi(jxtmp)-xjj
                     xgc(kxtmp)=xgc(kxtmp)+xi(kxtmp)-xkk
                     xgc(lxtmp)=xgc(lxtmp)+xi(lxtmp)-xll
c                     print *, 'GeoCons deriv; ',xi(ixtmp)-xii,
c     $                 xi(jxtmp)-xjj,xi(kxtmp)-xkk,xi(lxtmp)-xll,
c     $                 xii,xjj,xkk,xll
c                     print *,'Check; ',dtmp1,dtmp2,dtmp3,dtmp4,phi,
c     $                 dlambdag(icons),speng(icons),dotkjS
               end do 
            case default
               print *, 'Unknown constraint type in dfunc ', icons
            end select
         enddo
         xidotxigap=0.0d0
         xigapnorm=0.0d0
         xinorm=0.0d0
         do idim=1,nvar
            xidotxigap=xidotxigap+xi(idim)*xigap(idim)
            xigapnorm=xigapnorm+xigap(idim)*xigap(idim)
            xinorm=xinorm+xi(idim)*xi(idim)
         enddo
         gapnorm=xidotxigap*xidotxigap/xigapnorm
         perpnorm=xinorm-gapnorm
         gapnorm=sqrt(gapnorm)
         perpnorm=sqrt(perpnorm)

ctjm Now compute finite difference gradient and check
         if (ngeocons .ne. 0) then
            print *, 'Check Constraint Terms in gradient'
            do i=1,nvar
               ptmp(i)=p(i)
               pref(i)=p(i)
            enddo         
            fval1=funcgc(ptmp)
            print *, 'Function value: ',fval,fval1
            print *, 'State en: ',eni,enj
            do i=1,nvar
               stpnd=scalend(i)
               ptmp(i)=p(i)+stpnd
               fp=funcgc(ptmp)
               ptmp(i)=p(i)-stpnd
               fm=funcgc(ptmp)
               ptmp(i)=p(i)
               x1test(i)=(fp-fm)/(2.0d0*stpnd)
               print *, i,xgc(i),x1test(i)
            enddo       
         endif
         if (zcheckgrad) then
            print *, 'Check gradient'
            do i=1,nvar
               ptmp(i)=p(i)
               pref(i)=p(i)
            enddo         
            fval1=func(ptmp,gpval)
            print *, 'Function value: ',fval,fval1
            print *, 'State en: ',eni,enj
            xidotxigap=0.0d0
            xigapnorm=0.0d0
            xinorm=0.0d0
            do i=1,nvar
               stpnd=scalend(i)
               ptmp(i)=p(i)+stpnd
               fp=func(ptmp,gpp)
               fpi=envals(istate)
               fpj=envals(jstate)
               ptmp(i)=p(i)-stpnd
               fm=func(ptmp,gpm)
               fmi=envals(istate)
               fmj=envals(jstate)
               ptmp(i)=p(i)
               x1test(i)=(fp-fm)/(2.0d0*stpnd)
               xigap(i)=(gpp-gpm)/(2.0d0*stpnd)
               x1itest=(fpi-fmi)/(2.0d0*stpnd)
               x1jtest=(fpj-fmj)/(2.0d0*stpnd)
               xidotxigap=xidotxigap+x1test(i)*xigap(i)
               xigapnorm=xigapnorm+xigap(i)*xigap(i)
               xinorm=xinorm+x1test(i)*x1test(i)
               print *, i,xi(i),x1test(i)
               print *, '   IState:',x1itest,x1i(i)
               print *, '   JState:',x1jtest,x1j(i)
            enddo       
            gapnorm=xidotxigap*xidotxigap/xigapnorm
            perpnorm=xinorm-gapnorm
            gapnorm=sqrt(gapnorm)
            perpnorm=sqrt(perpnorm)
         endif
         istcallfunc=1
      else
         do i=1,nvar
            ptmp(i)=p(i)
            pref(i)=p(i)
         enddo
         fval=func(ptmp,gpval)
         xidotxigap=0.0d0
         xigapnorm=0.0d0
         xinorm=0.0d0
         if (zforward) then
            do i=1,nvar
               stpnd=scalend(i)
               ptmp(i)=p(i)+stpnd
               fp=func(ptmp,gpp)
               ptmp(i)=p(i)
               xi(i)=(fp-fval)/stpnd
               xigap(i)=(gpp-gpval)/stpnd
               xidotxigap=xidotxigap+xi(i)*xigap(i)
               xigapnorm=xigapnorm+xigap(i)*xigap(i)
               xinorm=xinorm+xi(i)*xi(i)
            enddo
            istcallfunc=1
         else
            do i=1,nvar
               stpnd=scalend(i)
               ptmp(i)=p(i)+stpnd
               fp=func(ptmp,gpp)
               ptmp(i)=p(i)-stpnd
               fm=func(ptmp,gpm)
               ptmp(i)=p(i)
               xi(i)=(fp-fm)/(2.0d0*stpnd)
               xigap(i)=(gpp-gpm)/(2.0d0*stpnd)
               xidotxigap=xidotxigap+xi(i)*xigap(i)
               xigapnorm=xigapnorm+xigap(i)*xigap(i)
               xinorm=xinorm+xi(i)*xi(i)
            enddo
            istcallfunc=1
         endif
         gapnorm=xidotxigap*xidotxigap/xigapnorm
         perpnorm=xinorm-gapnorm
         gapnorm=sqrt(gapnorm)
         perpnorm=sqrt(perpnorm)
      endif
      return
      end

      function funcgc(p)

      include 'implic'
      include 'efunc.cmn'

      real*8 funcgc, p(MaxDims)
      real*8 grad(MaxDims),grad2(MaxDims)
      real*8 en,gap,rdiff
      real*8 dlambda,violation
      integer icons
      real*8 EvalCons
      real*8,parameter::Pi=3.14159625
      external EvalCons

c      if (zrestart) then
c         read(37,*,END=3000) en,gap,func
c         return
c 3000    continue
c         close(37)
c         open(16,file='long.out',status='unknown')
c         zrestart=.false.
c      endif
      
c      call writexyz(p)

c      imode=0
c      call runmp(p,en,gap,grad,grad2,istate,imode)


      dlambda=dlambdagap
c      select case (nefunc)
c         case (1)
c            func=en
c         case (6)
c            func=envals(istate)+
c     $           dlambda*((gap*gap)/(gap+alpha))
c         case (7)
c            func=0.5d0*(envals(istate)+envals(jstate))+
c     $           dlambda*((gap*gap)/(gap+alpha))
c         case (8)
c            func=0.5d0*(envals(istate)+envals(jstate))+
c     $           0.5d0*dlambda*gap*gap
c         case (9)
c            func=0.5d0*(envals(istate)+envals(jstate))-
c     $           dlambda*gap+0.5d0*gap*spen*gap
c         case default
c            print *, 'Unknown nefunc in func',nefunc
c            stop
c      end select

ctjm  Add in geometric constraint terms
      funcgc=0
      do icons=1,ngeocons
         if (igeocons(icons) .eq. 2 .or. igeocons(icons) .eq. 3) then
            rdiff=EvalCons(p,icons)-EvalCons(pref,icons)
            if (istcallfunc .eq. 1 .and. rdiff .lt. -Pi) then
               violation=EvalCons(p,icons)+2*Pi-dgeocons(icons)
            else if(istcallfunc .eq. 1 .and. rdiff .gt. Pi) then
               violation=EvalCons(p,icons)-2*Pi-dgeocons(icons)
            else
               violation=EvalCons(p,icons)-dgeocons(icons)
            endif
         else
           violation=EvalCons(p,icons)-dgeocons(icons)
         endif
c         if (zseqpen) then
c            funcgc=funcgc+dlambdag(icons)*violation*violation
c         else
            funcgc=funcgc-dlambdag(icons)*violation+
     $           0.5d0*speng(icons)*violation*violation
c         endif
      enddo

c      open(16,file='long.out',status='unknown',position='append')
c      write(16,4000) en,gap,func,violation,(p(k),k=1,ndims)
c 4000 format(400f15.8)
c      close(16)

      end function

