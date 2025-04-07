!!****h* CIOpt/EvalCons.f
!  NAME
!    EvalCons.f
!*******
!!****f* EvalCons.f/EvalCons
!  NAME
!    EvalCons
!  DESCRIPTION
!    Evaluate geometric constraint terms in energy
!  SYNOPSIS
!    real*8 function EvalCons(geom,icons)
!***
      real*8 function EvalCons(geom,icons)
      include 'implic'
      include 'efunc.cmn'
      integer icons,i,j
      real*8 geom(MaxDims)
      real*8 pcoord(3,4),rij(3),rkj(3),rjk(3),rkl(3)
      real*8 BigA(3),BigB(3),BigS(3)
      real*8 dsumij,dsumkj,theta,dotijkj,dotlkkj
      real*8 dsum1,dsum2,dotab,dotkjS,phi
      real*8,parameter::Pi=3.14159625

      select case (igeocons(icons))

c     Bond Length
         case (1)
            do j=1,2
               do i=1,3
                  pcoord(i,j)=geom(iacons(j,icons)*3-3+i)
               enddo
            enddo
            dsumij=0
            do i=1,3
               rij(i)=pcoord(i,1)-pcoord(i,2)
               dsumij=dsumij+rij(i)*rij(i)
            enddo
            EvalCons=sqrt(dsumij)

c     Bond Angle
         case (2)
            do j=1,3
               do i=1,3
                  pcoord(i,j)=geom(iacons(j,icons)*3-3+i)
               enddo
            enddo
            dsumij=0
            dsumkj=0
            do i=1,3
               rij(i)=pcoord(i,1)-pcoord(i,2)
               rkj(i)=pcoord(i,3)-pcoord(i,2)
               dsumij=dsumij+rij(i)*rij(i)
               dsumkj=dsumkj+rkj(i)*rkj(i)
            enddo
            do i=1,3
               rij(i)=rij(i)/sqrt(dsumij)
               rkj(i)=rkj(i)/sqrt(dsumkj)
            enddo
            theta=0
            do i=1,3
               theta=theta+rij(i)*rkj(i)
            enddo
            EvalCons=acos(theta)

c     Dihedral Angle
         case (3)
            do j=1,4
               do i=1,3
                  pcoord(i,j)=geom(iacons(j,icons)*3-3+i)
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

            BigS(1)=BigA(2)*BigB(3)-BigA(3)*BigB(2)
            BigS(2)=BigA(3)*BigB(1)-BigA(1)*BigB(3)
            BigS(3)=BigA(1)*BigB(2)-BigA(2)*BigB(1)

            dsum1=0
            dsum2=0
          
            do i=1,3
               dsum1=dsum1+BigA(i)*BigA(i)
               dsum2=dsum2+BigB(i)*BigB(i)
            enddo
            do i=1,3
               BigA(i)=BigA(i)/sqrt(dsum1)
               BigB(i)=BigB(i)/sqrt(dsum2)
            enddo
            dotab=0
            dotkjS=0
            do i=1,3
               dotab=dotab+BigA(i)*BigB(i)
               dotkjS=dotkjS+(-rjk(i))*BigS(i)
            enddo
            phi=acos(dotab)
            if (dotkjS.lt.0) phi=-phi
            EvalCons=phi
         case default
            print *, 'Unknown constraint type in EvalCons ',icons
      end select
      return
      end
