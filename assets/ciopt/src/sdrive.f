C
C     ***********************
C     SIMPLE DRIVER FOR LBFGS
C     ***********************
C
C     Example of driver for LBFGS routine, using a
C     simple test problem. The solution point is at 
C     X=(1,...,1) and the optimal function value of 0.
C
C                          JORGE NOCEDAL
C                        *** July 1990 ***
C
      SUBROUTINE SDRIVE(X,N,ftolb,gtolb)
      include 'implic'
      include 'efunc.cmn'
      include 'dfphess.cmn'
      integer ndimo, msave, nwork,its, i, k, iatom
      real*8 fold, fret, ediff, test, gtolb, ftolb
      PARAMETER(NDIMO=2000,MSAVE=7,NWORK=NDIMO*(2*MSAVE +1)+2*MSAVE)
      DOUBLE PRECISION X(NDIMO),G(NDIMO),DIAG(NDIMO),W(NWORK)
      DOUBLE PRECISION F,EPS,XTOL,GTOL,T1,T2,STPMIN,STPMAX
      real*8 gapnorm, perpnorm
      INTEGER IPRINT(2),IFLAG,ICALL,N,M,MP,LP,J
      LOGICAL DIAGCO
C
C     The driver for LBFGS must always declare LB2 as EXTERNAL
C
      EXTERNAL LB2
      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
C
      M=5
      IPRINT(1)= 1
      IPRINT(2)= 0
C
C     We do not wish to provide the diagonal matrices Hk0, and 
C     therefore set DIAGCO to FALSE.
C
      DIAGCO= .false.
c      print *, 'Inverse Hessian Diagonal:'
c      do i=1,N
c         print *, hessini(i)
c         DIAG(i)=hessini(i)
c      enddo
      XTOL= 1.0D-16
      ICALL=0
      IFLAG=0
C
 20   CONTINUE
      Fold=F
      call dfunc(X,G,N,F,gapnorm,perpnorm)

      its=icall
      fret=F
      ediff=Fold-F
      test=sqrt(gapnorm*gapnorm+perpnorm*perpnorm)

      write(6,1081) its,fret,test
 1081 format(" Iteration ",i5," E=",f16.8," Cvg=",f16.8)
      open(46,file='iter.log',position='append')
      if (jstate .le. 0) then
         write(46,1082) its,fret,ediff,test,gapnorm,perpnorm,
     $        envals(istate)
      else
         if (kstate .le. 0) then
            if (nefunc.eq.11) then
               write(46,1082) its,fret,rmsd,ediff,test,gapnorm,
     $              perpnorm,envals(istate),envals(jstate)
            else
               write(46,1082) its,fret,ediff,test,gapnorm,perpnorm,
     $              envals(istate),envals(jstate)
            endif
         else
            write(46,1082) its,fret,ediff,test,gapnorm,perpnorm,
     $           envals(istate),envals(jstate),envals(kstate)
         endif
      endif
      close(46)
 1082 format(i5,100(1x,f15.8))
      
      if (.not. znoncart) then
         write(6,*) 'Current Geometry'
         do iatom=1,natoms
            write(6,1004) element(iatom),(x(k),k=iatom*3-2,iatom*3)
         enddo
      endif
      
 1004 format(1x,a1,3f16.10)

      if((((((gapnorm/dlambdagap).lt.gtolb .or. nefunc.eq.1 ) .and. 
     $     perpnorm.lt.gtolb ) .or.(nefunc.eq.11.and. 
     $     sqrt(gapnorm*gapnorm+perpnorm*perpnorm).lt.gtolb)).and. 
     $     (abs(ediff)/(dlambdagap+1.0d0)).lt.ftolb )) then
         print *, 'BFGS Converged'
         print *, 'Gap Direction Gradient Norm / dlambdagap = ',
     $        (gapnorm/dlambdagap)
         print *, 'Perpendicular Direction Gradient Norm = ',perpnorm
         print *, 'Gradient Convergence Threshold: ',gtolb
         print *, 'Energy Change = ',ediff
         print *, 'Convergence Threshold: ',ftolb
         goto 50
      endif
      if ( abs(ediff).lt.1.0d-10 .and.  test.ge.gtolb ) then
         print *, 'Warning: Convergence because ediff is zero.'
         print *, 'Gradient tolerance not met!'
         goto 50
      endif
        
      CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)
C      IF(IFLAG.LE.0) GO TO 50
      ICALL=ICALL + 1
C     We allow at most 300 evaluations of F and G
      IF(ICALL.GT.300) GO TO 50
      GO TO 20
  50  CONTINUE
      END
C
C     ** LAST LINE OF SIMPLE DRIVER (SDRIVE) **

