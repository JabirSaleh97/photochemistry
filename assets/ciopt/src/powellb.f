      subroutine powellb(p,xi,tol,iter,fret)

      include 'implic'
      include 'efunc.cmn'
      integer iter
      real*8 p(MaxDims),xi(MaxDims,MaxDims),tol,fret,fp,fpp,fbig
      real*8 p0(MaxDims),xii(MaxDims),estep
      integer i,j,ibig
      
      fp=0.0d0
      fret=0.0d0
      print *, "beginning powell optimization"

      do iter=1,200

         if (.not. znoncart) then
            write(6,*) 'Current Geometry'
            do i=1,natoms
               write(6,1004) element(i),(p(j),j=i*3-2,i*3)
            enddo
         endif
         
 1004    format(1x,a1,3f16.10)

         do i=1,ndims
            p0(i)=p(i)
         enddo
         
         fbig=0.0d0
         do i=1,ndims
            do j=1,ndims
               xii(j)=xi(i,j)
            enddo
            fpp=fret
            print *, "entering brentb for line search #",i
            call brentb(p,xii,ndims,tol,fret)
            if ((fpp-fret).gt.fbig) then
               fbig=(fpp-fret)
               ibig=i
            endif
         enddo
         
         do i=ibig,ndims-1
            do j=1,ndims
               xi(i,j)=xi(i+1,j)
            enddo
         enddo
            
         do j=1,ndims
            xi(ndims,j)=p(j)-p0(j)
         enddo

         do j=1,ndims
            xii(j)=xi(ndims,j)
         enddo
         print *, "entering brentb for conjugate direction line search"
         call brentb(p,xii,ndims,tol,fret)

         estep=abs(fp-fret)
         fp=fret

         open(46,file='iter.log',position='append')
         if (jstate .le. 0) then
            write(46,1082) iter,fret,estep,
     $           envals(istate)
         else
            if (kstate .le. 0) then
               write(46,1082) iter,fret,estep,
     $              envals(istate),envals(jstate)
            else
               write(46,1082) iter,fret,estep,
     $              envals(istate),envals(jstate),
     $              envals(kstate)
            endif
         endif
         close(46)
 1082    format(i5,100(1x,f15.8))
         if (estep.lt.tol) then
            print *, 'Powell Converged: DeltaE=',abs(fp-fret)
            print *, 'Convergence Threshold: ',tol
            return
         endif
      enddo
      end subroutine

      subroutine brentb(p,xii,ndims,tol,fret)
      integer ndims
      real*8 p(ndims),xii(ndims),fret,tol
      real*8 gold, a, b, x, u, v, w
      real*8 fa,fb,fx,fu,fv,fw,gappart,tmp
      real*8 xadiff,xbdiff,fxfadiff,fxfbdiff,numer,denom,xmean
      real*8 func1d
      real*8 cgold
      parameter (cgold=.381966)
      external func1d
      integer i,j,iter
      parameter(gold=1.618034d0,eps=1.0d-10)

      print *,"bracketing minimum"

c     first bracket the minimum
      a=0
      x=1
      fa=func1d(a,p,xii,ndims)
      fx=func1d(x,p,xii,ndims)
      if (fx.gt.fa) then
         tmp=a
         a=x
         x=tmp
         tmp=fa
         fa=fx
         fx=tmp
      endif

 1000 continue
      
      b=x+gold*(x-a)
      fb=func1d(b,p,xii,ndims)
      print *,"a, x, b, fa, fx, fb:"
      print *, a,x,b
      print *, fa, fx, fb      
      if (fx.gt.fb) then
         a=x
         x=b
         fa=fx
         fx=fb
         goto 1000
      endif

      if (b.lt.a) then
         tmp=a
         a=b
         b=tmp
      endif

      print *,"bracketting done"
c now that we have bracketed the minimum, do the optimization

      print *,"performing line search"

      fret=fx

      if (b.lt.a) then
         tmp=a
         a=b
         b=tmp
      endif
      
      if (fa.lt.fb) then
         v=b
         w=a
         fv=fb
         fw=fa
      else
         v=a
         w=b
         fv=fa
         fw=fb
      endif
      
      nlastdist=0
      lastdist=0
      step=0

      do iter=1,200
         print *,"a, x, b, fa, fx, fb:"
         print *, a, x, b
         print *, fa, fx, fb
         if (((max(fb,fa)-fx).lt.tol).or.((b-a).lt.eps)) then
            do i=1,ndims
               p(i)=p(i)+xii(i)*x
            enddo
            print *, "line search complete"
            return
         endif

         xmean=.5d0*(a+b)
         xvdiff=x-v
         xwdiff=x-w
         fxfvdiff=fx-fv
         fxfwdiff=fx-fw
         numer=xvdiff*xvdiff*fxfwdiff-xwdiff*xwdiff*fxfvdiff
         denom=xvdiff*fxfwdiff-xwdiff*fxfvdiff
         
         nlastdist=lastdist
         lastdist=abs(step)
         step=-.50d0*numer/denom
         u=x+step

         
         if ((denom.eq.0.0d0).or.(abs(step).gt..5*nlastdist).or.
     $        (u.ge.b).or.(u.le.a)) then
            if (x.gt.xmean) then
               step=cgold*(a-x)
            else
               step=cgold*(b-x)
            endif
            u=x+step
         endif
         
         fu=func1d(u,p,xii,ndims)
         print *, 'in line search u, fu:'
         print *, u,fu

         if (fu.le.fx) then
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
            fret=fx
         else
            if (u.lt.x) then
               a=u
               fa=fu
            elseif (u.gt.x) then
               b=u
               fb=fu
            endif
            if (fu.lt.fw) then
               v=w
               fv=fw
               w=u
               fw=fu
            endif
         endif
      enddo
      end subroutine
      
      function func1d(x,p,xii,ndims)
      real*8 func1d
      integer ndims
      real*8 x, p(ndims), xii(ndims)
      real*8 ptmp(ndims),gappart
      real*8 func
      external func
      integer i
      
      do i=1,ndims
         ptmp(i)=p(i)+x*xii(i)
      enddo
      func1d=func(ptmp,gappart)
      end function
