      subroutine runmp(geom,en,gap,grad1,grad2,grad3,ixstate,imode)
c
c     geom  - input geometry parsed by WriteInput
c     imode=0; energy only
c     imode=1; energy+gradient for ixstate
c     imode=2; energy+gradients for istate and jstate
c     imode=3; energy+gradients for istate, jstate, and kstate
c
      include 'implic'
      include 'efunc.cmn'

      integer imode,ixstate,j
      integer iswapra(MaxStates)
      real*8 geom(MaxDims),e(20),en,gap,dtemp
      real*8 grad1(MaxDims),grad2(MaxDims),grad3(MaxDims)
      real*8 pold(MaxDims),dlambdaold,spenold
      real*8 g1old(MaxDims),g2old(MaxDims),g3old(MaxDims)
      real*8 envold(20)
      character*256 cbuffer,cbuf1,ctmpw
      real*8 ptmp(3)
      real*8 tmp
      integer i,idiff,ixold,istcall
      save istcall,dlambdaold,spenold,ixold,pold,g1old,g2old,
     $     g3old,envold

      include 'tjmfunc.cmn'

      print *, 'Begin Run MolPro... imode = ',imode
      call flush(6)

      if (istcall .eq. 0) then
         istcall=1
         do i=1,ndims
            pold(i)=geom(i)
         enddo
         dlambdaold=dlambdagap
         spenold=spen
         ixold=ixstate
      else
         idiff=0
         do i=1,ndims
            if (pold(i) .ne. geom(i)) idiff=1
         enddo
         if (dlambdaold .ne. dlambdagap) idiff=1
         if (spenold .ne. spen) idiff=1
         if (ixold .ne. ixstate) idiff=1
         if (imode .ne. 0) idiff=1
         if (idiff .eq. 0) then
            print *, 'Skip calc!!!'
            do i=1,ndims
               grad1(i)=g1old(i)
               grad2(i)=g2old(i)
               grad3(i)=g3old(i)
            enddo
            en=envold(ixstate)
            gap=envold(istate)-envold(jstate)
            return
         endif
c         print *, 'Why we calculated:'
c         print *, dlambdaold,dlambdagap
c         print *, spenold,spen
c         print *, ixold,ixstate
c         print *, imode
c         do i=1,ndims
c            print *, geom(i),pold(i)
c         enddo
c         print *, 'And thats the story'
c         print *, ' '
         do i=1,ndims
            pold(i)=geom(i)
         enddo
         spenold=spen
         dlambdaold=dlambdagap
         ixold=ixstate
      endif

      open(13,file=trim(cinpdeck),status='unknown')
      close(13,status='delete')
      if (imode .eq. 0) then
         ctmpw=ctmpwrite
      else
         ctmpw=ctmpwriteg
      endif
      call WriteInput(geom,natoms*3,ixstate,ctmpw,cinpdeck)

      nescall=nescall+1
      call system(trim(crunstr))
      if (nescall .lt. 10) then
         write(cbuf1,'(i1)') nescall
      elseif (nescall .lt. 100) then
         write(cbuf1,'(i2)') nescall
      elseif (nescall .lt. 1000) then
         write(cbuf1,'(i3)') nescall
      else
         cbuf1='xxx'
      endif
      if (zdetails) then
         cbuffer='cp '//trim(coutfile)//' details/'//
     $        trim(coutfile)//'.'//trim(cbuf1)
         call system(trim(cbuffer))
      endif

c      print *, 'Read energies...'
      call ReadOutput(ptmp,nstates,ctmpread,coutfile)
      do i=1,nstates
         e(i)=ptmp(i)
      enddo
      if (ZExEnEv) then
         do i=2,nstates
            e(i)=e(1)+e(i)/27.21d0
         enddo
      endif
      
      do i=1,nstates
         iswapra(i)=i
      enddo
      do i=1,nstates
         do j=i+1,nstates
            if (e(i).gt.e(j)) then
               print *, '*** Root flipping detected.  Swapping: ',
     $              i,j,e(i),e(j)
               iswapra(i)=j
               iswapra(j)=i
               tmp=e(i)
               e(i)=e(j)
               e(j)=tmp
            endif
         enddo
      enddo
      do i=1,nstates
         envals(i)=e(i)
         envold(i)=e(i)
      enddo
      open(65,file='mplog.out',position='APPEND')
      write(65,*) (envals(i),i=1,nstates)
      close(65)
      en=e(ixstate)
      gap=e(istate)-e(jstate)


      print *, imode
      if (zangrad .and. imode .ne. 0) then
         call ReadOutput(grad1,3*natoms,ctmpgread,coutfile)
         do i=1,ndims
            g1old(i)=grad1(i)
         enddo
         print *, 'Read gradient for istate'
         do i=1,ndims
            print *, i, grad1(i)
         enddo
         call flush(6)
         if (imode .ge. 2) then
            call ReadOutput(grad2,3*natoms,ctmpg2read,coutfile)
            do i=1,ndims
               g2old(i)=grad2(i)
            enddo
            print *, 'Read gradient for jstate, multigrad',ndims
            do i=1,ndims
               print *, i, grad1(i),grad2(i)
            enddo
            if (imode.eq.3) then
               call ReadOutput(grad3,3*natoms,ctmpg3read,coutfile)
               do i=1,ndims
                  g3old(i)=grad3(i)
               enddo
               print *, 'Read gradient for kstate, multigrad'
               do i=1,ndims
                  print *, i, grad1(i),grad2(i),grad3(i)
               enddo
            endif
         endif
c     Check if we had to swap roots - the only case where this does not cause an abort 
c     at present is when the swapping was b/t istate and jstate and we are calculating
c     both gradients at once.  In order to make the code work in cases where the 
c     swapping involves some state other than istate/jstate or when we don't have the required
c     gradient to swap, we would have to add code here to calculate the required gradient.
         if (iswapra(istate) .ne. istate .and. imode .lt. 2) then
            print *, 'Had to swap energies but did ',
     $           'not have required gradients'
            print *, 'Aborting...'
            stop
         endif
         if (iswapra(istate) .ne. istate) then
            if (iswapra(istate) .eq. jstate .and.
     $           iswapra(jstate) .eq. istate) then
               print *, 'States i and j swapped. ',
     $              'Swapping gradients...'
               print *, 'This code is not tested ',
     $              'for root swapping with analytic ',
     $              'gradients. You were warned...'
               do i=1,ndims
                  dtemp=grad1(i)
                  grad1(i)=grad2(i)
                  grad2(i)=dtemp
               enddo
            else
               print *, 'Had to swap energies but did ',
     $              'not have required gradients'
               print *, 'Aborting...'
               stop       
            endif
         endif
      endif      

c      do i=1,3*natoms
c         print *, i, grad1(i)
c      enddo
c      stop

      print *, '...Run MolPro successful'
      call flush(6)
      return
 9000 continue
      write(6,*) 'error reading electronic structure output'
      stop
      end 
