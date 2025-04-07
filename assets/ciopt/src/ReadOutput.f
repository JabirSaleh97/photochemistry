!!****h* CIOpt/ReadOutput.f
!  NAME
!    ReadOutput.f
!*******
!!****f* ReadOutput.f/ReadOutput
!  NAME
!    ReadOutput
!  DESCRIPTION
!    Read an output file from electronic structure code using templates
!  SYNOPSIS
!    subroutine ReadOutput(par,n,cfnamein,cfnameout)
!***
      subroutine ReadOutput(par,n,cfnamein,cfnameout)
c
c     tjm 5/05
c     Attempt to write a template-driven output file reader
c     This is like the template-driven input deck writer in WriteInput
c
c     Keys to building template:
c     ^iiistring - find iiith occurrence of string after start of line
c     *iiistring - find iiith occurrence of string anywhere in a line 
c                  Note that we only count one occurrence per line
c     @iii - Skip iii lines
c     &%iiffffjjjkk%llggggmmmnn - Using format string ffff of length ii, 
c                  read par(j) starting at position kk, then do same for 
c                  par(m) with format string gggg of length ll starting at 
c                  position nn
c     !string%iiffffjjjkk%llggggmmmnn - Like above, but will search for line
c                  beginning with string 
c     Next one not yet implemented...
c     %string%%iiistring%%jjj - find string, then first number 
c                  goes into par(i), string and second number into par(j)
c
c     cfnamein - name of template file
c     cfnameout - name of output file
c
      include 'implic'
      integer n,ioIn,ioOut
      integer itimes,itmp,ntimes,ifound,ich,icomp,i
      integer lenfmt,ivar,ipos,lstkey
      real*8 par(n)
      parameter (ioIn=10)
      parameter (ioOut=11)
      character*256 cfnamein,cfnameout,cbuffer
      character(len=256) :: cbufout
      character*256 cfmt,ctmp,ckeystr

c      print *,cfnamein,cfnameout

      open(ioIn,file=trim(cfnamein))
      open(ioOut,file=trim(cfnameout))

 1000 format(a256)
 2000 read(ioIn,1000,end=5000) cbuffer

      if (cbuffer(1:1) .eq. '^') then
c         print *, 'Found key ^'
         read(cbuffer(2:4),'(i3)') itimes
         ctmp=cbuffer(5:len(cbuffer))//' '
         itmp=len_trim(cbuffer)-4
         ntimes=0
 3000    read(ioOut,1000,end=7000) cbufout
         if (trim(cbufout(1:itmp)) .eq. trim(ctmp)) ntimes=ntimes+1
         if (ntimes .lt. itimes) goto 3000
c         print *, 'Found ',itimes,' occs of x',trim(ctmp),'x'
         goto 2000
      endif

      if (cbuffer(1:1) .eq. '*') then
c         print *, 'Found key *'
         read(cbuffer(2:4),'(i3)') itimes
         ctmp=cbuffer(5:len(cbuffer))//' '
         itmp=len_trim(cbuffer)-4
         ntimes=0
 3001    read(ioOut,1000,end=7000) cbufout
         ifound=0
         ich=0
 3002    ich=ich+1
         icomp=min(len_trim(cbufout)-ich+1,len_trim(ctmp))
         if (trim(cbufout(ich:ich+icomp-1)) .eq. trim(ctmp(1:icomp)))
     $        ifound=1
         if (ifound .eq. 1) then
            ntimes=ntimes+1
            goto 3003
         endif
         if (ich .lt. len_trim(cbufout)) goto 3002
 3003    if (ntimes .lt. itimes) goto 3001
c         print *, 'Found ',itimes,' lines with x*',trim(ctmp),'x'
         goto 2000
      endif

      if (cbuffer(1:1) .eq. '&') then
c         print *, 'Found key &'
 3005    read(ioOut,1000,end=7000) cbufout
         i=2
 3006    read(cbuffer(i+1:i+2),'(i2)') lenfmt
         read(cbuffer(i+3:i+3+lenfmt-1),*) cfmt
         cfmt=cfmt//' '
         read(cbuffer(i+3+lenfmt:i+5+lenfmt),'(i3)') ivar
         if ((ivar .gt. n) .or. (ivar .lt. 1)) then
c            print *, 'ivar out of bounds: ',ivar,n
            stop
         endif
         read(cbuffer(i+6+lenfmt:i+7+lenfmt),'(i2)') ipos
         read(cbufout(ipos:),cfmt) par(ivar)
c         print *, 'read ',ivar,' variable as: ',par(ivar)
         if (cbuffer(i+8+lenfmt:i+8+lenfmt) .eq. '%') then
            i=i+8+lenfmt
            goto 3006
         endif
         goto 2000
      endif

      if (cbuffer(1:1) .eq. '!') then
c         print *, 'Found key !'
         i=2
         lstkey=0
 3007    if (cbuffer(i:i) .eq. '%') then
            lstkey=i-1
         else
            i=i+1
            if (i .gt. len_trim(cbuffer)) then
c               print *, 'Error interpreting key: ',trim(cbuffer)
               stop
            endif
            goto 3007
         endif
         if (lstkey .eq. 0) then
c            print *, 'Error interpreting key: ',trim(cbuffer)
            stop
         endif
         ckeystr=cbuffer(2:lstkey)//' '
 3008    read(ioOut,1000,end=7000) cbufout 
         if (cbufout(1:lstkey-1) .ne. ckeystr(1:lstkey-1)) goto 3008
         i=lstkey+1
 3009    read(cbuffer(i+1:i+2),'(i2)') lenfmt
         read(cbuffer(i+3:i+3+lenfmt-1),*) cfmt
         cfmt=cfmt//' '
         read(cbuffer(i+3+lenfmt:i+5+lenfmt),'(i3)') ivar
         if ((ivar .gt. n) .or. (ivar .lt. 1)) then
c            print *, 'ivar out of bounds: ',ivar,n
            stop
         endif
         read(cbuffer(i+6+lenfmt:i+7+lenfmt),'(i2)') ipos
         read(cbufout(ipos:),cfmt) par(ivar)
c         print *, 'read ',ivar,' variable as: ',par(ivar)
         if (cbuffer(i+8+lenfmt:i+8+lenfmt) .eq. '%') then
            i=i+8+lenfmt
            goto 3009
         endif
         goto 2000
      endif

      if (cbuffer(1:1) .eq. '@') then
c         print *, 'Found key @'
         read(cbuffer(2:4),'(i3)') itimes
c         print *, 'read ',itimes,' lines'
         do i=1,itimes
            read(ioOut,1000,end=7000) cbufout
         enddo
c         print *, 'done'
         goto 2000
      endif

      print *, 'Ignoring Unrecognized key: ',cbuffer(1:1)
      goto 2000

 7000 continue                  ! EOF reached on output file read
      print *, 'Premature EOF on output file'
      stop

 5000 continue                  ! EOF reached on input file
      close(ioOut)
      close(ioIn)
      return
      end
