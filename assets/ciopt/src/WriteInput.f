!!****h* CIOpt/WriteInput.f
!  NAME
!    WriteInput.f
!*******
!!****f* WriteInput.f/WriteInput
!  NAME
!    WriteInput
!  DESCRIPTION
!    Write an input file for electronic structure code using templates
!  SYNOPSIS
!    subroutine WriteInput(par,n,ixstate,cfnamein,cfnameout)
!***
      subroutine WriteInput(par,n,ixstate,cfnamein,cfnameout)
c
c     tjm 5/05
c     Attempt to write a template-driven input deck writer
c     The idea is to read in a template file and substitute all
c     occurrences of %%iii by par(i) where par is an input array
c     We also allow output formats to be specified by writing %%$jjccciii
c     where jj is number of characters in format string given by ccc and
c     output number should be par(i)
c     Some examples (a(i)=i in following):
c
c     Input Line: 
c     This is the first variable: %%001 and this is the second: %%002
c     Output Line:
c     This is the first variable: 1.000000 and this is the second: 2.000000
c     Input Line: 
c     First %%$6(f5.3)001 and second %%002
c     Output Line:
c     First 01.000 and second 2.000000
c
c     Certain special variables can also be substituted in:
c     %#STATE  = ixstate
c     %#ISTATE = istate
c     %#JSTATE = jstate
c
c     There is a parallel routine ReadOutput which needs to agree with 
c     this one in broad terms
c
      include 'implic'
      include 'efunc.cmn'

      integer ioIn,ioOut,n,ixstate,i,lstout,lenfmt,ivar
      real*8 par(n)
      parameter (ioIn=10)
      parameter (ioOut=11)
      character*256 cfnamein,cfnameout,cbuffer
      character(len=256) :: cbufout
      character*256 cfmt,ctmp

      open(ioIn,file=trim(cfnamein))
      open(ioOut,file=trim(cfnameout))

      do i=1,256
         cbufout(i:i)=' '
      enddo

 1000 format(a256)
 2000 read(ioIn,1000,end=5000) cbuffer
      i=1
      lstout=1
      if (len_trim(cbuffer) .eq. 0) then
         write(ioOut,'(a1)') ' ' 
         goto 2000
      endif
 2010 cfmt=' '
      if (cbuffer(i:i+1) .eq. '%%') then
         if (cbuffer(i+2:i+2) .eq. '$') then
            read(cbuffer(i+3:i+5),'(i3)') lenfmt
            read(cbuffer(i+6:i+6+lenfmt),*) cfmt
            i=i+6+lenfmt
            read(cbuffer(i+1:i+3),'(i3)') ivar
            if (ivar .gt. n) then
               print *, 'Variable out of range'
               stop
            endif
            cfmt=cfmt//' '
            write(ctmp,fmt=cfmt) par(ivar)
            cbufout=cbufout(1:lstout-1)//trim(ctmp)
            lstout=lstout+lenfmt+1
            i=i+3
         else
            read(cbuffer(i+2:i+4),'(i3)') ivar
            if (ivar .gt. n) then
               print *, 'Variable out of range'
               stop
            endif
            write(ctmp,'(f15.8) ') par(ivar)
            cbufout=cbufout(1:lstout-1)//trim(ctmp)
            lstout=lstout+15
            i=i+4
         endif
      elseif (cbuffer(i:i+1) .eq. '%#') then
         select case (cbuffer(i+2:i+7))
            case ('STATE ')
               write(ctmp,'(i3)') ixstate
            case ('ISTATE')
               write(ctmp,'(i3)') istate
               print *, 'istate = ',istate
               call flush(6)
            case ('JSTATE')
               write(ctmp,'(i3)') jstate
               print *, 'jstate = ',jstate
               call flush(6)
            case ('KSTATE')
               write(ctmp,'(i3)') kstate
               print *, 'jstate = ',kstate
               call flush(6)
            case default
               print *, 'Unknown variable for substitution: ',
     $              cbuffer(i+2:i+7)
               stop
         end select
         cbufout=cbufout(1:lstout-1)//trim(ctmp)
         lstout=lstout+3
         i=i+7
      else
         cbufout(lstout:lstout)=cbuffer(i:i)
         lstout=lstout+1
      endif
      i=i+1
      if (i .le. len_trim(cbuffer)) goto 2010 ! Not done with line
      write(ioOut,'(256A)') trim(cbufout)
      do i=1,len(cbufout)
         cbufout(i:i)=' '
      enddo
      goto 2000                 ! Not done with input file yet

 5000 continue                  ! EOF reached on input file
      close(ioOut)
      close(ioIn)
      return
      end
