      subroutine writexyz(p)

      include 'implic'
      include 'efunc.cmn'

      integer i,k
      real*8 p(MaxDims)
      character*256 cbuf1,cfname

      include 'tjmfunc.cmn'

c we only want to write xyzs if people want tree-killer output
      if (.not.zdetails) return

c     We will not create xyz files if coordinates are not Cartesian
      if (znoncart) return

      if (nescall .lt. 10) then
         write(cbuf1,'(i1)') nescall
      elseif (nescall .lt. 100) then
         write(cbuf1,'(i2)') nescall
      elseif (nescall .lt. 1000) then
         write(cbuf1,'(i3)') nescall
      else
         cbuf1='xxx'
      endif
      cfname='details/'//trim(cbuf1)//'.xyz'
      open(17,file=trim(cfname),status='unknown',form='formatted')
      write(17,*) natoms
      write(17,*) ' '
      do i=1,natoms
         write(17,1004) element(i),(p(k),k=i*3-2,i*3)
      enddo
 1004 format(1x,a1,3f16.10)
      close(17)
      return
      end subroutine
