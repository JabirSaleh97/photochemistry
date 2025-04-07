!!****h* CIOpt/ciopt.f
! NAME
!    CIOpt - Optimize conical intersections
! DESCRIPTION
!     Optimize conical intersection geometries without nonadiabatic coupling vectors
!     This code will also work without analytic gradients. It can optimize stationary points or
!     conical intersections.  It is meant to drive an electronic structure code and is set up
!     to work with many different codes using a flexible template-driven scheme for parsing
!     output files and writing input decks.  The algorithm used is described in 
!     Levine, Ko, Quenneville, and Martinez, Mol. Phys. v104 p1053 (2006).
! TODO
!     Spruce up multiplier penalty method as described in Fletcher.  Apply it to 
!     CASPT2 and MSCASPT2 and publish a short note. 
! COPYRIGHT 
!   Ben Levine and Todd J. Martinez, 2003-2005
!   This code may not be reproduced in whole or in part without express written permission
!   from the copyright holders.
! AUTHOR 
!   Ben Levine and Todd J. Martinez
!*******
      program ciopt

      include 'implic'
      include 'efunc.cmn'
      include 'dfphess.cmn'

      character*20 cversion
      real*8 p(MaxDims),xi(MaxDims,MaxDims)
      real*8 grad(MaxDims),grad2(MaxDims),grad3(MaxDims)
      real*8 gviol(MaxCons)
      integer iter,i,j,k,itop,imode,icons,igapok,iconv
      integer icleanup,itmp
      real*8 tol,fret,step,en,gap,gapold,tolramp,tolmin,cigap
      real*8 tolfirstsigma,gtolfirstsigma
      real*8 gtolramp,gtolmax
      logical ztolramp,zgtolramp,ztmp,zmat,zsmartsigma
      real*8 dviol,spenmax,dval,gtol,dlambdagapmax
      real*8 EvalCons,tmp,hessinguess
      external EvalCons

      data istate/0/
      data jstate/0/
      data kstate/0/
      data step/1.0d-2/
      data stepnd/0.01d0/
      data stepmax/0.02d0/
      data ndims/0/
      data alpha/0.025/
      data dlambdagap/3.5001d0/
      data dlambdagapmax/100.01d0/
      data spen/1.0d0/
      data spenmax/10000.0d0/
      data tol/1.001d-6/
      data gtol/5.001d-3/
      data tolramp/0.5/
      data tolmin/1.0d-5/
      data cigap/1.001d-3/
      data alpharmsd/.1d0/
      data nefunc/7/ 
      data nopt/3/
      data zmultigrad/.true./
      data zrestart/.false./
      data zangrad/.true./
      data znoncart/.false./
      data zforward/.true./
      data ztolramp/.false./
      data istcallfunc/0/
      data crunstr/
     $     '/usr/localmqm/SVNBin/molpro2002.6/bin/molpro tmp.com'/
      data ctmpread/'template.read'/
      data ctmpgread/'template.readg'/
      data ctmpg2read/'template.readg2'/
      data ctmpg3read/'template.readg3'/
      data ctmpwrite/'template.write'/
      data ctmpwriteg/'template.writeg'/
      data cinpdeck/'tmp.com'/
      data coutfile/'tmp.out'/
      data zdetails/.false./
      data rmsdweights/MAXATOMS*0.0d0/
      data zibf/.false./
      data ribf/0.04001d0/
      data Emaxibf/1.001d10/
      data dlambdarmsd/1.001d0/
      data ZExEnEv/.false./
      data zmat/.true./
      data hessinguess/.1001d0/
      data zsmartsigma/.true./
      data tolfirstsigma/1.00e-5/
      data gtolfirstsigma/1.00d-2/

      namelist /control/ dlambdagap,alpha,nefunc,istate,jstate,kstate,
     $     natoms,tol,nstates,ndims,step,nopt,zangrad,znoncart,stepnd,
     $     spen,cigap,zforward,ztolramp,tolramp,tolmin,
     $     zrestart,crunstr,ctmpread,ctmpgread,ctmpg2read,
     $     ctmpg3read,ctmpwrite,zmultigrad,coutfile,cinpdeck,
     $     ngeocons,igeocons,iacons,dgeocons,dlambdag,speng,
     $     dgcthresh,iskip,zcheckgrad,zdetails,gtol,dlambdagapmax,
     $     rmsdweights,rmsdgeo,alpharmsd,zibf,ribf,Emaxibf,
     $     dlambdarmsd,ZExEnEv,stepmax,zmat,hessinguess,zsmartsigma,
     $     tolfirstsigma,gtolfirstsigma

c------------------------------------------------------------------------------------------------------
c     Namelist variables
c------------------------------------------------------------------------------------------------------
c     dlambdagap - Fixed lagrange multiplier (not optimized) for nefunc={2,3}; Initial guess value in cases
c                  when the multiplier is optimized (intersection search w/ zlagrnge=.true.)
c     nefunc     - Choose energy functional to optimize
c                         EAve = 0.5(E(istate)+E(jstate))
c                         gap = E(istate)-E(jstate)
c                     1 - Efunc=E(istate)
c                     6 - Efunc=E(istate)+dlambdagap*((gap*gap)/(gap+alpha)) is minimized first.  Then
c                            we switch to nefunc=9.  I.e., this is just a way to get a good guess for
c                            the multiplier penalty method 
c                     7 - Eq. 14 of TDDFT paper
c                            Efunc=EAve+dlambdagap*((gap*gap)/(gap+alpha))
c                     8 - Sequential penalty method (dlambdagap adjusted automatically)
c                            Efunc=EAve+0.5*dlambdagap*gap*gap
c                     9 - Multiplier penalty method (spen and dlambdagap adjusted automatically)
c                            Efunc=EAve-dlambdagap*gap+0.5*spen*gap*gap
c                     10 - The equivalent of 7 for the optimization of three state intersections
c                           Efunc=EAveijk+dlambdagap/3*(((gapij*gapij)/(gapij+alpha))+((gapik*gapik)/(gapik+alpha))+((gapjk*gapjk)/(gapjk+alpha)))
c     alpha      - See Eq 14 of TDDFT paper, only used for nefunc=7
c     istate     - The state of interest for optimization (upper state if CI optimization)
c     jstate     - Lower state in CI optimization (defaults to istate-1)
c     natoms     - Number of atoms in molecule
c     tol
c     cigap      - Largest gap we will accept as a conical intersection for sequential and multiplier
c                  penalty methods
c     nstates    - Number of states which will be solved for by electronic structure program
c     ndims      - Number of variables to optimize.  If set to zero (default), the code assumes
c                  you are optimizing in Cartesian coordinates and the number of optimized variables
c                  is 3*natoms.  If set to non-zero value, input format for initial geometry is changed
c     step       - Initial step sizes for variable optimization (used in powell, not sure about cg)
c     stepnd     - Step size for numerical differentiation
c     iskip(i)   - Only used with Powell.  Set iskip(i)=1 to skip line minimization for ith DOF
c     nopt       - 1 = Use Powell
c                - 2 = Use conjugate gradient
c                - 3 = Use BFGS
c     spen       - Initial value of penalty in multiplier penalty method
c     ngeocons   - Number of geometry constraints (Cartesian coordinate optimization)
c     igeocons(i)- Type of ith constraint 
c                  1 - Bond length
c                  2 - Bond Angle
c                  3 - Dihedral Angle
c     iacons(i,j)- Atoms involved in jth geometric constraint
c     dgeocons(i)- Value for jth geometric constraint
c     dlambdag(i)- Lagrange multiplier term for ith geometric constraint
c     speng(i)   - Penalty term for ith geometric constraint
c     zrestart   - If true, we will attempt to restart optimization
c     zangrad    - If true, analytic gradients can be computed and read
c     znoncart   - If true, non-Cartesian coordinates are being used and xyz files will not be generated
c     zforward   - If true, forward differences will be used; default is central differences
c     ztolramp   - If true and using penalty method, tolerance is made tighter each time penalty is increased
c     zmultigrad - If true, we can read two gradients from a single calculation.  In this case, the 
c                  first gradient is read with ctmpgread and the second with ctmpg2read
c     zcheckgrad - if true, analytic gradients will be checked against finite differences.  The code will
c                  does not analyze the agreement or lack thereof, but rather just compares the analytic and
c                  numerical results in stdout 
c     crunstr    - String to run electronic structure code
c     ctmpread   - file name for template when reading energy from output
c     ctmpgread  - file name for template when reading analytic gradient from output
c     ctmpg2read - file name for template when reading second gradient from output
c     ctmpg2read - file name for template when reading third gradient from output
c     ctmpwrite  - file name for template when writing input decks (energy only)
c     ctmpwriteg - file name for template when writing input decks (energy+gradient)
c     cinpdeck   - file name to use for input deck
c     coutfile   - file name of output from electronic structure code
c------------------------------------------------------------------------------------------------------
c     Further notes about input file (Control.dat)
c     After the namelist, the code expects to read in geometry information 
c     If znoncart is false (default) this is done in an xyz-like format with element name, x, y, z 
c     Otherwise, we just read the initial values of optimization variables, one per line with optional
c     trailing comments
c------------------------------------------------------------------------------------------------------

      print *, '                *** CIOpt *** '
      print *, 'Intersection Optimization and Energy Minimization'
      print *, 'without analytic derivatives or nonadiabatic coupling'
      print *, '      - Ben Levine and Todd Martinez - '
      print *, '            Copyright 2004-2005        '
      call svn_version(cversion)
      print *, ' SVN Version: ',trim(cversion)
      print *, ' '

      open(11,file='Control.dat',status='old',form='formatted')
      read(11,control)

      if (nopt.eq.2) then
         print *, 'Conjugate Gradient (nopt=2) is no longer available'
         stop
      endif

      if (istate .eq. 0) then
         print *, 'You must set istate'
         stop
      endif
      if (jstate .eq. 0.and.nefunc.ne.1) jstate=istate-1
      if (kstate .eq. 0.and.nefunc.eq.10) kstate=istate-2

      if (zdetails) then
         call system('mkdir details')
         call system('rm details/*')
      endif

      do i=1,MaxDims
         do j=1,MaxDims
            xi(i,j)=0.0d0
         enddo
      enddo

      if (zrestart) then
         open(37,file='long.out',status='old',form='formatted')
      endif

      if (ngeocons .gt. MaxCons) then
         print *, 'Too many geometric constraints: ',ngeocons
         print *, 'Maximum is: ',MaxCons
         stop
      endif

ctjm  This is underhanded trickery to make the user unaware of defaults
ctjm  which change depending on method...  
      if (nopt .ne. 3) then
         if (ztolramp) then
            if (tol .eq. 1.001d-06) tol=0.01
         endif
         if (tol .eq. 1.001d-06) then
            if (nefunc.eq.11) then
               tol=1.0d-5
            else
               tol=1.0d-5
            endif
         endif
      else
         if (tol .eq. 1.001d-06) then
            if (ztolramp) then
               tol=0.2
               tolmin=0.01
c            else
c               tol=0.01
            endif
         endif
      endif
      if (gtol.eq.5.001d-3) then
         if (nefunc.eq.11) then
            gtol=5.0d-3
         else 
            gtol=5.0d-3
         endif
      endif
      if (nefunc .eq. 9) then
         if (dlambdagap .eq. 3.5001d0) dlambdagap=0.0d0
      endif
      if (dlambdagap .eq. 3.5001d0) then
         if (nefunc.eq.7.or.nefunc.eq.12) dlambdagap=3.5d0
         if (nefunc.eq.8) dlambdagap=200.0d0
      endif
      if (cigap.eq.1.001d-3) then
         if (nefunc.eq.10) then
            cigap=4.0d-3
         else
            cigap=1.0d-3
         endif
      endif
      if (dlambdagapmax .eq. 100.01d0) then
         if (nefunc.eq.7) dlambdagapmax=100.0d0
         if (nefunc.eq.8) dlambdagapmax=5000.0d0
      endif
      icleanup=0                ! Set to 1 when cleanup optimization is finished
      bfgseconv=1.d-3           ! Energy convergence threshold in BFGS

      if (hessinguess.eq..1001d0) then
         if (nefunc.eq.11) then
            hessinguess=1.0d0
         else
            hessinguess=.1d0
         endif
      endif

      if (.not. znoncart) then         
         print *, 'Optimizing in Cartesian coordinates'
         ndims=natoms*3
         do i=1,natoms
            read(11,*) element(i),(p(k),k=i*3-2,i*3)
         enddo
         do i=1,3*natoms
            scalefac(i)=step
            scalend(i)=stepnd
            hessini(i)=hessinguess
         enddo
      else
         print *, 'Optimizing in user-defined coordinates'
         do i=1,ndims
            read(11,*) p(i),scalefac(i),scalend(i),hessini(i)
         enddo
         write(6,1090)
 1090    format('Coord   ScaleFac(Powell)   FiniteDiffStep   InvHess')
         do i=1,ndims
            if (scalefac(i) .eq. 0.0d0) then
               if (zmat) then
                  print *, "Step sizes chosen assuming z-matrix"
                  if (i.eq.2. .or. MOD(i,3).eq.1) then
                     scalefac(i)=step
                  else
                     scalefac(i)=step*10
                  endif
               else
                  scalefac(i)=step
               endif
               if (scalend(i) .eq. 0.0d0) scalend(i)=scalefac(i)
               if (hessini(i) .eq. 0.0d0) hessini(i)=5000.0d0*
     $              scalefac(i)
            else
               if (scalend(i) .eq. 0.0d0) scalend(i)=step
               if (hessini(i) .eq. 0.0d0) hessini(i)=5000.0d0*step
            endif
            write(6,1091) i,scalefac(i),scalend(i),hessini(i)
         enddo
 1091    format(i5,3(1x,f14.8))
         print *, '  '
      endif

c      if (dlambdarmsd.eq.1.001d0) dlambdarmsd = 

      if (nefunc.eq.11.or.nefunc.eq.12) then
         ztmp=.false.
         do i=1,natoms*3
            if (abs(rmsdgeo(i)).gt.1.0d-7) ztmp=.true.
         enddo
         if (.not.ztmp) then
            print *, 'you must define rmsdgeo with nefunc=11'
            stop
         endif
         itmp=0
         do i=1,natoms
         print *,rmsdweights(i)
            if (rmsdweights(i).gt.0.0d0) itmp=itmp+1
         enddo
         print *,itmp
         if (itmp.eq.0) then
            do i=1,natoms
               if (element(i).eq.'H') then
                  rmsdweights(i)=1.0d0
               elseif (element(i).eq.'C') then
                  rmsdweights(i)=12.0d0
               elseif (element(i).eq.'N') then
                  rmsdweights(i)=14.0d0
               elseif (element(i).eq.'O') then
                  rmsdweights(i)=16.0d0
               elseif (element(i).eq.'S') then
                  rmsdweights(i)=32.0d0
               else
                  print *, 'Unrecognized atom type.'
                  print *, 'to use nefunc=11 you must use only H, C,',
     $                 'N, and O, or input rmsdweights'
                  stop
               endif
            enddo            
         elseif (itmp.ne.ndims) then
            print *, 'some rmsdweights are zero.  Exiting'
            stop
         endif
         if (dlambdagap .eq. 3.5001d0.and.nefunc.eq.11) then
            tmp=0.0d0
            do i=1,natoms
               tmp=tmp+rmsdweights(i)
            enddo
            tmp=sqrt(tmp)
            dlambdagap=10.0d0*tmp
         endif
         if (dlambdagapmax .eq. 100.01d0.and.nefunc.eq.11) then
            tmp=0.0d0
            do i=1,natoms
               tmp=tmp+rmsdweights(i)
            enddo
            tmp=sqrt(tmp)
            dlambdagapmax=1000d0*tmp
         endif
         if (dlambdarmsd.eq.1.001d0.and.nefunc.eq.12) then
            tmp=0.0d0
            do i=1,natoms
               tmp=tmp+rmsdweights(i)
            enddo
            tmp=sqrt(tmp)
            dlambdarmsd=.5d0/tmp
         endif
      endif


      select case (nefunc)
         case (1)
            print *, 'Energy minimization'
         case(6)
            print *, 'Excited state minimization followed',
     $           ' by intersection search'
            print *, 'Ben doesn''t like nefunc=6, so he turned it off'
            print *, 'You should probably use nefunc=7 instead'
            stop
         case (7)
            print *, 'Intersection search with smoothed gap'
            zseqpen=.true.
         case (8)
            print *, 'Sequential penalty intersection search'
            zseqpen=.true.
            print *, 'nefunc=8 is an interesting choice...'
            print *, 'You might prefer nefunc=7, as it tends to',
     $           ' be more efficient'
         case (9)
            print *, 'Multiplier penalty intersection search'
            zseqpen=.false.
            print *, 'Ben doesn''t like nefunc=9, so he turned it off'
            print *, 'You should probably use nefunc=7 instead'
c            stop
         case (10)
            print *, '3-state intersection search with smoothed gap'
            zseqpen=.true.
         case (11)
            print *, 'Closest intersection search with smoothed gap'
            zseqpen=.true.
         case (12)
            print *, 'MECI subject to rmsd constraint with smoothed gap'
            zseqpen=.true.
      end select

      if (ngeocons .ne. 0) then
c         if (nefunc .ne. 9) then
c            print *, 'Geometric constraints only with ',
c     $           'nefunc=9 for now'
c            stop
c         endif
         print *, 'Enforcing geometric constraints'
         do i=1,ngeocons
            if (dgcthresh(i) .eq. 0.0d0) then
               select case (igeocons(i))
                  case (1)
                     dgcthresh(i)=0.01
                  case (2)
                     dgcthresh(i)=0.017
                  case (3)
                     dgcthresh(i)=0.017
               end select
            endif
         enddo
         open(47,file='cons.log')
         close(47,status='delete')
         open(47,file='cons.log')
         write(47,1891)
         close(47)
 1891    format('#Iter  Cons     Current        Target')
         write(6,1892)
 1892    format('Cons   Type Atom1 Atom2 Atom3 Atom4   Initial   ',
     $        ' Target    Lambda      Spen    Thresh')
         do i=1,ngeocons
            itop=igeocons(i)+1
            do j=1,itop
               if (iacons(j,i) .eq. 0) then
                  print *, 'Error in geometric constraints'
                  stop
               endif
            enddo
            dval=EvalCons(p,i)
c            if (zseqpen) then
c               dlambdag(i)=1.0
c            else
Cremove               dlambdag(i)=0.0
Cremove               speng(i)=1.0
c            endif
            write(6,1893) i,igeocons(i),(iacons(j,i),j=1,4),dval,
     $           dgeocons(i),dlambdag(i),speng(i),dgcthresh(i)
         enddo
      endif
 1893 format(6(i5,1x),5(f9.5,1x))

c     Set up directions for Powell.  It would be nice to allow for different initial step
c     sizes along different directions, even better to automatically determine this.      
c     If we are working in zmatrix coordinates we put the angles first.
c     This seems to improve convergence.
      if (znoncart.and.zmat) then
         j=1
         k=ndims
         do i=1,ndims
            if (i.eq.2. .or. MOD(i,3).eq.1) then
               xi(k,i)=scalefac(i)
               k=k-1
            else
               xi(j,i)=scalefac(i)
               j=j+1
            endif
         enddo
      else
         do i=1,ndims
            xi(i,i)=scalefac(i)
            print *, xi(i,i)
         enddo
      endif

c     Done reading input file
      close(11)

      open(45,file='linmin.out')
      close(45,status='DELETE')
      open(45,file='mplog.out')
      close(45,status='DELETE')
      if (.not.zrestart) then
         open(16,file='long.out',status='unknown')
      endif
      close(16,status='delete')

      open(46,file='iter.log')
      close(46,status='delete')
      open(46,file='iter.log')
      select case (nopt)
         case (1)
            print *, '*** Using Powell Optimization Scheme ***'
            write(46,1083)
 1083       format('# Iter         F             DeltaF      ',
     $           '       Ei          Ej')
         case (2)
            print *, '* Using Conj Gradient Optimization Scheme *'
            if (nefunc.ne.10) then
               write(46,1085)
            else
               write(46,1086)
            endif
         case (3)
            print *, '*** Using BFGS Optimization Scheme ***'
            if (nefunc.ne.10.and.nefunc.ne.11.and.nefunc.ne.12) then
               write(46,1085)
 1085          format('# Iter         F            Fstep       ',
     $              '     |grad|        |grad(gap)|     ',
     $              '|grad(perp)|     Ei              Ej ')
            elseif (nefunc.eq.11.or.nefunc.eq.10) then
               write(46,1087)
 1087          format('# Iter         F            rmsd        ',
     $              '    Fstep       ',
     $              '     |grad|        |grad(gap)|     ',
     $              '|grad(perp)|     Ei              Ej ')
            else
               write(46,1086)
 1086          format('# Iter         F            Fstep       ',
     $              '     |grad|        |grad(gap)|     |grad(perp)|',
     $              '     Ei              Ej              Ek')
            endif
         case default
            print *, 'Unknown Optimization Scheme requested: ',nopt
            stop
      end select
      close(46)
      if (zangrad) then
         print *, 'Analytic gradients will be used'
      else
         if (zforward) then 
            print *, 'Gradients by forward difference'
         else
            print *, 'Gradients by central difference'
         endif
      endif

ctjm  Case where we optimize with smoothed functional and then switch over to multiplier
ctjm  penalty scheme. The proposed advantage is that we emphasize minimization of 
ctjm  excited state energy (while avoiding derivative discontinuities near intersections)
ctjm  before we start using the average energy...
      if (nefunc .eq. 6) then
         print *, 'Minimize excited state w/smoothed functional...'
         select case (nopt)
            case (1)
c               call powell(p,xi,ndims,MaxDims,tol,iter,fret)
               call powellb(p,xi,tol,iter,fret)
c            case (2)
c               call frprmn(p,ndims,tol,gtol,iter,fret)
            case (3)
c               call dfpmin(p,ndims,tol,gtol,iter,fret)
               call sdrive(p,ndims,tol,gtol)
            case default
               print *, 'Unknown Optimization Scheme requested: ',nopt
               stop
         end select            
         dlambdagap=-0.5d0
         spen=10.0
         if (ztolramp) then
            if (nopt .ne. 3) then
               tol=0.001
            else
               tol=0.1
            endif
         endif
         nefunc=9
ctjm     Ensure that we reset the hessian inverse if using BFGS
         if (nopt .eq. 3) istcalldfp=0
ctjm     Now we should check that we are not already converged...
         gap=envals(istate)-envals(jstate)
         if (gap .lt. cigap) then
            print *, ' '
            print *, 'Ei=',envals(istate),' Ej=',envals(jstate)
            print *, 'Gap=',gap,' CI Threshold=',cigap
            print *, ' '
            print *, 'Gap low enough to qualify result of excited'
            print *, 'state minimization as a conical intersection'
            print *, 'No refinement of intersection will be sought.'
            print *, 'If you are paranoid, you should restart this'
            print *, 'calculation using a sequential or multiplier'
            print *, 'penalty method without preminimization of '
            print *, 'the excited state minimum using the smoothed'
            print *, 'function.'
            goto 9800
         endif
      endif

ctjm  If doing multiplier penalty we need to know the initial value of the gap
      if (nefunc .eq. 9) then
         imode=0
         call runmp(p,en,gap,grad,grad2,grad3,istate,imode)
         gapold=gap
         print *, 'Multiplier penalty, gap = ',gap
      endif

      do icons=1,ngeocons
         gviol(icons)=EvalCons(p,icons)-dgeocons(icons)
      enddo
      
cbgl try to inteligently choose sigma
      if (zsmartsigma) then
         print *, 'Running optimization to ',
     $        'choose dlambdagap'
         call sdrive(p,ndims,tolfirstsigma,gtolfirstsigma)
         gap=abs(envals(istate)-envals(jstate))
         tmp=cigap/alpha+1.0d0
         dlambdagap=-1.5*tmp*tmp/(dZsmart*(tmp*tmp-1))
c         print *, alpha, tmp, dZsmart
c         print *, 'gaps',gap,alpha*(sqrt(dlambdagap*dZsmart/
c     $              (1+dlambdagap*dZsmart))-1.0d0)
         print *, 'done with first optimization'
         print *, 'dlambdagap = ', dlambdagap
      endif
         
 5000 continue
      select case (nopt)
         case (1)
c            call powell(p,xi,ndims,MaxDims,tol,iter,fret)
            call powellb(p,xi,tol,iter,fret)
c         case (2)
c            call frprmn(p,ndims,tol,gtol,iter,fret)
         case (3)
c            call dfpmin(p,ndims,tol,gtol,iter,fret)
            call sdrive(p,ndims,tol,gtol)
         case default
            print *, 'Unknown Optimization Scheme requested: ',nopt
            stop
      end select

      call writexyz(p)
      imode=0
      call runmp(p,en,gap,grad,grad2,grad3,istate,imode)
      write(6,2000) en,gap
 2000 format('FINAL GEOMETRY, E=',f15.8,', gap=',f15.8)

      if (.not. znoncart) then
         do i=1,natoms
            write(6,1004) element(i),(p(k),k=i*3-2,i*3)
         enddo
      else
         do i=1,ndims
            write(6,*) i,p(i)
         enddo
      endif
 1004 format(1x,a1,3f16.10)

ctjm  Reset the hessian inverse matrix if we are ramping a penalty
      istcalldfp=0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c new code for geometric constraints
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         iconv=1
         do icons=1,ngeocons
            dviol=EvalCons(p,icons)-dgeocons(icons)
            if (abs(dviol) .gt. dgcthresh(icons)) iconv=0
         enddo
         if (iconv .eq. 0) then
            do icons=1,ngeocons
               dviol=EvalCons(p,icons)-dgeocons(icons)
               if (dviol > 0.25*gviol(icons) .and.
     $              abs(dviol) .gt. dgcthresh(icons)) then
                  speng(icons)=4.0d0*speng(icons)
                  print *, 'Increase geom penalty ',icons,
     $                 gviol(icons),dviol,speng(icons)
                  gviol(icons)=dviol
               else
                  dlambdag(icons)=dlambdag(icons)-speng(icons)*dviol
                  print *, 'New linear geom penalty term ',icons,
     $                 gviol(icons),dviol,dlambdag(icons)
                  gviol(icons)=dviol
               endif
            enddo
            write(6,1892)
            do i=1,ngeocons
               write(6,1893) i,igeocons(i),(iacons(j,i),j=1,4),dval,
     $              dgeocons(i),dlambdag(i),speng(i),dgcthresh(i)
            enddo
            goto 5000
         endif
         

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ctjm  Now for sequential penalty.  Need to go through the minimizer again with an adjusted penalty.

      if (nefunc .eq. 8 .or. nefunc .eq. 7 .or. nefunc .eq. 10 .or. 
     $        nefunc .eq. 11 .or. nefunc .eq. 12 ) then
         gap=abs(envals(istate)-envals(jstate))
         if (nefunc .eq. 10) gap=abs(envals(istate)-envals(kstate))
         if (gap .gt. cigap) then
            if (dlambdagap.ge.(dlambdagapmax-.01d0)) then
               print *, 'This is as small as the gap gets.'
               print *, 'Increase dlambdagapmax if you wish, but'
               print *, 'the quality of the MECI may be degraded.'
               print *, 'Exiting'
               stop
            endif
            dlambdagap=dlambdagap*4.0d0
            if (zibf) ribf=ribf/10.0d0
            if (dlambdagap.gt.dlambdagapmax) dlambdagap=dlambdagapmax
            print *, 'Increasing penalty in seq penalty method',
     $           dlambdagap
            if (ztolramp) tol=max(tolmin,tolramp*tol)
            print *, 'New tolerance: ',tol
            print *, 'New gradient tolerance (gtol): ',gtol
            if (spen .gt. spenmax) then
               print *, 'Penalty term too large.  Aborting...'
               stop
            endif
            goto 5000 
         endif
         print *, 'Sequential penalty converged, gap=',gap
      endif

ctjm  Multiplier penalty
      if (nefunc .eq. 9) then
c     Are constraints satisfied to desired tolerances?
         igapok=0
         gap=abs(envals(istate)-envals(jstate))
         iconv=1
         if (gap .gt. cigap) iconv=0
         do icons=1,ngeocons
            dviol=EvalCons(p,icons)-dgeocons(icons)
            if (abs(dviol) .gt. dgcthresh(icons)) iconv=0
         enddo
         if (iconv .eq. 0) then
c     First, deal with energy gap constraint
c            if (gap > 0.25*gapold .and.
c     $           gap .gt. cigap) then
c     Need to determine if we have improved "sufficiently"
            if (gap > 0.8*gapold .and.
     $           gap .gt. cigap) then
               spen=10.0d0*spen
               gapold=gap
               print *, 'Increase gap penalty ',spen
               if (ztolramp) tol=max(tolmin,tolramp*tol)
               print *, 'New tolerance: ',tol
               if (nopt .eq. 3) istcalldfp=2
            else
               dlambdagap=dlambdagap-spen*gap
               gapold=gap
               print *, 'New linear gap penalty term ',dlambdagap
               if (ztolramp) tol=max(tolmin,tolramp*tol)
               print *, 'New tolerance: ',tol
               if (nopt .eq. 3) istcalldfp=2
            endif
c     Now deal with geometric constraints
            do icons=1,ngeocons
               dviol=EvalCons(p,icons)-dgeocons(icons)
               if (dviol > 0.25*gviol(icons) .and.
     $              abs(dviol) .gt. dgcthresh(icons)) then
                  speng(icons)=10.0d0*speng(icons)
                  print *, 'Increase geom penalty ',icons,
     $                 gviol(icons),dviol,speng(icons)
                  gviol(icons)=dviol
               else
                  dlambdag(icons)=dlambdag(icons)-speng(icons)*dviol
                  print *, 'New linear geom penalty term ',icons,
     $                 gviol(icons),dviol,dlambdag(icons)
                  gviol(icons)=dviol
               endif
            enddo
            if (spen .gt. spenmax) then
               print *, 'Gap penalty term too large.  Aborting...'
               stop
            endif
            goto 5000
         else
ctjm     We force one cleanup iteration with relatively tight convergence criteria
            if (icleanup .eq. 0) then
               icleanup=1
               istcalldfp=1
               bfgseconv=1.d-6
               if (ztolramp) tol=tolmin
               if (nopt .eq. 3) istcalldfp=2
               print *, ' Begin Cleanup Optimization Cycle...'
               goto 5000
            else
               write(6,1933) 
 1933       format('Cons     Lambda      Sigma             Deviation',
     $              '      Threshold')
 1934          format(i5,4(f14.8,1x))
 1935          format('GAP  ',4(f14.8,1x))
               write(6,1935) dlambdagap,spen,gap,cigap
               do icons=1,ngeocons
                  dviol=EvalCons(p,icons)-dgeocons(icons)
                  write(6,1934) icons,dlambdag(icons),speng(icons),
     $                 dviol,dgcthresh(icons)
               enddo
            endif
         endif
      endif

 9800 continue
      end program
