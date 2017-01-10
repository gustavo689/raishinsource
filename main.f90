!***********************************************************************
!      RAISHIN code: 3D General Relativistic MHD (3D MPI) version 
!      main program
!      written by Y. Mizuno
!      ver. 140514
!***********************************************************************
!
program main

  use pram !- for input parameter for simulation structure -!
  implicit none
  include 'mpif.h' !- for MPI -!
!
!======================================================================@
!    Difinition for variables
!======================================================================@
!
!- uu: converved variables (cell-center) -!
!- uu(1)=\rho*\gamma -!
!- uu(2)=momentum density of i(x1)-component -!
!- uu(3)=momentum density of j(x2)-component -!
!- uu(4)=momentum density of k(x3)-component -!
!- uu(5)=energy density -!
!- uu(6)=ink (advection scalar), not used -!
!- uu(7)=magnetic field of i-component -!
!- uu(8)=magnetic field of j-component -!
!- uu(9)=magnetic field of k-component -!
!
!- uo: 1 step before data of uu -!
!- u0: initial value of uu -!
!- uusr: store data for restart -!
!- uh, us, un: 1st, 2nd, 3rd step value of uu -!
!- uu*r, uu*l: conserved variables (cell-surface in i,j,k-direction) -!
!
  real(8), allocatable :: uu(:,:,:,:), uo(:,:,:,:), &
          u0(:,:,:,:), uh(:,:,:,:), us(:,:,:,:), un(:,:,:,:), &
          uupr(:,:,:,:), uusr(:,:,:,:), uuir(:,:,:,:), uuil(:,:,:,:), &
          uujr(:,:,:,:), uujl(:,:,:,:), uukr(:,:,:,:), uukl(:,:,:,:)
!
!- uri: primitive variables (cell-center) -!
!- uri(1)=density -!
!- uri(2)=4-velocity of i(x1)-component -!
!- uri(3)=4-velocity of j(x2)-component -!
!- uri(4)=4-velocity of k(x3)-component -!
!- uri(5)=gas pressure -!
!- uri(6)=ink (scalar), not used -!
!- uri(7)=magnetic field of i-component -!
!- uri(8)=magnetic field of j-component -!
!- uri(9)=magnetic field of k-component -!
!
!- urio: 1 step before data of uri -!
!- uri0: initial value of uri -!
!- urisr: store data for restart -!
!- uri*r, uri*l: premitive variables (cell-surface in i,j,k-direction) -!
!
  real(8), allocatable :: uri(:,:,:,:), urio(:,:,:,:), &
          uri0(:,:,:,:), urip(:,:,:,:), urisr(:,:,:,:), &
          uriir(:,:,:,:), urijr(:,:,:,:), urikr(:,:,:,:), &
          uriil(:,:,:,:), urijl(:,:,:,:), urikl(:,:,:,:)
  real(8), allocatable :: uria(:,:,:,:), urib(:,:,:,:)
!
!- Metric terms: gcon: contravariant 4-metric g^ij -!
!                gcov: covariant 4-metric     g_ij -!
!
  real(8), allocatable ::  gcon(:,:,:,:,:), gcov(:,:,:,:,:), &
          gconi(:,:,:,:,:), gconj(:,:,:,:,:), gconk(:,:,:,:,:), &
          gcovi(:,:,:,:,:), gcovj(:,:,:,:,:), gcovk(:,:,:,:,:)
!
!- Metric terms: alpha: lapse function -!
!                beta^i: shift vector  -!
!                detg: determinant g \sqrt{-g} -!
!
  real(8), allocatable ::  detg(:,:,:), detgi(:,:,:), &
                           detgj(:,:,:), detgk(:,:,:)    
!
!- Metric terms: christ: Christfel symbol -!
!
  real(8), allocatable :: christ(:,:,:,:,:,:)
!
!- ww: numerical flux (cell-center) -!
!- wwo: 1step before numerical flux (cell-center) -!
!- ww*r, ww*l: numerical flux (at cell-boundary of i,j,k directions) -!
!
  real(8), allocatable :: ww(:,:,:,:,:), wwo(:,:,:,:,:), &
          wwir(:,:,:,:), wwjr(:,:,:,:), wwkr(:,:,:,:), & 
          wwil(:,:,:,:), wwjl(:,:,:,:), wwkl(:,:,:,:)
!
!- sf: source term -!
!- tenr: energy-momentum tensor -!
!
  real(8), allocatable :: sf(:,:,:,:), tenr(:,:,:,:,:)
!
!- cmax, cmin: maximum and minimum wave speed -!
!
  real(8), allocatable :: cmaxi(:,:,:), cmini(:,:,:), &
          cmaxj(:,:,:), cminj(:,:,:), cmaxk(:,:,:), cmink(:,:,:)
!
!- Primitive variables -!
!- de: density, pr: gas pressure -!
!- util1,util2,util3: 4-velocity -!
!- b1,b2,b3: magnetic field -!
!
  real(8), allocatable :: b1(:,:,:), b2(:,:,:), b3(:,:,:), &
          util1(:,:,:), util2(:,:,:), util3(:,:,:), &
          v1(:,:,:), v2(:,:,:), v3(:,:,:), &
          de(:,:,:), pr(:,:,:)
!
!- Grid: x1: cell-center, x1a: cell-boundary -!
!        dx1: x1(i+1)-x1(i), dx1b=x1a(i+1)-x1a(i) -!
!        akap1=dt/dx1, akap1b=dt/dx1b -!
!        for modified Kerr-Schild coordinates
!        x1b: cell-center position for MKS
!        x1ab: cell-boundary position for MKS  
!
  real(8) :: x1(imax), x2(jmax), x3(kmax), &
             x1a(imax), x2a(jmax), x3a(kmax), &
             x1b(imax), x2b(jmax), x3b(kmax), &
             x1ab(imax), x2ab(jmax), x3ab(kmax), &             
             dx1(imax), dx2(jmax), dx3(kmax), &
             dx1b(imax-1), dx2b(jmax-1), dx3b(kmax-1), &
             dx1ab(imax-1), dx2ab(jmax-1), dx3ab(kmax-1), &             
             akap1(imax), akap2(jmax), akap3(kmax), &
             akap1b(imax-1), akap2b(jmax-1), akap3b(kmax-1)
!
!- for Kolmogorov-like Power-law Spectrum -!
!
  real(8) :: wkn(nkmax), pkn1d(nkmax), pkn2d(nkmax), &
             pkn3d(nkmax), thetan(nkmax), phin(nkmax), &
             thetan1(nkmax), phin1(nkmax)
!
  character*256 :: filename, filename1 !- filenames for outputs -!
  integer :: i, j, k, m, n 
  integer :: is, ie, is1, ie1 !- start & end grid for each cpus in i-direc -!
  integer :: js, je, js1, je1 !- start & end grid for each cpus in j-direc -!
  integer :: ks, ke, ks1, ke1 !- start & end grid for each cpus in k-direc -!
  integer :: iup, idown, jup, jdown, kup, kdown !- up and down cpu number -!
  integer :: icputable(-1:iprocs,-1:jprocs,-1:kprocs) !- cpu table -!
  real(8) :: dt, dtshot, dt1, dt2, dt3, dtcfl, dtcflg !- time steps -!
  real(8) :: dx1a, dx2a, dx3a, akap1a, akap2a, akap3a 
  real(8) :: time, timeh, tnext, timesr, timeo, time1 !- time -!
  integer :: istop, ierror                            !- error flags-!
  integer :: it0, it, ihi, its, itsr, ihisr, isr, iit !- iterations -!
  integer :: ig, nd
  real(8) :: cputime !- cputime -!
  integer :: nm0, nm1, nm2, ider
  integer :: npe, myrank, merr, myranki, myrankj, myrankk !- MPI cpu number -!
!
!======================================================================@
!- for MPI setup -!
  call mpi_init(merr)
  call mpi_comm_size(mpi_comm_world,npe   ,merr)
  call mpi_comm_rank(mpi_comm_world,myrank,merr)
!
!- check consistency of total number of cpu -!
!
  if(npe .ne. iprocs*jprocs*kprocs) then
    if(myrank .eq. 0) then
      write(*,*) 'Error: for ncpu \= iprocs*kprocs'
    endif
    call mpi_finalize(merr)
    stop
  endif 
!
!- check parameter range for grids -!
!
  call makecputable(myrank,icputable,myranki,myrankj,myrankk)
  call pararange(1,imax,iprocs,myranki,is,ie)
  call pararange(1,jmax,jprocs,myrankj,js,je)
  call pararange(1,kmax,kprocs,myrankk,ks,ke)
!
  if(myrank .eq. 0) then
    write(*,*) 'npe,iprocs,kprocs=', npe, iprocs, jprocs, kprocs
  endif
  write(*,*) 'myrank,myranki,myrankj,myrankk=', myrank, myranki, &
             myrankj, myrankk
!      
!- Buffer grid number -!
!
  nm0=1
  nm1=1
!- from reconstruction schemes -!      
  if(irec .eq. 1 .or. irec .eq. 2 .or. irec .eq. 3 .or. &
     irec .eq. 11 .or. irec .eq. 12 .or. irec .eq. 13) then
    nm0=2
  elseif(irec .eq. 4 .or. irec .eq. 6 .or. irec .eq. 7 .or. &
         irec .eq. 8 .or. irec .eq. 9 .or. irec .eq. 10) then
    nm0=3
  elseif(irec .eq. 5) then
    nm0=4
  endif
!- from CT schemes -!      
  if(ict .eq. 0) then
    nm1=nm0
  elseif(ict .eq. 1) then
    nm1=nm0+1
  elseif(ict .eq. 2) then
    nm1=nm0+2
  elseif(ict .eq. 3) then
    nm1=nm0+3
  endif
!- from interporations (not used)-!      
  if(irec .le. 3) then
    ider=0
    nm2=nm1
  elseif(irec .eq. 4 .or. irec .eq. 6 .or. irec .eq. 7 .or. irec .eq. 8) then
    ider=1
    nm2=nm1+1
  elseif(irec .eq. 5) then
    ider=2
    nm2=nm1+2
  endif
!- adjust calculation grid number with buffer for each cpus -!
!- i-direction -!
  if(myranki .eq. 0) then
    is1=1
  else
    is1=is-nm1
  endif
  if(myranki .eq. iprocs-1) then
    ie1=ie
  else
    ie1=ie+nm1
  endif
!- j-direction -!
  if(myrankj .eq. 0) then
    js1=1
  else
    js1=js-nm1
  endif
  if(myrankj .eq. jprocs-1) then
    je1=je
  else
    je1=je+nm1
  endif
!- k-direction -!
  if(myrankk .eq. 0) then
    ks1=1
  else
    ks1=ks-nm1
  endif
  if(myrankk .eq. kprocs-1) then
    ke1=ke
  else
    ke1=ke+nm1
  endif
!- check adjusted grid number is smaller than minimum grid number or not -!
  if(ie1-is1+1 .le. 2*nm1 .or. je1-js1+1 .le. 2*nm1 &
      .or. ke1-ks1+1 .le. 2*nm1) then
    write(*,*) 'Error: for grid size < minimum number of grid'
    write(*,*) 'myrank, ith grid #, jth grid #, kth grid # =',myrank, &
               ie1-is1+1, je1-js1+1, ke1-ks1+1, 2*nm1 
    call mpi_finalize(merr)
    stop        
  endif
!
  write(*,*) 'is,ie,js,je,ks,ke=',is,ie,js,je,ks,ke 
!
!======================================================================@
!-  Allocate variables -!
!======================================================================@
!
  allocate( uu(nv,is1:ie1,js1:je1,ks1:ke1), &
            uo(nv,is1:ie1,js1:je1,ks1:ke1), &
            u0(nv,is1:ie1,js1:je1,ks1:ke1), &
            uh(nv,is1:ie1,js1:je1,ks1:ke1), &
            us(nv,is1:ie1,js1:je1,ks1:ke1), &
            un(nv,is1:ie1,js1:je1,ks1:ke1), &
            uupr(nv,is1:ie1,js1:je1,ks1:ke1), &
            uusr(nv,is1:ie1,js1:je1,ks1:ke1), &
            uuir(nv,is1:ie1,js1:je1,ks1:ke1), &
            uuil(nv,is1:ie1,js1:je1,ks1:ke1), &
            uujr(nv,is1:ie1,js1:je1,ks1:ke1), &
            uujl(nv,is1:ie1,js1:je1,ks1:ke1), &
            uukr(nv,is1:ie1,js1:je1,ks1:ke1), &
            uukl(nv,is1:ie1,js1:je1,ks1:ke1), stat=merr )

  allocate( uri(nv,is1:ie1,js1:je1,ks1:ke1), &
            urio(nv,is1:ie1,js1:je1,ks1:ke1), &
            uri0(nv,is1:ie1,js1:je1,ks1:ke1), &
            urip(nv,is1:ie1,js1:je1,ks1:ke1), &
            urisr(nv,is1:ie1,js1:je1,ks1:ke1), &
            uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
            uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
            urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
            urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
            urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
            urikl(nv,is1:ie1,js1:je1,ks1:ke1), &
            uria(nv,is1:ie1,js1:je1,ks1:ke1),  &
            urib(nv+1,is1:ie1,js1:je1,ks1:ke1), stat=merr )
  
  allocate( gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
            gconi(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
            gconj(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
            gconk(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
            gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
            gcovi(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
            gcovj(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
            gcovk(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
            detg(is1:ie1,js1:je1,ks1:ke1), &
            detgi(is1:ie1,js1:je1,ks1:ke1), &
            detgj(is1:ie1,js1:je1,ks1:ke1), &
            detgk(is1:ie1,js1:je1,ks1:ke1), stat=merr)

  allocate( christ(0:3,0:3,0:3,is1:ie1,js1:je1,ks1:ke1), stat=merr)

  allocate( ww(3,nv,is1:ie1,js1:je1,ks1:ke1), &
            wwo(3,nv,is1:ie1,js1:je1,ks1:ke1), &
            wwir(nv,is1:ie1,js1:je1,ks1:ke1), &
            wwil(nv,is1:ie1,js1:je1,ks1:ke1), &
            wwjr(nv,is1:ie1,js1:je1,ks1:ke1), &
            wwjl(nv,is1:ie1,js1:je1,ks1:ke1), &
            wwkr(nv,is1:ie1,js1:je1,ks1:ke1), &
            wwkl(nv,is1:ie1,js1:je1,ks1:ke1), stat=merr )

  allocate( sf(2:5,is1:ie1,js1:je1,ks1:ke1), &
            tenr(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), stat=merr )

  allocate( cmaxi(is1:ie1,js1:je1,ks1:ke1), &
            cmini(is1:ie1,js1:je1,ks1:ke1), &
            cmaxj(is1:ie1,js1:je1,ks1:ke1), &
            cminj(is1:ie1,js1:je1,ks1:ke1), &
            cmaxk(is1:ie1,js1:je1,ks1:ke1), &
            cmink(is1:ie1,js1:je1,ks1:ke1), stat=merr )

  allocate( b1(is1:ie1,js1:je1,ks1:ke1), b2(is1:ie1,js1:je1,ks1:ke1), &
            b3(is1:ie1,js1:je1,ks1:ke1), &
            util1(is1:ie1,js1:je1,ks1:ke1), util2(is1:ie1,js1:je1,ks1:ke1), &
            util3(is1:ie1,js1:je1,ks1:ke1), &
            v1(is1:ie1,js1:je1,ks1:ke1), v2(is1:ie1,js1:je1,ks1:ke1), &
            v3(is1:ie1,js1:je1,ks1:ke1), &
            de(is1:ie1,js1:je1,ks1:ke1), pr(is1:ie1,js1:je1,ks1:ke1), &
            stat=merr )
!
!======================================================================@
!-    file open for restart -!
!======================================================================@
!- file is outputed from each cpus "restart000.outdat" -!
!    open( unit=1, file='ipram', status='old')
  write(filename1,991) myrank
  open( unit=9, file=filename1, status='unknown', &
        form='unformatted')
!
 991 format ('restart',i3.3,'.outdat')
!
!======================================================================@
!- Initialization of varibles -!
!======================================================================@
!
!- Fixed parameter -!
!
  istop=0
  ierror=0
!
  nd=0
!
  dt=0.01d0
  dtshot=tmax/float(nshot)
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
        sf(2,i,j,k)=0.d0
        sf(3,i,j,k)=0.d0
        sf(4,i,j,k)=0.d0
        sf(5,i,j,k)=0.d0
      enddo
    enddo
  enddo
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
        do m=1,nv
           ww(1,m,i,j,k)=0.d0
           ww(2,m,i,j,k)=0.d0
           ww(3,m,i,j,k)=0.d0
          
           uri(m,i,j,k)=0.d0
           urip(m,i,j,k)=0.d0
           uriir(m,i,j,k)=0.d0
           urijr(m,i,j,k)=0.d0
           urikr(m,i,j,k)=0.d0
           uriil(m,i,j,k)=0.d0
           urijl(m,i,j,k)=0.d0
           urikl(m,i,j,k)=0.d0
         
           uupr(m,i,j,k)=0.d0
           uuir(m,i,j,k)=0.d0
           uujr(m,i,j,k)=0.d0
           uukr(m,i,j,k)=0.d0
           uuil(m,i,j,k)=0.d0
           uujl(m,i,j,k)=0.d0
           uukl(m,i,j,k)=0.d0
         
           wwir(m,i,j,k)=0.d0
           wwjr(m,i,j,k)=0.d0
           wwkr(m,i,j,k)=0.d0
           wwil(m,i,j,k)=0.d0
           wwjl(m,i,j,k)=0.d0
           wwkl(m,i,j,k)=0.d0         
        enddo
      enddo
    enddo
  enddo
!
!======================================================================@
!- Print list of parameters -!
!======================================================================@
!
  if(myrank .eq. 0) then
    write(6,*) 'job name: ',jname,',  date: ',date
    write(6,*) 'Restart : ',icres
    write(6,*) '         List of Parameters '
    write(6,*) 'dtime =',dt,',  tmax =', tmax
    write(6,*) 'itmax =',itmax
    write(6,*) 'imax  =',imax,',  jmax =',jmax,',  kmax =',kmax
    write(6,*) 'npe =',npe
    write(6,*) 'xmin  =',xmin,',  xmax =',xmax
    write(6,*) 'ymin  =',ymin,',  ymax =',ymax
    write(6,*) 'zmin  =',zmin,',  zmax =',zmax
    write(6,*) 'ix1 =',ix1,',  ix2 =',ix2,',  ix3=',ix3
    write(6,*) 'model =',model,',  metric =',metric
    write(6,*) 'ieos =',ieos
    write(6,*) 'c0 =',c0
    write(6,*) 'gamma =',gam,',  gam0 =',gam0
    write(6,*) 'akm =',akm
    write(6,*) 'irec =',irec,', icha =',icha
    write(6,*) 'ihll =',ihll,',  ict = ',ict
    write(6,*) 'iwvec =',iwvec
  endif
!
!======================================================================@
!- Grid generation -!
!======================================================================@
!
  call grid(x1,x2,x3,dx1,dx2,dx3,dx1a,dx2a,dx3a, &
            dx1b,dx2b,dx3b,x1a,x2a,x3a,x1b,x2b,x3b,x1ab,x2ab,x3ab)
!
!======================================================================@
!- Initial setting for Metric (coordinates) -!
!======================================================================@
!
  call coord1(gcov,gcovi,gcovj,gcovk,gcon,gconi,gconj,gconk, &
              detg,detgi,detgj,detgk,x1,x2,x3,x1a,x2a,&
              x3a,x1b,x2b,x3b,is1,ie1,js1,je1,ks1,ke1)
!
!- calculation of christfel symbols -!
  call calchrist(christ,x1,x2,x3,x1b,x2b,x3b,is1,ie1,js1,je1,ks1,ke1)
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
!- Researrt simulation or Not -!
!
  if( icres .eq. 0 ) then !if icres = 0 => not restarted simulation
!
!- Initial Model Setup
!
    if(model .eq. 1) then !-1D Test -!
      call md1dtest(uri,is1,ie1,js1,je1,ks1,ke1)
    elseif(model .eq. 2) then !- 1D shock-tube Test -!
      call mdshoctub(uri,is1,ie1,js1,je1,ks1,ke1)
    elseif(model .eq. 3) then !- 1D accretion Test -!
      call md1dacc(uri,is1,ie1,js1,je1,ks1,ke1)
    elseif(model .eq. 4) then
      call mdeqcor1(uri,x1,is1,ie1,js1,je1,ks1,ke1)
    elseif(model .eq. 5) then
      call mdeqcor2(uri,x1,is1,ie1,js1,je1,ks1,ke1)
    elseif(model .eq. 6) then !- Bondi accretion Test -!
      call mdffcor1(uri,gcov,gcon,detg,x1,x2,x3,x1b,x2b,x3b,&
                    is1,ie1,js1,je1,ks1,ke1)
    elseif(model .eq. 7) then !- Fishbone-Moncrief Torus model -!
      call mdtorus1(uri,gcov,gcon,detg,x1,x2,x3,x1b,x2b,x3b,dx1,dx2,dx3, &
                    myrank,is1,ie1,js1,je1,ks1,ke1)
    elseif(model .eq. 8) then !- Power-law Rotating Torus model -!
      call mdtorus2(uri,gcov,gcon,detg,x1,x2,x3,x1b,x2b,x3b,dx1,dx2,dx3, &
                    myrank,is1,ie1,js1,je1,ks1,ke1)
    elseif(model .eq. 9) then !- constant-l Torus model -!
      call mdtorus3(uri,gcov,gcon,detg,x1,x2,x3,x1b,x2b,x3b,dx1,dx2,dx3, &
                    myrank,is1,ie1,js1,je1,ks1,ke1)
    elseif(model .eq. 12) then !- constant-l magnetized Torus model -!
      call mdtorus4(uri,gcov,gcon,x1,x2,x3,x1b,x2b,x3b,dx1,dx2,dx3, &
                    is1,ie1,js1,je1,ks1,ke1)
    elseif(model .eq. 14) then !- tilted constant-l Torus model -!
      call mdtildisk1(uri,gcov,gcon,detg,x1,x2,x3,x1b,x2b,x3b,dx1,dx2,dx3, &
                      myrank,is1,ie1,js1,je1,ks1,ke1)
    elseif(model .eq. 15) then !- tilted constant-l Torus model v2-!
      call mdtildisk2(uri,gcov,gcon,detg,x1,x2,x3,x1b,x2b,x3b,dx1,dx2,dx3, &
                      myrank,is1,ie1,js1,je1,ks1,ke1)  
    elseif(model .eq. 10) then !- Recoiling BH model -!
      call mdrecoilbh(uri,gcov,gcon,detg,x1,x2,x3,x1b,x2b,x3b,dx1,dx2,dx3, &
           myrank,is1,ie1,js1,je1,ks1,ke1)
!      call mdrecoilbh1a(uri,gcov,gcon,detg,x1,x2,x3,x1b,x2b,x3b,dx1,dx2,dx3, &
!                    myrank,is1,ie1,js1,je1,ks1,ke1)
    elseif(model .eq. 13) then !- Recoiling BH in magnetized Torus model -!
      call mdrecoilbh2(uri,gcov,gcon,x1,x2,x3,x1b,x2b,x3b,dx1,dx2,dx3, &
                       is1,ie1,js1,je1,ks1,ke1)
    elseif(model .eq. 20) then !- monopole field in BH magnetoshere -!
      call mdmagsph1(uri,gcov,gcon,detg,x1,x2,x3,is1,ie1,js1,je1,ks1,ke1)  
    elseif(model .eq. 31) then !- wind-wind collision model -!
      call mdwindcol(uri,gcov,gcon,x1,x2,x3,myrank,is1,ie1,js1,je1,ks1,ke1)
    else
      stop
    endif
!
!- Primitive Variables Set (uri) -!
!
!    do k=ks1,ke1
!      do j=js1,je1
!        do i=is1,ie1
!
!          uri(1,i,j,k)=de(i,j,k) !- density -!
!          uri(2,i,j,k)=util1(i,j,k) !- 4-velocity in i-component -!
!          uri(3,i,j,k)=util2(i,j,k) !- 4-velocity in j-component -!
!          uri(4,i,j,k)=util3(i,j,k) !- 4-velocity in k-component -!
!!          uri(2,i,j,k)=v1(i,j,k) !- 3-velocity in i-component -!
!!          uri(3,i,j,k)=v2(i,j,k) !- 3-velocity in j-component -!
!!          uri(4,i,j,k)=v3(i,j,k) !- 3-velocity in k-component -!
!          uri(5,i,j,k)=pr(i,j,k) ! gas pressure -!
!          uri(6,i,j,k)=0.d0 !- tracer (not used) -!
!          uri(7,i,j,k)=b1(i,j,k) !- magnetic field in i-component -!
!          uri(8,i,j,k)=b2(i,j,k) !- magnetic field in j-component -!
!          uri(9,i,j,k)=b3(i,j,k) !- magnetic field in k-component -!
!
!        enddo
!      enddo
!    enddo
!
!- calcualtion of Conserved Variables Set (uu) -!
!
    call caluu1(uri,uu,gcov,gcon,detg,x1,is1,ie1,js1,je1,ks1,ke1)
!
!- copy of variables as initial value -!
!
    call ident(uu,uo,is1,ie1,js1,je1,ks1,ke1)
    call ident(uu,u0,is1,ie1,js1,je1,ks1,ke1)
    call ident(uu,uusr,is1,ie1,js1,je1,ks1,ke1)
    call ident(uri,urio,is1,ie1,js1,je1,ks1,ke1)
    call ident(uri,uri0,is1,ie1,js1,je1,ks1,ke1)
    call ident(uri,urisr,is1,ie1,js1,je1,ks1,ke1)
!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
!- Initialize for time advance loop start -!
!
    it0=0
    ihi=0
    time=0.d0
    timeh=0.d0
    tnext=0.d0
  endif
!
!- come back point if icres = 1 (for restart) -!
!
!======================================================================@
!- Startup for restart simulation -!
!======================================================================@
!
!-  If restart is 'on' skip the initial setup and move this line -!
!
  if( icres .ne. 0 ) then
!    call restar(uu,u0,uri,uri0,it0,ihi,timeh,ks1,ke1)
    call restar2(uu,u0,uri,uri0,it0,ihi,timeh,nd,npe,myrank,&
                 is1,ie1,js1,je1,ks1,ke1)
!
!- next stored time -! 
!!      tnext=timeh+dtshot
    tnext=(int(timeh/dtshot)+1)*dtshot
!- copy of variables -!
    call ident(uu,uusr,is1,ie1,js1,je1,ks1,ke1)
    call ident(uri,urisr,is1,ie1,js1,je1,ks1,ke1)
  endif
!
!**********************************************************************@
!- !! Time Advance Step Loop Start !! -!
!**********************************************************************@
! 
  its=(itmax-it0+1)/nshot
  if( its .lt. 1 ) then
    its = 1
  endif
  itsr=it0
  ihisr=ihi
  isr=0
  timesr=timeh
!
  do it=it0,itmax !- iteration for numerical step -!
!
!-----------------------------------------------------------------------
!- CPU Time Check -!
!-----------------------------------------------------------------------
!
!    call kclock(ig)
    if(myrank .eq. 0) then
      call cpu_time(cputime)
    endif
!
    call mpi_barrier(mpi_comm_world,merr)
    call mpi_bcast(cputime,1,mpi_double_precision,0,mpi_comm_world,merr)
!
!    if( ig .ge. icpu ) then
    if(cputime .ge. float(icpu)) then
      if(myrank .eq. 0) then
        write(6,*) ' CPU Time Over : cputime =',cputime,'(sec)'
       endif
       istop=1
    endif
!
    iit=it-it0
!
!x   if( mod(it,its) .eq. 0  .or. istop.eq.1 .or. it.eq.itmax ) then
!
!
!----------------------------------------------------------------------
!- Simulation time check -!
!----------------------------------------------------------------------
!
    if( dtmin.ge.0.d0 ) then 
      timeo=time
      time=timeh
    else
      time = dt*float(it)
    endif
!
    if( time.gt.tmax ) then !- simulation time > final simulation time -!
      if(myrank .eq. 0) then
        write(6,*) '  '
        write(6,*) 'Congraturations! Run Succeeded! : time =' &
                    ,time,', tmax',tmax
        write(6,*) '  '
      endif
      istop = 1
    endif
!
!----------------------------------------------------------------------
!- Data store and output -!
!----------------------------------------------------------------------
!
    call mpi_barrier(mpi_comm_world,merr)
!
!- Data output criterion: 1. simulation time > output simulation time    -!
!-                        2. time step number > maximum time step number -!
!-                        3. simulation has numerical error or false signal -!
! 
    if( time.ge.tnext .or. it.eq.itmax .or. & 
        ierror.ne.0 .or. istop.eq.1 ) then
!
      if( ierror.eq.1 ) then
        write(6,*) 'stop by fault; time =',time
      endif
!
      ihi=ihi+it-it0+1
!- data stored number -!
      if(time .eq. 0.0d0) then
       nd=0
      else
       nd=nd+1
      endif
!
!- Data output for restart simulation -!
!
      if( (isskip.eq.0 .or. istop.eq.1) .and. ierror.eq.0 ) then
        call store1(uu,u0,uri,uri0,it,ihi,time,nd,npe,myrank,&
                    is1,ie1,js1,je1,ks1,ke1)
      elseif( it.eq.itmax .and. ierror.eq.0 ) then
        call store1(uu,u0,uri,uri0,it,ihi,time,nd,npe,myrank,&
                    is1,ie1,js1,je1,ks1,ke1)
      elseif( mod(isr,isskip).eq.0 .and. ierror.eq.0 ) then
        call store1(uusr,u0,urisr,uri0,itsr,ihisr,timesr,nd,npe,myrank, &
                    is1,ie1,js1,je1,ks1,ke1)
!
!- copy of output time information -!
!
        itsr=it
        ihisr=ihi
        timesr=time
        call ident(uu,uusr,is1,ie1,js1,je1,ks1,ke1)
        call ident(uri,urisr,is1,ie1,js1,je1,ks1,ke1)
      endif
!
!- Data output of Simulation results -!
!
      if(myrank .eq. 0) then
        write(6,*) ' data output for analysis'
        write(6,*) 'time=',time
!        write(6,*) 'dtshot=',dtshot
!        write(6,*) 'tnext=',tnext
        write(6,*) 'cputime=',cputime
      endif
!
!- convert 4-velocity ?-!
      if(metric .eq. 103 .or. metric .eq. 203) then !- Kerr BL -!
!- convert 4-velocity to 3-velocity
!        call cal4to3vel(uri,uria,gcov,is1,ie1,js1,je1,ks1,ke1)
!- no convertion: 4-velocity spacetial component u_tilda 
        call ident(uri,uria,is1,ie1,js1,je1,ks1,ke1)         
!- convert to 4-velocity (for radiation)
        call cal4to4vel(uri,urib,gcov,gcon,is1,ie1,js1,je1,ks1,ke1)
!        
      elseif(metric .eq. 303) then !- Kerr KS -!
!- convert from 4-vel (tilda) in KS to 4-vel (tilda) in BL
!         call transks2bl(uri,uria,gcov,gcon,x1,x2,x3,is1,ie1,js1,je1,ks1,ke1)
!- convert to 4-velocity (for radiation)
        call cal4to4vel(uri,urib,gcov,gcon,is1,ie1,js1,je1,ks1,ke1)
!- no convertion: 4-velocity spacetial component u_tilda 
        call ident(uri,uria,is1,ie1,js1,je1,ks1,ke1)           
      else
        call ident(uri,uria,is1,ie1,js1,je1,ks1,ke1)
      endif
!    
!- Output for simulation data
      if(metric .eq. 403) then
        call output1(uria,detg,x1b,x2b,x3b,time,nm1,nd,npe,myrank,&
                     myranki,myrankj,myrankk,is1,ie1,js1,je1,ks1,ke1)
      else
        call output1(uria,detg,x1,x2,x3,time,nm1,nd,npe,myrank,&
                     myranki,myrankj,myrankk,is1,ie1,js1,je1,ks1,ke1)
      endif  
!- Output for radiation
!      call output2(urib,x1,x2,x3,time,nm1,nd,npe,myrank,myranki,myrankj,&
!                   myrankk,is1,ie1,js1,je1,ks1,ke1)
!- output 4-velocity -!
!      call output1(uri,detg,x1,x2,x3,time,nm1,nd,npe,myrank,myranki,myrankj,&
!                   myrankk,is1,ie1,js1,je1,ks1,ke1)
!
!- next output time -!
      isr=isr+1
      tnext=tnext+dtshot
!
    endif
!
    if( ierror .ne. 0 .or. istop .eq. 1 .or. it .eq. itmax)  then
      goto 100
    endif
!
!=======================================================================
!- Start next time step calculation -!
!=======================================================================
!
!-----------------------------------------------------------------------
!- Copy of variables (before 1st step)
!-----------------------------------------------------------------------
!
    call ident(uu,uo,is1,ie1,js1,je1,ks1,ke1)
    call ident(uri,urio,is1,ie1,js1,je1,ks1,ke1)
!
!!======================================================================@
!- Calcuration of dt for next time step -!
!======================================================================@
!
    if( dtmin.ge.0.d0 ) then
!      if(metric .eq. 403) then
        call cdtcfl4(uri,gcov,gcon,dtcfl,dx1b,dx2b,dx3b,dt1,dt2,dt3, &
                     nm1,is1,ie1,js1,je1,ks1,ke1)
!      else
!        call cdtcfl4(uri,gcov,gcon,dtcfl,dx1,dx2,dx3,dt1,dt2,dt3, &
!                     nm1,is1,ie1,js1,je1,ks1,ke1)
!      endif
!
!- MPI
!
      call mpi_allreduce(dtcfl,dtcflg,1,mpi_double_precision,mpi_min, &
                         mpi_comm_world,merr)
      dt=cfl*dtcflg
!- dt/dx in i-direction -!
      do i=1,imax-1
        akap1(i)=dt/dx1(i)
        akap1b(i)=dt/dx1b(i)
      enddo
      akap1a=dt/dx1a
!- dt/dx in j-direction -! 
      do j=1,jmax-1
        akap2(j)=dt/dx2(j)
        akap2b(j)=dt/dx2b(j)
      enddo
      akap2a=dt/dx2a
!- dt/dx in k-direction -!
      do k=1,kmax-1
        akap3(k)=dt/dx3(k)
        akap3b(k)=dt/dx3b(k)
      enddo
      akap3a=dt/dx3a
!
      timeh=time+dt
    else
      timeh= dt*(float(it)+1.0)
    endif
       
    if( dt.lt.dtmin ) then
!      write(6,*) ' dt < dtmin: dt =',dt,', dtmin =',dtmin
      istop = 1
    endif
!
!-----------------------------------------------------------------------
!
    if(myrank .eq. 0) then
      write(6,620) 'Simulation running: it = ',it, &
                 ',  dt =',dt,', time=',time
620   format(1h ,a27,i6, 2(a7,1pe11.3))
    endif
!
!**********************************************************************@
!- First Step -!
!**********************************************************************@
!
!- Reconstraction Step -!
!
    call rec2(uri,x1,x2,x3,x1a,x2a,x3a,dx1,dx2,dx3,dx1b,dx2b,dx3b, &
              x1b,x2b,x3b,x1ab,x2ab,x3ab, & 
              uriir,urijr,urikr,uriil,urijl,urikl, &
              uuir,uujr,uukr,uuil,uujl,uukl, &
              wwir,wwjr,wwkr,wwil,wwjl,wwkl, &
              gcovi,gcovj,gcovk,gconi,gconj,gconk, &
              detgi,detgj,detgk,nm0,is1,ie1,js1,je1,ks1,ke1)
!
!- Calculation of Characteristics -!
!
    call calcha(uriir,urijr,urikr,uriil,urijl,urikl, &
                gcovi,gcovj,gcovk,gconi,gconj,gconk, &
                cmaxi,cmini,cmaxj,cminj,cmaxk,cmink,nm0, &
                is1,ie1,js1,je1,ks1,ke1)
!
!- Calculation of Numerical flux -!
!
    call hll(ww,uuir,uujr,uukr,uuil,uujl,uukl, &
             wwir,wwjr,wwkr,wwil,wwjl,wwkl, &
             cmaxi,cmini,cmaxj,cminj,cmaxk,cmink,nm0, &
             is1,ie1,js1,je1,ks1,ke1)
!
!- Constrained Transport -!
!
    call ct(ww,wwo,uri,nm0,nm1,is1,ie1,js1,je1,ks1,ke1)
!
!    call deriv(ww,ider,nm1,nm2,is1,ie1,js1,je1,ks1,ke1)
!
!- Calcuration of Source term -!
!
    call caltenr(uri,tenr,gcov,gcon,is1,ie1,js1,je1,ks1,ke1)
    call calsf(sf,tenr,christ,detg,is1,ie1,js1,je1,ks1,ke1)
!
!- Time advance of first step -!
!
    call rk2fst(uh,ww,uu,nm1,akap1b,akap2b,akap3b,is1,ie1,js1,je1,ks1,ke1)
!
!- Add Source term -!
!
    call rk2adsff(uh,sf,dt,nm1,is1,ie1,js1,je1,ks1,ke1)
!
!- Boundary condition for conserved variables -!
!
!y   call bnd4(uh,uo,x1,nm1,is1,ie1,js1,je1,ks1,ke1, &
!y             myranki,myrankj,myrankk,icputable)
!
!- Recovery step -!
!
    call recov(uh,uri,gcov,gcon,detg,x1,x2,x3,nm1, &
               is1,ie1,js1,je1,ks1,ke1)      
!    
!- Boundary condition for primitive variables -!
!
    call bnd4(uri,urio,x1,nm1,is1,ie1,js1,je1,ks1,ke1, &
              myranki,myrankj,myrankk,icputable) 
!
!- Inflow/outflow check from radial inner/outer boundary -!
!   
    if(model .ge. 7 .and. model .le. 20) then
       call infchk(uri,gcov,gcon,nm1,is1,ie1,js1,je1,ks1,ke1, &
                  myranki)
    endif     
!
!- Check and Correction (Atmosphere treatment) -!
!
    if(model .ge. 7 .and. model .le. 15) then
      call atmosphere1a(uh,uri,gcov,gcon,detg,x1,is1,ie1,js1,je1,ks1,ke1)
    endif    
!
!- Check and Correction (pressure and density) -!
!
    call pminmax2b(uh,uri,gcov,gcon,detg,x1,is1,ie1,js1,je1,ks1,ke1)
!
!
!- Calculation of conserved variables at boundary region -!
!
!    call calconv(uh,uri,gcov,gcon,detg,x1,nm1, &
!                 is1,ie1,js1,je1,ks1,ke1,myranki,myrankj,myrankk)
!
!- Artificial Damping (primitive variables) -!
!
!    if(adamp.gt.0.d0 .and. rdamp.gt.0.d0) then
!      call damp4(uh,uri,uri0,gcov,gcon,detg, &
!                 x1,x3,is1,ie1,js1,je1,ks1,ke1)
!    endif
!    if(adamp.gt.0.d0) then
!      call damp4a(uh,uri,uri0,gcov,gcon,detg, &
!                  x1,x3,is1,ie1,js1,je1,ks1,ke1)
!    endif
!    
!- calcualtion of Conserved Variables Set (uu) -!
!
    call caluu1(uri,uh,gcov,gcon,detg,x1,is1,ie1,js1,je1,ks1,ke1)
!
!- MPI data exchange -!                                                        
!                                                                              
    call mpi_barrier(mpi_comm_world,merr)
    call mpiex(uri,nm1,is1,ie1,js1,je1,ks1,ke1,myranki,myrankj,&
               myrankk,icputable)
!    call mpiex(uh,nm1,is1,ie1,js1,je1,ks1,ke1,myranki,myrankj,&
!               myrankk,icputable)  
!
!- Copy of variables (uh => uu)-!
!
    call ident(uh,uu,is1,ie1,js1,je1,ks1,ke1)
!
!**********************************************************************@
!- Second Step -!
!**********************************************************************@
!
!- Reconstraction Step -!
!
    call rec2(uri,x1,x2,x3,x1a,x2a,x3a,dx1,dx2,dx3,dx1b,dx2b,dx3b, &
              x1b,x2b,x3b,x1ab,x2ab,x3ab, & 
              uriir,urijr,urikr,uriil,urijl,urikl, &
              uuir,uujr,uukr,uuil,uujl,uukl, &
              wwir,wwjr,wwkr,wwil,wwjl,wwkl, &
              gcovi,gcovj,gcovk,gconi,gconj,gconk, &
              detgi,detgj,detgk,nm0,is1,ie1,js1,je1,ks1,ke1)
!
!- Calculation of Characteristics -!
!
    call calcha(uriir,urijr,urikr,uriil,urijl,urikl, &
                gcovi,gcovj,gcovk,gconi,gconj,gconk, &
                cmaxi,cmini,cmaxj,cminj,cmaxk,cmink,nm0, &
                is1,ie1,js1,je1,ks1,ke1)
!
!- Calculation of Numerical flux -!
!
    call hll(ww,uuir,uujr,uukr,uuil,uujl,uukl, &
             wwir,wwjr,wwkr,wwil,wwjl,wwkl, &
             cmaxi,cmini,cmaxj,cminj,cmaxk,cmink,nm0, &
             is1,ie1,js1,je1,ks1,ke1)
!
!- Constrained Transport -!
!
    call ct(ww,wwo,uri,nm0,nm1,is1,ie1,js1,je1,ks1,ke1)
!
!    call deriv(ww,ider,nm1,nm2,is1,ie1,js1,je1,ks1,ke1)
!
!- Calcuration of source term -!
!
    call caltenr(uri,tenr,gcov,gcon,is1,ie1,js1,je1,ks1,ke1)
    call calsf(sf,tenr,christ,detg,is1,ie1,js1,je1,ks1,ke1)
!
!- Time advance of second step -!
!
    if(irkt .eq. 2) then !for 2nd-order RK scheme
      call rk2snd(us,ww,uo,uh,nm1,akap1b,akap2b,akap3b,&
                  is1,ie1,js1,je1,ks1,ke1)
    elseif(irkt .eq. 3) then !for 3rd-order RK scheme
      call rk3snd(us,ww,uo,uh,nm1,akap1b,akap2b,akap3b,&
                  is1,ie1,js1,je1,ks1,ke1)
    endif
!
!- Add the Source Term -!
!
    if(irkt .eq. 2) then
      call rk2adsfs(us,sf,dt,nm1,is1,ie1,js1,je1,ks1,ke1)
    elseif(irkt .eq. 3) then
      call rk3adsfs(us,sf,dt,nm1,is1,ie1,js1,je1,ks1,ke1)
    endif
!
!- Boundary Condition for conserved variables -!
!
!y   call bnd4(us,uo,x1,nm1,is1,ie1,js1,je1,ks1,ke1,&
!y             myranki,myrankj,myrankk,icputable)
!
!- Recovery step -!
!
    call recov(us,uri,gcov,gcon,detg,x1,x2,x3,nm1, &
               is1,ie1,js1,je1,ks1,ke1)  
!
!- Boundary Condition for primitive variables -!
!
    call bnd4(uri,urio,x1,nm1,is1,ie1,js1,je1,ks1,ke1,&
              myranki,myrankj,myrankk,icputable)
!
!- Inflow/outflow check from radial inner/outer boundary -!
!
    if(model .ge. 7 .and. model .le. 20) then
       call infchk(uri,gcov,gcon,nm1,is1,ie1,js1,je1,ks1,ke1, &
                 myranki)
    endif    
!
!- Check and Correction (Atmosphere treatment) -!
!
    if(model .ge. 7 .and. model .le. 15) then
      call atmosphere1a(us,uri,gcov,gcon,detg,x1,is1,ie1,js1,je1,ks1,ke1)
    endif     
!
!- Check and Correction (pressure and density) -!
!
    call pminmax2b(us,uri,gcov,gcon,detg,x1,is1,ie1,js1,je1,ks1,ke1)
!
!- Calculation of conserved variables at boundary region -!
!
!    call calconv(us,uri,gcov,gcon,detg,x1,nm1, &
!                 is1,ie1,js1,je1,ks1,ke1,myranki,myrankj,myrankk)
!
!- Artificial Damping (primitive variables) -!
!
!    if(adamp.gt.0.d0 .and. rdamp.gt.0.d0) then
!      call damp4(us,uri,uri0,gcov,gcon,detg, &
!                   x1,x3,is1,ie1,js1,je1,ks1,ke1)
!    endif
!    if(adamp.gt.0.d0) then
!      call damp4a(us,uri,uri0,gcov,gcon,detg, &
!                    x1,x3,is1,ie1,js1,je1,ks1,ke1)
!    endif
!    
!- calcualtion of Conserved Variables Set (uu) -!
!
    call caluu1(uri,us,gcov,gcon,detg,x1,is1,ie1,js1,je1,ks1,ke1)
!                                                                              
!- MPI data exchange -!                                                        
!                                                                              
    call mpi_barrier(mpi_comm_world,merr)
    call mpiex(uri,nm1,is1,ie1,js1,je1,ks1,ke1,myranki,myrankj,&
               myrankk,icputable)
!    call mpiex(us,nm1,is1,ie1,js1,je1,ks1,ke1,myranki,myrankj,& 
!               myrankk,icputable)    
!
!- Copy of variables (us => uu)-!
!
    call ident(us,uu,is1,ie1,js1,je1,ks1,ke1)
!
!**********************************************************************@
!- Additional calculation for third Step of time evolution
!**********************************************************************@
!ccc
    if(irkt .eq. 3) then
!ccc
!
!- Reconstraction Step -!
!
      call rec2(uri,x1,x2,x3,x1a,x2a,x3a,dx1,dx2,dx3,dx1b,dx2b,dx3b, &
                x1b,x2b,x3b,x1ab,x2ab,x3ab, &       
                uriir,urijr,urikr,uriil,urijl,urikl, &
                uuir,uujr,uukr,uuil,uujl,uukl, &
                wwir,wwjr,wwkr,wwil,wwjl,wwkl, &
                gcovi,gcovj,gcovk,gconi,gconj,gconk, &
                detgi,detgj,detgk,nm0,is1,ie1,js1,je1,ks1,ke1)
!
!- Calculation of Characteristics -!
!
      call calcha(uriir,urijr,urikr,uriil,urijl,urikl, &
                  gcovi,gcovj,gcovk,gconi,gconj,gconk, &
                  cmaxi,cmini,cmaxj,cminj,cmaxk,cmink,nm0, &
                  is1,ie1,js1,je1,ks1,ke1)
!
!- Calculation of Numerical flux -!
!
      call hll(ww,uuir,uujr,uukr,uuil,uujl,uukl, &
               wwir,wwjr,wwkr,wwil,wwjl,wwkl, &
               cmaxi,cmini,cmaxj,cminj,cmaxk,cmink,nm0, &
               is1,ie1,js1,je1,ks1,ke1)
!
!- Constrained Transport -!
!
      call ct(ww,wwo,uri,nm0,nm1,is1,ie1,js1,je1,ks1,ke1)
!
!      call deriv(ww,ider,nm1,nm2,is1,ie1,js1,je1,ks1,ke1)
!
!- Calcuration of source term -!
!
      call caltenr(uri,tenr,gcov,gcon,is1,ie1,js1,je1,ks1,ke1)
      call calsf(sf,tenr,christ,detg,is1,ie1,js1,je1,ks1,ke1)
!
!- Time advance of 3rd step -!
!
      call rk3trd(un,ww,uo,us,nm1,akap1b,akap2b,akap3b,&
                  is1,ie1,js1,je1,ks1,ke1)
!
!- Add the Source Term -!
!
      call rk3adsft(un,sf,dt,nm1,is1,ie1,js1,je1,ks1,ke1)
!
!- Boundary Condition for conserved variables -!
!
!y     call bnd4(un,uo,x1,nm1,is1,ie1,js1,je1,ks1,ke1,&
!y               myranki,myrankj,myrankk,icputable)
!
!- Recovery Step -!
!
      call recov(un,uri,gcov,gcon,detg,x1,x2,x3,nm1, &
                 is1,ie1,js1,je1,ks1,ke1) 
!    
!- Boundary Condition for primitive variables -!
!
      call bnd4(uri,urio,x1,nm1,is1,ie1,js1,je1,ks1,ke1, &
                myranki,myrankj,myrankk,icputable)
!
!- Inflow/outflow check from radial inner/outer boundary -!
!
      if(model .ge. 7 .and. model .le. 20 ) then
       call infchk(uri,gcov,gcon,nm1,is1,ie1,js1,je1,ks1,ke1, &
                   myranki)
      endif     
!
!- Check and Correction (Atmosphere treatment) -!
!
      if(model .ge. 7 .and. model .le. 15) then
        call atmosphere1a(un,uri,gcov,gcon,detg,x1,is1,ie1,js1,je1,ks1,ke1)
      endif       
!
!- Check and Correction (density, pressure) -!
!
      call pminmax2b(un,uri,gcov,gcon,detg,x1,is1,ie1,js1,je1,ks1,ke1)
!
!- Calculation of conserved variables at boundary region -!
!
!      call calconv(un,uri,gcov,gcon,detg,x1,nm1, &
!                   is1,ie1,js1,je1,ks1,ke1,myranki,myrankj,myrankk)
!
!- Artificial Damping (primitive variables) -!
!
!      if(adamp.gt.0.d0 .and. rdamp.gt.0.d0) then
!        call damp4(un,uri,uri0,gcov,gcon,detg, &
!                   x1,x3,is1,ie1,js1,je1,ks1,ke1)
!      endif
!      if(adamp.gt.0.d0) then
!        call damp4a(un,uri,uri0,gcov,gcon,detg, &
!                    x1,x3,is1,ie1,js1,je1,ks1,ke1)
!      endif
!    
!- calcualtion of Conserved Variables Set (uu) -!
!
    call caluu1(uri,un,gcov,gcon,detg,x1,is1,ie1,js1,je1,ks1,ke1)
!                                                                              
!- MPI data exchange -!                                      
!                                                              
      call mpi_barrier(mpi_comm_world,merr)
      call mpiex(uri,nm1,is1,ie1,js1,je1,ks1,ke1,myranki,myrankj,&
                 myrankk,icputable)
!      call mpiex(un,nm1,is1,ie1,js1,je1,ks1,ke1,myranki,myrankj,&
!                 myrankk,icputable)  
!
!- Copy of variables (un => uu)-!
!
      call ident(un,uu,is1,ie1,js1,je1,ks1,ke1)
!ccc
    endif
!
  enddo
!
!**********************************************************************@
!- !! Time Step Loop End !! -!
!**********************************************************************@
!
100 continue

  call mpi_barrier(mpi_comm_world,merr)
  call mpi_finalize(merr)

  if(myrank.eq. 0) then
    write(6,*) '== end =='
  endif

  rewind(8)
!
!======================================================================@
!- CPU Time for End of Job -!
!======================================================================@
!
!y  it0=it
!y  call kclock(ig)
!y  write(6,*) 'CPU Time at End of Job :',ig,'sec'
!   .
!======================================================================@
  stop
end program main
