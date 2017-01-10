!***********************************************************************
!                     for I / O 
!***********************************************************************
!--------------------------------------------------------------------
subroutine store1(uu,u0,uri,uri0,it0,ih0,time,nd,npe,myrank,&
                  is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv 
  implicit none
  include 'mpif.h'
!
  integer :: i, j, k, m, n
  integer :: is1,ie1,js1,je1,ks1,ke1
!
  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1), &
             u0(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1), &
             uri0(nv,is1:ie1,js1:je1,ks1:ke1)
!
  integer :: it0, ih0, nd
  real(8) :: time
!
  integer :: npe, myrank, merr
!
!  integer :: istop  
!
!
  if(myrank .eq. 0) then
   write(6,*) ' data store for restart '
  endif
!
!     Check before store for restart                         1998.04.09
!
!  do k=1,kmax
!    do j=1,jmax
!      do i=1,imax
!        if( uu(1,i,j,k).le.0.d0 ) then
!          write(4,*) 'Warning in store for restart: density =< 0'
!          write(4,*) 'density =',uu(1,i,j,k),'at i=',i,',j=',j,',k=',k
!          if( istop.eq.1 ) then
!            write(4,*) 'Return before store'
!            return
!          endif
!         endif
!         if( uu(5,i,j,k).lt.0.d0 ) then
!           write(4,*) 'Warning in store for restart: pressure < 0'
!           write(4,*) 'pressure=',uu(5,i,j,k),'at i=',i,',j=',j,',k=',k
!           if( istop.eq.1 ) then
!             write(4,*) 'Return before store'
!             return
!           endif
!         endif
!       enddo
!     enddo
!   enddo
!
  call mpi_barrier(mpi_comm_world,merr)
!
! data write
  rewind(9)
  write(9) myrank, nd
  write(9) it0,ih0
  write(9) time
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
        write(9) uu(1,i,j,k),uu(2,i,j,k),uu(3,i,j,k),uu(4,i,j,k), &
                 uu(5,i,j,k),uu(6,i,j,k),uu(7,i,j,k),uu(8,i,j,k), &
                 uu(9,i,j,k)
        write(9) u0(1,i,j,k),u0(2,i,j,k),u0(3,i,j,k),u0(4,i,j,k), &
                 u0(5,i,j,k),u0(6,i,j,k),u0(7,i,j,k),u0(8,i,j,k), &
                 u0(9,i,j,k)
        write(9) uri(1,i,j,k),uri(2,i,j,k),uri(3,i,j,k),uri(4,i,j,k), &
                 uri(5,i,j,k),uri(6,i,j,k),uri(7,i,j,k),uri(8,i,j,k), &
                 uri(9,i,j,k)
        write(9) uri0(1,i,j,k),uri0(2,i,j,k),uri0(3,i,j,k),uri0(4,i,j,k), &
                 uri0(5,i,j,k),uri0(6,i,j,k),uri0(7,i,j,k),uri0(8,i,j,k), &
                 uri0(9,i,j,k)
!        write(9) (uu(m,i,j,k),m=1,nv)
!        write(9) (u0(m,i,j,k),m=1,nv)
!        write(9) (uri(m,i,j,k),m=1,nv)
!        write(9) (uri0(m,i,j,k),m=1,nv)
      enddo
    enddo
  enddo
!
  return
end subroutine store1
!
!--------------------------------------------------------------------
subroutine output1(uri,detg,x1,x2,x3,time,nm1,nd,npe,myrank,myranki,myrankj,&
                   myrankk,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
! - Output simulation data for data analysis
! - Output format
! - x1, x2, x3, \rho, u_tilde^i, p, ink, B^i, \determinant g
!  
  use pram, only : imax, jmax, kmax, nv, iprocs, jprocs, kprocs 
  implicit none
  include 'mpif.h'  
!
  integer :: i, j, k, nm1
  integer :: is, ie, is1, ie1
  integer :: js, je, js1, je1 
  integer :: ks, ke, ks1, ke1
  character*256 :: filename
!
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: x1(imax), x2(jmax), x3(kmax)
  real(8) :: time
!
  integer :: npe, myrank, myranki, myrankj, myrankk, merr, nd
!!! 
  if(myranki .eq. 0) then
    is=is1
  else
    is=is1+nm1
  endif
  if(myranki .eq. iprocs-1) then
    ie=ie1
  else
    ie=ie1-nm1
  endif
!
  if(myrankj .eq. 0) then
    js=js1
  else
    js=js1+nm1
  endif
  if(myrankj .eq. jprocs-1) then
    je=je1
  else
    je=je1-nm1
  endif
!
  if(myrankk .eq. 0) then
    ks=ks1
  else
    ks=ks1+nm1
  endif
  if(myrankk .eq. kprocs-1) then
    ke=ke1
  else
    ke=ke1-nm1
  endif
!
!  do k=ks,ke
!    do j=js,je
!      do i=is,ie
!        if(i .eq. 4 .and. j .eq. 4) then
!          write(*,*) 'k,bx,by,bx=',k,uri(7,i,j,k),uri(8,i,j,k),uri(9,i,j,k)
!        endif
!      enddo
!    enddo
!  enddo

!
! data output for data analysis
!
  write(filename,990) myrank, nd
  open( unit=8,file=filename,form='unformatted',status='unknown')
!
  call mpi_barrier(mpi_comm_world,merr)
!    do k=1, kmax
!      write(*,*) 'k, ro=', k, uri1(1,4,4,k)
!    enddo
  write(8) time, myrank, npe
  write(8) is,ie,js,je,ks,ke
  do i=is,ie
    do j=js,je
      do k=ks,ke
        write(8) x1(i),x2(j),x3(k),uri(1,i,j,k),uri(2,i,j,k), &
                 uri(3,i,j,k),uri(4,i,j,k),uri(5,i,j,k),uri(6,i,j,k), &
                 uri(7,i,j,k),uri(8,i,j,k),uri(9,i,j,k),detg(i,j,k)
      enddo
    enddo
  enddo
!
  close(8)
990 format('structr',i4.4,'-',i4.4,'.outdat')
!
  return
end subroutine output1
!
!--------------------------------------------------------------------
subroutine output2(urib,x1,x2,x3,time,nm1,nd,npe,myrank,myranki,myrankj,&
                   myrankk,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
! - Output simulation data for radiation calculation
! - Output format
! - x1, x2, x3, \rho, u^\alpha, p, bb^\alpha  
  
  use pram, only : imax, jmax, kmax, nv, iprocs, jprocs, kprocs 
  implicit none
  include 'mpif.h'  
!
  integer :: i, j, k, nm1
  integer :: is, ie, is1, ie1
  integer :: js, je, js1, je1 
  integer :: ks, ke, ks1, ke1
  character*256 :: filename2
!
  real(8) :: urib(nv+1,is1:ie1,js1:je1,ks1:ke1)
!- [1]: density, [2-5]: 4-velocity,
!- [6]: gas pressure, [7-10]: 4-magnetic field  
!
  real(8) :: x1(imax), x2(jmax), x3(kmax)
  real(8) :: time
!
  integer :: npe, myrank, myranki, myrankj, myrankk, merr, nd
!!! 
  if(myranki .eq. 0) then
    is=is1
  else
    is=is1+nm1
  endif
  if(myranki .eq. iprocs-1) then
    ie=ie1
  else
    ie=ie1-nm1
  endif
!
  if(myrankj .eq. 0) then
    js=js1
  else
    js=js1+nm1
  endif
  if(myrankj .eq. jprocs-1) then
    je=je1
  else
    je=je1-nm1
  endif
!
  if(myrankk .eq. 0) then
    ks=ks1
  else
    ks=ks1+nm1
  endif
  if(myrankk .eq. kprocs-1) then
    ke=ke1
  else
    ke=ke1-nm1
  endif
!
! data output for data analysis
!
  write(filename2,990) myrank, nd
  open( unit=10,file=filename2,form='unformatted',status='unknown')
!
  call mpi_barrier(mpi_comm_world,merr)
!
  write(10) time, myrank, npe
  write(10) is,ie,js,je,ks,ke
  do i=is,ie
    do j=js,je
      do k=ks,ke
        write(10) x1(i),x2(j),x3(k),urib(1,i,j,k),urib(2,i,j,k), &
                  urib(3,i,j,k),urib(4,i,j,k),urib(5,i,j,k),urib(6,i,j,k), &
                  urib(7,i,j,k),urib(8,i,j,k),urib(9,i,j,k),urib(10,i,j,k)
      enddo
    enddo
  enddo
!
  close(10)
990 format('structr1',i4.4,'-',i4.4,'.outdat')
!
  return
end subroutine output2
!
!--------------------------------------------------------------------
subroutine restar2(uu,u0,uri,uri0,it0,ih0,time,nd,npe,myrank,&
                   is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv 
  implicit none
  include 'mpif.h'
!
  integer :: i, j, k, m, is1, ie1, js1, je1, ks1, ke1
!
  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1), &
             u0(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1), &
             uri0(nv,is1:ie1,js1:je1,ks1:ke1)
!
  integer :: it0, ih0, nd
  real(8) :: time
  integer :: npe, myrank, myrank1, merr
!
  read(9) myrank1, nd
  read(9) it0,ih0
  read(9) time
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
        read(9) uu(1,i,j,k),uu(2,i,j,k),uu(3,i,j,k),uu(4,i,j,k), &
                uu(5,i,j,k),uu(6,i,j,k),uu(7,i,j,k),uu(8,i,j,k), &
                uu(9,i,j,k)
        read(9) u0(1,i,j,k),u0(2,i,j,k),u0(3,i,j,k),u0(4,i,j,k), &
                u0(5,i,j,k),u0(6,i,j,k),u0(7,i,j,k),u0(8,i,j,k), &
                u0(9,i,j,k)
        read(9) uri(1,i,j,k),uri(2,i,j,k),uri(3,i,j,k),uri(4,i,j,k), &
                uri(5,i,j,k),uri(6,i,j,k),uri(7,i,j,k),uri(8,i,j,k), &
                uri(9,i,j,k)
        read(9) uri0(1,i,j,k),uri0(2,i,j,k),uri0(3,i,j,k),uri0(4,i,j,k), &
                uri0(5,i,j,k),uri0(6,i,j,k),uri0(7,i,j,k),uri0(8,i,j,k), &
                uri0(9,i,j,k)

!        read(9) (uu(m,i,j,k),m=1,nv)
!        read(9) (u0(m,i,j,k),m=1,nv)
!        read(9) (uri(m,i,j,k),m=1,nv)
!        read(9) (uri0(m,i,j,k),m=1,nv)
      enddo
    enddo
  enddo
!
!     Check before store for restart                         1998.04.09
!
  if(myrank .ne. myrank1) then
    write(*,*) 'Error: for myrank \= myrank1'
    call mpi_finalize(merr)
  endif
!
!  do k=1,kmax
!    do j=1,jmax
!      do i=1,imax
!        if( uu(1,i,j,k).le.0.d0 ) then
!          write(4,*) 'Warning in restar: density =< 0'
!          write(4,*) 'density =',uu(1,i,j,k),'at i=',i,',j=',j,',k=',k
!          write(4,*) 'stop before restart in restar'
!          stop
!        endif
!        if( uu(5,i,j,k).lt.0.d0 ) then
!          write(4,*) 'Warning in store for restart: pressure < 0'
!          write(4,*) 'pressure=',uu(5,i,j,k),'at i=',i,',j=',j,',k=',k
!          write(4,*) 'stop before restart in restar'
!          stop
!        endif
!      enddo
!    enddo
!  enddo
!
  return
end subroutine restar2
!
!***********************************************************************
subroutine gofend(ifile)
!***********************************************************************
  implicit none
!
  integer :: ifile, ii
  real(8) :: xx

  ii=0
!
 100  continue
      read(ifile,end=200) xx
      ii=ii+1
      goto 100
!
 200  continue
      if( ii.eq.0 ) rewind(ifile)
!
  return
end subroutine gofend
