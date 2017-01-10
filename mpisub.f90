!***********************************************************************
!     for Calculation of Source Terms
!***********************************************************************
!
!---------------------------------------------------------------------
subroutine pararange(n1,n2,npes,myrank,ks,ke)
!---------------------------------------------------------------------
  implicit none
  integer :: n1, n2, npes, myrank, ks, ke, iwork1, iwork2
  
  iwork1=(n2-n1+1)/npes
  iwork2=mod(n2-n1+1,npes)

  ks=myrank*iwork1+n1+min(myrank,iwork2)
  ke=ks+iwork1-1
  if(iwork2 .gt. myrank) then
    ke=ke+1
  endif

  return

end subroutine pararange
!
!---------------------------------------------------------------------
subroutine makecputable(myrank,icputable,myranki,myrankj,myrankk)
!---------------------------------------------------------------------
  use pram, only : iprocs, jprocs, kprocs
  implicit none
  include 'mpif.h'
!
  integer :: i, j, k
  integer :: icputable(-1:iprocs,-1:jprocs,-1:kprocs)
  integer :: myrank, irank, myranki, myrankj, myrankk
  
  do k=-1,kprocs
    do j=-1,jprocs
      do i=-1,iprocs
        icputable(i,j,k)=mpi_proc_null
      enddo
    enddo
  enddo
  irank=0
  do k=0,kprocs-1
    do j=0,jprocs-1
      do i=0,iprocs-1
        icputable(i,j,k)=irank
        if(myrank .eq. irank) then
          myranki=i
          myrankj=j
          myrankk=k
        endif
        irank=irank+1
      enddo
    enddo
  enddo

  return

end subroutine makecputable
!
!---------------------------------------------------------------------
subroutine mpiex(uu,nm1,is1,ie1,js1,je1,ks1,ke1,myranki,myrankj, &
                 myrankk,icputable)
!--------------------------------------------------------------------- 
! MPI data exchange
!
  use pram, only : imax, jmax, kmax, nv, iprocs, jprocs, kprocs
  implicit none
  include 'mpif.h'

  integer :: i, j, k, n, nm1
  integer :: is, ie, is1, ie1
  integer :: js, je, js1, je1 
  integer :: ks, ke, ks1, ke1
 
  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  integer :: mmx, merr, iup, idown, jup, jdown, kup, kdown
  integer :: myranki, myrankj, myrankk
  integer :: icputable(-1:iprocs,-1:jprocs,-1:kprocs)
  integer :: mstatus(mpi_status_size)

  real(8), allocatable :: bufsnd(:,:,:,:), bufrcv(:,:,:,:)
!
!----------------------------------------------------------------------
!  from PE(myranki) to PE(myranki-1) for new da(ix)
!----------------------------------------------------------------------
!
  if(iprocs .ge. 2) then
    iup  =icputable(myranki+1, myrankj, myrankk)
    idown=icputable(myranki-1, myrankj, myrankk)
    if (myranki.eq.(iprocs-1)) iup=mpi_proc_null
    if (myranki.eq.0) idown=mpi_proc_null
!
!    if(myrankj .eq. 0) then
!      js=1
!    else
!      js=js1+nm1
!    endif
!    if(myrankj .eq. jprocs-1) then
!      je=je1
!    else
!      je=je1-nm1
!    endif
!!
!    if(myrankk .eq. 0) then
!      ks=1
!    else
!      ks=ks1+nm1
!    endif
!    if(myrankk .eq. kprocs-1) then
!      ke=ke1
!    else
!      ke=ke1-nm1
!    endif

    allocate(bufsnd(nv,nm1,js1:je1,ks1:ke1), &
             bufrcv(nv,nm1,js1:je1,ks1:ke1),stat=merr)
  
    mmx=nv*nm1*(je1-js1+1)*(ke1-ks1+1)
  
    do k=ks1,ke1
      do j=js1,je1
        do i=1,nm1
          do n=1,nv
            bufsnd(n,i,j,k)=uu(n,is1+nm1-1+i,j,k)
          enddo
        enddo
      enddo
    enddo

    call mpi_barrier(mpi_comm_world,merr)
    call mpi_sendrecv(bufsnd,mmx,mpi_double_precision,idown ,0, &
                      bufrcv,mmx,mpi_double_precision,iup   ,0, &
                      mpi_comm_world,mstatus,merr)

    if (myranki.ne.iprocs-1) then
      do k=ks1,ke1
        do j=js1,je1
          do i=1,nm1
            do n=1,nv
              uu(n,ie1-nm1+i,j,k)=bufrcv(n,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
  
    deallocate(bufsnd,bufrcv,stat=merr)
!
!----------------------------------------------------------------------
!  from PE(myrankk) to PE(myrankk+1) for new da(1)
!----------------------------------------------------------------------
!
    allocate(bufsnd(nv,nm1,js1:je1,ks1:ke1), &
             bufrcv(nv,nm1,js1:je1,ks1:ke1),stat=merr)
    mmx=nv*nm1*(je1-js1+1)*(ke1-ks1+1)

    do k=ks1,ke1
      do j=js1,je1
        do i=1,nm1
          do n=1,nv
            bufsnd(n,i,j,k)=uu(n,ie1-2*nm1+i,j,k)
          enddo
        enddo
      enddo
    enddo

    call mpi_barrier(mpi_comm_world,merr)
    call mpi_sendrecv(bufsnd,mmx,mpi_double_precision,iup  ,1, &
                      bufrcv,mmx,mpi_double_precision,idown,1, &
                      mpi_comm_world,mstatus,merr)

    if (myranki.ne.0) then
      do k=ks1,ke1
        do j=js1,je1
          do i=1,nm1
            do n=1,nv
              uu(n,is1-1+i,j,k)=bufrcv(n,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
!
    deallocate(bufsnd,bufrcv,stat=merr)
!
  endif
!
!----------------------------------------------------------------------
!  from PE(myrankj) to PE(myrankj-1) for new da(ix)
!----------------------------------------------------------------------
  if(jprocs .ge. 2) then
!
    jup  =icputable(myranki, myrankj+1, myrankk)
    jdown=icputable(myranki, myrankj-1, myrankk)
    if (myrankj.eq.(jprocs-1)) jup=mpi_proc_null
    if (myrankj.eq.0) jdown=mpi_proc_null    

!    if(myranki .eq. 0) then
!      is=1
!    else
!      is=is1+nm1
!    endif
!    if(myranki .eq. iprocs-1) then
!      ie=ie1
!    else
!      ie=ie1-nm1
!    endif
!!
!    if(myrankk .eq. 0) then
!      ks=1
!    else
!      ks=ks1+nm1
!    endif
!    if(myrankk .eq. kprocs-1) then
!      ke=ke1
!    else
!      ke=ke1-nm1
!    endif
!
    allocate(bufsnd(nv,is1:ie1,nm1,ks1:ke1), &
             bufrcv(nv,is1:ie1,nm1,ks1:ke1),stat=merr)
    mmx=nv*(ie1-is1+1)*nm1*(ke1-ks1+1)
  
    do j=1,nm1
      do k=ks1,ke1
        do i=is1,ie1
          do n=1,nv
            bufsnd(n,i,j,k)=uu(n,i,js1+nm1-1+j,k)
          enddo
        enddo
      enddo
    enddo

    call mpi_barrier(mpi_comm_world,merr)
    call mpi_sendrecv(bufsnd,mmx,mpi_double_precision,jdown ,2, &
                      bufrcv,mmx,mpi_double_precision,jup   ,2, &
                      mpi_comm_world,mstatus,merr)

    if (myrankj.ne.jprocs-1) then
      do j=1,nm1
        do k=ks1,ke1
          do i=is1,ie1
            do n=1,nv
              uu(n,i,je1-nm1+j,k)=bufrcv(n,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
!
    deallocate(bufsnd,bufrcv,stat=merr)
!
!----------------------------------------------------------------------
!  from PE(myrankj) to PE(myrankj+1) for new da(1)
!----------------------------------------------------------------------
!
    allocate(bufsnd(nv,is1:ie1,nm1,ks1:ke1), &
             bufrcv(nv,is1:ie1,nm1,ks1:ke1),stat=merr)
    mmx=nv*(ie1-is1+1)*nm1*(ke1-ks1+1)

    do j=1,nm1
      do k=ks1,ke1
        do i=is1,ie1
          do n=1,nv
            bufsnd(n,i,j,k)=uu(n,i,je1-2*nm1+j,k)
          enddo
        enddo
      enddo
    enddo

    call mpi_barrier(mpi_comm_world,merr)
    call mpi_sendrecv(bufsnd,mmx,mpi_double_precision,jup  ,3, &
                      bufrcv,mmx,mpi_double_precision,jdown,3, &
                      mpi_comm_world,mstatus,merr)

    if (myrankj.ne.0) then
      do j=1,nm1
        do k=ks1,ke1
          do i=is1,ie1
            do n=1,nv
              uu(n,i,js1-1+j,k)=bufrcv(n,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
!
    deallocate(bufsnd,bufrcv,stat=merr)
!
  endif
!
!----------------------------------------------------------------------
!  from PE(myrankk) to PE(myrankk-1) for new da(ix)
!----------------------------------------------------------------------
  if(kprocs .ge. 2) then
    kup  =icputable(myranki, myrankj, myrankk+1)
    kdown=icputable(myranki, myrankj, myrankk-1)
    if (myrankk.eq.(kprocs-1)) kup=mpi_proc_null
    if (myrankk.eq.0) kdown=mpi_proc_null
!
!    if(myranki .eq. 0) then
!      is=1
!    else
!      is=is1+nm1
!    endif
!    if(myranki .eq. iprocs-1) then
!      ie=ie1
!    else
!      ie=ie1-nm1
!    endif
!!
!    if(myrankj .eq. 0) then
!      js=1
!    else
!      js=js1+nm1
!    endif
!    if(myrankj .eq. jprocs-1) then
!      je=je1
!    else
!      je=je1-nm1
!    endif
!
    allocate(bufsnd(nv,is1:ie1,js1:je1,nm1), &
             bufrcv(nv,is1:ie1,js1:je1,nm1),stat=merr)
    mmx=nv*(ie1-is1+1)*(je1-js1+1)*nm1
  
    do k=1,nm1
      do j=js1,je1
        do i=is1,ie1
          do n=1,nv
            bufsnd(n,i,j,k)=uu(n,i,j,ks1+nm1-1+k)
          enddo
        enddo
      enddo
    enddo

    call mpi_barrier(mpi_comm_world,merr)
    call mpi_sendrecv(bufsnd,mmx,mpi_double_precision,kdown ,4, &
                      bufrcv,mmx,mpi_double_precision,kup   ,4, &
                      mpi_comm_world,mstatus,merr)

    if (myrankk.ne.kprocs-1) then
      do k=1,nm1
        do j=js1,je1
          do i=is1,ie1
            do n=1,nv
              uu(n,i,j,ke1-nm1+k)=bufrcv(n,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
!
    deallocate(bufsnd,bufrcv,stat=merr)
!
!----------------------------------------------------------------------
!  from PE(myrankk) to PE(myrankk+1) for new da(1)
!----------------------------------------------------------------------
!
    allocate(bufsnd(nv,is1:ie1,js1:je1,nm1), &
             bufrcv(nv,is1:ie1,js1:je1,nm1),stat=merr)
    mmx=nv*(ie1-is1+1)*(je1-js1+1)*nm1

    do k=1,nm1
      do j=js1,je1
        do i=is1,ie1
          do n=1,nv
            bufsnd(n,i,j,k)=uu(n,i,j,ke1-2*nm1+k)
          enddo
        enddo
      enddo
    enddo

    call mpi_barrier(mpi_comm_world,merr)
    call mpi_sendrecv(bufsnd,mmx,mpi_double_precision,kup  ,5, &
                      bufrcv,mmx,mpi_double_precision,kdown,5, &
                      mpi_comm_world,mstatus,merr)

    if (myrankk.ne.0) then
      do k=1,nm1
        do j=js1,je1
          do i=is1,ie1
            do n=1,nv
              uu(n,i,j,ks1-1+k)=bufrcv(n,i,j,k)
            enddo
          enddo
        enddo
      enddo
    endif
!
    deallocate(bufsnd,bufrcv,stat=merr)
!
  endif
!
  return
end subroutine mpiex
