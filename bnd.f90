!***********************************************************************
!     for BOUNDARY CONDITION
!***********************************************************************
subroutine bnd4(uu,uo,x1,nm1,is1,ie1,js1,je1,ks1,ke1, &
                myranki,myrankj,myrankk,icputable)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv, iprocs, jprocs, kprocs
  implicit none
  
  integer :: nm1, is1, ie1, js1, je1, ks1, ke1
  integer :: icputable(-1:iprocs,-1:jprocs,-1:kprocs)
  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1), &
             uo(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: x1(imax)
  integer :: myranki, myrankj, myrankk
!
  call bndpefb4(uu,nm1,is1,ie1,js1,je1,ks1,ke1, &
                myranki,myrankj,myrankk,icputable)
  call bndryfb4(uu,uo,x1,nm1,is1,ie1,js1,je1,ks1,ke1, &
                myranki,myrankj,myrankk)
!
  return
end subroutine bnd4
!
!--------------------------------------------------------------------
subroutine bndpefb4(uu,nm1,is1,ie1,js1,je1,ks1,ke1, &
                    myranki,myrankj,myrankk,icputable)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv, iprocs, jprocs, kprocs, &
                   iboux1in, iboux1ot, iboux2in, iboux2ot, iboux3in, iboux3ot  
  implicit none
  include 'mpif.h'
  
  integer :: i, j, k, n, m
  integer :: nm1, nm2, ibnd0, ibnd1, is1, ie1, js1, je1, ks1, ke1
  integer :: icputable(-1:iprocs,-1:jprocs,-1:kprocs)
  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)

  integer :: myranki, myrankj, myrankk
  integer :: mmx, merr, msnd, mrcv
  integer :: mstatus(mpi_status_size), mreq
  real(8), allocatable :: bufsnd(:,:,:), bufrcv(:,:,:)
!
  nm2=2*nm1-1
!
  ibnd0=1+nm1
  ibnd1=imax-nm1
  call mpi_barrier(mpi_comm_world,merr)
!
  if(iprocs .eq. 1) then

    do k=ks1+nm1,ke1-nm1
      do j=js1+nm1,je1-nm1
        do i=1,nm1
          do m=1,nv 
           if( iboux1in(m) .eq. 1 .or. iboux1ot(m) .eq. 1 ) then
!             uu(m,n,j,k)=uu(m,imax-nm1-n,j,k)
!             uu(m,imax-n,j,k)=uu(m,nm2-n,j,k)

!             uu(m,n,j,k)=uu(m,imax-nm2+n,j,k)
!             uu(m,imax-n,j,k)=uu(m,nm2-n,j,k)

              uu(m,ibnd0-i,j,k)=uu(m,ibnd1+1-i,j,k)
              uu(m,ibnd1+i,j,k)=uu(m,ibnd0-1+i,j,k)

            endif
          enddo
        enddo
      enddo
    enddo
!
  else
    do m=1,nv
      if( iboux1in(m) .eq. 1 ) then
!
        if(myranki .eq. 0 .or. myranki .eq. iprocs-1) then
          allocate(bufsnd(nm1,js1:je1,ks1:ke1), &
                   bufrcv(nm1,js1:je1,ks1:ke1), stat=merr)
          mmx=nm1*(je1-js1+1)*(ke1-ks1+1)
          msnd=icputable(0,myrankj,myrankk)
          mrcv=icputable(iprocs-1,myrankj,myrankk)

          if(myranki .eq. iprocs-1) then
            do k=ks1,ke1
              do j=js1,je1
                do i=1,nm1
                  bufsnd(i,j,k)=uu(m,ibnd1+1-i,j,k)
                enddo
              enddo
            enddo
            call mpi_isend(bufsnd,mmx,mpi_double_precision, &
                       msnd,0,mpi_comm_world,mreq,merr)
            call mpi_wait(mreq,mstatus,merr)
          endif

          if(myranki .eq. 0) then
            call mpi_irecv(bufrcv,mmx,mpi_double_precision, &
                       mrcv,0,mpi_comm_world,mreq,merr)
            call mpi_wait(mreq,mstatus,merr)          
            do k=ks1,ke1
              do j=js1,je1
                do i=1,nm1             
                  uu(m,ibnd0-i,j,k)=bufrcv(i,j,k) 
                enddo
              enddo
            enddo
          endif
!
        endif
        deallocate(bufsnd, bufrcv, stat=merr)
      endif
!
      if(iboux1ot(m) .eq. 1 ) then
!
        if(myranki .eq. 0 .or. myranki .eq. iprocs-1) then 
          allocate(bufsnd(nm1,js1:je1,ks1:ke1), &
                   bufrcv(nm1,js1:je1,ks1:ke1), stat=merr)
          mmx=nm1*(je1-js1+1)*(ke1-ks1+1)
          msnd=icputable(iprocs-1,myrankj,myrankk)
          mrcv=icputable(0,myrankj,myrankk)
        
          if(myranki .eq. 0) then
            do k=ks1,ke1
              do j=js1,je1
                do i=1,nm1
                  bufsnd(i,j,k)=uu(m,ibnd0-1+i,j,k)
                enddo
              enddo
            enddo
            call mpi_isend(bufsnd,mmx,mpi_double_precision, &
                       msnd,0,mpi_comm_world,mreq,merr)
            call mpi_wait(mreq,mstatus,merr)          
          endif
!
          if(myranki .eq. iprocs-1) then
            call mpi_irecv(bufrcv,mmx,mpi_double_precision, &
                           mrcv,0,mpi_comm_world,mreq,merr)
            call mpi_wait(mreq,mstatus,merr)             
            do k=ks1,ke1
              do j=js1,je1
                do i=1,nm1             
                  uu(m,ibnd1+i,j,k)=bufrcv(i,j,k) 
                enddo
              enddo
            enddo
          endif
!
        endif  
        deallocate(bufsnd, bufrcv, stat=merr)
      endif
!
    enddo
!      
  endif
!
!
  ibnd0=1+nm1
  ibnd1=jmax-nm1
  call mpi_barrier(mpi_comm_world,merr)
!
  if(jprocs .eq. 1) then

    do k=ks1+nm1,ke1-nm1
      do j=1,nm1
        do i=is1+nm1,ie1-nm1
          do m=1,nv 
           if( iboux2in(m) .eq. 1 .or. iboux2ot(m) .eq. 1 ) then
!             uu(m,i,n,k)=uu(m,i,jmax-nm1-n,k)
!             uu(m,i,jmax-n,k)=uu(m,i,nm2-n,k)

!             uu(m,i,n,k)=uu(m,i,jmax-nm2+n,k)
!             uu(m,i,jmax-n,k)=uu(m,i,nm2-n,k)

              uu(m,i,ibnd0-j,k)=uu(m,i,ibnd1+1-j,k)
              uu(m,i,ibnd1+j,k)=uu(m,i,ibnd0-1+j,k)

            endif
          enddo
        enddo
      enddo
    enddo
!
  else
    do m=1,nv
      if( iboux2in(m) .eq. 1 ) then
!
        if(myrankj .eq. 0 .or. myrankj .eq. jprocs-1) then
          allocate(bufsnd(is1:ie1,nm1,ks1:ke1), &
                   bufrcv(is1:ie1,nm1,ks1:ke1), stat=merr)
          mmx=(ie1-is1+1)*nm1*(ke1-ks1+1)
          msnd=icputable(myranki,0,myrankk)
          mrcv=icputable(myranki,jprocs-1,myrankk)

          if(myrankj .eq. jprocs-1) then
            do k=ks1,ke1
              do j=1,nm1
                do i=is1,ie1
                  bufsnd(i,j,k)=uu(m,i,ibnd1+1-j,k)
                enddo
              enddo
            enddo
            call mpi_isend(bufsnd,mmx,mpi_double_precision, &
                       msnd,0,mpi_comm_world,mreq,merr)
            call mpi_wait(mreq,mstatus,merr)
          endif

          if(myrankj .eq. 0) then
            call mpi_irecv(bufrcv,mmx,mpi_double_precision, &
                       mrcv,0,mpi_comm_world,mreq,merr)
            call mpi_wait(mreq,mstatus,merr)          
            do k=ks1,ke1
              do j=1,nm1
                do i=is1,ie1             
                  uu(m,i,ibnd0-j,k)=bufrcv(i,j,k) 
                enddo
              enddo
            enddo
          endif
!
        endif
        deallocate(bufsnd, bufrcv, stat=merr)
      endif
!
      if(iboux2ot(m) .eq. 1 ) then
!
        if(myrankj .eq. 0 .or. myrankj .eq. jprocs-1) then 
          allocate(bufsnd(is1:ie1,nm1,ks1:ke1), &
                   bufrcv(is1:ie1,nm1,ks1:ke1), stat=merr)
          mmx=(ie1-is1+1)*nm1*(ke1-ks1+1)
          msnd=icputable(myranki,jprocs-1,myrankk)
          mrcv=icputable(myranki,0,myrankk)
        
          if(myrankj .eq. 0) then
            do k=ks1,ke1
              do j=1,nm1
                do i=is1,ie1
                  bufsnd(i,j,k)=uu(m,i,ibnd0-1+j,k)
                enddo
              enddo
            enddo
            call mpi_isend(bufsnd,mmx,mpi_double_precision, &
                       msnd,0,mpi_comm_world,mreq,merr)
            call mpi_wait(mreq,mstatus,merr)          
          endif
!
          if(myrankj .eq. jprocs-1) then
            call mpi_irecv(bufrcv,mmx,mpi_double_precision, &
                           mrcv,0,mpi_comm_world,mreq,merr)
            call mpi_wait(mreq,mstatus,merr)             
            do k=ks1,ke1
              do j=1,nm1
                do i=is1,ie1             
                  uu(m,i,ibnd1+j,k)=bufrcv(i,j,k) 
                enddo
              enddo
            enddo
          endif
!
        endif  
        deallocate(bufsnd, bufrcv, stat=merr)
      endif
!
    enddo
!      
  endif
!
!
  ibnd0=1+nm1
  ibnd1=kmax-nm1
!
  call mpi_barrier(mpi_comm_world,merr)
!
  if(kprocs .eq. 1) then
!
    do k=1,nm1
      do j=js1+nm1,je1-nm1
        do i=is1+nm1,ie1-nm1
          do m=1,nv
            if( iboux3in(m) .eq. 1 .or. iboux3ot(m) .eq. 1 ) then
!             uu(m,i,j,n)=uu(m,i,j,kmax-nm1-n)
!             uu(m,i,j,kmax-n)=uu(m,i,j,nm2-n)

!             uu(m,i,j,n)=uu(m,i,j,kmax-nm2+n)
!             uu(m,i,j,kmax-n)=uu(m,i,j,nm2-n)
 
              uu(m,i,j,ibnd0-k)=uu(m,i,j,ibnd1+1-k)
              uu(m,i,j,ibnd1+k)=uu(m,i,j,ibnd0-1+k)

            endif
          enddo
        enddo
      enddo
    enddo
!
  else
!
    do m=1,nv
      if( iboux3in(m) .eq. 1 ) then
  
        if(myrankk .eq. 0 .or. myrankk .eq. kprocs-1) then
          allocate(bufsnd(is1:ie1,js1:je1,nm1), &
                   bufrcv(is1:ie1,js1:je1,nm1), stat=merr)
          mmx=(ie1-is1+1)*(je1-js1+1)*nm1
          msnd=icputable(myranki,myrankj,0)
          mrcv=icputable(myranki,myrankj,kprocs-1)

          if(myrankk .eq. kprocs-1) then
            do k=1,nm1
              do j=js1,je1
                do i=is1,ie1
                  bufsnd(i,j,k)=uu(m,i,j,ibnd1+1-k)
                enddo
              enddo
            enddo
            call mpi_isend(bufsnd,mmx,mpi_double_precision, &
                       msnd,0,mpi_comm_world,mreq,merr)
            call mpi_wait(mreq,mstatus,merr)
          endif

          if(myrankk .eq. 0) then
            call mpi_irecv(bufrcv,mmx,mpi_double_precision, &
                       mrcv,0,mpi_comm_world,mreq,merr)
            call mpi_wait(mreq,mstatus,merr)          
            do k=1,nm1
              do j=js1,je1
                do i=is1,ie1             
                  uu(m,i,j,ibnd0-k)=bufrcv(i,j,k) 
                enddo
              enddo
            enddo
          endif
!
        endif
        deallocate(bufsnd, bufrcv, stat=merr)
      endif
!
      if(iboux3ot(m) .eq. 1 ) then
!
        if(myrankk .eq. 0 .or. myrankk .eq. kprocs-1) then 
          allocate(bufsnd(is1:ie1,js1:je1,nm1), &
                   bufrcv(is1:ie1,js1:je1,nm1), stat=merr)
          mmx=(ie1-is1+1)*(je1-js1+1)*nm1
          msnd=icputable(myranki,myrankj,kprocs-1)
          mrcv=icputable(myranki,myrankj,0)
        
          if(myrankk .eq. 0) then
            do k=1,nm1
              do j=js1,je1
                do i=is1,ie1
                  bufsnd(i,j,k)=uu(m,i,j,ibnd0-1+k)
                enddo
              enddo
            enddo
            call mpi_isend(bufsnd,mmx,mpi_double_precision, &
                       msnd,0,mpi_comm_world,mreq,merr)
            call mpi_wait(mreq,mstatus,merr)          
          endif
!
          if(myrankk .eq. kprocs-1) then
            call mpi_irecv(bufrcv,mmx,mpi_double_precision, &
                           mrcv,0,mpi_comm_world,mreq,merr)
            call mpi_wait(mreq,mstatus,merr)             
            do k=1,nm1
              do j=js1,je1
                do i=is1,ie1             
                  uu(m,i,j,ibnd1+k)=bufrcv(i,j,k) 
                enddo
              enddo
            enddo
          endif
!
        endif  
        deallocate(bufsnd, bufrcv, stat=merr)
      endif
!
    enddo
!
  endif
!
  return
end subroutine bndpefb4
!
!--------------------------------------------------------------------
subroutine bndryfb4(uu,uo,x1,nm1,is1,ie1,js1,je1,ks1,ke1, &
                    myranki,myrankj,myrankk)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv, iprocs, jprocs, kprocs, &
                   iboux1in, iboux1ot, & 
                   iboux2in, iboux2ot, iboux3in, iboux3ot  
  implicit none
  
  integer :: i, j, k, n, m
  integer :: nm1, is1, ie1, js1, je1, ks1, ke1
  integer :: myranki, myrankj, myrankk
  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1), &
             uo(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: x1(imax)
  real(8) :: rb
!
! ibou
! 1   :    periodic boundary
! 2   :    fixed boundary
! 3   :    Neumann boundary
! 4   :    free (extrapolate) boundary
! 5   :    copy boundary, (u_-1=u_0, u_-2=u_0)       
! 7   :    radiative boundary without eigenvalue
! 12   :    antisymmetric boundary condition (u_-1=-u_0 & u_-2=-u_1)
! 13   :    regid wall boundary condition (u_0=0.0)
! 15   :    special reflecting boundary 
! 14  :    symmetric boundary (u_-1=u_0 & u_-2=u_1)
!           for jet propagation (z-direction only)
! 17   :    special radiative boundary without eigenvalue
!           for jet propagation (z-direction only)
! 19   :    special boundary of u_0=u_1
!           for jet propagation (z-direction only)
!
!     x-direction
!
  do m = 1, nv
!
    if( myranki .eq. 0 .or. iprocs .eq. 1) then
!
      if( iboux1in(m) .eq. 1 ) then
!
      elseif( iboux1in(m) .eq. 2 ) then

        do k = ks1, ke1
          do j = js1, je1
            do n = 1, nm1           
              uu(m,n,j,k)=uo(m,n,j,k)
            enddo
          enddo
        enddo
 
      elseif( iboux1in(m) .eq. 3 ) then
 
        do k = ks1, ke1
          do j = js1, je1
            do n = nm1, 1, -1
              uu(m,n,j,k)=(4.0*uu(m,n+1,j,k)-uu(m,n+2,j,k))/3.0
            enddo
          enddo
        enddo
 
      elseif( iboux1in(m) .eq. 4 ) then
 
        do k = ks1, ke1
          do j = js1, je1
            do n = nm1, 1, -1
              uu(m,n,j,k)=2.0*uu(m,n+1,j,k)-uu(m,n+2,j,k)
            enddo
          enddo
        enddo
 
      elseif( iboux1in(m) .eq. 5 ) then
 
        do k = ks1, ke1
          do j = js1, je1
            do n = nm1, 1, -1
              uu(m,n,j,k)=uu(m,1+nm1,j,k)
            enddo
          enddo
        enddo
 
      elseif( iboux1in(m) .eq. 12 ) then
 
        do k = ks1, ke1
          do j = js1, je1
            do n = 1, nm1
!              uu(m,n,j,k)=-uu(m,nm1+n,j,k)
              uu(m,n,j,k)=-uu(m,2*nm1-n+1,j,k)
            enddo
          enddo
        enddo
 
      elseif( iboux1in(m) .eq. 13 ) then
    
        do k = ks1, ke1
          do j = js1, je1
            do n = 1, nm1
              uu(m,n,j,k)=0.0d0
            enddo
          enddo
        enddo

      elseif( iboux1in(m) .eq. 14 ) then
 
        do k = ks1, ke1
          do j = js1, je1
            do n = 1, nm1
              uu(m,n,j,k)=uu(m,2*nm1-n+1,j,k)
            enddo
          enddo
        enddo

      elseif( iboux1in(m) .eq. 7 ) then
  
        do k = ks1, ke1
          do j = js1, je1
            do n = nm1, 1, -1
              uu(m,n,j,k)=uo(m,n,j,k)+uu(m,n+1,j,k)-uo(m,n+1,j,k)
            enddo
          enddo
        enddo
  
      else
!      write(4,*) 'inner boundary setting error for x1 axis'
!      write(4,*) 'stop in bndry'
      endif
!
    endif
!
    if(myranki .eq. iprocs-1 .or. iprocs .eq. 1) then

      if( iboux1ot(m) .eq. 1 ) then
!
      elseif(iboux1ot(m) .eq. 2 ) then

        do k = ks1, ke1
          do j = js1, je1
            do n = 0, nm1-1
              uu(m,imax-n,j,k)=uo(m,imax-n,j,k)
            enddo
          enddo
        enddo
 
      elseif(iboux1ot(m) .eq. 3 ) then

        do k = ks1, ke1
          do j = js1, je1
            do n = nm1-1, 0, -1
              uu(m,imax-n,j,k)=(4.0*uu(m,imax-n-1,j,k)-uu(m,imax-n-2,j,k))/3.0
            enddo
          enddo
        enddo
 
      elseif(iboux1ot(m) .eq. 4 ) then

        do k = ks1, ke1
          do j = js1, je1
            do n = nm1-1, 0, -1
              uu(m,imax-n,j,k)=2.0*uu(m,imax-n-1,j,k)-uu(m,imax-n-2,j,k)
            enddo
          enddo
        enddo
 
      elseif(iboux1ot(m) .eq. 5 ) then 
        do k = ks1, ke1
          do j = js1, je1
            do n = nm1-1, 0, -1
              uu(m,imax-n,j,k)=uu(m,imax-nm1,j,k)
            enddo
          enddo
        enddo
 
      elseif(iboux1ot(m) .eq. 12 ) then
 
        do k = ks1, ke1
          do j = js1, je1
            do n = 0, nm1-1
!              uu(m,imax-n,j,k)=-uu(m,imax-nm1-n,j,k)
              uu(m,imax-n,j,k)=-uu(m,imax-2*nm1+1+n,j,k)
            enddo
          enddo
        enddo

      elseif(iboux1ot(m) .eq. 13 ) then

        do k = ks1, ke1
          do j = js1, je1
            do n = 0, nm1-1
              uu(m,imax-n,j,k)=0.d0
            enddo
          enddo
        enddo

      elseif( iboux1ot(m) .eq. 14 ) then
 
        do k = ks1, ke1
          do j = js1, je1
            do n = 0, nm1-1
              uu(m,imax-n,j,k)=uu(m,imax-2*nm1+1+n,j,k)
            enddo
          enddo
        enddo

      elseif(iboux1ot(m) .eq. 7 ) then

        do k = ks1, ke1
          do j = js1, je1
            do n = nm1-1, 0, -1
             uu(m,imax-n,j,k)=uo(m,imax-n,j,k) &
                                 +uu(m,imax-n-1,j,k)-uo(m,imax-n-1,j,k)
            enddo
          enddo
        enddo
!
      elseif(iboux1ot(m) .eq. 18) then
         
        do k = ks1, ke1
          do j = js1, je1
            do n = nm1-1, 0, -1      
              if(uu(2,imax-nm1,j,k) .lt. 0.d0) then
                uu(m,imax-n,j,k)=0.d0
              else
                uu(m,imax-n,j,k)=uu(m,imax-nm1,j,k)
              endif   
            enddo
          enddo
        enddo
!       
      else
!        write(4,*) 'outer boundary setting error for x1 axis'
!        write(4,*) 'stop in bndry'
      endif
!
    endif
!
  enddo
!
!     y-direction
!

  do m = 1, nv

!
    if( myrankj .eq. 0 .or. jprocs .eq. 1) then
!
      if( iboux2in(m) .eq. 1 ) then
!
      elseif( iboux2in(m) .eq. 2 ) then

        do k = ks1, ke1
          do i = is1, ie1
            do n = 1, nm1
              uu(m,i,n,k)=uo(m,i,n,k)
            enddo
          enddo
        enddo
 
      elseif( iboux2in(m) .eq. 3 ) then

        do k = ks1, ke1
          do i = is1, ie1
            do n = nm1, 1, -1
              uu(m,i,n,k)=(4.0*uu(m,i,n+1,k)-uu(m,i,n+2,k))/3.0
            enddo
          enddo
        enddo
 
      elseif( iboux2in(m) .eq. 4 ) then

        do k = ks1, ke1
          do i = is1, ie1
            do n = nm1, 1, -1
              uu(m,i,n,k)=2.0*uu(m,i,n+1,k)-uu(m,i,n+2,k)
            enddo
          enddo
        enddo
 
      elseif( iboux2in(m) .eq. 5 ) then

        do k = ks1, ke1
          do i = is1, ie1
            do n = nm1, 1, -1
              uu(m,i,n,k)=uu(m,i,1+nm1,k)
            enddo
          enddo
        enddo
 
      elseif( iboux2in(m) .eq. 12 ) then

        do k = ks1, ke1
          do i = is1, ie1
            do n = 1, nm1
!              uu(m,i,n,k)=-uu(m,i,nm1+n+1,k)
              uu(m,i,n,k)=-uu(m,i,2*nm1-n+1,k)
            enddo
          enddo
        enddo

      elseif( iboux2in(m) .eq. 13 ) then

        do k = ks1, ke1
          do i = is1, ie1
            do n = 1, nm1
              uu(m,i,n,k)=0.0d0
            enddo
          enddo
        enddo

      elseif( iboux2in(m) .eq. 14 ) then
 
        do k = ks1, ke1
          do i = is1, ie1
            do n = 1, nm1
              uu(m,i,n,k)=uu(m,i,2*nm1-n+1,k)
            enddo
          enddo
        enddo


      elseif( iboux2in(m) .eq. 7 ) then

        do k = ks1, ke1
          do i = is1, ie1
            do n = nm1, 1, -1
              uu(m,i,n,k)=uo(m,i,n,k)+uu(m,i,n+1,k)-uo(m,i,n+1,k)
            enddo
          enddo
        enddo
        
      else
!        write(4,*) 'inner boundary setting error for x2 axis'
!        write(4,*) 'stop in bndry'
      endif
!
    endif
!
    if(myrankj .eq. jprocs-1 .or. jprocs .eq. 1) then

      if( iboux2ot(m) .eq. 1 ) then
!
      elseif( iboux2ot(m) .eq. 2 ) then

        do k = ks1, ke1
          do i = is1, ie1
            do n = 0, nm1-1
              uu(m,i,jmax-n,k)=uo(m,i,jmax-n,k)
            enddo
          enddo
        enddo
 
      elseif( iboux2ot(m) .eq. 3 ) then

        do k = ks1, ke1
          do i = is1, ie1
            do n = nm1-1, 0, -1
              uu(m,i,jmax-n,k)=(4.0*uu(m,i,jmax-n-1,k)-uu(m,i,jmax-n-2,k))/3.0
            enddo
          enddo
        enddo
 
      elseif(iboux2ot(m) .eq. 4 ) then

        do k = ks1, ke1
          do i = is1, ie1
            do n = nm1-1, 0, -1
              uu(m,i,jmax-n,k)=2.0*uu(m,i,jmax-n-1,k)-uu(m,i,jmax-n-2,k)
            enddo
          enddo
        enddo
 
      elseif(iboux2ot(m) .eq. 5 ) then

        do k = ks1, ke1
          do i = is1, ie1
            do n = nm1-1, 0, -1
              uu(m,i,jmax-n,k)=uu(m,i,jmax-nm1,k)
            enddo
          enddo
        enddo
 
      elseif(iboux2ot(m) .eq. 12 ) then

        do k = ks1, ke1
          do i = is1, ie1
            do n = 0, nm1-1
!              uu(m,i,jmax-n,k)=-uu(m,i,jmax-nm1-n,k)
              uu(m,i,jmax-n,k)=-uu(m,i,jmax-2*nm1+1+n,k)
            enddo
          enddo
        enddo

      elseif(iboux2ot(m) .eq. 13 ) then

        do k = ks1, ke1
          do i = is1, ie1
            do n = 0, nm1-1
              uu(m,i,jmax-n,k)=0.0d0
            enddo
          enddo
        enddo

      elseif( iboux2ot(m) .eq. 14 ) then
 
        do k = ks1, ke1
          do i = is1, ie1
            do n = 0, nm1-1
              uu(m,i,jmax-n,k)=uu(m,i,jmax-2*nm1+1+n,k)
            enddo
          enddo
        enddo

      elseif(iboux2ot(m) .eq. 7 ) then

        do k = ks1, ke1
          do i = is1, ie1
            do n = nm1-1, 0, -1
              uu(m,i,jmax-n,k)=uo(m,i,jmax-n,k)&
                                +uu(m,i,jmax-n-1,k)-uo(m,i,jmax-n-1,k)
            enddo
          enddo
        enddo
 
      else
!        write(4,*) 'outer boundary setting error for x2 axis'
!        write(4,*) 'stop in bndry'
      endif
!
    endif
!
  enddo
!
!     z-direction
!
  do m = 1, nv
!
    if(myrankk .eq. 0 .or. kprocs .eq. 1) then
!
      if( iboux3in(m) .eq. 1 ) then

      elseif(iboux3in(m) .eq. 2 ) then

        do j = js1, je1
          do i = is1, ie1
            do n = 1, nm1
              uu(m,i,j,n)=uo(m,i,j,n)
            enddo
          enddo
        enddo
 
      elseif(iboux3in(m) .eq. 3 ) then

        do j = js1, je1
          do i = is1, ie1
            do n = nm1, 1, -1
              uu(m,i,j,n)=(4.0*uu(m,i,j,n+1)-uu(m,i,j,n+2))/3.0
            enddo
          enddo
        enddo
 
      elseif(iboux3in(m) .eq. 4 ) then

        do j = js1, je1
          do i = is1, ie1
            do n = nm1, 1, -1
              uu(m,i,j,n)=2.0*uu(m,i,j,n+1)-uu(m,i,j,n+2)
            enddo
          enddo
        enddo
 
      elseif(iboux3in(m) .eq. 5 ) then

        do j = js1, je1
          do i = is1, ie1
            do n = nm1, 1, -1
              uu(m,i,j,n)=uu(m,i,j,1+nm1)
            enddo
          enddo
        enddo
 
      elseif(iboux3in(m) .eq. 12 ) then

        do j = js1, je1
          do i = is1, ie1
            do n = 1, nm1
!              uu(m,i,j,n)=-uu(m,i,j,nm1+n)
              uu(m,i,j,n)=-uu(m,i,j,2*nm1-n+1)
            enddo
          enddo
        enddo

      elseif(iboux3in(m) .eq. 13 ) then

        do j = js1, je1
          do i = is1, ie1
            do n = 1, nm1
              uu(m,i,j,n)=0.d0
            enddo
          enddo
        enddo

      elseif( iboux3in(m) .eq. 14 ) then
 
        do j = js1, je1
          do i = is1, ie1
            do n = 1, nm1
              uu(m,i,j,n)=uu(m,i,j,2*nm1-n+1)
            enddo
          enddo
        enddo

      elseif(iboux3in(m) .eq. 7 ) then

        do j = js1, je1
          do i = is1, ie1
            do n = nm1, 1, -1
              uu(m,i,j,n)=uo(m,i,j,n)+uu(m,i,j,n+1)-uo(m,i,j,n+1)
            enddo
          enddo
        enddo

      elseif( iboux3in(m) .eq. 15 ) then

        rb=1.0d0
        do j= js1, je1
          do i= is1, ie1
            do n = nm1, 1, -1
          
              if (x1(i) .le. rb) then
                uu(m,i,j,n)=uo(m,i,j,n)
              else
                uu(m,i,j,n)=uu(m,i,j,1+nm1)
              endif
            enddo
          enddo
        enddo

      elseif( iboux3in(m) .eq. 17 ) then
!
! rb : jet beam radius
!
       rb=1.d0
        do j= js1, je1
          do i= is1, ie1
            do n = nm1, 1, -1
          
              if (x1(i) .le. rb) then
                uu(m,i,j,n)=uo(m,i,j,n)
              else
                uu(m,i,j,n)=uo(m,i,j,n)+uu(m,i,j,n+1)-uo(m,i,j,n+1)
              endif
            enddo
          enddo
        enddo

      elseif( iboux3in(m) .eq. 19 ) then

        rb=1.d0

        do j= js1, je1
          do i= is1, ie1
            do n = 1, nm1

              if (x1(i) .le. rb) then
                uu(m,i,j,n)=uo(m,i,j,n)
              else
                uu(m,i,j,n)=-uu(m,i,j,nm1+n)
              endif
            enddo
          enddo
        enddo
!
      else
!        write(4,*) 'inner boundary setting error for x3 axis'
!        write(4,*) 'stop in bndry'
      endif
!
    endif
!
    if(myrankk .eq. kprocs-1 .or. kprocs .eq. 1) then
!
      if( iboux3ot(m) .eq. 1 ) then
!
      elseif( iboux3ot(m) .eq. 2 ) then

        do j = js1, je1
          do i = is1, ie1
            do n = 0, nm1-1
              uu(m,i,j,kmax-n)=uo(m,i,j,kmax-n)
            enddo
          enddo
        enddo
 
      elseif( iboux3ot(m) .eq. 3 ) then

        do j = js1, je1
          do i = is1, ie1
            do n = nm1-1, 0, -1
              uu(m,i,j,kmax-n)=(4.0*uu(m,i,j,kmax-n-1)-uu(m,i,j,kmax-n-2))/3.0
            enddo
          enddo
        enddo
 
      elseif(iboux3ot(m) .eq. 4 ) then

        do j = js1, je1
          do i = is1, ie1
            do n = nm1-1, 0, -1
              uu(m,i,j,kmax-n)=2.0*uu(m,i,j,kmax-n-1)-uu(m,i,j,kmax-n-2)
            enddo
          enddo
        enddo
 
      elseif(iboux3ot(m) .eq. 5 ) then

        do j = js1, je1
          do i = is1, ie1
            do n = nm1-1, 0, -1
              uu(m,i,j,kmax-n)=uu(m,i,j,kmax-nm1)
            enddo
          enddo
        enddo
 
      elseif(iboux3ot(m) .eq. 12 ) then

        do j = js1, je1
          do i = is1, ie1
            do n = 0, nm1-1
!              uu(m,i,j,kmax-n)=-uu(m,i,j,kmax-nm1-n)
              uu(m,i,j,kmax-n)=-uu(m,i,j,kmax-2*nm1+1+n)
            enddo
          enddo
        enddo

      elseif(iboux3ot(m) .eq. 13 ) then

        do j = js1, je1
          do i = is1, ie1
            do n = 0, nm1-1
              uu(m,i,j,kmax-n)=0.d0
            enddo
          enddo
        enddo
!
      elseif( iboux3ot(m) .eq. 14 ) then
 
        do j = js1, je1
          do i = is1, ie1
            do n = 0, nm1-1
              uu(m,i,j,kmax-n)=uu(m,i,j,kmax-2*nm1+1+n)
            enddo
          enddo
        enddo
!
      elseif(iboux3ot(m) .eq. 7 ) then

        do j = js1, je1
          do i = is1, ie1
            do n = nm1-1, 0, -1
              uu(m,i,j,kmax-n)=uo(m,i,j,kmax-n)+uu(m,i,j,kmax-n-1) &
                             -uo(m,i,j,kmax-n-1)
            enddo
          enddo
        enddo
 
      else
!        write(4,*) 'outer boundary setting error for x3 axis'
!        write(4,*) 'stop in bndry'
      endif
!
    endif
!
  enddo
!
  return
end subroutine bndryfb4
! 
!--------------------------------------------------------------------
subroutine calconv(uu,uri,gcov,gcon,detg,x1,nm1, &
                   is1,ie1,js1,je1,ks1,ke1,myranki,myrankj,myrankk)
!--------------------------------------------------------------------
!
  use pram, only : imax, jmax, kmax, nv, iprocs, jprocs, kprocs, &
                   gam, c0, ieos, iflag
  implicit none

  integer :: i, j, k, n, m, nn, nm1, n2
  integer :: is1, ie1, js1, je1, ks1, ke1, myranki, myrankj, myrankk 
  integer :: merr

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
 
  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)

  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: delta(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:)
  real(8), allocatable :: util(:), bcon(:), vcon(:)  

  real(8) :: x1(imax)

  real(8) :: de,  pr, roh, roe, gfl, bbsq
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), delta(0:3,0:3), stat=merr)
  allocate(util(1:3), bcon(1:3), vcon(1:3), stat=merr)
  allocate(ucov(0:3), ucon(0:3), bbcov(0:3), bbcon(0:3), stat=merr)
!
!=====================================================================
!
!!!!
  do k=ks1+nm1,ke1-nm1
    do i=is1+nm1,ie1-nm1
      do n=1,nm1

        if(myrankj .eq. 0 .or. jprocs .eq. 1) then
!
!- copy of primitive variables -!
          de=uri(1,i,n,k)
          util(1)=uri(2,i,n,k)
          util(2)=uri(3,i,n,k)
          util(3)=uri(4,i,n,k)
!          vcon(1)=uri(2,i,n,k)
!          vcon(2)=uri(3,i,n,k)
!          vcon(3)=uri(4,i,n,k)
          pr=uri(5,i,n,k)
          bcon(1)=uri(7,i,n,k)
          bcon(2)=uri(8,i,n,k)
          bcon(3)=uri(9,i,n,k)
!- calculation of relativistic enthalpy -!  
          if(ieos .eq. 0 .or. ieos .eq. 3) then
            roh=de+(gam/(gam-1.0))*pr
          elseif(ieos .eq. 1) then
            roh=(5./2.)*pr+sqrt((9./4.)*pr**2+de**2)
          elseif(ieos .eq. 2) then
            roe=(3./2.)*(pr+((3.*pr**2) &
                /(2.0*de+sqrt(2.*pr**2+4.*de**2)) ))
            roh=de+roe+pr
          endif
!- copy of metric terms -!
          do m=0,3
            do nn=0,3
              gcov1(m,nn)=gcov(m,nn,i,n,k)
              gcon1(m,nn)=gcon(m,nn,i,n,k)
            enddo
          enddo 
!
!- cal of lorentz factor -!        
          call calgfl(util,gfl,gcov1) 
!- cal of contravariant 4-velocity -!
          call calucon(util,ucon,gfl,gcon1)
!- cal of covariant 4-velocity -! 
          call lower(ucon,ucov,gcov1)
!- cal of contravariant 4-magnetic field -!
          call calbbcon(bcon,bbcon,ucov,ucon)
!- cal of covariant 4-magnetic field -! 
          call lower(bbcon,bbcov,gcov1) 
!- cal of 4-magnetic field square -! 
          call calbbsq(bbcon,bbcov,bbsq)
!- cal of delta function -!
          call caldelta(delta,gcov1,gcon1) 
!
!- calculation of conserved variables -!
!         
          uu(1,i,n,k)=de*ucon(0) !- relativistic mass -!
          uu(2,i,n,k)=(roh+bbsq)*ucon(0)*ucov(1)+(pr+0.5*bbsq)*delta(0,1) &
                     -bbcon(0)*bbcov(1) !- momentum density in i-direc -!
          uu(3,i,n,k)=(roh+bbsq)*ucon(0)*ucov(2)+(pr+0.5*bbsq)*delta(0,2) &
                     -bbcon(0)*bbcov(2) !- momentum density in j-direc -!
          uu(4,i,n,k)=(roh+bbsq)*ucon(0)*ucov(3)+(pr+0.5*bbsq)*delta(0,3) &
                     -bbcon(0)*bbcov(3) !- momentum density in k-direc -!
          uu(5,i,n,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
                     -bbcon(0)*bbcov(0)-de*ucon(0) !- total energy -!
!          uu(5,i,n,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
!                     -bbcon(0)*bbcov(0)+de*ucon(0) !- total energy -!
!          uu(5,i,n,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
!                     -bbcon(0)*bbcov(0) !- total energy -!
           
          if( iflag(6).le.1 ) then
            uu(6,i,n,k)=gfl*de*x1(i)
          elseif( iflag(6).le.3 ) then
            uu(6,i,n,k)=gfl*pr/de**(gam-1.0)
          elseif( iflag(6).eq.4 ) then
            uu(6,i,n,k)=gfl*de*log(pr/de**gam)
          endif         
           
          uu(7,i,n,k)=bcon(1) !- contravariant B-field in i-direc -!
          uu(8,i,n,k)=bcon(2) !- contravariant B-field in j-direc -!
          uu(9,i,n,k)=bcon(3) !- contravariant B-field in k-direc -!
!
          do m=1,9
            uu(m,i,n,k)=detg(i,n,k)*uu(m,i,n,k)
          enddo
!
        endif
!!!!
        if(myrankj .eq. jprocs-1 .or. jprocs .eq. 1) then
          n2=n-1
!
!- copy of primitive variables -!
          de=uri(1,i,jmax-n2,k)
          util(1)=uri(2,i,jmax-n2,k)
          util(2)=uri(3,i,jmax-n2,k)
          util(3)=uri(4,i,jmax-n2,k)
!          vcon(1)=uri(2,i,jmax-n2,k)
!          vcon(2)=uri(3,i,jmax-n2,k)
!          vcon(3)=uri(4,i,jmax-n2,k)
          pr=uri(5,i,jmax-n2,k)
          bcon(1)=uri(7,i,jmax-n2,k)
          bcon(2)=uri(8,i,jmax-n2,k)
          bcon(3)=uri(9,i,jmax-n2,k)
!- calculation of relativistic enthalpy -!  
          if(ieos .eq. 0 .or. ieos .eq. 3) then
            roh=de+(gam/(gam-1.0))*pr
          elseif(ieos .eq. 1) then
            roh=(5./2.)*pr+sqrt((9./4.)*pr**2+de**2)
          elseif(ieos .eq. 2) then
            roe=(3./2.)*(pr+((3.*pr**2) &
                  /(2.0*de+sqrt(2.*pr**2+4.*de**2)) ))
            roh=de+roe+pr
          endif
!- copy of metric terms -!
          do m=0,3
            do nn=0,3
              gcov1(m,nn)=gcov(m,nn,i,jmax-n2,k)
              gcon1(m,nn)=gcon(m,nn,i,jmax-n2,k)
            enddo
          enddo
!
!- cal of lorentz factor -!        
          call calgfl(util,gfl,gcov1) 
!- cal of contravariant 4-velocity -!
          call calucon(util,ucon,gfl,gcon1)
!- cal of covariant 4-velocity -! 
          call lower(ucon,ucov,gcov1)
!- cal of contravariant 4-magnetic field -!
          call calbbcon(bcon,bbcon,ucov,ucon)
!- cal of covariant 4-magnetic field -! 
          call lower(bbcon,bbcov,gcov1) 
!- cal of 4-magnetic field square -! 
          call calbbsq(bbcon,bbcov,bbsq)
!- cal of delta function -!
          call caldelta(delta,gcov1,gcon1) 
!
!- calculation of conserved variables -!
!
          uu(1,i,jmax-n2,k)=de*ucon(0) !- relativistic mass -!
          uu(2,i,jmax-n2,k)=(roh+bbsq)*ucon(0)*ucov(1) &
                           +(pr+0.5*bbsq)*delta(0,1) &
                           -bbcon(0)*bbcov(1) !- momentum density in i-direc -!
          uu(3,i,jmax-n2,k)=(roh+bbsq)*ucon(0)*ucov(2) &
                           +(pr+0.5*bbsq)*delta(0,2) &
                           -bbcon(0)*bbcov(2) !- momentum density in j-direc -!
          uu(4,i,jmax-n2,k)=(roh+bbsq)*ucon(0)*ucov(3) &
                           +(pr+0.5*bbsq)*delta(0,3) &
                           -bbcon(0)*bbcov(3) !- momentum density in k-direc -!
          uu(5,i,jmax-n2,k)=(roh+bbsq)*ucon(0)*ucov(0) &
                           +(pr+0.5*bbsq)*delta(0,0) &
                           -bbcon(0)*bbcov(0)-de*ucon(0) !- total energy -!
!          uu(5,i,jmax-n2,k)=(roh+bbsq)*ucon(0)*ucov(0) &
!                           +(pr+0.5*bbsq)*delta(0,0) &
!                           -bbcon(0)*bbcov(0)+de*ucon(0) !- total energy -!
!          uu(5,i,jmax-n2,k)=(roh+bbsq)*ucon(0)*ucov(0) &
!                           +(pr+0.5*bbsq)*delta(0,0) &
!                           -bbcon(0)*bbcov(0) !- total energy -!
           
          if( iflag(6).le.1 ) then
            uu(6,i,jmax-n2,k)=gfl*de*x1(i)
          elseif( iflag(6).le.3 ) then
            uu(6,i,jmax-n2,k)=gfl*pr/de**(gam-1.0)
          elseif( iflag(6).eq.4 ) then
            uu(6,i,jmax-n2,k)=gfl*de*log(pr/de**gam)
          endif         
           
          uu(7,i,jmax-n2,k)=bcon(1) !- contravariant B-field in i-direc -!
          uu(8,i,jmax-n2,k)=bcon(2) !- contravariant B-field in j-direc -!
          uu(9,i,jmax-n2,k)=bcon(3) !- contravariant B-field in k-direc -!
!
          do m=1,9
            uu(m,i,jmax-n2,k)=detg(i,jmax-n2,k)*uu(m,i,jmax-n2,k)
          enddo
!
        endif
      enddo
    enddo
  enddo
!
!!
!
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
      do n=1,nm1
        
        if(myranki .eq. 0 .or. iprocs .eq. 1) then
!
!- copy of primitive variables -!
          de=uri(1,n,j,k)
          util(1)=uri(2,n,j,k)
          util(2)=uri(3,n,j,k)
          util(3)=uri(4,n,j,k)
!          vcon(1)=uri(2,n,j,k)
!          vcon(2)=uri(3,n,j,k)
!          vcon(3)=uri(4,n,j,k)
          pr=uri(5,n,j,k)
          bcon(1)=uri(7,n,j,k)
          bcon(2)=uri(8,n,j,k)
          bcon(3)=uri(9,n,j,k)
!- calculation of relativistic enthalpy -!  
          if(ieos .eq. 0 .or. ieos .eq. 3) then
            roh=de+(gam/(gam-1.0))*pr
          elseif(ieos .eq. 1) then
            roh=(5./2.)*pr+sqrt((9./4.)*pr**2+de**2)
          elseif(ieos .eq. 2) then
            roe=(3./2.)*(pr+((3.*pr**2) &
                /(2.0*de+sqrt(2.*pr**2+4.*de**2)) ))
            roh=de+roe+pr
          endif
!- copy of metric terms -!
          do m=0,3
            do nn=0,3
              gcov1(m,nn)=gcov(m,nn,n,j,k)
              gcon1(m,nn)=gcon(m,nn,n,j,k)
            enddo
          enddo
!
!- cal of lorentz factor -!        
          call calgfl(util,gfl,gcov1) 
!- cal of contravariant 4-velocity -!
          call calucon(util,ucon,gfl,gcon1)
!- cal of covariant 4-velocity -! 
          call lower(ucon,ucov,gcov1)
!- cal of contravariant 4-magnetic field -!
          call calbbcon(bcon,bbcon,ucov,ucon)
!- cal of covariant 4-magnetic field -! 
          call lower(bbcon,bbcov,gcov1) 
!- cal of 4-magnetic field square -! 
          call calbbsq(bbcon,bbcov,bbsq)
!- cal of delta function -!
          call caldelta(delta,gcov1,gcon1) 
!
!- calculation of conserved variables -!
!         
          uu(1,n,j,k)=de*ucon(0) !- relativistic mass -!
          uu(2,n,j,k)=(roh+bbsq)*ucon(0)*ucov(1)+(pr+0.5*bbsq)*delta(0,1) &
                     -bbcon(0)*bbcov(1) !- momentum density in i-direc -!
          uu(3,n,j,k)=(roh+bbsq)*ucon(0)*ucov(2)+(pr+0.5*bbsq)*delta(0,2) &
                     -bbcon(0)*bbcov(2) !- momentum density in j-direc -!
          uu(4,n,j,k)=(roh+bbsq)*ucon(0)*ucov(3)+(pr+0.5*bbsq)*delta(0,3) &
                     -bbcon(0)*bbcov(3) !- momentum density in k-direc -!
          uu(5,n,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
                     -bbcon(0)*bbcov(0)-de*ucon(0) !- total energy -!
!          uu(5,n,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
!                     -bbcon(0)*bbcov(0)+de*ucon(0) !- total energy -!
!          uu(5,n,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
!                     -bbcon(0)*bbcov(0) !- total energy -!
           
          if( iflag(6).le.1 ) then
            uu(6,n,j,k)=gfl*de*x1(n)
          elseif( iflag(6).le.3 ) then
            uu(6,n,j,k)=gfl*pr/de**(gam-1.0)
          elseif( iflag(6).eq.4 ) then
            uu(6,n,j,k)=gfl*de*log(pr/de**gam)
          endif         
           
          uu(7,n,j,k)=bcon(1) !- contravariant B-field in i-direc -!
          uu(8,n,j,k)=bcon(2) !- contravariant B-field in j-direc -!
          uu(9,n,j,k)=bcon(3) !- contravariant B-field in k-direc -!
!
          do m=1,9
            uu(m,n,j,k)=detg(n,j,k)*uu(m,n,j,k)
          enddo
!
        endif
!!!!
        if(myranki .eq. iprocs-1 .or. iprocs .eq. 1) then
          n2=n-1
!
!- copy of primitive variables -!
          de=uri(1,imax-n2,j,k)
          util(1)=uri(2,imax-n2,j,k)
          util(2)=uri(3,imax-n2,j,k)
          util(3)=uri(4,imax-n2,j,k)
!          vcon(1)=uri(2,imax-n2,j,k)
!          vcon(2)=uri(3,imax-n2,j,k)
!          vcon(3)=uri(4,imax-n2,j,k)
          pr=uri(5,imax-n2,j,k)
          bcon(1)=uri(7,imax-n2,j,k)
          bcon(2)=uri(8,imax-n2,j,k)
          bcon(3)=uri(9,imax-n2,j,k)
!- calculation of relativistic enthalpy -!  
          if(ieos .eq. 0 .or. ieos .eq. 3) then
            roh=de+(gam/(gam-1.0))*pr
          elseif(ieos .eq. 1) then
            roh=(5./2.)*pr+sqrt((9./4.)*pr**2+de**2)
          elseif(ieos .eq. 2) then
            roe=(3./2.)*(pr+((3.*pr**2) &
                  /(2.0*de+sqrt(2.*pr**2+4.*de**2)) ))
            roh=de+roe+pr
          endif
!- copy of metric terms -!
          do m=0,3
            do nn=0,3
              gcov1(m,nn)=gcov(m,nn,imax-n2,j,k)
              gcon1(m,nn)=gcon(m,nn,imax-n2,j,k)
            enddo
          enddo
!
!- cal of lorentz factor -!        
          call calgfl(util,gfl,gcov1) 
!- cal of contravariant 4-velocity -!
          call calucon(util,ucon,gfl,gcon1)
!- cal of covariant 4-velocity -! 
          call lower(ucon,ucov,gcov1)
!- cal of contravariant 4-magnetic field -!
          call calbbcon(bcon,bbcon,ucov,ucon)
!- cal of covariant 4-magnetic field -! 
          call lower(bbcon,bbcov,gcov1) 
!- cal of 4-magnetic field square -! 
          call calbbsq(bbcon,bbcov,bbsq)
!- cal of delta function -!
          call caldelta(delta,gcov1,gcon1) 
!
!- calculation of conserved variables -!
!
          uu(1,imax-n2,j,k)=de*ucon(0) !- relativistic mass -!
          uu(2,imax-n2,j,k)=(roh+bbsq)*ucon(0)*ucov(1) &
                           +(pr+0.5*bbsq)*delta(0,1) &
                           -bbcon(0)*bbcov(1) !- momentum density in i-direc -!
          uu(3,imax-n2,j,k)=(roh+bbsq)*ucon(0)*ucov(2) &
                           +(pr+0.5*bbsq)*delta(0,2) &
                           -bbcon(0)*bbcov(2) !- momentum density in j-direc -!
          uu(4,imax-n2,j,k)=(roh+bbsq)*ucon(0)*ucov(3) &
                           +(pr+0.5*bbsq)*delta(0,3) &
                           -bbcon(0)*bbcov(3) !- momentum density in k-direc -!
          uu(5,imax-n2,j,k)=(roh+bbsq)*ucon(0)*ucov(0) &
                           +(pr+0.5*bbsq)*delta(0,0) &
                           -bbcon(0)*bbcov(0)-de*ucon(0) !- total energy -!
!          uu(5,imax-n2,j,k)=(roh+bbsq)*ucon(0)*ucov(0) &
!                           +(pr+0.5*bbsq)*delta(0,0) &
!                           -bbcon(0)*bbcov(0)+de*ucon(0) !- total energy -!
!          uu(5,imax-n2,j,k)=(roh+bbsq)*ucon(0)*ucov(0) &
!                           +(pr+0.5*bbsq)*delta(0,0) &
!                           -bbcon(0)*bbcov(0) !- total energy -!
           
          if( iflag(6).le.1 ) then
            uu(6,imax-n2,j,k)=gfl*de*x1(imax-n2)
          elseif( iflag(6).le.3 ) then
            uu(6,imax-n2,j,k)=gfl*pr/de**(gam-1.0)
          elseif( iflag(6).eq.4 ) then
            uu(6,imax-n2,j,k)=gfl*de*log(pr/de**gam)
          endif         
           
          uu(7,imax-n2,j,k)=bcon(1) !- contravariant B-field in i-direc -!
          uu(8,imax-n2,j,k)=bcon(2) !- contravariant B-field in j-direc -!
          uu(9,imax-n2,j,k)=bcon(3) !- contravariant B-field in k-direc -!
!
          do m=1,9
            uu(m,imax-n2,j,k)=detg(imax-n2,j,k)*uu(m,imax-n2,j,k)
          enddo
!
        endif
      enddo
    enddo
  enddo
!
!!!
!
  do j=js1+nm1,je1-nm1
    do i=is1+nm1,ie1-nm1
      do n=1,nm1

        if(myrankk .eq. 0 .or. kprocs .eq. 1) then
!- copy of primitive variables -!
          de=uri(1,i,j,n)
          util(1)=uri(2,i,j,n)
          util(2)=uri(3,i,j,n)
          util(3)=uri(4,i,j,n)
!          vcon(1)=uri(2,i,j,n)
!          vcon(2)=uri(3,i,j,n)
!          vcon(3)=uri(4,i,j,n)
          pr=uri(5,i,j,n)
          bcon(1)=uri(7,i,j,n)
          bcon(2)=uri(8,i,j,n)
          bcon(3)=uri(9,i,j,n)
!- calculation of relativistic enthalpy -!           
          if(ieos .eq. 0 .or. ieos .eq. 3) then
            roh=de+(gam/(gam-1.0))*pr
          elseif(ieos .eq. 1) then
            roh=(5./2.)*pr+sqrt((9./4.)*pr**2+de**2)
          elseif(ieos .eq. 2) then
            roe=(3./2.)*(pr+((3.*pr**2) &
                /(2.0*de+sqrt(2.*pr**2+4.*de**2)) ))
            roh=de+roe+pr
          endif
!- copy of metric terms -!
          do m=0,3
            do nn=0,3
              gcov1(m,nn)=gcov(m,nn,i,j,n)
              gcon1(m,nn)=gcon(m,nn,i,j,n)
            enddo
          enddo
!
!- cal of lorentz factor -!        
          call calgfl(util,gfl,gcov1) 
!- cal of contravariant 4-velocity -!
          call calucon(util,ucon,gfl,gcon1)
!- cal of covariant 4-velocity -! 
          call lower(ucon,ucov,gcov1)
!- cal of contravariant 4-magnetic field -!
          call calbbcon(bcon,bbcon,ucov,ucon)
!- cal of covariant 4-magnetic field -! 
          call lower(bbcon,bbcov,gcov1) 
!- cal of 4-magnetic field square -! 
          call calbbsq(bbcon,bbcov,bbsq)
!- cal of delta function -!
          call caldelta(delta,gcov1,gcon1) 
!
!- calculation of conserved variables -!
!         
          uu(1,i,j,n)=de*ucon(0) !- relativistic mass -!
          uu(2,i,j,n)=(roh+bbsq)*ucon(0)*ucov(1)+(pr+0.5*bbsq)*delta(0,1) &
                     -bbcon(0)*bbcov(1) !- momentum density in i-direc -!
          uu(3,i,j,n)=(roh+bbsq)*ucon(0)*ucov(2)+(pr+0.5*bbsq)*delta(0,2) &
                     -bbcon(0)*bbcov(2) !- momentum density in j-direc -!
          uu(4,i,j,n)=(roh+bbsq)*ucon(0)*ucov(3)+(pr+0.5*bbsq)*delta(0,3) &
                     -bbcon(0)*bbcov(3) !- momentum density in k-direc -!
          uu(5,i,j,n)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
                     -bbcon(0)*bbcov(0)-de*ucon(0) !- total energy -!
!          uu(5,i,j,n)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
!                     -bbcon(0)*bbcov(0)+de*ucon(0) !- total energy -!
!          uu(5,i,j,n)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
!                     -bbcon(0)*bbcov(0) !- total energy -!
           
          if( iflag(6).le.1 ) then
            uu(6,i,j,n)=gfl*de*x1(i)
          elseif( iflag(6).le.3 ) then
            uu(6,i,j,n)=gfl*pr/de**(gam-1.0)
          elseif( iflag(6).eq.4 ) then
            uu(6,i,j,n)=gfl*de*log(pr/de**gam)
          endif         
           
          uu(7,i,j,n)=bcon(1) !- contravariant B-field in i-direc -!
          uu(8,i,j,n)=bcon(2) !- contravariant B-field in j-direc -!
          uu(9,i,j,n)=bcon(3) !- contravariant B-field in k-direc -!
!
          do m=1,9
            uu(m,i,j,n)=detg(i,j,n)*uu(m,i,j,n)
          enddo
!
        endif
!!!!
        if(myrankk .eq. kprocs-1 .or. kprocs .eq. 1) then
          n2=n-1
!
!- copy of primitive variables -! 
          de=uri(1,i,j,kmax-n2)
          util(1)=uri(2,i,j,kmax-n2)
          util(2)=uri(3,i,j,kmax-n2)
          util(3)=uri(4,i,j,kmax-n2)
!          vcon(1)=uri(2,i,j,kmax-n2)
!          vcon(2)=uri(3,i,j,kmax-n2)
!          vcon(3)=uri(4,i,j,kmax-n2)
          pr=uri(5,i,j,kmax-n2)
          bcon(1)=uri(7,i,j,kmax-n2)
          bcon(2)=uri(8,i,j,kmax-n2)
          bcon(3)=uri(9,i,j,kmax-n2)

          if(ieos .eq. 0 .or. ieos .eq. 3) then
            roh=de+(gam/(gam-1.0))*pr
          elseif(ieos .eq. 1) then
           roh=(5./2.)*pr+sqrt((9./4.)*pr**2+de**2)
          elseif(ieos .eq. 2) then
            roe=(3./2.)*(pr+((3.*pr**2) &
                /(2.0*de+sqrt(2.*pr**2+4.*de**2)) ))
            roh=de+roe+pr
          endif
!- copy of metric terms -!
          do m=0,3
            do nn=0,3
              gcov1(m,nn)=gcov(m,nn,i,j,kmax-n2)
              gcon1(m,nn)=gcon(m,nn,i,j,kmax-n2)
            enddo
          enddo
!
!- cal of lorentz factor -!        
          call calgfl(util,gfl,gcov1) 
!- cal of contravariant 4-velocity -!
          call calucon(util,ucon,gfl,gcon1)
!- cal of covariant 4-velocity -! 
          call lower(ucon,ucov,gcov1)
!- cal of contravariant 4-magnetic field -!
          call calbbcon(bcon,bbcon,ucov,ucon)
!- cal of covariant 4-magnetic field -! 
          call lower(bbcon,bbcov,gcov1) 
!- cal of 4-magnetic field square -! 
          call calbbsq(bbcon,bbcov,bbsq)
!- cal of delta function -!
          call caldelta(delta,gcov1,gcon1) 
!
!- calculation of conserved variables -!
! 
          uu(1,i,j,kmax-n2)=de*ucon(0) !- relativistic mass -!
          uu(2,i,j,kmax-n2)=(roh+bbsq)*ucon(0)*ucov(1) &
                           +(pr+0.5*bbsq)*delta(0,1) &
                           -bbcon(0)*bbcov(1) !- momentum density in i-direc -!
          uu(3,i,j,kmax-n2)=(roh+bbsq)*ucon(0)*ucov(2) &
                           +(pr+0.5*bbsq)*delta(0,2) &
                           -bbcon(0)*bbcov(2) !- momentum density in j-direc -!
          uu(4,i,j,kmax-n2)=(roh+bbsq)*ucon(0)*ucov(3) &
                           +(pr+0.5*bbsq)*delta(0,3) &
                           -bbcon(0)*bbcov(3) !- momentum density in k-direc -!
          uu(5,i,j,kmax-n2)=(roh+bbsq)*ucon(0)*ucov(0) &
                           +(pr+0.5*bbsq)*delta(0,0) &
                           -bbcon(0)*bbcov(0)-de*ucon(0) !- total energy -!
!          uu(5,i,j,kmax-n2)=(roh+bbsq)*ucon(0)*ucov(0) &
!                           +(pr+0.5*bbsq)*delta(0,0) &
!                           -bbcon(0)*bbcov(0)+de*ucon(0) !- total energy -!
!          uu(5,i,j,kmax-n2)=(roh+bbsq)*ucon(0)*ucov(0) &
!                           +(pr+0.5*bbsq)*delta(0,0) &
!                           -bbcon(0)*bbcov(0) !- total energy -!
           
          if( iflag(6).le.1 ) then
            uu(6,i,j,kmax-n2)=gfl*de*x1(i)
          elseif( iflag(6).le.3 ) then
            uu(6,i,j,kmax-n2)=gfl*pr/de**(gam-1.0)
          elseif( iflag(6).eq.4 ) then
            uu(6,i,j,kmax-n2)=gfl*de*log(pr/de**gam)
          endif         
           
          uu(7,i,j,kmax-n2)=bcon(1) !- contravariant B-field in i-direc -!
          uu(8,i,j,kmax-n2)=bcon(2) !- contravariant B-field in j-direc -!
          uu(9,i,j,kmax-n2)=bcon(3) !- contravariant B-field in k-direc -!     
!
          do m=1,9
            uu(m,i,j,kmax-n2)=detg(i,j,kmax-n2)*uu(m,i,j,kmax-n2)
          enddo
!
        endif
      enddo
    enddo
  enddo
!
  deallocate( gcov1, gcon1, delta, util, bcon, vcon, ucov, ucon, &
              bbcov, bbcon, stat=merr)
!
  return
end subroutine calconv
!
!--------------------------------------------------------------------
subroutine infchk(uri,gcov,gcon,nm1,is1,ie1,js1,je1,ks1,ke1, &
                  myranki)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv, iprocs, metric

  implicit none 

  integer :: i, j, k, n, m, nn, nm1, n2
  integer :: is1, ie1, js1, je1, ks1, ke1, myranki
  integer :: merr

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: util(:), ucon(:), vel(:), beta1(:)  
  real(8) :: alpha1, gfl, gflmax, vsq
!
!-allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), stat=merr)
  allocate(util(1:3), ucon(0:3), vel(1:3), beta1(1:3), stat=merr)

!--------------------------------------------------------------------
!- Parameter -!
!
  gflmax=10.d0
!
!- inflow check from inner & outer radial boundary -!
   

  if(myranki .eq. 0 .or. iprocs .eq. 1) then
     
    do k=ks1,ke1
      do j=js1,je1
        do n=1,nm1   
!
!- copy of primitive variables -!
          util(1)=uri(2,n,j,k)
          util(2)=uri(3,n,j,k)
          util(3)=uri(4,n,j,k)
!
!- copy of metric terms -!
          do m=0,3
            do nn=0,3
              gcov1(m,nn)=gcov(m,nn,n,j,k)
              gcon1(m,nn)=gcon(m,nn,n,j,k)
            enddo
          enddo
!
!- cal of lorentz factor -!        
          call calgfl(util,gfl,gcov1) 
!- cal of contravariant 4-velocity -!
          call calucon(util,ucon,gfl,gcon1)

!- 4-velocity=> 3-velocity -!
          vel(1)=util(1)*(1./gfl)
          vel(2)=util(2)*(1./gfl)
          vel(3)=util(3)*(1./gfl)
!         
          if(vel(1) .gt. 0.d0) then
!          if(ucon(1) .gt. 0.d0) then   
!             
!          write(*,*) 'check 1e in infchk part1 start' 
!

            alpha1=1./sqrt(-gcon1(0,0))
!
            beta1(1)=alpha1*alpha1*gcon1(0,1)
            beta1(2)=alpha1*alpha1*gcon1(0,2)
            beta1(3)=alpha1*alpha1*gcon1(0,3)
!
!- reset radial 4-velocity -!            
!            vel(1)=beta1(1)/alpha1
            vel(1)=0.d0
!            
!- new calculation of Lorentz factor -!
            vsq=0.d0
            do m=1,3
              do nn=1,3
                vsq=vsq+gcov1(m,nn)*vel(m)*vel(nn)     
              enddo
            enddo
!
!            if(abs(vsq) .le. 1.d-10) then
!              vsq=1.d-10
!            endif

            if(vsq .ge. 1.d0) then
              vsq=1.-1./(gflmax*gflmax)
            endif
!
!- new set of 4-velocity -! 
            gfl=1./sqrt(1.-vsq)
!            uri(2,n,j,k)=vel(1)*gfl
!            uri(3,n,j,k)=vel(2)*gfl
!            uri(4,n,j,k)=vel(3)*gfl
!
!            write(*,*) 'check 1e in infchk part1 OK'
!            
          endif
       !
        enddo
      enddo       
    enddo
!   
  endif
!!!!
  if(myranki .eq. iprocs-1 .or. iprocs .eq. 1) then
     
    do k=ks1,ke1
      do j=js1,je1
        do n=0,nm1-1
            
!          n2=n-1
!
!- copy of primitive variables -!
          util(1)=uri(2,imax-n,j,k)
          util(2)=uri(3,imax-n,j,k)
          util(3)=uri(4,imax-n,j,k)
!
!- copy of metric terms -!
          do m=0,3
            do nn=0,3
              gcov1(m,nn)=gcov(m,nn,imax-n,j,k)
              gcon1(m,nn)=gcon(m,nn,imax-n,j,k)
            enddo
          enddo
!
!- cal of lorentz factor -!        
          call calgfl(util,gfl,gcov1) 
!- cal of contravariant 4-velocity -!
          call calucon(util,ucon,gfl,gcon1)
!
!- 4-velocity=> 3-velocity -!
          vel(1)=util(1)*(1./gfl)
          vel(2)=util(2)*(1./gfl)
          vel(3)=util(3)*(1./gfl)
          
          if(vel(1) .lt. 0.d0) then
!          if(ucon(1) .lt. 0.d0) then   
!
!            write(*,*) 'check 1f in infchk part2 start'                
!
            alpha1=1./sqrt(-gcon1(0,0))
!
            beta1(1)=alpha1*alpha1*gcon1(0,1)
            beta1(2)=alpha1*alpha1*gcon1(0,2)
            beta1(3)=alpha1*alpha1*gcon1(0,3)
!
!- reset radial 4-velocity -!            
            vel(1)=beta1(1)*(1./alpha1)
!            vel(1)=0.d0
!            
!- new calculation of Lorentz factor -!
            vsq=0.d0
            do m=1,3
              do nn=1,3
                vsq=vsq+gcov1(m,nn)*vel(m)*vel(nn)     
              enddo
            enddo
!
!            if(abs(vsq) .le. 1.d-10) then
!              vsq=1.d-10
!            endif
            if(vsq .ge. 1.d0) then
              vsq=1.-1./(gflmax*gflmax)
            endif
!
!- new set of 4-velocity -! 
            gfl=1./sqrt(1.-vsq)
            uri(2,imax-n,j,k)=vel(1)*gfl
            uri(3,imax-n,j,k)=vel(2)*gfl
            uri(4,imax-n,j,k)=vel(3)*gfl
!
!            write(*,*) 'check 1f in infchk part2 OK', n, k               
          endif
!
        enddo
      enddo
    enddo
!
  endif
! 
  deallocate(gcov1, gcon1, util, ucon, vel, beta1, stat=merr)
!
  return
end subroutine infchk

