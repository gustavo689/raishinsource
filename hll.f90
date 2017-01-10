!---------------------------------------------------------------------@
subroutine hll(ww,uuir,uujr,uukr,uuil,uujl,uukl, &
               wwir,wwjr,wwkr,wwil,wwjl,wwkl, &
               cmaxi,cmini,cmaxj,cminj,cmaxk,cmink,nm0, &
               is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Calculate numerical flux by Approximate Riemann Solver 
!
  use pram, only : imax, jmax, kmax, nv, ihll
  implicit none

  integer :: nm0, is1, ie1, js1, je1, ks1, ke1

  real(8) :: uuir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uuil(nv,is1:ie1,js1:je1,ks1:ke1), &
             uujr(nv,is1:ie1,js1:je1,ks1:ke1), &
             uujl(nv,is1:ie1,js1:je1,ks1:ke1), &
             uukr(nv,is1:ie1,js1:je1,ks1:ke1), &
             uukl(nv,is1:ie1,js1:je1,ks1:ke1)
     
  real(8) :: ww(3,nv,is1:ie1,js1:je1,ks1:ke1), &
             wwir(nv,is1:ie1,js1:je1,ks1:ke1), &
             wwil(nv,is1:ie1,js1:je1,ks1:ke1), &
             wwjr(nv,is1:ie1,js1:je1,ks1:ke1), &
             wwjl(nv,is1:ie1,js1:je1,ks1:ke1), &
             wwkr(nv,is1:ie1,js1:je1,ks1:ke1), &
             wwkl(nv,is1:ie1,js1:je1,ks1:ke1)
      
  real(8) :: cmaxi(is1:ie1,js1:je1,ks1:ke1), cmini(is1:ie1,js1:je1,ks1:ke1), &
             cmaxj(is1:ie1,js1:je1,ks1:ke1), cminj(is1:ie1,js1:je1,ks1:ke1), &
             cmaxk(is1:ie1,js1:je1,ks1:ke1), cmink(is1:ie1,js1:je1,ks1:ke1)
!
!-------------------
      
  if(ihll .eq. 1) then
    call hlle(ww,uuir,uujr,uukr,uuil,uujl,uukl, &
              wwir,wwjr,wwkr,wwil,wwjl,wwkl, &
              cmaxi,cmini,cmaxj,cminj,cmaxk,cmink,nm0, &
              is1,ie1,js1,je1,ks1,ke1)
  endif  

  return
end subroutine hll  
!
!---------------------------------------------------------------------@
subroutine hlle(ww,uuir,uujr,uukr,uuil,uujl,uukl, &
                wwir,wwjr,wwkr,wwil,wwjl,wwkl, &
                cmaxi,cmini,cmaxj,cminj,cmaxk,cmink,nm0, &
                is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Calculate numerical flux by HLLE method
!
  use pram, only : imax, jmax, kmax, nv
  implicit none

  integer :: i, j, k, n
  integer :: nm0, is1, ie1, js1, je1, ks1, ke1

  real(8) :: uuir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uuil(nv,is1:ie1,js1:je1,ks1:ke1), & 
             uujr(nv,is1:ie1,js1:je1,ks1:ke1), &
             uujl(nv,is1:ie1,js1:je1,ks1:ke1), &
             uukr(nv,is1:ie1,js1:je1,ks1:ke1), &
             uukl(nv,is1:ie1,js1:je1,ks1:ke1)
     
  real(8) :: ww(3,nv,is1:ie1,js1:je1,ks1:ke1), &
             wwir(nv,is1:ie1,js1:je1,ks1:ke1), &
             wwil(nv,is1:ie1,js1:je1,ks1:ke1), &
             wwjr(nv,is1:ie1,js1:je1,ks1:ke1), &
             wwjl(nv,is1:ie1,js1:je1,ks1:ke1), &
             wwkr(nv,is1:ie1,js1:je1,ks1:ke1), &
             wwkl(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: cmaxi(is1:ie1,js1:je1,ks1:ke1), cmini(is1:ie1,js1:je1,ks1:ke1), &
             cmaxj(is1:ie1,js1:je1,ks1:ke1), cminj(is1:ie1,js1:je1,ks1:ke1), &
             cmaxk(is1:ie1,js1:je1,ks1:ke1), cmink(is1:ie1,js1:je1,ks1:ke1)
  real(8) :: tmp1i, tmp1j, tmp1k, tmp2i, tmp2j, tmp2k
  real(8) :: cmaxi1, cmaxj1, cmaxk1
  
  real(8), parameter :: small=1.d-12  
!
!=====================================================================@

  do k=ks1+nm0-1,ke1-nm0
    do j=js1+nm0-1,je1-nm0
      do i=is1+nm0-1,ie1-nm0
        do n=1,nv
!
!!- i-th direction          
!
          if(cmini(i,j,k) .gt. 0.d0) then
!
            ww(1,n,i,j,k)=wwil(n,i,j,k)
!
          elseif(cmini(i,j,k) .le. 0.d0 .and. cmaxi(i,j,k) .gt. 0.d0) then
!            if(n .le. 6) then
!- HLL flux
!              tmp1i=(cmaxi(i,j,k)*wwil(n,i,j,k) &
!                    -cmini(i,j,k)*wwir(n,i,j,k) &
!                    +cmaxi(i,j,k)*cmini(i,j,k) &
!                    *(uuir(n,i,j,k)-uuil(n,i,j,k)))
!  
!              ww(1,n,i,j,k)=tmp1i*(1./(cmaxi(i,j,k)-cmini(i,j,k)+small))
!- LF flux
!            elseif(n .ge. 7 .and. n .le. 9) then
!
              cmaxi1=max(abs(cmaxi(i,j,k)), abs(cmini(i,j,k)))
              ww(1,n,i,j,k)=0.5*( wwil(n,i,j,k)+wwir(n,i,j,k) &
                            -cmaxi1*(uuir(n,i,j,k)-uuil(n,i,j,k)) )
!            endif
              
          elseif(cmaxi(i,j,k) .le. 0.d0) then
            ww(1,n,i,j,k)=wwir(n,i,j,k)
          endif
!             
          if(abs(ww(1,n,i,j,k)) .lt. 1.0d-40) then
            ww(1,n,i,j,k)=0.d0
          endif
!         
!!- j-th direction       
!
          if(cminj(i,j,k) .gt. 0.d0) then
!
             ww(2,n,i,j,k)=wwjl(n,i,j,k)
!
          elseif(cminj(i,j,k) .le. 0.d0 .and. cmaxj(i,j,k) .gt. 0.d0) then
!
!            if(n .le. 6) then               
!- HLL flux
!              tmp1j=(cmaxj(i,j,k)*wwjl(n,i,j,k) &
!                    -cminj(i,j,k)*wwjr(n,i,j,k) &
!                    +cmaxj(i,j,k)*cminj(i,j,k) &
!                    *(uujr(n,i,j,k)-uujl(n,i,j,k)))              
!
!               ww(2,n,i,j,k)=tmp1j*(1./(cmaxj(i,j,k)-cminj(i,j,k)+small))   
!
!            elseif(n .ge. 7 .and. n .le. 9) then
!- LF flux
              cmaxj1=max(abs(cmaxj(i,j,k)), abs(cminj(i,j,k)))
              ww(2,n,i,j,k)=0.5*( wwjl(n,i,j,k)+wwjr(n,i,j,k) &
                            -cmaxj1*(uujr(n,i,j,k)-uujl(n,i,j,k)) )

!            endif
!           
          elseif(cmaxj(i,j,k) .le. 0.d0) then
            ww(2,n,i,j,k)=wwjr(n,i,j,k)
          endif
!          
          if(abs(ww(2,n,i,j,k)) .lt. 1.0d-40) then
            ww(2,n,i,j,k)=0.d0
          endif
!
!!- k-th direction
!                      
          if(cmink(i,j,k) .gt. 0.d0) then
!
            ww(3,n,i,j,k)=wwkl(n,i,j,k)
!
          elseif(cmink(i,j,k) .le. 0.d0 .and. cmaxk(i,j,k) .gt. 0.d0) then
!- HLL flux
!            if(n .le. 6) then
!              tmp1k=(cmaxk(i,j,k)*wwkl(n,i,j,k) &
!                    -cmink(i,j,k)*wwkr(n,i,j,k) &
!                    +cmaxk(i,j,k)*cmink(i,j,k) &
!                    *(uukr(n,i,j,k)-uukl(n,i,j,k)))   
!
!               ww(3,n,i,j,k)=tmp1k*(1./(cmaxk(i,j,k)-cmink(i,j,k)+small)) 
!
!            elseif(n .ge. 7 .and. n .le. 9) then
!- LF flux
              cmaxk1=max(abs(cmaxk(i,j,k)), abs(cmink(i,j,k)))
              ww(3,n,i,j,k)=0.5*( wwkl(n,i,j,k)+wwkr(n,i,j,k) &
                            -cmaxk1*(uukr(n,i,j,k)-uukl(n,i,j,k)) )
!            endif

          elseif(cmaxk(i,j,k) .le. 0.d0) then
            ww(3,n,i,j,k)=wwkr(n,i,j,k)
          endif
!
          if(abs(ww(3,n,i,j,k)) .lt. 1.0d-40) then
           ww(3,n,i,j,k)=0.d0
          endif
!         
        enddo
      enddo
    enddo
  enddo

  return
end subroutine hlle
!
