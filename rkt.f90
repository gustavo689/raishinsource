!----------------------------------------------------------------------
subroutine rk2fst(uh,ww,uu,nm1,akap1b,akap2b,akap3b, &
                  is1,ie1,js1,je1,ks1,ke1)
!----------------------------------------------------------------------
!
!      First step of 2nd order Runge-Kutta timeadvance step
!
!      Variables
!       uu: conserved variables
!       ww: numerical flux
!
  use pram, only : imax, jmax, kmax, nv
  implicit none
!
  integer :: i, j, k, m, nm1, is1, ie1, js1, je1, ks1, ke1

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1), &
             uh(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: ww(3,nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: akap1b(imax-1),akap2b(jmax-1),akap3b(kmax-1)

!  real(8) akap1a, akap2a, akap3a
      
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
      do i=is1+nm1,ie1-nm1
        do m=1,nv

          uh(m,i,j,k)=uu(m,i,j,k) &
           -akap1b(i)*(ww(1,m,i,j,k)-ww(1,m,i-1,j,k)) &
           -akap2b(j)*(ww(2,m,i,j,k)-ww(2,m,i,j-1,k)) &
           -akap3b(k)*(ww(3,m,i,j,k)-ww(3,m,i,j,k-1)) 

        enddo
      enddo
    enddo
  enddo
! 
  return
end subroutine rk2fst
!
!----------------------------------------------------------------------
subroutine rk2snd(uu,ww,uo,uh,nm1,akap1b,akap2b,akap3b, &
                  is1,ie1,js1,je1,ks1,ke1)
!----------------------------------------------------------------------
!
!      Second step of 2nd order Runge-Kutta timeadvance step
!
!      Variables
!       uu: conserved variables(new), 
!       uh: conserved variables(half-step)
!       uo: conserved variables(previous time)
!       ww: numerical flux
!
  use pram, only : imax, jmax, kmax, nv
  implicit none
!
  integer :: i, j, k, m, nm1, is1, ie1, js1, je1, ks1, ke1

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1), &
             uh(nv,is1:ie1,js1:je1,ks1:ke1), &
             uo(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: ww(3,nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: akap1b(imax-1),akap2b(jmax-1),akap3b(kmax-1)

!  real(8) akap1a, akap2a, akap3a
!
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
      do i=is1+nm1,ie1-nm1
        do m=1,nv
!
          uu(m,i,j,k)=(uo(m,i,j,k)+uh(m,i,j,k))*0.5+ &
           (-akap1b(i)*(ww(1,m,i,j,k)-ww(1,m,i-1,j,k)) &
            -akap2b(j)*(ww(2,m,i,j,k)-ww(2,m,i,j-1,k)) &
            -akap3b(k)*(ww(3,m,i,j,k)-ww(3,m,i,j,k-1)))*0.5
!     
        enddo
      enddo
    enddo
  enddo 
!
  return
end subroutine rk2snd
!
!----------------------------------------------------------------------
subroutine rk3snd(uu,ww,uo,uh,nm1,akap1b,akap2b,akap3b,&
                  is1,ie1,js1,je1,ks1,ke1)
!----------------------------------------------------------------------
!
!      Second step of 3rd order Runge-Kutta timeadvance step
!
!      Variables
!       uu: conserved variables(new), 
!       uh: conserved variables(half-step)
!       uo: conserved variables(previous time)
!       ww: numerical flux
!
  use pram, only : imax, jmax, kmax, nv
  implicit none
!
  integer :: i, j, k, m, nm1, is1, ie1, js1, je1, ks1, ke1

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1), &
             uh(nv,is1:ie1,js1:je1,ks1:ke1), &
             uo(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: ww(3,nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: akap1b(imax-1),akap2b(jmax-1),akap3b(kmax-1)

!  real(8) akap1a, akap2a, akap3a
!
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
      do i=is1+nm1,ie1-nm1
        do m=1,nv
!         
          uu(m,i,j,k)=(3.0*uo(m,i,j,k)+uh(m,i,j,k))*0.25+ &
            (-akap1b(i)*(ww(1,m,i,j,k)-ww(1,m,i-1,j,k)) &
             -akap2b(j)*(ww(2,m,i,j,k)-ww(2,m,i,j-1,k)) &
             -akap3b(k)*(ww(3,m,i,j,k)-ww(3,m,i,j,k-1)))*0.25
!         
        enddo
      enddo
    enddo
  enddo 
!
  return
end subroutine rk3snd
!
!----------------------------------------------------------------------
subroutine rk3trd(uu,ww,uo,us,nm1,akap1b,akap2b,akap3b,&
                  is1,ie1,js1,je1,ks1,ke1)
!----------------------------------------------------------------------
!
!      Third step of 3rd order Runge-Kutta timeadvance step
!
!      Variables
!       uu: conserved variables(new), 
!       us: conserved variables(second-step)
!       uo: conserved variables(previous time)
!       ww: numerical flux
!
  use pram, only : imax, jmax, kmax, nv
  implicit none
!
  integer :: i, j, k, m, nm1, is1, ie1, js1, je1, ks1, ke1

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1), &
             us(nv,is1:ie1,js1:je1,ks1:ke1), &
             uo(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: ww(3,nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: akap1b(imax-1),akap2b(jmax-1),akap3b(kmax-1)

!  real(8) akap1a, akap2a, akap3a
!
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
      do i=is1+nm1,ie1-nm1
        do m=1,nv

          uu(m,i,j,k)=(uo(m,i,j,k)+2.*us(m,i,j,k))*(1./3.)+ &
           (-akap1b(i)*(ww(1,m,i,j,k)-ww(1,m,i-1,j,k)) &
            -akap2b(j)*(ww(2,m,i,j,k)-ww(2,m,i,j-1,k)) &
            -akap3b(k)*(ww(3,m,i,j,k)-ww(3,m,i,j,k-1)))*(2./3.)
!    
        enddo
      enddo
    enddo
  enddo

  return
end subroutine rk3trd
!
!----------------------------------------------------------------------
subroutine rk2adsff(uu,sf,dt,nm1,is1,ie1,js1,je1,ks1,ke1)
!----------------------------------------------------------------------
!
!      Add source term to conservec variables for first step
!
!      Variables
!       uh: conserved variables
!       sf, sou: source term
!
  use pram, only : imax, jmax, kmax, nv
  implicit none
!
  integer :: i, j, k, nm1, is1, ie1, js1, je1, ks1, ke1
  
  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: sf(2:5,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: dt
!
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
      do i=is1+nm1,ie1-nm1
         
        uu(2,i,j,k)=uu(2,i,j,k)+sf(2,i,j,k)*dt
        uu(3,i,j,k)=uu(3,i,j,k)+sf(3,i,j,k)*dt
        uu(4,i,j,k)=uu(4,i,j,k)+sf(4,i,j,k)*dt
        uu(5,i,j,k)=uu(5,i,j,k)+sf(5,i,j,k)*dt
          
      enddo
    enddo
  enddo
!
  return
end subroutine rk2adsff
!
!----------------------------------------------------------------------
subroutine rk2adsfs(uu,sf,dt,nm1,is1,ie1,js1,je1,ks1,ke1)
!----------------------------------------------------------------------
!
!      Add source term to conservec variables for second step
!
!      Variables
!       uh: conserved variables
!       sf, sou: source term
!
  use pram, only : imax, jmax, kmax, nv
  implicit none
!
  integer :: i, j, k, nm1, is1, ie1, js1, je1, ks1, ke1
  
  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: sf(2:5,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: dt
!
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
      do i=is1+nm1,ie1-nm1
          
        uu(2,i,j,k)=uu(2,i,j,k)+sf(2,i,j,k)*dt*0.5
        uu(3,i,j,k)=uu(3,i,j,k)+sf(3,i,j,k)*dt*0.5
        uu(4,i,j,k)=uu(4,i,j,k)+sf(4,i,j,k)*dt*0.5
        uu(5,i,j,k)=uu(5,i,j,k)+sf(5,i,j,k)*dt*0.5

      enddo
    enddo
  enddo
!
  return
end subroutine rk2adsfs
!
!----------------------------------------------------------------------
subroutine rk3adsfs(uu,sf,dt,nm1,is1,ie1,js1,je1,ks1,ke1)
!----------------------------------------------------------------------
!
!      Add source term to conservec variables for second step
!
!      Variables
!       uh: conserved variables
!       sf, sou: source term
!
  use pram, only : imax, jmax, kmax, nv
  implicit none
!
  integer :: i, j, k, nm1, is1, ie1, js1, je1, ks1, ke1
  
  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: sf(2:5,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: dt
!
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
      do i=is1+nm1,ie1-nm1
          
        uu(2,i,j,k)=uu(2,i,j,k)+sf(2,i,j,k)*dt*0.25
        uu(3,i,j,k)=uu(3,i,j,k)+sf(3,i,j,k)*dt*0.25
        uu(4,i,j,k)=uu(4,i,j,k)+sf(4,i,j,k)*dt*0.25
        uu(5,i,j,k)=uu(5,i,j,k)+sf(5,i,j,k)*dt*0.25
     
      enddo
    enddo
  enddo
!
  return
end subroutine rk3adsfs
!
!----------------------------------------------------------------------
subroutine rk3adsft(uu,sf,dt,nm1,is1,ie1,js1,je1,ks1,ke1)
!----------------------------------------------------------------------
!
!      Add source term to conservec variables for third step
!
!      Variables
!       uh: conserved variables
!       sf, sou: source term
!
  use pram, only : imax, jmax, kmax, nv
  implicit none
!
  integer :: i, j, k, nm1, is1, ie1, js1, je1, ks1, ke1
  
  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: sf(2:5,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: dt
!
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
      do i=is1+nm1,ie1-nm1
          
        uu(2,i,j,k)=uu(2,i,j,k)+sf(2,i,j,k)*dt*(2./3.)
        uu(3,i,j,k)=uu(3,i,j,k)+sf(3,i,j,k)*dt*(2./3.)
        uu(4,i,j,k)=uu(4,i,j,k)+sf(4,i,j,k)*dt*(2./3.)
        uu(5,i,j,k)=uu(5,i,j,k)+sf(5,i,j,k)*dt*(2./3.)
     
      enddo
    enddo
  enddo
!
  return
end subroutine rk3adsft
!
