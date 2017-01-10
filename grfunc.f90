!--------------------------------------------------------------------
subroutine raise1(ucov,ucon,gcona,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax
  implicit none

  integer ::  i, j, k, is1, ie1, js1, je1, ks1, ke1
  real(8) :: ucov(0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: ucon(0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcona(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
        ucon(0,i,j,k)=gcona(0,0,i,j,k)*ucov(0,i,j,k) &
                     +gcona(0,1,i,j,k)*ucov(1,i,j,k) &
                     +gcona(0,2,i,j,k)*ucov(2,i,j,k) &
                     +gcona(0,3,i,j,k)*ucov(3,i,j,k)

        ucon(1,i,j,k)=gcona(1,0,i,j,k)*ucov(0,i,j,k) &
                     +gcona(1,1,i,j,k)*ucov(1,i,j,k) &
                     +gcona(1,2,i,j,k)*ucov(2,i,j,k) &
                     +gcona(1,3,i,j,k)*ucov(3,i,j,k)

        ucon(2,i,j,k)=gcona(2,0,i,j,k)*ucov(0,i,j,k) &
                     +gcona(2,1,i,j,k)*ucov(1,i,j,k) &
                     +gcona(2,2,i,j,k)*ucov(2,i,j,k) &
                     +gcona(2,3,i,j,k)*ucov(3,i,j,k)

        ucon(3,i,j,k)=gcona(3,0,i,j,k)*ucov(0,i,j,k) &
                     +gcona(3,1,i,j,k)*ucov(1,i,j,k) &
                     +gcona(3,2,i,j,k)*ucov(2,i,j,k) &
                     +gcona(3,3,i,j,k)*ucov(3,i,j,k)
      enddo
    enddo
  enddo
!
  return
end subroutine raise1
!  
!--------------------------------------------------------------------
subroutine lower1(ucon,ucov,gcova,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax
  implicit none

  integer ::  i, j, k, is1, ie1, js1, je1, ks1, ke1
  real(8) :: ucon(0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: ucov(0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcova(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
        ucov(0,i,j,k)=gcova(0,0,i,j,k)*ucon(0,i,j,k) &
                     +gcova(0,1,i,j,k)*ucon(1,i,j,k) &
                     +gcova(0,2,i,j,k)*ucon(2,i,j,k) &
                     +gcova(0,3,i,j,k)*ucon(3,i,j,k)

        ucov(1,i,j,k)=gcova(1,0,i,j,k)*ucon(0,i,j,k) &
                     +gcova(1,1,i,j,k)*ucon(1,i,j,k) &
                     +gcova(1,2,i,j,k)*ucon(2,i,j,k) &
                     +gcova(1,3,i,j,k)*ucon(3,i,j,k)

        ucov(2,i,j,k)=gcova(2,0,i,j,k)*ucon(0,i,j,k) &
                     +gcova(2,1,i,j,k)*ucon(1,i,j,k) &
                     +gcova(2,2,i,j,k)*ucon(2,i,j,k) &
                     +gcova(2,3,i,j,k)*ucon(3,i,j,k)

        ucov(3,i,j,k)=gcova(3,0,i,j,k)*ucon(0,i,j,k) &
                     +gcova(3,1,i,j,k)*ucon(1,i,j,k) &
                     +gcova(3,2,i,j,k)*ucon(2,i,j,k) &
                     +gcova(3,3,i,j,k)*ucon(3,i,j,k)
      enddo
    enddo
  enddo
!
  return
end subroutine lower1  

