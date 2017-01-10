!***********************************************************************
!      Metric for GRMHD
!***********************************************************************
!--------------------------------------------------------------------
subroutine coord1(gcov,gcovi,gcovj,gcovk,gcon,gconi,gconj,gconk,&
                  detg,detgi,detgj,detgk,x1,x2,x3,x1a,x2a,x3a, &
                  x1b,x2b,x3b,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, c0, metric, akm, tilang
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1
  integer :: merr

  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), & !- covariant metric -!
             gcovi(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gcovj(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gcovk(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), & !- contravariant metric -!
             gconi(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gconj(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gconk(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8), allocatable :: gcov1(:,:), gcov1i(:,:), gcov1k(:,:)  
  real(8), allocatable :: gcon1(:,:), gcon1i(:,:), gcon1k(:,:)
  real(8), allocatable :: gcon2(:,:)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1), & !- determinant g -!
             detgi(is1:ie1,js1:je1,ks1:ke1), &
             detgj(is1:ie1,js1:je1,ks1:ke1), &
             detgk(is1:ie1,js1:je1,ks1:ke1)
 
  real(8) :: x1(imax), x2(jmax), x3(kmax), & !- grid position -!
             x1a(imax), x2a(jmax), x3a(kmax), &
             x1b(imax), x2b(jmax), x3b(kmax)

  real(8) :: x1aa, x1ai, x2aa, x3aa, x3ak, x3t
  real(8) :: tmp1, tmp1a, rbh1
!
  allocate(gcov1(0:3,0:3), gcov1i(0:3,0:3), gcov1k(0:3,0:3), &
           gcon1(0:3,0:3), gcon1i(0:3,0:3), gcon1k(0:3,0:3), &
           gcon2(0:3,0:3), stat=merr )
!
  rbh1=1.+sqrt(1.-akm*akm)
  write(*,*) "rbh=", rbh1
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
          x1aa=x1(i)
          x1ai=x1a(i)
          x2aa=x2(j)
          x3aa=x3(k)
          x3ak=x3a(k)
          x3t=x3b(k)
!
        if(metric .eq. 1) then !- Minkowski cartesian spacetime -!
          call mincov1(gcov1,x1aa,x2aa,x3aa)
          call mincov1(gcov1i,x1ai,x2aa,x3aa)
          call mincov1(gcov1k,x1aa,x2aa,x3ak)
          call mincon1(gcon1,x1aa,x2aa,x3aa)
          call mincon1(gcon1i,x1ai,x2aa,x3aa)
          call mincon1(gcon1k,x1aa,x2aa,x3ak)
        elseif(metric .eq. 3) then !- Minkowski spherical spacetime -!
          call sphcov1(gcov1,x1aa,x2aa,x3aa)
          call sphcov1(gcov1i,x1ai,x2aa,x3aa)
          call sphcov1(gcov1k,x1aa,x2aa,x3ak)
          call sphcon1(gcon1,x1aa,x2aa,x3aa)
          call sphcon1(gcon1i,x1ai,x2aa,x3aa)
          call sphcon1(gcon1k,x1aa,x2aa,x3ak)
        elseif(metric .eq. 103) then !- BL Schwarzschild spacetime -!
          call schcov1(gcov1,x1aa,x2aa,x3aa)
          call schcov1(gcov1i,x1ai,x2aa,x3aa)
          call schcov1(gcov1k,x1aa,x2aa,x3ak)
          call schcon1(gcon1,x1aa,x2aa,x3aa)
          call schcon1(gcon1i,x1ai,x2aa,x3aa)
          call schcon1(gcon1k,x1aa,x2aa,x3ak)             
        elseif(metric .eq. 203) then !- BL Kerr spacetime -!
          call kercov1(gcov1,x1aa,x2aa,x3aa)
          call kercov1(gcov1i,x1ai,x2aa,x3aa)
          call kercov1(gcov1k,x1aa,x2aa,x3ak)
          call kercon1(gcon1,x1aa,x2aa,x3aa)
          call kercon1(gcon1i,x1ai,x2aa,x3aa)
          call kercon1(gcon1k,x1aa,x2aa,x3ak)
          call invert_matrix2(gcov1,gcon2,4)
        elseif(metric .eq. 303) then !- Kerr-Schild spacetime -!
          call kshcov1(gcov1,x1aa,x2aa,x3aa)
          call kshcov1(gcov1i,x1ai,x2aa,x3aa)
          call kshcov1(gcov1k,x1aa,x2aa,x3ak)
          call kshcon1(gcon1,x1aa,x2aa,x3aa)
          call kshcon1(gcon1i,x1ai,x2aa,x3aa)
          call kshcon1(gcon1k,x1aa,x2aa,x3ak)
          call invert_matrix2(gcov1,gcon2,4)
        elseif(metric .eq. 403) then !- Modified Kerr-Schild spacetime -!
          call mkshcov1(gcov1,x1aa,x2aa,x3aa,x3t)
          call mkshcov1(gcov1i,x1ai,x2aa,x3aa,x3t)
          call mkshcov1(gcov1k,x1aa,x2aa,x3ak,x3t)
          call invert_matrix2(gcov1,gcon1,4)
          call invert_matrix2(gcov1i,gcon1i,4)
          call invert_matrix2(gcov1k,gcon1k,4)
        elseif(metric .eq. 503) then !- tilted Kerr-Schild spacetime -!
          call tkshcov1(gcov1,x1aa,x2aa,x3aa)
          call tkshcov1(gcov1i,x1ai,x2aa,x3aa)
          call tkshcov1(gcov1k,x1aa,x2aa,x3ak)
          call invert_matrix2(gcov1,gcon1,4)
          call invert_matrix2(gcov1i,gcon1i,4)
          call invert_matrix2(gcov1k,gcon1k,4)          
          call kshcon1(gcon2,x1aa,x2aa,x3aa)
       endif
!
        do n=0,3
          do m=0,3
            gcov(m,n,i,j,k)=gcov1(m,n)
            gcovi(m,n,i,j,k)=gcov1i(m,n)
            gcovj(m,n,i,j,k)=gcov1(m,n)
            gcovk(m,n,i,j,k)=gcov1k(m,n)
!
            gcon(m,n,i,j,k)=gcon1(m,n)
            gconi(m,n,i,j,k)=gcon1i(m,n)
            gconj(m,n,i,j,k)=gcon1(m,n)
            gconk(m,n,i,j,k)=gcon1k(m,n)
!
!            if(i .eq. 64 .and. j .eq. 4 .and. k .eq. 25) then
!              write(*,*) 'm,n,gcon1(m,n), gcon2=', m, n, gcon1(m,n), gcon2(m,n)
!
!         endif
           
          enddo
        enddo
!      
      enddo
    enddo
  enddo
!
!- Calculation of determinant g  -!
  call caldetg(detg,gcov,is1,ie1,js1,je1,ks1,ke1)
  call caldetg(detgi,gcovi,is1,ie1,js1,je1,ks1,ke1)
  call caldetg(detgj,gcovj,is1,ie1,js1,je1,ks1,ke1)
  call caldetg(detgk,gcovk,is1,ie1,js1,je1,ks1,ke1)
!
  deallocate(gcov1, gcov1i, gcov1k, gcon1, gcon1i, gcon1k, gcon2, &
             stat=merr)
!
  return
end subroutine coord1
!
!--------------------------------------------------------------------
subroutine mincov1(gcov1,x1aa,x2aa,x3aa)
!--------------------------------------------------------------------
! Calculation of covariant metric for Minkowski spacetime
!
  implicit none

  integer :: m, n

  real(8) ::  gcov1(0:3,0:3) !- covariant metric -!
  real(8) :: x1aa, x2aa, x3aa !- grid positions -!
!
  do n=0,3
    do m=0,3
      gcov1(m,n)=0.d0
    enddo
  enddo

  gcov1(0,0)=-1.0 !- g_tt -!
  gcov1(1,1)=1.0  !- g_xx -!
  gcov1(2,2)=1.0  !- g_yy -!
  gcov1(3,3)=1.0  !- g_zz -!   
!
  return
end subroutine mincov1
!
!--------------------------------------------------------------------
subroutine mincon1(gcon1,x1aa,x2aa,x3aa)
!--------------------------------------------------------------------
! Calculation of contravariant metric for Minkowski spacetime
!
  implicit none

  integer :: m, n

  real(8) ::  gcon1(0:3,0:3) !- contravariant metric -!
  real(8) :: x1aa, x2aa, x3aa !- positions -!
!
  do n=0,3
    do m=0,3
      gcon1(m,n)=0.d0
    enddo
  enddo

  gcon1(0,0)=-1.0 !- g^tt -!
  gcon1(1,1)=1.0  !- g^xx -!
  gcon1(2,2)=1.0  !- g^yy -!
  gcon1(3,3)=1.0  !- g^zz -!             
!
  return
end subroutine mincon1
!
!--------------------------------------------------------------------
subroutine sphcov1(gcov1,x1aa,x2aa,x3aa)
!--------------------------------------------------------------------
! Calculation of covariant metric for Spherical coordinates
!
  use pram, only : pi
  implicit none

  integer :: m, n

  real(8) ::  gcov1(0:3,0:3) !- covariant metric -!
  real(8) :: x1aa, x2aa, x3aa !- positions -!
!
  real(8) :: rr, theta, del, sig, aa
!
  rr=x1aa
  if(x3aa .gt. pi) then
    theta=2.*pi-x3aa
  else
    theta=x3aa
  endif
!
  del=rr**2
  sig=rr**2
  aa=rr**4 
!
  do n=0,3
    do m=0,3
      gcov1(m,n)=0.d0
    enddo
  enddo

  gcov1(0,0)=-1.0      !- g_tt -!
  gcov1(1,1)=1.0               !- g_rr -!
  gcov1(2,2)=(aa/sig)*sin(theta)**2 !- g_phi phi -!
  gcov1(3,3)=sig                   !- g_th th -!
!
  return
end subroutine sphcov1
!
!--------------------------------------------------------------------
subroutine sphcon1(gcon1,x1aa,x2aa,x3aa)
!--------------------------------------------------------------------
! Calculation of contravariant metric for Spherical coordinates
!
  use pram, only : pi
  implicit none

  integer :: m, n

  real(8) ::  gcon1(0:3,0:3) !- contravariant metric -!
  real(8) :: x1aa, x2aa, x3aa !- positions -!
!
  real(8) :: rr, theta, del, sig, aa
!
  rr=x1aa
  if(x3aa .gt. pi) then
    theta=2.*pi-x3aa
  else
    theta=x3aa
  endif
!
  del=rr**2
  sig=rr**2
  aa=rr**4 
!
  do n=0,3
    do m=0,3
      gcon1(m,n)=0.d0
    enddo
  enddo

  gcon1(0,0)=-1.0                    !- g^tt -!
  gcon1(1,1)=1.0                          !- g^rr -!
  gcon1(2,2)=1.0/(del*sin(theta)**2) !- g^phi phi -!
  gcon1(3,3)=1./sig                           !- g^th th -!            
!
  return
end subroutine sphcon1
!
!--------------------------------------------------------------------
subroutine schcov1(gcov1,x1aa,x2aa,x3aa)
!--------------------------------------------------------------------
! Calculation of covariant metric for Schwarzschild spacetime
!
  use pram, only : pi, rbh
  implicit none

  integer :: m, n

  real(8) ::  gcov1(0:3,0:3) !- convariant metric -!
  real(8) :: x1aa, x2aa, x3aa !- positions -!
!
  real(8) :: rr, theta, del, sig, aa
!
  rr=x1aa
  if(x3aa .gt. pi) then
    theta=2.*pi-x3aa
  else
    theta=x3aa
  endif
!
  del=(rr**2)-2.*rbh*rr
  sig=rr**2
  aa=rr**4 
!
  if(del .le. 0.d0) then
    write(*,*) "del < 0"
  endif
!
  do n=0,3
    do m=0,3
      gcov1(m,n)=0.d0
    enddo
  enddo

  gcov1(0,0)=-(sig-2.*rbh*rr)/sig      !- g_tt -!
  gcov1(1,1)=sig/del               !- g_rr -!
  gcov1(2,2)=(aa/sig)*sin(theta)**2 !- g_phi phi -!
  gcov1(3,3)=sig                   !- g_th th -!
!
  return
end subroutine schcov1
!
!--------------------------------------------------------------------
subroutine schcon1(gcon1,x1aa,x2aa,x3aa)
!--------------------------------------------------------------------
! Calculation of contravariant metric for Schwarzschild spacetime
!
  use pram, only : pi, rbh
  implicit none

  integer :: m, n

  real(8) ::  gcon1(0:3,0:3) !- contravariant metric -!
  real(8) :: x1aa, x2aa, x3aa !- positions -!
!
  real(8) :: rr, theta, del, sig, aa
!
  rr=x1aa
  if(x3aa .gt. pi) then
    theta=2.*pi-x3aa
  else
    theta=x3aa
  endif
!
  del=(rr**2)-2.*rbh*rr
  sig=rr**2
  aa=rr**4 
!
  if(del .le. 0.d0) then
    write(*,*) "rr < 2.*rbh"
  endif
!
  do n=0,3
    do m=0,3
      gcon1(m,n)=0.d0
    enddo
  enddo

  gcon1(0,0)=-aa/(sig*del)                    !- g^tt -!
  gcon1(1,1)=del/sig                          !- g^rr -!
  gcon1(2,2)=(sig-2.*rbh*rr)/(del*sig*sin(theta)**2) !- g^phi phi -!
  gcon1(3,3)=1./sig                           !- g^th th -!             
!
  return
end subroutine schcon1
!
!--------------------------------------------------------------------
subroutine kercov1(gcov1,x1aa,x2aa,x3aa)
!--------------------------------------------------------------------
! Calculation of covariant metric for Kerr spacetime
!
  use pram, only : akm, pi
  implicit none

  integer :: m, n

  real(8) ::  gcov1(0:3,0:3) !- convariant metric -!
  real(8) :: x1aa, x2aa, x3aa !- positions -!
!
  real(8) :: rr, theta, del, sig, aa
  real(8) :: sth, s2, cth, a2, r2
  real(8), parameter :: small=1.d-12
!
  rr=x1aa
  if(x3aa .gt. pi) then
    theta=2.*pi-x3aa
  else
    theta=x3aa
  endif
!
  sth=sin(theta)
  if(abs(sth) .lt. small) then
    if(sth .gt. 0.d0) then
      sth=small
    elseif(sth .le. 0.d0) then
      sth=-small
    endif   
  endif
 
  s2=sth*sth
  
  cth=cos(theta)
  if(abs(cth) .lt. small) then
    if(cth .gt. 0.d0) then
      cth=small
    elseif(cth .le. 0.d0) then
      cth=-small
    endif   
  endif
 
  a2=akm*akm
  r2=rr*rr
  del=1.-(2./rr)+(a2*(1./r2))
  sig=1.+a2*cth*cth*(1./r2)

!  del=(rr**2)-2.*rr+(akm**2)
!  sig=rr**2+(akm**2)*(cos(theta)**2)
!  aa=(((rr**2)+(akm**2))**2)-(akm**2)*del*(sin(theta)**2) 
!
  if(del .le. 0.d0) then
    write(*,*) "del < 0"
  endif
!
  do n=0,3
    do m=0,3
      gcov1(m,n)=0.d0
    enddo
  enddo

!  gcov1(0,0)=-(sig-2.*rr)/sig      !- g_tt -!
!  gcov1(1,1)=sig/del               !- g_rr -!
!  gcov1(2,2)=(aa/sig)*sin(theta)**2 !- g_phi phi -!
!  gcov1(0,2)=-(2.*akm*rr*(sin(theta)**2))/sig !- g_t phi -!
!  gcov1(2,0)=gcov1(0,2)              !- g_phi t -!
!  gcov1(3,3)=sig                   !- g_th th -!
!

  gcov1(0,0)=-(1.-2.*(1./(rr*sig)))      !- g_tt -!
  gcov1(0,2)=-2.*akm*s2*(1./(rr*sig)) !- g_t phi -!
!
  gcov1(1,1)= sig/del              !- g_rr -!
!
  gcov1(2,0)=gcov1(0,2)              !- g_phi t -!
  gcov1(2,2)= r2*s2*(1.+(a2*(1./r2))+2.*a2*s2*(1./(r2*rr*sig)))  !- g_phi phi -!
!
  gcov1(3,3)=r2*sig                   !- g_th th -!
!
  return
end subroutine kercov1
!
!--------------------------------------------------------------------
subroutine kercon1(gcon1,x1aa,x2aa,x3aa)
!--------------------------------------------------------------------
! Calculation of contravariant metric for Kerr spacetime
! 
  use pram, only : akm, pi
  implicit none

  integer :: m, n

  real(8) ::  gcon1(0:3,0:3) !- contravariant metric -!
  real(8) :: x1aa, x2aa, x3aa !- positions -!
!
  real(8) :: rr, theta, del, sig, aa
  real(8) :: sth, cth, a2, r2, r3
  real(8), parameter :: small=1.d-12
!
  rr=x1aa
  if(x3aa .gt. pi) then
    theta=2.*pi-x3aa
  else
    theta=x3aa
  endif
!
  sth=sin(theta)
  if(abs(sth) .lt. small) then
    if(sth .gt. 0.d0) then
      sth=small
    elseif(sth .le. 0.d0) then
      sth=-small
    endif   
  endif
  
  cth=cos(theta)
  if(abs(cth) .lt. small) then
    if(cth .gt. 0.d0) then
      cth=small
    elseif(cth .le. 0.d0) then
      cth=-small
    endif   
  endif

  a2=akm*akm
  r2=rr*rr
  r3=r2*rr 
  del=1.-(2./rr)+(a2*(1./r2))
  sig=1.+a2*cth*cth*(1./r2)

!  del=(rr**2)-2.*rr+(akm**2)
!  sig=rr**2+(akm**2)*(cos(theta)**2)
!  aa=(((rr**2)+(akm**2))**2)-(akm**2)*del*(sin(theta)**2) 
!
  if(del .le. 0.d0) then
    write(*,*) "del < 0"
  endif
!
  do n=0,3
    do m=0,3
      gcon1(m,n)=0.d0
    enddo
  enddo

!  gcon1(0,0)=-aa/(sig*del)               !- g^tt -!
!  gcon1(1,1)=del/sig                          !- g^rr -!
!  gcon1(2,2)=(sig-2.*rr)/(del*sig*sin(theta)**2) !- g^phi phi -!
!  gcon1(0,2)=-(2.*akm*rr)/(sig*del)           !- g^t phi -!
!  gcon1(2,0)=gcon1(0,2)                  !- g^phi t -!
!  gcon1(3,3)=1./sig 
!
  gcon1(0,0)=-1.-2.*(1.+a2*(1./r2))/(rr*del*sig)               !- g^tt -!
  gcon1(0,2)=-2.*akm*(1./(r3*del*sig))           !- g^t phi -!
!
  gcon1(1,1)=del/sig                          !- g^rr -!
!
  gcon1(2,0)=gcon1(0,2)                  !- g^phi t -!
  gcon1(2,2)=(1.-2./(rr*sig))/(r2*sth*sth*del) !- g^phi phi -!
!
  gcon1(3,3)=1./(r2*sig) 

  return
end subroutine kercon1
!
!--------------------------------------------------------------------
subroutine kshcov1(gcov1,x1aa,x2aa,x3aa)
!--------------------------------------------------------------------
! Calculation of covariant metric for Kerr-Schild spacetime
!
  use pram, only : akm, pi
  implicit none

  integer :: m, n

  real(8) ::  gcov1(0:3,0:3) !- convariant metric -!
  real(8) :: x1aa, x2aa, x3aa !- positions -!
!
  real(8) :: rr, theta, del, sig, aa
  real(8) :: sth, s2, cth, a2, r2
  real(8), parameter :: small=1.d-12
!
  rr=x1aa
  if(x3aa .gt. pi) then
    theta=2.*pi-x3aa
  else
    theta=x3aa
  endif
!
  sth=sin(theta)
  if(abs(sth) .lt. small) then
    if(sth .gt. 0.d0) then
      sth=small
    elseif(sth .le. 0.d0) then
      sth=-small
    endif   
  endif
  s2=sth*sth

  cth=cos(theta)
  if(abs(cth) .lt. small) then
    if(cth .gt. 0.d0) then
      cth=small
    elseif(cth .le. 0.d0) then
      cth=-small
    endif   
  endif
 
  a2=akm*akm
  r2=rr*rr
  del=r2-(2.*rr)+a2
  sig=r2+a2*cth*cth
  aa=((r2+a2)**2)-a2*del*s2

!  del=(rr**2)-2.*rr+(akm**2)
!  sig=rr**2+(akm**2)*(cos(theta)**2)
!  aa=(((rr**2)+(akm**2))**2)-(akm**2)*del*(sin(theta)**2) 
!
!  if(del .le. 0.d0) then
!    write(*,*) "del < 0"
!  endif
!
  do n=0,3
    do m=0,3
      gcov1(m,n)=0.d0
    enddo
  enddo

  gcov1(0,0)=-1.+2.*rr/sig      !- g_tt -!
  gcov1(0,1)=2.*rr/sig           !- g_tr -!
  gcov1(0,2)=-2.*akm*rr*s2*(1./sig) !- g_t phi -!
!
  gcov1(1,0)=gcov1(0,1)        !- g_rt -!
  gcov1(1,1)=1.+2.*rr/sig      !- g_rr -!
  gcov1(1,2)=-1.*akm*s2*(1.+2.*rr/sig) !- g_r phi -!
! 
  gcov1(2,0)=gcov1(0,2)        !- g_phi t -!
  gcov1(2,1)=gcov1(1,2)        !- g_phi r -!
  gcov1(2,2)=aa*s2*(1./sig)         !- g_phi phi -!
!  gcov1(2,2)=s2*(sig+a2*s2*(1.+2.*rr/sig))
!
  gcov1(3,3)=sig                !- g_th th -!
!
  return
end subroutine kshcov1
!
!--------------------------------------------------------------------
subroutine kshcon1(gcon1,x1aa,x2aa,x3aa)
!--------------------------------------------------------------------
! Calculation of contravariant metric for Kerr-Schild spacetime
!
  use pram, only : akm, pi
  implicit none

  integer :: m, n

  real(8) ::  gcon1(0:3,0:3) !- contravariant metric -!
  real(8) :: x1aa, x2aa, x3aa !- positions -!
!
  real(8) :: rr, theta, del, sig
  real(8) :: sth, s2, cth, a2, r2
  real(8), parameter :: small=1.d-12
!
  rr=x1aa
  if(x3aa .gt. pi) then
    theta=2.*pi-x3aa
  else
    theta=x3aa
  endif
!
  sth=sin(theta)
  if(abs(sth) .lt. small) then
    if(sth .gt. 0.d0) then
      sth=small
    elseif(sth .le. 0.d0) then
      sth=-small
    endif   
  endif
  s2=sth*sth

  cth=cos(theta)
  if(abs(cth) .lt. small) then
    if(cth .gt. 0.d0) then
      cth=small
    elseif(cth .le. 0.d0) then
      cth=-small
    endif   
  endif
 
  a2=akm*akm
  r2=rr*rr
  del=r2-(2.*rr)+a2
  sig=r2+a2*cth*cth

!  del=(rr**2)-2.*rr+(akm**2)
!  sig=rr**2+(akm**2)*(cos(theta)**2)
!  aa=(((rr**2)+(akm**2))**2)-(akm**2)*del*(sin(theta)**2) 
!
!  if(del .le. 0.d0) then
!    write(*,*) "del < 0"
!  endif
!
  do n=0,3
    do m=0,3
      gcon1(m,n)=0.d0
    enddo
  enddo

  gcon1(0,0)=-1.-2.*rr/sig      !- g^tt -!
  gcon1(0,1)=2.*rr/sig           !- g^tr -!
!
  gcon1(1,0)=gcon1(0,1)        !- g^rt -!
  gcon1(1,1)=del/sig           !- g^rr -!
  gcon1(1,2)=akm*(1./sig)            !- g^r phi -!
! 
  gcon1(2,1)=gcon1(1,2)        !- g^phi r -!
  gcon1(2,2)=1./(sig*s2)         !- g^phi phi -!
!
  gcon1(3,3)=1./sig                !- g_th th -!
!
  return
end subroutine kshcon1
!
!--------------------------------------------------------------------
subroutine mkshcov1(gcov1,x1aa,x2aa,x3aa,x3t)
!--------------------------------------------------------------------
! Calculation of covariant metric for Modified Kerr-Schild spacetime
!
  use pram, only : akm, pi, ix1, ix2, ix3, R0, hslope, &
                   zmax, zmin, imax, jmax, kmax
  implicit none

  integer :: m, n

  real(8) ::  gcov1(0:3,0:3) !- convariant metric -!
  real(8) :: x1aa, x2aa, x3aa !- positions -!
!
  real(8) :: rr, theta, del, sig, aa
  real(8) :: sth, s2, cth, a2, r2, &
             tfac, rfac, hfac, pfac, x3t
  real(8), parameter :: small=1.d-12  
!
  rr=x1aa
  if(x3aa .gt. pi) then
    theta=2.*pi-x3aa
  else
    theta=x3aa
  endif
!
  sth=sin(theta)
  if(abs(sth) .lt. small) then
    if(sth .gt. 0.d0) then
      sth=small
    elseif(sth .le. 0.d0) then
      sth=-small
    endif   
  endif
  s2=sth*sth

  cth=cos(theta)
  if(abs(cth) .lt. small) then
    if(cth .gt. 0.d0) then
      cth=small
    elseif(cth .le. 0.d0) then
      cth=-small
    endif   
  endif

  a2=akm*akm
  r2=rr*rr
  del=r2-(2.*rr)+a2
  sig=r2+a2*cth*cth
  aa=((r2+a2)**2)-a2*del*s2

!- Modification parameter -!
!
  tfac=1.d0

  if(ix1 .eq. 1) then
    rfac=1.d0
  elseif(ix1 .eq. 2) then
    rfac=rr
  endif

  pfac=1.d0

  if(ix3 .eq. 1 ) then
    hfac=1.d0
  elseif(ix3 .eq. 2) then
    hfac=1.+(1.-hslope)*cos(2.*x3t)
!    write(*,*) 'hfac=', hfac 
  endif

!  del=(rr**2)-2.*rr+(akm**2)
!  sig=rr**2+(akm**2)*(cos(theta)**2)
!  aa=(((rr**2)+(akm**2))**2)-(akm**2)*del*(sin(theta)**2) 
!
!  if(del .le. 0.d0) then
!    write(*,*) "del < 0"
!  endif
!
  do n=0,3
    do m=0,3
      gcov1(m,n)=0.d0
    enddo
  enddo

  gcov1(0,0)=(-1.+2.*rr/sig) *tfac*tfac      !- g_tt -!
  gcov1(0,1)=(2.*rr/sig) *tfac*rfac           !- g_tr -!
  gcov1(0,2)=(-2.*akm*rr*s2*(1./sig)) *tfac*pfac !- g_t phi -!
!
  gcov1(1,0)=gcov1(0,1)        !- g_rt -!
  gcov1(1,1)=(1.+2.*rr/sig) *rfac*rfac      !- g_rr -!
  gcov1(1,2)=(-1.*akm*s2*(1.+2.*rr/sig)) *rfac*pfac !- g_r phi -!
! 
  gcov1(2,0)=gcov1(0,2)        !- g_phi t -!
  gcov1(2,1)=gcov1(1,2)        !- g_phi r -!
  gcov1(2,2)=(aa*s2*(1./sig)) *pfac*pfac         !- g_phi phi -!
!  gcov1(2,2)=(s2*(sig+a2*s2*(1.+2.*rr/sig))) *pfac*pfac
!
  gcov1(3,3)=sig *hfac*hfac     !- g_th th -!
!
  return
end subroutine mkshcov1
!
!--------------------------------------------------------------------
subroutine tkshcov1(gcov1,x1aa,x2aa,x3aa)
!--------------------------------------------------------------------
! Calculation of covariant metric for tilted Kerr-Schild spacetime
!
  use pram, only : akm, pi
  implicit none

  integer :: m, n

  real(8) ::  gcov1(0:3,0:3) !- convariant metric -!
  real(8) :: x1aa, x2aa, x3aa !- positions -!
!
  real(8) :: rr, phi, theta
  real(8) :: sph1, cph1, sth1, cth1, pp1, pp2, tt1, tt2
  real(8) :: del, sig, aa, s2, a2, r2
!
!- set position in tilted grid -!
  rr=x1aa
! 
  if(x2aa .le. 0.d0) then
    phi=2.*pi-x2aa
  elseif(x2aa .gt. 2.*pi) then
    phi=x2aa-2.*pi
  else
    phi=x2aa
  endif
!
  if(x3aa .le. 0.d0) then
    theta=2.*pi-x3aa
  elseif(x3aa .gt. 2.*pi) then
    theta=x3aa-2.*pi
  else
    theta=x3aa
  endif
!
  call calrotposder(phi,theta,sph1,cph1,sth1,cth1,&
                    pp1,pp2,tt1,tt2)
!
!  if (rr .gt. 3.d0 .and. rr .le. 3.1d0 .and. &
!       theta .gt. 1.5608d0 .and. theta .le. 1.5808d0) then
!     write(*,*) "r, theta, p1, p2, t1, t2=", rr, theta, pp1, pp2, tt1, tt2
!  endif   
  
  s2=sth1*sth1
  a2=akm*akm
  r2=rr*rr
  del=r2-(2.*rr)+a2
  sig=r2+a2*cth1*cth1
  aa=(((r2+a2)**2)-a2*del*s2)*s2/sig
!
  do n=0,3
    do m=0,3
      gcov1(m,n)=0.d0
    enddo
  enddo

  gcov1(0,0)=-1.+2.*rr/sig      !- g_tt -!
  gcov1(0,1)=2.*rr/sig           !- g_tr -!
  gcov1(0,2)=-2.*akm*rr*s2*pp2*(1./sig) !- g_t phi -!
  gcov1(0,3)=-2.*akm*rr*s2*pp1*(1./sig) !- g_t th -!
!
  gcov1(1,0)=gcov1(0,1)        !- g_rt -!
  gcov1(1,1)=1.+2.*rr/sig      !- g_rr -!
  gcov1(1,2)=-1.*akm*s2*(1.+2.*rr/sig)*pp2 !- g_r phi -!
  gcov1(1,3)=-1.*akm*s2*(1.+2.*rr/sig)*pp1 !- g_r th -!
! 
  gcov1(2,0)=gcov1(0,2)        !- g_phi t -!
  gcov1(2,1)=gcov1(1,2)        !- g_phi r -!
  gcov1(2,2)=(sig*tt2*tt2+aa*pp2*pp2) !- g_phi phi -!
  gcov1(2,3)=(sig*tt1*tt2+aa*pp1*pp2)  !- g_phi th -!
!
  gcov1(3,0)=gcov1(0,3)         !- g_th t -!
  gcov1(3,1)=gcov1(1,3)         !- g_th r -!
  gcov1(3,2)=gcov1(2,3)         !- g_th phi -!
  gcov1(3,3)=(sig*tt1*tt1+aa*pp1*pp1) !- g_th th -!
!
  return
end subroutine tkshcov1
!
!--------------------------------------------------------------------
subroutine calrotposder(phi,theta,sph1,cph1,sth1,cth1,&
                        pp1,pp2,tt1,tt2)
!--------------------------------------------------------------------
!- Calculation of rotated position (sin theta & phi, cos theta & phi) -!
!- and derivatives -!
!
  use pram, only: pi, tilang
  implicit none

  real(8) :: phi, theta
  real(8) :: phip, phim, thetap, thetam, delta
!  
  real(8) :: sph, cph, sth, cth, sph1, cph1, sth1, cth1, stil, ctil
  real(8) :: sph1pp, cph1pp, sth1pp, cth1pp,sph1pm, cph1pm, sth1pm, cth1pm 
  real(8) :: sph1tp, cph1tp, sth1tp, cth1tp,sph1tm, cph1tm, sth1tm, cth1tm
!
  real(8) :: dsph1dth, dsph1dph, dcth1dth, dcth1dph
  real(8) :: pp1, pp2, tt1, tt2
  !
!- Parameter
  delta=1.0d-5 !- small position difference -!
!
!- set \pm positions -!
  phip=phi+delta
  phim=phi-delta
  thetap=theta+delta
  thetam=theta-delta
!
!- calculation of tilted position -!
  call calrotpos(phi,theta,sph1,cph1,sth1,cth1)
  call calrotpos(phip,theta,sph1pp,cph1pp,sth1pp,cth1pp)  
  call calrotpos(phim,theta,sph1pm,cph1pm,sth1pm,cth1pm) 
  call calrotpos(phi,thetap,sph1tp,cph1tp,sth1tp,cth1tp)  
  call calrotpos(phi,thetam,sph1tm,cph1tm,sth1tm,cth1tm)
!
  dsph1dth=(sph1tp-sph1tm)/(thetap-thetam)
  dsph1dph=(sph1pp-sph1pm)/(phip-phim)
  dcth1dth=(cth1tp-cth1tm)/(thetap-thetam)
  dcth1dph=(cth1pp-cth1pm)/(phip-phim)
!
  pp1=(1./cph1)*dsph1dth
  pp2=(1./cph1)*dsph1dph
  tt1=(-1./sth1)*dcth1dth
  tt2=(-1./sth1)*dcth1dph
!  tt1=(-1./sth1)*dcth1dth
!  tt2=(-1./sth1)*dcth1dph

!  if (theta .gt. 1.5608d0 .and. theta .le. 1.5808d0) then
!     write(*,*) "theta=", theta
!     write(*,*) "sth, cth, sth1, cth1=", sin(theta), cos(theta), sth1, cth1
!     write(*,*) "sph, cph, sph1, cph1=", sin(phi), cos(phi), sph1, cph1     
!     write(*,*) "p1, p2, t1, t2=", pp1, pp2, tt1, tt2
!
!  endif
  
  return
end subroutine calrotposder
!
!--------------------------------------------------------------------
subroutine calrotpos(phi,theta,sph1,cph1,sth1,cth1)
!--------------------------------------------------------------------

  use pram, only: pi, tilang
  implicit none
  
  real(8) :: phi, theta, thetan
  real(8) :: sph, cph, sth, cth, stil, ctil
  real(8) :: sph1, cph1, sth1, cth1
!
!- parameter -!
  thetan=theta+tilang
!- cal sin & cos of grid pisition -!

  sph=sin(phi)
  cph=cos(phi)
  sth=sin(theta)
  cth=cos(theta)
  stil=sin(tilang)
  ctil=cos(tilang)
!
  cth1=(ctil*cth-stil*sth*cph)
!
  if(thetan .gt. pi) then
    sth1=-sqrt(1.-cth1*cth1)
  else
    sth1=sqrt(1.-cth1*cth1)
  endif
! 
  sph1=(sth*sph)*(1./sth1)
  cph1=(ctil*sth*cph+stil*cth)*(1./sth1)
!  
  return
end subroutine calrotpos
!
!--------------------------------------------------------------------
subroutine cal3met(gcov1,gcon1,gcov3a,gcon3a)
!--------------------------------------------------------------------
!- Calculation of spatial 3-metric -!
!
  implicit none

  integer :: m, n, merr

  real(8) ::  gcov1(0:3,0:3) !- convariant 4-metric -!
  real(8) ::  gcon1(0:3,0:3) !- contravariant 4-metric -!
  real(8) ::  gcov3a(1:3,1:3) !- convariant 3-metric -!
  real(8) ::  gcon3a(1:3,1:3) !- contravariant 3-metric -!

  real(8), allocatable :: oncov(:), oncon(:)
  real(8), allocatable :: beta1(:)
  real(8) :: alpha1 
!
  allocate(oncov(0:3), oncon(0:3), beta1(1:3), stat=merr)
!
  alpha1=1./sqrt(-gcon1(0,0))
  beta1(1)=alpha1*alpha1*gcon1(0,1)
  beta1(2)=alpha1*alpha1*gcon1(0,2)
  beta1(3)=alpha1*alpha1*gcon1(0,3)   
!
  oncon(0)=1./alpha1
  oncon(1)=beta1(1)/alpha1
  oncon(2)=beta1(2)/alpha1
  oncon(3)=beta1(3)/alpha1
!
  oncov(0)=-alpha1
  oncov(1)=0.d0
  oncov(2)=0.d0
  oncov(3)=0.d0
!
  do m=1,3
    do n=1,3
      gcov3a(m,n)=gcov1(m,n)+(oncov(m)*oncov(n))
      gcon3a(m,n)=gcon1(m,n)+(oncon(m)*oncon(n))
    enddo
  enddo
!
  deallocate(oncov, oncon, beta1, stat=merr)
  return
end subroutine cal3met
!
!--------------------------------------------------------------------
subroutine caldetg(detg,gcov,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
! Calcualtion of determinant 
!
  use pram, only : imax, jmax, kmax, c0
  implicit none
!
  integer :: i, j, k, is1, ie1, js1, je1, ks1, ke1

  real(8) ::  gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1) !- covariant metric -!
  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1) !- determinant g -!
  real(8) :: aa, bb, cc, pp, qq, rr, xx, yy, zz, tmp1, tmp2, tmp3, tmp4 
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1 
!
!  Calculation determinant of 3x3 
        aa=gcov(1,1,i,j,k)
        bb=gcov(1,2,i,j,k)
        cc=gcov(1,3,i,j,k)
        pp=gcov(2,1,i,j,k)
        qq=gcov(2,2,i,j,k)
        rr=gcov(2,3,i,j,k)
        xx=gcov(3,1,i,j,k)
        yy=gcov(3,2,i,j,k)
        zz=gcov(3,3,i,j,k)
        tmp1=aa*qq*zz+bb*rr*xx+cc*pp*yy-cc*qq*xx-bb*pp*zz-aa*rr*yy
!
        aa=gcov(0,1,i,j,k)
        bb=gcov(0,2,i,j,k)
        cc=gcov(0,3,i,j,k)
        pp=gcov(2,1,i,j,k)
        qq=gcov(2,2,i,j,k)
        rr=gcov(2,3,i,j,k)
        xx=gcov(3,1,i,j,k)
        yy=gcov(3,2,i,j,k)
        zz=gcov(3,3,i,j,k)
        tmp2=aa*qq*zz+bb*rr*xx+cc*pp*yy-cc*qq*xx-bb*pp*zz-aa*rr*yy
!
        aa=gcov(0,1,i,j,k)
        bb=gcov(0,2,i,j,k)
        cc=gcov(0,3,i,j,k)
        pp=gcov(1,1,i,j,k)
        qq=gcov(1,2,i,j,k)
        rr=gcov(1,3,i,j,k)
        xx=gcov(3,1,i,j,k)
        yy=gcov(3,2,i,j,k)
        zz=gcov(3,3,i,j,k)
        tmp3=aa*qq*zz+bb*rr*xx+cc*pp*yy-cc*qq*xx-bb*pp*zz-aa*rr*yy
!
        aa=gcov(0,1,i,j,k)
        bb=gcov(0,2,i,j,k)
        cc=gcov(0,3,i,j,k)
        pp=gcov(1,1,i,j,k)
        qq=gcov(1,2,i,j,k)
        rr=gcov(1,3,i,j,k)
        xx=gcov(2,1,i,j,k)
        yy=gcov(2,2,i,j,k)
        zz=gcov(2,3,i,j,k)
        tmp4=aa*qq*zz+bb*rr*xx+cc*pp*yy-cc*qq*xx-bb*pp*zz-aa*rr*yy
!
        detg(i,j,k)=sqrt(-1.*(gcov(0,0,i,j,k)*tmp1-gcov(1,0,i,j,k)*tmp2 &
                   +gcov(2,0,i,j,k)*tmp3-gcov(3,0,i,j,k)*tmp4))
!
      enddo
    enddo
  enddo
!
  return
end subroutine caldetg
!
!--------------------------------------------------------------------
subroutine calchrist(christ,x1,x2,x3,x1b,x2b,x3b,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
! Calculation of christfel symbols
!
  use pram, only : imax, jmax, kmax, c0, metric, ix1, ix2, ix3, &
                   R0, hslope, pi, tilang 
  implicit none
!
  integer :: i, j, k, l, m, n, is1, ie1, js1, je1, ks1, ke1
  integer :: merr

  real(8) :: christ(0:3,0:3,0:3,is1:ie1,js1:je1,ks1:ke1) ! - Christfel symbol -!
!- christ(l,m,n)=\Gamma^l_{mn} -!
!
  real(8), allocatable :: dgcov(:,:,:) !derivative of covariant metric -!
!- dgcov(l,m,n)=g_{lm,n} = dg_{lm}/dn -!
!
!- covariant metric with small \pm different positions -!
  real(8), allocatable :: gcov1ia(:,:), gcov1ib(:,:), & 
           gcov1ja(:,:), gcov1jb(:,:), gcov1ka(:,:), gcov1kb(:,:), &
           gcov1(:,:), gcon1(:,:)
!
  real(8), allocatable :: tmp2(:,:,:)   
!
  real(8) :: x1(imax), x2(jmax), x3(kmax) !- positions -!
  real(8) :: x1b(imax), x2b(jmax), x3b(kmax) !- positions for Mod. KS metric-!
!
  real(8) :: delta, tmp1i, tmp1j, tmp1k, tmp2t, tmp2i, tmp2j, tmp2k, &
             x1aa, x1ia, x1ib, x2aa, x2ja, x2jb, x3aa, x3ka, x3kb, &
             x1iab, x1ibb, x3kab, x3kbb, x3t
  real(8) :: tmp3a, tmp3b, tmp3c, tmp3d
  real(8), parameter :: small=1.d-12
!
! -------------------------------------------------------
  allocate(gcov1ia(0:3,0:3), gcov1ib(0:3,0:3), gcov1ja(0:3,0:3), &
           gcov1jb(0:3,0:3), gcov1ka(0:3,0:3), gcov1kb(0:3,0:3), &
           dgcov(0:3,0:3,0:3), gcov1(0:3,0:3), gcon1(0:3,0:3), &
           tmp2(0:3,0:3,0:3), stat=merr)
!
! ========================================================
! Parameter
!
  delta=1.0d-5 !- small position difference -!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1 
! - Initialize -!
        do n=0,3
          do m=0,3
            do l=0,3
              christ(n,m,l,i,j,k)=0.d0
            enddo
          enddo
        enddo
!- set \pm positions -! 
        x1aa=x1(i)
        x1ia=x1(i)+delta
        x1ib=x1(i)-delta
        x2aa=x2(j)
        x2ja=x2(j)+delta
        x2jb=x2(j)-delta
        x3aa=x3(k)
        x3ka=x3(k)+delta
        x3kb=x3(k)-delta
!
        if(metric .eq. 403) then
          if(ix1 .eq. 1) then
            x1iab=x1b(i)+delta
            x1ibb=x1b(i)-delta
            x1ia=x1iab
            x1ib=x1ibb              
          elseif(ix1 .eq. 2) then
            x1iab=x1b(i)+delta ! s
            x1ibb=x1b(i)-delta
            x1ia=exp(x1iab)    ! r_KS=exp(s) 
            x1ib=exp(x1ibb)
          endif

          if(ix3 .eq. 1) then
            x3kab=x3b(k)+delta
            x3kbb=x3b(k)-delta
            x3ka=x3kab
            x3kb=x3kbb
            x3t=x3b(k)             
          elseif(ix3 .eq. 2) then 
            x3kab=x3b(k)+delta ! theta_mKS
            x3kbb=x3b(k)-delta
            x3ka=x3kab+0.5*(1.-hslope)*sin(2.*x3kab) ! theta_KS
            x3kb=x3kbb+0.5*(1.-hslope)*sin(2.*x3kbb)
            x3t=x3b(k)
          endif
        endif

!- calculation of metrics at each positions -!
        if(metric .eq. 1) then !- Minkowski cartesian spacetime -!
          call mincov1(gcov1ia,x1ia,x2aa,x3aa)
          call mincov1(gcov1ib,x1ib,x2aa,x3aa)
          call mincov1(gcov1ja,x1aa,x2ja,x3aa)
          call mincov1(gcov1jb,x1aa,x2jb,x3aa)
          call mincov1(gcov1ka,x1aa,x2aa,x3ka)
          call mincov1(gcov1kb,x1aa,x2aa,x3kb)
          call mincon1(gcon1,x1aa,x2aa,x3aa)
       elseif(metric .eq. 3) then !- Minkowski spherical spacetime -!
          call sphcov1(gcov1ia,x1ia,x2aa,x3aa)
          call sphcov1(gcov1ib,x1ib,x2aa,x3aa)
          call sphcov1(gcov1ja,x1aa,x2ja,x3aa)
          call sphcov1(gcov1jb,x1aa,x2jb,x3aa)
          call sphcov1(gcov1ka,x1aa,x2aa,x3ka)
          call sphcov1(gcov1kb,x1aa,x2aa,x3kb)
          call sphcon1(gcon1,x1aa,x2aa,x3aa)
        elseif(metric .eq. 103) then !- BL schwarzschild spacetime -!
          call schcov1(gcov1ia,x1ia,x2aa,x3aa)
          call schcov1(gcov1ib,x1ib,x2aa,x3aa)
          call schcov1(gcov1ja,x1aa,x2ja,x3aa)
          call schcov1(gcov1jb,x1aa,x2jb,x3aa)
          call schcov1(gcov1ka,x1aa,x2aa,x3ka)
          call schcov1(gcov1kb,x1aa,x2aa,x3kb)
          call schcon1(gcon1,x1aa,x2aa,x3aa)
        elseif(metric .eq. 203) then !- BL Kerr spacetime -!
          call kercov1(gcov1ia,x1ia,x2aa,x3aa)
          call kercov1(gcov1ib,x1ib,x2aa,x3aa)
          call kercov1(gcov1ja,x1aa,x2ja,x3aa)
          call kercov1(gcov1jb,x1aa,x2jb,x3aa)
          call kercov1(gcov1ka,x1aa,x2aa,x3ka)
          call kercov1(gcov1kb,x1aa,x2aa,x3kb)
          call kercon1(gcon1,x1aa,x2aa,x3aa)
        elseif(metric .eq. 303) then !- Kerr-schild spacetime -!
          call kshcov1(gcov1ia,x1ia,x2aa,x3aa)
          call kshcov1(gcov1ib,x1ib,x2aa,x3aa)
          call kshcov1(gcov1ja,x1aa,x2ja,x3aa)
          call kshcov1(gcov1jb,x1aa,x2jb,x3aa)
          call kshcov1(gcov1ka,x1aa,x2aa,x3ka)
          call kshcov1(gcov1kb,x1aa,x2aa,x3kb)
          call kshcon1(gcon1,x1aa,x2aa,x3aa)
        elseif(metric .eq. 403) then !- Modified Kerr-schild spacetime -!
          call mkshcov1(gcov1ia,x1ia,x2aa,x3aa,x3t)
          call mkshcov1(gcov1ib,x1ib,x2aa,x3aa,x3t)
          call mkshcov1(gcov1ja,x1aa,x2ja,x3aa,x3t)
          call mkshcov1(gcov1jb,x1aa,x2jb,x3aa,x3t)
          call mkshcov1(gcov1ka,x1aa,x2aa,x3ka,x3kab)
          call mkshcov1(gcov1kb,x1aa,x2aa,x3kb,x3kbb)
          call mkshcov1(gcov1,x1aa,x2aa,x3aa,x3t)
          call invert_matrix2(gcov1,gcon1,4)
        elseif(metric .eq. 503) then !- tilted Kerr-Schild spacetime -!
          call tkshcov1(gcov1ia,x1ia,x2aa,x3aa)
          call tkshcov1(gcov1ib,x1ib,x2aa,x3aa)
          call tkshcov1(gcov1ja,x1aa,x2ja,x3aa)
          call tkshcov1(gcov1jb,x1aa,x2jb,x3aa)
          call tkshcov1(gcov1ka,x1aa,x2aa,x3ka)
          call tkshcov1(gcov1kb,x1aa,x2aa,x3kb)
          call tkshcov1(gcov1,x1aa,x2aa,x3aa)
          call invert_matrix2(gcov1,gcon1,4)          
       endif
!
!- calcualtion of derivative of covariant metrics -!
        do m=0,3
          do l=0,3
            dgcov(l,m,0)=0.d0
!- i-direction -!
            tmp1i=gcov1ia(l,m)-gcov1ib(l,m)
            if(metric .eq. 403) then
              dgcov(l,m,1)=tmp1i/(x1iab-x1ibb+small) ! derivative s
            else
              dgcov(l,m,1)=tmp1i/(x1ia-x1ib+small)
            endif
!- j-direction -!
            tmp1j=gcov1ja(l,m)-gcov1jb(l,m)
            dgcov(l,m,2)=tmp1j/(x2ja-x2jb+small)
!- k-direction -!
            tmp1k=gcov1ka(l,m)-gcov1kb(l,m)
            if(metric .eq. 403) then
              dgcov(l,m,3)=tmp1k/(x3kab-x3kbb+small) ! derivative theta_mKS
            else                
              dgcov(l,m,3)=tmp1k/(x3ka-x3kb+small)
            endif          
          enddo
        enddo
!
!- calculation of christfel symbols -!
        do n=0,3
          do m=0,3
            do l=0,3
!- tmp2(l,m,n)=0.5*(g_{lm,n}+g_{ln,m}-g_{mn,l}) -!
!- dgcov(l,m,n)=g_{lm,n} -!
!
              tmp2(l,m,n)=0.5d0*(dgcov(m,l,n)+dgcov(n,l,m)-dgcov(n,m,l)) ! from HARM
!              tmp2(l,m,n)=0.5d0*(dgcov(l,m,n)+dgcov(l,n,m)-dgcov(m,n,l))  ! from my calculation
!
            enddo
          enddo
        enddo

        do n=0,3
          do m=0,3
            do l=0,3
! from HARM
              christ(l,m,n,i,j,k)=gcon1(l,0)*tmp2(0,m,n) &
                                 +gcon1(l,1)*tmp2(1,m,n) &
                                 +gcon1(l,2)*tmp2(2,m,n) &
                                 +gcon1(l,3)*tmp2(3,m,n)
!
            enddo
          enddo
        enddo
!
!       do m=0,3
!         do l=0,3
!           christ(l,m,0,i,j,k)=gcon1(0,l)*(tmp2(0,m,0)+tmp2(0,m,1) &
!                               +tmp2(0,m,2)+tmp2(0,m,3))
!           christ(l,m,1,i,j,k)=gcon1(1,l)*(tmp2(1,m,0)+tmp2(1,m,1) &
!                               +tmp2(1,m,2)+tmp2(1,m,3))
!           christ(l,m,2,i,j,k)=gcon1(2,l)*(tmp2(2,m,0)+tmp2(2,m,1) &
!                               +tmp2(2,m,2)+tmp2(2,m,3))
!           christ(l,m,3,i,j,k)=gcon1(3,l)*(tmp2(3,m,0)+tmp2(3,m,1) &
!                               +tmp2(3,m,2)+tmp2(3,m,3))
!         enddo
!       enddo


!        if(i .eq. 128 .and. j .eq. 4 .and. k .eq. 4) then
!          write(*,*) 'dgcov(1,1,1)', dgcov(1,1,1)
!          write(*,*) 'christ(1,1,1)', christ(1,1,1,i,j,k)
!        endif
!       
      enddo
    enddo
  enddo
!
  deallocate(gcov1ia, gcov1ib, gcov1ja, gcov1jb, gcov1ka, gcov1kb, &
             dgcov, gcov1, gcon1, tmp2, stat=merr)
  return
end subroutine calchrist
!
!--------------------------------------------------------------------
subroutine invert_matrix(A_org,A_inv,N)
!--------------------------------------------------------------------
! Invertion of Matrix for calculation of contravariant metric 
!
  implicit none
!  
  integer :: i, j, k, N
  real(8) :: A(N,N), b(N)
  real(8) :: A_org(N,N), A_inv(N,N)
  real(8) :: col(N), E(N,N)
  integer :: p(N)
  integer :: flag
  real(8) :: U(N,N),L(N,N)
  real(8) :: tmp, d
  character :: str*10

!-  initialize -!
  do i=1,N 
    b(i)=1.d0
    p(i)=i
  enddo
  flag=1

!- copy original matrix -!
  do j=1,N
    do i=1,N
      A(i,j)=A_org(i,j)
    enddo
  enddo 

!- LU decompose -!
  d=1.d0
  call lu_decompose(A,N,p,d,flag)
!
  if(flag .eq. 0) then
    goto 300
  endif

!- 
  call lu_solve(A,N,p,b)
  str = 'A after LU'
!
!- Calculate invert matrix -!
  do i=1,N

    do j=1,N
      col(j)=0.d0
    enddo

    
    col(i)=1.d0

    call lu_solve(A,N,p,col)
    do j=1,N
      A_inv(i,j)=col(j)
    enddo
  enddo
  str = 'A inverse'
!
  do j=1,N
    do i=1,N
      E(i,j)=0.d0
      do k=1,N
        E(i,j) = E(i,j) + A_org(k,j)*A_inv(i,k)
      enddo
    enddo
  enddo
  str = 'A-1*A = E'
!
300 continue
!
  if(flag .eq. 0) then
    write(6,*) "A is singular matrix"
  endif
!
  return
end subroutine invert_matrix
!
!--------------------------------------------------------------------
subroutine lu_decompose(A,N,p,d,flag)
!--------------------------------------------------------------------
! Calculation of LU decomposition
!
  implicit none
!
  real(8) :: Eps
  integer :: i, j, t, N, mj

  real(8) :: tmp
  integer :: tmp2
  real(8) :: A(N,N)
  integer :: p(N)
  integer :: flag
  real(8) :: d 
!
  Eps = 1d-10
      
  do t=1,N
!####  PIVOT #####
    mj=t
    do j=t+1,N
      if(abs(A(t,j)) .gt. abs(A(t,mj)) ) then
        mj=j
      endif
    enddo
!
    if(mj .ne. t) then
      do i=1,N
        tmp=A(i,mj)
        A(i,mj)=A(i,t)
        A(i,t)=tmp
      enddo
      tmp2=p(mj)
      p(mj)=p(t)
      p(t)=tmp2
      d=-d
    endif
!####  END of PIVOT #####

    if(abs(A(t,t)) .lt. Eps) then
      write(6,*) "A is singular matrix"
      flag = 0
      return
    endif
        
    do j=t+1,N
      A(t,j)=A(t,j)/A(t,t)
      do i=t+1,N
        A(i,j)=A(i,j)-A(t,j)*A(i,t)
      enddo
    enddo
         
  enddo
      
  flag=1
  return 
end subroutine lu_decompose
!
!--------------------------------------------------------------------
subroutine lu_solve(A,N,p,b)
!--------------------------------------------------------------------
! Solving LU decompose
!
  implicit none
!   
  integer ::i,j,N
  real(8) :: A(N,N), b(N)
  integer :: p(N)
  real(8) :: sum, tmp
!
  do j=1,N
    if(j .ne. p(j) .and. j .lt. p(j)) then
      tmp=b(j)
      b(j)=b(p(j))
      b(p(j))=tmp
    endif
  enddo

!- Forward -!
  do j=1,N
    sum=b(j)
    do i=1,j-1
      sum=sum-A(i,j)*b(i)
    enddo
    b(j)=sum
  enddo

!- Backward -!
  do j=N,1,-1
    sum=b(j)
    do i=N,j+1,-1
      sum=sum-A(i,j)*b(i)
    enddo
    b(j)=sum/A(j,j)
  enddo
!
  return
end subroutine lu_solve
!
!--------------------------------------------------------------------
subroutine invert_matrix2(a,axn,n)
!--------------------------------------------------------------------
!
  implicit none
!
  integer :: n, indx(1:n)
  integer :: i, j
  real(8) :: an(1:n,1:n), a(0:n-1,0:n-1), ax(1:n,1:n), axn(0:n-1,0:n-1)
  real(8) :: y(1:n,1:n), b(1:n)           
  real(8) :: d
!  
!- input matrix A copy to matrix An

  do i=1,n
    do j=1,n 
      an(i,j)=a(i-1,j-1)
    enddo
  enddo  
!
!- set up identity matrix
  do i=1,n
    do j=1,n
      if(i .eq. j) then  
        y(i,j)=1.d0
      else
        y(i,j)=0.d0
      endif   
    enddo
  enddo  
!
  call ludcmp(an,n,indx,d)
!
  do j=1,n
    do i=1,n 
      b(i)=y(i,j)
    enddo 
!   
    call lubksb(an,n,indx,b)
!
    do i=1,n
      ax(i,j)=b(i)
    enddo         
  enddo
!
  do j=1,n
    do i=1,n
      axn(i-1,j-1)=ax(i,j)
    enddo
  enddo
!        
  return
end subroutine invert_matrix2
!
!--------------------------------------------------------------------
subroutine ludcmp(a,n,indx,d)
!--------------------------------------------------------------------
!- input :: a(n,n), n
!- output :: a(n,n), indx, d  
  
  implicit none
!
  integer :: n, indx(1:n)
  integer, parameter :: nmax=500
  real(8) :: d, a(1:n,1:n)
  real(8), parameter :: tiny=1.e-20
!  
  integer :: i, imax, j, k
  real(8) :: aamax, dum, sum, vv(nmax)
!
!- parameter
  d=1.d0
!
  do i=1,n
!
    aamax=0.d0
!
    do j=1,n 
      if (abs(a(i,j)) .gt. aamax) then
        aamax=abs(a(i,j))
      endif
    enddo
!
    if (aamax .eq. 0.d0) then
      write(*,*) "sigular matrix in ludcmp"
    endif   
    vv(i)=1./aamax
!
  enddo   
!  
  do j=1,n
!
     do i=0,j-1
      sum=a(i,j)
      do k=1,i-1
        sum=sum-a(i,k)*a(k,j) 
      enddo
      a(i,j)=sum
    enddo
!
    aamax=0.d0
!
    do i=j,n
      sum=a(i,j)
      do k=1,j-1
         sum=sum-a(i,k)*a(k,j)
      enddo
      a(i,j)=sum
      dum=vv(i)*abs(sum)
      if(dum.ge.aamax) then
        imax=i
        aamax=dum
      endif  
    enddo
!
    if(j .ne. imax) then
      do k=1,n
        dum=a(imax,k)
        a(imax,k)=a(j,k)
        a(j,k)=dum
      enddo
      d=-d
      vv(imax)=vv(j)
    endif
!
    indx(j)=imax
!
    if(a(j,j) .eq. 0.d0) then
      a(j,j)=tiny  
    endif
!
    if(j .ne. n) then
      dum=1./a(j,j)
      do i=j+1,n
        a(i,j)=a(i,j)*dum
      enddo   
    endif   
!    
  enddo
!
  return
end subroutine ludcmp  
!
!--------------------------------------------------------------------
subroutine lubksb(a,n,indx,b)
!--------------------------------------------------------------------
!
  implicit none
!
  integer :: n, indx(n)
  real(8) :: a(1:n,1:n), b(1:n)
  integer :: i,ii,j,ll
  real(8) :: sum
!
  ii=0
!  
  do i=1,n
    ll=indx(i)
    sum=b(ll)
    b(ll)=b(i)
!
    if(ii .ne. 0) then
       
      do j=ii, i-1
        sum=sum-a(i,j)*b(j)
      enddo

    else if (sum .ne. 0) then
        ii=i  
    endif
!
    b(i)=sum
!
  enddo
!
  do i=n,1,-1
    sum=b(i)
    do j=i+1, n
      sum=sum-a(i,j)*b(j) 
    enddo   
!
    b(i)=sum/a(i,i)
!    
  enddo
!
  return
end subroutine lubksb
