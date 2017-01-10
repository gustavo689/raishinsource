!***********************************************************************
!     for Calculation of Source Terms
!***********************************************************************
!
!--------------------------------------------------------------------
subroutine caltenr(uri,tenr,gcov,gcon,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
! Calculation of energy-momentum tensor T^\mu_\nu 
! for source term (cell-center value)
!
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos  
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1
  integer :: merr
       
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: tenr(0:3,0:3,is1:ie1,js1:je1,ks1:ke1) !tenr(m,n)=T^m_n

  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: delta(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:)
  real(8), allocatable :: util(:), bcon(:), vcon(:)  

  real(8) :: de,  pr, roh, roe, gfl, bbsq
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), delta(0:3,0:3), stat=merr)
  allocate(util(1:3), bcon(1:3), vcon(1:3), stat=merr)
  allocate(ucov(0:3), ucon(0:3), bbcov(0:3), bbcon(0:3), stat=merr)
!      
!=====================================================================

  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!
!- copy of primitive variables -!
        de=uri(1,i,j,k)
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
!        vcon(1)=uri(2,i,j,k)
!        vcon(2)=uri(3,i,j,k)
!        vcon(3)=uri(4,i,j,k)
        pr=uri(5,i,j,k)
        bcon(1)=uri(7,i,j,k)
        bcon(2)=uri(8,i,j,k)
        bcon(3)=uri(9,i,j,k)
!- calculation of relativistic enthalpy -!        
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          roh=de+(gam/(gam-1.0))*pr
        elseif(ieos .eq. 1) then
          roh=(5./2.)*pr+sqrt((9./4.)*pr**2+de**2)
        elseif(ieos .eq. 2) then
          roe=(3./2.)*(pr+((3.*pr**2)/(2.0*de+sqrt(2.*pr**2+4.*de**2)) ))
          roh=de+roe+pr
        endif
 !- copy of metric terms -!
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
            gcon1(m,n)=gcon(m,n,i,j,k)
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
!- calculation of energy-momentum tensor -!
!
        do n=0,3
          do m=0,3
            tenr(m,n,i,j,k)=(roh+bbsq)*ucon(m)*ucov(n) &
                           +(pr+0.5*bbsq)*delta(m,n) &
                           -bbcon(m)*bbcov(n)
          enddo
        enddo      
!
      enddo
    enddo
  enddo
!
!
  deallocate( gcov1, gcon1, delta, util, bcon, vcon, ucov, ucon, &
              bbcov, bbcon, stat=merr)
!      
  return
end subroutine caltenr

!---------------------------------------------------------------------
subroutine calsf(sf,tenr,christ,detg,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv, gam, c0
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1  
!
  real(8) :: sf(2:5,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: tenr(0:3,0:3,is1:ie1,js1:je1,ks1:ke1) ! tenr(m,n)=T^m_n
  real(8) :: christ(0:3,0:3,0:3,is1:ie1,js1:je1,ks1:ke1) 
             ! christ(l,m,n)=lambda^l_m n
  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)
  real(8) :: tmp1a, tmp1b, tmp1c, tmp1d, tmp2a, tmp2b, tmp2c, tmp2d
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1

        tmp1a=0.d0
        tmp1b=0.d0
        tmp1c=0.d0
        tmp1d=0.d0

        do n=0,3
          do m=0,3
!
            tmp1a=tmp1a+tenr(m,n,i,j,k)*christ(n,1,m,i,j,k)
            tmp1b=tmp1b+tenr(m,n,i,j,k)*christ(n,2,m,i,j,k)
            tmp1c=tmp1c+tenr(m,n,i,j,k)*christ(n,3,m,i,j,k)
            tmp1d=tmp1d+tenr(m,n,i,j,k)*christ(n,0,m,i,j,k)
!
          enddo
        enddo
!
        sf(2,i,j,k)=detg(i,j,k)*tmp1a
        sf(3,i,j,k)=detg(i,j,k)*tmp1b
        sf(4,i,j,k)=detg(i,j,k)*tmp1c
        sf(5,i,j,k)=detg(i,j,k)*tmp1d
!        
      enddo
    enddo
  enddo
!
  return
end subroutine calsf
