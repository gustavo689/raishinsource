!----------------------------------------------------------------------
subroutine atmosphere1(uu,uri,gcov,gcon,detg,x1,nm1,is1,ie1,js1,je1,ks1,ke1)
!----------------------------------------------------------------------
!     floor model for density, pressure, and velocity for vacuum atmosphere  
!
  use pram, only : imax, jmax, kmax, nv, pmin, pmax, dmin, gam, kpol, ieos
  implicit none
!
  integer :: i, j, k, m, n, nm1, is1, ie1, js1, je1, ks1, ke1, merr

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: x1(imax)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: util(:), vel(:), beta1(:), bcon(:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:)
!
  real(8) :: alpha1, pr, de, demin, prmin, vrmin, dthe, dthe1, gfl, gflmax, vsq
  real(8) :: rr, rin, epsilon, rde, bsq, bbsq, roh, roe 
  integer :: isetp1
!
!-allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), stat=merr)
  allocate(util(1:3), vel(1:3), beta1(1:3), bcon(1:3),stat=merr)
  allocate(ucov(0:3), ucon(0:3), bbcov(0:3), bbcon(0:3), stat=merr)

!--------------------------------------------------------------------!
!- Parameter -!
  epsilon=1.d-2
  gflmax=10.d0
  isetp1=0
!
!
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
      do i=is1+nm1,ie1-nm1
         
        rr=x1(i)
        rin=25.d0
!
!- copy of primitive variables and conserved variables-! 
        de=uri(1,i,j,k)
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
        pr=uri(5,i,j,k)
        bcon(1)=uri(7,i,j,k)
        bcon(2)=uri(8,i,j,k)
        bcon(3)=uri(9,i,j,k)
!
!- copy of metric terms -!
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo
!
        alpha1=1./sqrt(-gcon1(0,0))
!
        beta1(1)=alpha1*alpha1*gcon1(0,1)
        beta1(2)=alpha1*alpha1*gcon1(0,2)
        beta1(3)=alpha1*alpha1*gcon1(0,3)
!
        rde=alpha1*uu(1,i,j,k)/detg(i,j,k) !- rde= de*gamma -!
!
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

!- 4-velocity=> 3-velocity -!
          vel(1)=util(1)*(1./gfl)
          vel(2)=util(2)*(1./gfl)
          vel(3)=util(3)*(1./gfl)
!
!- set atmosphere  value -!         
        demin=dmin*(rr)**(-3./2.)          
        prmin=(gam-1.)*pmin*(rr)**(-5./2.)
!          
!         demin=dmin
!         prmin=pmin 

!         demin=dmin*exp(-3.*rr/rin)
!        prmin=(gam-1.)*pmin
!        prmin=kpol*(demin**gam)

!        vrmin=vel(1)
!        vrmin=0.1d0*vel(1)   
!         vrmin=beta1(1)*(1./alpha1)
        vrmin=0.d0
        
!- set density limitter   
!        dthe=5.0*dmin*(1.d0+epsilon)
!        dthe=2.0*(1.d0+epsilon)*dmin*rr**(-4./2.)
        dthe=1.5*(1.+epsilon)*demin
        dthe1=1.0*(1.+epsilon)*prmin
!
!- cal bsq
        bsq=bcon(1)*bcon(1)+bcon(2)*bcon(2)+bcon(3)*bcon(3)
!
!- check the density & pressure limit -!
!
!        if( de .lt. dthe .and. bsq .eq. 0.d0) then        
        if( de .lt. dthe .or. pr .lt. dthe1) then   
!        if( de .lt. dthe) then
!          de=demin*rr**(-4./2.)
!          pr=prmin*rr**(-5./2.)
          de=demin
          pr=prmin 
!          vel(1)=vrmin
!          vel(2)=0.d0
!          vel(3)=0.d0
!          bcon(1)=0.d0
!          bcon(2)=0.d0
!          bcon(3)=0.d0
          
          isetp1=1
!
!- new calculation of Lorentz factor -!
          vsq=0.d0
          do m=1,3
            do n=1,3
              vsq=vsq+gcov1(m,n)*vel(m)*vel(n)     
            enddo
          enddo
!
!          if(abs(vsq) .le. 1.d-10) then
!            vsq=1.d-10
!          endif
          if(vsq .ge. 1.d0) then
            vsq=1.-1./(gflmax*gflmax)
          endif
!
          gfl=1./sqrt(1.-vsq)
!          
          util(1)=vel(1)*gfl
          util(2)=vel(2)*gfl
          util(3)=vel(3)*gfl

        endif
          
!- set of primitive variables

        uri(1,i,j,k)=de
!
        uri(2,i,j,k)=util(1)
        uri(3,i,j,k)=util(2)
        uri(4,i,j,k)=util(3)
        uri(5,i,j,k)=pr
        uri(6,i,j,k)=0.d0
        uri(7,i,j,k)=bcon(1)
        uri(8,i,j,k)=bcon(2)
        uri(9,i,j,k)=bcon(3)
!
!- new set of conservative variables -!

        if(ieos .eq. 0 .or. ieos .eq. 3) then
          roh=de+(gam/(gam-1.0))*pr
        elseif(ieos .eq. 1) then
          roh=(5./2.)*pr+sqrt((9./4.)*pr**2+de**2)
        elseif(ieos .eq. 2) then
          roe=(3./2.)*(pr+((3.*pr**2)/(2.0*de+sqrt(2.*pr**2+4.*de**2)) ))
          roh=de+roe+pr
        endif
          
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

!          uu(1,i,j,k)=de*ucon(0) !- relativistic mass -!
!          uu(2,i,j,k)=(roh+bbsq)*ucon(0)*ucov(1)-bbcon(0)*bbcov(1)
!          uu(3,i,j,k)=(roh+bbsq)*ucon(0)*ucov(2)-bbcon(0)*bbcov(2)
!          uu(4,i,j,k)=(roh+bbsq)*ucon(0)*ucov(3)-bbcon(0)*bbcov(3)
!!          uu(5,i,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq) &
!!                     -bbcon(0)*bbcov(0)-de*ucon(0)
!          uu(5,i,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq) &
!                     -bbcon(0)*bbcov(0)+de*ucon(0)
!!          uu(5,i,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq) &
!!                     -bbcon(0)*bbcov(0)
!
!          do n=1,5
!            uu(n,i,j,k)=detg(i,j,k)*uu(n,i,j,k)
!          enddo
!
      enddo
    enddo
 enddo
!
  if(isetp1 .ne. 0) then
    call caluu3a(uri,uu,gcov,gcon,detg,x1,nm1,is1,ie1,js1,je1,ks1,ke1)
  endif
!
  deallocate(gcov1, gcon1, util, vel, beta1, bcon, &
             ucov, ucon, bbcov, bbcon, stat=merr)
!
  return
end subroutine atmosphere1
!
!----------------------------------------------------------------------
subroutine atmosphere1a(uu,uri,gcov,gcon,detg,x1,is1,ie1,js1,je1,ks1,ke1)
!----------------------------------------------------------------------
!     floor model for density, pressure, and velocity for vacuum atmosphere  
!
  use pram, only : imax, jmax, kmax, nv, pmin, pmax, dmin, gam, kpol, ieos
  implicit none
!
  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1, merr

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: x1(imax)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: util(:), vel(:), beta1(:), bcon(:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:)
!
  real(8) :: alpha1, pr, de, demin, prmin, vrmin, dthe, dthe1, gfl, gflmax, vsq
  real(8) :: rr, rin, epsilon, rde, bsq, bbsq, roh, roe 
  integer :: isetp1
!
!-allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), stat=merr)
  allocate(util(1:3), vel(1:3), beta1(1:3), bcon(1:3),stat=merr)
  allocate(ucov(0:3), ucon(0:3), bbcov(0:3), bbcon(0:3), stat=merr)

!--------------------------------------------------------------------!
!- Parameter -!
  epsilon=1.d-2
  gflmax=10.d0
  isetp1=0
!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
         
        rr=x1(i)
        rin=25.d0
!
!- copy of primitive variables and conserved variables-! 
        de=uri(1,i,j,k)
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
        pr=uri(5,i,j,k)
        bcon(1)=uri(7,i,j,k)
        bcon(2)=uri(8,i,j,k)
        bcon(3)=uri(9,i,j,k)
!
!- copy of metric terms -!
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo
!
        alpha1=1./sqrt(-gcon1(0,0))
!
        beta1(1)=alpha1*alpha1*gcon1(0,1)
        beta1(2)=alpha1*alpha1*gcon1(0,2)
        beta1(3)=alpha1*alpha1*gcon1(0,3)
!
        rde=alpha1*uu(1,i,j,k)/detg(i,j,k) !- rde= de*gamma -!
!
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

!- 4-velocity=> 3-velocity -!
          vel(1)=util(1)*(1./gfl)
          vel(2)=util(2)*(1./gfl)
          vel(3)=util(3)*(1./gfl)
!
!- set atmosphere  value -!         
          demin=dmin*((rr)**(-3./2.))
          prmin=(gam-1.)*pmin*((rr)**(-5./2.))
!          
!         demin=dmin
!         prmin=pmin 

!         demin=dmin*exp(-3.*rr/rin)
!        prmin=(gam-1.)*pmin
!        prmin=kpol*(demin**gam)
!
!        vrmin=vel(1)
!        vrmin=0.1d0*vel(1)   
!        vrmin=beta1(1)*(1./alpha1)
        vrmin=0.d0
        
!- set density limitter   
!        dthe=5.0*dmin*(1.d0+epsilon)
!        dthe=2.0*(1.d0+epsilon)*dmin*rr**(-4./2.)
        dthe=1.0*(1.+epsilon)*demin
        dthe1=1.0*(1.+epsilon)*prmin
!
!- cal bsq
        bsq=bcon(1)*bcon(1)+bcon(2)*bcon(2)+bcon(3)*bcon(3)
!
!- check the density & pressure limit -!
!
!        if( de .lt. dthe .and. bsq .le. 1.d-20) then        
        if( de .lt. dthe .or. pr .lt. dthe1 .and. bsq .le. 1.d-20) then   
!        if( de .lt. dthe .or. pr .lt. dthe1) then
!          de=demin*rr**(-4./2.)
!          pr=prmin*rr**(-5./2.)
          de=demin
          pr=prmin 
          vel(1)=vrmin
          vel(2)=0.d0
          vel(3)=0.d0
!          bcon(1)=0.d0
!          bcon(2)=0.d0
!          bcon(3)=0.d0
          
          isetp1=1
!
!- new calculation of Lorentz factor -!
          vsq=0.d0
          do m=1,3
            do n=1,3
              vsq=vsq+gcov1(m,n)*vel(m)*vel(n)     
            enddo
          enddo
!
!          if(abs(vsq) .le. 1.d-10) then
!            vsq=1.d-10
!          endif
          if(vsq .ge. 1.d0) then
            vsq=1.-1./(gflmax*gflmax)
          endif
!
          gfl=1./sqrt(1.-vsq)
!          
          util(1)=vel(1)*gfl
          util(2)=vel(2)*gfl
          util(3)=vel(3)*gfl

        endif
          
!- set of primitive variables

        uri(1,i,j,k)=de
!
        uri(2,i,j,k)=util(1)
        uri(3,i,j,k)=util(2)
        uri(4,i,j,k)=util(3)
        uri(5,i,j,k)=pr
        uri(6,i,j,k)=0.d0
        uri(7,i,j,k)=bcon(1)
        uri(8,i,j,k)=bcon(2)
        uri(9,i,j,k)=bcon(3)
!
!- new set of conservative variables -!

        if(ieos .eq. 0 .or. ieos .eq. 3) then
          roh=de+(gam/(gam-1.0))*pr
        elseif(ieos .eq. 1) then
          roh=(5./2.)*pr+sqrt((9./4.)*pr**2+de**2)
        elseif(ieos .eq. 2) then
          roe=(3./2.)*(pr+((3.*pr**2)/(2.0*de+sqrt(2.*pr**2+4.*de**2)) ))
          roh=de+roe+pr
        endif
          
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

!          uu(1,i,j,k)=de*ucon(0) !- relativistic mass -!
!          uu(2,i,j,k)=(roh+bbsq)*ucon(0)*ucov(1)-bbcon(0)*bbcov(1)
!          uu(3,i,j,k)=(roh+bbsq)*ucon(0)*ucov(2)-bbcon(0)*bbcov(2)
!          uu(4,i,j,k)=(roh+bbsq)*ucon(0)*ucov(3)-bbcon(0)*bbcov(3)
!!          uu(5,i,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq) &
!!                     -bbcon(0)*bbcov(0)-de*ucon(0)
!          uu(5,i,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq) &
!                     -bbcon(0)*bbcov(0)+de*ucon(0)
!!          uu(5,i,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq) &
!!                     -bbcon(0)*bbcov(0)
!
!          do n=1,5
!            uu(n,i,j,k)=detg(i,j,k)*uu(n,i,j,k)
!          enddo
!
      enddo
    enddo
 enddo
!
!  if(isetp1 .ne. 0) then
!    call caluu1(uri,uu,gcov,gcon,detg,x1,is1,ie1,js1,je1,ks1,ke1)
!  endif
!
  deallocate(gcov1, gcon1, util, vel, beta1, bcon, &
             ucov, ucon, bbcov, bbcon, stat=merr)
!
  return
end subroutine atmosphere1a
!
!----------------------------------------------------------------------
subroutine pminmax1a(uu,uri,gcov,gcon,detg,x1,nm1,is1,ie1,js1,je1,ks1,ke1)
!----------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv, pmin, pmax, dmin
  implicit none
!
  integer :: i, j, k, nm1, is1, ie1, js1, je1, ks1, ke1

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)

  real(8) :: x1(imax)
  
  real(8) :: pr, de
  integer :: isetp
!
  isetp=0
!
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
      do i=is1+nm1,ie1-nm1

        isetp=0
         
        pr=uri(5,i,j,k)
        de=uri(1,i,j,k)
         
        if( pr.lt.pmin .or. pr.gt.pmax .or. de.lt.dmin) then
!
          if( pr.lt.pmin ) then
            pr=pmin
            isetp=1
          endif
!
          if( pr.gt.pmax ) then
            pr=pmax
            isetp=2
          endif

          if( de.lt.dmin ) then
            de=dmin
            isetp=3
          endif
!
          uri(5,i,j,k)=pr
          uri(1,i,j,k)=de
!
        endif

      enddo
    enddo
  enddo  
!     
  if(isetp .ne. 0) then
    call caluu3a(uri,uu,gcov,gcon,detg,x1,nm1,is1,ie1,js1,je1,ks1,ke1)
  endif
!
  if( isetp .eq. 1 ) then
    write(6,*) 'Correction: pmin in pminmax'
    write(6,*) 'pmin =', pmin
  elseif( isetp .eq. 2) then
    write(6,*) 'Correction: pmax in pminmax'
    write(6,*) 'pmax =', pmax
  elseif( isetp .eq. 3) then
    write(6,*) 'Correction: dmin in pminmax'
    write(6,*) 'dmin =', dmin
  endif
!
  return
end subroutine pminmax1a
!
!----------------------------------------------------------------------
subroutine pminmax2a(uu,uri,gcov,gcon,detg,x1,nm1,is1,ie1,js1,je1,ks1,ke1)
!----------------------------------------------------------------------
!     floor model for density and pressure                       
!
  use pram, only : imax, jmax, kmax, nv, pmin, pmax, dmin, &
                   gam, kpol, ieos
  implicit none
!
  integer :: i, j, k, nm1, is1, ie1, js1, je1, ks1, ke1
  integer :: m, n, merr

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:)
  real(8), allocatable :: util(:), bcon(:), bcon1(:), bcov1(:)
  
  real(8) :: x1(imax)
!
  real(8) :: pr, de, demin, uumin, prmin, rr, rin, rcen, roh
!  integer, allocatable :: isetp(:,:,:)
  integer :: isetp1
  real(8) :: gfl, bbsq, bsq, bbromax, bbprmax, prromax, bbro, bbpr, prro
  real(8) :: tmp1a, tmp1b
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), stat=merr)
  allocate(util(1:3), bcon(1:3), stat=merr)
  allocate(ucov(0:3), ucon(0:3), bbcov(0:3), bbcon(0:3), &
           bcon1(0:3), bcov1(0:3), stat=merr)
! 
!  allocate(isetp(is1:ie1,js1:je1,ks1:ke1), stat=merr)
!
!
!  do k=ks1+nm1,ke1-nm1
!    do j=js1+nm1,je1-nm1
!      do i=is1+nm1,ie1-nm1  
!        isetp(i,j,k)=0
!      enddo
!    enddo
!  enddo
  isetp1=0 
!
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
      do i=is1+nm1,ie1-nm1
         
        rr=x1(i)
        rin=1.d0
        rcen=25.d0

!- copy of primitive variables and conserved variables-! 
        de=uri(1,i,j,k)
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)        
        pr=uri(5,i,j,k)
        bcon(1)=uri(7,i,j,k)
        bcon(2)=uri(8,i,j,k)
        bcon(3)=uri(9,i,j,k)
!        
        bcon1(0)=0.d0
        bcon1(1)=bcon(1)
        bcon1(2)=bcon(2)
        bcon1(3)=bcon(3)
!       
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
!        
!- check roh*gamma^2-p < flc*B^2 
        call lower(bcon1,bcov1,gcov1)
        bsq=bcon1(1)*bcov1(1)+bcon1(2)*bcov1(2)+bcon1(3)*bcov1(3)

95      continue
        
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          roh=de+(gam/(gam-1.0))*pr
        elseif(ieos .eq. 1) then
          roh=(5./2.)*pr+sqrt((9./4.)*pr**2+de**2)
        endif
         
        tmp1a=roh*gfl*gfl-pr
        tmp1b=abs(0.05*bsq)
        if(tmp1a .le. tmp1b) then
          de=1.05*de
          pr=1.05*pr
!          write(*,*) "roh*gf^2-p < fluc*B^2 at ", i, j, k
          go to 95
        endif 
        
!- set min & max values -!
        
        demin=dmin*(rr)**(-3./2.)
        prmin=(gam-1.)*pmin*(rr)**(-5./2.)

!        demin=dmin
!        prmin=(gam-1.)*pmin

!        demin=dmin*exp(-3.*rr/rcen)
!        prmin=kpol*(demin**gam)
!        prmin=pmin
        
!
        bbromax=1.d2
        bbprmax=1.d2
        prromax=1.d9
!
        if(bbsq .lt. 1.d-20) then
          bbro=0.d0
          bbpr=0.d0
        else
          bbro=bbsq*(1./de)
          bbpr=bbsq*(1./pr)
        endif
        prro=pr*(1./de)       
!        
!        if( pr .lt. prmin .or. pr .gt. pmax .or. de .lt. demin) then
!           
        if( pr .lt. prmin .or. pr .gt. pmax .or. de .lt. demin .or. &
             bbro .gt. bbromax .or. bbpr .gt. bbprmax &
             .or. prro .gt. prromax) then
!
          if( pr.lt.prmin ) then
            uri(5,i,j,k)=prmin
!            isetp(i,j,k)=1
            isetp1=1
         endif
!
          if( pr.gt.pmax ) then
            uri(5,i,j,k)=pmax
!            isetp(i,j,k)=2
            isetp1=1
          endif

          if( de.lt.demin ) then          
            uri(1,i,j,k)=demin
!            isetp(i,j,k)=3
            isetp1=1
          endif          
!
          if( bbro .gt. bbromax ) then
            uri(1,i,j,k)=bbsq/bbromax
            isetp1=1 
          endif
!
          if( bbpr .gt. bbprmax ) then
            uri(5,i,j,k)=bbsq/bbprmax
            isetp1=1 
          endif
!
          if(prro .gt. prromax ) then
            uri(1,i,j,k)=pr/prromax 
            isetp1=1
          endif
!         
!          uri(5,i,j,k)=pr
!          uri(1,i,j,k)=de
          
!
        endif       
!
      enddo
    enddo
 enddo
!
  if(isetp1 .ne. 0) then
    call caluu3a(uri,uu,gcov,gcon,detg,x1,nm1,is1,ie1,js1,je1,ks1,ke1)
  endif
!
!  do k=ks1+nm1,ke1-nm1
!    do j=js1+nm1,je1-nm1
!      do i=is1+nm1,ie1-nm1
!        if( isetp(i,j,k) .eq. 1) then
!          write(6,*) 'Correction: pressure min'
!          write(6,*) 'pmin =', pmin, i, j, k            
!        elseif( isetp(i,j,k) .eq. 2) then
!          write(6,*) 'Correction: pressure max'
!          write(6,*) 'pmax =', pmax, i, j, k
!        elseif( isetp(i,j,k) .eq. 3) then
!          write(6,*) 'Correction: density min'
!          write(6,*) 'dmin =', dmin, i, j, k
!        endif
!      enddo
!    enddo
!  enddo  
!
!  deallocate(isetp, stat=merr)
  deallocate(gcov1, gcon1, util, bcon, ucov, ucon, bbcov, bbcon, &
             bcon1, bcov1, stat=merr)
!
  return
end subroutine pminmax2a
!
!----------------------------------------------------------------------
subroutine pminmax2b(uu,uri,gcov,gcon,detg,x1,is1,ie1,js1,je1,ks1,ke1)
!----------------------------------------------------------------------
!     floor model for density and pressure                       
!
  use pram, only : imax, jmax, kmax, nv, pmin, pmax, dmin, &
                   gam, kpol, ieos
  implicit none
!
  integer :: i, j, k, is1, ie1, js1, je1, ks1, ke1
  integer :: m, n, merr

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:), gcov3a(:,:), gcon3a(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:)
  real(8), allocatable :: util(:), bcon(:), bcon1(:), bcov1(:)
  real(8), allocatable :: vel(:)
  
  real(8) :: x1(imax)
!
  real(8) :: pr, de, demin, uumin, prmin, rr, rin, rcen, roh
!  integer, allocatable :: isetp(:,:,:)
  integer :: isetp1
  real(8) :: gfl, bbsq, bsq, bbromax, bbprmax, prromax, bbro, bbpr, prro
  real(8) :: tmp1a, tmp1b, vsq, gflmax
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), &
           gcov3a(1:3,1:3), gcon3a(1:3,1:3), stat=merr)
  allocate(util(1:3), bcon(1:3), stat=merr)
  allocate(ucov(0:3), ucon(0:3), bbcov(0:3), bbcon(0:3), &
           bcon1(0:3), bcov1(0:3), vel(1:3), stat=merr)
! 
!  allocate(isetp(is1:ie1,js1:je1,ks1:ke1), stat=merr)
!
!
!  do k=ks1+nm1,ke1-nm1
!    do j=js1+nm1,je1-nm1
!      do i=is1+nm1,ie1-nm1  
!        isetp(i,j,k)=0
!      enddo
!    enddo
!  enddo
  isetp1=0 
  gflmax=20.d0
  !
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
         
        rr=x1(i)
        rin=1.d0
        rcen=25.d0

!- copy of primitive variables and conserved variables-! 
        de=uri(1,i,j,k)
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)        
        pr=uri(5,i,j,k)
        bcon(1)=uri(7,i,j,k)
        bcon(2)=uri(8,i,j,k)
        bcon(3)=uri(9,i,j,k)
!        
        bcon1(0)=0.d0
        bcon1(1)=bcon(1)
        bcon1(2)=bcon(2)
        bcon1(3)=bcon(3)
!       
!- copy of metric terms -!
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo
!
        call cal3met(gcov1,gcon1,gcov3a,gcon3a)
        
!- Check Lorentz factor -!
        call calgfl(util,gfl,gcov1) 
!
        vel(1)=util(1)*(1./gfl)
        vel(2)=util(2)*(1./gfl)
        vel(3)=util(3)*(1./gfl)        
!
        vsq=0.d0
        do m=1,3
          do n=1,3
            vsq=vsq+gcov3a(m,n)*vel(m)*vel(n)     
          enddo
        enddo        
!
        if(vsq .ge. 1.d0 .or. gfl .gt. gflmax) then
          vsq=1.-1./(gflmax*gflmax)
        endif
!
        gfl=1./sqrt(1.-vsq)
!          
        util(1)=vel(1)*gfl
        util(2)=vel(2)*gfl
        util(3)=vel(3)*gfl        
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
!        
!- check roh*gamma^2-p < flc*B^2 
        call lower(bcon1,bcov1,gcov1)
        bsq=bcon1(1)*bcov1(1)+bcon1(2)*bcov1(2)+bcon1(3)*bcov1(3)

95      continue
        
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          roh=de+(gam/(gam-1.0))*pr
        elseif(ieos .eq. 1) then
          roh=(5./2.)*pr+sqrt((9./4.)*pr**2+de**2)
        endif
         
        tmp1a=roh*gfl*gfl-pr
        tmp1b=abs(0.05*bsq)
        if(tmp1a .le. tmp1b) then
          de=1.05*de
          pr=1.05*pr
!          write(*,*) "roh*gf^2-p < fluc*B^2 at ", i, j, k
          go to 95
        endif 
!
        if(de .ne. uri(1,i,j,k)) then
          uri(1,i,j,k)=de
          isetp1=1
        elseif(pr .ne. uri(5,i,j,k)) then
          uri(5,i,j,k)=pr
          isetp1=1
        endif
!        
!- set min & max values -!
        
        demin=dmin*(rr**(-3./2.))
        prmin=(gam-1.)*pmin*(rr**(-5./2.))

!        demin=dmin
!        prmin=(gam-1.)*pmin

!        demin=dmin*exp(-3.*rr/rcen)
!        prmin=kpol*(demin**gam)
!        prmin=pmin
!
        bbromax=50.d0
        bbprmax=2.d3
        prromax=1.d6
!
        if( pr.lt.prmin ) then
          uri(5,i,j,k)=prmin
!          isetp(i,j,k)=1
          isetp1=1
        endif
!
        if( pr.gt.pmax ) then
          uri(5,i,j,k)=pmax
!          isetp(i,j,k)=2
          isetp1=1
        endif

        if( de.lt.demin ) then          
          uri(1,i,j,k)=demin
!          isetp(i,j,k)=3
          isetp1=1
        endif   
!
        prro=uri(5,i,j,k)*(1./uri(1,i,j,k))
!
        if(prro .gt. prromax ) then
          uri(1,i,j,k)=uri(5,i,j,k)/prromax 
          isetp1=1
        endif
!
!- check b^2/ro, b^2/pg
!        
!        if(bbsq .lt. 1.d-20) then
!          bbro=0.d0
!          bbpr=0.d0
!        else
          bbro=bbsq*(1./uri(1,i,j,k))
          bbpr=bbsq*(1./uri(5,i,j,k))
!        endif
!        
!        if( pr .lt. prmin .or. pr .gt. pmax .or. de .lt. demin) then
!           
!        if( pr .lt. prmin .or. pr .gt. pmax .or. de .lt. demin .or. &
!             bbro .gt. bbromax .or. bbpr .gt. bbprmax &
!             .or. prro .gt. prromax) then
       
!
        if( bbro .gt. bbromax ) then
          uri(1,i,j,k)=bbsq/bbromax
          isetp1=1 
        endif
!
        if( bbpr .gt. bbprmax ) then
          uri(5,i,j,k)=bbsq/bbprmax
          isetp1=1 
        endif
!         
!          uri(5,i,j,k)=pr
!          uri(1,i,j,k)=de
          
!
!        endif       
!
      enddo
    enddo
 enddo
!
!  if(isetp1 .ne. 0) then
!    call caluu1(uri,uu,gcov,gcon,detg,x1,is1,ie1,js1,je1,ks1,ke1)
!  endif
!
!  do k=ks1,ke1
!    do j=js1,je1
!      do i=is1,ie1
!        if( isetp(i,j,k) .eq. 1) then
!          write(6,*) 'Correction: pressure min'
!          write(6,*) 'pmin =', pmin, i, j, k            
!        elseif( isetp(i,j,k) .eq. 2) then
!          write(6,*) 'Correction: pressure max'
!          write(6,*) 'pmax =', pmax, i, j, k
!        elseif( isetp(i,j,k) .eq. 3) then
!          write(6,*) 'Correction: density min'
!          write(6,*) 'dmin =', dmin, i, j, k
!        endif
!      enddo
!    enddo
!  enddo  
!
!  deallocate(isetp, stat=merr)
  deallocate(gcov1, gcon1, gcov3a, gcon3a, util, bcon, ucov, ucon, &
             bbcov, bbcon, bcon1, bcov1, vel, stat=merr)
!
  return
end subroutine pminmax2b
!
!***********************************************************************
!                   ARTIFICIAL DAMP  
!***********************************************************************
!-----------------------------------------------------------------------
subroutine damp4(uu,uri,uri0,gcov,gcon,detg,x1,x3,is1,ie1,js1,je1,ks1,ke1)
!-----------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv, metric, adamp, rdamp, &
                   xmin, xmax
  implicit none
  
  integer :: i, j, k, is1, ie1, js1, je1, ks1, ke1

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri0(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)

  real(8) :: x1(imax), x3(kmax)
  real(8) :: rr, adt, de, util1, util2, util3, pr, b1, b2, b3
!
  real(8) :: drdamp, dx1a

  drdamp=0.4d0
  dx1a=(xmax-xmin)/float(imax)

  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
        if(rdamp .gt. 0.d0) then
          if(metric.eq.3 .or. metric.eq.103 .or. metric.eq.203) then
            rr=x1(i)
          else
            rr=sqrt(x1(i)**2+x3(k)**2)
          endif
            adt=adamp*0.5*(1.0+tanh((rdamp-rr)*(0.3/dx1a)))
!            adt=adamp*0.5*(1.0+tanh((rdamp-rr)/drdamp))

            de= uri(1,i,j,k)-adt*(uri(1,i,j,k)-uri0(1,i,j,k))
            util1= uri(2,i,j,k)-adt*(uri(2,i,j,k)-uri0(2,i,j,k))
            util2= uri(3,i,j,k)-adt*(uri(3,i,j,k)-uri0(3,i,j,k))
            util3= uri(4,i,j,k)-adt*(uri(4,i,j,k)-uri0(4,i,j,k))
            pr= uri(5,i,j,k)-adt*(uri(5,i,j,k)-uri0(5,i,j,k))
!            b1= uri(7,i,j,k)-adt*(uri(7,i,j,k)-uri0(7,i,j,k))
!            b2= uri(8,i,j,k)-adt*(uri(8,i,j,k)-uri0(8,i,j,k))
!            b3= uri(9,i,j,k)-adt*(uri(9,i,j,k)-uri0(9,i,j,k))
            b1=uri(7,i,j,k)
            b2=uri(8,i,j,k)
            b3=uri(9,i,j,k)
          
            uri(1,i,j,k)=de
            uri(2,i,j,k)=util1
            uri(3,i,j,k)=util2
            uri(4,i,j,k)=util3
            uri(5,i,j,k)=pr
            uri(7,i,j,k)=b1
            uri(8,i,j,k)=b2
            uri(9,i,j,k)=b3
          
        endif
      enddo
    enddo
  enddo
      
  if(rdamp .gt. 0.d0) then
    call caluu2a(uri,uu,gcov,gcon,detg,is1,ie1,js1,je1,ks1,ke1)
  endif
!
  return
end subroutine damp4
!
!-----------------------------------------------------------------------
subroutine damp4a(uu,uri,uri0,gcov,gcon,detg,x1,x3,is1,ie1,js1,je1,ks1,ke1)
!-----------------------------------------------------------------------
!    set damping zone in z-direction (for jet propagation simulation)
!
  use pram, only : imax, jmax, kmax, nv, metric, adamp, pi
  implicit none
  
  integer :: i, j, k, is1, ie1, js1, je1, ks1, ke1

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri0(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: alpha(is1:ie1,js1:je1,ks1:ke1)
  real(8) :: beta(0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)

  real(8) :: x1(imax), x3(kmax)
  real(8) :: zz, zz1, zz2, adt, de, util1, util2, util3, pr, b1, b2, b3, tmp1
!
  zz1=-0.95d0
  zz2=-0.025d0

  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
        if(adamp .gt. 0.d0) then
          if(metric.eq.3 .or. metric.eq.103 .or. metric.eq.203) then
            zz=x1(i)*cos(x3(k))
          else
            zz=x3(k)
          endif
            tmp1=(zz-zz1)/(zz2-zz1)
            adt=adamp*0.5*(cos(tmp1*pi)+1.0)

            de= uri(1,i,j,k)-adt*(uri(1,i,j,k)-uri0(1,i,j,k))
            util1= uri(2,i,j,k)-adt*(uri(2,i,j,k)-uri0(2,i,j,k))
            util2= uri(3,i,j,k)-adt*(uri(3,i,j,k)-uri0(3,i,j,k))
            util3= uri(4,i,j,k)-adt*(uri(4,i,j,k)-uri0(4,i,j,k))
            pr= uri(5,i,j,k)-adt*(uri(5,i,j,k)-uri0(5,i,j,k))
!            b1= uri(7,i,j,k)-adt*(uri(7,i,j,k)-uri0(7,i,j,k))
!            b2= uri(8,i,j,k)-adt*(uri(8,i,j,k)-uri0(8,i,j,k))
!            b3= uri(9,i,j,k)-adt*(uri(9,i,j,k)-uri0(9,i,j,k))
            b1=uri(7,i,j,k)
            b2=uri(8,i,j,k)
            b3=uri(9,i,j,k)
          
            uri(1,i,j,k)=de
            uri(2,i,j,k)=util1
            uri(3,i,j,k)=util2
            uri(4,i,j,k)=util3
            uri(5,i,j,k)=pr
            uri(7,i,j,k)=b1
            uri(8,i,j,k)=b2
            uri(9,i,j,k)=b3
          
        endif
      enddo
    enddo
  enddo
      
  if(adamp .gt. 0.d0) then      
    call caluu2a(uri,uu,gcov,gcon,detg,is1,ie1,js1,je1,ks1,ke1)
  endif
  
  return
end subroutine damp4a
!
