!--------------------------------------------------------------------
subroutine md1dtest(uri,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
!         model = 1; 1D Test for debugging

  use pram, only : imax, jmax, kmax, kmax, nv, gam, c0
  implicit none

  integer :: i, j, k, is1, ie1, js1, je1, ks1, ke1
  integer :: merr

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: de, pr
  real(8), allocatable :: util(:), bcon(:)
!
!-allocate variables -!
  allocate(util(1:3), bcon(1:3), stat=merr)  

!--------------------------------------------------------------------
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
  
!        gam=2.d0
          de=1.d0
          pr=0.01d0
          util(1)=0.d0
          util(2)=0.d0
          util(3)=0.d0
          bcon(1)=0.d0
          bcon(2)=0.d0 
          bcon(3)=0.0d0
!          
!- set primitive variables -!
!          
          uri(1,i,j,k)=de
          uri(2,i,j,k)=util(1)
          uri(3,i,j,k)=util(2)
          uri(4,i,j,k)=util(3)
          uri(5,i,j,k)=pr
          uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
          uri(7,i,j,k)=bcon(1)
          uri(8,i,j,k)=bcon(2)          
          uri(9,i,j,k)=bcon(3)
!
      enddo
    enddo
  enddo
!
  deallocate(util, bcon, stat=merr)
!  
  return
end subroutine md1dtest
!
!--------------------------------------------------------------------
subroutine mdshoctub(uri,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
!         model = 2; 1D Shock Tube Test Problem

  use pram, only : imax, jmax, kmax, kmax, nv, gam, c0
  implicit none

  integer :: i, j, k, is1, ie1, js1, je1, ks1, ke1
  integer :: merr

    real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: de, pr
  real(8), allocatable :: util(:), bcon(:)
!
!-allocate variables -!
  allocate(util(1:3), bcon(1:3), stat=merr)   

!--------------------------------------------------------------------
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
  
!        gam=2.d0
        if(i .le. imax/2) then
          de=1.d0
          pr=1.d0
          util(1)=0.d0
          util(2)=0.d0
          util(3)=0.d0
          bcon(1)=0.5d0
          bcon(2)=1.0d0 
          bcon(3)=0.0d0
        else
          de=0.125d0
          pr=0.1d0
          util(1)=0.d0
          util(2)=0.d0
          util(3)=0.d0
          bcon(1)=0.5d0
          bcon(2)=-1.0d0
          bcon(3)=0.0d0
        endif
!          
!- set primitive variables -!
!          
          uri(1,i,j,k)=de
          uri(2,i,j,k)=util(1)
          uri(3,i,j,k)=util(2)
          uri(4,i,j,k)=util(3)
          uri(5,i,j,k)=pr
          uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
          uri(7,i,j,k)=bcon(1)
          uri(8,i,j,k)=bcon(2)          
          uri(9,i,j,k)=bcon(3)
!
      enddo
    enddo
  enddo
!
  deallocate(util, bcon, stat=merr)
!  
  return
end subroutine mdshoctub
!--------------------------------------------------------------------
subroutine md1dacc(uri,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
!         model = 1; 1D GRHD Test for debugging (free-fall)

  use pram, only : imax, jmax, kmax, kmax, nv, gam, c0
  implicit none

  integer :: i, j, k, is1, ie1, js1, je1, ks1, ke1
  integer :: merr
!
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: de, pr
  real(8), allocatable :: util(:), bcon(:)
!
!-allocate variables -!
  allocate(util(1:3), bcon(1:3), stat=merr)  
!
!--------------------------------------------------------------------
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
  
!        gam=2.d0
          de=1.d-1
          pr=1.d-3
          util(1)=0.d0
          util(2)=0.d0
          util(3)=0.d0
          bcon(1)=0.d0
          bcon(2)=0.d0 
          bcon(3)=0.0d0
!          
!- set primitive variables -!
!          
          uri(1,i,j,k)=de
          uri(2,i,j,k)=util(1)
          uri(3,i,j,k)=util(2)
          uri(4,i,j,k)=util(3)
          uri(5,i,j,k)=pr
          uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
          uri(7,i,j,k)=bcon(1)
          uri(8,i,j,k)=bcon(2)          
          uri(9,i,j,k)=bcon(3)
!
      enddo
    enddo
  enddo

  deallocate(util, bcon, stat=merr)  
  
  return
end subroutine md1dacc
!
!--------------------------------------------------------------------
subroutine mdeqcor1(uri,x1,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
!         model = 1; 1D GRHD Test for debugging (free-fall)

  use pram, only : imax, jmax, kmax, kmax, nv, gam
  implicit none

  integer :: i, j, k, is1, ie1, js1, je1, ks1, ke1
  integer :: merr

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: de, pr
  real(8), allocatable :: util(:), bcon(:)
  real(8) :: x1(imax)
  real(8) :: rdcen1, ck1, ee0, dd0, cssq, rr
!
!-allocate variables -!
  allocate(util(1:3), bcon(1:3), stat=merr)  
!
!--------------------------------------------------------------------
! Parameter
!
  rdcen1=10.0d0
  ck1=0.1d0
  ee0=0.1d0
  dd0=1.d0

  cssq=gam*(gam-1.0)*ee0/dd0
!
!--------------------------------------------------------------------
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1  

        rr=x1(i)
        de=dd0*exp(-3.0*rr/rdcen1)

!        pr=ck1*de**gam
        pr=ck1*de**((gam+1.0)/gam)
!        pr=cssq/gam*de
!        pr=(gam-1.0)*e0*de

        util(1)=0.d0
        util(2)=0.d0
        util(3)=0.d0

        bcon(1)=0.d0
        bcon(2)=0.d0
        bcon(3)=0.d0
!
!          
!- set primitive variables -!
!          
        uri(1,i,j,k)=de
        uri(2,i,j,k)=util(1)
        uri(3,i,j,k)=util(2)
        uri(4,i,j,k)=util(3)
        uri(5,i,j,k)=pr
        uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
        uri(7,i,j,k)=bcon(1)
        uri(8,i,j,k)=bcon(2)          
        uri(9,i,j,k)=bcon(3)
!        
      enddo
    enddo
  enddo
!
  deallocate(util, bcon, stat=merr)
!
  return
end subroutine mdeqcor1
!
!--------------------------------------------------------------------
subroutine mdeqcor2(uri,x1,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
!         model = 1; 1D GRHD Test for debugging (free-fall)

  use pram, only : imax, jmax, kmax, kmax, nv, gam
  implicit none

  integer :: i, j, k, is1, ie1, js1, je1, ks1, ke1
  integer :: merr
!  
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: de, pr
  real(8), allocatable :: util(:), bcon(:)
  real(8) :: x1(imax)
!
  real(8) :: dd0, rr, ee0, pp0, rin
!
!-allocate variables -!
  allocate(util(1:3), bcon(1:3), stat=merr)
!  
!--------------------------------------------------------------------
! Parameter
!
  dd0=0.0001d0
  pp0=0.01d0*dd0
  rin=1.d0
!
!--------------------------------------------------------------------
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1  

        rr=x1(i)
        de=dd0*(rr/rin)**(-1.5d0)
!        ee0=0.1*dd0*rr**(-5./2.)
        ee0=pp0*(rr/rin)**(-2.5d0)
!        pr=(gam-1.0)*ee0
        pr=(gam-1.d0)*pp0*(rr/rin)**(-2.5d0)

!        util(1)=0.d0
        util(1)=-1.d0/(rr*rr)
        util(2)=0.d0
        util(3)=0.d0

        bcon(1)=0.d0
        bcon(2)=0.d0
        bcon(3)=0.d0
!          
!- set primitive variables -!
!          
        uri(1,i,j,k)=de
        uri(2,i,j,k)=util(1)
        uri(3,i,j,k)=util(2)
        uri(4,i,j,k)=util(3)
        uri(5,i,j,k)=pr
        uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
        uri(7,i,j,k)=bcon(1)
        uri(8,i,j,k)=bcon(2)          
        uri(9,i,j,k)=bcon(3)
!        
      enddo
    enddo
  enddo
!
  deallocate(util, bcon, stat=merr)
!    
  return
end subroutine mdeqcor2
!
!--------------------------------------------------------------------
subroutine mdmagsph1(uri,gcov,gcon,detg,x1,x2,x3,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
!   monopole magnetosphere of BH (Komissarov 2004)
!   model=20  
!
  use pram, only : imax, jmax, kmax, nv, gam, c0, rbh, metric, akm, pi, kpol
  implicit none

  integer :: i, j, k, m, n, is1 ,ie1, js1, je1, ks1, ke1
  integer :: merr, iflg1, irs, ichk

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)

  real(8) :: x1(imax), x2(jmax), x3(kmax)

  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: bcon(:), bcon1(:), bcov1(:), ucon(:), ucov(:), &
                          bbcov(:), bbcon(:)
  real(8), allocatable :: util(:), beta1(:)

  real(8) :: rr, th, alpha1, b0, flc, bsq, bbsq, gfl, de, pr, hrel

!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), stat=merr)
  allocate(bcon(1:3), bcon1(0:3), bcov1(0:3), util(1:3), beta1(1:3), stat=merr)
  allocate(ucov(0:3), ucon(0:3), bbcov(0:3), bbcon(0:3), stat=merr)
!
! 
!--------------------------------------------------------------------
!- Parameter -!
!
  b0=1.d0
  flc=0.08d0
  hrel=1.5d0
!  
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1 

        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo
       
        alpha1=1./sqrt(-gcon1(0,0))
        beta1(1)=alpha1*alpha1*gcon1(0,1)
        beta1(2)=alpha1*alpha1*gcon1(0,2)
        beta1(3)=alpha1*alpha1*gcon1(0,3)
        
        rr=x1(i)
        th=x3(k)
!
!        util(1)=-beta1(1)/alpha1
!        util(2)=-beta1(2)/alpha1
!        util(3)=-beta1(3)/alpha1
        util(1)=0.d0
        util(2)=0.d0
        util(3)=0.d0
!
        bcon(1)=b0*sin(th)/(detg(i,j,k))
!        bcon(1)=b0*(1./(detg(i,j,k)*rr**2))        
!        bcon(1)=0.d0
        bcon(2)=0.d0
        bcon(3)=0.d0
!        
        bcon1(0)=0.d0
        bcon1(1)=bcon(1)
        bcon1(2)=bcon(2)
        bcon1(3)=bcon(3)
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
!- cal of covariant 3-magnetic field -! !
        call lower(bcon1,bcov1,gcov1)
!
!- cal of B-field square (3-vector) -!
        bsq=bcon1(1)*bcov1(1)+bcon1(2)*bcov1(2)+bcon1(3)*bcov1(3)
!        
        de=flc*bsq
!        pr=kpol*de**gam
        pr=flc*bsq

!        if(hrel/alpha1 .le. 1.d0) then
!          write(*,*) 'hrel/alpha < 1 at', rr 
!          stop 
!        endif   

!        pr=1.d-2/rr
!        pr=0.5d0*de**gam
!        pr=((gam-1.)/gam)*de*((hrel/alpha1)-1.)        
!       de=(gam*pr/(gam-1.))*(1./((hrel/alpha1)-1.))
        !          
!- set primitive variables -!
!          
        uri(1,i,j,k)=de
        uri(2,i,j,k)=util(1)
        uri(3,i,j,k)=util(2)
        uri(4,i,j,k)=util(3)
        uri(5,i,j,k)=pr
        uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
        uri(7,i,j,k)=bcon(1)
        uri(8,i,j,k)=bcon(2)          
        uri(9,i,j,k)=bcon(3)
!
     enddo
    enddo
  enddo
!
  deallocate(gcov1, gcon1, bcon, bcon1, bcov1, util, ucov, ucon, &
             bbcov, bbcon, beta1, stat=merr)
 
  return
end subroutine mdmagsph1
!
!--------------------------------------------------------------------
subroutine mdffcor1(uri,gcov,gcon,detg,x1,x2,x3,x1b,x2b,x3b, &
                    is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
!   Bondi Accretion Model based on Hawley, Smarr, Wilson (1984)
!   + angular momentum (Proga et al. 2003)
!
  use pram, only : imax, jmax, kmax, nv, gam, c0, rbh, metric, akm, pi
  implicit none

  integer :: i, j, k, m, n, is1 ,ie1, js1, je1, ks1, ke1
  integer :: merr, iflg1, irs, ichk

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)
!  
  real(8) :: de, pr

  real(8) :: x1(imax), x2(jmax), x3(kmax) 
  real(8) :: x1b(imax), x2b(jmax), x3b(kmax)

  real(8), allocatable :: util(:), bcon(:)
  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:)
  real(8), allocatable :: ffur(:), ffte(:)
  real(8), allocatable :: utiln(:)

  real(8) :: dd0, vp, rc, gam1, urcsq, urc, vvcsq, tec, cnst1, cnst2, &
             rr, th, ur, te, al0, ffunc, uphi, b0, tmp1, tmp1a, tmp1b, &
             gfl, bbsq, beta, beta_act, bnorm, beta_min1, bsq_max, &
             sig_act, sig_max1
  real(8) :: x1aa, x2aa, x3aa, x3t
  real(8) :: del,sig, gamrr, alpha1, det1   
!
!-allocate variables -!
  allocate(util(1:3), bcon(1:3), stat=merr)
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), stat=merr)
  allocate(ucov(0:3), ucon(0:3), bbcov(0:3), bbcon(0:3), stat=merr)
  allocate(ffur(imax), ffte(imax), utiln(1:3), stat=merr)
!
! 
!--------------------------------------------------------------------
!  Parameter
!
  gam1=1.d0/(gam-1.d0)
!
  dd0=1.0d0 !- normalized density -!
  vp=1.d0 !- velocity correction -!
  rc=8.d0 !- sonic point radius -!
  al0=0.0d0 !- angular momentum paramter -!
!  b0=0.5084d0 !- magnetic field strength parameter -!
  b0=2.0d0
!
  beta_min1=1.d10
  beta_act=1.d10
  bsq_max=0.d0
  sig_max1=1.d-10
  sig_act=1.d-10
!  
!
  urcsq=1./(2.*rc) !- (u^r_c)^2
  urc=-sqrt(urcsq) !
  vvcsq=(urcsq)*(1./(1.-3.*urcsq))
!  tec=gam1*vvcsq/((1.+gam1)*(1.-gam1*vvcsq))
  tec=gam1*vvcsq*(1./((1.+gam1)-gam1*(1.+gam1)*vvcsq))
! 
  cnst1=(tec**gam1)*urc*(rc**2)
!  cnst2=(1.-(2.*rbh/rc)+urcsq)*(1.+(1.+gam1)*tec)**2
  tmp1a=cnst1**2
  tmp1b=(rc**4)*(tec**(2.*gam1))
  tmp1=tmp1a/tmp1b
  cnst2=(1.-(2./rc)+tmp1)*((1.+(1.+gam1)*tec)**2)
!
  iflg1=0
  irs=1
!      
  do i=1,imax
    if(iflg1 .eq. 0) then
      if(x1(i) .ge. rc) then
        iflg1=1
        irs=i
      endif
    endif
  enddo
!
!======================================================================@
!     Calculation of alpha and ur by Newton-Raphson
!======================================================================@
!
  do i=1,imax
    ffur(i)=urc
    ffte(i)=tec
  enddo
!
! ==== rr > rcs ====  
!
  do i=irs,imax,1
  
    rr=x1(i)
    
    if(i .gt. irs) then
      ur=ffur(i-1)+0.01d0
      te=ffte(i-1)
    else
      ur=urc
      te=tec
    endif
       
    call calws1(rr,ur,te,cnst1,cnst2)
    
    ffur(i)=ur
    ffte(i)=te

    
  enddo
!
! ==== rr < rcs ====
!
  do i=irs-1,1,-1

    rr=x1(i)
    th=x3(k)

    if(i .eq. irs-1) then
      ur=urc
      te=tec
    else
      ur=ffur(i+1)-0.01d0
      te=ffte(i+1)
    endif

    call calws1(rr,ur,te,cnst1,cnst2)

    ffur(i)=ur
    ffte(i)=te

  enddo
!
  ichk=0
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1 
 
        rr=x1(i)
        th=x3(k)
!
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo
!
        alpha1=1./sqrt(-gcon1(0,0))
!
        del=rr*rr-2.*rr
        sig=rr*rr
        gamrr=sig/del
        det1=-rr*rr*sin(th)*sin(th)
!
!- set Bondi-flow structure -!
!        
        de=dd0*ffte(i)**gam1
        pr=de*ffte(i)
        util(1)=vp*ffur(i)
        util(2)=0.d0
        util(3)=0.d0
!- set B-field structure -!
!       
!        bcon(1)=0.d0
!        bcon(1)=b0/rr**2
!        bcon(1)=alpha1*b0*sin(th)/detg(i,j,k)
!        bcon(1)=(sqrt(b0*2.*dd0*tec**(gam1+1.))*(rc/rr)**2)/detg(i,j,k)        
        bcon(1)=b0*(1./(detg(i,j,k)*rr**2))
!        bcon(1)=b0*(1./(det1*rr**2))
        bcon(2)=0.d0
        bcon(3)=0.d0

!
!- Angular momentum -!
!
!        ffunc=1.- abs(cos(th))
!        ffunc=1.- (cos(th)**10.0)
!        if(th .le. theta0 .or. th .ge. pi-theta0)
!          ffunc=0.d0
!        else
!          ffunc=1.d0
!        endif
!
!        uphi=al0*ffunc/rr 
!        util(2)=uphi                
!
!- transform from BL coordinate to KS coordinate -!        
        if(metric .eq. 303 .or. metric .eq. 403) then
!
          x1aa=x1(i)
          x2aa=x2(j)
          x3aa=x3(k)   
          x3t=x3b(k)
!
          call transbl2ks(util,utiln,x1aa,x2aa,x3aa,x3t)
!
          util(1)=utiln(1)
          util(2)=utiln(2)
          util(3)=utiln(3)
!
        endif
!          
!- set primitive variables -!
!          
        uri(1,i,j,k)=de
        uri(2,i,j,k)=util(1)
        uri(3,i,j,k)=util(2)
        uri(4,i,j,k)=util(3)
        uri(5,i,j,k)=pr
        uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
        uri(7,i,j,k)=bcon(1)
        uri(8,i,j,k)=bcon(2)          
        uri(9,i,j,k)=bcon(3)
!       
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo
!
        call calgfl(util,gfl,gcov1)
        call calucon(util,ucon,gfl,gcon1)
        call lower(ucon,ucov,gcov1)
        call calbbcon(bcon,bbcon,ucov,ucon)  
        call lower(bbcon,bbcov,gcov1)
        call calbbsq(bbcon,bbcov,bbsq)
!
        if(bbsq .ne. 0.d0) then
           
          beta_act=2.*pr*(1./bbsq)
          sig_act=bbsq*(1./de) 
        endif
        if(beta_act .lt. beta_min1) then
          beta_min1=beta_act
        endif
        if(sig_act .gt. sig_max1) then
          sig_max1=sig_act
        endif   
!       
        if(bbsq .gt. bsq_max) then
          bsq_max=bbsq
        endif
!        
!        if(j .eq. 4 .and. k .eq. 4) then
!          write(*,*) 'i, pr=', i, uri(5,i,j,k)
!        endif   

      enddo
    enddo
  enddo
!
  write(*,*) "beta_min, sig_max=", beta_min1, sig_max1
!
  deallocate(gcov1, gcon1, ucov, ucon, bbcov, bbcon, &
             ffur,ffte, util, bcon, utiln, stat=merr)
 
  return
end subroutine mdffcor1
!
!-----------------------------------------------------------------------
subroutine calws1(rr,ur,te,cnst1,cnst2)
!----------------------------------------------------------------------- 
  use pram, only : gam, c0, rbh, iter
  implicit none

  integer :: nnn
  real(8) :: rr, ur, te, cnst1, cnst2, gam1
  real(8) :: f1, g1, dfx, dfy, dgx, dgy, det, dx1, dy1, delta, deltaend
  real(8) :: ursq, rrsq, tmp3a, tmp3b, xx_old, yy_old, tmp1a, tmp1b 

!     Parameter for Bondi Accretion

  delta=1.d0
  deltaend=1.d-8
  gam1=1.d0/(gam-1.d0)

  do nnn=0,iter
    rrsq=rr*rr
    ursq=ur*ur

    f1=ur*(te**gam1)*rrsq-cnst1
    tmp3a=(1.+(1.+gam1)*te)**2
    tmp3b=1.-(2.*rbh/rr)+ursq
    tmp1a=cnst1**2
    tmp1b=(rr**4)*te**(2.*gam1)
!    tmp3b=1.-(2.*rbh/rr)+tmp1a/tmp1b
    g1=tmp3a*tmp3b-cnst2
      
    dfx=(te**gam1)*rrsq
    dfy=ur*gam1*(te**(gam1-1.))*rrsq
    dgx=2.*tmp3a*ur
!    dgx=0.d0
    dgy=2.*(1.+gam1)*(1.+(1.+gam1)*te)*tmp3b
!    dgy=2.*(1.+gam1)*(1.+(1.+gam1)*te)*tmp3b &
!       -(2.*gam1*(cnst1**2)/(rr**4))*tmp3a*te**(-2.*gam1-1.)

    det=dfx*dgy-dfy*dgx

!      write(*,*) 'f1, g1=', f1, g1
!      write(*,*) 'dfx,dfy,dgx,dgy=', dfx, dfy, dgx, dgy

    if(det .eq. 0.d0) then
      write(6,*) ' >> Jacobian =0 in grbondi'
      write(6,*) ' >> delta, iter =', delta, nnn,', at x1:',rr
 !       stop
    endif

    if(det .eq. 0.d0) then
      dx1=0.d0
      dy1=0.d0
      delta=0.d0
    else        
      dx1=(-dgy*f1+dfy*g1)/det
      dy1=(dgx*f1-dfx*g1)/det

      delta=abs((dx1+dy1)/2.0)
    endif

    if(delta .lt. deltaend) then
      goto 101
    endif

    xx_old=ur
    yy_old=te

    ur=ur+dx1
    te=te+dy1

  enddo

  if(delta .gt. deltaend) then 
    write(6,*) ' >> Not convergence in grbondi'
    write(6,*) ' >> delta, iter =', delta, nnn,', at x1:',rr
    stop
  endif
!
!----------------------------------------------------
!     Calculation of primitive variables
!----------------------------------------------------
!
  101 continue 

  return
end subroutine calws1
!
!--------------------------------------------------------------------
subroutine callfish(rmax,clm)
!--------------------------------------------------------------------
  use pram, only : akm
  implicit none

  real(8) :: rmax, clm
  real(8) :: tmp1a, tmp1b, tmp1c, tmp1d, tmp1e, tmp1f, &
             tmp2a, tmp2b, tmp3a, tmp3b
  tmp1a=((akm**2)-2.*akm*sqrt(rmax)+ (rmax**2))
  tmp1b=(-2.*akm*rmax)*tmp1a
  tmp1c=sqrt(2.*akm*sqrt(rmax)+(-3.+rmax)*rmax)
  tmp1d=(akm+(-2.+rmax)*sqrt(rmax))
  tmp1e=((rmax**3) + (akm**2)*(2.+rmax))
  tmp1f=sqrt(1.+(2.*akm/rmax**(1.5))-(3./rmax))

  tmp2a=(rmax**3)*sqrt(2.*akm*sqrt(rmax)+(-3.+rmax)*rmax)
  tmp2b=((akm**2) + (-2.+rmax)*rmax)

  tmp3a=tmp1a*((tmp1b/tmp1c)+tmp1d*tmp1e/tmp1f)
  tmp3b=tmp2a*tmp2b

  clm=tmp3a/tmp3b

  return
end subroutine callfish
!
!--------------------------------------------------------------------
subroutine getmaxrho(rin,rmax,rhomax,clm,kpol1)
!--------------------------------------------------------------------
  use pram, only : akm, gam
  implicit none
!
  real(8) :: rin, rmax, clm, delin, aain, sigin, delmax, aamax, sigmax
  real(8) :: tmp1a, tmp1c, tmp1d, tmp1e, tmp1f, &
             tmp2a, tmp2c, tmp2d, tmp2e, tmp2f
  real(8) :: kpol1, alnh, hm1, rhomax

!- value at rin -!
!  
  delin=rin*rin-2.d0*rin+akm*akm
  aain=(rin*rin+akm*akm)*(rin*rin+akm*akm)-delin*akm*akm
  sigin=rin*rin 
!
!- value at rmax -!
!
  delmax=rmax*rmax-2.d0*rmax+akm*akm
  aamax=(rmax*rmax+akm*akm)*(rmax*rmax+akm*akm)-delmax*akm*akm
  sigmax=rmax*rmax 
!
!
  tmp1a=4.d0*(clm*clm*sigmax*sigmax)*delmax*(1.d0/(aamax*aamax))
  tmp1c=sigmax*delmax*(1.d0/aamax)
  tmp1d=0.5d0*log((1.d0+sqrt(1.+tmp1a))*(1.d0/tmp1c))
  tmp1e=0.5d0*sqrt(1.d0+tmp1a)
  tmp1f=2.d0*akm*rmax*clm*(1.d0/aamax)
!
  tmp2a=4.d0*(clm*clm*sigin*sigin)*delin*(1.d0/(aain*aain))
  tmp2c=sigin*delin*(1.d0/aain)
  tmp2d=0.5d0*log((1.d0+sqrt(1.+tmp2a))*(1.d0/tmp2c))
  tmp2e=0.5d0*sqrt(1.d0+tmp2a)
  tmp2f=2.d0*akm*rin*clm*(1.d0/aain)
!
  alnh=tmp1d-tmp1e-tmp1f-(tmp2d-tmp2e-tmp2f)
  hm1=exp(alnh) -1.0d0

!- maxmum density
  rhomax=(hm1*(gam-1.0d0)/(kpol1*gam))**(1.0d0/(gam-1.0d0))
!  
  return
end subroutine getmaxrho  
!--------------------------------------------------------------------
subroutine mdtorus1(uri,gcov,gcon,detg, &
                    x1,x2,x3,x1b,x2b,x3b,dx1,dx2,dx3, &
                    myrank,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
!
  use pram, only : imax, jmax, kmax, nv, gam, metric, akm, pi, kpol, &
                   dmin, pmin
  implicit none
  include 'mpif.h'

  integer :: i, j, k, m, n, is1 ,ie1, js1, je1, ks1, ke1
  integer :: merr, myrank, jh
   
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:), beta1(:)
  real(8), allocatable :: util(:), utiln(:), bcon(:)
!
  real(8) :: x1(imax), x2(jmax), x3(kmax) 
  real(8) :: x1b(imax), x2b(jmax), x3b(kmax) 
  real(8) :: dx1(imax), dx2(jmax), dx3(kmax)
!
  real(8),allocatable :: alnh(:,:,:), pom(:,:,:)
!
  real(8) :: de, pr
  real(8) :: rin, rmax, kappa, beta, rhomin, prmin, clm, hm1, &
             expm2chi, uphi1, rhomax1, rhomax1a, &
             prmax1, prmax1a, prmax2, prmax2a, rin1,thedge, thedge1, &
             rhomin1, prmin1
  real(8) :: rr, th, sth, cth, del, aa, sig, rr1, th1, &
             thin, sthin, cthin, delin, aain, sigin
  real(8) :: x1aa, x2aa, x3aa, x3t
  real(8) :: tmp1a, tmp1b, tmp1c, tmp1d, tmp1e, tmp1f, &
             tmp2a, tmp2b, tmp2c, tmp2d, tmp2e, tmp2f, &
             tmp3a, tmp3b, tmp3c, tmp3d, tmp_rand
  real(8) :: bsq_max, bsq_max1, bsq_max2, bsq_max2a, beta_act, &
             bnorm, pom1, rhoav, gfl, bbsq, &
             pbeta, pbeta_min, pbeta_min1, pbeta_min2, pbeta_min2a
  real(8) :: tmp4a, tmp4b, tmp4c, discr, alpha1, omg, uut
  real(8) :: dx1a, dx3a, small
  real(8) :: ff1, ff2, qq1, fc1, pow1, pow2, ss1, tt1, kpol1
! 
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), stat=merr)
  allocate(util(1:3), utiln(1:3), bcon(1:3), ucov(0:3), ucon(0:3), &
           bbcov(0:3), bbcon(0:3), beta1(1:3), stat=merr)
  allocate(alnh(is1:ie1,js1:je1,ks1:ke1), pom(is1:ie1,js1:je1,ks1:ke1), &
           stat=merr)

!--------------------------------------------------------------------
!- Parameter -!

!  rin=15.d0
  rin=6.0d0
!  rmax=34.d0
  rmax=12.d0
!  rmax=15.0d0
  thedge=0.02*pi
  thedge1=0.98*pi

  kappa=1.0d-3
  kpol1=1.d-3
  beta=100.d0
!
  pow1=1.d0
  pow2=2.d0
  fc1=0.2d0
!
  rin1=1.d0 
!  rhomin=5.d-4
  rhomin=dmin
  prmin=pmin
!
  rhomax1=0.d0
  prmax1=0.d0 
  prmax2=0.d0
!  
  bsq_max=0.d0
  bsq_max1=0.d0
  bsq_max2=0.d0
  bsq_max2a=0.d0
!  
  pbeta_min=1.d5
  pbeta_min1=1.d5
  pbeta_min2=1.d5
  pbeta_min2a=1.d5
  !
  small=1.d-12
!  
!- calculate Angular momentum l -!
  call callfish(rmax,clm)
  if(myrank .eq. 0) then
    write(*,*) 'rmax, clm=', rmax, clm
  endif
!
!- get maximum density in the torus -!
!  call getmaxrho(rin,rmax,rhomax1a,clm,kappa) 
!  if(myrank .eq. 0) then
!    write(*,*) 'rhomax=', rhomax1a
!  endif  
!
!- initially set minimum plasma everywhere -!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1

        do m=0,3
          do n=0,3
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
        rr=x1(i)
!
        de=rhomin*((rr)**(-1.5d0))
!        de=rhomin*((rr/rin1))**(-1.0d0)
!        pr=prmin*((rr)**(-2.5d0))
        pr=prmin*(gam-1.)*((rr)**(-2.5d0))
!        pr=kpol*(de**gam)
!        pr=1.0d-3*(gam-1.)*de
!
!        de=rhomin
!        pr=(gam-1.)*prmin
!
        util(1)=0.d0
!        util(1)=-1.d0/(rr*rr)
!        util(1)=beta1(1)/alpha1
        util(2)=0.d0
        util(3)=0.d0
        bcon(1)=0.d0
        bcon(2)=0.d0 
        bcon(3)=0.d0
!          
!- set primitive variables -!
!          
        uri(1,i,j,k)=de
        uri(2,i,j,k)=util(1)
        uri(3,i,j,k)=util(2)
        uri(4,i,j,k)=util(3)
        uri(5,i,j,k)=pr
        uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
        uri(7,i,j,k)=bcon(1)
        uri(8,i,j,k)=bcon(2)          
        uri(9,i,j,k)=bcon(3)
!
        alnh(i,j,k)=0.d0
!        
      enddo
    enddo
  enddo
!
!- calculate ln h -!
!
!- value at inner disk radius
!  
  thin=pi/2.
!  sthin=sin(thin)
  sthin=1.d0
!  cthin=cos(thin)
  cthin=0.d0
  delin=rin*rin-2.d0*rin+akm*akm
  aain=(rin*rin+akm*akm)*(rin*rin+akm*akm)-delin*akm*akm
  sigin=rin*rin 
!
  tmp2a=4.d0*(clm*clm*sigin*sigin)*delin*(1.d0/(aain*aain))
!  tmp2b=aain*sthin*aain*sthin
  tmp2c=sigin*delin*(1.d0/aain)
  tmp2d=0.5d0*log((1.d0+sqrt(1.+tmp2a))*(1.d0/tmp2c))
  tmp2e=0.5d0*sqrt(1.d0+tmp2a)
  tmp2f=2.d0*akm*rin*clm*(1.d0/aain)
  
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!
        rr=x1(i)
        th=x3(k)
        sth=sin(th)
        if(abs(sth) .lt. small) then
          if(sth .ge. 0.d0) then
            sth=small
          elseif(sth .lt. 0.d0) then
            sth=-small 
          endif   
        endif
!
        cth=cos(th) 
!
        del=(rr*rr)-(2.d0*rr)+(akm*akm)
        aa=(rr*rr+akm*akm)*(rr*rr+akm*akm)-del*akm*akm*sth*sth
        sig=rr*rr+akm*akm*cth*cth
!
        if(rr .ge. rin .and. th .ge. thedge .and. th .le. thedge1) then
!
 
           
          tmp1a=4.d0*(clm*clm*sig*sig)*del*(1.d0/(aa*sth*aa*sth))
!          tmp1b=aa*sth*aa*sth
          tmp1c=sig*del*(1.d0/aa)
          tmp1d=0.5d0*log((1.d0+sqrt(1.+tmp1a))*(1.d0/tmp1c))

          tmp1e=0.5d0*sqrt(1.d0+tmp1a)
          tmp1f=2.d0*akm*rr*clm*(1.d0/aa)

          alnh(i,j,k)=tmp1d-tmp1e-tmp1f-(tmp2d-tmp2e-tmp2f)

          if(alnh(i,j,k) .le. 0.d0) then
            alnh(i,j,k)=0.d0 
          endif   
!
        else
          alnh(i,j,k)=0.d0
        endif
!
!- set torus structure -!
!
        hm1=exp(alnh(i,j,k)) -1.d0
        
        if(hm1 .gt. 0.d0 .and. rr .ge. rin &
           .and. th .ge. thedge .and. th .le. thedge1) then    
!    
!!          hm1=exp(alnh(i,j,k)) -1.
!          de=((hm1*(gam-1.d0)/(kpol1*gam))**(1.d0/(gam-1.d0))) /rhomax1a
!          pr=(kpol1*(de**gam)) /rhomax1a
          
          de=(hm1*(gam-1.d0)/(kpol1*gam))**(1.d0/(gam-1.d0))
          pr=kpol1*(de**gam)
!
!          de=rhomin*((rr/rin1)**(-1.5d0))
!          de=rhomin*((rr/rin1))**(-1.0d0)
!          pr=(gam-1.d0)*prmin*((rr/rin1)**(-2.5d0))

          util(1)=0.d0
          util(3)=0.d0
!
          expm2chi=sig*sig*del*(1.d0/(aa*aa*sth*sth))
          uphi1=sqrt((-1.d0+sqrt(1.d0+4.d0*clm*clm*expm2chi))/2.d0)
          tmp3a=2.d0*akm*rr*sqrt(1.+uphi1*uphi1)
          tmp3b=sqrt(aa*sig*del)
          tmp3c=sqrt(sig*(1.d0/aa))*uphi1*(1.d0/sth)
          tmp3d=tmp3a*(1.d0/tmp3b)+tmp3c
!
!          omg=-(gcov1(0,2)+gcov1(0,0)*clm)/(gcov1(2,2)+gcov1(0,2)*clm)         
!          uut=sqrt(-1./(gcov1(0,0)+2.*omg*gcov1(0,2)+omg*omg*gcov1(2,2)))

          
          util(2)=tmp3d
!          util(2)=omg*uut
!          util(2)=1.03d0*tmp3d   
!          util(2)=0.d0
!
!- check density & pressure max -!
!
          if(de .gt. rhomax1) then
            rhomax1=de
          endif
          if(pr .gt. prmax1) then
            prmax1=pr
          endif
!
!- transform velocity from BL to KS coordinates -!
!
          if(metric .eq. 303 .or. metric .eq. 403) then
            x1aa=x1(i)
            x2aa=x2(j)
            x3aa=x3(k)   
            x3t=x3b(k)
            call transbl2ks(util,utiln,x1aa,x2aa,x3aa,x3t)
            util(1)=utiln(1)
            util(2)=utiln(2)
            util(3)=utiln(3)
!
          elseif(metric .eq. 203 .and. akm .ne. 0.d0) then

            do m=0,3
              do n=0,3
                gcov1(m,n)=gcov(m,n,i,j,k)
                gcon1(m,n)=gcon(m,n,i,j,k)
              enddo
            enddo

            alpha1=1./sqrt(-gcon1(0,0))
!
            beta1(1)=alpha1*alpha1*gcon1(0,1)
            beta1(2)=alpha1*alpha1*gcon1(0,2)
            beta1(3)=alpha1*alpha1*gcon1(0,3)
!
            ucon(1)=util(1)
            ucon(2)=util(2)
            ucon(3)=util(3)
!
            tmp4a=gcov1(0,0)
            tmp4b=2.d0*(gcov1(0,1)*ucon(1)+gcov1(0,2)*ucon(2) &
                       +gcov1(0,3)*ucon(3))
            tmp4c=1.d0+gcov1(1,1)*ucon(1)*ucon(1)+gcov1(2,2)*ucon(2)*ucon(2) &
                  +gcov1(3,3)*ucon(3)*ucon(3) &
                  +2.d0*(gcov1(1,2)*ucon(1)*ucon(2) &
                  +gcov1(1,3)*ucon(1)*ucon(3)+gcov1(2,3)*ucon(2)*ucon(3))   
            discr=tmp4b*tmp4b-4.d0*tmp4a*tmp4c
            ucon(0)=(-tmp4b-sqrt(discr))*(1.d0/(2.d0*tmp4a))
! 
            gfl=ucon(0)*alpha1
!
            util(1)=ucon(1)+gfl*beta1(1)*(1.d0/alpha1)
            util(2)=ucon(2)+gfl*beta1(2)*(1.d0/alpha1)
            util(3)=ucon(3)+gfl*beta1(3)*(1.d0/alpha1)   
 
          endif
!          
!- set primitive variables -!
!          
          uri(1,i,j,k)=de
          uri(2,i,j,k)=util(1)
          uri(3,i,j,k)=util(2)
          uri(4,i,j,k)=util(3)
          uri(5,i,j,k)=pr
          uri(6,i,j,k)=pr*(1.d0/(de**(gam-1.)))        
!          
        else
!          
!          de=rhomin*((rr/rin1)**(-1.5d0))
!          de=rhomin*((rr/rin1))**(-1.0d0)
!          pr=(gam-1.d0)*prmin*((rr/rin1)**(-2.5d0))
!          pr=1.0d-3*de(i,j,k)**gam
!
!          de=rhomin
!          pr=prmin
!
!          util(1)=0.d0
!          util(1)=-1.d0/(rr*rr)
!          util(2)=0.d0
!          util(3)=0.d0
!
!          uri(1,i,j,k)=de
!          uri(5,i,j,k)=pr
!
        endif
!        
      enddo
    enddo
  enddo
!
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_allreduce(rhomax1,rhomax1a,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
  call mpi_allreduce(prmax1,prmax1a,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
!
  if(myrank .eq. 0) then
    write(*,*) "rhomax1a, prmax1a=", rhomax1a, prmax1a
  endif
!
!
!- Normalize the densities and pressure -!
!
  rhomax1=0.d0
  prmax1=0.d0
  
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1 
        rr=x1(i)
        th=x3(k)

!        rhomin1=rhomin
!        prmin1=(gam-1.)*prmin
        
        rhomin1=rhomin*(rr**(-3.d0/2.d0))
        prmin1=prmin*(gam-1.)*(rr**(-5.d0/2.d0))
!        prmin1=prmin*(rr**(-5.d0/2.d0))
!        prmin1=kpol*(rhomin1**gam)
        
        hm1=exp(alnh(i,j,k)) -1.d0 

        if(hm1 .gt. 0.d0 .and. rr .ge. rin  &
             .and. th .ge. thedge .and. th .le. thedge1) then
!           
          de=uri(1,i,j,k)*(1.d0/rhomax1a)
          pr=uri(5,i,j,k)*(1.d0/rhomax1a)

!          de=uri(1,i,j,k)
!          pr=uri(5,i,j,k)
!!
!!!          de=uri(1,i,j,k)/(1.d-1*rhomax1a)
!!!          pr=uri(5,i,j,k)/(1.d-1*rhomax1a)
!!
!- Put perturbation in pressure -!
!
          call random_number(tmp_rand)
!!          write(*,*) "tmp_rand=", tmp_rand-0.5
          pr=pr*(1.d0+0.05d0*(tmp_rand-0.5d0))


          if(de .le. rhomin1) then
            de=rhomin1
            uri(2,i,j,k)=0.d0
            uri(3,i,j,k)=0.d0
            uri(4,i,j,k)=0.d0
            pr=prmin1
!            uri(6,i,j,k)=pr*(1./(de**(gam-1.)))     
          endif
          if(pr .le. prmin1) then
            de=rhomin1
            uri(2,i,j,k)=0.d0
            uri(3,i,j,k)=0.d0
            uri(4,i,j,k)=0.d0
            pr=prmin1
!            uri(6,i,j,k)=pr*(1./(de**(gam-1.)))        
          endif
!
!-check maximum density & pressure -!
          if(de .gt. rhomax1) then
            rhomax1=de
          endif
          if(pr .gt. prmax1) then
           prmax1=pr 
          endif

          uri(1,i,j,k)=de
          uri(5,i,j,k)=pr
          uri(6,i,j,k)=pr*(1.d0/(de**(gam-1.d0)))
          
!!        else
!!
!!          de=rhomin1
!!          de=rhomin*((rr/rin1))**(-1.0d0)
!!          pr=prmin1
!!          pr=1.0d-3*de(i,j,k)**gam 
         
        endif
!
      enddo
    enddo
  enddo
 !
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_allreduce(rhomax1,rhomax1a,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
  call mpi_allreduce(prmax1,prmax1a,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
  
   prmax2a=prmax1a/rhomax1a
   rhomax1a=1.d0
!
  if(myrank .eq. 0) then
    write(*,*) 'rhomax1a, prmax1a=', rhomax1a, prmax1a
  endif
!
!
!- Set vector potential at cell-center (r-theta frame)-!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
        pom(i,j,k)=0.d0
      enddo
    enddo
  enddo
!
!  do k=ks1+1,ke1-1
!    do j=js1,je1
!      do i=is1+1,ie1-1
!
  do k=ks1+1,ke1
    do j=js1,je1
      do i=is1+1,ie1
!  do k=ks1,ke1
!    do j=js1,je1
!      do i=is1,ie1  

        rr1=x1(i)
        th1=x3(k)
!        rr1=(x1(i+1)+x1(i))/2.
!        th1=(x3(k+1)+x3(k))/2.
!        pom1=0.5*rr1*sin(th1)

!        rhoav=0.25*(uri(1,i-1,j,k)+uri(1,i,j,k-1) &
!                    +uri(1,i+1,j,k)+uri(1,i,j,k+1))
!
! A_phi(i,j,k)=A_phi(i-1/2,j,k-1/2)
!
        rhoav=0.25*(uri(1,i,j,k)+uri(1,i-1,j,k) &
                    +uri(1,i,j,k-1)+uri(1,i-1,j,k-1))
!        rhoav=uri(1,i,j,k)

!- single poloidal loop        
!- A_phi \propto rho_{disk}
        if(rhoav .ne. 0.d0) then
          pom1=((rhoav*(1./rhomax1a)) - 0.2d0)
        else
          pom1=0.d0
        endif
!
!- Large poloidal loop (flipping)
!        
!        qq1=((uri(5,i,j,k)/prmax1)-fc1)
!        if(qq1 .le. 0.d0) then
!          ff1=0.d0
!        else
!          ff1=abs(qq1)*(abs(rr1*sin(th1)))**2
!        endif
!        ff2=sin(log(rr1/ss1)/tt1)  
!!        ff2=1.d0
!        pom1=ff1*ff2
        
        if(pom1 .gt. 0.d0) then
          pom(i,j,k)=pom1
        endif
!
      enddo
    enddo
  enddo
!
!- put boundary data for magnetic field -!
!
  do j=js1,je1
    do i=is1,ie1
      pom(i,j,ks1)=pom(i,j,ks1+1)
!      pom(i,j,ke1)=pom(i,j,ke1-1)
    enddo
  enddo
!
  do k=ks1,ke1
    do j=js1,je1
      pom(is1,j,k)=pom(is1+1,j,k)
!      pom(ie1,j,k)=pom(ie1-1,j,k)
    enddo
  enddo
!
!- Set-up of Magnetic field -!
  do k=ks1+1,ke1-1
    do j=js1,je1
      do i=is1+1,ie1-1
!        
!!        bcon(1)=-(pom(i,j,k)-pom(i,j,k+1)+pom(i+1,j,k)-pom(i+1,j,k+1)) &
!!                   /(2.*(x3(k+1)-x3(k))*detg(i,j,k))
!!
!!        bcon(3)=(pom(i,j,k)+pom(i,j,k+1)-pom(i+1,j,k)-pom(i+1,j,k+1)) &
!!                   /(2.*(x1(i+1)-x1(i))*detg(i,j,k))
!
        tmp1a=0.5d0*(pom(i+1,j,k+1)+pom(i,j,k+1)-pom(i+1,j,k)-pom(i,j,k))
!        tmp1a=pom(i,j,k+1)-pom(i,j,k) !- GR way ?
!        tmp1a=sin(x3(k+1))*pom(i,j,k+1)-sin(x3(k))*pom(i,j,k) !- Newtonian

!        dx3a=x3(k+1)-x3(k)
        if(metric .eq. 403) then
          dx3a=0.5d0*(x3b(k+1)+x3b(k))-0.5*(x3b(k)+x3b(k-1))           
        else   
          dx3a=0.5d0*(x3(k+1)+x3(k))-0.5*(x3(k)+x3(k-1))
        endif
       
        if(tmp1a .ne. 0.d0) then
          bcon(1)=-tmp1a*(1.d0/(dx3a*detg(i,j,k))) !- GR way ?
!          bcon(1)=-tmp1a*(1.d0/(x1(i)*sin(x3(k))*dx3a)) !- Newtonian
!          bcon(1)=-tmp1a*(1.d0/(x1(i)*dx3a)) 
        else
          bcon(1)=0.d0
        endif
!!
        tmp1b=0.5d0*(pom(i+1,j,k+1)+pom(i+1,j,k)-pom(i,j,k+1)-pom(i,j,k))
!        tmp1b=pom(i+1,j,k)-pom(i,j,k) !- GR way ?
!        tmp1b=x1(i+1)*pom(i+1,j,k)-x1(i)*pom(i,j,k) !- Newtonian

!        dx1a=x1(i+1)-x1(i)
        if(metric .eq. 403) then
          dx1a=0.5d0*(x1b(i+1)+x1b(i))-0.5d0*(x1b(i)+x1b(i-1))
        else   
          dx1a=0.5d0*(x1(i+1)+x1(i))-0.5d0*(x1(i)+x1(i-1))
        endif   

        if(tmp1b .ne. 0.d0) then
          bcon(3)=tmp1b*(1.d0/(dx1a*detg(i,j,k))) !- GR way ?
!          bcon(3)=tmp1b*(1./(x1(i)*dx1a)) !- Newtonian
!          bcon(3)=tmp1b*(1./dx1a)  
        else
          bcon(3)=0.d0
        endif
!
!!        if(k .eq. ks1) then
!!          bcon(1)=(pom(i,j,k+1)-pom(i,j,k))/((x3(k+1)-x3(k))*detg(i,j,k))
!!        elseif(k .eq. ke1) then
!!          bcon(1)=(pom(i,j,k)-pom(i,j,k-1))/((x3(k)-x3(k-1))*detg(i,j,k))
!!        else         
!!          bcon(1)=(pom(i,j,k+1)-pom(i,j,k-1))/((x3(k+1)-x3(k-1))*detg(i,j,k))
!!        endif
!
!!        if(i .eq. is1) then
!!          bcon(3)=-(pom(i+1,j,k)-pom(i,j,k))/((x1(i+1)-x1(i))*detg(i,j,k))
!!        elseif(i .eq. ie1) then
!!          bcon(3)=-(pom(i,j,k)-pom(i-1,j,k))/((x1(i)-x1(i-1))*detg(i,j,k))
!!        else         
!!          bcon(3)=-(pom(i+1,j,k)-pom(i-1,j,k))/((x1(i+1)-x1(i-1))*detg(i,j,k))
! !       endif
!
        uri(7,i,j,k)=bcon(1)
        uri(9,i,j,k)=bcon(3)
!
      enddo
    enddo
  enddo
!
!- put boundary data for magnetic field -!
!
  do j=js1,je1
    do i=is1,ie1
      uri(7,i,j,ks1)=uri(7,i,j,ks1+1)
      uri(7,i,j,ke1)=uri(7,i,j,ke1-1)
!
      uri(9,i,j,ks1)=uri(9,i,j,ks1+1)
      uri(9,i,j,ke1)=uri(9,i,j,ke1-1)
    enddo
  enddo
!
  do k=ks1,ke1
    do j=js1,je1
      uri(7,is1,j,k)=uri(7,is1+1,j,k)
      uri(7,ie1,j,k)=uri(7,ie1-1,j,k)
!
      uri(9,is1,j,k)=uri(9,is1+1,j,k)
      uri(9,ie1,j,k)=uri(9,ie1-1,j,k)
    enddo
  enddo
!
!- calculation maximum magnetic field strength -!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
        pr=uri(5,i,j,k)
        bcon(1)=uri(7,i,j,k)
        bcon(2)=uri(8,i,j,k)
        bcon(3)=uri(9,i,j,k)
!       
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo
!
        call calgfl(util,gfl,gcov1)
        call calucon(util,ucon,gfl,gcon1)
        call lower(ucon,ucov,gcov1)
        call calbbcon(bcon,bbcon,ucov,ucon)
        call lower(bbcon,bbcov,gcov1)
        call calbbsq(bbcon,bbcov,bbsq)
!
        if(bbsq .ne. 0.d0) then        
          pbeta=2.*pr/bbsq
        else
          pbeta=1.d5  
        endif

        if(pbeta .lt. pbeta_min) then
          pbeta_min=pbeta           
        endif
       
        if(bbsq .gt. bsq_max) then
          bsq_max=bbsq
        endif
!
      enddo
    enddo
  enddo    
!
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_allreduce(bsq_max,bsq_max1,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
  call mpi_allreduce(pbeta_min,pbeta_min1,1,mpi_double_precision,mpi_min, &
                     mpi_comm_world,merr)
!
!  if(bsq_max1 .eq. 0.d0) then
!    beta_act=0.d0
!  else
!    beta_act=prmax1a/(0.5*bsq_max1)
!  endif   
  beta_act=pbeta_min1 
!
  if(beta .eq. 0.d0) then
    bnorm=0.d0
  else
    bnorm=sqrt(beta_act/beta)
  endif
!
  if(myrank .eq. 0) then
    write(*,*) "beta_act, bnorm=", beta_act, bnorm
  endif
!
!- Normalize magnetic field strength -!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!
        bcon(1)=uri(7,i,j,k)*bnorm            
        bcon(2)=uri(8,i,j,k)*bnorm  
        bcon(3)=uri(9,i,j,k)*bnorm  
!
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
        pr=uri(5,i,j,k)
!
        uri(7,i,j,k)=bcon(1)
        uri(8,i,j,k)=bcon(2)
        uri(9,i,j,k)=bcon(3)
!       
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo        
!
        call calgfl(util,gfl,gcov1)
        call calucon(util,ucon,gfl,gcon1)
        call lower(ucon,ucov,gcov1)
        call calbbcon(bcon,bbcon,ucov,ucon)
        call lower(bbcon,bbcov,gcov1)
        call calbbsq(bbcon,bbcov,bbsq)
!
        if(bbsq .ne. 0.d0) then        
          pbeta=2.*pr/bbsq
        else
          pbeta=1.d5  
        endif

        if(pbeta .lt. pbeta_min2) then
          pbeta_min2=pbeta           
        endif

        if(bbsq .gt. bsq_max2) then
          bsq_max2=bbsq
        endif
!
      enddo
    enddo
  enddo 
!
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_allreduce(bsq_max2,bsq_max2a,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
  call mpi_allreduce(pbeta_min2,pbeta_min2a,1,mpi_double_precision,mpi_min, &
                     mpi_comm_world,merr)
  
!  if(bsq_max2a .gt. 0.d0) then
!    beta_act=prmax1a/(0.5*bsq_max2a)
!  else
!    beta_act=0.d0
!  endif   
   beta_act=pbeta_min2a
!
  if(myrank .eq. 0) then
    write(*,*) "beta_act=", beta_act
  endif    
!
  deallocate(gcov1, gcon1, util, utiln, bcon, ucov, ucon, &
             bbcov, bbcon, beta1, alnh, pom, stat=merr)
!
  return
end subroutine mdtorus1
!
!--------------------------------------------------------------------
subroutine mdtorus2(uri,gcov,gcon,detg, &
                    x1,x2,x3,x1b,x2b,x3b,dx1,dx2,dx3, &
                    myrank,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
!  power-law rotating torus model  (nearly Keplerian)
!  Devilliers et al. 2003, ApJ, 599, 1238
!  model=8  
!
  use pram, only : imax, jmax, kmax, nv, gam, metric, akm, pi, iter, kpol, &
                   dmin, pmin
  implicit none
  include 'mpif.h'

  integer :: i, j, k, m, n, is1 ,ie1, js1, je1, ks1, ke1
  integer :: merr, myrank, jh, nnn

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  
  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:), gcon1in(:,:), gcov1in(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:), beta1(:)
  real(8), allocatable :: util(:), utiln(:), bcon(:)
!
  real(8) :: x1(imax), x2(jmax), x3(kmax) 
  real(8) :: x1b(imax), x2b(jmax), x3b(kmax) 
  real(8) :: dx1(imax), dx2(jmax), dx3(kmax)
!
  real(8),allocatable :: ah(:,:,:), pom(:,:,:)
!
  real(8) :: de, pr  
  real(8) :: kappa, beta, clmin, clm, cqq, cqq1, eta, ckk, calp, calp2, &
             amdain, amda, utin, ut, uut, flin, fl, delta, deltaend, f, df, &
             rhomin, prmin, rhomin1, prmin1, hm1, omg, thedge, thedge1, rcen, &
             rhomax1, rhomax1a, prmax1, prmax1a, prmax2, prmax2a 
  real(8) :: rr, th, sth, cth, del, aa, sig, rr1, th1, &
             rin, thin, sthin, cthin, delin, aain, sigin, rin1
  real(8) :: x1aa, x2aa, x3aa, x3t
  real(8) :: tmp_rand
  real(8) :: bsq_max, bsq_max1, bsq_max2, bsq_max2a, beta_act, &
             bnorm, pom1, rhoav, gfl, bbsq
  real(8) :: tmp1a, tmp1b, dx1a, dx3a
  real(8) :: tmp4a, tmp4b, tmp4c, discr, alpha1
  real(8) :: tmp1, tmp5a1, tmp5a2, tmp5b1, tmp5b2, tmp5b3, tmp5b4, tmp5b5
  real(8) :: small
! 
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), &
           gcon1in(0:3,0:3), gcov1in(0:3,0:3), stat=merr)
  allocate(util(1:3), utiln(1:3), bcon(1:3), ucov(0:3), ucon(0:3), &
           bbcov(0:3), bbcon(0:3), beta1(1:3), stat=merr)
  allocate(ah(is1:ie1,js1:je1,ks1:ke1), pom(is1:ie1,js1:je1,ks1:ke1), &
           stat=merr)

!--------------------------------------------------------------------
!- Parameter -!

  rin=15.d0
  thin=pi/2.d0
  kappa=1.0d-2
!  beta=1.d2
!  beta=50.d0
  beta=100.d0
!  clmin=4.66d0 !a=0.0
!  clmin=4.61d0 !a=0.5
  clmin=4.57d0 !a=0.9 or 0.998
  cqq=1.68d0
!  
  thedge=0.2d0*pi
  thedge1=0.8d0*pi
  rcen=25.d0
!  
  rin1=1.d0
  rhomin=dmin 
!  prmin=(gam-1.)*pmin
  prmin=pmin
  !
  rhomax1=0.d0
  rhomax1a=0.d0
  prmax1=0.d0
  prmax1a=0.d0
  prmax2=0.d0
!
  bsq_max=0.d0
  bsq_max1=0.d0
  bsq_max2=0.d0
  bsq_max2a=0.d0
!
  small=1.d-12  
!
!- initally set minimum plasma everywhere -!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1

        do m=0,3
          do n=0,3
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo

        alpha1=1./sqrt(-gcon1(0,0))
!
        beta1(1)=alpha1*alpha1*gcon1(0,1)
        beta1(2)=alpha1*alpha1*gcon1(0,2)
        beta1(3)=alpha1*alpha1*gcon1(0,3)

        
        rr=x1(i)

!        de=rhomin
!        de=rhomin*exp(-3.*rr/rcen)
        de=rhomin*(rr)**(-3./2.)

        pr=(gam-1.d0)*prmin*((rr)**(-5./2.))
!        pr=kpol*(de**gam)
        pr=prmin

        util(1)=0.d0
!        util(1)=(gcon1(0,1)/gcon1(0,0))*(1.-(1./rr)**4)-beta1(1)/alpha1
!        util(1)=-beta1(1)*(1./alpha1)
        util(2)=0.d0
        util(3)=0.d0
        bcon(1)=0.d0
        bcon(2)=0.d0 
        bcon(3)=0.d0
!          
!- set primitive variables -!
        
        uri(1,i,j,k)=de
        uri(2,i,j,k)=util(1)
        uri(3,i,j,k)=util(2)
        uri(4,i,j,k)=util(3)
        uri(5,i,j,k)=pr
        uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
        uri(7,i,j,k)=bcon(1)
        uri(8,i,j,k)=bcon(2)          
        uri(9,i,j,k)=bcon(3)
!       
      enddo
    enddo
  enddo
!
!- calculation of constant paramter from lin  -!
!- calculation of BL Kerr metric at rin, thin -!
!  
!  sthin=sin(thin)
  sthin=1.d0
!  cthin=cos(thin)
  cthin=0.d0
  delin=rin*rin-2.*rin+akm*akm
  aain=(rin*rin+akm*akm)*(rin*rin+akm*akm)-delin*akm*akm*sthin*sthin
  sigin=rin*rin+akm*akm*cthin*cthin  
!
  do n=0,3
    do m=0,3
      gcon1in(m,n)=0.d0
      gcov1in(m,n)=0.d0 
    enddo
  enddo

  gcon1in(0,0)=-aain*(1./(sigin*delin))               !- g^tt -!
  gcon1in(2,2)=(sigin-2.*rin)*(1./(delin*sigin*sthin*sthin)) !- g^phi phi -!
  gcon1in(0,2)=-(2.*akm*rin)*(1./(sigin*delin))           !- g^t phi -!
!
!  gcov1in(0,0)=-(sigin-2.*rin)*(1./sigin)               !- g_tt -!
!  gcov1in(2,2)=(aain*(1./sigin))*sthin*sthin !- g_phi phi -!
!  gcov1in(0,2)=-(2.*akm*rin*sthin*sthin)*(1./sigin)  !- g_t phi -!  
!
!- calculation of lambda at rin -!
  tmp1b=(clmin*gcon1in(0,0)-clmin*clmin*gcon1in(0,2)) &
        /(gcon1in(0,2)-clmin*gcon1in(2,2))
  amdain=sqrt(abs(tmp1b))
!
!- calculation of ut at rin -!
!  tmp1a=gcov1in(0,2)*gcov1in(0,2)-gcov1in(0,0)*gcov1in(2,2)
!  tmp1b=gcov1in(0,0)*clm*clm+2.*clm*gcov1in(0,2)+gcov1in(2,2)
!  if(tmp1b .gt. 0.d0 .and. tmp1a .gt. 0.d0) then
!    utin=-sqrt(tmp1a/tmp1b)
!  else
!    utin=-1.d0  
!  endif
!
  tmp1a=abs(gcon1in(0,0)-2.*clmin*gcon1in(0,2)+clmin*clmin*gcon1in(2,2))
  utin=-1./sqrt(tmp1a) 
! 
!- calculation of eta, k, alpha -!
  eta=clmin*(1./(amdain**(2.-cqq)))
  ckk=eta**(-2./(cqq-2.))
  calp=cqq/(cqq-2.)
  calp2=calp+1.
!- calculation of f(lin) -!
  flin=abs(1.-ckk*clmin**calp2)**(1./calp2)
!
  if(myrank .eq. 0) then
     write(*,*) 'lin, utin, amdain, flin, eta, ckk=', &
                 clmin, utin, amdain, flin, eta, ckk
  endif
!
!- calculate h -!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!
!- initialize -!
        ah(i,j,k)=1.d0
!
        rr=x1(i)
        th=x3(k)

        sth=sin(th)
        if(abs(sth) .lt. small) then
          if(sth .ge. 0.d0) then
            sth=small
          elseif(sth .lt. 0.d0) then
            sth=-small 
          endif   
       endif
       
        cth=cos(th) 
!
!        if(rr .ge. rin) then         
        if(rr .ge. rin .and. th .ge. thedge .and. th .le. thedge1) then
!- calculation of BL Kerr metric -! 
          del=(rr*rr)-(2.*rr)+(akm*akm)
          aa=(rr*rr+akm*akm)*(rr*rr+akm*akm)-del*akm*akm*sth*sth
          sig=rr*rr+akm*akm*cth*cth
!
          do n=0,3
             do m=0,3
              gcov1(m,n)=0.d0  
              gcon1(m,n)=0.d0
            enddo
          enddo
!
          gcov1(0,0)=-(sig-2.*rr)*(1./sig)               !- g^tt -!
          gcov1(1,1)=sig*(1./del)                          !- g^rr -!
          gcov1(2,2)=(aa*(1./sig))*sth*sth             !- g^phi phi -!
          gcov1(0,2)=-(2.*akm*rr*sth*sth)*(1./sig)           !- g^t phi -!
          gcov1(2,0)=gcov1(0,2)                  !- g^phi t -!
          gcov1(3,3)=sig                      !- g^th th -!    
          
          gcon1(0,0)=-aa*(1./(sig*del))               !- g^tt -!
          gcon1(1,1)=del*(1./sig)                          !- g^rr -!
          gcon1(2,2)=(sig-2.*rr)*(1./(del*sig*sth*sth)) !- g^phi phi -!
          gcon1(0,2)=-(2.*akm*rr)*(1./(sig*del))           !- g^t phi -!       
          gcon1(2,0)=gcon1(0,2)                  !- g^phi t -!
          gcon1(3,3)=1./sig 
!
!- calculation of l -!
          if(akm .eq. 0.d0) then
            tmp1b=-gcon1(0,0)*(1./gcon1(2,2)) 
            amda=sqrt(abs(tmp1b))
            clm=eta*amda**(2.-cqq)
            fl=abs(1.-ckk*clm**calp2)**(1./calp2)
          else
            delta=1.d0
            deltaend=1.d-5
            f=1.d0
            df=1.d0
            nnn=0
!- initial value of l for Newton-Raphson method -!
            tmp1b=(clmin*(gcon1(0,0)-clmin*gcon1(0,2)) &
                  *(1./(gcon1(0,2)-clmin*gcon1(2,2))))
            amda=sqrt(abs(tmp1b))
            clm=eta*amda**(2.-cqq)   
            cqq1=0.5*(2.-cqq)          
!
!- Newton Raphson iteration method -! 
  95        continue
!
            if(abs(delta) .gt. deltaend) then
!
!              tmp5a=clm*((gcon1(0,2)-clm*gcon1(2,2))**cqq1)

              tmp1b=(clm*gcon1(0,0)-gcon1(0,2)*clm*clm) &
                    *(1./(gcon1(0,2)-clm*gcon1(2,2))) 
              tmp5a1=sqrt(abs(tmp1b))
              tmp5a2=eta*tmp5a1**(2.-cqq)
              f=tmp5a2-clm
!
              tmp5b1=(gcon1(0,0)-2.*clm*gcon1(0,2))*(gcon1(0,2)-clm*gcon1(2,2))
              tmp5b2=(clm*gcon1(0,0)-clm*clm*gcon1(0,2))*(-gcon1(2,2))
              tmp5b3=(gcon1(0,2)-clm*gcon1(2,2))**2
              tmp5b4=(tmp5b1-tmp5b2)*(1./tmp5b3)
              tmp5b5=tmp5a1**(-cqq)
              df=0.5*(2.-cqq)*eta*tmp5b4*tmp5b5-1.d0
!
              if(f .eq. 0.d0 .or. df .eq. 0.d0) then
                write(6,*) 'Bad numerical divergence in thick torus' 
                write(6,*) 'f, df =', f, df,'at', i, j, k
              else
                delta=f/df
              endif 
! 
              clm=clm-delta
              nnn=nnn+1
              if(nnn .gt. 1000) then
                write(6,*) 'iteration is over 1000 in thick torus'    
              endif

              goto 95
            endif
 
            fl=abs(1.-ckk*clm**calp2)**(1./calp2)          
          endif
!         
          tmp1=abs(gcon1(0,0)-2.*clm*gcon1(0,2)+clm*clm*gcon1(2,2))
          ut=-1./sqrt(tmp1) 

          ah(i,j,k)=(utin*flin)*(1./(ut*fl))
!
!          if(j .eq. 4 .and. i .eq. 60) then
!             write(*,*) 'ah=', ah(i,j,k), 180.*x3(k)/pi
!          endif

        else
          ah(i,j,k)=1.d0
        endif
!
!- set torus structure -!
!
        if(ah(i,j,k) .gt. 1.d0 .and. rr .ge. rin  &
           .and. th .ge. thedge .and. th .le. thedge1) then  
!        if(ah(i,j,k) .gt. 1.d0 .and. rr .ge. rin) then    
!    
          hm1=max((ah(i,j,k) -1.)*(gam-1.)/gam, 0.d0)
          de=(hm1*(1./kpol))**(1./(gam-1.))
          pr=kpol*(hm1*(1./kpol))**(gam/(gam-1.))
!
!          util(1)=0.d0
!          util(3)=0.d0
!
!
!- calculation of omega & contravariant u^t
!          
          omg=-(gcov1(0,2)+gcov1(0,0)*clm)*(1./(gcov1(2,2)+gcov1(0,2)*clm))
          uut=sqrt(-1./(gcov1(0,0)+2.*omg*gcov1(0,2)+omg*omg*gcov1(2,2)))
!          uut=-1./(ut*(1.-omg*clm))
!
          ucon(0)=uut
          ucon(1)=0.d0 
          ucon(2)=(omg*uut)
!          ucon(2)=0.d0
          ucon(3)=0.d0
!
!          ucov(0)=ut
!          ucov(1)=0.d0
!          ucov(2)=-clm*ut
!          ucov(3)=0.d0
!
!- cal spacial 4-velocity (in BL coord.)-!
!
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
!          call upper(ucov,ucon,gcon1)
!
          gfl=ucon(0)*alpha1

          util(1)=ucon(1)+gfl*beta1(1)*(1./alpha1)
          util(2)=ucon(2)+gfl*beta1(2)*(1./alpha1)
          util(3)=ucon(3)+gfl*beta1(3)*(1./alpha1)             
!
!- Put perturbation in pressure -!
!
         call random_number(tmp_rand)
!!          write(*,*) "tmp_rand=", tmp_rand-0.5
          pr=pr*(1.+0.05*(tmp_rand-0.5))
!
!- check density & pressure max -!
!
          if(de .gt. rhomax1) then
            rhomax1=de
          endif
          if(pr .gt. prmax1) then
            prmax1=pr
          endif
!
!- transform velocity from BL to KS coordinates -!
!
!          if(metric .eq. 303 .or. metric .eq. 403) then
!            x1aa=x1(i)
!            x2aa=x2(j)
!            x3aa=x3(k)   
!            x3t=x3b(k)
!            call transbl2ks(util,utiln,x1aa,x2aa,x3aa,x3t)
!            util(1)=utiln(1)
!            util(2)=utiln(2)
!            util(3)=utiln(3)
!
!          elseif(metric .eq. 203 .and. akm .ne. 0.d0) then
!
!            do m=0,3
!              do n=0,3
!                gcov1(m,n)=gcov(m,n,i,j,k)
!                gcon1(m,n)=gcon(m,n,i,j,k)
!              enddo
!            enddo
!
!            alpha1=1./sqrt(-gcon1(0,0))
!
!            beta1(1)=alpha1*alpha1*gcon1(0,1)
!            beta1(2)=alpha1*alpha1*gcon1(0,2)
!            beta1(3)=alpha1*alpha1*gcon1(0,3)
!
!            ucon(1)=util(1)
!            ucon(2)=util(2)
!            ucon(3)=util(3)
!
!            tmp4a=gcov1(0,0)
!            tmp4b=2.*(gcov1(0,1)*ucon(1)+gcov1(0,2)*ucon(2)+gcov1(0,3)*ucon(3))
!            tmp4c=1.+gcov1(1,1)*ucon(1)*ucon(1)+gcov1(2,2)*ucon(2)*ucon(2) &   ! 
!                  +gcov1(3,3)*ucon(3)*ucon(3) &
!                  +2.*(gcov1(1,2)*ucon(1)*ucon(2)+gcov1(1,3)*ucon(1)*ucon(3) & ! 
!                  +gcov1(2,3)*ucon(2)*ucon(3))   
!            discr=tmp4b*tmp4b-4.*tmp4a*tmp4c
!            ucon(0)=(-tmp4b-sqrt(discr))/(2.*tmp4a)
! 
!            gfl=ucon(0)*alpha1
!
!            util(1)=ucon(1)+gfl*beta1(1)/alpha1
!            util(2)=ucon(2)+gfl*beta1(2)/alpha1
!            util(3)=ucon(3)+gfl*beta1(3)/alpha1  
!
!          endif
!
!          
!- set primitive variables -!
!          
          uri(1,i,j,k)=de
          uri(2,i,j,k)=util(1)
          uri(3,i,j,k)=util(2)
          uri(4,i,j,k)=util(3)
          uri(5,i,j,k)=pr
          uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
!          
        endif
!
      enddo
    enddo
  enddo
!
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_allreduce(rhomax1,rhomax1a,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
  call mpi_allreduce(prmax1,prmax1a,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)

  if(myrank .eq. 0) then
    write(*,*) "rhomax1a, prmax1a=", rhomax1a, prmax1a
  endif
!
!
!- Normalize the densities and pressure -!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1 
        rr=x1(i)
        th=x3(k)
        
!        rhomin1=rhomin
!        rhomin1=rhomin*exp(-3.*rr/rcen)
        rhomin1=rhomin*(rr)**(-3./2.)

        prmin1=(gam-1.d0)*prmin*((rr)**(-5./2.))
!        prmin1=kpol*(rhomin1**gam)
!        prmin1=prmin

        if(ah(i,j,k) .gt. 1.d0 .and. rr .ge. rin .and. &
           th .ge. thedge .and. th .le. thedge1) then
!        if(ah(i,j,k) .gt. 1.d0 .and. rr .ge. rin) then
          de=uri(1,i,j,k)*(1./rhomax1a)
          pr=uri(5,i,j,k)*(1./rhomax1a)
!
          if(de .le. rhomin1) then
            de=rhomin1
          endif
          if(pr .le. prmin1) then
            pr=prmin1
          endif
!
          uri(1,i,j,k)=de
          uri(5,i,j,k)=pr
          uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
        endif
!
      enddo
    enddo
  enddo
!
   prmax2a=prmax1a/rhomax1a
   rhomax1a=1.d0
!
  if(myrank .eq. 0) then
    write(*,*) 'rhomax1a, prmax2a=', rhomax1a, prmax2a
  endif
!
!
!- Set vector potential at cell-center (r-theta frame)-!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
        pom(i,j,k)=0.d0
      enddo
    enddo
  enddo
!
!  do k=ks1,ke1
!    do j=js1,je1
!      do i=is1,ie1

  do k=ks1+1,ke1
    do j=js1,je1
      do i=is1+1,ie1


        rr=x1(i)
        th=x3(k)
         
!        if(ah(i,j,k) .gt. 1.d0 .and. rr .ge. rin .and. &
!           th .ge. thedge .and. th .le. thedge1) then

!
!        rr1=x1(i)
!        th1=x3(k)
!        rr1=(x1(i+1)+x1(i))/2.
!        th1=(x3(k+1)+x3(k))/2.
!        pom1=0.5*rr1*sin(th1)

!!        rhoav=0.25*(uri(1,i-1,j,k)+uri(1,i,j,k-1) &
!!                     +uri(1,i+1,j,k)+uri(1,i,j,k+1))
!
! A_phi(i,j,k)=A_phi(i-1/2,j,k-1/2)
!
        rhoav=0.25*(uri(1,i,j,k)+uri(1,i-1,j,k) &
                    +uri(1,i,j,k-1)+uri(1,i-1,j,k-1))
!          rhoav=uri(1,i,j,k)
!
          if(rhoav .ne. 0.d0) then
            pom1=(rhoav*(1./rhomax1a)) - 0.3d0
          else
            pom1=0.d0
          endif
!
          if(pom1 .gt. 0.d0) then
            pom(i,j,k)=pom1
          endif
!
!        endif        
      enddo
    enddo
  enddo
!
!- put boundary data for magnetic field -!
!
  do j=js1,je1
    do i=is1,ie1
      pom(i,j,ks1)=pom(i,j,ks1+1)
!      pom(i,j,ke1)=pom(i,j,ke1-1)
    enddo
  enddo
!
  do k=ks1,ke1
    do j=js1,je1
      pom(is1,j,k)=pom(is1+1,j,k)
!      pom(ie1,j,k)=pom(ie1-1,j,k)
    enddo
  enddo
!
!- Set-up of Magnetic field -!
  do k=ks1+1,ke1-1
    do j=js1,je1
      do i=is1+1,ie1-1
!        
!!        bcon(1)=-(pom(i,j,k)-pom(i,j,k+1)+pom(i+1,j,k)-pom(i+1,j,k+1)) &
!!                   /(2.*(x3(k+1)-x3(k))*detg(i,j,k))
!!
!!        bcon(3)=(pom(i,j,k)+pom(i,j,k+1)-pom(i+1,j,k)-pom(i+1,j,k+1)) &
!!                   /(2.*(x1(i+1)-x1(i))*detg(i,j,k))
!
        tmp1a=0.5*(pom(i+1,j,k+1)+pom(i,j,k+1)-pom(i+1,j,k)-pom(i,j,k))
!        tmp1a=pom(i,j,k+1)-pom(i,j,k) !- GR way ?
!        tmp1a=sin(x3(k+1))*pom(i,j,k+1)-sin(x3(k))*pom(i,j,k) !- Newtonian
!        dx3a=x3(k+1)-x3(k) 
        dx3a=0.5*(x3(k+1)+x3(k))-0.5*(x3(k)+x3(k-1))
        if(tmp1a .ne. 0.d0) then
          bcon(1)=-tmp1a*(1./(dx3a*detg(i,j,k))) !- GR way ?
!          bcon(1)=-tmp1a*(1./(x1(i)*sin(x3(k))*dx3a)) !- Newtonian
!          bcon(1)=-tmp1a*(1./(x1(i)*dx3a)) 
        else
          bcon(1)=0.d0
        endif
!!
        tmp1b=0.5*(pom(i+1,j,k+1)+pom(i+1,j,k)-pom(i,j,k+1)-pom(i,j,k))
!        tmp1b=pom(i+1,j,k)-pom(i,j,k) !- GR way ?
!        tmp1b=x1(i+1)*pom(i+1,j,k)-x1(i)*pom(i,j,k) !- Newtonian
!        dx1a=x1(i+1)-x1(i)
        dx1a=0.5*(x1(i+1)+x1(i))-0.5*(x1(i)+x1(i-1))
        if(tmp1b .ne. 0.d0) then
          bcon(3)=tmp1b*(1./(dx1a*detg(i,j,k))) !- GR way ?
!          bcon(3)=tmp1b*(1./(x1(i)*dx1a)) !- Newtonian
!          bcon(3)=tmp1b*(1./dx1a)            
        else
          bcon(3)=0.d0
        endif
!
!!        if(k .eq. ks1) then
!!          bcon(1)=(pom(i,j,k+1)-pom(i,j,k))/((x3(k+1)-x3(k))*detg(i,j,k))
!!        elseif(k .eq. ke1) then
!!          bcon(1)=(pom(i,j,k)-pom(i,j,k-1))/((x3(k)-x3(k-1))*detg(i,j,k))
!!        else         
!!          bcon(1)=(pom(i,j,k+1)-pom(i,j,k-1))/((x3(k+1)-x3(k-1))*detg(i,j,k))
!!        endif
!
!!        if(i .eq. is1) then
!!          bcon(3)=-(pom(i+1,j,k)-pom(i,j,k))/((x1(i+1)-x1(i))*detg(i,j,k))
!!        elseif(i .eq. ie1) then
!!          bcon(3)=-(pom(i,j,k)-pom(i-1,j,k))/((x1(i)-x1(i-1))*detg(i,j,k))
!!        else         
!!          bcon(3)=-(pom(i+1,j,k)-pom(i-1,j,k))/((x1(i+1)-x1(i-1))*detg(i,j,k))
! !       endif
!
        uri(7,i,j,k)=bcon(1)
        uri(9,i,j,k)=bcon(3)
!
      enddo
    enddo
  enddo
!
!- put boundary data for magnetic field -!
!
  do j=js1,je1
    do i=is1,ie1
      uri(7,i,j,ks1)=uri(7,i,j,ks1+1)
      uri(7,i,j,ke1)=uri(7,i,j,ke1-1)
!
      uri(9,i,j,ks1)=uri(9,i,j,ks1+1)
      uri(9,i,j,ke1)=uri(9,i,j,ke1-1)
    enddo
  enddo
!
  do k=ks1,ke1
    do j=js1,je1
      uri(7,is1,j,k)=uri(7,is1+1,j,k)
      uri(7,ie1,j,k)=uri(7,ie1-1,j,k)
!
      uri(9,is1,j,k)=uri(9,is1+1,j,k)
      uri(9,ie1,j,k)=uri(9,ie1-1,j,k)
    enddo
  enddo
!
!- calculation maximum magnetic field strength -!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
        bcon(1)=uri(7,i,j,k)
        bcon(2)=uri(8,i,j,k)
        bcon(3)=uri(9,i,j,k)
!       
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
          enddo
        enddo
!
        call calgfl(util,gfl,gcov1)
        call calucon(util,ucon,gfl,gcon1)
        call lower(ucon,ucov,gcov1)
        call calbbcon(bcon,bbcon,ucov,ucon)
        call lower(bbcon,bbcov,gcov1)
        call calbbsq(bbcon,bbcov,bbsq)
!
        if(bbsq .gt. bsq_max) then
          bsq_max=bbsq
        endif
!
      enddo
    enddo
  enddo    
!
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_allreduce(bsq_max,bsq_max1,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
!
  if(bsq_max1 .eq. 0.d0) then
    beta_act=0.d0 
  else   
    beta_act=prmax2a*(1./(0.5*bsq_max1))
  endif   
!
  if(beta .eq. 0.d0) then
    bnorm=0.d0
  else
    bnorm=sqrt(beta_act*(1./beta))
  endif
!
  if(myrank .eq. 0) then
    write(*,*) "beta_act, bnorm=", beta_act, bnorm
  endif
!
!- Normalize magnetic field strength -!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!!
        bcon(1)=uri(7,i,j,k)*bnorm            
        bcon(2)=uri(8,i,j,k)*bnorm  
        bcon(3)=uri(9,i,j,k)*bnorm  
!
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
!
        uri(7,i,j,k)=bcon(1)
        uri(8,i,j,k)=bcon(2)
        uri(9,i,j,k)=bcon(3)
!
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
          enddo
        enddo        
!
        call calgfl(util,gfl,gcov1)
        call calucon(util,ucon,gfl,gcon1)
        call lower(ucon,ucov,gcov1)
        call calbbcon(bcon,bbcon,ucov,ucon)
        call lower(bbcon,bbcov,gcov1)
        call calbbsq(bbcon,bbcov,bbsq)
!
        if(bbsq .gt. bsq_max2) then
          bsq_max2=bbsq
        endif
!
      enddo
    enddo
  enddo 
!
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_allreduce(bsq_max2,bsq_max2a,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
  if(bsq_max2a .eq. 0.d0) then
    beta_act=0.d0 
  else   
    beta_act=prmax2a*(1./(0.5*bsq_max2a))
   endif  
!
  if(myrank .eq. 0) then
    write(*,*) "beta_act=", beta_act
  endif    
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
        rr=x1(i)
        th=x3(k)*180.d0/pi 
!        if(x1(i) .ge. 35.d0 .and. x1(i) .le. 35.5d0 .and. j .eq. 4) then
!          if((uri(7,i,j,k)) .gt. 0.d0 .or. abs(uri(9,i,j,k)) .gt. 0.d0) then
!          write(*,*) 'rr, th, br, bph, bth=', rr, th, uri(7,i,j,k), uri(8,i,j,k), uri(9,i,j,k)
!          endif           
!        endif    

      enddo
    enddo  
  enddo
    
  deallocate(gcov1, gcon1, gcon1in, gcov1in, util, utiln, bcon, ucov, ucon, &
             bbcov, bbcon, beta1, ah, pom, stat=merr)
!
  return
end subroutine mdtorus2
!
!--------------------------------------------------------------------
subroutine mdtorus3(uri,gcov,gcon,detg, &
                    x1,x2,x3,x1b,x2b,x3b,dx1,dx2,dx3, &
                    myrank,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
!  constant-l rotating torus model
!  Font & Daigne 2002, MNRAS, 334, 383
!  model=9  
!
  use pram, only : imax, jmax, kmax, nv, gam, metric, akm, pi, iter, &
                   kpol, ieos, dmin, pmin
  implicit none
  include 'mpif.h'

  integer :: i, j, k, m, n, is1 ,ie1, js1, je1, ks1, ke1
  integer :: merr, myrank, jh, nnn

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  
  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:), gcov1in(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:), beta1(:)
  real(8), allocatable :: util(:), utiln(:), bcon(:)
!
  real(8) :: x1(imax), x2(jmax), x3(kmax) 
  real(8) :: x1b(imax), x2b(jmax), x3b(kmax) 
  real(8) :: dx1(imax), dx2(jmax), dx3(kmax)
!
  real(8),allocatable :: ah(:,:,:), pom(:,:,:)
!
  real(8) :: de, pr  
  real(8) :: kappa, beta, clm, utin, ut, uut, &
             rhomin, prmin, rhomin1, prmin1, hm1, omg, thedge, thedge1, rcen, &
             rhomax1, rhomax1a, prmax1, prmax1a, prmax2, prmax2a 
  real(8) :: rr, th, sth, cth, del, aa, sig, rr1, th1, &
             rin, thin, sthin, cthin, delin, aain, sigin, rin1
  real(8) :: x1aa, x2aa, x3aa, x3t
  real(8) :: tmp_rand
  real(8) :: bsq_max, bsq_max1, bsq_max2, bsq_max2a, beta_act, &
             bnorm, pom1, rhoav, gfl, bbsq, &
             pbeta, pbeta_min, pbeta_min1, pbeta_min2, pbeta_min2a      
  real(8) :: tmp4a, tmp4b, tmp4c, discr, alpha1
  real(8) :: tmp1a, tmp1b, tmp2a, tmp2b, dx1a, dx3a
  real(8) :: ucovph
  real(8) :: small
! 
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), gcov1in(0:3,0:3), stat=merr)
  allocate(util(1:3), utiln(1:3), bcon(1:3), ucov(0:3), ucon(0:3), &
           bbcov(0:3), bbcon(0:3), beta1(1:3), stat=merr)
  allocate(ah(is1:ie1,js1:je1,ks1:ke1), pom(is1:ie1,js1:je1,ks1:ke1), &
           stat=merr)

!--------------------------------------------------------------------
!- Parameter -!

  rin=9.34d0
!  rin=9.0d0
  thin=pi/2.d0
  kappa=1.0d-3
!  beta=1.d2
!  beta=100.d0
  beta=0.d0
  clm=4.5d0 ! a=0.0
!  clm=4.35d0 ! a=0.5
!  clm=4.25d0 ! a=0.9 (0.998)
!  
  thedge=0.1d0*pi
  thedge1=0.9d0*pi
  rcen=25.d0
!  
  rin1=1.d0
  rhomin=dmin
!  rhomin=1.d-5 
!  prmin=0.5*rhomin
!  if(ieos .eq. 3) then
!    prmin=kpol*(rhomin**gam)
!  else
    prmin=pmin
!    prmin=kpol*(rhomin**gam)
!  endif
!
  rhomax1=0.d0
  rhomax1a=0.d0
  prmax1=0.d0 
  prmax1a=0.d0
  prmax2=0.d0
!  
  bsq_max=0.d0
  bsq_max1=0.d0
  bsq_max2=0.d0
  bsq_max2a=0.d0
!
  pbeta_min=1.d5
  pbeta_min1=1.d5
  pbeta_min2=1.d5
  pbeta_min2a=1.d5
 !
  small=1.d-12  
!
!- initally set minimum plasma everywhere -!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1

        do m=0,3
          do n=0,3
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
        
        rr=x1(i)

!        de=rhomin
!        de=rhomin*exp(-3.*rr/rcen)
        de=rhomin*rr**(-3./2.)

        pr=(gam-1.d0)*prmin*rr**(-5./2.)
!        pr=kpol*(de**gam)
!        pr=(gam-1.)*prmin
        
!        pr=(gam-1.)*prmin*rr**(-5./2.)
        
        util(1)=0.d0
!        util(1)=(gcon1(0,1)/gcon1(0,0))*(1.-(1./rr)**4)
!        util(1)=beta1(1)*(1./alpha1)        
        util(2)=0.d0
        util(3)=0.d0
        bcon(1)=0.d0
        bcon(2)=0.d0 
        bcon(3)=0.d0
!          
!- set primitive variables -!
!          
        uri(1,i,j,k)=de
        uri(2,i,j,k)=util(1)
        uri(3,i,j,k)=util(2)
        uri(4,i,j,k)=util(3)
        uri(5,i,j,k)=pr
        uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
        uri(7,i,j,k)=bcon(1)
        uri(8,i,j,k)=bcon(2)          
        uri(9,i,j,k)=bcon(3)
!       
      enddo
    enddo
  enddo
!  
!- calculation of BL Kerr metric at rin, thin -!
!  
!  sthin=sin(thin)
  sthin=1.d0
!  cthin=cos(thin)
  cthin=0.d0
  delin=rin*rin-2.*rin+akm*akm
  aain=(rin*rin+akm*akm)*(rin*rin+akm*akm)-delin*akm*akm*sthin*sthin
  sigin=rin*rin+akm*akm*cthin*cthin  
!
  do n=0,3
    do m=0,3
      gcov1in(m,n)=0.d0
    enddo
  enddo

  gcov1in(0,0)=-(sigin-2.*rin)*(1./sigin)               !- g_tt -!
  gcov1in(2,2)=(aain*(1./sigin))*sthin*sthin !- g_phi phi -!
  gcov1in(0,2)=-(2.*akm*rin*sthin*sthin)*(1./sigin)  !- g_t phi -!
! 
!- calculation of ut at rin -!
  tmp1a=gcov1in(0,2)*gcov1in(0,2)-gcov1in(0,0)*gcov1in(2,2)
  tmp1b=gcov1in(0,0)*clm*clm+2.*clm*gcov1in(0,2)+gcov1in(2,2)
  if(tmp1b .gt. 0.d0 .and. tmp1a .gt. 0.d0) then
    utin=-sqrt(tmp1a*(1./tmp1b))
  else
    utin=-1.d0  
  endif
!
  if(myrank .eq. 0) then
    write(*,*) 'lin, utin=', clm, utin
  endif
!
!- calculate h -!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!          
!- initialize -!
!         
        ah(i,j,k)=1.d0
!         
        rr=x1(i)
        th=x3(k)
        
        sth=sin(th)
        if(abs(sth) .lt. small) then
          if(sth .ge. 0.d0) then
            sth=small
          elseif(sth .lt. 0.d0) then
            sth=-small 
          endif   
        endif

        cth=cos(th) 
!
!        if(rr .ge. rin) then         
        if(rr .ge. rin .and. th .ge. thedge .and. th .le. thedge1) then
!- calculation of BL Kerr metric -! 
          del=(rr*rr)-(2.*rr)+(akm*akm)
          aa=(rr*rr+akm*akm)*(rr*rr+akm*akm)-del*akm*akm*sth*sth
          sig=rr*rr+akm*akm*cth*cth
!
          do n=0,3
            do m=0,3
              gcov1(m,n)=0.d0
              gcon1(m,n)=0.d0 
            enddo
          enddo
!
          gcov1(0,0)=-(sig-2.*rr)*(1./sig)               !- g^tt -!
          gcov1(1,1)=sig*(1./del)                          !- g^rr -!
          gcov1(2,2)=(aa*(1./sig))*sth*sth             !- g^phi phi -!
          gcov1(0,2)=-(2.*akm*rr*sth*sth)*(1./sig)           !- g^t phi -!
          gcov1(2,0)=gcov1(0,2)                  !- g^phi t -!
          gcov1(3,3)=sig                      !- g^th th -!        
!
          gcon1(0,0)=-aa*(1./(sig*del))               !- g^tt -!
          gcon1(1,1)=del*(1./sig)                          !- g^rr -!
          gcon1(2,2)=(sig-2.*rr)*(1./(del*sig*sth*sth)) !- g^phi phi -!
          gcon1(0,2)=-(2.*akm*rr)*(1./(sig*del))           !- g^t phi -!
          gcon1(2,0)=gcon1(0,2)                  !- g^phi t -!
          gcon1(3,3)=1./sig 
!          
!- calculation of u_t at r,theta -!
!          
          tmp2a=gcov1(0,2)*gcov1(0,2)-gcov1(0,0)*gcov1(2,2)
          tmp2b=gcov1(0,0)*clm*clm+2.*clm*gcov1(0,2)+gcov1(2,2)
          if(tmp2b .gt. 0.d0 .and. tmp2a .gt. 0.d0) then
            ut=-sqrt(tmp2a*(1./tmp2b))
          else
            ut=-1.d0  
          endif        

!
!- calculation of h = u_t,in/u_t
!          
          ah(i,j,k)=utin*(1./ut)
!
        else
          ah(i,j,k)=1.d0
        endif
!
!- set torus structure -!
!
        if(ah(i,j,k) .gt. 1.d0 .and. rr .ge. rin  &
           .and. th .ge. thedge .and. th .le. thedge1) then    
!         if((ah(i,j,k)-1.d0) .gt. 1.d-7 .and. rr .ge. rin) then
!    
          hm1=max((ah(i,j,k) -1.)*(gam-1.)/gam, 0.d0)
          de=(hm1*(1./kpol))**(1./(gam-1.))
          pr=kpol*(hm1*(1./kpol))**(gam/(gam-1.))
!
          ucon(1)=0.d0
          ucon(3)=0.d0
!
!- calculation of omega & contravariant u^t
!          
          omg=-(gcov1(0,2)+gcov1(0,0)*clm)*(1./(gcov1(2,2)+gcov1(0,2)*clm))
          uut=sqrt(-1./(gcov1(0,0)+2.*omg*gcov1(0,2)+omg*omg*gcov1(2,2)))
!          uut=-1./(ut*(1.-omg*clm))
!
          ucovph=-omg*ut

          ucon(0)=uut
          ucon(2)=(omg*uut)
!          ucon(2)=gcon1(0,2)*ut+gcon1(2,2)*ucovph
!          ucon(2)=0.d0
!
!- cal spacial 4-velocity -!
!
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
          gfl=ucon(0)*alpha1

          util(1)=ucon(1)+gfl*beta1(1)*(1./alpha1)
          util(2)=ucon(2)+gfl*beta1(2)*(1./alpha1)
          util(3)=ucon(3)+gfl*beta1(3)*(1./alpha1)             
!
!- Put perturbation in pressure -!
!
         call random_number(tmp_rand)
!!          write(*,*) "tmp_rand=", tmp_rand-0.5
!           pr=pr*(1.+0.05*(tmp_rand-0.5))
!           util(1)=util(1)+0.01*(tmp_rand-0.5)             
!
!- check density & pressure max -!
!
          if(de .gt. rhomax1) then
            rhomax1=de
!            write(*,*) 'de,i,k,rr,th=', de,i,k,x1(i),x3(k)
          endif
          if(pr .gt. prmax1) then
            prmax1=pr
          endif
!
!- transform velocity from BL to KS coordinates (do not need?)-!
!
!          if(metric .eq. 303 .or. metric .eq. 403) then
!            x1aa=x1(i)
!            x2aa=x2(j)
!            x3aa=x3(k)   
!            x3t=x3b(k)
!            call transbl2ks(util,utiln,x1aa,x2aa,x3aa,x3t)
!            util(1)=utiln(1)
!            util(2)=utiln(2)
!            util(3)=utiln(3)
!
!          elseif(metric .eq. 203 .and. akm .ne. 0.d0) then
!
!            do m=0,3
!              do n=0,3
!                gcov1(m,n)=gcov(m,n,i,j,k)
!                gcon1(m,n)=gcon(m,n,i,j,k)
!              enddo
!            enddo
!
!            alpha1=1./sqrt(-gcon1(0,0))
!
!            beta1(1)=alpha1*alpha1*gcon1(0,1)
!            beta1(2)=alpha1*alpha1*gcon1(0,2)
!            beta1(3)=alpha1*alpha1*gcon1(0,3)
!
!            ucon(1)=util(1)
!            ucon(2)=util(2)
!            ucon(3)=util(3)
!
!            tmp4a=gcov1(0,0)
!            tmp4b=2.*(gcov1(0,1)*ucon(1)+gcov1(0,2)*ucon(2)+gcov1(0,3)*ucon(3))
!            tmp4c=1.+gcov1(1,1)*ucon(1)*ucon(1)+gcov1(2,2)*ucon(2)*ucon(2) &
!                  +gcov1(3,3)*ucon(3)*ucon(3) &
!                  +2.*(gcov1(1,2)*ucon(1)*ucon(2)+gcov1(1,3)*ucon(1)*ucon(3) &
!                  +gcov1(2,3)*ucon(2)*ucon(3))   
!            discr=tmp4b*tmp4b-4.*tmp4a*tmp4c
!            ucon(0)=(-tmp4b-sqrt(discr))/(2.*tmp4a)
! 
!            gfl=ucon(0)*alpha1
!
!            util(1)=ucon(1)+gfl*beta1(1)/alpha1
!            util(2)=ucon(2)+gfl*beta1(2)/alpha1
!            util(3)=ucon(3)+gfl*beta1(3)/alpha1  
!
!          endif
!
!          
!- set primitive variables -!
!          
          uri(1,i,j,k)=de
          uri(2,i,j,k)=util(1)
          uri(3,i,j,k)=util(2)
          uri(4,i,j,k)=util(3)
          uri(5,i,j,k)=pr
          uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
!          
        endif
!
      enddo
    enddo
  enddo
!
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_allreduce(rhomax1,rhomax1a,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
  call mpi_allreduce(prmax1,prmax1a,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
!
  if(myrank .eq. 0) then
    write(*,*) "rhomax1a, prmax1a=", rhomax1a, prmax1a
  endif
!
!
!- Normalize the densities and pressure -!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1 
        rr=x1(i)
        th=x3(k)
        
!        rhomin1=rhomin
!        rhomin1=rhomin*exp(-3.*rr/rcen)
        rhomin1=rhomin*rr**(-3./2.)

!        prmin1=(gam-1.d0)*prmin*((rr/rin1)**(-2.5d0))
!        prmin1=kpol*(de**gam)
!        prmin1=(gam-1.)*prmin
        prmin1=(gam-1.)*prmin*rr**(-5./2.)

        if(ah(i,j,k) .gt. 1.d0 .and. rr .ge. rin &
           .and. th .ge. thedge .and. th .le. thedge1) then
!        if((ah(i,j,k)-1.d0) .gt. 1.d-7 .and. rr .ge. rin) then
!          de=uri(1,i,j,k)/rhomax1a
!          pr=uri(5,i,j,k)/rhomax1a
!
          de=uri(1,i,j,k)/1.d0
          pr=uri(5,i,j,k)/1.d0          
!
          if(de .le. rhomin1) then
            de=rhomin1
          endif
          if(pr .le. prmin1) then
            pr=prmin1
          endif
!
          uri(1,i,j,k)=de
          uri(5,i,j,k)=pr
          uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
        endif
!
      enddo
    enddo
  enddo
!

!   prmax2a=prmax1a/rhomax1a
   prmax2a=prmax1a
!   rhomax1a=1.d0
!
  if(myrank .eq. 0) then
    write(*,*) 'rhomax1a, prmax2a=', rhomax1a, prmax2a
  endif
!
!
!- Set vector potential at cell-center (r-theta frame)-!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
        pom(i,j,k)=0.d0
      enddo
    enddo
  enddo
!
!  do k=ks1,ke1
!    do j=js1,je1
!      do i=is1,ie1
!  
  do k=ks1+1,ke1
    do j=js1,je1
      do i=is1+1,ie1

!  do k=ks1,ke1
!    do j=js1,je1
!      do i=is1,ie1

!!
        rr1=x1(i)
        th1=x3(k)
!!        rr1=(x1(i+1)+x1(i))/2.
!!        th1=(x3(k+1)+x3(k))/2.
!!        pom1=0.5*rr1*sin(th1)

!        rhoav=0.25*(uri(1,i-1,j,k)+uri(1,i,j,k-1) &
!                     +uri(1,i+1,j,k)+uri(1,i,j,k+1))
!
! A_phi(i,j,k)=A_phi(i-1/2,j,k-1/2)
!
        rhoav=0.25*(uri(1,i,j,k)+uri(1,i-1,j,k) &
                    +uri(1,i,j,k-1)+uri(1,i-1,j,k-1))
!        rhoav=uri(1,i,j,k)
!
        if(rhoav .ne. 0.d0) then
          pom1=(rhoav*(1./rhomax1a)) - 0.2d0
        else
          pom1=0.d0
        endif
!
!
        if(pom1 .gt. 0.d0) then
          pom(i,j,k)=pom1
        endif
!
      enddo
    enddo
  enddo
!
!
!- put boundary data for magnetic field -!
!
  do j=js1,je1
    do i=is1,ie1
      pom(i,j,ks1)=pom(i,j,ks1+1)
!      pom(i,j,ke1)=pom(i,j,ke1-1)
    enddo
  enddo
!
  do k=ks1,ke1
    do j=js1,je1
      pom(is1,j,k)=pom(is1+1,j,k)
!      pom(ie1,j,k)=pom(ie1-1,j,k)
    enddo
  enddo
!
!- Set-up of Magnetic field -!
  do k=ks1,ke1-1
    do j=js1,je1
      do i=is1,ie1-1
!        
!!        bcon(1)=-(pom(i,j,k)-pom(i,j,k+1)+pom(i+1,j,k)-pom(i+1,j,k+1)) &
!!                   /(2.*(x3(k+1)-x3(k))*detg(i,j,k))
!!
!!        bcon(3)=(pom(i,j,k)+pom(i,j,k+1)-pom(i+1,j,k)-pom(i+1,j,k+1)) &
!!                   /(2.*(x1(i+1)-x1(i))*detg(i,j,k))
!
        tmp1a=0.5*(pom(i+1,j,k+1)+pom(i,j,k+1)-pom(i+1,j,k)-pom(i,j,k))
!        tmp1a=pom(i,j,k+1)-pom(i,j,k) !- GR way ?
!        tmp1a=sin(x3(k+1))*pom(i,j,k+1)-sin(x3(k))*pom(i,j,k) !- Newtonian
!        dx3a=x3(k+1)-x3(k) 
        dx3a=0.5*(x3(k+1)+x3(k))-0.5*(x3(k)+x3(k-1))
        if(tmp1a .ne. 0.d0) then
          bcon(1)=-tmp1a*(1./(dx3a*detg(i,j,k))) !- GR way ?
!          bcon(1)=-tmp1a*(1./(x1(i)*sin(x3(k))*dx3a)) !- Newtonian
!          bcon(1)=-tmp1a*(1./(x1(i)*dx3a)) 
        else
          bcon(1)=0.d0
        endif
!!
        tmp1b=0.5*(pom(i+1,j,k+1)+pom(i+1,j,k)-pom(i,j,k+1)-pom(i,j,k))
!        tmp1b=pom(i+1,j,k)-pom(i,j,k) !- GR way ?
!        tmp1b=x1(i+1)*pom(i+1,j,k)-x1(i)*pom(i,j,k) !- Newtonian
        dx1a=0.5*(x1(i+1)+x1(i))-0.5*(x1(i)+x1(i-1))
!        dx1a=x1(i+1)-x1(i)
        if(tmp1b .ne. 0.d0) then
          bcon(3)=tmp1b*(1./(dx1a*detg(i,j,k))) !- GR way ?
!          bcon(3)=tmp1b*(1./(x1(i)*dx1a)) !- Newtonian
!          bcon(3)=tmp1b*(1./dx1a)  
        else
          bcon(3)=0.d0
        endif
!
!!        if(k .eq. ks1) then
!!          bcon(1)=(pom(i,j,k+1)-pom(i,j,k))/((x3(k+1)-x3(k))*detg(i,j,k))
!!        elseif(k .eq. ke1) then
!!          bcon(1)=(pom(i,j,k)-pom(i,j,k-1))/((x3(k)-x3(k-1))*detg(i,j,k))
!!        else         
!!          bcon(1)=(pom(i,j,k+1)-pom(i,j,k-1))/((x3(k+1)-x3(k-1))*detg(i,j,k))
!!        endif
!
!!        if(i .eq. is1) then
!!          bcon(3)=-(pom(i+1,j,k)-pom(i,j,k))/((x1(i+1)-x1(i))*detg(i,j,k))
!!        elseif(i .eq. ie1) then
!!          bcon(3)=-(pom(i,j,k)-pom(i-1,j,k))/((x1(i)-x1(i-1))*detg(i,j,k))
!!        else         
!!          bcon(3)=-(pom(i+1,j,k)-pom(i-1,j,k))/((x1(i+1)-x1(i-1))*detg(i,j,k))
! !       endif
!
        uri(7,i,j,k)=bcon(1)
        uri(9,i,j,k)=bcon(3)
!
      enddo
    enddo
  enddo
!
!- put boundary data for magnetic field -!
!
  do j=js1,je1
    do i=is1,ie1
      uri(7,i,j,ks1)=uri(7,i,j,ks1+1)
      uri(7,i,j,ke1)=uri(7,i,j,ke1-1)
!
      uri(9,i,j,ks1)=uri(9,i,j,ks1+1)
      uri(9,i,j,ke1)=uri(9,i,j,ke1-1)
    enddo
  enddo
!
  do k=ks1,ke1
    do j=js1,je1
      uri(7,is1,j,k)=uri(7,is1+1,j,k)
      uri(7,ie1,j,k)=uri(7,ie1-1,j,k)
!
      uri(9,is1,j,k)=uri(9,is1+1,j,k)
      uri(9,ie1,j,k)=uri(9,ie1-1,j,k)
    enddo
  enddo
!
!- calculation maximum magnetic field strength -!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
        bcon(1)=uri(7,i,j,k)
        bcon(2)=uri(8,i,j,k)
        bcon(3)=uri(9,i,j,k)
!       
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
            gcon1(m,n)=gcon(m,n,i,j,k)             
          enddo
        enddo
!
        call calgfl(util,gfl,gcov1)
        call calucon(util,ucon,gfl,gcon1)
        call lower(ucon,ucov,gcov1)
        call calbbcon(bcon,bbcon,ucov,ucon)
        call lower(bbcon,bbcov,gcov1)
        call calbbsq(bbcon,bbcov,bbsq)
!
        if(bbsq .ne. 0.d0) then        
          pbeta=2.*pr/bbsq
        else
          pbeta=1.d5  
        endif

        if(pbeta .lt. pbeta_min) then
          pbeta_min=pbeta           
        endif
        if(bbsq .gt. bsq_max) then
          bsq_max=bbsq
        endif
!
      enddo
    enddo
  enddo    
!
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_allreduce(bsq_max,bsq_max1,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
  call mpi_allreduce(pbeta_min,pbeta_min1,1,mpi_double_precision,mpi_min, &
                     mpi_comm_world,merr)
!
!  if(bsq_max1 .eq. 0.d0) then
!    beta_act=0.d0 
!  else   
!    beta_act=prmax2a*(1./(0.5*bsq_max1))
!  endif   
  beta_act=pbeta_min1 
 !
  if(beta .eq. 0.d0) then
    bnorm=0.d0
  else
    bnorm=sqrt(beta_act*(1./beta))
  endif
!
  if(myrank .eq. 0) then
    write(*,*) "beta_act, bnorm=", beta_act, bnorm
  endif
!
!- Normalize magnetic field strength -!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!!
        bcon(1)=uri(7,i,j,k)*bnorm            
        bcon(2)=uri(8,i,j,k)*bnorm  
        bcon(3)=uri(9,i,j,k)*bnorm  
!
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
!
        uri(7,i,j,k)=bcon(1)
        uri(8,i,j,k)=bcon(2)
        uri(9,i,j,k)=bcon(3)
!
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo        
!
        call calgfl(util,gfl,gcov1)
        call calucon(util,ucon,gfl,gcon1)
        call lower(ucon,ucov,gcov1)
        call calbbcon(bcon,bbcon,ucov,ucon)
        call lower(bbcon,bbcov,gcov1)
        call calbbsq(bbcon,bbcov,bbsq)
!
        if(bbsq .ne. 0.d0) then        
          pbeta=2.*pr/bbsq
        else
          pbeta=1.d5  
        endif

        if(pbeta .lt. pbeta_min2) then
          pbeta_min2=pbeta           
        endif
       
        if(bbsq .gt. bsq_max2) then
          bsq_max2=bbsq
        endif
!
      enddo
    enddo
  enddo 
!
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_allreduce(bsq_max2,bsq_max2a,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
  call mpi_allreduce(pbeta_min2,pbeta_min2a,1,mpi_double_precision,mpi_min, &
                     mpi_comm_world,merr)  

!  if(bsq_max2a .eq. 0.d0) then
!    beta_act=0.d0 
!  else   
!    beta_act=prmax2a*(1./(0.5*bsq_max2a))
!  endif  
  beta_act=pbeta_min2a
   !
  if(myrank .eq. 0) then
    write(*,*) "beta_act=", beta_act
 endif    
!
  deallocate(gcov1, gcon1, gcov1in, util, utiln, bcon, ucov, ucon, &
             bbcov, bbcon, beta1, ah, pom, stat=merr)
!
  return
end subroutine mdtorus3
!
!--------------------------------------------------------------------
subroutine mdtorus4(uri,gcov,gcon, &
                    x1,x2,x3,x1b,x2b,x3b,dx1,dx2,dx3, &
                    is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
!  constant-l rotating torus model with toroidal magnetic field 
!  Komissarov 2006, MNRAS, 368, 993
!  model=12  
!
  use pram, only : imax, jmax, kmax, nv, gam, metric, akm, pi, iter, kpol, &
                   dmin, pmin
  implicit none
  include 'mpif.h'

  integer :: i, j, k, m, n, is1 ,ie1, js1, je1, ks1, ke1
  integer :: merr

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  
  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:), gcov1cen(:,:), gcov1in(:,:)
  real(8), allocatable :: util(:), utiln(:), bcon(:), bbcon(:)
  real(8), allocatable :: ucov(:), ucon(:), beta1(:)
!
  real(8) :: x1(imax), x2(jmax), x3(kmax) 
  real(8) :: x1b(imax), x2b(jmax), x3b(kmax) 
  real(8) :: dx1(imax), dx2(jmax), dx3(kmax)
!
  real(8) :: de, pr, pm, ut, wp, ah, hm1
  real(8) :: rin, thin, utin, win,win1,  &
             rcen, thcen, ahcen, utcen, wcen, wcen1, clm, betacen, &
             rhocen, thedge, thedge1, rin1, rhomin, prmin, polkp, polkm
  real(8) :: rr, th, sth, cth, del, aa, sig, &
             sthin, cthin, delin, aain, sigin, & 
             sthcen, cthcen, delcen, aacen, sigcen 
  real(8) :: bb, aa1, bbin, aa1in, bbcen, aa1cen, omg, uut, ucovph
  real(8) :: tmp1a, tmp1b, tmp1c, tmp1d, tmp_rand
  real(8) :: x1aa, x2aa, x3aa, x3t
  real(8) :: tmp4a, tmp4b, tmp4c, discr, alpha1, gfl
  real(8) :: rmb, rms, zz1, zz2
  real(8) :: rhomax1, rhomax1a, prmax1, prmax1a, prmax2
  real(8) :: small
! 
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), gcov1in(0:3,0:3), &
           gcov1cen(0:3,0:3), stat=merr)
  allocate(util(1:3), utiln(1:3), bcon(1:3), bbcon(0:3), &
           ucov(0:3), ucon(0:3), beta1(1:3), stat=merr)

!--------------------------------------------------------------------
!- Parameter -!
!- we assume polytropic index for gas pressure and mag pressure is same -!  

  zz1=1.+((1.-akm*akm)**(1./3.))*(((1.+akm)**(1./3.))+((1.-akm)**(1./3.)))
  zz2=sqrt((3.*akm*akm+zz1*zz1))
  rmb=2.*(1.-0.5*akm+sqrt(1.-akm))
  rms=3.+zz2-sqrt((3.-zz1)*(3.+zz1+2.*zz2))

  write(*,*) "rmb, rms=", rmb, rms
  
!  rin=rms     !- radial position at inner radius of torus
  rin=9.34d0
!  rin=3.0d0
  thin=pi/2.d0  !- \theta position at center of torus  
  win1=-0.030d0    !- W at surface of torus
!  rcen=5.0d0    !- radial position at center of torus
  rcen=25.d0
  thcen=pi/2.d0  !- \theta position at center of torus
  ahcen=1.5d0     !- specific enthalpy at center of torus
  wcen1=-0.103d0   !- W at center of torus (calculated from rcen)
  rhocen=1.d0     !- density at center of torus
!  clm=2.8d0      !- angular momentum
!  clm=4.25d0
  clm=4.35d0
!  clm=4.5d0
!  betacen=1.0d0  !- plasma beta at center of torus (\beta=p_gas/p_mag)
  betacen=1.d0
  !
  thedge=0.1*pi
  thedge1=0.9*pi
!  
  rin1=1.d0
  rhomin=dmin    !- minimum density at atmosphere (corona) 
  prmin=pmin !-minimum pressure at atmosphere (corona)
!
  rhomax1=0.d0
  rhomax1a=0.d0
  prmax1=0.d0 
  prmax1a=0.d0
  prmax2=0.d0
!  
  small=1.d-12
!
!- initally set minimum plasma everywhere -!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1

        do m=0,3
          do n=0,3
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo
        alpha1=1./sqrt(-gcon1(0,0))
!
        beta1(1)=alpha1*alpha1*gcon1(0,1)
        beta1(2)=alpha1*alpha1*gcon1(0,2)
        beta1(3)=alpha1*alpha1*gcon1(0,3)
        
        rr=x1(i)

        de=rhomin
!        de=rhomin*exp(-3.*rr/rcen)
!        de=rhomin*((rr/rin1)**(-1.5d0))

!        pr=(gam-1.d0)*prmin*((rr/rin1)**(-2.5d0))
!        pr=kpol*(de**gam)
!
        pr=prmin
!        util(1)=0.d0
!        util(1)=(gcon1(0,1)/gcon1(0,0))*(1.-(1./rr)**4)
        util(1)=beta1(1)*(1./alpha1)
        util(2)=0.d0
        util(3)=0.d0
        bcon(1)=0.d0
        bcon(2)=0.d0 
        bcon(3)=0.d0
!          
!- set primitive variables -!
!          
        uri(1,i,j,k)=de
        uri(2,i,j,k)=util(1)
        uri(3,i,j,k)=util(2)
        uri(4,i,j,k)=util(3)
        uri(5,i,j,k)=pr
        uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
        uri(7,i,j,k)=bcon(1)
        uri(8,i,j,k)=bcon(2)          
        uri(9,i,j,k)=bcon(3)
!       
      enddo
    enddo
  enddo
!  
!- calculation of BL Kerr metric at rin, thin -!
! 
!  sthin=sin(thin)
  sthin=1.d0
!  cthin=cos(thin)
  cthin=0.d0
  delin=rin*rin-2.*rin+akm*akm
  aain=(rin*rin+akm*akm)*(rin*rin+akm*akm)-delin*akm*akm*sthin*sthin
  sigin=rin*rin+akm*akm*cthin*cthin  
!
  do n=0,3
    do m=0,3
      gcov1in(m,n)=0.d0
    enddo
  enddo
!
!-Kerr metric in BL & KS coordinates -!
  gcov1in(0,0)=-(sigin-2.*rin)*(1./sigin)    !- g_tt -!!
  gcov1in(2,2)=(aain*(1./sigin))*sthin*sthin !- g_phi phi -!
  gcov1in(0,2)=-(2.*akm*rin*sthin*sthin)*(1./sigin)  !- g_t phi -!
!          
!- calculation of u_t at r_in,theta_in -!
!          
  bbin=gcov1in(0,2)*gcov1in(0,2)-gcov1in(0,0)*gcov1in(2,2)
  aa1in=gcov1in(0,0)*clm*clm+2.*clm*gcov1in(0,2)+gcov1in(2,2)
          if(aa1in .gt. 0.d0) then
            utin=-sqrt(bbin*(1./aa1in))
!            utin=bbin*(1./aa1in)
          else
            utin=-1.d0  
          endif        
!
!- calculation of W from u_t at disc center
!          
          win=log(abs(utin))
!          win=0.5*log(abs(utin))
!          win=win1
          write(*,*) 'win=', win           
!  
!- calculation of BL Kerr metric at rcen, thcen -!
! 
!  sthcen=sin(thcen)
  sthcen=1.d0
!  cthcen=cos(thcen)
  cthcen=0.d0
  delcen=rcen*rcen-2.*rcen+akm*akm
  aacen=(rcen*rcen+akm*akm)*(rcen*rcen+akm*akm)-delcen*akm*akm*sthcen*sthcen
  sigcen=rcen*rcen+akm*akm*cthcen*cthcen  
!
  do n=0,3
    do m=0,3
      gcov1cen(m,n)=0.d0
    enddo
 enddo
 
!-Kerr metric in BL & KS coordinates -!
!
  gcov1cen(0,0)=-(sigcen-2.*rcen)*(1./sigcen)               !- g_tt -!
  gcov1cen(2,2)=(aacen*(1./sigcen))*sthcen*sthcen !- g_phi phi -!
  gcov1cen(0,2)=-(2.*akm*rcen*sthcen*sthcen)*(1./sigcen)  !- g_t phi -!
!          
  bbcen=gcov1cen(0,2)*gcov1cen(0,2)-gcov1cen(0,0)*gcov1cen(2,2)
  aa1cen=gcov1cen(0,0)*clm*clm+2.*clm*gcov1cen(0,2)+gcov1cen(2,2)
          if(aa1cen .gt. 0.d0) then
            utcen=-sqrt(bbcen/aa1cen)
!            utcen=-(bbcen*(1./aa1cen))
          else
            utcen=-1.d0  
          endif        
!
!- calculation of W from u_t at disc center
!          
          wcen=log(abs(utcen))
!          wcen=0.5*log(abs(utcen))
!          wcen=wcen1
!
          write(*,*) 'wcen=', wcen          
!
!- calculate plokp and polkm
!- (polytropic constant for gas and mag pressure)
!
            
!  tmp1c=exp((win-wcen)/(1.+1./betacen))-1.        
!  tmp1d=rhocen**(gam-1.)
!  polkp=((gam-1.)/gam)*tmp1c/tmp1d
!        
  if(betacen .eq. 0.d0) then
    polkp=1.d-3
    polkm=0.d0     
  else   
    tmp1c=(betacen*(1./(1.+betacen)))        
    polkp=(win-wcen)*((gam-1.)/gam)*tmp1c*(ahcen**(1.-gam))
    polkm=(bbcen**(1.-gam))*polkp*(1./betacen)    
  endif  

  write(*,*) 'polkp, polkm=', polkp, polkm
  
!
!- calculate w at each grid point -!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!          
!- initialize -!
!         
        wp=1.d0
!          
        rr=x1(i)
        th=x3(k)
!        
        sth=sin(th)
        if(abs(sth) .lt. small) then
          if(sth .ge. 0.d0) then
            sth=small
          elseif(sth .lt. 0.d0) then
            sth=-small 
          endif   
        endif
!       
        cth=cos(th) 
!
!        if(rr .ge. rin) then         
        if(rr .ge. rin .and. th .ge. thedge .and. th .le. thedge1) then
!- calculation of BL Kerr metric -! 
          del=(rr*rr)-(2.*rr)+(akm*akm)
          aa=(rr*rr+akm*akm)*(rr*rr+akm*akm)-del*akm*akm*sth*sth
          sig=rr*rr+akm*akm*cth*cth
!
          do n=0,3
            do m=0,3
              gcov1(m,n)=0.d0
              gcon1(m,n)=0.d0 
            enddo
          enddo
!         
!- Kerr metric in BL coordinates -!
          gcov1(0,0)=-(sig-2.*rr)*(1./sig)               !- g^tt -!
          gcov1(1,1)=sig*(1./del)                          !- g^rr -!
          gcov1(2,2)=(aa*(1./sig))*sth*sth             !- g^phi phi -!
          gcov1(0,2)=-(2.*akm*rr*sth*sth)*(1./sig)           !- g^t phi -!
          gcov1(2,0)=gcov1(0,2)                  !- g^phi t -!
          gcov1(3,3)=sig                      !- g^th th -!        
!
          gcon1(0,0)=-aa*(1./(sig*del))               !- g^tt -!
          gcon1(1,1)=del*(1./sig)                          !- g^rr -!
          gcon1(2,2)=(sig-2.*rr)*(1./(del*sig*sth*sth)) !- g^phi phi -!
          gcon1(0,2)=-(2.*akm*rr)*(1./(sig*del))           !- g^t phi -!
          gcon1(2,0)=gcon1(0,2)                  !- g^phi t -!
          gcon1(3,3)=1./sig 
!
!- Kerr metric in KS coordinates -!
!          gcov1(0,0)=-1.+2.*rr/sig      !- g_tt -!
!          gcov1(0,1)=2.*rr/sig           !- g_tr -!
!          gcov1(0,2)=-2.*akm*rr*sth*sth/sig !- g_t phi -!
!          gcov1(1,0)=gcov1(0,1)        !- g_rt -!
!          gcov1(1,1)=1.+2.*rr/sig      !- g_rr -!
!          gcov1(1,2)=-1.*akm*sth*sth*(1.+2.*rr/sig) !- g_r phi -! 
!          gcov1(2,0)=gcov1(0,2)        !- g_phi t -!
!          gcov1(2,1)=gcov1(1,2)        !- g_phi r -!
!          gcov1(2,2)=aa*sth*sth/sig         !- g_phi phi -!
!          gcov1(3,3)=sig                !- g_th th -!
!
!          gcon1(0,0)=-1.-2.*rr/sig      !- g^tt -!
!          gcon1(0,1)=2.*rr/sig           !- g^tr -!
!          gcon1(1,0)=gcon1(0,1)        !- g^rt -!
!          gcon1(1,1)=del/sig           !- g^rr -!
!          gcon1(1,2)=akm/sig            !- g^r phi -!
!          gcon1(2,1)=gcon1(1,2)        !- g^phi r -!
!          gcon1(2,2)=1./(sig*sth*sth)         !- g^phi phi -!
!          gcon1(3,3)=1./sig                !- g^th th -!
!          
!- calculation of u_t at r,theta -!
!          
          bb=gcov1(0,2)*gcov1(0,2)-gcov1(0,0)*gcov1(2,2)
          aa1=gcov1(0,0)*clm*clm+2.*clm*gcov1(0,2)+gcov1(2,2)
          if(aa1 .gt. 0.d0) then
            ut=-sqrt(bb*(1./aa1))
!            ut=-(bb*(1./aa1))
          else
            ut=-1.d0  
          endif        
!
!- calculation of W from u_t
!          
          wp=log(abs(ut))
!          wp=0.5*log(abs(ut))
!
!- calculation of h
!
          if(wp .le. win) then
!- from Vincent et al.             
!            tmp1a=(wp-win)/(wcen-win)
!            tmp1b=(polkp+polkm*bbcen**(gam-1.))/(polkp+polkm*bb**(gam-1.))
!            ah=ahcen*(tmp1a*tmp1b)**(1./(gam-1.))          
!
!- from Komissarov 
            tmp1a=(win-wp)
!            tmp1b=1.+(polkm/polkp)*bb**(gam-1.)
            tmp1b=(gam/(gam-1.))*(polkp+polkm*bb**(gam-1.))
!            tmp1b=1.d0
            if(tmp1a .gt. 0.d0) then
!              ah=exp(tmp1a/tmp1b)            
              ah=(tmp1a*(1./tmp1b))**(1./(gam-1.))
            else
              ah=-1.d0
            endif
!
!- for ideal EoS used?
!            tmp1a=exp(win-wp)-1.
!            tmp1b=(gam/(gam-1.))*(polkp+polkm*bb**(gam-1.))
!            if(tmp1a .gt. 0.d0) then
!!              ah=exp(tmp1a/tmp1b)            
!              ah=(tmp1a*(1./tmp1b))**(1./(gam-1.))
!            else
!              ah=-1.d0
!            endif
!
          else
            ah=-1.d0
          endif            
!          
        else
          ah=-1.d0
        endif
!
!- set torus structure -!
!
!        if(ah .gt. 1.d-4 .and. rr .ge. rin) then  &
!           .and. th .ge. thedge .and. th .le. thedge1) then    
        if(ah .gt. 1.d-4) then
!
!          hm1=max((ah-1.)*(gam-1.)/gam, 0.d0) 
!          de=(hm1/polkp)**(1./(gam-1.)) 
!          pr=polkp*(hm1/polkp)**(gam/(gam-1.))
!          pm=polkm*(bb**(gam-1.))*(pr/polkp)
!
!- Komissarov
          pr=polkp*(ah**gam)
          pm=polkm*(bb**(gam-1.))*(ah**gam)
!          de=ah-pr
          de=(pr/polkp)**(1./gam)
!          de=ah-((gam)/(gam-1.))*pr
          
!          de=ah
!          pr=polkp*(ah**gam)
!          pm=polkm*(bb**(gam-1.))*(1.+(gam/(gam-1.))*de**(gam-1.))**gam

          if(de .lt. rhomin) then
            de=rhomin
          endif
          if(pr .lt. prmin) then
            pr=prmin
          endif   
          
          ucon(1)=0.d0
          ucon(3)=0.d0
!
!- calculation of omega & contravariant u^t
!          
          omg=-(gcov1(0,2)+gcov1(0,0)*clm)*(1./(gcov1(2,2)+gcov1(0,2)*clm))
          uut=sqrt(-1./(gcov1(0,0)+2.*omg*gcov1(0,2)+omg*omg*gcov1(2,2)))
!          uut=-1./(ut*(1.-omg*clm))
!
          ucovph=-omg*ut

          ucon(0)=uut
          ucon(2)=(omg*uut)
!          ucon(2)=gcon1(0,2)*ut+gcon1(2,2)*ucovph
!          ucon(2)=0.d0
!
!
          if(betacen .eq. 0.d0) then
            bbcon(1)=0.d0
            bbcon(2)=0.d0
            bbcon(3)=0.d0
            bbcon(0)=0.d0
          else   
            bbcon(1)=0.d0
            if(aa1 .gt. 0.d0) then 
              bbcon(2)=sqrt(2.*pm*(1./aa1))
            else
              bbcon(2)=0.d0
            endif   
            bbcon(3)=0.d0            
            bbcon(0)=clm*bbcon(2)
!            
!            bcon(1)=0.d0
!            bcon(2)=sqrt(2.*pm)
!            bcon(3)=0.d0
          endif
!
!- cal 3-magnetic field -!
!          
          bcon(1)=ucon(0)*bbcon(1)-bbcon(0)*ucon(1)
          bcon(2)=ucon(0)*bbcon(2)-bbcon(0)*ucon(2)
          bcon(3)=ucon(0)*bbcon(3)-bbcon(0)*ucon(3)
!
!          bcon(1)=gfl*bbcon(1)-bbcon(0)*ucon(1)
!          bcon(2)=gfl*bbcon(2)-bbcon(0)*ucon(2)
!          bcon(3)=gfl*bbcon(3)-bbcon(0)*ucon(3)
!
!- cal spacial 4-velocity -!
!
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
          gfl=ucon(0)*alpha1

          util(1)=ucon(1)+gfl*beta1(1)*(1./alpha1)
          util(2)=ucon(2)+gfl*beta1(2)*(1./alpha1)
          util(3)=ucon(3)+gfl*beta1(3)*(1./alpha1)             
!
!- Put perturbation in pressure -!
!
!         call random_number(tmp_rand)
!!          write(*,*) "tmp_rand=", tmp_rand-0.5
!           pr=pr*(1.+0.05*(tmp_rand-0.5))
!           util(1)=util(1)+0.01*(tmp_rand-0.5)             
!
!- calculation of bcon -!
!
!          do m=0,3
!            do n=0,3
!              gcov1(m,n)=gcov(m,n,i,j,k)
!              gcon1(m,n)=gcon(m,n,i,j,k)
!            enddo
!          enddo
!
!          alpha1=1./sqrt(-gcon1(0,0))
!                
!          call calgfl(util,gfl,gcov1)
!          call calucon(util,ucon,gfl,gcon1)
!
!          bcon(1)=ucon(0)*bbcon(1)-bbcon(0)*ucon(1)
!          bcon(2)=ucon(0)*bbcon(2)-bbcon(0)*ucon(2)
!          bcon(3)=ucon(0)*bbcon(3)-bbcon(0)*ucon(3)
!          
!          bcon(1)=alpha1*ucon(0)*bbcon(1)-alpha1*bbcon(0)*ucon(1)
!          bcon(2)=alpha1*ucon(0)*bbcon(2)-alpha1*bbcon(0)*ucon(2)
!          bcon(3)=alpha1*ucon(0)*bbcon(3)-alpha1*bbcon(0)*ucon(3)
!          
!          if(j .eq. 4 .and. k .eq. 50) then
!            write(*,*) 'rr, B\phi=', x1(i), bcon(2)
!          endif
!          
!- check density & pressure max -!
!
          if(de .gt. rhomax1) then
            rhomax1=de
!            write(*,*) 'de,i,k,rr,th=', de,i,k,x1(i),x3(k)
          endif
          if(pr .gt. prmax1) then
            prmax1=pr
          endif
!         
!- set primitive variables -!
!          
          uri(1,i,j,k)=de
          uri(2,i,j,k)=util(1)
          uri(3,i,j,k)=util(2)
          uri(4,i,j,k)=util(3)
          uri(5,i,j,k)=pr
          uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
          uri(7,i,j,k)=bcon(1)
          uri(8,i,j,k)=bcon(2)
          uri(9,i,j,k)=bcon(3)
!          
        endif
!
      enddo
    enddo
  enddo
!
  deallocate(gcov1, gcon1, gcov1cen, gcov1in, util, utiln, bcon, bbcon, &
             ucov, ucon, beta1, stat=merr)
!
  return
end subroutine mdtorus4
!
!--------------------------------------------------------------------
subroutine mdtildisk1(uri,gcov,gcon,detg, &
                      x1,x2,x3,x1b,x2b,x3b,dx1,dx2,dx3, &
                      myrank,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
!  tilted constant-l rotating torus model
!  Based on Font & Daigne 2002, MNRAS, 334, 383
!  and Fragile et al. 2007, ApJ, 668, 417
!  model=14  
!
  use pram, only : imax, jmax, kmax, nv, gam, metric, akm, pi, iter, kpol
  implicit none
  include 'mpif.h'

  integer :: i, j, k, m, n, is1 ,ie1, js1, je1, ks1, ke1
  integer :: merr, myrank, jh, nnn

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  
  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: gcov1(:,:), gcov1bl(:,:), gcon1(:,:), gcon1bl(:,:), &
                          gcon1ks(:,:), gcov1in(:,:)
  real(8), allocatable :: ucov(:), ucon(:), ucon1(:), &
                          bbcov(:), bbcon(:), beta1ks(:)
  real(8), allocatable :: util(:), utiln(:), bcon(:)
  real(8), allocatable :: trans(:,:), tmp1(:)
!
  real(8) :: x1(imax), x2(jmax), x3(kmax) 
  real(8) :: x1b(imax), x2b(jmax), x3b(kmax) 
  real(8) :: dx1(imax), dx2(jmax), dx3(kmax)
!
  real(8),allocatable :: ah(:,:,:), pom(:,:,:)
!
  real(8) :: de, pr  
  real(8) :: kappa, beta, clm, utin, ut, uut, &
             rhomin, prmin, hm1, omg, rcen, &
             rhomax1, rhomax1a, prmax1, prmax1a, prmax2, prmax2a 
  real(8) :: rr, theta, phi, th, sth1, cth1, sph1, cph1, &
             del, aa, sig, rr1, th1, &
             rin, thin, sthin, cthin, delin, aain, sigin, rin1
  real(8) :: x1aa, x2aa, x3aa, x3t
  real(8) :: tmp_rand
  real(8) :: bsq_max, bsq_max1, bsq_max2, bsq_max2a, beta_act, &
             bnorm, pom1, rhoav, gfl, bbsq
  real(8) :: tmp4a, tmp4b, tmp4c, discr, alpha1ks
  real(8) :: tmp1a, tmp1b, tmp2a, tmp2b
  real(8) :: ucovph
! 
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcov1bl(0:3,0:3), gcon1(0:3,0:3), gcon1bl(0:3,0:3), &
           gcon1ks(0:3,0:3), gcov1in(0:3,0:3), stat=merr)
  allocate(util(1:3), utiln(1:3), bcon(1:3), ucov(0:3), ucon(0:3), &
           ucon1(0:3), &
           bbcov(0:3), bbcon(0:3), beta1ks(1:3), trans(0:3,0:3), tmp1(0:3), &
           stat=merr)
  allocate(ah(is1:ie1,js1:je1,ks1:ke1), pom(is1:ie1,js1:je1,ks1:ke1), &
           stat=merr)

!--------------------------------------------------------------------
!- Parameter -!

  rin=9.34d0
!  rin=9.0d0
!  rin=4.d0
  thin=pi/2.d0
  kappa=1.0d-3
  beta=0.d0
!  beta=100.d0
!  beta=0.d0
!  clm=3.385d0
!  clm=4.5d0 ! a=0.0
  clm=4.35d0 ! a=0.5
!  clm=4.25d0 ! a=0.9
!
  rcen=25.d0
!  
  rin1=1.d0
  rhomin=1.d-5
  prmin=1.d-3*rhomin
!
  rhomax1=0.d0
  rhomax1a=0.d0
  prmax1=0.d0 
  prmax1a=0.d0
  prmax2=0.d0
!  
  bsq_max=0.d0
  bsq_max1=0.d0
  bsq_max2=0.d0
  bsq_max2a=0.d0
!
!- initally set minimum plasma everywhere -!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1

        do m=0,3
          do n=0,3
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo

        rr=x1(i)

        de=rhomin
!        de=1.d-3*exp(-3.*rr/rcen)
!        de=rhomin*((rr/rin1)**(-1.5d0))

!        pr=(gam-1.d0)*prmin*((rr/rin1)**(-2.5d0))
!        pr=kpol*(de**gam)
!
        pr=prmin
        util(1)=0.d0
!        util(1)=(gcon1(0,1)/gcon1(0,0))*(1.-(1./rr)**4)
        util(2)=0.d0
        util(3)=0.d0
        bcon(1)=0.d0
        bcon(2)=0.d0 
        bcon(3)=0.d0
!          
!- set primitive variables -!
!          
        uri(1,i,j,k)=de
        uri(2,i,j,k)=util(1)
        uri(3,i,j,k)=util(2)
        uri(4,i,j,k)=util(3)
        uri(5,i,j,k)=pr
        uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
        uri(7,i,j,k)=bcon(1)
        uri(8,i,j,k)=bcon(2)          
        uri(9,i,j,k)=bcon(3)
!       
      enddo
    enddo
  enddo
!  
!- calculation of BL Kerr metric at rin, thin -!
!  
  sthin=sin(thin)
  cthin=cos(thin)
  delin=rin*rin-2.*rin+akm*akm
  aain=(rin*rin+akm*akm)*(rin*rin+akm*akm)-delin*akm*akm*sthin*sthin
  sigin=rin*rin+akm*akm*cthin*cthin  
!
  do n=0,3
    do m=0,3
      gcov1in(m,n)=0.d0
    enddo
  enddo

  gcov1in(0,0)=-(sigin-2.*rin)/sigin               !- g_tt -!
  gcov1in(1,1)=sigin/delin                          !- g_rr -!
  gcov1in(2,2)=(aain/sigin)*sthin*sthin !- g_phi phi -!
  gcov1in(0,2)=-(2.*akm*rin*sthin*sthin)/sigin  !- g_t phi -!
  gcov1in(2,0)=gcov1in(0,2)                  !- g_phi t -!
  gcov1in(3,3)=sigin                      !- g_th th -!
! 
!- calculation of ut at rin -!
  tmp1a=gcov1in(0,2)*gcov1in(0,2)-gcov1in(0,0)*gcov1in(2,2)
  tmp1b=gcov1in(0,0)*clm*clm+2.*clm*gcov1in(0,2)+gcov1in(2,2)
  if(tmp1b .gt. 0.d0) then
    utin=-sqrt(tmp1a/tmp1b)
  else
    utin=-1.d0  
  endif
!
  if(myrank .eq. 0) then
    write(*,*) 'lin, win=', clm, log(abs(utin))
  endif
!
!- calculate h -!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!
!- set tilted position         
        rr=x1(i)
!         
        if(x2(j) .le. 0.d0) then
          phi=2.*pi-x2(j)
        elseif(x2(j) .gt. 2.*pi) then
          phi=x2(j)-2.*pi
        else
          phi=x2(j)
        endif
!       
        if(x3(k) .le. 0.d0) then
          theta=2.*pi-x3(k)
        elseif(x3(k) .gt. 2.*pi) then
          theta=x3(k)-2.*pi
        else
          theta=x3(k)
        endif
!
        call calrotpos(phi,theta,sph1,cph1,sth1,cth1)
!
        if(rr .ge. rin) then         
!- calculation of BL Kerr metric -! 
          del=(rr*rr)-(2.*rr)+(akm*akm)
          aa=(rr*rr+akm*akm)*(rr*rr+akm*akm)-del*akm*akm*sth1*sth1
          sig=rr*rr+akm*akm*cth1*cth1
!
          do n=0,3
            do m=0,3
              gcov1bl(m,n)=0.d0
              gcon1bl(m,n)=0.d0 
            enddo
          enddo
!
          gcov1bl(0,0)=-(sig-2.*rr)/sig               !- g^tt -!
          gcov1bl(1,1)=sig/del                          !- g^rr -!
          gcov1bl(2,2)=(aa/sig)*sth1*sth1             !- g^phi phi -!
          gcov1bl(0,2)=-(2.*akm*rr*sth1*sth1)/sig           !- g^t phi -!
          gcov1bl(2,0)=gcov1(0,2)                  !- g^phi t -!
          gcov1bl(3,3)=sig                      !- g^th th -!        
!
          gcon1bl(0,0)=-aa/(sig*del)               !- g^tt -!
          gcon1bl(1,1)=del/sig                          !- g^rr -!
          gcon1bl(2,2)=(sig-2.*rr)/(del*sig*sth1*sth1) !- g^phi phi -!
          gcon1bl(0,2)=-(2.*akm*rr)/(sig*del)           !- g^t phi -!
          gcon1bl(2,0)=gcon1bl(0,2)                  !- g^phi t -!
          gcon1bl(3,3)=1./sig 
!          
!- calculation of u_t at r,theta -!
!          
          tmp2a=gcov1bl(0,2)*gcov1bl(0,2)-gcov1bl(0,0)*gcov1bl(2,2)
          tmp2b=gcov1bl(0,0)*clm*clm+2.*clm*gcov1bl(0,2)+gcov1bl(2,2)
          if(tmp2b .gt. 0.d0) then
            ut=-sqrt(tmp2a/tmp2b)
          else
            ut=-1.d0  
          endif        

!
!- calculation of h = u_t,in/u_t
!          
          ah(i,j,k)=utin/ut
!
        else
          ah(i,j,k)=-1.d0
        endif
!
!- set torus structure -!
!
        if(ah(i,j,k) .gt. 1.d0 .and. rr .ge. rin) then    
!    
          hm1=max((ah(i,j,k) -1.)*(gam-1.)/gam, 0.d0)
          de=(hm1/kpol)**(1./(gam-1.))
          pr=kpol*(hm1/kpol)**(gam/(gam-1.))
!
!         if(j .eq. 4 .and. k .eq. 32) then
!           write(*,*) 'ah, de=', ah(i,j,k), de(i,j,k), i
!         endif
!
          util(1)=0.d0
          util(3)=0.d0
!
!- calculation of omega & contravariant u^t
!          
          omg=-(gcov1bl(0,2)+gcov1bl(0,0)*clm)/(gcov1bl(2,2)+gcov1bl(0,2)*clm)
          uut=sqrt(-1./(gcov1bl(0,0)+2.*omg*gcov1bl(0,2)+omg*omg*gcov1bl(2,2)))
!          uut=-1./(ut*(1.-omg*clm))
!
          ucovph=-omg*ut

          util(2)=(omg*uut)
!          util(2)=gcon1bl(0,2)*ut+gcon1bl(2,2)*ucovph
!          util(2)=0.d0
!
!- Put perturbation in pressure -!
!
         call random_number(tmp_rand)
!!          write(*,*) "tmp_rand=", tmp_rand-0.5
!           pr=pr*(1.+0.05*(tmp_rand-0.5))
!           util(1)=util(1)+0.01*(tmp_rand-0.5)             
!
!- check density & pressure max -!
!
          if(de .gt. rhomax1) then
            rhomax1=de
!            write(*,*) 'de,i,k,rr,th=', de,i,k,x1(i),x3(k)
          endif
          if(pr .gt. prmax1) then
            prmax1=pr
          endif
!
!- transform velocity from BL to KS coordinates -!
!- set 4-velocity in BL coordinate
          ucon(1)=util(1)
          ucon(2)=util(2)
          ucon(3)=util(3)
!
          tmp4a=gcov1bl(0,0)
          tmp4b=2.*(gcov1bl(0,1)*ucon(1)+gcov1bl(0,2)*ucon(2) &
                +gcov1bl(0,3)*ucon(3))
          tmp4c=1.+gcov1bl(1,1)*ucon(1)*ucon(1) &
                +gcov1bl(2,2)*ucon(2)*ucon(2) &
                +gcov1bl(3,3)*ucon(3)*ucon(3) &
                +2.*(gcov1bl(1,2)*ucon(1)*ucon(2) &
                +gcov1bl(1,3)*ucon(1)*ucon(3) &
                +gcov1bl(2,3)*ucon(2)*ucon(3))   
          discr=tmp4b*tmp4b-4.*tmp4a*tmp4c
          ucon(0)=(-tmp4b-sqrt(discr))/(2.*tmp4a)                    
!
!- make transform matrix -!
          do m=0,3
            do n=0,3
              trans(m,n)=0.d0
            enddo
          enddo
          trans(0,0)=1.d0
          trans(0,1)=2.*rr/(rr*rr -(2.*rr) +akm*akm)
          trans(1,1)=1.d0 
          trans(2,2)=1.d0
          trans(2,1)=akm/(rr*rr -(2.*rr) +akm*akm)
          trans(3,3)=1.d0          
!
!- initialize -!
          do m=0,3
           tmp1(m)=0.d0
           ucon1(m)=0.d0
          enddo
!
          do n=0,3
            tmp1(0)=tmp1(0)+trans(0,n)*ucon(n)
            tmp1(1)=tmp1(1)+trans(1,n)*ucon(n)
            tmp1(2)=tmp1(2)+trans(2,n)*ucon(n)
            tmp1(3)=tmp1(3)+trans(3,n)*ucon(n)
          enddo
!
!- contravariant velocity in KS coordinates -!
          do m=0,3
            ucon1(m)=tmp1(m)
          enddo
!
!- get Kerr-Schild metric (contravariant metric) -!          
          do n=0,3
            do m=0,3
              gcon1ks(m,n)=0.d0
            enddo
          enddo
!
          gcon1ks(0,0)=-1.-2.*rr/sig      !- g^tt -!
          gcon1ks(0,1)=2.*rr/sig           !- g^tr -!
!
          gcon1ks(1,0)=gcon1ks(0,1)        !- g^rt -!
          gcon1ks(1,1)=del/sig           !- g^rr -!
          gcon1ks(1,2)=akm/sig            !- g^r phi -!
! 
          gcon1ks(2,1)=gcon1ks(1,2)        !- g^phi r -!
          gcon1ks(2,2)=1./(sig*sth1*sth1)    !- g^phi phi -!
!
          gcon1ks(3,3)=1./sig                !- g_th th -!          
!
          alpha1ks=1./sqrt(-gcon1ks(0,0))
          gfl=ucon1(0)*alpha1ks
!!
          beta1ks(1)=alpha1ks*alpha1ks*gcon1ks(0,1)
          beta1ks(2)=alpha1ks*alpha1ks*gcon1ks(0,2)
          beta1ks(3)=alpha1ks*alpha1ks*gcon1ks(0,3)
!!          
          util(1)=ucon1(1)+gfl*beta1ks(1)/alpha1ks
          util(2)=ucon1(2)+gfl*beta1ks(2)/alpha1ks
          util(3)=ucon1(3)+gfl*beta1ks(3)/alpha1ks          
!
!          
!- set primitive variables -!
!          
          uri(1,i,j,k)=de
          uri(2,i,j,k)=util(1)
          uri(3,i,j,k)=util(2)
          uri(4,i,j,k)=util(3)
          uri(5,i,j,k)=pr
          uri(6,i,j,k)=pr*(1./(de**(gam-1.))) 
        endif
!
      enddo
    enddo
  enddo
!
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_allreduce(rhomax1,rhomax1a,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
  call mpi_allreduce(prmax1,prmax1a,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
!
  if(myrank .eq. 0) then
    write(*,*) "rhomax1a, prmax1a=", rhomax1a, prmax1a
  endif
!
!
!- Normalize the densities and pressure -!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1 
        rr=x1(i)
        if(ah(i,j,k) .gt. 1.d0 .and. rr .ge. rin ) then
          de=uri(1,i,j,k)/rhomax1a
          pr=uri(5,i,j,k)/rhomax1a
!
          if(de .le. rhomin) then
            de=rhomin
          endif
          if(pr .le. prmin) then
            pr=prmin
          endif
!
          uri(1,i,j,k)=de
          uri(5,i,j,k)=pr
          uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
        endif
!
      enddo
    enddo
  enddo
!
   prmax2a=prmax1a/rhomax1a
   rhomax1a=1.d0
!
  if(myrank .eq. 0) then
    write(*,*) 'rhomax1a, prmax2a=', rhomax1a, prmax2a
  endif
!
!- Set vector potential at cell-center (r-theta frame)-!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
        pom(i,j,k)=0.d0
      enddo
    enddo
  enddo
!
!  do k=ks1+1,ke1-1
!    do j=js1,je1
!      do i=is1+1,ie1-1
!
!!        rr1=x1(i)
!!        th1=x3(k)
!!        rr1=(x1(i+1)+x1(i))/2.
!!        th1=(x3(k+1)+x3(k))/2.
!!        pom1=0.5*rr1*sin(th1)

!!!        rhoav=0.25*(uri(1,i-1,j,k)+uri(1,i,j,k-1) &
!!!                     +uri(1,i+1,j,k)+uri(1,i,j,k+1))
!         rhoav=0.25*(uri(1,i,j,k)+uri(1,i-1,j,k) &
!                     +uri(1,i,j,k-1)+uri(1,i-1,j,k-1))
!        rhoav=uri(1,i,j,k)
!
!        if(rhoav .ne. 0.d0) then
!          pom1=rhoav/rhomax1a - 0.3
!        else
!          pom1=0.d0
!        endif
!
!         if(j .eq. 4 .and. k .eq. 64) then
!           write(*,*) 'pom1=', pom1, i, j, k
!         endif
!
!        if(pom1 .gt. 0.d0) then
!          pom(i,j,k)=pom1
!        endif
!
!      enddo
!    enddo
!  enddo
!
!- Set-up of Magnetic field -!
!  do k=ks1+1,ke1-1
!    do j=js1,je1
!      do i=is1+1,ie1-1
!        
!!        bcon(1)=-(pom(i,j,k)-pom(i,j,k+1)+pom(i+1,j,k)-pom(i+1,j,k+1)) &
!!                   /(2.*(x3(k+1)-x3(k))*detg(i,j,k))
!!
!!        bcon(3)=(pom(i,j,k)+pom(i,j,k+1)-pom(i+1,j,k)-pom(i+1,j,k+1)) &
!!                   /(2.*(x1(i+1)-x1(i))*detg(i,j,k))
!
!        tmp1a=pom(i,j,k)-pom(i,j,k+1)+pom(i+1,j,k)-pom(i+1,j,k+1)
!        if(tmp1a .ne. 0.d0) then
!          bcon(1)=-tmp1a/(2.*dx3(k)*detg(i,j,k))
!        else
!          bcon(1)=0.d0
!        endif
!!
!        tmp1b=pom(i,j,k)+pom(i,j,k+1)-pom(i+1,j,k)-pom(i+1,j,k+1)
!        if(tmp1b .ne. 0.d0) then
!          bcon(3)=tmp1b/(2.*dx1(i)*detg(i,j,k))
!        else
!          bcon(3)=0.d0
!        endif
!
!!        if(k .eq. ks1) then
!!          bcon(1)=(pom(i,j,k+1)-pom(i,j,k))/((x3(k+1)-x3(k))*detg(i,j,k))
!!        elseif(k .eq. ke1) then
!!          bcon(1)=(pom(i,j,k)-pom(i,j,k-1))/((x3(k)-x3(k-1))*detg(i,j,k))
!!        else         
!!          bcon(1)=(pom(i,j,k+1)-pom(i,j,k-1))/((x3(k+1)-x3(k-1))*detg(i,j,k))
!!        endif
!
!!        if(i .eq. is1) then
!!          bcon(3)=-(pom(i+1,j,k)-pom(i,j,k))/((x1(i+1)-x1(i))*detg(i,j,k))
!!        elseif(i .eq. ie1) then
!!          bcon(3)=-(pom(i,j,k)-pom(i-1,j,k))/((x1(i)-x1(i-1))*detg(i,j,k))
!!        else         
!!          bcon(3)=-(pom(i+1,j,k)-pom(i-1,j,k))/((x1(i+1)-x1(i-1))*detg(i,j,k))
!!        endif
!!
!        uri(7,i,j,k)=bcon(1)
!        uri(9,i,j,k)=bcon(3)
!
!      enddo
!    enddo
!  enddo
!
!- put boundary data for magnetic field -!
!
!  do j=js1,je1
!    do i=is1,ie1
!      uri(7,i,j,ks1)=uri(7,i,j,ks1+1)
!      uri(7,i,j,ke1)=uri(7,i,j,ke1-1)
!
!      uri(9,i,j,ks1)=uri(9,i,j,ks1+1)
!      uri(9,i,j,ke1)=uri(9,i,j,ke1-1)
!    enddo
!  enddo
!
!  do k=ks1,ke1
!    do j=js1,je1
!      uri(7,is1,j,k)=uri(7,is1+1,j,k)
!      uri(7,ie1,j,k)=uri(7,ie1-1,j,k)
!
!      uri(9,is1,j,k)=uri(9,is1+1,j,k)
!      uri(9,ie1,j,k)=uri(9,ie1-1,j,k)
!    enddo
!  enddo
!
!- calculation maximum magnetic field strength -!
!  do k=ks1,ke1
!    do j=js1,je1
!      do i=is1,ie1
!
!        util(1)=uri(2,i,j,k)
!        util(2)=uri(3,i,j,k)
!        util(3)=uri(4,i,j,k)
!        bcon(1)=uri(7,i,j,k)
!        bcon(2)=uri(8,i,j,k)
!        bcon(3)=uri(9,i,j,k)
!       
!        do m=0,3
!          do n=0,3
!            gcov1(m,n)=gcov(m,n,i,j,k)
!          enddo
!        enddo
!
!        call calgfl(util,gfl,gcov1)
!        call calucon(util,ucon,gfl,gcon1)
!        call lower(ucon,ucov,gcov1)
!        call calbbcon(bcon,bbcon,ucov,ucon)
!        call lower(bbcon,bbcov,gcov1)
!        call calbbsq(bbcon,bbcov,bbsq)
!
!        if(bbsq .gt. bsq_max) then
!          bsq_max=bbsq
!        endif
!
!      enddo
!    enddo
!  enddo    
!
!  call mpi_barrier(mpi_comm_world,merr)
!  call mpi_allreduce(bsq_max,bsq_max1,1,mpi_double_precision,mpi_max, &
!                     mpi_comm_world,merr)
!
!  beta_act=prmax2a/(0.5*bsq_max1)
!
!  if(beta .eq. 0.d0) then
!    bnorm=0.d0
!  else
!    bnorm=sqrt(beta_act/beta)
!  endif
!
!  if(myrank .eq. 0) then
!    write(*,*) "beta_act, bnorm=", beta_act, bnorm
!  endif
!
!- Normalize magnetic field strength -!
!
!  do k=ks1,ke1
!    do j=js1,je1
!      do i=is1,ie1
!
!        bcon(1)=uri(7,i,j,k)*bnorm            
!        bcon(2)=uri(8,i,j,k)*bnorm  
!        bcon(3)=uri(9,i,j,k)*bnorm  
!
!        util(1)=uri(2,i,j,k)
!        util(2)=uri(3,i,j,k)
!        util(3)=uri(4,i,j,k)
!
!        uri(7,i,j,k)=bcon(1)
!        uri(8,i,j,k)=bcon(2)
!        uri(9,i,j,k)=bcon(3)
!
!        do m=0,3
!          do n=0,3
!            gcov1(m,n)=gcov(m,n,i,j,k)
!          enddo
!        enddo        
!
!        call calgfl(util,gfl,gcov1)
!        call calucon(util,ucon,gfl,gcon1)
!        call lower(ucon,ucov,gcov1)
!        call calbbcon(bcon,bbcon,ucov,ucon)
!        call lower(bbcon,bbcov,gcov1)
!        call calbbsq(bbcon,bbcov,bbsq)
!
!        if(bbsq .gt. bsq_max2) then
!          bsq_max2=bbsq
!        endif
!
!      enddo
!    enddo
!  enddo 
!
!  call mpi_barrier(mpi_comm_world,merr)
!  call mpi_allreduce(bsq_max2,bsq_max2a,1,mpi_double_precision,mpi_max, &
!                     mpi_comm_world,merr)
!  beta_act=prmax2a/(0.5*bsq_max2a)
!
!  if(myrank .eq. 0) then
!    write(*,*) "beta_act=", beta_act
!  endif    
!
  deallocate(gcov1, gcov1bl, gcon1, gcon1bl, gcon1ks, gcov1in, &
             util, utiln, bcon, ucov, ucon, ucon1, &
             bbcov, bbcon, beta1ks, trans, tmp1, ah, pom, stat=merr)
!
  return
end subroutine mdtildisk1
!
!--------------------------------------------------------------------
subroutine mdtildisk2(uri,gcov,gcon,detg, &
                      x1,x2,x3,x1b,x2b,x3b,dx1,dx2,dx3, &
                      myrank,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
!  constant-l rotating tilted torus model
!  Font & Daigne 2002, MNRAS, 334, 383
!
  use pram, only : imax, jmax, kmax, nv, gam, metric, akm, pi, iter, kpol, &
                    tilang
  implicit none
  include 'mpif.h'

  integer :: i, j, k, m, n, is1 ,ie1, js1, je1, ks1, ke1
  integer :: merr, myrank, jh, nnn

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  
  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:), gcov1in(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:), beta1(:)
  real(8), allocatable :: util(:), utiln(:), bcon(:)
!
  real(8) :: x1(imax), x2(jmax), x3(kmax) 
  real(8) :: x1b(imax), x2b(jmax), x3b(kmax) 
  real(8) :: dx1(imax), dx2(jmax), dx3(kmax)
!
  real(8),allocatable :: ah(:,:,:), pom(:,:,:)
!
  real(8) :: de, pr  
  real(8) :: kappa, beta, clm, utin, ut, uut, &
             rhomin, prmin, hm1, omg, thedge, thedge1, rcen, &
             rhomax1, rhomax1a, prmax1, prmax1a, prmax2, prmax2a 
  real(8) :: rr, ph, th, th1, sph, cph, sth, cth, sth1, cth1, &
             del, aa, sig, rr1, &
             rin, thin, sthin, cthin, delin, aain, sigin, rin1
  real(8) :: x1aa, x2aa, x3aa, x3t
  real(8) :: tmp_rand
  real(8) :: bsq_max, bsq_max1, bsq_max2, bsq_max2a, beta_act, &
             bnorm, pom1, rhoav, gfl, bbsq
  real(8) :: tmp4a, tmp4b, tmp4c, discr, alpha1
  real(8) :: tmp1a, tmp1b, tmp2a, tmp2b
  real(8) :: ucovph
  real(8) :: vx, vy, vz, vxn, vyn, vzn
! 
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), gcov1in(0:3,0:3), stat=merr)
  allocate(util(1:3), utiln(1:3), bcon(1:3), ucov(0:3), ucon(0:3), &
           bbcov(0:3), bbcon(0:3), beta1(1:3), stat=merr)
  allocate(ah(is1:ie1,js1:je1,ks1:ke1), pom(is1:ie1,js1:je1,ks1:ke1), &
           stat=merr)

!--------------------------------------------------------------------
!- Parameter -!

  rin=9.34d0
!  rin=9.0d0
  thin=pi/2.d0
  kappa=1.0d-3
  beta=0.d0
!  beta=100.d0
!  beta=0.d0
!  clm=4.5d0 ! a=0.0
  clm=4.35d0 ! a=0.5
!  clm=4.25d0 ! a=0.9
!  
  thedge=0.1d0
  thedge1=pi-0.1d0
  rcen=25.d0
!  
  rin1=1.d0
  rhomin=1.d-5 
  prmin=1.d-3*rhomin
!
  rhomax1=0.d0
  rhomax1a=0.d0
  prmax1=0.d0 
  prmax1a=0.d0
  prmax2=0.d0
!  
  bsq_max=0.d0
  bsq_max1=0.d0
  bsq_max2=0.d0
  bsq_max2a=0.d0
!
!- initally set minimum plasma everywhere -!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1

        do m=0,3
          do n=0,3
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo

        rr=x1(i)

        de=rhomin
!        de=1.d-3*exp(-3.*rr/rcen)
!        de=rhomin*((rr/rin1)**(-1.5d0))

!        pr=(gam-1.d0)*prmin*((rr/rin1)**(-2.5d0))
        pr=kpol*(de**gam)
!
!        pr=prmin
        util(1)=0.d0
!        util(1)=(gcon1(0,1)/gcon1(0,0))*(1.-(1./rr)**4)
        util(2)=0.d0
        util(3)=0.d0
        bcon(1)=0.d0
        bcon(2)=0.d0 
        bcon(3)=0.d0
!          
!- set primitive variables -!
!          
        uri(1,i,j,k)=de
        uri(2,i,j,k)=util(1)
        uri(3,i,j,k)=util(2)
        uri(4,i,j,k)=util(3)
        uri(5,i,j,k)=pr
        uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
        uri(7,i,j,k)=bcon(1)
        uri(8,i,j,k)=bcon(2)          
        uri(9,i,j,k)=bcon(3)
!       
      enddo
    enddo
  enddo
!  
!- calculation of BL Kerr metric at rin, thin -!
!  
  sthin=sin(thin)
  cthin=cos(thin)
  delin=rin*rin-2.*rin+akm*akm
  aain=(rin*rin+akm*akm)*(rin*rin+akm*akm)-delin*akm*akm*sthin*sthin
  sigin=rin*rin+akm*akm*cthin*cthin  
!
  do n=0,3
    do m=0,3
      gcov1in(m,n)=0.d0
    enddo
  enddo

  gcov1in(0,0)=-(sigin-2.*rin)/sigin               !- g_tt -!
  gcov1in(1,1)=sigin/delin                          !- g_rr -!
  gcov1in(2,2)=(aain/sigin)*sthin*sthin !- g_phi phi -!
  gcov1in(0,2)=-(2.*akm*rin*sthin*sthin)/sigin  !- g_t phi -!
  gcov1in(2,0)=gcov1in(0,2)                  !- g_phi t -!
  gcov1in(3,3)=sigin                      !- g_th th -!
! 
!- calculation of ut at rin -!
  tmp1a=gcov1in(0,2)*gcov1in(0,2)-gcov1in(0,0)*gcov1in(2,2)
  tmp1b=gcov1in(0,0)*clm*clm+2.*clm*gcov1in(0,2)+gcov1in(2,2)
  if(tmp1b .gt. 0.d0) then
    utin=-sqrt(tmp1a/tmp1b)
  else
    utin=-1.d0  
  endif
!
  if(myrank .eq. 0) then
    write(*,*) 'lin, utin=', clm, utin
  endif
!
!- calculate h -!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!
        rr=x1(i)
        ph=x2(j)
        th=x3(k) 
        th1=x3(k)-tilang

        if(th1 .le. 0.d0) then
          th1=2.*pi-th1
        elseif(th .gt. 2.*pi) then
          th1=th1-2.*pi
        endif   

        sph=sin(ph)
        cph=cos(ph)
        sth=sin(th)
        cth=cos(th) 
        sth1=sin(th1)
        cth1=cos(th1) 
!
        if(rr .ge. rin) then         
!- calculation of BL Kerr metric -! 
          del=(rr*rr)-(2.*rr)+(akm*akm)
          aa=(rr*rr+akm*akm)*(rr*rr+akm*akm)-del*akm*akm*sth1*sth1
          sig=rr*rr+akm*akm*cth1*cth1
!
          do n=0,3
            do m=0,3
              gcov1(m,n)=0.d0
              gcon1(m,n)=0.d0 
            enddo
          enddo
!
          gcov1(0,0)=-(sig-2.*rr)/sig               !- g^tt -!
          gcov1(1,1)=sig/del                          !- g^rr -!
          gcov1(2,2)=(aa/sig)*sth1*sth1             !- g^phi phi -!
          gcov1(0,2)=-(2.*akm*rr*sth1*sth1)/sig           !- g^t phi -!
          gcov1(2,0)=gcov1(0,2)                  !- g^phi t -!
          gcov1(3,3)=sig                      !- g^th th -!        
!
          gcon1(0,0)=-aa/(sig*del)               !- g^tt -!
          gcon1(1,1)=del/sig                          !- g^rr -!
          gcon1(2,2)=(sig-2.*rr)/(del*sig*sth1*sth1) !- g^phi phi -!
          gcon1(0,2)=-(2.*akm*rr)/(sig*del)           !- g^t phi -!
          gcon1(2,0)=gcon1(0,2)                  !- g^phi t -!
          gcon1(3,3)=1./sig 
!          
!- calculation of u_t at r,theta -!
!          
          tmp2a=gcov1(0,2)*gcov1(0,2)-gcov1(0,0)*gcov1(2,2)
          tmp2b=gcov1(0,0)*clm*clm+2.*clm*gcov1(0,2)+gcov1(2,2)
          if(tmp2b .gt. 0.d0) then
            ut=-sqrt(tmp2a/tmp2b)
          else
            ut=-1.d0  
          endif        

!
!- calculation of h = u_t,in/u_t
!          
          ah(i,j,k)=utin/ut
!
        else
          ah(i,j,k)=-1.d0
        endif
!
!- set torus structure -!
!
        if(ah(i,j,k) .gt. 1.d0 .and. rr .ge. rin) then    
!    
          hm1=max((ah(i,j,k) -1.)*(gam-1.)/gam, 0.d0)
          de=(hm1/kpol)**(1./(gam-1.))
          pr=kpol*(hm1/kpol)**(gam/(gam-1.))
!
!         if(j .eq. 4 .and. k .eq. 32) then
!           write(*,*) 'ah, de=', ah(i,j,k), de(i,j,k), i
!         endif
!
          util(1)=0.d0
          util(3)=0.d0
!
!- calculation of omega & contravariant u^t
!          
          omg=-(gcov1(0,2)+gcov1(0,0)*clm)/(gcov1(2,2)+gcov1(0,2)*clm)          
          uut=sqrt(-1./(gcov1(0,0)+2.*omg*gcov1(0,2)+omg*omg*gcov1(2,2)))
!          uut=-1./(ut*(1.-omg*clm))
!
          ucovph=-omg*ut

          util(2)=(omg*uut)
!          util(2)=gcon1(0,2)*ut+gcon1(2,2)*ucovph
!          util(2)=0.d0
!          
!- rotate (shifted) velicity
!
          vx=sth1*cph*util(1)+cth1*cph*util(3)-cth1*util(2)
          vy=sth1*sph*util(1)+cth1*sph*util(3)+cth1*util(2)
          vz=cth1*util(1)-sth1*util(3)
!          
          vxn=cos(tilang)*vx-sin(tilang)*vz
          vyn=vy
          vzn=sin(tilang)*vx+cos(tilang)*vz
!
          util(1)=sth*cph*vxn+sth*sph*vyn+cth*vzn
          util(2)=cth*cph*vxn+cth*sph*vyn-sth*vzn
          util(3)=-sph*vxn+cph*vyn
          
!
!- Put perturbation in pressure -!
!
!         call random_number(tmp_rand)
!!          write(*,*) "tmp_rand=", tmp_rand-0.5
!           pr=pr*(1.+0.05*(tmp_rand-0.5))
!           util(1)=util(1)+0.01*(tmp_rand-0.5)             
!
!- check density & pressure max -!
!
          if(de .gt. rhomax1) then
            rhomax1=de
!            write(*,*) 'de,i,k,rr,th=', de,i,k,x1(i),x3(k)
          endif
          if(pr .gt. prmax1) then
            prmax1=pr
          endif
!
!- transform velocity from BL to KS coordinates -!
!
          if(metric .eq. 303 .or. metric .eq. 403) then
            x1aa=x1(i)
            x2aa=x2(j)
!            x3aa=x3(k)+tilang   
            x3aa=x3(k)
            x3t=x3b(k)
            call transbl2ks(util,utiln,x1aa,x2aa,x3aa,x3t)
            util(1)=utiln(1)
            util(2)=utiln(2)
            util(3)=utiln(3)
!
          elseif(metric .eq. 203 .and. akm .ne. 0.d0) then
!
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
            ucon(1)=util(1)
            ucon(2)=util(2)
            ucon(3)=util(3)
!
            tmp4a=gcov1(0,0)
            tmp4b=2.*(gcov1(0,1)*ucon(1)+gcov1(0,2)*ucon(2)+gcov1(0,3)*ucon(3))
            tmp4c=1.+gcov1(1,1)*ucon(1)*ucon(1)+gcov1(2,2)*ucon(2)*ucon(2) &    
                  +gcov1(3,3)*ucon(3)*ucon(3) &
                  +2.*(gcov1(1,2)*ucon(1)*ucon(2)+gcov1(1,3)*ucon(1)*ucon(3) &  
                  +gcov1(2,3)*ucon(2)*ucon(3))   
            discr=tmp4b*tmp4b-4.*tmp4a*tmp4c
            ucon(0)=(-tmp4b-sqrt(discr))/(2.*tmp4a)
! 
            gfl=ucon(0)*alpha1
!
            util(1)=ucon(1)+gfl*beta1(1)/alpha1
            util(2)=ucon(2)+gfl*beta1(2)/alpha1
            util(3)=ucon(3)+gfl*beta1(3)/alpha1  
!
          endif
!
!          
!- set primitive variables -!
!          
          uri(1,i,j,k)=de
          uri(2,i,j,k)=util(1)
          uri(3,i,j,k)=util(2)
          uri(4,i,j,k)=util(3)
          uri(5,i,j,k)=pr
          uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
!          
        endif
!
      enddo
    enddo
  enddo
!
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_allreduce(rhomax1,rhomax1a,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
  call mpi_allreduce(prmax1,prmax1a,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
!
  if(myrank .eq. 0) then
    write(*,*) "rhomax1a, prmax1a=", rhomax1a, prmax1a
  endif
!
!
!- Normalize the densities and pressure -!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1 
        rr=x1(i)
        th=x3(k)
        if(ah(i,j,k) .gt. 1.d0 .and. rr .ge. rin) then
          de=uri(1,i,j,k)/rhomax1a
          pr=uri(5,i,j,k)/rhomax1a
!
          if(de .le. rhomin) then
            de=rhomin
          endif
          if(pr .le. prmin) then
            pr=prmin
          endif
!
          uri(1,i,j,k)=de
          uri(5,i,j,k)=pr
          uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
        endif
!
      enddo
    enddo
  enddo
!
   prmax2a=prmax1a/rhomax1a
   rhomax1a=1.d0
!
  if(myrank .eq. 0) then
    write(*,*) 'rhomax1a, prmax2a=', rhomax1a, prmax2a
  endif
!
!
!- Set vector potential at cell-center (r-theta frame)-!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
        pom(i,j,k)=0.d0
      enddo
    enddo
  enddo
!
!  do k=ks1+1,ke1-1
!    do j=js1,je1
!      do i=is1+1,ie1-1
!
!!        rr1=x1(i)
!!        th1=x3(k)
!!        rr1=(x1(i+1)+x1(i))/2.
!!        th1=(x3(k+1)+x3(k))/2.
!!        pom1=0.5*rr1*sin(th1)

!!!        rhoav=0.25*(uri(1,i-1,j,k)+uri(1,i,j,k-1) &
!!!                     +uri(1,i+1,j,k)+uri(1,i,j,k+1))
!         rhoav=0.25*(uri(1,i,j,k)+uri(1,i-1,j,k) &
!                     +uri(1,i,j,k-1)+uri(1,i-1,j,k-1))
!        rhoav=uri(1,i,j,k)
!
!        if(rhoav .ne. 0.d0) then
!          pom1=rhoav/rhomax1a - 0.3
!        else
!          pom1=0.d0
!        endif
!
!
!        if(pom1 .gt. 0.d0) then
!          pom(i,j,k)=pom1
!        endif
!
!      enddo
!    enddo
!  enddo
!
!- Set-up of Magnetic field -!
!  do k=ks1+1,ke1-1
!    do j=js1,je1
!      do i=is1+1,ie1-1
!        
!!        bcon(1)=-(pom(i,j,k)-pom(i,j,k+1)+pom(i+1,j,k)-pom(i+1,j,k+1)) &
!!                   /(2.*(x3(k+1)-x3(k))*detg(i,j,k))
!!
!!        bcon(3)=(pom(i,j,k)+pom(i,j,k+1)-pom(i+1,j,k)-pom(i+1,j,k+1)) &
!!                   /(2.*(x1(i+1)-x1(i))*detg(i,j,k))
!
!        tmp1a=pom(i,j,k)-pom(i,j,k+1)+pom(i+1,j,k)-pom(i+1,j,k+1)
!        if(tmp1a .ne. 0.d0) then
!          bcon(1)=-tmp1a/(2.*dx3(k)*detg(i,j,k))
!        else
!          bcon(1)=0.d0
!        endif
!!
!        tmp1b=pom(i,j,k)+pom(i,j,k+1)-pom(i+1,j,k)-pom(i+1,j,k+1)
!        if(tmp1b .ne. 0.d0) then
!          bcon(3)=tmp1b/(2.*dx1(i)*detg(i,j,k))
!        else
!          bcon(3)=0.d0
!        endif
!!
!!        if(k .eq. ks1) then
!!          bcon(1)=(pom(i,j,k+1)-pom(i,j,k))/((x3(k+1)-x3(k))*detg(i,j,k))
!!        elseif(k .eq. ke1) then
!!          bcon(1)=(pom(i,j,k)-pom(i,j,k-1))/((x3(k)-x3(k-1))*detg(i,j,k))
!!        else         
!!          bcon(1)=(pom(i,j,k+1)-pom(i,j,k-1))/((x3(k+1)-x3(k-1))*detg(i,j,k))
!!        endif
!
!!        if(i .eq. is1) then
!!          bcon(3)=-(pom(i+1,j,k)-pom(i,j,k))/((x1(i+1)-x1(i))*detg(i,j,k))
!!        elseif(i .eq. ie1) then
!!          bcon(3)=-(pom(i,j,k)-pom(i-1,j,k))/((x1(i)-x1(i-1))*detg(i,j,k))
!!        else         
!!          bcon(3)=-(pom(i+1,j,k)-pom(i-1,j,k))/((x1(i+1)-x1(i-1))*detg(i,j,k))
!!        endif
!
!        uri(7,i,j,k)=bcon(1)
!        uri(9,i,j,k)=bcon(3)
!
!      enddo
!    enddo
!  enddo
!
!- put boundary data for magnetic field -!
!
!  do j=js1,je1
!    do i=is1,ie1
!      uri(7,i,j,ks1)=uri(7,i,j,ks1+1)
!      uri(7,i,j,ke1)=uri(7,i,j,ke1-1)
!
!      uri(9,i,j,ks1)=uri(9,i,j,ks1+1)
!      uri(9,i,j,ke1)=uri(9,i,j,ke1-1)
!    enddo
!  enddo
!
!  do k=ks1,ke1
!    do j=js1,je1
!      uri(7,is1,j,k)=uri(7,is1+1,j,k)
!      uri(7,ie1,j,k)=uri(7,ie1-1,j,k)
!
!      uri(9,is1,j,k)=uri(9,is1+1,j,k)
!      uri(9,ie1,j,k)=uri(9,ie1-1,j,k)
!    enddo
!  enddo
!
!- calculation maximum magnetic field strength -!
!  do k=ks1,ke1
!    do j=js1,je1
!      do i=is1,ie1
!
!        util(1)=uri(2,i,j,k)
!        util(2)=uri(3,i,j,k)
!        util(3)=uri(4,i,j,k)
!        bcon(1)=uri(7,i,j,k)
!        bcon(2)=uri(8,i,j,k)
!        bcon(3)=uri(9,i,j,k)
!       
!        do m=0,3
!          do n=0,3
!            gcov1(m,n)=gcov(m,n,i,j,k)
!          enddo
!        enddo
!
!        call calgfl(util,gfl,gcov1)
!        call calucon(util,ucon,gfl,gcon1)
!        call lower(ucon,ucov,gcov1)
!        call calbbcon(bcon,bbcon,ucov,ucon)
!        call lower(bbcon,bbcov,gcov1)
!        call calbbsq(bbcon,bbcov,bbsq)
!
!        if(bbsq .gt. bsq_max) then
!          bsq_max=bbsq
!        endif
!
!      enddo
!    enddo
!  enddo    
!
!  call mpi_barrier(mpi_comm_world,merr)
!  call mpi_allreduce(bsq_max,bsq_max1,1,mpi_double_precision,mpi_max, &
!                     mpi_comm_world,merr)
!
!  beta_act=prmax2a/(0.5*bsq_max1)
!
!  if(beta .eq. 0.d0) then
!    bnorm=0.d0
!  else
!    bnorm=sqrt(beta_act/beta)
!  endif
!
!  if(myrank .eq. 0) then
!    write(*,*) "beta_act, bnorm=", beta_act, bnorm
!  endif
!
!- Normalize magnetic field strength -!
!
!  do k=ks1,ke1
!    do j=js1,je1
!      do i=is1,ie1
!
!        bcon(1)=uri(7,i,j,k)*bnorm            
!        bcon(2)=uri(8,i,j,k)*bnorm  
!        bcon(3)=uri(9,i,j,k)*bnorm  
!
!        util(1)=uri(2,i,j,k)
!        util(2)=uri(3,i,j,k)
!        util(3)=uri(4,i,j,k)
!
!        uri(7,i,j,k)=bcon(1)
!        uri(8,i,j,k)=bcon(2)
!        uri(9,i,j,k)=bcon(3)
!
!        do m=0,3
!          do n=0,3
!            gcov1(m,n)=gcov(m,n,i,j,k)
!          enddo
!        enddo        
!
!        call calgfl(util,gfl,gcov1)
!        call calucon(util,ucon,gfl,gcon1)
!        call lower(ucon,ucov,gcov1)
!        call calbbcon(bcon,bbcon,ucov,ucon)
!        call lower(bbcon,bbcov,gcov1)
!        call calbbsq(bbcon,bbcov,bbsq)
!
!        if(bbsq .gt. bsq_max2) then
!          bsq_max2=bbsq
!        endif
!
!      enddo
!    enddo
!  enddo 
!
!  call mpi_barrier(mpi_comm_world,merr)
!  call mpi_allreduce(bsq_max2,bsq_max2a,1,mpi_double_precision,mpi_max, &
!                     mpi_comm_world,merr)
!  beta_act=prmax2a/(0.5*bsq_max2a)
!
!  if(myrank .eq. 0) then
!    write(*,*) "beta_act=", beta_act
!  endif    
!
  deallocate(gcov1, gcon1, gcov1in, util, utiln, bcon, ucov, ucon, &
             bbcov, bbcon, beta1, ah, pom, stat=merr)
!
  return
end subroutine mdtildisk2
!
!--------------------------------------------------------------------
subroutine mdrecoilbh(uri,gcov,gcon,detg, &
                      x1,x2,x3,x1b,x2b,x3b,dx1,dx2,dx3, &
                      myrank,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
!  Recoiling BHs, Zanotti et al. (2010)
!  model=10 in pram.f90
!
  use pram, only : imax, jmax, kmax, nv, gam, metric, akm, pi, iter, kpol
  implicit none
  include 'mpif.h'

  integer :: i, j, k, m, n, is1 ,ie1, js1, je1, ks1, ke1
  integer :: merr, myrank, jh, nnn

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  
  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:), gcov1in(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:), beta1(:)
  real(8), allocatable :: util(:), utiln(:), bcon(:)
!
  real(8) :: x1(imax), x2(jmax), x3(kmax) 
  real(8) :: x1b(imax), x2b(jmax), x3b(kmax) 
  real(8) :: dx1(imax), dx2(jmax), dx3(kmax)
!
  real(8),allocatable :: ah(:,:,:), pom(:,:,:)
!
  real(8) :: de, pr  
  real(8) :: kappa, beta, clm, utin, ut, uut, &
             rhomin, prmin, hm1, omg, thedge, thedge1, &
             rhomax1, rhomax1a, prmax1, prmax1a, prmax2, prmax2a 
  real(8) :: rr, ph, th, sth, cth, del, aa, sig, rr1, th1, &
             rin, thin, sthin, cthin, delin, aain, sigin, rin1
  real(8) :: x1aa, x2aa, x3aa, x3t
  real(8) :: tmp_rand
  real(8) :: bsq_max, bsq_max1, bsq_max2, bsq_max2a, beta_act, &
             bnorm, pom1, rhoav, gfl, bbsq
  real(8) :: tmp4a, tmp4b, tmp4c, discr, alpha1
  real(8) :: tmp1a, tmp1b, tmp2a, tmp2b
  real(8) :: ucovph, small, vrecoil
! 
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), gcov1in(0:3,0:3), stat=merr)
  allocate(util(1:3), utiln(1:3), bcon(1:3), ucov(0:3), ucon(0:3), &
           bbcov(0:3), bbcon(0:3), beta1(1:3), stat=merr)
  allocate(ah(is1:ie1,js1:je1,ks1:ke1), pom(is1:ie1,js1:je1,ks1:ke1), &
           stat=merr)

!--------------------------------------------------------------------
!- Parameter -!

  rin=40.0d0
!  rin=9.0d0
  thin=pi/2.d0
  kappa=1.0d-3
  beta=0.d0
!  beta=100.d0
!  beta=0.d0
!  clm=4.5d0
  clm=8.0d0
!
  thedge=0.1*pi
  thedge1=0.9*pi
!  
  rin1=1.d0
  rhomin=1.d-6 
  prmin=1.d-3*rhomin
!
  rhomax1=0.d0
  rhomax1a=0.d0
  prmax1=0.d0
  prmax1a=0.d0
  prmax2=0.d0
!
  bsq_max=0.d0
  bsq_max1=0.d0
  bsq_max2=0.d0
  bsq_max2a=0.d0  
  
  vrecoil=1.d-3

  small=1.d-12
!  
!- initally set minimum plasma everywhere -!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1

        do m=0,3
          do n=0,3
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo

        rr=x1(i)

        de=rhomin
        pr=prmin
        util(1)=0.d0
!        util(1)=(gcon1(0,1)/gcon1(0,0))*(1.-(1./rr)**4)
        util(2)=0.d0
        util(3)=0.d0
        bcon(1)=0.d0
        bcon(2)=0.d0 
        bcon(3)=0.d0
!          
!- set primitive variables -!
!          
        uri(1,i,j,k)=de
        uri(2,i,j,k)=util(1)
        uri(3,i,j,k)=util(2)
        uri(4,i,j,k)=util(3)
        uri(5,i,j,k)=pr
        uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
        uri(7,i,j,k)=bcon(1)
        uri(8,i,j,k)=bcon(2)          
        uri(9,i,j,k)=bcon(3)
!       
      enddo
    enddo
  enddo
!  
!- calculation of BL Kerr metric at rin, thin -!
!  
!  sthin=sin(thin)
  sthin=1.d0
!  cthin=cos(thin)
  cthin=0.d0
  delin=rin*rin-2.*rin+akm*akm
  aain=(rin*rin+akm*akm)*(rin*rin+akm*akm)-delin*akm*akm*sthin*sthin
  sigin=rin*rin+akm*akm*cthin*cthin  
!
  do n=0,3
    do m=0,3
      gcov1in(m,n)=0.d0
    enddo
  enddo

  gcov1in(0,0)=-(sigin-2.*rin)*(1./sigin)               !- g_tt -!
  gcov1in(2,2)=(aain*(1./sigin))*sthin*sthin !- g_phi phi -!
  gcov1in(0,2)=-(2.*akm*rin*sthin*sthin)*(1./sigin)  !- g_t phi -!
! 
!- calculation of ut at rin -!
  tmp1a=gcov1in(0,2)*gcov1in(0,2)-gcov1in(0,0)*gcov1in(2,2)
  tmp1b=gcov1in(0,0)*clm*clm+2.*clm*gcov1in(0,2)+gcov1in(2,2)
  if(tmp1b .gt. 0.d0) then
    utin=-sqrt(tmp1a*(1./tmp1b))
  else
    utin=-1.d0  
  endif
!
  if(myrank .eq. 0) then
    write(*,*) 'lin, utin=', clm, utin
  endif
!
!- calculate h -!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!          
!- initialize -!
!         
        ah(i,j,k)=1.d0
        ucon(0)=0.d0
        ucon(1)=0.d0
        ucon(2)=0.d0
        ucon(3)=0.d0
        util(1)=0.d0
        util(2)=0.d0
        util(3)=0.d0
!
        rr=x1(i)
         
        if(x2(j) .lt. 0.d0) then
          ph=2.*pi-x2(j)
        elseif(x2(j) .gt. 2.*pi) then
          ph=x2(j)-2.*pi
        else   
          ph=x2(j) 
        endif

!        if(i .eq. 10 .and. k .eq. 4) then
!          write(*,*) "j, x2, ph=", j, x2(j), ph
!        endif   
        
        th=x3(k)
        sth=sin(th)
        if(abs(sth) .lt. small) then
          if(sth .ge. 0.d0) then
            sth=small
          elseif(sth .lt. 0.d0) then
            sth=-small 
          endif   
        endif   
        cth=cos(th) 
!
        if(rr .ge. rin .and. th .ge. thedge .and. th .le. thedge1) then         
!- calculation of BL Kerr metric -! 
          del=(rr*rr)-(2.*rr)+(akm*akm)
          aa=(rr*rr+akm*akm)*(rr*rr+akm*akm)-del*akm*akm*sth*sth
          sig=rr*rr+akm*akm*cth*cth
!
          do n=0,3
            do m=0,3
              gcov1(m,n)=0.d0
              gcon1(m,n)=0.d0 
            enddo
          enddo
!
          gcov1(0,0)=-(sig-2.*rr)*(1./sig)               !- g^tt -!
          gcov1(1,1)=sig*(1./del)                          !- g^rr -!
          gcov1(2,2)=(aa*(1./sig))*sth*sth             !- g^phi phi -!
          gcov1(0,2)=-(2.*akm*rr*sth*sth)*(1./sig)           !- g^t phi -!
          gcov1(2,0)=gcov1(0,2)                  !- g^phi t -!
          gcov1(3,3)=sig                      !- g^th th -!        
!
          gcon1(0,0)=-aa*(1./(sig*del))               !- g^tt -!
          gcon1(1,1)=del*(1./sig)                          !- g^rr -!
          gcon1(2,2)=(sig-2.*rr)*(1./(del*sig*sth*sth)) !- g^phi phi -!
          gcon1(0,2)=-(2.*akm*rr)*(1./(sig*del))           !- g^t phi -!
          gcon1(2,0)=gcon1(0,2)                  !- g^phi t -!
          gcon1(3,3)=1./sig 
!          
!- calculation of u_t at r,theta -!
!          
          tmp2a=gcov1(0,2)*gcov1(0,2)-gcov1(0,0)*gcov1(2,2)
          tmp2b=gcov1(0,0)*clm*clm+2.*clm*gcov1(0,2)+gcov1(2,2)
          if(tmp2b .gt. 0.d0) then
            ut=-sqrt(tmp2a*(1./tmp2b))
          else
            ut=-1.d0  
          endif        

!
!- calculation of h = u_t,in/u_t
!          
          ah(i,j,k)=utin*(1./ut)
!
        else
          ah(i,j,k)=1.d0
        endif
!
!- set torus structure -!
!
        if(ah(i,j,k) .gt. 1.d0 .and. rr .ge. rin &
           .and. th .ge. thedge .and. th .le. thedge1) then    
!    
          hm1=max((ah(i,j,k) -1.)*(gam-1.)/gam, 1.d-10)
          de=(hm1*(1./kpol))**(1./(gam-1.))
          pr=kpol*(hm1*(1./kpol))**(gam/(gam-1.))
!
!         if(j .eq. 4 .and. k .eq. 32) then
!           write(*,*) 'ah, de=', ah(i,j,k), de(i,j,k), i
!         endif
!
          ucon(1)=0.d0
          ucon(3)=0.d0
!
!- calculation of omega & contravariant u^t
!          
          omg=-(gcov1(0,2)+gcov1(0,0)*clm)*(1./(gcov1(2,2)+gcov1(0,2)*clm))
          uut=sqrt(-1./(gcov1(0,0)+2.*omg*gcov1(0,2)+omg*omg*gcov1(2,2)))
!          uut=-1./(ut*(1.-omg*clm))
!
          ucovph=-omg*ut

          ucon(0)=uut
          ucon(2)=(omg*uut)
!          util(2)=gcon1(0,2)*ut+gcon1(2,2)*ucovph
!          util(2)=0.d0
!
!- cal spacial 4-velocity -!
!
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
          gfl=ucon(0)*alpha1

          util(1)=ucon(1)+gfl*beta1(1)*(1./alpha1)
          util(2)=ucon(2)+gfl*beta1(2)*(1./alpha1)
          util(3)=ucon(3)+gfl*beta1(3)*(1./alpha1)          
!
!- Put perturbation in pressure -!
!
!         call random_number(tmp_rand)
!!          write(*,*) "tmp_rand=", tmp_rand-0.5
!          pr=pr*(1.+0.05*(tmp_rand-0.5))
!
!- Put recoiling speed -!
!
          util(1)=util(1)+vrecoil*cos(ph)*sth*(1./sqrt(gcov1(1,1)))
          util(2)=util(2)-vrecoil*sin(ph)*(1./sqrt(gcov1(2,2)))
          util(3)=util(3)+vrecoil*cos(ph)*cth*(1./sqrt(gcov1(3,3)))
!          
!- check density & pressure max -!
!
          if(de .gt. rhomax1) then
            rhomax1=de
!            write(*,*) 'de,i,k,rr,th=', de,i,k,x1(i),x3(k)
          endif
          if(pr .gt. prmax1) then
            prmax1=pr
          endif         
!
!- transform velocity from BL to KS coordinates -!
!
!          if(metric .eq. 303) then
!            x1aa=x1(i)
!            x2aa=x2(j)
!            x3aa=x3(k)   
!            x3t=x3b(k)
!            call transbl2ks(util,utiln,x1aa,x2aa,x3aa,x3t)
!            util(1)=utiln(1)
!            util(2)=utiln(2)
!            util(3)=utiln(3)
!
!          elseif(metric .eq. 203 .and. akm .ne. 0.d0) then
!
!            do m=0,3
!              do n=0,3
!                gcov1(m,n)=gcov(m,n,i,j,k)
!                gcon1(m,n)=gcon(m,n,i,j,k)
!              enddo
!            enddo
!
!            alpha1=1./sqrt(-gcon1(0,0))
!
!            beta1(1)=alpha1*alpha1*gcon1(0,1)
!            beta1(2)=alpha1*alpha1*gcon1(0,2)
!            beta1(3)=alpha1*alpha1*gcon1(0,3)
!
!            ucon(1)=util(1)
!            ucon(2)=util(2)
!            ucon(3)=util(3)
!
!            tmp4a=gcov1(0,0)
!            tmp4b=2.*(gcov1(0,1)*ucon(1)+gcov1(0,2)*ucon(2)+gcov1(0,3)*ucon(3))
!            tmp4c=1.+gcov1(1,1)*ucon(1)*ucon(1)+gcov1(2,2)*ucon(2)*ucon(2) &   ! 
!                  +gcov1(3,3)*ucon(3)*ucon(3) &
!                  +2.*(gcov1(1,2)*ucon(1)*ucon(2)+gcov1(1,3)*ucon(1)*ucon(3) & ! 
!                  +gcov1(2,3)*ucon(2)*ucon(3))   
!            discr=tmp4b*tmp4b-4.*tmp4a*tmp4c
!            ucon(0)=(-tmp4b-sqrt(discr))/(2.*tmp4a)
! 
!            gfl=ucon(0)*alpha1
!
!            util(1)=ucon(1)+gfl*beta1(1)/alpha1
!            util(2)=ucon(2)+gfl*beta1(2)/alpha1
!            util(3)=ucon(3)+gfl*beta1(3)/alpha1  
!
!          endif
!
!          
!- set primitive variables -!
!          
          uri(1,i,j,k)=de
          uri(2,i,j,k)=util(1)
          uri(3,i,j,k)=util(2)
          uri(4,i,j,k)=util(3)
          uri(5,i,j,k)=pr
          uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
!          
        endif
!
      enddo
    enddo
  enddo
!
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_allreduce(rhomax1,rhomax1a,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
  call mpi_allreduce(prmax1,prmax1a,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)

  if(myrank .eq. 0) then
    write(*,*) "rhomax1a, prmax1a=", rhomax1a, prmax1a
  endif
!
!
!- Normalize the densities and pressure -!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1 
        rr=x1(i)
        th=x3(k)
        if(ah(i,j,k) .gt. 1.d0 .and. rr .ge. rin &
           .and. th .ge. thedge .and. th .le. thedge1) then
!        if(ah(i,j,k) .gt. 1.d0 .and. rr .ge. rin) then
!          de=uri(1,i,j,k)/rhomax1a
!          pr=uri(5,i,j,k)/rhomax1a
!          
          de=uri(1,i,j,k)/1.d0
          pr=uri(5,i,j,k)/1.d0   
!
          if(de .le. rhomin) then
            de=rhomin
          endif
          if(pr .le. prmin) then
            pr=prmin
          endif
!
          uri(1,i,j,k)=de
          uri(5,i,j,k)=pr
          uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
        endif
!
      enddo
    enddo
  enddo
!
!   prmax2a=prmax1a/rhomax1a
   prmax2a=prmax1a
   rhomax1a=1.d0
!
  if(myrank .eq. 0) then
    write(*,*) 'rhomax1a, prmax2a=', rhomax1a, prmax2a
  endif   
!
!- put recoiling BH speed -!
!  
!  do k=ks1,ke1
!    do j=js1,je1
!      do i=is1,ie1
!          
!        uri(2,i,j,k)=uri(2,i,j,k)+vrecoil
!         
!      enddo
!    enddo  
!  enddo
!
!
!
!- Set vector potential at cell-center (r-theta frame)-!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
        pom(i,j,k)=0.d0
      enddo
    enddo
  enddo
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
  
!  do k=ks1+1,ke1-1
!    do j=js1,je1
!      do i=is1+1,ie1-1
!!
!!        rr1=x1(i)
!!        th1=x3(k)
!!        rr1=(x1(i+1)+x1(i))/2.
!!        th1=(x3(k+1)+x3(k))/2.
!!        pom1=0.5*rr1*sin(th1)

!        rhoav=0.25*(uri(1,i-1,j,k)+uri(1,i,j,k-1) &
!                     +uri(1,i+1,j,k)+uri(1,i,j,k+1))
!         rhoav=0.25*(uri(1,i,j,k)+uri(1,i-1,j,k) &
!                     +uri(1,i,j,k-1)+uri(1,i-1,j,k-1))
        rhoav=uri(1,i,j,k)
!
        if(rhoav .ne. 0.d0) then
          pom1=(rhoav*(1./rhomax1a)) - 0.3d0
        else
          pom1=0.d0
        endif
!
!
        if(pom1 .gt. 0.d0) then
          pom(i,j,k)=pom1
        endif
!
      enddo
    enddo
  enddo
!
!
!- put boundary data for magnetic field -!
!
!  do j=js1,je1
!    do i=is1,ie1
!      pom(i,j,ks1)=pom(i,j,ks1+1)
!      pom(i,j,ke1)=pom(i,j,ke1-1)
!    enddo
!  enddo
!
!  do k=ks1,ke1
!    do j=js1,je1
!      pom(is1,j,k)=pom(is1+1,j,k)
!      pom(ie1,j,k)=pom(ie1-1,j,k)
!    enddo
!  enddo
!

  !
!- Set-up of Magnetic field -!
  do k=ks1+1,ke1-1
    do j=js1,je1
      do i=is1+1,ie1-1
!        
!!        bcon(1)=-(pom(i,j,k)-pom(i,j,k+1)+pom(i+1,j,k)-pom(i+1,j,k+1)) &
!!                   /(2.*(x3(k+1)-x3(k))*detg(i,j,k))
!!
!!        bcon(3)=(pom(i,j,k)+pom(i,j,k+1)-pom(i+1,j,k)-pom(i+1,j,k+1)) &
!!                   /(2.*(x1(i+1)-x1(i))*detg(i,j,k))
!
!        tmp1a=pom(i,j,k)-pom(i,j,k+1)+pom(i+1,j,k)-pom(i+1,j,k+1)
        tmp1a=pom(i,j,k+1)-pom(i,j,k-1)
        if(tmp1a .ne. 0.d0) then
!          bcon(1)=-tmp1a/(2.*dx3(k)*detg(i,j,k))
          bcon(1)=-tmp1a*(1./((x3(k+1)-x3(k-1))*detg(i,j,k))) 
        else
          bcon(1)=0.d0
        endif
!!
!        tmp1b=pom(i,j,k)+pom(i,j,k+1)-pom(i+1,j,k)-pom(i+1,j,k+1)
        tmp1b=pom(i+1,j,k)-pom(i-1,j,k)
        if(tmp1b .ne. 0.d0) then
!          bcon(3)=tmp1b/(2.*dx1(i)*detg(i,j,k))
          bcon(3)=tmp1b*(1./((x1(i+1)-x1(i-1))*detg(i,j,k))) 
        else
          bcon(3)=0.d0
        endif
!
!!        if(k .eq. ks1) then
!!          bcon(1)=(pom(i,j,k+1)-pom(i,j,k))/((x3(k+1)-x3(k))*detg(i,j,k))
!!        elseif(k .eq. ke1) then
!!          bcon(1)=(pom(i,j,k)-pom(i,j,k-1))/((x3(k)-x3(k-1))*detg(i,j,k))
!!        else         
!!          bcon(1)=(pom(i,j,k+1)-pom(i,j,k-1))/((x3(k+1)-x3(k-1))*detg(i,j,k))
!!        endif
!
!!        if(i .eq. is1) then
!!          bcon(3)=-(pom(i+1,j,k)-pom(i,j,k))/((x1(i+1)-x1(i))*detg(i,j,k))
!!        elseif(i .eq. ie1) then
!!          bcon(3)=-(pom(i,j,k)-pom(i-1,j,k))/((x1(i)-x1(i-1))*detg(i,j,k))
!!        else         
!!          bcon(3)=-(pom(i+1,j,k)-pom(i-1,j,k))/((x1(i+1)-x1(i-1))*detg(i,j,k))
! !       endif
!
        uri(7,i,j,k)=bcon(1)
        uri(9,i,j,k)=bcon(3)
!
      enddo
    enddo
  enddo
!
!- put boundary data for magnetic field -!
!
  do j=js1,je1
    do i=is1,ie1
      uri(7,i,j,ks1)=uri(7,i,j,ks1+1)
      uri(7,i,j,ke1)=uri(7,i,j,ke1-1)
!
      uri(9,i,j,ks1)=uri(9,i,j,ks1+1)
      uri(9,i,j,ke1)=uri(9,i,j,ke1-1)
    enddo
  enddo
!
  do k=ks1,ke1
    do j=js1,je1
      uri(7,is1,j,k)=uri(7,is1+1,j,k)
      uri(7,ie1,j,k)=uri(7,ie1-1,j,k)
!
      uri(9,is1,j,k)=uri(9,is1+1,j,k)
      uri(9,ie1,j,k)=uri(9,ie1-1,j,k)
    enddo
  enddo
!
!- calculation maximum magnetic field strength -!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
        bcon(1)=uri(7,i,j,k)
        bcon(2)=uri(8,i,j,k)
        bcon(3)=uri(9,i,j,k)
!       
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
          enddo
        enddo
!
        call calgfl(util,gfl,gcov1)
        call calucon(util,ucon,gfl,gcon1)
        call lower(ucon,ucov,gcov1)
        call calbbcon(bcon,bbcon,ucov,ucon)
        call lower(bbcon,bbcov,gcov1)
        call calbbsq(bbcon,bbcov,bbsq)
!
        if(bbsq .gt. bsq_max) then
          bsq_max=bbsq
        endif
!
      enddo
    enddo
  enddo    
!
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_allreduce(bsq_max,bsq_max1,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
!
  if(bsq_max1 .eq. 0.d0) then
    beta_act=0.d0 
  else   
    beta_act=prmax2a*(1./(0.5*bsq_max1))
  endif   
!
  if(beta .eq. 0.d0) then
    bnorm=0.d0
  else
    bnorm=sqrt(beta_act*(1./beta))
  endif
!
  if(myrank .eq. 0) then
    write(*,*) "beta_act, bnorm=", beta_act, bnorm
  endif
!
!- Normalize magnetic field strength -!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!!
        bcon(1)=uri(7,i,j,k)*bnorm            
        bcon(2)=uri(8,i,j,k)*bnorm  
        bcon(3)=uri(9,i,j,k)*bnorm  
!
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
!
        uri(7,i,j,k)=bcon(1)
        uri(8,i,j,k)=bcon(2)
        uri(9,i,j,k)=bcon(3)
!
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
          enddo
        enddo        
!
        call calgfl(util,gfl,gcov1)
        call calucon(util,ucon,gfl,gcon1)
        call lower(ucon,ucov,gcov1)
        call calbbcon(bcon,bbcon,ucov,ucon)
        call lower(bbcon,bbcov,gcov1)
        call calbbsq(bbcon,bbcov,bbsq)
!
        if(bbsq .gt. bsq_max2) then
          bsq_max2=bbsq
        endif
!
      enddo
    enddo
  enddo 
!
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_allreduce(bsq_max2,bsq_max2a,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
  if(bsq_max2a .eq. 0.d0) then
    beta_act=0.d0 
  else   
    beta_act=prmax2a*(1./(0.5*bsq_max2a))
   endif  
!
  if(myrank .eq. 0) then
    write(*,*) "beta_act=", beta_act
  endif    
!  
  deallocate(gcov1, gcon1, gcov1in, util, utiln, bcon, ucov, ucon, &
             bbcov, bbcon, beta1, ah, pom, stat=merr)
!
  return
end subroutine mdrecoilbh
!
!--------------------------------------------------------------------
subroutine mdrecoilbh1a(uri,gcov,gcon,detg, &
                    x1,x2,x3,x1b,x2b,x3b,dx1,dx2,dx3, &
                    myrank,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
!  Recoiling BHs, Zanotti et al. (2010)
!  model=10 in pram.f90
!
  use pram, only : imax, jmax, kmax, nv, gam, metric, akm, pi, iter, kpol
  implicit none
  include 'mpif.h'

  integer :: i, j, k, m, n, is1 ,ie1, js1, je1, ks1, ke1
  integer :: merr, myrank, jh, nnn

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  
  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:), gcov1in(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:), beta1(:)
  real(8), allocatable :: util(:), utiln(:), bcon(:)
!
  real(8) :: x1(imax), x2(jmax), x3(kmax) 
  real(8) :: x1b(imax), x2b(jmax), x3b(kmax) 
  real(8) :: dx1(imax), dx2(jmax), dx3(kmax)
!
  real(8),allocatable :: ah(:,:,:), pom(:,:,:)
!
  real(8) :: de, pr  
  real(8) :: kappa, beta, clm, utin, ut, uut, &
             rhomin, prmin, hm1, omg, thedge, thedge1, rcen, &
             rhomax1, rhomax1a, prmax1, prmax1a, prmax2, prmax2a 
  real(8) :: rr, th, sth, cth, del, aa, sig, rr1, th1, &
             rin, thin, sthin, cthin, delin, aain, sigin, rin1
  real(8) :: x1aa, x2aa, x3aa, x3t
  real(8) :: tmp_rand
  real(8) :: bsq_max, bsq_max1, bsq_max2, bsq_max2a, beta_act, &
             bnorm, pom1, rhoav, gfl, bbsq
  real(8) :: tmp4a, tmp4b, tmp4c, discr, alpha1
  real(8) :: tmp1a, tmp1b, tmp2a, tmp2b
  real(8) :: ucovph, vrecoil
  real(8) :: small
! 
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), gcov1in(0:3,0:3), stat=merr)
  allocate(util(1:3), utiln(1:3), bcon(1:3), ucov(0:3), ucon(0:3), &
           bbcov(0:3), bbcon(0:3), beta1(1:3), stat=merr)
  allocate(ah(is1:ie1,js1:je1,ks1:ke1), pom(is1:ie1,js1:je1,ks1:ke1), &
           stat=merr)

!--------------------------------------------------------------------
!- Parameter -!

  rin=40.d0
!  rin=9.0d0
  thin=pi/2.d0
  kappa=1.0d-3
  beta=0.d0
!  beta=50.d0
!  beta=0.d0
  clm=8.0d0
!  
  thedge=0.1d0*pi
  thedge1=0.9d0*pi
  rcen=25.d0
!  
  rin1=1.d0
!  rhomin=1.d-4
  rhomin=1.d-5 
!  prmin=0.5*rhomin
  prmin=1.d-3*rhomin
!
  rhomax1=0.d0
  rhomax1a=0.d0
  prmax1=0.d0 
  prmax1a=0.d0
  prmax2=0.d0
!  
  bsq_max=0.d0
  bsq_max1=0.d0
  bsq_max2=0.d0
  bsq_max2a=0.d0
!
  small=1.d-12  
!
!- initally set minimum plasma everywhere -!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1

        do m=0,3
          do n=0,3
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
        
        rr=x1(i)

        de=rhomin
!        de=rhomin*exp(-3.*rr/rcen)
!        de=rhomin*rr**(-3./2.)

!        pr=(gam-1.d0)*prmin*((rr/rin1)**(-2.5d0))
!        pr=kpol*(de**gam)
        pr=prmin
!        pr=prmin*rr**(-5./2.)
        
        util(1)=0.d0
!        util(1)=(gcon1(0,1)/gcon1(0,0))*(1.-(1./rr)**4)
!        util(1)=beta1(1)*(1./alpha1)        
        util(2)=0.d0
        util(3)=0.d0
        bcon(1)=0.d0
        bcon(2)=0.d0 
        bcon(3)=0.d0
!          
!- set primitive variables -!
!          
        uri(1,i,j,k)=de
        uri(2,i,j,k)=util(1)
        uri(3,i,j,k)=util(2)
        uri(4,i,j,k)=util(3)
        uri(5,i,j,k)=pr
        uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
        uri(7,i,j,k)=bcon(1)
        uri(8,i,j,k)=bcon(2)          
        uri(9,i,j,k)=bcon(3)
!       
      enddo
    enddo
  enddo
!  
!- calculation of BL Kerr metric at rin, thin -!
!  
!  sthin=sin(thin)
  sthin=1.d0
!  cthin=cos(thin)
  cthin=0.d0
  delin=rin*rin-2.*rin+akm*akm
  aain=(rin*rin+akm*akm)*(rin*rin+akm*akm)-delin*akm*akm*sthin*sthin
  sigin=rin*rin+akm*akm*cthin*cthin  
!
  do n=0,3
    do m=0,3
      gcov1in(m,n)=0.d0
    enddo
  enddo

  gcov1in(0,0)=-(sigin-2.*rin)*(1./sigin)               !- g_tt -!
  gcov1in(2,2)=(aain*(1./sigin))*sthin*sthin !- g_phi phi -!
  gcov1in(0,2)=-(2.*akm*rin*sthin*sthin)*(1./sigin)  !- g_t phi -!
! 
!- calculation of ut at rin -!
  tmp1a=gcov1in(0,2)*gcov1in(0,2)-gcov1in(0,0)*gcov1in(2,2)
  tmp1b=gcov1in(0,0)*clm*clm+2.*clm*gcov1in(0,2)+gcov1in(2,2)
  if(tmp1b .gt. 0.d0 .and. tmp1a .gt. 0.d0) then
    utin=-sqrt(tmp1a*(1./tmp1b))
  else
    utin=-1.d0  
  endif
!
  if(myrank .eq. 0) then
    write(*,*) 'lin, utin=', clm, utin
  endif
!
!- calculate h -!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!          
!- initialize -!
!         
        ah(i,j,k)=1.d0
!         
        rr=x1(i)
        th=x3(k)
        
        sth=sin(th)
        if(abs(sth) .lt. small) then
          if(sth .ge. 0.d0) then
            sth=small
          elseif(sth .lt. 0.d0) then
            sth=-small 
          endif   
        endif

        cth=cos(th) 
!
!        if(rr .ge. rin) then         
        if(rr .ge. rin .and. th .ge. thedge .and. th .le. thedge1) then
!- calculation of BL Kerr metric -! 
          del=(rr*rr)-(2.*rr)+(akm*akm)
          aa=(rr*rr+akm*akm)*(rr*rr+akm*akm)-del*akm*akm*sth*sth
          sig=rr*rr+akm*akm*cth*cth
!
          do n=0,3
            do m=0,3
              gcov1(m,n)=0.d0
              gcon1(m,n)=0.d0 
            enddo
          enddo
!
          gcov1(0,0)=-(sig-2.*rr)*(1./sig)               !- g^tt -!
          gcov1(1,1)=sig*(1./del)                          !- g^rr -!
          gcov1(2,2)=(aa*(1./sig))*sth*sth             !- g^phi phi -!
          gcov1(0,2)=-(2.*akm*rr*sth*sth)*(1./sig)           !- g^t phi -!
          gcov1(2,0)=gcov1(0,2)                  !- g^phi t -!
          gcov1(3,3)=sig                      !- g^th th -!        
!
          gcon1(0,0)=-aa*(1./(sig*del))               !- g^tt -!
          gcon1(1,1)=del*(1./sig)                          !- g^rr -!
          gcon1(2,2)=(sig-2.*rr)*(1./(del*sig*sth*sth)) !- g^phi phi -!
          gcon1(0,2)=-(2.*akm*rr)*(1./(sig*del))           !- g^t phi -!
          gcon1(2,0)=gcon1(0,2)                  !- g^phi t -!
          gcon1(3,3)=1./sig 
!          
!- calculation of u_t at r,theta -!
!          
          tmp2a=gcov1(0,2)*gcov1(0,2)-gcov1(0,0)*gcov1(2,2)
          tmp2b=gcov1(0,0)*clm*clm+2.*clm*gcov1(0,2)+gcov1(2,2)
          if(tmp2b .gt. 0.d0 .and. tmp2a .gt. 0.d0) then
            ut=-sqrt(tmp2a*(1./tmp2b))
          else
            ut=-1.d0  
          endif        

!
!- calculation of h = u_t,in/u_t
!          
          ah(i,j,k)=utin*(1./ut)
!
        else
          ah(i,j,k)=1.d0
        endif
!
!- set torus structure -!
!
        if(ah(i,j,k) .gt. 1.d0 .and. rr .ge. rin  &
           .and. th .ge. thedge .and. th .le. thedge1) then    
!         if((ah(i,j,k)-1.d0) .gt. 1.d-7 .and. rr .ge. rin) then
!    
          hm1=max((ah(i,j,k) -1.)*(gam-1.)/gam, 0.d0)
          de=(hm1*(1./kpol))**(1./(gam-1.))
          pr=kpol*(hm1*(1./kpol))**(gam/(gam-1.))
!
          ucon(1)=0.d0
          ucon(3)=0.d0
!
!- calculation of omega & contravariant u^t
!          
          omg=-(gcov1(0,2)+gcov1(0,0)*clm)*(1./(gcov1(2,2)+gcov1(0,2)*clm))
          uut=sqrt(-1./(gcov1(0,0)+2.*omg*gcov1(0,2)+omg*omg*gcov1(2,2)))
!          uut=-1./(ut*(1.-omg*clm))
!
          ucovph=-omg*ut

          ucon(0)=uut
          ucon(2)=(omg*uut)
!          ucon(2)=gcon1(0,2)*ut+gcon1(2,2)*ucovph
!          ucon(2)=0.d0
!
!- cal spacial 4-velocity -!
!
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
          gfl=ucon(0)*alpha1

          util(1)=ucon(1)+gfl*beta1(1)*(1./alpha1)
          util(2)=ucon(2)+gfl*beta1(2)*(1./alpha1)
          util(3)=ucon(3)+gfl*beta1(3)*(1./alpha1)             
!
!- Put perturbation in pressure -!
!
         call random_number(tmp_rand)
!!          write(*,*) "tmp_rand=", tmp_rand-0.5
!           pr=pr*(1.+0.001*(tmp_rand-0.5))
!           util(1)=util(1)+0.01*(tmp_rand-0.5)             
!
!- check density & pressure max -!
!
          if(de .gt. rhomax1) then
            rhomax1=de
!            write(*,*) 'de,i,k,rr,th=', de,i,k,x1(i),x3(k)
          endif
          if(pr .gt. prmax1) then
            prmax1=pr
          endif
!
!- transform velocity from BL to KS coordinates (do not need?)-!
!
!          if(metric .eq. 303 .or. metric .eq. 403) then
!            x1aa=x1(i)
!            x2aa=x2(j)
!            x3aa=x3(k)   
!            x3t=x3b(k)
!            call transbl2ks(util,utiln,x1aa,x2aa,x3aa,x3t)
!            util(1)=utiln(1)
!            util(2)=utiln(2)
!            util(3)=utiln(3)
!
!          elseif(metric .eq. 203 .and. akm .ne. 0.d0) then
!
!            do m=0,3
!              do n=0,3
!                gcov1(m,n)=gcov(m,n,i,j,k)
!                gcon1(m,n)=gcon(m,n,i,j,k)
!              enddo
!            enddo
!
!            alpha1=1./sqrt(-gcon1(0,0))
!
!            beta1(1)=alpha1*alpha1*gcon1(0,1)
!            beta1(2)=alpha1*alpha1*gcon1(0,2)
!            beta1(3)=alpha1*alpha1*gcon1(0,3)
!
!            ucon(1)=util(1)
!            ucon(2)=util(2)
!            ucon(3)=util(3)
!
!            tmp4a=gcov1(0,0)
!            tmp4b=2.*(gcov1(0,1)*ucon(1)+gcov1(0,2)*ucon(2)+gcov1(0,3)*ucon(3))
!            tmp4c=1.+gcov1(1,1)*ucon(1)*ucon(1)+gcov1(2,2)*ucon(2)*ucon(2) &    
!                  +gcov1(3,3)*ucon(3)*ucon(3) &
!                  +2.*(gcov1(1,2)*ucon(1)*ucon(2)+gcov1(1,3)*ucon(1)*ucon(3) &  
!                  +gcov1(2,3)*ucon(2)*ucon(3))   
!            discr=tmp4b*tmp4b-4.*tmp4a*tmp4c
!            ucon(0)=(-tmp4b-sqrt(discr))/(2.*tmp4a)
! 
!            gfl=ucon(0)*alpha1
!
!            util(1)=ucon(1)+gfl*beta1(1)/alpha1
!            util(2)=ucon(2)+gfl*beta1(2)/alpha1
!            util(3)=ucon(3)+gfl*beta1(3)/alpha1  
!
!          endif
!
!          
!- set primitive variables -!
!          
          uri(1,i,j,k)=de
          uri(2,i,j,k)=util(1)
          uri(3,i,j,k)=util(2)
          uri(4,i,j,k)=util(3)
          uri(5,i,j,k)=pr
          uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
!          
        endif
!
      enddo
    enddo
  enddo
!
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_allreduce(rhomax1,rhomax1a,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
  call mpi_allreduce(prmax1,prmax1a,1,mpi_double_precision,mpi_max, &
                     mpi_comm_world,merr)
!
  if(myrank .eq. 0) then
    write(*,*) "rhomax1a, prmax1a=", rhomax1a, prmax1a
  endif
!
!
!- Normalize the densities and pressure -!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1 
        rr=x1(i)
        th=x3(k)
        if(ah(i,j,k) .gt. 1.d0 .and. rr .ge. rin &
           .and. th .ge. thedge .and. th .le. thedge1) then
!        if((ah(i,j,k)-1.d0) .gt. 1.d-7 .and. rr .ge. rin) then
!          de=uri(1,i,j,k)/rhomax1a
!          pr=uri(5,i,j,k)/rhomax1a
!
          de=uri(1,i,j,k)/1.d0
          pr=uri(5,i,j,k)/1.d0          
!
          if(de .le. rhomin) then
            de=rhomin
          endif
          if(pr .le. prmin) then
            pr=prmin
          endif
!
          uri(1,i,j,k)=de
          uri(5,i,j,k)=pr
          uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
        endif
!
      enddo
    enddo
  enddo
!

!   prmax2a=prmax1a/rhomax1a
   prmax2a=prmax1a
   rhomax1a=1.d0
!
  if(myrank .eq. 0) then
    write(*,*) 'rhomax1a, prmax2a=', rhomax1a, prmax2a
  endif
!
!
!- Set vector potential at cell-center (r-theta frame)-!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
        pom(i,j,k)=0.d0
      enddo
    enddo
  enddo
!
!  do k=ks1,ke1
!    do j=js1,je1
!      do i=is1,ie1
  
  do k=ks1+1,ke1-1
    do j=js1,je1
      do i=is1+1,ie1-1
!!
!!        rr1=x1(i)
!!        th1=x3(k)
!!        rr1=(x1(i+1)+x1(i))/2.
!!        th1=(x3(k+1)+x3(k))/2.
!!        pom1=0.5*rr1*sin(th1)

!        rhoav=0.25*(uri(1,i-1,j,k)+uri(1,i,j,k-1) &
!                     +uri(1,i+1,j,k)+uri(1,i,j,k+1))
!!         rhoav=0.25*(uri(1,i,j,k)+uri(1,i-1,j,k) &
!!                     +uri(1,i,j,k-1)+uri(1,i-1,j,k-1))
!!        rhoav=uri(1,i,j,k)
!!
!        if(rhoav .ne. 0.d0) then
!          pom1=(rhoav*(1./rhomax1a)) - 0.3d0
!        else
!          pom1=0.d0
!        endif
!
!
!        if(pom1 .gt. 0.d0) then
!          pom(i,j,k)=pom1
!        endif
!
      enddo
    enddo
  enddo
!
!- Set-up of Magnetic field -!
!  do k=ks1+1,ke1-1
!    do j=js1,je1
!      do i=is1+1,ie1-1
!        
!!        bcon(1)=-(pom(i,j,k)-pom(i,j,k+1)+pom(i+1,j,k)-pom(i+1,j,k+1)) &
!!                   /(2.*(x3(k+1)-x3(k))*detg(i,j,k))
!!
!!        bcon(3)=(pom(i,j,k)+pom(i,j,k+1)-pom(i+1,j,k)-pom(i+1,j,k+1)) &
!!                   /(2.*(x1(i+1)-x1(i))*detg(i,j,k))
!
!!        tmp1a=pom(i,j,k)-pom(i,j,k+1)+pom(i+1,j,k)-pom(i+1,j,k+1)
!        tmp1a=pom(i,j,k+1)-pom(i,j,k-1)
!        if(tmp1a .ne. 0.d0) then
!!          bcon(1)=-tmp1a/(2.*dx3(k)*detg(i,j,k))
!          bcon(1)=-tmp1a*(1./((x3(k+1)-x3(k-1))*detg(i,j,k))) 
!        else
!          bcon(1)=0.d0
!        endif
!!
!!        tmp1b=pom(i,j,k)+pom(i,j,k+1)-pom(i+1,j,k)-pom(i+1,j,k+1)
!        tmp1b=pom(i+1,j,k)-pom(i-1,j,k)
!        if(tmp1b .ne. 0.d0) then
!!          bcon(3)=tmp1b/(2.*dx1(i)*detg(i,j,k))
!          bcon(3)=tmp1b*(1./((x1(i+1)-x1(i-1))*detg(i,j,k))) 
!        else
!          bcon(3)=0.d0
!        endif
!
!!        if(k .eq. ks1) then
!!          bcon(1)=(pom(i,j,k+1)-pom(i,j,k))/((x3(k+1)-x3(k))*detg(i,j,k))
!!        elseif(k .eq. ke1) then
!!          bcon(1)=(pom(i,j,k)-pom(i,j,k-1))/((x3(k)-x3(k-1))*detg(i,j,k))
!!        else         
!!          bcon(1)=(pom(i,j,k+1)-pom(i,j,k-1))/((x3(k+1)-x3(k-1))*detg(i,j,k))
!!        endif
!
!!        if(i .eq. is1) then
!!          bcon(3)=-(pom(i+1,j,k)-pom(i,j,k))/((x1(i+1)-x1(i))*detg(i,j,k))
!!        elseif(i .eq. ie1) then
!!          bcon(3)=-(pom(i,j,k)-pom(i-1,j,k))/((x1(i)-x1(i-1))*detg(i,j,k))
!!        else         
!!          bcon(3)=-(pom(i+1,j,k)-pom(i-1,j,k))/((x1(i+1)-x1(i-1))*detg(i,j,k))
! !       endif
!
!        uri(7,i,j,k)=bcon(1)
!        uri(9,i,j,k)=bcon(3)
!
!      enddo
!    enddo
!  enddo
!
!- put boundary data for magnetic field -!
!
!  do j=js1,je1
!    do i=is1,ie1
!      uri(7,i,j,ks1)=uri(7,i,j,ks1+1)
!      uri(7,i,j,ke1)=uri(7,i,j,ke1-1)
!
!      uri(9,i,j,ks1)=uri(9,i,j,ks1+1)
!      uri(9,i,j,ke1)=uri(9,i,j,ke1-1)
!    enddo
!  enddo
!
!  do k=ks1,ke1
!    do j=js1,je1
!      uri(7,is1,j,k)=uri(7,is1+1,j,k)
!      uri(7,ie1,j,k)=uri(7,ie1-1,j,k)
!
!      uri(9,is1,j,k)=uri(9,is1+1,j,k)
!      uri(9,ie1,j,k)=uri(9,ie1-1,j,k)
!    enddo
!  enddo
!
!- calculation maximum magnetic field strength -!
!  do k=ks1,ke1
!    do j=js1,je1
!      do i=is1,ie1
!
!        util(1)=uri(2,i,j,k)
!        util(2)=uri(3,i,j,k)
!        util(3)=uri(4,i,j,k)
!        bcon(1)=uri(7,i,j,k)
!        bcon(2)=uri(8,i,j,k)
!        bcon(3)=uri(9,i,j,k)
!       
!        do m=0,3
!          do n=0,3
!            gcov1(m,n)=gcov(m,n,i,j,k)
!          enddo
!        enddo
!
!        call calgfl(util,gfl,gcov1)
!        call calucon(util,ucon,gfl,gcon1)
!        call lower(ucon,ucov,gcov1)
!        call calbbcon(bcon,bbcon,ucov,ucon)
!        call lower(bbcon,bbcov,gcov1)
!        call calbbsq(bbcon,bbcov,bbsq)
!
!        if(bbsq .gt. bsq_max) then
!          bsq_max=bbsq
!        endif
!
!      enddo
!    enddo
!  enddo    
!
!  call mpi_barrier(mpi_comm_world,merr)
!  call mpi_allreduce(bsq_max,bsq_max1,1,mpi_double_precision,mpi_max, &
!                     mpi_comm_world,merr)
!
!  if(bsq_max1 .eq. 0.d0) then
!    beta_act=0.d0 
!  else   
!    beta_act=prmax2a*(1./(0.5*bsq_max1))
!  endif   
!
!  if(beta .eq. 0.d0) then
!    bnorm=0.d0
!  else
!    bnorm=sqrt(beta_act*(1./beta))
!  endif
!
!  if(myrank .eq. 0) then
!    write(*,*) "beta_act, bnorm=", beta_act, bnorm
!  endif
!
!- Normalize magnetic field strength -!
!
!  do k=ks1,ke1
!    do j=js1,je1
!      do i=is1,ie1
!!
!        bcon(1)=uri(7,i,j,k)*bnorm            
!        bcon(2)=uri(8,i,j,k)*bnorm  
!        bcon(3)=uri(9,i,j,k)*bnorm  
!
!        util(1)=uri(2,i,j,k)
!        util(2)=uri(3,i,j,k)
!        util(3)=uri(4,i,j,k)
!
!        uri(7,i,j,k)=bcon(1)
!        uri(8,i,j,k)=bcon(2)
!        uri(9,i,j,k)=bcon(3)
!
!        do m=0,3
!          do n=0,3
!            gcov1(m,n)=gcov(m,n,i,j,k)
!          enddo
!        enddo        
!
!        call calgfl(util,gfl,gcov1)
!        call calucon(util,ucon,gfl,gcon1)
!        call lower(ucon,ucov,gcov1)
!        call calbbcon(bcon,bbcon,ucov,ucon)
!        call lower(bbcon,bbcov,gcov1)
!        call calbbsq(bbcon,bbcov,bbsq)
!
!        if(bbsq .gt. bsq_max2) then
!          bsq_max2=bbsq
!        endif
!
!      enddo
!    enddo
!  enddo 
!
!  call mpi_barrier(mpi_comm_world,merr)
!  call mpi_allreduce(bsq_max2,bsq_max2a,1,mpi_double_precision,mpi_max, &
!                     mpi_comm_world,merr)
!  if(bsq_max2a .eq. 0.d0) then
!    beta_act=0.d0 
!  else   
!    beta_act=prmax2a*(1./(0.5*bsq_max2a))
!   endif  
!
!  if(myrank .eq. 0) then
!    write(*,*) "beta_act=", beta_act
!  endif    
!
  deallocate(gcov1, gcon1, gcov1in, util, utiln, bcon, ucov, ucon, &
             bbcov, bbcon, beta1, ah, pom, stat=merr)
!
  return
end subroutine mdrecoilbh1a
!
!--------------------------------------------------------------------
subroutine mdrecoilbh2(uri,gcov,gcon, &
                       x1,x2,x3,x1b,x2b,x3b,dx1,dx2,dx3, &
                       is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
!  constant-l rotating torus model with toroidal magnetic field 
!  Komissarov 2006, MNRAS, 368, 993 + recoiling velocity
!  model=13  
!
  use pram, only : imax, jmax, kmax, nv, gam, metric, akm, pi, iter, kpol
  implicit none
  include 'mpif.h'

  integer :: i, j, k, m, n, is1 ,ie1, js1, je1, ks1, ke1
  integer :: merr

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  
  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:), gcov1cen(:,:), gcov1in(:,:)
  real(8), allocatable :: util(:), utiln(:), bcon(:), bbcon(:)
  real(8), allocatable :: ucov(:), ucon(:), beta1(:)
!
  real(8) :: x1(imax), x2(jmax), x3(kmax) 
  real(8) :: x1b(imax), x2b(jmax), x3b(kmax) 
  real(8) :: dx1(imax), dx2(jmax), dx3(kmax)
!
  real(8) :: de, pr, pm, ut, wp, ah, hm1
  real(8) :: rin, thin, utin, win,win1,  &
             rcen, thcen, ahcen, utcen, wcen, wcen1, clm, betacen, &
             rhocen, thedge, thedge1, rin1, rhomin, prmin, polkp, polkm
  real(8) :: rr, ph, th, sth, cth, del, aa, sig, &
             sthin, cthin, delin, aain, sigin, & 
             sthcen, cthcen, delcen, aacen, sigcen 
  real(8) :: bb, aa1, bbin, aa1in, bbcen, aa1cen, omg, uut, ucovph
  real(8) :: tmp1a, tmp1b, tmp1c, tmp1d, tmp_rand
  real(8) :: x1aa, x2aa, x3aa, x3t
  real(8) :: tmp4a, tmp4b, tmp4c, discr, alpha1, gfl
  real(8) :: rmb, rms, zz1, zz2
  real(8) :: rhomax1, rhomax1a, prmax1, prmax1a, prmax2
  real(8) :: small, vrecoil
! 
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), gcov1in(0:3,0:3), &
           gcov1cen(0:3,0:3), stat=merr)
  allocate(util(1:3), utiln(1:3), bcon(1:3), bbcon(0:3), &
           ucov(0:3), ucon(0:3), beta1(1:3), stat=merr)

!--------------------------------------------------------------------
!- Parameter -!
!- we assume polytropic index for gas pressure and mag pressure is same -!  

  zz1=1.+((1.-akm*akm)**(1./3.))*(((1.+akm)**(1./3.))+((1.-akm)**(1./3.)))
  zz2=sqrt((3.*akm*akm+zz1*zz1))
  rmb=2.*(1.-0.5*akm+sqrt(1.-akm))
  rms=3.+zz2-sqrt((3.-zz1)*(3.+zz1+2.*zz2))

  write(*,*) "rmb, rms=", rmb, rms
  
!  rin=rms     !- radial position at inner radius of torus
  rin=40.0d0
  thin=pi/2.d0  !- \theta position at center of torus  
  win1=-0.053d0    !- W at surface of torus
!  rcen=3.40d0    !- radial position at center of torus
  rcen=70.d0
  thcen=pi/2.d0  !- \theta position at center of torus
  ahcen=1.5d0     !- specific enthalpy at center of torus
  wcen1=-0.136d0   !- W at center of torus (calculated from rcen)
  rhocen=1.d0     !- density at center of torus
!  clm=2.6d0      !- angular momentum
  clm=8.0d0
  !  clm=4.5d0
  !  betacen=1.0d0  !- plasma beta at center of torus (\beta=p_gas/p_mag)
  betacen=1.d0
  !
  vrecoil=1.d-3
  !
  thedge=0.1*pi
  thedge1=0.9*pi
!  
  rin1=1.d0
  rhomin=1.d-5    !- minimum density at atmosphere (corona) 
  prmin=1.d-3*rhomin !-minimum pressure at atmosphere (corona)
!
  rhomax1=0.d0
  rhomax1a=0.d0
  prmax1=0.d0 
  prmax1a=0.d0
  prmax2=0.d0
!  
  small=1.d-12
!
  
!- initally set minimum plasma everywhere -!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1

        do m=0,3
          do n=0,3
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo
!        
       alpha1=1./sqrt(-gcon1(0,0))
!
        beta1(1)=alpha1*alpha1*gcon1(0,1)
        beta1(2)=alpha1*alpha1*gcon1(0,2)
        beta1(3)=alpha1*alpha1*gcon1(0,3)
        
        rr=x1(i)

        de=rhomin
!        de=rhomin*exp(-3.*rr/rcen)
!        de=rhomin*((rr/rin1)**(-1.5d0))

!        pr=(gam-1.d0)*prmin*((rr/rin1)**(-2.5d0))
!        pr=kpol*(de**gam)
!
        pr=prmin
        util(1)=0.d0
!        util(1)=(gcon1(0,1)/gcon1(0,0))*(1.-(1./rr)**4)
        util(2)=0.d0
        util(3)=0.d0
        bcon(1)=0.d0
        bcon(2)=0.d0 
        bcon(3)=0.d0
!          
!- set primitive variables -!
!          
        uri(1,i,j,k)=de
        uri(2,i,j,k)=util(1)
        uri(3,i,j,k)=util(2)
        uri(4,i,j,k)=util(3)
        uri(5,i,j,k)=pr
        uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
        uri(7,i,j,k)=bcon(1)
        uri(8,i,j,k)=bcon(2)          
        uri(9,i,j,k)=bcon(3)
!       
      enddo
    enddo
  enddo
!  
!- calculation of BL Kerr metric at rin, thin -!
!
!  sthin=sin(thin)
  sthin=1.d0
!  cthin=cos(thin)
  cthin=0.d0
  delin=rin*rin-2.*rin+akm*akm
  aain=(rin*rin+akm*akm)*(rin*rin+akm*akm)-delin*akm*akm*sthin*sthin
  sigin=rin*rin+akm*akm*cthin*cthin  
!
  do n=0,3
    do m=0,3
      gcov1in(m,n)=0.d0
    enddo
  enddo

  gcov1in(0,0)=-(sigin-2.*rin)*(1./sigin)               !- g_tt -!
  gcov1in(2,2)=(aain*(1./sigin))*sthin*sthin !- g_phi phi -!
  gcov1in(0,2)=-(2.*akm*rin*sthin*sthin)*(1./sigin)  !- g_t phi -!  
!          
!- calculation of u_t at r_in,theta_in -!
!          
  bbin=gcov1in(0,2)*gcov1in(0,2)-gcov1in(0,0)*gcov1in(2,2)
  aa1in=gcov1in(0,0)*clm*clm+2.*clm*gcov1in(0,2)+gcov1in(2,2)
          if(aa1in .gt. 0.d0) then
            utin=-sqrt(bbin*(1./aa1in))
          else
            utin=-1.d0  
          endif        
!
!- calculation of W from u_t at disc center
!          
          win=log(abs(utin))
          write(*,*) 'win=', win           
!  
!- calculation of BL Kerr metric at rcen, thcen -!
! 
!  sthcen=sin(thcen)
  sthcen=1.d0
!  cthcen=cos(thcen)
  cthcen=0.d0  
  delcen=rcen*rcen-2.*rcen+akm*akm
  aacen=(rcen*rcen+akm*akm)*(rcen*rcen+akm*akm)-delcen*akm*akm*sthcen*sthcen
  sigcen=rcen*rcen+akm*akm*cthcen*cthcen  
!
  do n=0,3
    do m=0,3
      gcov1cen(m,n)=0.d0
    enddo
  enddo

  gcov1cen(0,0)=-(sigcen-2.*rcen)*(1./sigcen)               !- g_tt -!
  gcov1cen(2,2)=(aacen*(1./sigcen))*sthcen*sthcen !- g_phi phi -!
  gcov1cen(0,2)=-(2.*akm*rcen*sthcen*sthcen)*(1./sigcen)  !- g_t phi -!
!          
!- calculation of u_t at r_cen,theta_cen -!
!          
  bbcen=gcov1cen(0,2)*gcov1cen(0,2)-gcov1cen(0,0)*gcov1cen(2,2)
  aa1cen=gcov1cen(0,0)*clm*clm+2.*clm*gcov1cen(0,2)+gcov1cen(2,2)
          if(aa1cen .gt. 0.d0) then
            utcen=-sqrt(bbcen*(1./aa1cen))
          else
            utcen=-1.d0  
          endif        
!
!- calculation of W from u_t at disc center
!          
          wcen=log(abs(utcen))
!
          write(*,*) 'wcen=', wcen          
!
!- calculate plokp and polkm
!- (polytropic constant for gas and mag pressure)
!             

!  tmp1c=exp((win-wcen)/(1.+1./betacen))-1.        
!  tmp1d=rhocen**(gam-1.)
!  polkp=((gam-1.)/gam)*tmp1c/tmp1d
!        
  if(betacen .eq. 0.d0) then
    polkp=1.d-3
    polkm=0.d0     
  else   
    tmp1c=(betacen*(1./(1.+betacen)))        
    polkp=(win-wcen)*((gam-1.)/gam)*tmp1c*(ahcen**(1.-gam))
    polkm=(bbcen**(1.-gam))*polkp*(1./betacen)    
  endif  
  
  write(*,*) 'polkp, polkm=', polkp, polkm
  
!
!- calculate w at each grid point -!
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!          
!- initialize -!
!         
        wp=1.d0
!
        rr=x1(i)
        th=x3(k)
        
        if(x2(j) .lt. 0.d0) then
          ph=2.*pi-x2(j)
        elseif(x2(j) .gt. 2.*pi) then
          ph=x2(j)-2.*pi
        else   
          ph=x2(j) 
        endif
!       
        sth=sin(th)
        if(abs(sth) .lt. small) then
          if(sth .ge. 0.d0) then
            sth=small
          elseif(sth .lt. 0.d0) then
            sth=-small 
          endif   
        endif        
!
        cth=cos(th) 
!
!        if(rr .ge. rin) then         
        if(rr .ge. rin .and. th .ge. thedge .and. th .le. thedge1) then
!- calculation of BL Kerr metric -! 
          del=(rr*rr)-(2.*rr)+(akm*akm)
          aa=(rr*rr+akm*akm)*(rr*rr+akm*akm)-del*akm*akm*sth*sth
          sig=rr*rr+akm*akm*cth*cth
!
          do n=0,3
            do m=0,3
              gcov1(m,n)=0.d0
              gcon1(m,n)=0.d0 
            enddo
          enddo
!
!- Kerr metric in BL coordinates -!
          gcov1(0,0)=-(sig-2.*rr)*(1./sig)               !- g^tt -!
          gcov1(1,1)=sig*(1./del)                          !- g^rr -!
          gcov1(2,2)=(aa*(1./sig))*sth*sth             !- g^phi phi -!
          gcov1(0,2)=-(2.*akm*rr*sth*sth)*(1./sig)           !- g^t phi -!
          gcov1(2,0)=gcov1(0,2)                  !- g^phi t -!
          gcov1(3,3)=sig                      !- g^th th -!        
!
          gcon1(0,0)=-aa*(1./(sig*del))               !- g^tt -!
          gcon1(1,1)=del*(1./sig)                          !- g^rr -!
          gcon1(2,2)=(sig-2.*rr)*(1./(del*sig*sth*sth)) !- g^phi phi -!
          gcon1(0,2)=-(2.*akm*rr)*(1./(sig*del))           !- g^t phi -!
          gcon1(2,0)=gcon1(0,2)                  !- g^phi t -!
          gcon1(3,3)=1./sig 
!
!- Kerr metric in KS coordinates -!
!          gcov1(0,0)=-1.+2.*rr/sig      !- g_tt -!
!          gcov1(0,1)=2.*rr/sig           !- g_tr -!
!          gcov1(0,2)=-2.*akm*rr*sth*sth/sig !- g_t phi -!
!          gcov1(1,0)=gcov1(0,1)        !- g_rt -!
!          gcov1(1,1)=1.+2.*rr/sig      !- g_rr -!
!          gcov1(1,2)=-1.*akm*sth*sth*(1.+2.*rr/sig) !- g_r phi -! 
!          gcov1(2,0)=gcov1(0,2)        !- g_phi t -!
!          gcov1(2,1)=gcov1(1,2)        !- g_phi r -!
!          gcov1(2,2)=aa*sth*sth/sig         !- g_phi phi -!
!          gcov1(3,3)=sig                !- g_th th -!
!
!          gcon1(0,0)=-1.-2.*rr/sig      !- g^tt -!
!          gcon1(0,1)=2.*rr/sig           !- g^tr -!
!          gcon1(1,0)=gcon1(0,1)        !- g^rt -!
!          gcon1(1,1)=del/sig           !- g^rr -!
!          gcon1(1,2)=akm/sig            !- g^r phi -!
!          gcon1(2,1)=gcon1(1,2)        !- g^phi r -!
!          gcon1(2,2)=1./(sig*sth*sth)         !- g^phi phi -!
!          gcon1(3,3)=1./sig                !- g^th th -!
!          
!- calculation of u_t at r,theta -!
!          
          bb=gcov1(0,2)*gcov1(0,2)-gcov1(0,0)*gcov1(2,2)
          aa1=gcov1(0,0)*clm*clm+2.*clm*gcov1(0,2)+gcov1(2,2)
          if(aa1 .gt. 0.d0) then
            ut=-sqrt(bb*(1./aa1))
          else
            ut=-1.d0  
          endif        
!
!- calculation of W from u_t
!          
          wp=log(abs(ut))
!
!- calculation of h
!
          if(wp .le. win) then
!            tmp1a=(wp-win)/(wcen-win)
!            tmp1b=(polkp+polkm*bbcen**(gam-1.))/(polkp+polkm*bb**(gam-1.))
!            ah=ahcen*(tmp1a*tmp1b)**(1./(gam-1.))          
!
!- from Komissarov 
            tmp1a=(win-wp)
!            tmp1b=1.+(polkm/polkp)*bb**(gam-1.)
            tmp1b=(gam/(gam-1.))*(polkp+polkm*bb**(gam-1.))
!            tmp1b=1.d0
            if(tmp1a .gt. 0.d0) then
!              ah=exp(tmp1a/tmp1b)            
              ah=(tmp1a*(1./tmp1b))**(1./(gam-1.))
            else
              ah=-1.d0
            endif   
          else
            ah=-1.d0
          endif            
!          
        else
          ah=-1.d0
        endif
!
!- set torus structure -!
!  
!        if(ah .gt. 1.d-4 .and. rr .ge. rin) then  &
!           .and. th .ge. thedge .and. th .le. thedge1) then    
        if(ah .gt. 0.d0) then
!
!          hm1=max((ah-1.)*(gam-1.)/gam, 0.d0) 
!          de=(hm1/polkp)**(1./(gam-1.)) 
!          pr=polkp*(hm1/polkp)**(gam/(gam-1.))
!          pm=polkm*(bb**(gam-1.))*(pr/polkp)
!
          pr=polkp*(ah**gam)
          pm=polkm*(bb**(gam-1.))*(ah**gam)
!          de=ah-pr
          de=(pr/polkp)**(1./gam)

          if(de .lt. rhomin) then
            de=rhomin
          endif
          if(pr .lt. prmin) then
            pr=prmin
         endif
          
          ucon(1)=0.d0
          ucon(3)=0.d0
!
!- calculation of omega & contravariant u^t
!          
          omg=-(gcov1(0,2)+gcov1(0,0)*clm)*(1./(gcov1(2,2)+gcov1(0,2)*clm))
          uut=sqrt(-1./(gcov1(0,0)+2.*omg*gcov1(0,2)+omg*omg*gcov1(2,2)))
!          uut=-1./(ut*(1.-omg*clm))
!
          ucovph=-omg*ut

          ucon(0)=uut
          ucon(2)=(omg*uut)
!          ucon(2)=gcon1(0,2)*ut+gcon1(2,2)*ucovph
!          ucon(2)=0.d0
!
!
          if(betacen .eq. 0.d0) then
            bbcon(1)=0.d0
            bbcon(2)=0.d0
            bbcon(3)=0.d0
            bbcon(0)=0.d0
          else   
            bbcon(1)=0.d0
            if(aa1 .gt. 0.d0) then 
              bbcon(2)=sqrt(2.*pm*(1./aa1))
            else
              bbcon(2)=0.d0
            endif   
            bbcon(3)=0.d0            
            bbcon(0)=clm*bbcon(2)
!            
!            bcon(1)=0.d0
!            bcon(2)=sqrt(2.*pm)
!            bcon(3)=0.d0
          endif
!
!
!- cal 3-magnetic field -!
!          
          bcon(1)=ucon(0)*bbcon(1)-bbcon(0)*ucon(1)
          bcon(2)=ucon(0)*bbcon(2)-bbcon(0)*ucon(2)
          bcon(3)=ucon(0)*bbcon(3)-bbcon(0)*ucon(3)
!
!          bcon(1)=gfl*bbcon(1)-bbcon(0)*ucon(1)
!          bcon(2)=gfl*bbcon(2)-bbcon(0)*ucon(2)
!          bcon(3)=gfl*bbcon(3)-bbcon(0)*ucon(3)
!
!- cal spacial 4-velocity -!
!
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
          gfl=ucon(0)*alpha1

          util(1)=ucon(1)+gfl*beta1(1)*(1./alpha1)
          util(2)=ucon(2)+gfl*beta1(2)*(1./alpha1)
          util(3)=ucon(3)+gfl*beta1(3)*(1./alpha1)             
!
!- Put perturbation in pressure -!
!
!         call random_number(tmp_rand)
!!          write(*,*) "tmp_rand=", tmp_rand-0.5
!           pr=pr*(1.+0.05*(tmp_rand-0.5))
!           util(1)=util(1)+0.01*(tmp_rand-0.5)             
!
!- Put recoiling speed -!
!
          util(1)=util(1)+vrecoil*cos(ph)*sin(th)/sqrt(gcov1(1,1))
          util(2)=util(2)-vrecoil*sin(ph)/sqrt(gcov1(2,2))
          util(3)=util(3)+vrecoil*cos(ph)*cos(th)/sqrt(gcov1(3,3))
!
!- calculation of bcon -!
!
!          do m=0,3
!            do n=0,3
!              gcov1(m,n)=gcov(m,n,i,j,k)
!              gcon1(m,n)=gcon(m,n,i,j,k)
!            enddo
!          enddo
!
!          alpha1=1./sqrt(-gcon1(0,0))
!                
!          call calgfl(util,gfl,gcov1)
!          call calucon(util,ucon,gfl,gcon1)
!
!          bcon(1)=ucon(0)*bbcon(1)-bbcon(0)*ucon(1)
!          bcon(2)=ucon(0)*bbcon(2)-bbcon(0)*ucon(2)
!          bcon(3)=ucon(0)*bbcon(3)-bbcon(0)*ucon(3)
!          
!          bcon(1)=alpha1*ucon(0)*bbcon(1)-alpha1*bbcon(0)*ucon(1)
!          bcon(2)=alpha1*ucon(0)*bbcon(2)-alpha1*bbcon(0)*ucon(2)
!          bcon(3)=alpha1*ucon(0)*bbcon(3)-alpha1*bbcon(0)*ucon(3)
!          
!          if(j .eq. 4 .and. k .eq. 50) then
!            write(*,*) 'rr, B\phi=', x1(i), bcon(2)
!          endif
!          
!- check density & pressure max -!
!
          if(de .gt. rhomax1) then
            rhomax1=de
!            write(*,*) 'de,i,k,rr,th=', de,i,k,x1(i),x3(k)
          endif
          if(pr .gt. prmax1) then
            prmax1=pr
          endif
!         
!- set primitive variables -!
!          
          uri(1,i,j,k)=de
          uri(2,i,j,k)=util(1)
          uri(3,i,j,k)=util(2)
          uri(4,i,j,k)=util(3)
          uri(5,i,j,k)=pr
          uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
          uri(7,i,j,k)=bcon(1)
          uri(8,i,j,k)=bcon(2)
          uri(9,i,j,k)=bcon(3)
!          
        endif
!
      enddo
    enddo
  enddo
!
  deallocate(gcov1, gcon1, gcov1cen, gcov1in, util, utiln, bcon, bbcon, &
             ucov, ucon, beta1, stat=merr)
!
  return
end subroutine mdrecoilbh2
!
!--------------------------------------------------------------------
subroutine mdwindcol(uri,gcov,gcon,x1,x2,x3, &
                     myrank,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
!         model = 11; 2D Wind Collision Problem

  use pram, only : nv, imax, jmax, kmax, gam, c0, pi
  implicit none
  include 'mpif.h' ! for MPI
  
  integer :: i, j, k, n, m, is1, ie1, js1, je1, ks1, ke1
  integer :: myrank, merr 

  integer, parameter :: nkmax=8

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1) 
!
  real(8) :: x1(imax), x2(jmax), x3(kmax)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: vel(:), bcon(:), phi1(:)
! 
  real(8) :: rr, th, de, pr, eps, rsh, rsh1, &
            tmp1, tmp2, tmp3, vsq, gfl, gflmax
  real(8) :: gamfw, gamsw, rofw, rosw, roatm, afw, asw, aatm, &
            eta, epsfw, epssw, epsatm, bfw, bsw, batm, prconst
!
!-allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), stat=merr)
  allocate(vel(1:3), bcon(1:3), phi1(1:nkmax), stat=merr)  
  
!--------------------------------------------------------------------
!- Parameter -!
!
  rsh=500.d0 !- collision radius -!
  rsh1=1000.d0 !- slow wind head radius -!
  gflmax=10.d0

  gamfw=4./3.
  gamsw=5./3.

  eta=2.d0 ! density slope for atmopshere
  
  rofw=0.01d0 ! density of fast-wind at collision radius
  rosw=1.d0   ! density of slow-wind at collision radius
  roatm=1.d-6 ! density of atmosphere at slow-wind head radius
  
  afw=rofw*rsh**(2.)
  asw=rosw*rsh**(2.)
  aatm=roatm*rsh1**(eta)
  
  epsfw=1.d0 ! internal energy of fast-wind at collision radius
  epssw=((gamfw-1.)/(gamsw-1.))*(rofw/rosw)*epsfw
            ! internal energy of slow-wind at collision raidus
  epsatm=1.d-3 ! internal energy of atmpshere at slow-wind head radius
  
  bfw=epsfw*rsh**(2.*(gamfw-1.))
  bsw=epssw*rsh**(2.*(gamsw-1.))
  batm=epsatm*rsh1**(eta*(gamsw-1.))

  write(*,*) 'batm=', batm
  prconst=roatm*epsatm*(gamsw-1.)
  
  if(myrank .eq. 0) then
    do n=1, nkmax
      call random_number(tmp3)
      phi1(n)=2.0*pi*tmp3
!      phi1(n)=0.d0       
    enddo
  endif
!
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_bcast(phi1,nkmax,mpi_double_precision,0,mpi_comm_world,merr)    
!    
!
! Initial setup
!  
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!         
        rr=x1(i)
        th=x3(k)
!        
        if(rr .le. rsh) then
!- fast wind region
!          de=afw*rr**(-2.)
!          eps=bfw*rr**(-2.*(gamfw-1.)) 
          de=rofw
          eps=epsfw 
          pr=(gamfw-1.)*de*eps
          vel(1)=0.4d0
          vel(2)=0.d0
          vel(3)=0.d0
          bcon(1)=0.d0
          bcon(2)=0.d0
          bcon(3)=0.d0

        elseif(rr .gt. rsh .and. rr .le. rsh1) then          
!- slow wind region           
!          de=asw*rr**(-2.)
!          eps=bsw*rr**(-2.*(gamsw-1.))
          de=rosw
          eps=epssw 
          pr=(gamsw-1.)*de*eps
          vel(1)=0.1d0
          vel(2)=0.d0
          vel(3)=0.d0
          bcon(1)=0.d0
          bcon(2)=0.d0
          bcon(3)=0.d0
!
        else
!        elseif(rr .gt. rsh1) then
!- atmosphere region                                                            
          de=aatm*rr**(-eta)
!          eps=batm*rr*(-eta*(gamsw-1.)) 
!          eps=batm*rr*(-0.5*eta)
!          de=roatm
          eps=epsatm 
!          eps=prconst/((gamsw-1.)*de)
          pr=(gamsw-1.)*de*eps
!          pr=1.d-5*(de**gamsw)
          vel(1)=0.0d0
          vel(2)=0.d0
          vel(3)=0.d0
          bcon(1)=0.d0
          bcon(2)=0.d0
          bcon(3)=0.d0
!
        endif
!
!- perturbation
!       
!- random number -!
!        call random_number(tmp1)
!        tmp2=2.0*tmp1-1.0
!- sinusoidal -!
!        do n=1, nkmax
!          tmp2=tmp2+sin(float(n)*th+phi1(n))
!        enddo   
!        tmp2=tmp2/float(nkmax)  
!       
!        if(rr .gt. 0.9d0*rsh .and. rr .lt. 1.1d0*rsh) then
!          vel(1)=vel(1)+0.05*tmp2 
!        endif
!      
!- convert 3-vel => 4-vel
!
        vsq=0.d0
        do m=1,3
          do n=1,3
            vsq=vsq+gcov1(m,n)*vel(m)*vel(n) 
          enddo
        enddo
!
        if(vsq .ge. 1.d0) then
          vsq=1.-1./(gflmax*gflmax)
        endif   
!- new set of primitive variables -!
        gfl=1./sqrt(1.-vsq)
        uri(1,i,j,k)=de
        uri(2,i,j,k)=vel(1)*gfl
        uri(3,i,j,k)=vel(2)*gfl
        uri(4,i,j,k)=vel(3)*gfl
        uri(5,i,j,k)=pr
        uri(6,i,j,k)=pr*(1./(de**(gam-1.)))
        uri(7,i,j,k)=bcon(1)
        uri(8,i,j,k)=bcon(2)
        uri(9,i,j,k)=bcon(3)
!        
      enddo
    enddo
  enddo
!
  deallocate(gcov1, gcon1, vel, bcon, phi1, stat=merr)
!
  return
end subroutine mdwindcol
!
!***********************************************************************
!           GENERAL FUNCTIONS
!***********************************************************************
!--------------------------------------------------------------------
function ffjump(xx,bi,ee,aa,dd,bb)
!--------------------------------------------------------------------
  implicit none
  real(8) :: xx, bi, ee, aa, dd, bb, ffjump
!
  if( xx .le. aa-bb ) then
    ffjump=bi
   elseif( xx .ge. aa+bb ) then
    ffjump=ee
   else
    ffjump=1.0+tanh((xx-aa)/dd)/tanh(bb/dd)
    ffjump=bi+0.5*ffjump*(ee-bi)
   endif
!
   return
end function ffjump
!

