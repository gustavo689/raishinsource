!--------------------------------------------------------------------
subroutine calgfl(util,gfl,gcov1)
!--------------------------------------------------------------------
!- Calculation of Lorentz factor -!
!
  implicit none

  integer :: m, n
  real(8) :: gcov1(0:3,0:3) !- covariant metric -!
  real(8) :: util(1:3) !- \tilde{u}^i = \gamma v^i -!
  real(8) :: utsq, gfl
!
  utsq=0.d0
!  vsq=0.d0
!
!- g_ij (4-metric) =\gamma_ij (3-metric)
  do m=1,3
    do n=1,3
      utsq=utsq+gcov1(m,n)*util(m)*util(n)     
    enddo
  enddo
!
  if(utsq .lt. 0.d0) then
    if(abs(utsq) .le. 1.0d-10) then
      utsq=abs(utsq)
    else
      utsq=1.d-10 !- set floot number -!
    endif
  endif
  if(abs(utsq) .gt. 5.d2) then
    utsq=5.0d2
  endif
!
  gfl=sqrt(1.+utsq) !- Lorentz factor -!
!  gfl=1./sqrt(1.-vsq)
!
!  util(1)=gfl*vcon(1)
!  util(2)=gfl*vcon(2)
!  util(3)=gfl*vcon(3)
!
  return
end subroutine calgfl
!
!--------------------------------------------------------------------
subroutine calucon(util,ucon,gfl,gcon1)
!--------------------------------------------------------------------
!- calculation of contravariant 4-velocity -!
!
  implicit none

!- \tilde{u}^i = \gamma v^i = u^i + \gamma \beta^i / \alpha -! 
!- beta: shift vector, alpha: lapse function, gfl: Lorentz factor -!
!
  integer :: merr  
  real(8) :: util(1:3)
  real(8) :: gcon1(0:3,0:3) !- contravariant metric -!
  real(8) :: ucon(0:3) !- contravariant 4-velocity -!
  real(8), allocatable :: beta1(:)
  
  real(8) :: alpha1, gfl
! 
!- allocate variables -!
  allocate(beta1(1:3), stat=merr)  
!
  alpha1=1./sqrt(-gcon1(0,0))
  beta1(1)=alpha1*alpha1*gcon1(0,1)
  beta1(2)=alpha1*alpha1*gcon1(0,2)
  beta1(3)=alpha1*alpha1*gcon1(0,3)
!
  ucon(0)=gfl*(1./alpha1)
  ucon(1)=util(1)-gfl*beta1(1)*(1./alpha1)
  ucon(2)=util(2)-gfl*beta1(2)*(1./alpha1)
  ucon(3)=util(3)-gfl*beta1(3)*(1./alpha1)
!
  deallocate(beta1, stat=merr)
!  
  return
end subroutine calucon
!
!--------------------------------------------------------------------
subroutine calbbcon(bcon,bbcon,ucov,ucon)
!--------------------------------------------------------------------
!-  Calcualtion of contravariant 4-magnetic field -! 
!
  implicit none

!- bbcon: contravariant 4-magnetic field (output) -!
!- ucov: covariant 4-velocity, ucon: contravariant 4-velocity (input) -!
!- bcon: contravariant 3-magnetic field (input) -!
!
  real(8) :: bbcon(0:3), ucov(0:3), ucon(0:3)
  real(8) :: bcon(1:3)

  bbcon(0)=bcon(1)*ucov(1)+bcon(2)*ucov(2)+bcon(3)*ucov(3)
  bbcon(1)=(bcon(1)+bbcon(0)*ucon(1))*(1./ucon(0))
  bbcon(2)=(bcon(2)+bbcon(0)*ucon(2))*(1./ucon(0))
  bbcon(3)=(bcon(3)+bbcon(0)*ucon(3))*(1./ucon(0))

  return
end subroutine calbbcon
!
!--------------------------------------------------------------------
subroutine calbbsq(bbcon,bbcov,bbsq)
!--------------------------------------------------------------------
!- Calculation of square of 4-magnetic field -!
!
  implicit none

!- bbcon: contravariant 4-magnetic field -!
!- bbcov: covariant 4-magnetic field -!

  real(8) :: bbcon(0:3), bbcov(0:3)
  real(8) :: bbsq

  bbsq=bbcon(0)*bbcov(0)+bbcon(1)*bbcov(1) &
      +bbcon(2)*bbcov(2)+bbcon(3)*bbcov(3)

  return
end subroutine calbbsq
!
!--------------------------------------------------------------------
subroutine caldelta(delta,gcov1,gcon1)
!--------------------------------------------------------------------
!- Calculation of Kronecker delta function -!
!- \delta^\mu_\nu=g^{\mu \tau} *g_{\tau \nu} -!
!
  implicit none

  integer :: m, n
  real(8) :: gcov1(0:3,0:3), gcon1(0:3,0:3), delta(0:3,0:3)

  do m=0,3
    do n=0,3
      delta(m,n)=0.d0
    enddo
  enddo

!  do m=0,3
!    do n=0,3
!      delta(0,n)=gcon1(0,m)*gcov1(m,n) 
!      delta(1,n)=gcon1(1,m)*gcov1(m,n) 
!      delta(2,n)=gcon1(2,m)*gcov1(m,n) 
!      delta(3,n)=gcon1(3,m)*gcov1(m,n) 
!    enddo
!  enddo
!
  do m=0,3
    do n=0,3
      if(m .eq. n) then
        delta(m,n)=1.d0
      else
        delta(m,n)=0.d0
      endif
    enddo
  enddo
!
  return
end subroutine caldelta
!
!--------------------------------------------------------------------
subroutine lower(vcon,vcov,gcov1)
!--------------------------------------------------------------------
!- contravariant 4-vector => covariant 4-vector -!
!
  implicit none
!
!- vcon: contravariant 4-vector (input)-!
!- vcov: covariant 4-vector (output)-!
!- gcov1: covariant metric (input)-! 
!
  real(8) :: vcon(0:3), vcov(0:3) 
  real(8) :: gcov1(0:3,0:3)

  vcov(0)=gcov1(0,0)*vcon(0) &
         +gcov1(0,1)*vcon(1) &
         +gcov1(0,2)*vcon(2) &
         +gcov1(0,3)*vcon(3)

  vcov(1)=gcov1(1,0)*vcon(0) &
         +gcov1(1,1)*vcon(1) &
         +gcov1(1,2)*vcon(2) &
         +gcov1(1,3)*vcon(3)

  vcov(2)=gcov1(2,0)*vcon(0) &
         +gcov1(2,1)*vcon(1) &
         +gcov1(2,2)*vcon(2) &
         +gcov1(2,3)*vcon(3)

  vcov(3)=gcov1(3,0)*vcon(0) &
         +gcov1(3,1)*vcon(1) &
         +gcov1(3,2)*vcon(2) &
         +gcov1(3,3)*vcon(3)

  return
end subroutine lower
!
!--------------------------------------------------------------------
subroutine upper(vcov,vcon,gcon1)
!--------------------------------------------------------------------
!- covariant 4-vector => contravariant 4-vector -!
!
  implicit none
!
!- vcov: covariant 4-vector (input)-!
!- vcon: contravariant 4-vector (output)-!
!- gcon1: contravariant metric (input)-! 
!
  real(8) :: vcov(0:3), vcon(0:3) 
  real(8) :: gcon1(0:3,0:3)

  vcon(0)=gcon1(0,0)*vcov(0) &
         +gcon1(0,1)*vcov(1) &
         +gcon1(0,2)*vcov(2) &
         +gcon1(0,3)*vcov(3)

  vcon(1)=gcon1(1,0)*vcov(0) &
         +gcon1(1,1)*vcov(1) &
         +gcon1(1,2)*vcov(2) &
         +gcon1(1,3)*vcov(3)

  vcon(2)=gcon1(2,0)*vcov(0) &
         +gcon1(2,1)*vcov(1) &
         +gcon1(2,2)*vcov(2) &
         +gcon1(2,3)*vcov(3)

  vcon(3)=gcon1(3,0)*vcov(0) &
         +gcon1(3,1)*vcov(1) &
         +gcon1(3,2)*vcov(2) &
         +gcon1(3,3)*vcov(3)

  return
end subroutine upper
!
!--------------------------------------------------------------------
subroutine lower1(vcon,vcov,gcov3a)
!--------------------------------------------------------------------
!- contravariant 4-vector => covariant 4-vector -!
!
  implicit none
!
!- vcon: contravariant 3-vector (input) -!
!- vcov: covariant 3-vector (output) -!
!- gcov3a: covariant 3-metric (input) -! 
!
  real(8) :: vcon(1:3), vcov(1:3) 
  real(8) :: gcov3a(1:3,1:3)


  vcov(1)=gcov3a(1,1)*vcon(1) &
         +gcov3a(1,2)*vcon(2) &
         +gcov3a(1,3)*vcon(3)

  vcov(2)=gcov3a(2,1)*vcon(1) &
         +gcov3a(2,2)*vcon(2) &
         +gcov3a(2,3)*vcon(3)

  vcov(3)=gcov3a(3,1)*vcon(1) &
         +gcov3a(3,2)*vcon(2) &
         +gcov3a(3,3)*vcon(3)

  return
end subroutine lower1
!
!--------------------------------------------------------------------
subroutine upper1(vcov,vcon,gcon3a)
!--------------------------------------------------------------------
!- covariant 3-vector => contravariant 3-vector -!
!
  implicit none
!
!- vcov: covariant 3-vector (input) -!
!- vcon: contravariant 3-vector (output) -!
!- gcon3a: contravariant 3-metric (input) -! 
!
  real(8) :: vcov(1:3), vcon(1:3) 
  real(8) :: gcon3a(1:3,1:3)

  vcon(1)=gcon3a(1,1)*vcov(1) &
         +gcon3a(1,2)*vcov(2) &
         +gcon3a(1,3)*vcov(3)

  vcon(2)=gcon3a(2,1)*vcov(1) &
         +gcon3a(2,2)*vcov(2) &
         +gcon3a(2,3)*vcov(3)

  vcon(3)=gcon3a(3,1)*vcov(1) &
         +gcon3a(3,2)*vcov(2) &
         +gcon3a(3,3)*vcov(3)

  return
end subroutine upper1
!
!--------------------------------------------------------------------
subroutine cal4to3vel(uri,uria,gcov,x1,x2,x3,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
!- Calculation of 3-velocity from 4-velocity -!
!- uri: primitive variables with \tilde{u} (input)
!- urib: primiritve variables with 3-velocity v^i (output)  
!
  use pram, only : nv, imax, jmax, kmax  
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1 
  integer :: merr, jh, kh

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1) ! prim variables with 4-vel -!
  real(8) :: uria(nv,is1:ie1,js1:je1,ks1:ke1) ! prim variables with 3-vel -!
  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: x1(imax),x2(jmax),x3(kmax)
  
  real(8), allocatable :: gcov1(:,:) !- covariant metric -!
  real(8), allocatable :: util(:) !- \tilde{u}^i = \gamma v^i -!
  real(8) :: gfl
!
  allocate(gcov1(0:3,0:3), util(1:3), stat=merr)

!=====================================================================
!
  jh=(jmax+1)/2
  kh=((kmax+1)/2)  

!  write(*,*) "x3=",x3(kh)
  
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
         
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
!
!
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
          enddo
        enddo
!        
        call calgfl(util,gfl,gcov1)
!
        uria(1,i,j,k)=uri(1,i,j,k)
        uria(2,i,j,k)=util(1)*(1./gfl)
        uria(3,i,j,k)=util(2)*(1./gfl)
        uria(4,i,j,k)=util(3)*(1./gfl)
        uria(5,i,j,k)=uri(5,i,j,k)
        uria(6,i,j,k)=uri(6,i,j,k)
        uria(7,i,j,k)=uri(7,i,j,k)
        uria(8,i,j,k)=uri(8,i,j,k)
        uria(9,i,j,k)=uri(9,i,j,k)
!
!        if(j .eq. jh .and. k .eq. kh) then
!          write(*,*) "r, v^\phi=", x1(i), uria(3,i,j,k) 
!        endif
!        
      enddo
    enddo
  enddo
!
  deallocate(gcov1, util, stat=merr)
!  
  return
end subroutine cal4to3vel
!
!--------------------------------------------------------------------
subroutine transbl2ks(util,utiln,x1aa,x2aa,x3aa,x3t)
!--------------------------------------------------------------------
!- convert BL 4-velocity to 4-velocity in Kerr-Schild coordinates -!
!
!- util: \tilde{u}_BL (input)-!
!- utiln: \tilde{u}_KS (output)-! 
!
!
  use pram, only : metric, akm, ix1, ix2, ix3, R0, hslope, pi
  implicit none

  integer :: m, n, merr
! 
  real(8) :: util(1:3), utiln(1:3) 
!
  real(8), allocatable :: ucon(:), ucon1(:)
  real(8), allocatable :: gcov1(:,:), gcon1(:,:), beta1(:)
  real(8), allocatable :: tmp1(:), trans(:,:)
!
  real(8) :: x1aa, x2aa, x3aa, rr, aa, bb, cc, discr, alpha1, gfl, x3t, &
             tmp2, rbh1
!

  allocate(ucon(0:3), ucon1(0:3), gcov1(0:3,0:3), gcon1(0:3,0:3), &
           beta1(1:3), tmp1(0:3), trans(0:3,0:3), stat=merr)
!
  rbh1=1.+sqrt(1.-akm*akm)  
  rr=x1aa

  if(rr .gt. rbh1) then
!
!- get BL kerr metric (covariant metric) -!
    call kercov1(gcov1,x1aa,x2aa,x3aa)
!
!- contravariant velocity in BL coordinates -!
    ucon(0)=0.d0
    ucon(1)=util(1)
    ucon(2)=util(2)
    ucon(3)=util(3)
!
    aa=gcov1(0,0)
    bb=2.*(gcov1(0,1)*ucon(1)+gcov1(0,2)*ucon(2)+gcov1(0,3)*ucon(3))
    cc=1.+gcov1(1,1)*ucon(1)*ucon(1)+gcov1(2,2)*ucon(2)*ucon(2) &    
          +gcov1(3,3)*ucon(3)*ucon(3) &
          +2.*(gcov1(1,2)*ucon(1)*ucon(2)+gcov1(1,3)*ucon(1)*ucon(3) &  
          +gcov1(2,3)*ucon(2)*ucon(3))   
    discr=bb*bb-4.*aa*cc
    ucon(0)=(-bb-sqrt(discr))*(1./(2.*aa))
!
!- make transform matrix -!
    do m=0,3
      do n=0,3
        trans(m,n)=0.d0
      enddo
    enddo

    tmp2=rr*rr -(2.*rr) +akm*akm

    if(tmp2 .gt. 0.d0) then
      trans(0,0)=1.d0
      trans(0,1)=2.*rr*(1./tmp2)
      trans(1,1)=1.d0 
      trans(2,2)=1.d0
      trans(2,1)=akm*(1./tmp2)
      trans(3,3)=1.d0
    endif  
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
!- KS => modified KS coords -!
    if(metric .eq. 403) then
      if(ix1 .eq. 2) then
        ucon1(1)=ucon1(1)*(1./(rr))
      endif
      if(ix3 .eq. 2) then
        ucon1(3)=ucon1(3)*(1./(1.+(1.-hslope)*cos(2.*x3t)))
      endif
    endif
!
!- get Kerr-Schild metric (contravariant metric) -!
    if(metric .eq. 303) then
      call kshcon1(gcon1,x1aa,x2aa,x3aa)
    elseif(metric .eq. 403) then
      call mkshcov1(gcov1,x1aa,x2aa,x3aa,x3t)
      call invert_matrix(gcov1,gcon1,4)
    endif
! 
    alpha1=1./sqrt(-gcon1(0,0))
    gfl=ucon1(0)*alpha1
!
    beta1(1)=alpha1*alpha1*gcon1(0,1)
    beta1(2)=alpha1*alpha1*gcon1(0,2)
    beta1(3)=alpha1*alpha1*gcon1(0,3)
!
    utiln(1)=ucon1(1)+gfl*beta1(1)*(1./alpha1)
    utiln(2)=ucon1(2)+gfl*beta1(2)*(1./alpha1)
    utiln(3)=ucon1(3)+gfl*beta1(3)*(1./alpha1)

  else
!
!- get Kerr-Schild metric (contravariant metric) -!
!
    if(metric .eq. 303) then
      call kshcon1(gcon1,x1aa,x2aa,x3aa)
    elseif(metric .eq. 403) then
      call mkshcov1(gcov1,x1aa,x2aa,x3aa,x3t)
      call invert_matrix(gcov1,gcon1,4)
    endif
! 
    alpha1=1./sqrt(-gcon1(0,0))
!
    beta1(1)=alpha1*alpha1*gcon1(0,1)
    beta1(2)=alpha1*alpha1*gcon1(0,2)
    beta1(3)=alpha1*alpha1*gcon1(0,3)
!
    utiln(1)=beta1(1)*(1./alpha1)
    utiln(2)=beta1(2)*(1./alpha1)
    utiln(3)=beta1(3)*(1./alpha1)
!
  endif
!  
  deallocate(ucon, ucon1, gcov1, gcon1, beta1, tmp1, trans, stat=merr)
!
  return
end subroutine transbl2ks
!--------------------------------------------------------------------
subroutine transks2bl(uri,uria,gcov,gcon,x1,x2,x3,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
!- convert 4-velocity in KS to 4-velocity in BL coordinates -!
!-  reading whole primitive variables -!
!
!- uri: primitive variables with 4-velocity in KS (input)-!
!- uria: primitive variables with 4-velocity in BL (output)-! 
!
!  
  use pram, only : imax, jmax, kmax, nv, akm  
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1 
  integer :: merr
       
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uria(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8), allocatable :: gcov1(:,:), gcon1(:,:), gcon1a(:,:), beta1(:)
  real(8), allocatable :: util(:), utiln(:)
  real(8), allocatable :: ucon(:), ucon1(:)
  real(8), allocatable :: tmp1(:), trans(:,:)

  real(8) :: x1(imax), x2(jmax), x3(kmax)      

  real(8) :: x1aa, x2aa, x3aa, rr, rbh1, alpha1, gfl, tmp2
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), gcon1a(0:3,0:3), beta1(1:3), &
           util(1:3), utiln(1:3), ucon(0:3), ucon1(0:3), &
           tmp1(0:3), trans(0:3,0:3), stat=merr)
!
!=====================================================================
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1

        rbh1=1.+sqrt(1.-akm*akm)
        rr=x1(i)
!
        if(rr .gt. rbh1) then
          util(1)=uri(2,i,j,k)
          util(2)=uri(3,i,j,k)
          util(3)=uri(4,i,j,k)
!
!          v1=uri(2,i,j,k)
!          v2=uri(3,i,j,k)
!          v3=uri(4,i,j,k)

!
          do m=0,3
            do n=0,3
              gcov1(m,n)=gcov(m,n,i,j,k)
              gcon1(m,n)=gcon(m,n,i,j,k)
            enddo
          enddo
!
!- calculation of contravariant 4-velocity in KS coordinates -!        
          call calgfl(util,gfl,gcov1)
          call calucon(util,ucon,gfl,gcon1)
!
!- make transform matrix -!
        !- make transform matrix -!
          do m=0,3
            do n=0,3
              trans(m,n)=0.d0
            enddo
          enddo

          tmp2=rr*rr -(2.*rr) +akm*akm
          
          if(tmp2 .gt. 0.d0) then
            trans(0,0)=1.d0
            trans(0,1)=-2.*rr*(1./tmp2)
            trans(1,1)=1.d0 
            trans(2,2)=1.d0
            trans(2,1)=-akm*(1./tmp2)
            trans(3,3)=1.d0
          endif  
!
!- initialize -!
          do m=0,3
            tmp1(m)=0.d0
            ucon1(m)=0.d0
          enddo
          do n=0,3
            tmp1(0)=tmp1(0)+trans(0,n)*ucon(n)
            tmp1(1)=tmp1(1)+trans(1,n)*ucon(n)
            tmp1(2)=tmp1(2)+trans(2,n)*ucon(n)
            tmp1(3)=tmp1(3)+trans(3,n)*ucon(n)
          enddo
!
!- contravariant velocity in BL coordinates -!
          do m=0,3
            ucon1(m)=tmp1(m)
          enddo
!
!- get Boyer-Lindquist metric (contravariant metric) -!
          x1aa=x1(i)
          x2aa=x2(j)
          x3aa=x3(k)
!
!          if(x1aa .le. rbh1) then
!            call kshcon1(gcon1a,x1aa,x2aa,x3aa)
!          else
            call kercon1(gcon1a,x1aa,x2aa,x3aa)
!          endif
!
          alpha1=1./sqrt(-gcon1a(0,0))
          gfl=ucon1(0)*alpha1
!
          beta1(1)=alpha1*alpha1*gcon1a(0,1)
          beta1(2)=alpha1*alpha1*gcon1a(0,2)
          beta1(3)=alpha1*alpha1*gcon1a(0,3)
!- calculation new u^tilda
! 
          utiln(1)=ucon1(1)+gfl*beta1(1)*(1./alpha1)
          utiln(2)=ucon1(2)+gfl*beta1(2)*(1./alpha1)
          utiln(3)=ucon1(3)+gfl*beta1(3)*(1./alpha1)
!
        else
          utiln(1)=0.d0
          utiln(2)=0.d0
          utiln(3)=0.d0
        endif
!
!- copy new primitive variables array -!
        uria(1,i,j,k)=uri(1,i,j,k)
!        uria(2,i,j,k)=utiln(1)*(1./gfl)
!        uria(3,i,j,k)=utiln(2)*(1./gfl)
!        uria(4,i,j,k)=utiln(3)*(1./gfl)
        uria(2,i,j,k)=utiln(1)
        uria(3,i,j,k)=utiln(2)
        uria(4,i,j,k)=utiln(3)
        uria(5,i,j,k)=uri(5,i,j,k)
        uria(6,i,j,k)=uri(6,i,j,k)
        uria(7,i,j,k)=uri(7,i,j,k)
        uria(8,i,j,k)=uri(8,i,j,k)
        uria(9,i,j,k)=uri(9,i,j,k)  
!
      enddo
    enddo
  enddo
!
  deallocate( gcov1, gcon1, gcon1a, beta1, util, utiln, ucon, ucon1, &
              tmp1, trans, stat=merr)
      
  return
end subroutine transks2bl
!
!--------------------------------------------------------------------
subroutine cal4to4vel(uri,urib,gcov,gcon,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
!- Calculation of 4-velocity & 4-magnetic field -!
!- (for output of radiation calculation) -!
!
!- uri: primitive variables with \tilde{u} & 3-magnetic field (input)
!- urib: primitive vairbales with 4-vel & 4-mag (output)  
!
  use pram, only : nv, imax, jmax, kmax  
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1 
  integer :: merr, jh, kh

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1) ! prim variables with 4-vel -!
  real(8) :: urib(nv+1,is1:ie1,js1:je1,ks1:ke1) ! prim variables with 4-vel, 4-Bfield -!
  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)  

  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcon(:)
  real(8), allocatable :: util(:), bcon(:) 
  real(8) :: gfl
!
!- allocate variables -!  
  allocate(gcov1(0:3,0:3),  gcon1(0:3,0:3), stat=merr)
  allocate(util(1:3), bcon(1:3), stat=merr)
  allocate(ucov(0:3), ucon(0:3), bbcon(0:3), stat=merr)

!=====================================================================
!
  jh=(jmax+1)/2
  kh=((kmax+1)/2)  
  
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
         
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
        bcon(1)=uri(7,i,j,k)
        bcon(2)=uri(8,i,j,k)
        bcon(3)=uri(9,i,j,k)
!
!
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
!
        urib(1,i,j,k)=uri(1,i,j,k)
        urib(2,i,j,k)=ucon(0)
        urib(3,i,j,k)=ucon(1)
        urib(4,i,j,k)=ucon(2)
        urib(5,i,j,k)=ucon(3)
        urib(6,i,j,k)=uri(5,i,j,k)
        urib(7,i,j,k)=bbcon(0)
        urib(8,i,j,k)=bbcon(1)
        urib(9,i,j,k)=bbcon(2)
        urib(10,i,j,k)=bbcon(3)
!      
      enddo
    enddo
  enddo
!
  deallocate(gcov1, gcon1, util, bcon, ucov, ucon, bbcon, stat=merr)
!  
  return
end subroutine cal4to4vel
!
