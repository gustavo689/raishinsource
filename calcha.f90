!---------------------------------------------------------------------
subroutine calcha(uriir,urijr,urikr,uriil,urijl,urikl, &
                  gcovi,gcovj,gcovk,gconi,gconj,gconk, &
                  cmaxi,cmini,cmaxj,cminj,cmaxk,cmink,nm0, &
                  is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------
!
!     Calculate characteristics at cell-interface
!
  use pram, only : imax, jmax, kmax, nv, icha
  implicit none

  integer :: i, j, k, is1, ie1, js1, je1, ks1, ke1
  integer :: nm0
  integer :: merr
!
  real(8) :: uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), & 
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: gcovi(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gcovj(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gcovk(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: gconi(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gconj(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gconk(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: cmaxi(is1:ie1,js1:je1,ks1:ke1), cmini(is1:ie1,js1:je1,ks1:ke1)
  real(8) :: cmaxj(is1:ie1,js1:je1,ks1:ke1), cminj(is1:ie1,js1:je1,ks1:ke1)
  real(8) :: cmaxk(is1:ie1,js1:je1,ks1:ke1), cmink(is1:ie1,js1:je1,ks1:ke1)
!      
  real(8), allocatable :: vfipr(:,:,:), vfimr(:,:,:), &
           vfjpr(:,:,:), vfjmr(:,:,:), vfkpr(:,:,:), vfkmr(:,:,:), &
           vfipl(:,:,:), vfiml(:,:,:), vfjpl(:,:,:), vfjml(:,:,:), &
           vfkpl(:,:,:), vfkml(:,:,:)

  real(8) :: cmaxi1, cmini1, cmaxj1, cminj1, cmaxk1, cmink1
!
  allocate( vfipr(is1:ie1,js1:je1,ks1:ke1), vfimr(is1:ie1,js1:je1,ks1:ke1), &
            vfjpr(is1:ie1,js1:je1,ks1:ke1), vfjmr(is1:ie1,js1:je1,ks1:ke1), &
            vfkpr(is1:ie1,js1:je1,ks1:ke1), vfkmr(is1:ie1,js1:je1,ks1:ke1), &
            vfipl(is1:ie1,js1:je1,ks1:ke1), vfiml(is1:ie1,js1:je1,ks1:ke1), &
            vfjpl(is1:ie1,js1:je1,ks1:ke1), vfjml(is1:ie1,js1:je1,ks1:ke1), &
            vfkpl(is1:ie1,js1:je1,ks1:ke1), vfkml(is1:ie1,js1:je1,ks1:ke1), &
            stat=merr)
!
!=====================================================================@
!
!- initialize -!
  do k=ks1, ke1
    do j=js1, je1
      do i=is1, ie1
!          
        vfipr(i,j,k)=0.d0
        vfimr(i,j,k)=0.d0
        vfjpr(i,j,k)=0.d0
        vfjmr(i,j,k)=0.d0
        vfkpr(i,j,k)=0.d0
        vfkmr(i,j,k)=0.d0
!        
          
        vfipl(i,j,k)=0.d0
        vfiml(i,j,k)=0.d0
        vfjpl(i,j,k)=0.d0
        vfjml(i,j,k)=0.d0
        vfkpl(i,j,k)=0.d0
        vfkml(i,j,k)=0.d0
!        
      enddo
    enddo
  enddo  

  
  if(icha .eq. 0 .or. icha .eq. 1 .or. icha .eq. 2) then
    call calcha4i(uriir,gcovi,gconi,vfipr,vfimr,nm0,is1,ie1,js1,je1,ks1,ke1)
    call calcha4j(urijr,gcovj,gconj,vfjpr,vfjmr,nm0,is1,ie1,js1,je1,ks1,ke1)
    call calcha4k(urikr,gcovk,gconk,vfkpr,vfkmr,nm0,is1,ie1,js1,je1,ks1,ke1)
!
    call calcha4i(uriil,gcovi,gconi,vfipl,vfiml,nm0,is1,ie1,js1,je1,ks1,ke1)
    call calcha4j(urijl,gcovj,gconj,vfjpl,vfjml,nm0,is1,ie1,js1,je1,ks1,ke1)
    call calcha4k(urikl,gcovk,gconk,vfkpl,vfkml,nm0,is1,ie1,js1,je1,ks1,ke1)
  elseif(icha .eq. 3) then
    call calcha5i(uriir,gcovi,gconi,vfipr,vfimr,nm0,is1,ie1,js1,je1,ks1,ke1)
    call calcha5j(urijr,gcovj,gconj,vfjpr,vfjmr,nm0,is1,ie1,js1,je1,ks1,ke1)
    call calcha5k(urikr,gcovk,gconk,vfkpr,vfkmr,nm0,is1,ie1,js1,je1,ks1,ke1)
!
    call calcha5i(uriil,gcovi,gconi,vfipl,vfiml,nm0,is1,ie1,js1,je1,ks1,ke1)
    call calcha5j(urijl,gcovj,gconj,vfjpl,vfjml,nm0,is1,ie1,js1,je1,ks1,ke1)
    call calcha5k(urikl,gcovk,gconk,vfkpl,vfkml,nm0,is1,ie1,js1,je1,ks1,ke1)    
  endif
!
  do k=ks1+nm0-1,ke1-nm0
    do j=js1+nm0-1,je1-nm0
      do i=is1+nm0-1,ie1-nm0

         cmaxi1=max(0.d0,vfipr(i,j,k),vfipl(i,j,k))
         cmini1=min(0.d0,vfimr(i,j,k),vfiml(i,j,k))

         cmaxj1=max(0.d0,vfjpr(i,j,k),vfjpl(i,j,k))
         cminj1=min(0.d0,vfjmr(i,j,k),vfjml(i,j,k))

         cmaxk1=max(0.d0,vfkpr(i,j,k),vfkpl(i,j,k))
         cmink1=min(0.d0,vfkpr(i,j,k),vfkml(i,j,k))

         cmaxi(i,j,k)=cmaxi1
         cmini(i,j,k)=cmini1
         
         cmaxj(i,j,k)=cmaxj1
         cminj(i,j,k)=cminj1
     
         cmaxk(i,j,k)=cmaxk1
         cmink(i,j,k)=cmink1

      enddo
    enddo
  enddo

  deallocate( vfipr, vfimr, vfjpr, vfjmr, vfkpr, vfkmr, &
              vfipl, vfiml, vfjpl, vfjml, vfkpl, vfkml, stat=merr )
!
  return
end subroutine calcha
!
!---------------------------------------------------------------------@
subroutine calcha4(uri,gcov,gcon,vfip,vfim,vfjp,vfjm,vfkp,vfkm, &
                   is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Calculate characteristics for HLLE method at cell-interface
!
  use pram, only : imax, jmax, kmax, nv, gam, c0, metric, ieos, icha
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1
  integer :: merr
!
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: vfip(is1:ie1,js1:je1,ks1:ke1), vfim(is1:ie1,js1:je1,ks1:ke1)
  real(8) :: vfjp(is1:ie1,js1:je1,ks1:ke1), vfjm(is1:ie1,js1:je1,ks1:ke1)
  real(8) :: vfkp(is1:ie1,js1:je1,ks1:ke1), vfkm(is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:), gcov3a(:,:), gcon3a(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:)
  real(8), allocatable :: util(:), bcon(:), beta1(:), vcon(:), vcov(:)  

  real(8) :: de, pr, roh, roe, alpha1, gfl, bbsq, vsq, vb, &
             cssq, vasq, omsq, f1i, f1j, f1k, r1, f2, f3i, f3j, f3k
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), &
           gcov3a(1:3,1:3), gcon3a(1:3,1:3), stat=merr)
  allocate(util(1:3), bcon(1:3), beta1(1:3), vcon(1:3), vcov(1:3), stat=merr)
  allocate(ucov(0:3), ucon(0:3), bbcov(0:3), bbcon(0:3), stat=merr)
!     
!=====================================================================@

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
!
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
        alpha1=1./sqrt(-gcon1(0,0))
        beta1(1)=alpha1*alpha1*gcon1(0,1)
        beta1(2)=alpha1*alpha1*gcon1(0,2)
        beta1(3)=alpha1*alpha1*gcon1(0,3)
!
        call cal3met(gcov1,gcon1,gcov3a,gcon3a)

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

        vcon(1)=util(1)*(1./gfl)
        vcon(2)=util(2)*(1./gfl)
        vcon(3)=util(3)*(1./gfl)
!
!- cal of covariant 3-velocity
        
        call lower1(vcon,vcov,gcov3a)

!        vsq=((gfl**2)-1.)/(gfl**2)
!        vsq=vcon(1)*vcon(1)+vcon(2)*vcon(2)+vcon(3)*vcon(3)
!        vsq=vcov(1)*vcon(1)+vcov(2)*vcon(2)+vcov(3)*vcon(3)

        vsq=0.d0
        do m=1,3
          do n=1,3
            vsq=vsq+gcov3a(m,n)*vcon(m)*vcon(n)     
          enddo
        enddo         
!
!        vb=vcon(1)*bcon(1)+vcon(2)*bcon(2)+vcon(3)*bcon(3)
        vb=vcov(1)*bcon(1)+vcov(2)*bcon(2)+vcov(3)*bcon(3)
        
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          cssq=gam*pr/roh
        elseif(ieos .eq. 1) then
          cssq=(pr/(3.0*roh))*((5.0*roh-8.0*pr)/(roh-pr))
        endif
         
        vasq=bbsq*(1./(roh+bbsq))
        omsq=vasq+cssq*(1.0-vasq)

        if(icha .eq. 0) then
!
!  --- Gammie et al. (2003)

          vfip(i,j,k)=vcon(1)+sqrt(omsq)
          vfjp(i,j,k)=vcon(2)+sqrt(omsq)
          vfkp(i,j,k)=vcon(3)+sqrt(omsq)
         
          vfim(i,j,k)=vcon(1)-sqrt(omsq)
          vfjm(i,j,k)=vcon(2)-sqrt(omsq)
          vfkm(i,j,k)=vcon(3)-sqrt(omsq)
          
        elseif(icha .eq. 1) then
!
!  --- Leismann et al. (2005)

          f1i=vcon(1)*(1.0-omsq)
          f1j=vcon(2)*(1.0-omsq)
          f1k=vcon(3)*(1.0-omsq)

          r1=cssq*(vb**2)*(1./((roh+bbsq)*gfl**2))
         
          f2=1.0-vsq*omsq-r1
         
!          f3i=((vsq-1.0)*omsq+r1)*((vsq-(vcon(1)**2))*omsq &
!             +(vcon(1)**2)-1.0+r1)
!          f3j=((vsq-1.0)*omsq+r1)*((vsq-(vcon(2)**2))*omsq &
!             +(vcon(2)**2)-1.0+r1)
!          f3k=((vsq-1.0)*omsq+r1)*((vsq-(vcon(3)**2))*omsq &
!             +(vcon(3)**2)-1.0+r1)
          
          f3i=((vsq-1.0)*omsq+r1)*((gcon3a(1,1)*vsq-(vcon(1)**2))*omsq &
             +(vcon(1)**2)-gcon3a(1,1)+r1)
          f3j=((vsq-1.0)*omsq+r1)*((gcon3a(2,2)*vsq-(vcon(2)**2))*omsq &
             +(vcon(2)**2)-gcon3a(2,2)+r1)
          f3k=((vsq-1.0)*omsq+r1)*((gcon3a(3,3)*vsq-(vcon(3)**2))*omsq &
             +(vcon(3)**2)-gcon3a(3,3)+r1)
          
          vfip(i,j,k)=(f1i+sqrt(f3i))*(1./f2)
          vfjp(i,j,k)=(f1j+sqrt(f3j))*(1./f2)
          vfkp(i,j,k)=(f1k+sqrt(f3k))*(1./f2)
         
          vfim(i,j,k)=(f1i-sqrt(f3i))*(1./f2)
          vfjm(i,j,k)=(f1j-sqrt(f3j))*(1./f2)
          vfkm(i,j,k)=(f1k-sqrt(f3k))*(1./f2)

        elseif(icha .eq. 2) then
!
! --- Del Zanna et al. (2007)
         
          f1i=vcon(1)*(1.0-omsq)
          f1j=vcon(2)*(1.0-omsq)
          f1k=vcon(3)*(1.0-omsq)
         
          f2=1.0-vsq*omsq
         
!          f3i=omsq*(1.0-vsq)*((1.0-vsq*omsq) &
!                             -(1.0-omsq)*vcon(1)*vcon(1))
!          f3j=omsq*(1.0-vsq)*((1.0-vsq*omsq) &
!                             -(1.0-omsq)*vcon(2)*vcon(2))
!          f3k=omsq*(1.0-vsq)*((1.0-vsq*omsq) &
!                             -(1.0-omsq)*vcon(3)*vcon(3))

          f3i=omsq*(1.0-vsq)*((1.0-vsq*omsq)*gcon3a(1,1) &
                             -(1.0-omsq)*vcon(1)*vcon(1))
          f3j=omsq*(1.0-vsq)*((1.0-vsq*omsq)*gcon3a(2,2) &
                             -(1.0-omsq)*vcon(2)*vcon(2))
          f3k=omsq*(1.0-vsq)*((1.0-vsq*omsq)*gcon3a(3,3) &
                             -(1.0-omsq)*vcon(3)*vcon(3))
          
          vfip(i,j,k)=(f1i+sqrt(f3i))*(1./f2)
          vfjp(i,j,k)=(f1j+sqrt(f3j))*(1./f2)
          vfkp(i,j,k)=(f1k+sqrt(f3k))*(1./f2)
         
          vfim(i,j,k)=(f1i-sqrt(f3i))*(1./f2)
          vfjm(i,j,k)=(f1j-sqrt(f3j))*(1./f2)
          vfkm(i,j,k)=(f1k-sqrt(f3k))*(1./f2)

        endif 
!
!- correction for GR -!
        vfip(i,j,k)=alpha1*vfip(i,j,k)-beta1(1)
        vfjp(i,j,k)=alpha1*vfjp(i,j,k)-beta1(2)
        vfkp(i,j,k)=alpha1*vfkp(i,j,k)-beta1(3)

        vfim(i,j,k)=alpha1*vfim(i,j,k)-beta1(1)
        vfjm(i,j,k)=alpha1*vfjm(i,j,k)-beta1(2)
        vfkm(i,j,k)=alpha1*vfkm(i,j,k)-beta1(3)
!- check wave speed exceed light speed -!
!
        if(abs(vfip(i,j,k)) .gt. 1.d0) then
          vfip(i,j,k)=1.d0
        endif
        if(abs(vfjp(i,j,k)) .gt. 1.d0) then
          vfjp(i,j,k)=1.d0
        endif
        if(abs(vfkp(i,j,k)) .gt. 1.d0) then
          vfkp(i,j,k)=1.d0
        endif
        if(abs(vfim(i,j,k)) .gt. 1.d0) then
          vfim(i,j,k)=-1.d0
        endif
        if(abs(vfjm(i,j,k)) .gt. 1.d0) then
          vfjm(i,j,k)=-1.d0
        endif
        if(abs(vfkm(i,j,k)) .gt. 1.d0) then
          vfkm(i,j,k)=-1.d0
        endif
         
      enddo
    enddo
  enddo
!
  deallocate( gcov1, gcon1, gcov3a, gcon3a, util, bcon, beta1, &
              vcov, vcon, ucov, ucon, bbcov, bbcon, stat=merr)
!
   return
end subroutine calcha4
!
!---------------------------------------------------------------------@
subroutine calcha4i(uri,gcov,gcon,vfip,vfim,nm0,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Calculate characteristics for HLLE method at cell-interface
!
  use pram, only : imax, jmax, kmax, nv, gam, c0, metric, ieos, icha
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1
  integer :: nm0, merr
!
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: vfip(is1:ie1,js1:je1,ks1:ke1), vfim(is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:), gcov3a(:,:), gcon3a(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:)
  real(8), allocatable :: util(:), bcon(:), beta1(:), vcon(:), vcov(:)  
!
  real(8) :: de, pr, roh, roe, alpha1, gfl, bbsq, vsq, vb, &
             cssq, vasq, omsq, f1i, r1, f2, f3i
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), &
           gcov3a(1:3,1:3), gcon3a(1:3,1:3), stat=merr)
  allocate(util(1:3), bcon(1:3), beta1(1:3), vcon(1:3), vcov(1:3), stat=merr)
  allocate(ucov(0:3), ucon(0:3), bbcov(0:3), bbcon(0:3), stat=merr)
!     
!=====================================================================@

  do k=ks1+nm0-1,ke1-nm0
    do j=js1+nm0-1,je1-nm0
      do i=is1+nm0-1,ie1-nm0
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
!
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
        alpha1=1./sqrt(-gcon1(0,0))
        beta1(1)=alpha1*alpha1*gcon1(0,1)
        beta1(2)=alpha1*alpha1*gcon1(0,2)
        beta1(3)=alpha1*alpha1*gcon1(0,3)
!
        call cal3met(gcov1,gcon1,gcov3a,gcon3a)
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
        vcon(1)=util(1)*(1./gfl)
        vcon(2)=util(2)*(1./gfl)
        vcon(3)=util(3)*(1./gfl)
!
!- cal of covariant 3-velocity
        
        call lower1(vcon,vcov,gcov3a)

!        vsq=((gfl**2)-1.)/(gfl**2)
!        vsq=vcon(1)*vcon(1)+vcon(2)*vcon(2)+vcon(3)*vcon(3)
!        vsq=vcov(1)*vcon(1)+vcov(2)*vcon(2)+vcov(3)*vcon(3)

        vsq=0.d0
        do m=1,3
          do n=1,3
            vsq=vsq+gcov3a(m,n)*vcon(m)*vcon(n)     
          enddo
        enddo 

!
!        vb=vcon(1)*bcon(1)+vcon(2)*bcon(2)+vcon(3)*bcon(3)
        vb=vcov(1)*bcon(1)+vcov(2)*bcon(2)+vcov(3)*bcon(3)
!         
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          cssq=gam*pr/roh
        elseif(ieos .eq. 1) then
          cssq=(pr/(3.0*roh))*((5.0*roh-8.0*pr)/(roh-pr))
        endif
         
        vasq=bbsq*(1./(roh+bbsq))
        omsq=vasq+cssq*(1.0-vasq)

        if(icha .eq. 0) then
!
!  --- Gammie et al. (2003)

          vfip(i,j,k)=vcon(1)+sqrt(omsq)
          vfim(i,j,k)=vcon(1)-sqrt(omsq)
          
        elseif(icha .eq. 1) then
!
!  --- Leismann et al. (2005)

          f1i=vcon(1)*(1.0-omsq)
          r1=cssq*(vb**2)*(1./((roh+bbsq)*gfl**2))
          f2=1.0-vsq*omsq-r1         
!          f3i=((vsq-1.0)*omsq+r1)*((vsq-vcon(1)**2)*omsq+vcon(1)**2-1.0+r1)
          f3i=((vsq-1.0)*omsq+r1)*((gcon3a(1,1)*vsq-vcon(1)**2)*omsq &
                                   +vcon(1)**2-gcon3a(1,1)+r1)

          vfip(i,j,k)=(f1i+sqrt(f3i))*(1./f2)         
          vfim(i,j,k)=(f1i-sqrt(f3i))*(1./f2)


        elseif(icha .eq. 2) then
!
! --- Del Zanna et al. (2007)
         
          f1i=vcon(1)*(1.0-omsq)         
          f2=1.0-vsq*omsq         

!          f3i=omsq*(1.0-vsq)*((1.0-vsq*omsq) &
!                             -(1.0-omsq)*vcon(1)*vcon(1))  

          f3i=omsq*(1.0-vsq)*((1.0-vsq*omsq)*gcon3a(1,1) &
                             -(1.0-omsq)*vcon(1)*vcon(1))
          
          vfip(i,j,k)=(f1i+sqrt(f3i))*(1./f2)         
          vfim(i,j,k)=(f1i-sqrt(f3i))*(1./f2)

        endif 
!
        vfip(i,j,k)=alpha1*vfip(i,j,k)-beta1(1)
        vfim(i,j,k)=alpha1*vfim(i,j,k)-beta1(1)
!
!- check wave speed exceed light speed -!

        if(abs(vfip(i,j,k)) .gt. 1.d0) then
          vfip(i,j,k)=1.d0
        endif
        if(abs(vfim(i,j,k)) .gt. 1.d0) then
          vfim(i,j,k)=-1.d0
        endif
!       
      enddo
    enddo
  enddo
!
  deallocate( gcov1, gcon1, gcov3a, gcon3a, util, bcon, beta1, &
              vcov, vcon, ucov, ucon, bbcov, bbcon, stat=merr)
!
   return
end subroutine calcha4i
!
!---------------------------------------------------------------------@
subroutine calcha4j(uri,gcov,gcon,vfjp,vfjm,nm0,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Calculate characteristics for HLLE method at cell-interface
!
  use pram, only : imax, jmax, kmax, nv, gam, c0, metric, ieos, icha
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1
  integer :: nm0, merr
!
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: vfjp(is1:ie1,js1:je1,ks1:ke1), vfjm(is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:), gcov3a(:,:), gcon3a(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:)
  real(8), allocatable :: util(:), bcon(:), beta1(:), vcon(:), vcov(:)  
!
  real(8) :: de, pr, roh, roe, alpha1, gfl, bbsq, vsq, vb, &
             cssq, vasq, omsq, f1j, r1, f2, f3j
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), &
           gcov3a(1:3,1:3), gcon3a(1:3,1:3), stat=merr)
  allocate(util(1:3), bcon(1:3), beta1(1:3), vcon(1:3), vcov(1:3), stat=merr)
  allocate(ucov(0:3), ucon(0:3), bbcov(0:3), bbcon(0:3), stat=merr)
!     
!=====================================================================@

  do k=ks1+nm0-1,ke1-nm0
    do j=js1+nm0-1,je1-nm0
      do i=is1+nm0-1,ie1-nm0
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
!
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
        alpha1=1./sqrt(-gcon1(0,0))
        beta1(1)=alpha1*alpha1*gcon1(0,1)
        beta1(2)=alpha1*alpha1*gcon1(0,2)
        beta1(3)=alpha1*alpha1*gcon1(0,3)

        call cal3met(gcov1,gcon1,gcov3a,gcon3a)
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

        vcon(1)=util(1)*(1./gfl)
        vcon(2)=util(2)*(1./gfl)
        vcon(3)=util(3)*(1./gfl)
!
!- cal of covariant 3-velocity
        
        call lower1(vcon,vcov,gcov3a)

!        vsq=((gfl**2)-1.)/(gfl**2)
!        vsq=vcon(1)*vcon(1)+vcon(2)*vcon(2)+vcon(3)*vcon(3)
!        vsq=vcov(1)*vcon(1)+vcov(2)*vcon(2)+vcov(3)*vcon(3)

        vsq=0.d0
        do m=1,3
          do n=1,3
            vsq=vsq+gcov3a(m,n)*vcon(m)*vcon(n)     
          enddo
        enddo
!
!        vb=vcon(1)*bcon(1)+vcon(2)*bcon(2)+vcon(3)*bcon(3)
        vb=vcov(1)*bcon(1)+vcov(2)*bcon(2)+vcov(3)*bcon(3)
         
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          cssq=gam*pr/roh
        elseif(ieos .eq. 1) then
          cssq=(pr/(3.0*roh))*((5.0*roh-8.0*pr)/(roh-pr))
        endif
         
        vasq=bbsq/(roh+bbsq)
        omsq=vasq+cssq*(1.0-vasq)

        if(icha .eq. 0) then
!
!  --- Gammie et al. (2003)

          vfjp(i,j,k)=vcon(2)+sqrt(omsq)
          vfjm(i,j,k)=vcon(2)-sqrt(omsq)
          
        elseif(icha .eq. 1) then
!
!  --- Leismann et al. (2005)

          f1j=vcon(2)*(1.0-omsq)
          r1=cssq*(vb**2)*(1./((roh+bbsq)*gfl**2))
          f2=1.0-vsq*omsq-r1
!          f3j=((vsq-1.0)*omsq+r1)*((vsq-vcon(2)**2)*omsq+vcon(2)**2-1.0+r1)
          f3j=((vsq-1.0)*omsq+r1)*((gcon3a(2,2)*vsq-vcon(2)**2)*omsq &
                                   +vcon(2)**2-gcon3a(2,2)+r1)
          
          vfjp(i,j,k)=(f1j+sqrt(f3j))*(1./f2)
          vfjm(i,j,k)=(f1j-sqrt(f3j))*(1./f2)

        elseif(icha .eq. 2) then
!
! --- Del Zanna et al. (2007)
         
          f1j=vcon(2)*(1.0-omsq)
          f2=1.0-vsq*omsq
!          f3j=omsq*(1.0-vsq)*((1.0-vsq*omsq) &
!                             -(1.0-omsq)*vcon(2)*vcon(2))         
          f3j=omsq*(1.0-vsq)*((1.0-vsq*omsq)*gcon3a(2,2) &
                             -(1.0-omsq)*vcon(2)*vcon(2))
          
          vfjp(i,j,k)=(f1j+sqrt(f3j))*(1./f2)
          vfjm(i,j,k)=(f1j-sqrt(f3j))*(1./f2)

        endif 
!
        vfjp(i,j,k)=alpha1*vfjp(i,j,k)-beta1(2)
        vfjm(i,j,k)=alpha1*vfjm(i,j,k)-beta1(2)
!
!- check wave speed exceed light speed -!

        if(abs(vfjp(i,j,k)) .gt. 1.d0) then
          vfjp(i,j,k)=1.d0
        endif
        if(abs(vfjm(i,j,k)) .gt. 1.d0) then
          vfjm(i,j,k)=-1.d0
        endif
!       
      enddo
    enddo
  enddo
!
  deallocate( gcov1, gcon1, gcov3a, gcon3a, util, bcon, beta1, &
              vcov, vcon, ucov, ucon, bbcov, bbcon, stat=merr)
!
   return
end subroutine calcha4j
!
!---------------------------------------------------------------------@
subroutine calcha4k(uri,gcov,gcon,vfkp,vfkm,nm0,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Calculate characteristics for HLLE method at cell-interface
!
  use pram, only : imax, jmax, kmax, nv, gam, c0, metric, ieos, icha
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1
  integer :: nm0, merr
!
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: vfkp(is1:ie1,js1:je1,ks1:ke1), vfkm(is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:), gcov3a(:,:), gcon3a(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:)
  real(8), allocatable :: util(:), bcon(:), beta1(:), vcon(:), vcov(:)  
!
  real(8) :: de, pr, roh, roe, alpha1, gfl, bbsq, vsq, vb, &
             cssq, vasq, omsq, f1k, r1, f2, f3k
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), &
           gcov3a(1:3,1:3), gcon3a(1:3,1:3), stat=merr)
  allocate(util(1:3), bcon(1:3), beta1(1:3), vcon(1:3), vcov(1:3), stat=merr)
  allocate(ucov(0:3), ucon(0:3), bbcov(0:3), bbcon(0:3), stat=merr)
!     
!=====================================================================@

  do k=ks1+nm0-1,ke1-nm0
    do j=js1+nm0-1,je1-nm0
      do i=is1+nm0-1,ie1-nm0
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
!
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
        alpha1=1./sqrt(-gcon1(0,0))
        beta1(1)=alpha1*alpha1*gcon1(0,1)
        beta1(2)=alpha1*alpha1*gcon1(0,2)
        beta1(3)=alpha1*alpha1*gcon1(0,3)

        call cal3met(gcov1,gcon1,gcov3a,gcon3a)
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

        vcon(1)=util(1)*(1./gfl)
        vcon(2)=util(2)*(1./gfl)
        vcon(3)=util(3)*(1./gfl)
!
!- cal of covariant 3-velocity
        
        call lower1(vcon,vcov,gcov3a)

!        vsq=((gfl**2)-1.)/(gfl**2)
!        vsq=vcon(1)*vcon(1)+vcon(2)*vcon(2)+vcon(3)*vcon(3)
!        vsq=vcov(1)*vcon(1)+vcov(2)*vcon(2)+vcov(3)*vcon(3)

        vsq=0.d0
        do m=1,3
          do n=1,3
            vsq=vsq+gcov3a(m,n)*vcon(m)*vcon(n)     
          enddo
        enddo 
!
!        vb=vcon(1)*bcon(1)+vcon(2)*bcon(2)+vcon(3)*bcon(3)
        vb=vcov(1)*bcon(1)+vcov(2)*bcon(2)+vcov(3)*bcon(3)
         
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          cssq=gam*pr*(1./roh)
        elseif(ieos .eq. 1) then
          cssq=(pr/(3.0*roh))*((5.0*roh-8.0*pr)/(roh-pr))
        endif
         
        vasq=bbsq*(1./(roh+bbsq))
        omsq=vasq+cssq*(1.0-vasq)

        if(icha .eq. 0) then
!
!  --- Gammie et al. (2003)

          vfkp(i,j,k)=vcon(3)+sqrt(omsq)
          vfkm(i,j,k)=vcon(3)-sqrt(omsq)
          
        elseif(icha .eq. 1) then
!
!  --- Leismann et al. (2005)

          f1k=vcon(3)*(1.0-omsq)
          r1=cssq*(vb**2)*(1./((roh+bbsq)*gfl**2))
          f2=1.0-vsq*omsq-r1
!          f3k=((vsq-1.0)*omsq+r1)*((vsq-vcon(3)**2)*omsq+vcon(3)**2-1.0+r1)
          f3k=((vsq-1.0)*omsq+r1)*((gcon3a(3,3)*vsq-vcon(3)**2)*omsq &
                                   +vcon(3)**2-gcon3a(3,3)+r1)
          
          vfkp(i,j,k)=(f1k+sqrt(f3k))*(1./f2)
          vfkm(i,j,k)=(f1k-sqrt(f3k))*(1./f2)

        elseif(icha .eq. 2) then
!
! --- Del Zanna et al. (2007)
         
          f1k=vcon(3)*(1.0-omsq)
          f2=1.0-vsq*omsq
!          f3k=omsq*(1.0-vsq)*((1.0-vsq*omsq) &
!                             -(1.0-omsq)*vcon(3)*vcon(3))
          f3k=omsq*(1.0-vsq)*((1.0-vsq*omsq)*gcon3a(3,3) &
                             -(1.0-omsq)*vcon(3)*vcon(3))
          
          vfkp(i,j,k)=(f1k+sqrt(f3k))*(1./f2)
          vfkm(i,j,k)=(f1k-sqrt(f3k))*(1./f2)

        endif 

        vfkp(i,j,k)=alpha1*vfkp(i,j,k)-beta1(3)
        vfkm(i,j,k)=alpha1*vfkm(i,j,k)-beta1(3)

!- check wave speed exceed light speed -!

        if(abs(vfkp(i,j,k)) .gt. 1.d0) then
          vfkp(i,j,k)=1.d0
        endif
        if(abs(vfkm(i,j,k)) .gt. 1.d0) then
          vfkm(i,j,k)=-1.d0
        endif
       
      enddo
    enddo
  enddo
!
  deallocate( gcov1, gcon1, gcov3a, gcon3a, util, bcon, beta1, &
              vcov, vcon, ucov, ucon, bbcov, bbcon, stat=merr)
!
   return
end subroutine calcha4k
!
!
!---------------------------------------------------------------------@
subroutine calcha5(uri,gcov,gcon,vfip,vfim,vfjp,vfjm,vfkp,vfkm, &
                   is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Calculate characteristics at cell-interface
!     (from Gammie et al. 2003)
!
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1
  integer :: merr
!
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: vfip(is1:ie1,js1:je1,ks1:ke1), vfim(is1:ie1,js1:je1,ks1:ke1)
  real(8) :: vfjp(is1:ie1,js1:je1,ks1:ke1), vfjm(is1:ie1,js1:je1,ks1:ke1)
  real(8) :: vfkp(is1:ie1,js1:je1,ks1:ke1), vfkm(is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:), gcov3a(:,:), gcon3a(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcon(:), bbcov(:)  
  real(8), allocatable :: util(:), vcon(:), vcona(:), vcova(:), bcon(:), &
                          beta1(:)

  real(8) :: de, pr, roh, roe, alpha1, alpsq, gfl, bbsq, vsq, &
             cssq, vasq, omsq, vnsq, vnsq1, f1, f2i, f2j, f2k, &
             f3, f4i, f4j, f4k
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), gcov3a(1:3,1:3), gcon3a(1:3,1:3), &
       stat=merr)
  allocate(ucov(0:3), ucon(0:3), bbcon(0:3), bbcov(0:3), stat=merr)
  allocate(util(1:3), vcon(1:3), vcona(1:3), vcova(1:3), bcon(1:3), &
           beta1(1:3), stat=merr)
!     
!=====================================================================@

  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!
!- copy of primitive variables -!
        de=uri(1,i,j,k)
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
        pr=uri(5,i,j,k)
        bcon(1)=uri(7,i,j,k)
        bcon(2)=uri(8,i,j,k)
        bcon(3)=uri(9,i,j,k)
!
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
        call cal3met(gcov1,gcon1,gcov3a,gcon3a)
!
        alpha1=1./sqrt(-gcon1(0,0))
        alpsq=alpha1*alpha1
!
        beta1(1)=alpsq*gcon1(0,1)
        beta1(2)=alpsq*gcon1(0,2)
        beta1(3)=alpsq*gcon1(0,3)
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
        vcon(1)=util(1)*(1./gfl)
        vcon(2)=util(2)*(1./gfl)
        vcon(3)=util(3)*(1./gfl)
!
        bcon(1)=bcon(1)*(1./alpha1)
        bcon(2)=bcon(2)*(1./alpha1)
        bcon(3)=bcon(3)*(1./alpha1)
!        
        vsq=0.d0
        do m=1,3
          do n=1,3
            vsq=vsq+gcov3a(m,n)*vcon(m)*vcon(n)     
          enddo
        enddo
!         
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          cssq=gam*pr/roh
        elseif(ieos .eq. 1) then
          cssq=(pr/(3.0*roh))*((5.0*roh-8.0*pr)/(roh-pr))
        endif
!         
        vasq=bbsq*(1./(roh+bbsq))
        omsq=vasq+cssq*(1.0-vasq)
!
!- v^i + beta^i
        vcona(1)=vcon(1)+beta1(1)
        vcona(2)=vcon(2)+beta1(2)
        vcona(3)=vcon(3)+beta1(3)
        call lower1(vcona,vcova,gcov3a)
!
        vnsq1=vcova(1)*vcona(1)+vcova(2)*vcona(2)+vcova(3)*vcona(3)
!
        vnsq=0.d0
        do m=1,3
          do n=1,3
            vnsq=vnsq+gcov3a(m,n)*vcona(m)*vcona(n)
          enddo
        enddo
!
        f1=alpsq-vnsq1*omsq
        f2i=vcon(1)*alpsq*(1.-omsq)-beta1(1)*omsq*(alpsq-vnsq)
        f2j=vcon(2)*alpsq*(1.-omsq)-beta1(2)*omsq*(alpsq-vnsq)
        f2k=vcon(3)*alpsq*(1.-omsq)-beta1(3)*omsq*(alpsq-vnsq)

        f3=alpha1*sqrt(omsq)
        f4i=(alpsq-vnsq)*(gcov3a(1,1)*(alpsq-vnsq*omsq) &
                                 -(1.-omsq)*vcona(1)*vcona(1))
        f4j=(alpsq-vnsq)*(gcov3a(2,2)*(alpsq-vnsq*omsq) &
                                 -(1.-omsq)*vcona(2)*vcona(2))
        f4k=(alpsq-vnsq)*(gcov3a(3,3)*(alpsq-vnsq*omsq) &
                                 -(1.-omsq)*vcona(3)*vcona(3))

        vfip(i,j,k)=(1./f1)*(f2i+f3*sqrt(f4i))
        vfjp(i,j,k)=(1./f1)*(f2j+f3*sqrt(f4j))
        vfkp(i,j,k)=(1./f1)*(f2k+f3*sqrt(f4k))
!
        vfim(i,j,k)=(1./f1)*(f2i-f3*sqrt(f4i))
        vfjm(i,j,k)=(1./f1)*(f2j-f3*sqrt(f4j))
        vfkm(i,j,k)=(1./f1)*(f2k-f3*sqrt(f4k))
!
!- check wave speed exceed light speed -!

          if(vfip(i,j,k) .gt. 1.d0) then
            vfip(i,j,k)=1.d0
          endif
          if(vfjp(i,j,k) .gt. 1.d0) then
            vfjp(i,j,k)=1.d0
          endif
          if(vfkp(i,j,k) .gt. 1.d0) then
            vfkp(i,j,k)=1.d0
          endif
          if(abs(vfim(i,j,k)) .gt. 1.d0) then
            vfim(i,j,k)=-1.d0
          endif
          if(abs(vfjm(i,j,k)) .gt. 1.d0) then
            vfjm(i,j,k)=-1.d0
          endif
          if(abs(vfkm(i,j,k)) .gt. 1.d0) then
            vfkm(i,j,k)=-1.d0
          endif
!
      enddo
    enddo
  enddo
!
  deallocate( gcov1, gcon1, gcov3a, gcon3a, ucov, ucon, bbcon, bbcov, &
              util, vcon, vcona, vcova, bcon, beta1, stat=merr)
!
   return
end subroutine calcha5
!
!---------------------------------------------------------------------@
subroutine calcha5i(uri,gcov,gcon,vfip,vfim,nm0,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Calculate characteristics at cell-interface
!     (from Gammie et al. 2003)
!
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1
  integer :: nm0, merr
!
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: vfip(is1:ie1,js1:je1,ks1:ke1), vfim(is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:), gcov3a(:,:), gcon3a(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcon(:), bbcov(:)  
  real(8), allocatable :: util(:), vcon(:), vcona(:), vcova(:), bcon(:), &
                          beta1(:)

  real(8) :: de, pr, roh, roe, alpha1, alpsq, gfl, bbsq, vsq, &
             cssq, vasq, omsq, vnsq, vnsq1, f1, f2i, f3, f4i
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), gcov3a(1:3,1:3), gcon3a(1:3,1:3),&
       stat=merr)
  allocate(ucov(0:3), ucon(0:3), bbcon(0:3), bbcov(0:3), stat=merr)
  allocate(util(1:3), vcon(1:3), vcona(1:3), vcova(1:3), bcon(1:3), &
           beta1(1:3), stat=merr)
!     
!=====================================================================@

  do k=ks1+nm0-1,ke1-nm0
    do j=js1+nm0-1,je1-nm0
      do i=is1+nm0-1,ie1-nm0
!
!- copy of primitive variables -!
        de=uri(1,i,j,k)
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
        pr=uri(5,i,j,k)
        bcon(1)=uri(7,i,j,k)
        bcon(2)=uri(8,i,j,k)
        bcon(3)=uri(9,i,j,k)
!
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
        call cal3met(gcov1,gcon1,gcov3a,gcon3a)         
!
        alpha1=1./sqrt(-gcon1(0,0))
        alpsq=alpha1*alpha1
!
        beta1(1)=alpsq*gcon1(0,1)
        beta1(2)=alpsq*gcon1(0,2)
        beta1(3)=alpsq*gcon1(0,3)
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
        vcon(1)=util(1)*(1./gfl)
        vcon(2)=util(2)*(1./gfl)
        vcon(3)=util(3)*(1./gfl)
!
        bcon(1)=bcon(1)*(1./alpha1)
        bcon(2)=bcon(2)*(1./alpha1)
        bcon(3)=bcon(3)*(1./alpha1)
!        
        vsq=0.d0
        do m=1,3
          do n=1,3
            vsq=vsq+gcov3a(m,n)*vcon(m)*vcon(n)     
          enddo
        enddo
!         
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          cssq=gam*pr*(1./roh)
        elseif(ieos .eq. 1) then
          cssq=(pr/(3.0*roh))*((5.0*roh-8.0*pr)/(roh-pr))
        endif
!         
        vasq=bbsq*(1./(roh+bbsq))
        omsq=vasq+cssq*(1.0-vasq)
!
!- v^i + beta^i
        vcona(1)=vcon(1)+beta1(1)
        vcona(2)=vcon(2)+beta1(2)
        vcona(3)=vcon(3)+beta1(3)
        call lower1(vcona,vcova,gcov3a)
!
        vnsq1=vcova(1)*vcona(1)+vcova(2)*vcona(2)+vcova(3)*vcona(3)
!
        vnsq=0.d0
        do m=1,3
          do n=1,3
            vnsq=vnsq+gcov3a(m,n)*vcona(m)*vcona(n)
          enddo
        enddo
!
        f1=alpsq-vnsq1*omsq
        f2i=vcon(1)*alpsq*(1.-omsq)-beta1(1)*omsq*(alpsq-vnsq)

        f3=alpha1*sqrt(omsq)
        f4i=(alpsq-vnsq)*(gcov3a(1,1)*(alpsq-vnsq*omsq) &
                                 -(1.-omsq)*vcona(1)*vcona(1))

        vfip(i,j,k)=(1./f1)*(f2i+f3*sqrt(f4i))
        vfim(i,j,k)=(1./f1)*(f2i-f3*sqrt(f4i))
!
!- check wave speed exceed light speed -!

          if(vfip(i,j,k) .gt. 1.d0) then
            vfip(i,j,k)=1.d0
          endif
          if(abs(vfim(i,j,k)) .gt. 1.d0) then
            vfim(i,j,k)=-1.d0
          endif
!
      enddo
    enddo
  enddo
!
  deallocate( gcov1, gcon1, gcov3a, gcon3a, ucov, ucon, bbcon, bbcov, &
              util, vcon, vcona, vcova, bcon, beta1, stat=merr)
!
   return
end subroutine calcha5i
!
!---------------------------------------------------------------------@
subroutine calcha5j(uri,gcov,gcon,vfjp,vfjm,nm0,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Calculate characteristics at cell-interface
!     (from Gammie et al. 2003)
!
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1
  integer :: nm0, merr
!
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: vfjp(is1:ie1,js1:je1,ks1:ke1), vfjm(is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:), gcov3a(:,:), gcon3a(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcon(:), bbcov(:)  
  real(8), allocatable :: util(:), vcon(:), vcona(:), vcova(:), bcon(:), &
                          beta1(:)

  real(8) :: de, pr, roh, roe, alpha1, alpsq, gfl, bbsq, vsq, &
             cssq, vasq, omsq, vnsq, vnsq1, f1, f2j, f3, f4j
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), gcov3a(1:3,1:3), gcon3a(1:3,1:3), &
       stat=merr)
  allocate(ucov(0:3), ucon(0:3), bbcon(0:3), bbcov(0:3), stat=merr)
  allocate(util(1:3), vcon(1:3), vcona(1:3), vcova(1:3), bcon(1:3), &
           beta1(1:3), stat=merr)
!     
!=====================================================================@

  do k=ks1+nm0-1,ke1-nm0
    do j=js1+nm0-1,je1-nm0
      do i=is1+nm0-1,ie1-nm0  
!
!- copy of primitive variables -!
        de=uri(1,i,j,k)
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
        pr=uri(5,i,j,k)
        bcon(1)=uri(7,i,j,k)
        bcon(2)=uri(8,i,j,k)
        bcon(3)=uri(9,i,j,k)
!
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
        call cal3met(gcov1,gcon1,gcov3a,gcon3a)
!
        alpha1=1./sqrt(-gcon1(0,0))
        alpsq=alpha1*alpha1
!
        beta1(1)=alpsq*gcon1(0,1)
        beta1(2)=alpsq*gcon1(0,2)
        beta1(3)=alpsq*gcon1(0,3)
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
        vcon(1)=util(1)*(1./gfl)
        vcon(2)=util(2)*(1./gfl)
        vcon(3)=util(3)*(1./gfl)
!
        bcon(1)=bcon(1)*(1./alpha1)
        bcon(2)=bcon(2)*(1./alpha1)
        bcon(3)=bcon(3)*(1./alpha1)
!        
        vsq=0.d0
        do m=1,3
          do n=1,3
            vsq=vsq+gcov3a(m,n)*vcon(m)*vcon(n)     
          enddo
        enddo
!         
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          cssq=gam*pr/roh
        elseif(ieos .eq. 1) then
          cssq=(pr/(3.0*roh))*((5.0*roh-8.0*pr)/(roh-pr))
        endif
!         
        vasq=bbsq*(1./(roh+bbsq))
        omsq=vasq+cssq*(1.0-vasq)
!
!- v^i + beta^i
        vcona(1)=vcon(1)+beta1(1)
        vcona(2)=vcon(2)+beta1(2)
        vcona(3)=vcon(3)+beta1(3)
        call lower1(vcona,vcova,gcov3a)
!
        vnsq1=vcova(1)*vcona(1)+vcova(2)*vcona(2)+vcova(3)*vcona(3)
!
        vnsq=0.d0
        do m=1,3
          do n=1,3
            vnsq=vnsq+gcov3a(m,n)*vcona(m)*vcona(n)
          enddo
        enddo
!
        f1=alpsq-vnsq1*omsq
        f2j=vcon(2)*alpsq*(1.-omsq)-beta1(2)*omsq*(alpsq-vnsq)

        f3=alpha1*sqrt(omsq)
        f4j=(alpsq-vnsq)*(gcov3a(2,2)*(alpsq-vnsq*omsq) &
                                 -(1.-omsq)*vcona(2)*vcona(2))

        vfjp(i,j,k)=(1./f1)*(f2j+f3*sqrt(f4j))
        vfjm(i,j,k)=(1./f1)*(f2j-f3*sqrt(f4j))
!
!- check wave speed exceed light speed -!

          if(vfjp(i,j,k) .gt. 1.d0) then
            vfjp(i,j,k)=1.d0
          endif
          if(abs(vfjm(i,j,k)) .gt. 1.d0) then
            vfjm(i,j,k)=-1.d0
          endif
!
      enddo
    enddo
  enddo
!
  deallocate( gcov1, gcon1, gcov3a, gcon3a, ucov, ucon, bbcon, bbcov, &
              util, vcon, vcona, vcova, bcon, beta1, stat=merr)
!
   return
end subroutine calcha5j
!
!---------------------------------------------------------------------@
subroutine calcha5k(uri,gcov,gcon,vfkp,vfkm,nm0,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Calculate characteristics at cell-interface
!     (from Gammie et al. 2003)
!
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1
  integer :: nm0, merr
!
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: vfkp(is1:ie1,js1:je1,ks1:ke1), vfkm(is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: gcov1(:,:), gcon1(:,:), gcov3a(:,:), gcon3a(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcon(:), bbcov(:)  
  real(8), allocatable :: util(:), vcon(:), vcona(:), vcova(:), bcon(:), &
                          beta1(:)

  real(8) :: de, pr, roh, roe, alpha1, alpsq, gfl, bbsq, vsq, &
             cssq, vasq, omsq, vnsq, vnsq1, f1, f2k, f3, f4k
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), gcov3a(1:3,1:3), gcon3a(1:3,1:3), &
       stat=merr)
  allocate(ucov(0:3), ucon(0:3), bbcon(0:3), bbcov(0:3), stat=merr)
  allocate(util(1:3), vcon(1:3), vcona(1:3), vcova(1:3), bcon(1:3), &
           beta1(1:3), stat=merr)
!     
!=====================================================================@
!
  do k=ks1+nm0-1,ke1-nm0
    do j=js1+nm0-1,je1-nm0
      do i=is1+nm0-1,ie1-nm0
!
!- copy of primitive variables -!
        de=uri(1,i,j,k)
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
        pr=uri(5,i,j,k)
        bcon(1)=uri(7,i,j,k)
        bcon(2)=uri(8,i,j,k)
        bcon(3)=uri(9,i,j,k)
!
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
        call cal3met(gcov1,gcon1,gcov3a,gcon3a)
!
        alpha1=1./sqrt(-gcon1(0,0))
        alpsq=alpha1*alpha1
!
        beta1(1)=alpsq*gcon1(0,1)
        beta1(2)=alpsq*gcon1(0,2)
        beta1(3)=alpsq*gcon1(0,3)
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
        vcon(1)=util(1)*(1./gfl)
        vcon(2)=util(2)*(1./gfl)
        vcon(3)=util(3)*(1./gfl)
!
        bcon(1)=bcon(1)*(1./alpha1)
        bcon(2)=bcon(2)*(1./alpha1)
        bcon(3)=bcon(3)*(1./alpha1)
!        
        vsq=0.d0
        do m=1,3
          do n=1,3
            vsq=vsq+gcov3a(m,n)*vcon(m)*vcon(n)     
          enddo
        enddo
!         
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          cssq=gam*pr/roh
        elseif(ieos .eq. 1) then
          cssq=(pr/(3.0*roh))*((5.0*roh-8.0*pr)/(roh-pr))
        endif
!         
        vasq=bbsq*(1./(roh+bbsq))
        omsq=vasq+cssq*(1.0-vasq)
!
!- v^i + beta^i
        vcona(1)=vcon(1)+beta1(1)
        vcona(2)=vcon(2)+beta1(2)
        vcona(3)=vcon(3)+beta1(3)
        call lower1(vcona,vcova,gcov3a)
!
        vnsq1=vcova(1)*vcona(1)+vcova(2)*vcona(2)+vcova(3)*vcona(3)
!
        vnsq=0.d0
        do m=1,3
          do n=1,3
            vnsq=vnsq+gcov3a(m,n)*vcona(m)*vcona(n)
          enddo
        enddo
!
        f1=alpsq-vnsq1*omsq
        f2k=vcon(3)*alpsq*(1.-omsq)-beta1(3)*omsq*(alpsq-vnsq)

        f3=alpha1*sqrt(omsq)
        f4k=(alpsq-vnsq)*(gcov3a(3,3)*(alpsq-vnsq*omsq) &
                                 -(1.-omsq)*vcona(3)*vcona(3))

        vfkp(i,j,k)=(1./f1)*(f2k+f3*sqrt(f4k))
        vfkm(i,j,k)=(1./f1)*(f2k-f3*sqrt(f4k))
!
!- check wave speed exceed light speed -!

          if(vfkp(i,j,k) .gt. 1.d0) then
            vfkp(i,j,k)=1.d0
          endif
          if(abs(vfkm(i,j,k)) .gt. 1.d0) then
            vfkm(i,j,k)=-1.d0
          endif
!
      enddo
    enddo
  enddo
!
  deallocate( gcov1, gcon1, gcov3a, gcon3a, ucov, ucon, bbcon, bbcov, &
              util, vcon, vcona, vcova, bcon, beta1, stat=merr)
!
   return
end subroutine calcha5k
!
