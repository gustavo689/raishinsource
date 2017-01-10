!--------------------------------------------------------------------
subroutine caluu1(uri,uu,gcov,gcon,detg,x1,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos, iflag  
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1 
  integer :: merr
       
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)

  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: delta(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:)
  real(8), allocatable :: util(:), bcon(:), vcon(:)  

  real(8) :: x1(imax)      

  real(8) :: de,  pr, roh, roe, gfl, bbsq
  real(8) :: alpha1
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), delta(0:3,0:3), stat=merr)
  allocate(util(1:3), bcon(1:3), vcon(1:3), stat=merr)
  allocate(ucov(0:3), ucon(0:3), bbcov(0:3), bbcon(0:3), stat=merr)
!
!=====================================================================
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!
!- copy of primitive variables -!
        de=uri(1,i,j,k)
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
!        v1=uri(2,i,j,k)
!        v2=uri(3,i,j,k)
!        v3=uri(4,i,j,k)
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
        uu(1,i,j,k)=de*ucon(0) !- relativistic mass -!
        uu(2,i,j,k)=(roh+bbsq)*ucon(0)*ucov(1)+(pr+0.5*bbsq)*delta(0,1) &
                   -bbcon(0)*bbcov(1) !- momentum density in i-drec. -!
        uu(3,i,j,k)=(roh+bbsq)*ucon(0)*ucov(2)+(pr+0.5*bbsq)*delta(0,2) &
                   -bbcon(0)*bbcov(2) !- momentum density in j-direc. -!
        uu(4,i,j,k)=(roh+bbsq)*ucon(0)*ucov(3)+(pr+0.5*bbsq)*delta(0,3) &
                   -bbcon(0)*bbcov(3) !- momentum density in k-direc. -!
!        uu(5,i,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
!                   -bbcon(0)*bbcov(0)-de*ucon(0) !- total energy -!
        uu(5,i,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
                   -bbcon(0)*bbcov(0)+de*ucon(0) !- total energy -!
!        uu(5,i,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
!                   -bbcon(0)*bbcov(0) !- total energy -!
           
        uu(6,i,j,k)=de*(pr/(de**(gam)))*ucon(0)         
           
        uu(7,i,j,k)=bcon(1) !- contravariant B-field in i-direc. -!
        uu(8,i,j,k)=bcon(2) !- contravariant B-field in j-direc. -!
        uu(9,i,j,k)=bcon(3) !- contravariant B-field in k-direc. -!
!
!- determinant g * U -!
        do n=1,9
          uu(n,i,j,k)=detg(i,j,k)*uu(n,i,j,k)  
        enddo
!
        alpha1=1./sqrt(-gcon1(0,0))
!        if(i .eq. 128 .and. j .eq. 4 .and. k .eq. 4) then
!           write(*,*) 'uu(1), uu(2)', uu(1,i,j,k), uu(2,i,j,k)
!           write(*,*) 'uu(3), uu(4)', uu(3,i,j,k), uu(4,i,j,k)
!           write(*,*) 'uu(5)', (alpha1*gcon1(0,0)*(uu(5,i,j,k)-uu(1,i,j,k)))
!        endif
!                
      enddo
    enddo
  enddo
!
  deallocate( gcov1, gcon1, delta, util, bcon, vcon, ucov, ucon, &
              bbcov, bbcon, stat=merr)
      
  return
end subroutine caluu1
!
!--------------------------------------------------------------------
subroutine caluu2a(uri,uu,gcov,gcon,detg,x1,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos, iflag  
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1 
  integer :: merr
       
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)

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
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1

        de=uri(1,i,j,k)
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
!        v1=uri(2,i,j,k)
!        v2=uri(3,i,j,k)
!        v3=uri(4,i,j,k)
        pr=uri(5,i,j,k)
        bcon(1)=uri(7,i,j,k)
        bcon(2)=uri(8,i,j,k)
        bcon(3)=uri(9,i,j,k)
!        
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          roh=de+(gam/(gam-1.0))*pr
        elseif(ieos .eq. 1) then
          roh=(5./2.)*pr+sqrt((9./4.)*pr**2+de**2)
        elseif(ieos .eq. 2) then
          roe=(3./2.)*(pr+((3.*pr**2)/(2.0*de+sqrt(2.*pr**2+4.*de**2)) ))
          roh=de+roe+pr
        endif
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
        call caldelta(delta,gcov1,gcon1)
         
        uu(1,i,j,k)=de*ucon(0)
        uu(2,i,j,k)=(roh+bbsq)*ucon(0)*ucov(1)+(pr+0.5*bbsq)*delta(0,1) &
                   -bbcon(0)*bbcov(1)
        uu(3,i,j,k)=(roh+bbsq)*ucon(0)*ucov(2)+(pr+0.5*bbsq)*delta(0,2) &
                   -bbcon(0)*bbcov(2)
        uu(4,i,j,k)=(roh+bbsq)*ucon(0)*ucov(3)+(pr+0.5*bbsq)*delta(0,3) &
                   -bbcon(0)*bbcov(3)
!        uu(5,i,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
!                   -bbcon(0)*bbcov(0)-de*ucon(0)
        uu(5,i,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
                   -bbcon(0)*bbcov(0)+de*ucon(0)
!        uu(5,i,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
!                   -bbcon(0)*bbcov(0)
!
        uu(6,i,j,k)=de*(pr/(de**(gam)))*ucon(0) 

        uu(7,i,j,k)=bcon(1)
        uu(8,i,j,k)=bcon(2)
        uu(9,i,j,k)=bcon(3)
!
        do n=1,9
          uu(n,i,j,k)=detg(i,j,k)*uu(n,i,j,k)
        enddo
!
      enddo
    enddo
  enddo
!
  deallocate( gcov1, gcon1, delta, util, bcon, vcon, ucov, ucon, &
              bbcov, bbcon, stat=merr)
      
  return
end subroutine caluu2a
!
!--------------------------------------------------------------------
subroutine caluu3(uri,uu,gcov,gcon,detg,x1,nm0,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos, iflag  
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1
  integer :: nm0 
  integer :: merr
       
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)

  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: delta(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:)
  real(8), allocatable :: util(:), bcon(:), vcon(:)  

  real(8) :: x1(imax)      

  real(8) :: de,  pr, roh, roe, gfl, bbsq
  real(8) :: alpha1
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), delta(0:3,0:3), stat=merr)
  allocate(util(1:3), bcon(1:3), vcon(1:3), stat=merr)
  allocate(ucov(0:3), ucon(0:3), bbcov(0:3), bbcon(0:3), stat=merr)
!
!=====================================================================
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
        uu(1,i,j,k)=de*ucon(0) !- relativistic mass -!
        uu(2,i,j,k)=(roh+bbsq)*ucon(0)*ucov(1)+(pr+0.5*bbsq)*delta(0,1) &
                   -bbcon(0)*bbcov(1) !- momentum density in i-direc -!
        uu(3,i,j,k)=(roh+bbsq)*ucon(0)*ucov(2)+(pr+0.5*bbsq)*delta(0,2) &
                   -bbcon(0)*bbcov(2) !- momentum density in j-direc -!
        uu(4,i,j,k)=(roh+bbsq)*ucon(0)*ucov(3)+(pr+0.5*bbsq)*delta(0,3) &
                   -bbcon(0)*bbcov(3) !- momentum density in k-direc -!
!        uu(5,i,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
!                   -bbcon(0)*bbcov(0)-de*ucon(0) !- total energy -!
        uu(5,i,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
                   -bbcon(0)*bbcov(0)+de*ucon(0) !- total energy -!
!        uu(5,i,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
!                   -bbcon(0)*bbcov(0) !- total energy -!

        uu(6,i,j,k)=de*(pr/(de**(gam)))*ucon(0) 
           
        uu(7,i,j,k)=bcon(1) !- contravariant B-field in i-direc -!
        uu(8,i,j,k)=bcon(2) !- contravariant B-field in j-direc -!
        uu(9,i,j,k)=bcon(3) !- contravariant B-field in k-direc -!
!
        do n=1,9
!- determinant g * U -!
          uu(n,i,j,k)=detg(i,j,k)*uu(n,i,j,k)
        enddo
!
        alpha1=1./sqrt(-gcon1(0,0))
!        
      enddo
    enddo
  enddo
!
  deallocate( gcov1, gcon1, delta, util, bcon, vcon, ucov, ucon, &
              bbcov, bbcon, stat=merr)
!      
  return
end subroutine caluu3
!
!--------------------------------------------------------------------
subroutine caluu3a(uri,uu,gcov,gcon,detg,x1,nm1,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos, iflag  
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1
  integer :: nm1 
  integer :: merr
       
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)

  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: delta(:,:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:)
  real(8), allocatable :: util(:), bcon(:), vcon(:)  

  real(8) :: x1(imax)      

  real(8) :: de,  pr, roh, roe, gfl, bbsq
  real(8) :: alpha1
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), delta(0:3,0:3), stat=merr)
  allocate(util(1:3), bcon(1:3), vcon(1:3), stat=merr)
  allocate(ucov(0:3), ucon(0:3), bbcov(0:3), bbcon(0:3), stat=merr)
!
!=====================================================================
!
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
      do i=is1+nm1,ie1-nm1
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
        uu(1,i,j,k)=de*ucon(0) !- relativistic mass -!
        uu(2,i,j,k)=(roh+bbsq)*ucon(0)*ucov(1)+(pr+0.5*bbsq)*delta(0,1) &
                   -bbcon(0)*bbcov(1) !- momentum density in i-direc -!
        uu(3,i,j,k)=(roh+bbsq)*ucon(0)*ucov(2)+(pr+0.5*bbsq)*delta(0,2) &
                   -bbcon(0)*bbcov(2) !- momentum density in j-direc -!
        uu(4,i,j,k)=(roh+bbsq)*ucon(0)*ucov(3)+(pr+0.5*bbsq)*delta(0,3) &
                   -bbcon(0)*bbcov(3) !- momentum density in k-direc -!
!        uu(5,i,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
!                   -bbcon(0)*bbcov(0)-de*ucon(0) !- total energy -!
        uu(5,i,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
                   -bbcon(0)*bbcov(0)+de*ucon(0) !- total energy -!
!        uu(5,i,j,k)=(roh+bbsq)*ucon(0)*ucov(0)+(pr+0.5*bbsq)*delta(0,0) &
!                   -bbcon(0)*bbcov(0) !- total energy -!
!
        uu(6,i,j,k)=de*(pr/(de**(gam)))*ucon(0)        
!           
        uu(7,i,j,k)=bcon(1) !- contravariant B-field in i-direc -!
        uu(8,i,j,k)=bcon(2) !- contravariant B-field in j-direc -!
        uu(9,i,j,k)=bcon(3) !- contravariant B-field in k-direc -!
!
        do n=1,9
!- determinant g * U -!
          uu(n,i,j,k)=detg(i,j,k)*uu(n,i,j,k)
        enddo
!
        alpha1=1./sqrt(-gcon1(0,0))
!        
      enddo
    enddo
  enddo
!
  deallocate( gcov1, gcon1, delta, util, bcon, vcon, ucov, ucon, &
              bbcov, bbcon, stat=merr)
!      
  return
end subroutine caluu3a
!
