! for Divergence-Free Magnetic Field
! Constrained Transport Schemes
!---------------------------------------------------------------------@
subroutine ct(ww,wwo,uri,nm0,nm1,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Calculate numerical flux for magnetic field by flux-CT method
!
  use pram, only : imax, jmax, kmax, nv, ict
  implicit none
  
  integer :: nm0, nm1, is1, ie1, js1, je1, ks1, ke1

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: ww(3,nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: wwo(3,nv,is1:ie1,js1:je1,ks1:ke1)
!
  if(ict .eq. 1) then
    call fluxct(ww,nm0,nm1,is1,ie1,js1,je1,ks1,ke1)
  elseif(ict .eq. 2) then
    call mflct(ww,wwo,nm0,nm1,is1,ie1,js1,je1,ks1,ke1)
  elseif(ict .eq. 3) then
    call upflct(ww,wwo,uri,nm0,nm1,is1,ie1,js1,je1,ks1,ke1)
  endif
!
  return
end subroutine ct
!
!---------------------------------------------------------------------@
subroutine fluxct(ww,nm0,nm1,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Calculate numerical flux for magnetic field by flux-CT method
!
  use pram, only : imax, jmax, kmax, nv
  implicit none
  
  integer :: i, j, k
  integer :: nm0, nm1, is1, ie1, js1, je1, ks1, ke1
  integer :: merr

  real(8) :: ww(3,nv,is1:ie1,js1:je1,ks1:ke1)
  
  real(8), allocatable :: fyx(:,:,:), fzx(:,:,:), &
           fxy(:,:,:), fzy(:,:,:), fxz(:,:,:), fyz(:,:,:)

  allocate( fyx(is1:ie1,js1:je1,ks1:ke1), fzx(is1:ie1,js1:je1,ks1:ke1), &
            fxy(is1:ie1,js1:je1,ks1:ke1), fzy(is1:ie1,js1:je1,ks1:ke1), &
            fxz(is1:ie1,js1:je1,ks1:ke1), fyz(is1:ie1,js1:je1,ks1:ke1), &
            stat=merr )
!
!=====================================================================@
!
  do k=ks1+nm0-1,ke1-nm0
    do j=js1+nm0-1,je1-nm0
      do i=is1+nm0-1,ie1-nm0

        fyx(i,j,k)=ww(2,7,i,j,k)
        fzx(i,j,k)=ww(3,7,i,j,k)
        fxy(i,j,k)=ww(1,8,i,j,k)
        fzy(i,j,k)=ww(3,8,i,j,k)
        fxz(i,j,k)=ww(1,9,i,j,k)
        fyz(i,j,k)=ww(2,9,i,j,k)

      enddo
    enddo
  enddo
! 
  do k=ks1+nm1-1,ke1-nm1
    do j=js1+nm1-1,je1-nm1
      do i=is1+nm1-1,ie1-nm1

!           fbyx=(1.0/8.0)*
!     &          (2.0*fyx(i,j,k)+fyx(i,j,k+1)+fyx(i,j,k-1)
!     &          -fxy(i,j,k)-fxy(i,j+1,k)
!     &          -fxy(i-1,j,k)-fxy(i-1,j+1,k))
!           fbzx=(1.0/8.0)*
!     &          (2.0*fzx(i,j,k)+fzx(i,j+1,k)+fzx(i,j-1,k)
!     &          -fxz(i,j,k)-fxz(i,j,k+1)
!     &          -fxz(i-1,j,k)-fxz(i-1,j,k+1))
!           fbxy=(1.0/8.0)*
!     &          (2.0*fxy(i,j,k)+fxy(i,j,k+1)+fxy(i,j,k-1)
!     &          -fyx(i,j,k)-fyx(i+1,j,k)
!     &          -fyx(i,j-1,k)-fyx(i+1,j-1,k))
!           fbzy=(1.0/8.0)*
!     &          (2.0*fzy(i,j,k)+fzy(i+1,j,k)+fzy(i-1,j,k)
!     &          -fyz(i,j,k)-fyz(i,j,k+1)
!     &          -fyz(i,j-1,k)-fyz(i,j-1,k+1))
!           fbxz=(1.0/8.0)*
!     &          (2.0*fxz(i,j,k)+fxz(i,j+1,k)+fxz(i,j-1,k)
!     &          -fzx(i,j,k)-fzx(i+1,j,k)
!     &          -fzx(i,j,k-1)-fzx(i+1,j,k-1))
!           fbyz=(1.0/8.0)*
!     &          (2.0*fyz(i,j,k)+fyz(i+1,j,k)+fyz(i-1,j,k)
!     &          -fzy(i,j,k)-fzy(i,j+1,k)
!     &          -fzy(i,j,k-1)-fzy(i,j+1,k-1))

        ww(2,7,i,j,k)=(1.0/8.0)*(2.0*fyx(i,j,k)+fyx(i+1,j,k)+fyx(i-1,j,k) &
                     -fxy(i,j,k)-fxy(i,j+1,k)-fxy(i-1,j,k)-fxy(i-1,j+1,k))
        ww(3,7,i,j,k)=(1.0/8.0)*(2.0*fzx(i,j,k)+fzx(i+1,j,k)+fzx(i-1,j,k) &
                     -fxz(i,j,k)-fxz(i,j,k+1)-fxz(i-1,j,k)-fxz(i-1,j,k+1))
        ww(1,8,i,j,k)=(1.0/8.0)*(2.0*fxy(i,j,k)+fxy(i,j+1,k)+fxy(i,j-1,k) &
                     -fyx(i,j,k)-fyx(i+1,j,k)-fyx(i,j-1,k)-fyx(i+1,j-1,k))
        ww(3,8,i,j,k)=(1.0/8.0)*(2.0*fzy(i,j,k)+fzy(i,j+1,k)+fzy(i,j-1,k) &
                     -fyz(i,j,k)-fyz(i,j,k+1)-fyz(i,j-1,k)-fyz(i,j-1,k+1))
        ww(1,9,i,j,k)=(1.0/8.0)*(2.0*fxz(i,j,k)+fxz(i,j,k+1)+fxz(i,j,k-1) &
                     -fzx(i,j,k)-fzx(i+1,j,k)-fzx(i,j,k-1)-fzx(i+1,j,k-1))
        ww(2,9,i,j,k)=(1.0/8.0)*(2.0*fyz(i,j,k)+fyz(i,j,k+1)+fyz(i,j,k-1) &
                     -fzy(i,j,k)-fzy(i,j+1,k)-fzy(i,j,k-1)-fzy(i,j+1,k-1))
  
      enddo
    enddo
  enddo

  deallocate(fyx, fzx, fxy, fzy, fxz, fyz, stat=merr)

  return
end subroutine fluxct
!
!---------------------------------------------------------------------@
subroutine mflct(ww,wwo,nm0,nm1,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Calculate numerical flux for magnetic field by flux-CT method
!
  use pram, only : imax, jmax, kmax, nv
  implicit none
  
  integer :: i, j, k
  integer :: nm0, nm1, is1, ie1, js1, je1, ks1, ke1
  integer :: merr

  real(8) :: ww(3,nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: wwo(3,nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: fyx(:,:,:), fzx(:,:,:), &
           fxy(:,:,:), fzy(:,:,:), fxz(:,:,:), fyz(:,:,:), &
           exc(:,:,:), eyc(:,:,:), ezc(:,:,:)
!
  real(8) :: fbyx, fbzx, fbxy, fbzy, fbxz, fbyz, &
             fhyx, fhzx, fhxy, fhzy, fhxz, fhyz

  allocate( fyx(is1:ie1,js1:je1,ks1:ke1), fzx(is1:ie1,js1:je1,ks1:ke1), &
            fxy(is1:ie1,js1:je1,ks1:ke1), fzy(is1:ie1,js1:je1,ks1:ke1), &
            fxz(is1:ie1,js1:je1,ks1:ke1), fyz(is1:ie1,js1:je1,ks1:ke1), &
            exc(is1:ie1,js1:je1,ks1:ke1), eyc(is1:ie1,js1:je1,ks1:ke1), &
            ezc(is1:ie1,js1:je1,ks1:ke1), stat=merr )
!
!=====================================================================@
!
  do k=ks1+nm0-1,ke1-nm0
    do j=js1+nm0-1,je1-nm0
      do i=is1+nm0-1,ie1-nm0

        fyx(i,j,k)=ww(2,7,i,j,k)
        fzx(i,j,k)=ww(3,7,i,j,k)
        fxy(i,j,k)=ww(1,8,i,j,k)
        fzy(i,j,k)=ww(3,8,i,j,k)
        fxz(i,j,k)=ww(1,9,i,j,k)
        fyz(i,j,k)=ww(2,9,i,j,k)
           
!        ezc(i,j,k)=wwo(2,7,i,j,k)
!        exc(i,j,k)=wwo(3,8,i,j,k)
!        eyc(i,j,k)=wwo(1,9,i,j,k)

      enddo
    enddo
  enddo
!
  do k=ks1+nm0,ke1-nm0
    do j=js1+nm0,je1-nm0
      do i=is1+nm0,ie1-nm0
           
        ezc(i,j,k)=0.25*(-fxy(i-1,j,k)-fxy(i,j,k)+fyx(i,j-1,k)+fyx(i,j,k))
        exc(i,j,k)=0.25*(-fyz(i,j-1,k)-fyz(i,j,k)+fzy(i,j,k-1)+fzy(i,j,k))
        eyc(i,j,k)=0.25*(-fzx(i,j,k-1)-fzx(i,j,k)+fxz(i-1,j,k)+fxz(i,j,k))
           
      enddo
    enddo
  enddo

! Memo: nm1=nm0+2

  do k=ks1+nm1-1,ke1-nm1
    do j=js1+nm1-1,je1-nm1
      do i=is1+nm1-1,ie1-nm1
!
! orginial flux ct *2
!
     
        fbyx=(1.0/4.0)*(2.0*fyx(i,j,k)+fyx(i+1,j,k)+fyx(i-1,j,k) &
             -fxy(i,j,k)-fxy(i,j+1,k)-fxy(i-1,j,k)-fxy(i-1,j+1,k))
        fbzx=(1.0/4.0)*(2.0*fzx(i,j,k)+fzx(i+1,j,k)+fzx(i-1,j,k) &
             -fxz(i,j,k)-fxz(i,j,k+1)-fxz(i-1,j,k)-fxz(i-1,j,k+1))
        fbxy=(1.0/4.0)*(2.0*fxy(i,j,k)+fxy(i,j+1,k)+fxy(i,j-1,k) &
             -fyx(i,j,k)-fyx(i+1,j,k)-fyx(i,j-1,k)-fyx(i+1,j-1,k))
        fbzy=(1.0/4.0)*(2.0*fzy(i,j,k)+fzy(i,j+1,k)+fzy(i,j-1,k) &
             -fyz(i,j,k)-fyz(i,j,k+1)-fyz(i,j-1,k)-fyz(i,j-1,k+1))
        fbxz=(1.0/4.0)*(2.0*fxz(i,j,k)+fxz(i,j,k+1)+fxz(i,j,k-1) &
             -fzx(i,j,k)-fzx(i+1,j,k)-fzx(i,j,k-1)-fzx(i+1,j,k-1))
        fbyz=(1.0/4.0)*(2.0*fyz(i,j,k)+fyz(i,j,k+1)+fyz(i,j,k-1) &
             -fzy(i,j,k)-fzy(i,j+1,k)-fzy(i,j,k-1)-fzy(i,j+1,k-1))
!
! Modified
!
        fhyx=-(1.0/8.0)*(2.*ezc(i,j,k)+2.*ezc(i,j+1,k)+ezc(i+1,j,k) &
             +ezc(i+1,j+1,k)+ezc(i-1,j,k)+ezc(i-1,j+1,k))
        fhzx=(1.0/8.0)*(2.*eyc(i,j,k)+2.*eyc(i,j,k+1)+eyc(i+1,j,k) &
             +eyc(i+1,j,k+1)+eyc(i-1,j,k)+eyc(i-1,j,k+1))
        fhxy=(1.0/8.0)*(2.*ezc(i,j,k)+2.*ezc(i+1,j,k)+ezc(i,j+1,k) &
             +ezc(i+1,j+1,k)+ezc(i,j-1,k)+ezc(i+1,j-1,k))
        fhzy=-(1.0/8.0)*(2.*exc(i,j,k)+2.*exc(i,j,k+1)+exc(i,j+1,k) &
             +exc(i,j+1,k+1)+exc(i,j-1,k)+exc(i,j-1,k+1))
        fhxz=-(1.0/8.0)*(2.*eyc(i,j,k)+2.*eyc(i+1,j,k)+eyc(i,j,k+1) &
             +eyc(i+1,j,k+1)+eyc(i,j,k-1)+eyc(i+1,j,k-1))
        fhyz=(1.0/8.0)*(2.*exc(i,j,k)+2.*exc(i,j+1,k)+exc(i,j,k+1) &
             +exc(i,j+1,k+1)+exc(i,j,k-1)+exc(i,j+1,k-1))

        ww(2,7,i,j,k)=fbyx+fhyx
        ww(3,7,i,j,k)=fbzx+fhzx
        ww(1,8,i,j,k)=fbxy+fhxy
        ww(3,8,i,j,k)=fbzy+fhzy
        ww(1,9,i,j,k)=fbxz+fhxz
        ww(2,9,i,j,k)=fbyz+fhyz

      enddo
    enddo
  enddo

  deallocate( fyx, fzx, fxy, fzy, fxz, fyz, exc, eyc, ezc, stat=merr )
!
  return
end subroutine mflct
!
!---------------------------------------------------------------------@
subroutine upflct(ww,wwo,uri,nm0,nm1,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Calculate numerical flux for magnetic field by flux-CT method
!
  use pram, only : imax, jmax, kmax, nv
  implicit none
  
  integer :: i, j, k
  integer :: nm0, nm1, nm2, nm3, is1, ie1, js1, je1, ks1, ke1
  integer :: merr

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: ww(3,nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: wwo(3,nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: fyx(:,:,:), fzx(:,:,:), &
           fxy(:,:,:), fzy(:,:,:), fxz(:,:,:), fyz(:,:,:)
!
  real(8), allocatable :: exc(:,:,:), eyc(:,:,:), ezc(:,:,:), &
           exb(:,:,:), eyb(:,:,:), ezb(:,:,:)
      
  real(8), allocatable :: fbyx(:,:,:), fbzx(:,:,:), &
           fbxy(:,:,:), fbzy(:,:,:), fbxz(:,:,:), fbyz(:,:,:)
      
  real(8), allocatable :: dyez(:,:,:), dxez(:,:,:), &
           dyex(:,:,:), dzex(:,:,:), dxey(:,:,:), dzey(:,:,:)
!
  real(8) :: ezh, ezhmi, ezhmj, exh, exhmj, exhmk, eyh, eyhmi, eyhmk
!
  allocate ( fyx(is1:ie1,js1:je1,ks1:ke1), fzx(is1:ie1,js1:je1,ks1:ke1), &
             fxy(is1:ie1,js1:je1,ks1:ke1), fzy(is1:ie1,js1:je1,ks1:ke1), &
             fxz(is1:ie1,js1:je1,ks1:ke1), fyz(is1:ie1,js1:je1,ks1:ke1), &
             exc(is1:ie1,js1:je1,ks1:ke1), eyc(is1:ie1,js1:je1,ks1:ke1), &
             ezc(is1:ie1,js1:je1,ks1:ke1), &
             exb(is1:ie1,js1:je1,ks1:ke1), eyb(is1:ie1,js1:je1,ks1:ke1), &
             ezb(is1:ie1,js1:je1,ks1:ke1), &
             fbyx(is1:ie1,js1:je1,ks1:ke1), fbzx(is1:ie1,js1:je1,ks1:ke1), &
             fbxy(is1:ie1,js1:je1,ks1:ke1), fbzy(is1:ie1,js1:je1,ks1:ke1), &
             fbxz(is1:ie1,js1:je1,ks1:ke1), fbyz(is1:ie1,js1:je1,ks1:ke1), &
             dyez(is1:ie1,js1:je1,ks1:ke1), dxez(is1:ie1,js1:je1,ks1:ke1), &
             dyex(is1:ie1,js1:je1,ks1:ke1), dzex(is1:ie1,js1:je1,ks1:ke1), &
             dxey(is1:ie1,js1:je1,ks1:ke1), dzey(is1:ie1,js1:je1,ks1:ke1), &
             stat=merr )
!
!=====================================================================@

  do k=ks1+nm0-1,ke1-nm0
    do j=js1+nm0-1,je1-nm0
      do i=is1+nm0-1,ie1-nm0

        fyx(i,j,k)=ww(2,7,i,j,k)
        fzx(i,j,k)=ww(3,7,i,j,k)
        fxy(i,j,k)=ww(1,8,i,j,k)
        fzy(i,j,k)=ww(3,8,i,j,k)
        fxz(i,j,k)=ww(1,9,i,j,k)
        fyz(i,j,k)=ww(2,9,i,j,k)
           
!       ezc(i,j,k)=wwo(2,7,i,j,k)
!       exc(i,j,k)=wwo(3,8,i,j,k)
!       eyc(i,j,k)=wwo(1,9,i,j,k)

      enddo
    enddo
  enddo
!
! cal E^z(i,j,k), E^z(i+1/2,j+1/2,k), E^z(i,j+1/2,k) ...
!
  nm2=nm0+1

  do k=ks1+nm2-1,ke1-nm2
    do j=js1+nm2-1,je1-nm2
      do i=is1+nm2-1,ie1-nm2
           
         ezc(i,j,k)=0.25*(-fxy(i-1,j,k)-fxy(i,j,k)+fyx(i,j-1,k)+fyx(i,j,k))
         exc(i,j,k)=0.25*(-fyz(i,j-1,k)-fyz(i,j,k)+fzy(i,j,k-1)+fzy(i,j,k))
         eyc(i,j,k)=0.25*(-fzx(i,j,k-1)-fzx(i,j,k)+fxz(i-1,j,k)+fxz(i,j,k))
           
         ezb(i,j,k)=0.25*(fyx(i,j,k)+fyx(i+1,j,k)-fxy(i,j,k)-fxy(i,j+1,k))
         exb(i,j,k)=0.25*(fzy(i,j,k)+fzy(i,j+1,k)-fyz(i,j,k)-fyz(i,j,k+1))
         eyb(i,j,k)=0.25*(fxz(i,j,k)+fxz(i+1,j,k)-fzx(i,j,k)-fzx(i,j+1,k))
     
         fbyx(i,j,k)=(1.0/8.0)*(2.0*fyx(i,j,k)+fyx(i,j,k+1)+fyx(i,j,k-1) &
                    -fxy(i,j,k)-fxy(i,j+1,k)-fxy(i-1,j,k)-fxy(i-1,j+1,k))
         fbzx(i,j,k)=(1.0/8.0)*(2.0*fzx(i,j,k)+fzx(i,j+1,k)+fzx(i,j-1,k) &
                    -fxz(i,j,k)-fxz(i,j,k+1)-fxz(i-1,j,k)-fxz(i-1,j,k+1))
         fbxy(i,j,k)=(1.0/8.0)*(2.0*fxy(i,j,k)+fxy(i,j,k+1)+fxy(i,j,k-1) &
                    -fyx(i,j,k)-fyx(i+1,j,k)-fyx(i,j-1,k)-fyx(i+1,j-1,k))
         fbzy(i,j,k)=(1.0/8.0)*(2.0*fzy(i,j,k)+fzy(i+1,j,k)+fzy(i-1,j,k) &
                    -fyz(i,j,k)-fyz(i,j,k+1)-fyz(i,j-1,k)-fyz(i,j-1,k+1))
         fbxz(i,j,k)=(1.0/8.0)*(2.0*fxz(i,j,k)+fxz(i,j+1,k)+fxz(i,j-1,k) &
                    -fzx(i,j,k)-fzx(i+1,j,k)-fzx(i,j,k-1)-fzx(i+1,j,k-1))
         fbyz(i,j,k)=(1.0/8.0)*(2.0*fyz(i,j,k)+fyz(i+1,j,k)+fyz(i-1,j,k) &
                    -fzy(i,j,k)-fzy(i,j+1,k)-fzy(i,j,k-1)-fzy(i,j+1,k-1))
     
      enddo
    enddo
  enddo
!
! cal dyE^z(i+1/2,j+1/4,k)
!
  nm3=nm0+2
     
  do k=ks1+nm2-1,ke1-nm3
    do j=js1+nm2-1,je1-nm3
      do i=is1+nm2-1,ie1-nm3

        if(uri(2,i,j,k) .gt. 0.0) then
          dyez(i,j,k)=fbyx(i,j,k)-ezc(i,j,k)
        elseif(uri(2,i,j,k) .lt. 0.0) then
          dyez(i,j,k)=fbyx(i+1,j,k)-ezc(i+1,j,k)
        else
          dyez(i,j,k)=0.5*(fbyx(i,j,k)+fbyx(i+1,j,k)-ezc(i,j,k)-ezc(i+1,j,k))
        endif
         
        if(uri(3,i,j,k) .gt. 0.0) then
          dxez(i,j,k)=fbxy(i,j,k)-ezc(i,j,k)
        elseif(uri(2,i,j,k) .lt. 0.0) then
          dxez(i,j,k)=fbxy(i,j+1,k)-ezc(i,j+1,k)
        else
          dxez(i,j,k)=0.5*(fbxy(i,j,k)+fbxy(i,j+1,k)-ezc(i,j,k)-ezc(i,j+1,k))
        endif
         
        if(uri(4,i,j,k) .gt. 0.0) then
          dyex(i,j,k)=fbyz(i,j,k)-exc(i,j,k)
        elseif(uri(2,i,j,k) .lt. 0.0) then
          dyex(i,j,k)=fbyz(i,j,k+1)-exc(i,j,k+1)
        else
          dyex(i,j,k)=0.5*(fbyz(i,j,k)+fbyz(i,j,k+1)-exc(i,j,k)-exc(i,j,k+1))
        endif
         
        if(uri(3,i,j,k) .gt. 0.0) then
          dzex(i,j,k)=fbzy(i,j,k)-exc(i,j,k)
        elseif(uri(2,i,j,k) .lt. 0.0) then
          dzex(i,j,k)=fbzy(i,j+1,k)-exc(i,j+1,k)
        else
          dzex(i,j,k)=0.5*(fbzy(i,j,k)+fbzy(i,j+1,k)-exc(i,j,k)-exc(i,j+1,k))
        endif
         
        if(uri(4,i,j,k) .gt. 0.0) then
          dxey(i,j,k)=fbxz(i,j,k)-eyc(i,j,k)
        elseif(uri(2,i,j,k) .lt. 0.0) then
          dxey(i,j,k)=fbxz(i,j,k+1)-eyc(i,j,k+1)
        else
          dxey(i,j,k)=0.5*(fbxz(i,j,k)+fbxz(i,j,k+1)-eyc(i,j,k)-eyc(i,j,k+1))
        endif
         
        if(uri(2,i,j,k) .gt. 0.0) then
          dzey(i,j,k)=fbzx(i,j,k)-eyc(i,j,k)
        elseif(uri(2,i,j,k) .lt. 0.0) then
          dzey(i,j,k)=fbzx(i+1,j,k)-eyc(i,j,k)
        else
          dzey(i,j,k)=0.5*(fbzx(i,j,k)+fbzx(i+1,j,k)-eyc(i,j,k)-eyc(i+1,j,k))
        endif
         
      enddo
    enddo
  enddo
!
! cal E^hat z(i+1/2,j+1/2,k), flux
!
!      nm1=nm0+3

  do k=ks1+nm1-1,ke1-nm1
    do j=js1+nm1-1,je1-nm1
      do i=is1+nm1-1,ie1-nm1
          
        ezh=ezb(i,j,k)+0.25*(dyez(i,j,k)-dyez(i,j+1,k)) &
           +0.25*(dxez(i,j,k)-dxez(i+1,j,k))
          
        ezhmi=ezb(i-1,j,k)+0.25*(dyez(i-1,j,k)-dyez(i-1,j+1,k)) &
             +0.25*(dxez(i-1,j,k)-dxez(i,j,k))
        ezhmj=ezb(i,j-1,k)+0.25*(dyez(i,j-1,k)-dyez(i,j,k)) &
             +0.25*(dxez(i,j-1,k)-dxez(i+1,j-1,k))
     
        exh=exb(i,j,k)+0.25*(dyex(i,j,k)-dyex(i,j+1,k)) &
           +0.25*(dzex(i,j,k)-dzex(i,j,k+1))
          
        exhmj=exb(i,j-1,k)+0.25*(dyex(i,j-1,k)-dyex(i,j,k)) &
             +0.25*(dzex(i,j-1,k)-dzex(i,j-1,k+1))
        exhmk=exb(i,j,k-1)+0.25*(dyex(i,j,k-1)-dyex(i,j+1,k-1)) &
             +0.25*(dzex(i,j,k-1)-dzex(i,j,k))
          
        eyh=eyb(i,j,k)+0.25*(dxey(i,j,k)-dxey(i+1,j,k)) &
           +0.25*(dzey(i,j,k)-dzey(i,j,k+1))
     
        eyhmi=eyb(i-1,j,k)+0.25*(dxey(i-1,j,k)-dxey(i,j,k)) &
             +0.25*(dzey(i-1,j,k)-dzey(i-1,j,k+1))
        eyhmk=eyb(i,j,k-1)+0.25*(dxey(i,j,k-1)-dxey(i+1,j,k-1)) &
             +0.25*(dzey(i,j,k-1)-dzey(i,j,k))
          
        ww(2,7,i,j,k)=0.5*(ezh+ezhmi)
        ww(3,7,i,j,k)=0.5*(-eyh-eyhmi)
        ww(1,8,i,j,k)=0.5*(-ezh-ezhmj)
        ww(3,8,i,j,k)=0.5*(exh+exhmj)
        ww(1,9,i,j,k)=0.5*(eyh+eyhmk)
        ww(2,9,i,j,k)=0.5*(-exh-exhmk)

      enddo
    enddo        
  enddo

  deallocate ( fyx, fzx, fxy, fzy, fxz, fyz, &
               exc, eyc, ezc, exb, eyb, ezb, &
               fbyx, fbzx, fbxy, fbzy, fbxz, fbyz, &
               dyez, dxez, dyex, dzex, dxey, dzey, stat=merr )

  return
end subroutine upflct
!
