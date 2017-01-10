!---------------------------------------------------------------------@
subroutine rec2(uri,x1,x2,x3,x1a,x2a,x3a,dx1,dx2,dx3,dx1b,dx2b,dx3b, &
                x1b,x2b,x3b,x1ab,x2ab,x3ab, & 
                uriir,urijr,urikr,uriil,urijl,urikl, &
                uuir,uujr,uukr,uuil,uujl,uukl, &
                wwir,wwjr,wwkr,wwil,wwjl,wwkl, &
                gcovi,gcovj,gcovk,gconi,gconj,gconk, &
                detgi,detgj,detgk,nm0,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Reconstruction step 
!     Calculate cell-interface variables(uuil, uuir, uril, urir, wwl, wwr) 
!     from cell-center variables (uu, uri)
!
  use pram, only : imax, jmax, kmax, nv, irec
  implicit none
!
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)
      
  real(8) :: uuir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uuil(nv,is1:ie1,js1:je1,ks1:ke1), &
             uujr(nv,is1:ie1,js1:je1,ks1:ke1), &
             uujl(nv,is1:ie1,js1:je1,ks1:ke1), &
             uukr(nv,is1:ie1,js1:je1,ks1:ke1), &
             uukl(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: wwir(nv,is1:ie1,js1:je1,ks1:ke1), &
             wwil(nv,is1:ie1,js1:je1,ks1:ke1), &
             wwjr(nv,is1:ie1,js1:je1,ks1:ke1), &
             wwjl(nv,is1:ie1,js1:je1,ks1:ke1), &
             wwkr(nv,is1:ie1,js1:je1,ks1:ke1), &
             wwkl(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: gcovi(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gcovj(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gcovk(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: gconi(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gconj(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gconk(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: detgi(is1:ie1,js1:je1,ks1:ke1), & 
             detgj(is1:ie1,js1:je1,ks1:ke1), & 
             detgk(is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: x1(imax), x2(jmax), x3(kmax), & !- cell center position -!
             x1a(imax), x2a(jmax), x3a(kmax), & !- cell boundary position -!
                                                !- x1a(i)=x1(i+1/2) -!
             x1b(imax), x2b(jmax), x3b(kmax), & !- cell center position for MKS-!
             x1ab(imax), x2ab(jmax), x3ab(kmax), & !- cell boundary position for MKS-!             
             dx1(imax), dx2(jmax), dx3(kmax), & !- x1(i+1)-x1(i) -!
             dx1b(2:imax-1), dx2b(2:jmax-1), dx3b(2:kmax-1) !-x1a(i+1)-x1a(i)-!

  integer :: nm0, is1, ie1, js1, je1, ks1, ke1
      
  if(irec .eq. 1) then
     call mclim(uri,x1,x2,x3,x1a,x2a,x3a,x1b,x2b,x3b,x1ab,x2ab,x3ab,&
                uriir,urijr,urikr,uriil,urijl,urikl,is1,ie1,js1,je1,ks1,ke1)
!     call mclim2(uri,x1,x2,x3,x1a,x2a,x3a,x1b,x2b,x3b,x1ab,x2ab,x3ab,&
!                uriir,urijr,urikr,uriil,urijl,urikl,is1,ie1,js1,je1,ks1,ke1)
!
  elseif(irec .eq. 2) then
    call minlim(uri,x1,x2,x3,x1a,x2a,x3a,x1b,x2b,x3b,x1ab,x2ab,x3ab,&
                uriir,urijr,urikr,uriil,urijl,urikl,is1,ie1,js1,je1,ks1,ke1)
!
  elseif(irec .eq. 3) then
    call muscl(uri,uriir,urijr,urikr,uriil,urijl,urikl,is1,ie1,js1,je1,ks1,ke1)
!
  elseif(irec .eq. 4) then
    call ceno(uri,x1,x2,x3,uriir,urijr,urikr,uriil,urijl,urikl, &
              is1,ie1,js1,je1,ks1,ke1)
!
  elseif(irec .eq. 5) then
    call ppm(uri,x1,x2,x3,dx1b,dx2b,dx3b, &
             uriir,urijr,urikr,uriil,urijl,urikl, &
             is1,ie1,js1,je1,ks1,ke1)
!
  elseif(irec .eq. 6) then
    call mp5(uri,uriir,urijr,urikr,uriil,urijl,urikl, &
             is1,ie1,js1,je1,ks1,ke1)
!
  elseif(irec .eq. 7) then
    call weno5(uri,uriir,urijr,urikr,uriil,urijl,urikl, &
               is1,ie1,js1,je1,ks1,ke1)
!
  elseif(irec .eq. 8) then
    call mpweno5(uri,uriir,urijr,urikr,uriil,urijl,urikl, &
                 is1,ie1,js1,je1,ks1,ke1)
!
  elseif(irec .eq. 9) then
    call weno5z(uri,uriir,urijr,urikr,uriil,urijl,urikl, &
                is1,ie1,js1,je1,ks1,ke1)
!
  elseif(irec .eq. 10) then
    call weno5m(uri,uriir,urijr,urikr,uriil,urijl,urikl, &
                is1,ie1,js1,je1,ks1,ke1)
!
  elseif(irec .eq. 11) then
!    call lim03(uri,x1,x2,x3,x1b,x2b,x3b,uriir,urijr,urikr,uriil,urijl,urikl, &
!               is1,ie1,js1,je1,ks1,ke1)

     call lim03a(uri,x1,x2,x3,x1a,x2a,x3a,x1b,x2b,x3b,x1ab,x2ab,x3ab,&
                uriir,urijr,urikr,uriil,urijl,urikl,is1,ie1,js1,je1,ks1,ke1)
!
  elseif(irec .eq. 12) then
    call mvllim(uri,x1,x2,x3,x1a,x2a,x3a,x1b,x2b,x3b,x1ab,x2ab,x3ab,&
                uriir,urijr,urikr,uriil,urijl,urikl,is1,ie1,js1,je1,ks1,ke1)
!    
  elseif(irec .eq. 13) then
    call koren(uri,x1,x2,x3,x1a,x2a,x3a,x1b,x2b,x3b,x1ab,x2ab,x3ab,&
                uriir,urijr,urikr,uriil,urijl,urikl,is1,ie1,js1,je1,ks1,ke1)
 endif
      
  call caluu3(uriir,uuir,gcovi,gconi,detgi,x1,nm0,is1,ie1,js1,je1,ks1,ke1)
  call caluu3(urijr,uujr,gcovj,gconj,detgj,x1,nm0,is1,ie1,js1,je1,ks1,ke1)
  call caluu3(urikr,uukr,gcovk,gconk,detgk,x1,nm0,is1,ie1,js1,je1,ks1,ke1)
  call caluu3(uriil,uuil,gcovi,gconi,detgi,x1,nm0,is1,ie1,js1,je1,ks1,ke1)
  call caluu3(urijl,uujl,gcovj,gconj,detgj,x1,nm0,is1,ie1,js1,je1,ks1,ke1)
  call caluu3(urikl,uukl,gcovk,gconk,detgk,x1,nm0,is1,ie1,js1,je1,ks1,ke1)

!     
  call calflx3i(uriir,uuir,wwir,gcovi,gconi,detgi, &
                  nm0,is1,ie1,js1,je1,ks1,ke1)
  call calflx3j(urijr,uujr,wwjr,gcovj,gconj,detgj, &
                  nm0,is1,ie1,js1,je1,ks1,ke1)
  call calflx3k(urikr,uukr,wwkr,gcovk,gconk,detgk, &
                  nm0,is1,ie1,js1,je1,ks1,ke1)
  call calflx3i(uriil,uuil,wwil,gcovi,gconi,detgi, &
                  nm0,is1,ie1,js1,je1,ks1,ke1)
  call calflx3j(urijl,uujl,wwjl,gcovj,gconj,detgj, &
                  nm0,is1,ie1,js1,je1,ks1,ke1)
  call calflx3k(urikl,uukl,wwkl,gcovk,gconk,detgk, &
                  nm0,is1,ie1,js1,je1,ks1,ke1)
!
  return
end subroutine rec2
!
!---------------------------------------------------------------------@
subroutine muscl(uri,uriir,urijr,urikr,uriil,urijl,urikl, &
                 is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Reconstruction step with MUSCL Method
!     Calculate cell-interface variables(uril, urir) 
!     from cell-center variables (uri)
!
  use pram, only : imax, jmax, kmax, nv
  implicit none
!
  integer :: i, j, k, n, is1, ie1, js1, je1, ks1, ke1
!
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: beta

!
!-----------------------------------------------------------------------

  beta=1.d0/3.d0

!=====================================================================@
  do k=ks1+1,ke1-2
    do j=js1+1,je1-2
      do i=is1+1,ie1-2
        do n=1,nv
!
           uriil(n,i,j,k)=uri(n,i,j,k) &
                         +0.25*((1.0+beta)*(uri(n,i+1,j,k)-uri(n,i,j,k)) &
                         +(1.0-beta)*(uri(n,i,j,k)-uri(n,i-1,j,k)))
           uriir(n,i,j,k)=uri(n,i+1,j,k) &
                         -0.25*((1.0+beta)*(uri(n,i+1,j,k)-uri(n,i,j,k)) &
                         +(1.0-beta)*(uri(n,i+2,j,k)-uri(n,i+1,j,k)))

           urijl(n,i,j,k)=uri(n,i,j,k) &
                         +0.25*((1.0+beta)*(uri(n,i,j+1,k)-uri(n,i,j,k)) &
                         +(1.0-beta)*(uri(n,i,j,k)-uri(n,i,j-1,k)))
           urijr(n,i,j,k)=uri(n,i,j+1,k) &
                         -0.25*((1.0+beta)*(uri(n,i,j+1,k)-uri(n,i,j,k)) &
                         +(1.0-beta)*(uri(n,i,j+2,k)-uri(n,i,j+1,k)))
    
           urikl(n,i,j,k)=uri(n,i,j,k) &
                         +0.25*((1.0+beta)*(uri(n,i,j,k+1)-uri(n,i,j,k)) &
                         +(1.0-beta)*(uri(n,i,j,k)-uri(n,i,j,k-1)))
           urikr(n,i,j,k)=uri(n,i,j,k+1) &
                         -0.25*((1.0+beta)*(uri(n,i,j,k+1)-uri(n,i,j,k)) &
                         +(1.0-beta)*(uri(n,i,j,k+2)-uri(n,i,j,k+1)))
!
        enddo
      enddo
    enddo
  enddo
!
  return
end subroutine muscl
!
!---------------------------------------------------------------------@
subroutine mclim(uri,x1,x2,x3,x1a,x2a,x3a,x1b,x2b,x3b,x1ab,x2ab,x3ab,&
                 uriir,urijr,urikr,uriil,urijl,urikl,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Reconstruction step with MC Limitter 
!     Calculate cell-interface variables(uril, urir, wwl, wwr) 
!     from cell-center variables (uri)
!
  use pram, only : imax, jmax, kmax, nv, metric
  implicit none
!
  integer :: i, j, k, n, is1, ie1, js1, je1, ks1, ke1
!
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: x1(imax), x2(jmax), x3(kmax), & !- cell center position
             x1a(imax), x2a(jmax), x3a(kmax), & !- cell boundary position
                                               !- x1a(i)=x1(i+1/2)
             x1b(imax), x2b(jmax), x3b(kmax), & !- cell center position for MKS-!
             x1ab(imax), x2ab(jmax), x3ab(kmax) !- cell boundary position for MKS-! 
!
  real(8) :: sp2i, sm2i, spm2i, dk2i, tmp2i, &
             sp2j, sm2j, spm2j, dk2j, tmp2j, &
             sp2k, sm2k, spm2k, dk2k, tmp2k, &
             x1p, x1c, x1m, x1ap, x1am, &
             x2p, x2c, x2m, x2ap, x2am, &
             x3p, x3c, x3m, x3ap, x3am
!
  real(8), parameter :: small=1.d-12 
!
!-----------------------------------------------------------------------

  do k=ks1+1,ke1-1
    do j=js1+1,je1-1
      do i=is1+1,ie1-1
        do n=1,nv
!
!!! i-th direction
!           
          if(metric .eq. 403) then
            x1p=x1b(i+1)
            x1c=x1b(i)
            x1m=x1b(i-1)
            x1ap=x1ab(i)
            x1am=x1ab(i-1)
          else
            x1p=x1(i+1)
            x1c=x1(i)
            x1m=x1(i-1)
            x1ap=x1a(i)
            x1am=x1a(i-1)
          endif   
!
          sp2i=(uri(n,i+1,j,k)-uri(n,i,j,k))*(1./(x1p-x1c+small))
          sm2i=(uri(n,i,j,k)-uri(n,i-1,j,k))*(1./(x1c-x1m+small))

          spm2i=sp2i*sm2i
!- check sign
          if(sp2i .lt. 0.d0) then
            tmp2i=-1.d0
          else
            tmp2i=1.d0
          endif
!- limiter
          dk2i=tmp2i*max(0.d0, min(2.d0*abs(sp2i), tmp2i*2.d0*sm2i, &
                         tmp2i*0.5d0*(sp2i+sm2i)))

!          uriil(n,i,j,k)=uri(n,i,j,k)+dk2i*(0.5*(x1p-x1c))
!          uriir(n,i-1,j,k)=uri(n,i,j,k)-dk2i*(0.5*(x1c-x1p))
!
          uriil(n,i,j,k)=uri(n,i,j,k)+dk2i*(x1ap-x1c)
          uriir(n,i-1,j,k)=uri(n,i,j,k)-dk2i*(x1c-x1am)
!
!!! j-th direction
!
          if(metric .eq. 403) then
            x2p=x2b(j+1)
            x2c=x2b(j)
            x2m=x2b(j-1)
            x2ap=x2ab(j)
            x2am=x2ab(j-1)
          else
            x2p=x2(j+1)
            x2c=x2(j)
            x2m=x2(j-1)
            x2ap=x2a(j)
            x2am=x2a(j-1)
          endif  

          sp2j=(uri(n,i,j+1,k)-uri(n,i,j,k))*(1./(x2p-x2c+small))
          sm2j=(uri(n,i,j,k)-uri(n,i,j-1,k))*(1./(x2c-x2m+small))

          spm2j=sp2j*sm2j
!- check sign
          if(sp2j .lt. 0.d0) then
            tmp2j=-1.d0
          else
            tmp2j=1.d0
          endif
!- limiter
          dk2j=tmp2j*max(0.d0, min(2.d0*abs(sp2j), tmp2j*2.d0*sm2j, &
                         tmp2j*0.5d0*(sp2j+sm2j)))
            
!          endif

!          urijl(n,i,j,k)=uri(n,i,j,k)+dk2j*(0.5*(x2p-x2c))
!          urijr(n,i,j-1,k)=uri(n,i,j,k)-dk2j*(0.5*(x2c-x2m))
!
          urijl(n,i,j,k)=uri(n,i,j,k)+dk2j*(x2ap-x2c)
          urijr(n,i,j-1,k)=uri(n,i,j,k)-dk2j*(x2c-x2am)
!
!!! k-th direction          
!
          if(metric .eq. 403) then
            x3p=x3b(k+1)
            x3c=x3b(k)
            x3m=x3b(k-1)
            x3ap=x3ab(k)
            x3am=x3ab(k-1)
          else
            x3p=x3(k+1)
            x3c=x3(k)
            x3m=x3(k-1)
            x3ap=x3a(k)
            x3am=x3a(k-1)
          endif
!         
          sp2k=(uri(n,i,j,k+1)-uri(n,i,j,k))*(1./(x3p-x3c+small))
          sm2k=(uri(n,i,j,k)-uri(n,i,j,k-1))*(1./(x3c-x3m+small))

          spm2k=sp2k*sm2k
!-check sign
          if(sp2k .lt. 0.d0) then
            tmp2k=-1.d0
          else
            tmp2k=1.d0
          endif
!- limiter
          dk2k=tmp2k*max(0.d0, min(2.d0*abs(sp2k), tmp2k*2.d0*sm2k, &
                         tmp2k*0.5d0*(sp2k+sm2k))) 
          
!          urikl(n,i,j,k)=uri(n,i,j,k)+dk2k*(0.5*(x3p-x3c))
!          urikr(n,i,j,k-1)=uri(n,i,j,k)-dk2k*(0.5*(x3c-x3m))
!
          urikl(n,i,j,k)=uri(n,i,j,k)+dk2k*(x3ap-x3c)
          urikr(n,i,j,k-1)=uri(n,i,j,k)-dk2k*(x3c-x3am)
          
        enddo
      enddo
    enddo
  enddo

  return
end subroutine mclim
!
!---------------------------------------------------------------------@
subroutine mclim2(uri,x1,x2,x3,x1a,x2a,x3a,x1b,x2b,x3b,x1ab,x2ab,x3ab,&
                 uriir,urijr,urikr,uriil,urijl,urikl,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Reconstruction step with MC Limitter 
!     Calculate cell-interface variables(uril, urir, wwl, wwr) 
!     from cell-center variables (uri)
!
  use pram, only : imax, jmax, kmax, nv, metric
  implicit none
!
  integer :: i, j, k, n, is1, ie1, js1, je1, ks1, ke1
!
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: x1(imax), x2(jmax), x3(kmax), & !- cell center position
             x1a(imax), x2a(jmax), x3a(kmax), & !- cell boundary position
                                               !- x1a(i)=x1(i+1/2)
             x1b(imax), x2b(jmax), x3b(kmax), & !- cell center position for MKS-!
             x1ab(imax), x2ab(jmax), x3ab(kmax) !- cell boundary position for MKS-! 
!
  real(8) :: sp2i, sm2i, spm2i, dk2i, tmp2i, &
             sp2j, sm2j, spm2j, dk2j, tmp2j, &
             sp2k, sm2k, spm2k, dk2k, tmp2k, &
             cfi, cbi, cfj, cbj, cfk, cbk, &
             x1p, x1c, x1m, x1ap, x1am, &
             x2p, x2c, x2m, x2ap, x2am, &
             x3p, x3c, x3m, x3ap, x3am, &
             tmp1ia, tmp1ib, tmp1ic, tmp1id, tmp1i, &
             tmp1ja, tmp1jb, tmp1jc, tmp1jd, tmp1j, &
             tmp1ka, tmp1kb, tmp1kc, tmp1kd, tmp1k
!
  real(8), parameter :: small=1.d-12 
!
!-----------------------------------------------------------------------

  do k=ks1+1,ke1-1
    do j=js1+1,je1-1
      do i=is1+1,ie1-1
        do n=1,nv
!
!!! i-th direction
!           
          if(metric .eq. 403) then
            x1p=x1b(i+1)
            x1c=x1b(i)
            x1m=x1b(i-1)
            x1ap=x1ab(i)
            x1am=x1ab(i-1)
          else
            x1p=x1(i+1)
            x1c=x1(i)
            x1m=x1(i-1)
            x1ap=x1a(i)
            x1am=x1a(i-1)
          endif
!         
          tmp1ia=x1p-x1c+small
          tmp1ib=x1c-x1m+small
          tmp1ic=x1ap-x1c+small
          tmp1id=x1c-x1am+small 
!
!- S+, S-          
          sp2i=(uri(n,i+1,j,k)-uri(n,i,j,k))*(1.d0/tmp1ia)   
          sm2i=(uri(n,i,j,k)-uri(n,i-1,j,k))*(1.d0/tmp1ib)  

          spm2i=sp2i*sm2i          
!- CF, CB  
          cfi=tmp1ia*(1.d0/tmp1ic)
          cbi=tmp1ib*(1.d0/tmp1id)
!
!- limiter
          if(sp2i .lt. 0.d0) then
            tmp2i=-1.d0
          else
            tmp2i=1.d0
          endif   
          
          dk2i=tmp2i*max(0.d0, min(abs(cfi*sp2i), tmp2i*cbi*sm2i, &
                         tmp2i*0.5d0*(sp2i+sm2i)))

!           uriil(n,i,j,k)=uri(n,i,j,k)+dk2i*(0.5*(x1p-x1c))
!           uriir(n,i-1,j,k)=uri(n,i,j,k)-dk2i*(0.5*(x1c-x1p))
!
          uriil(n,i,j,k)=uri(n,i,j,k)+dk2i*(x1ap-x1c)
          uriir(n,i-1,j,k)=uri(n,i,j,k)-dk2i*(x1c-x1am)
!
!!! j-th direction
!
          if(metric .eq. 403) then
            x2p=x2b(j+1)
            x2c=x2b(j)
            x2m=x2b(j-1)
            x2ap=x2ab(j)
            x2am=x2ab(j-1)
          else
            x2p=x2(j+1)
            x2c=x2(j)
            x2m=x2(j-1)
            x2ap=x2a(j)
            x2am=x2a(j-1)
          endif  

          tmp1ja=x2p-x2c+small
          tmp1jb=x2c-x2m+small
          tmp1jc=x2ap-x2c+small
          tmp1jd=x2c-x2am+small          
!
!          sp2j=(uri(n,i,j+1,k)-uri(n,i,j,k))/dx2(j+1)
!          sm2j=(uri(n,i,j,k)-uri(n,i,j-1,k))/dx2(j)
!  
          sp2j=(uri(n,i,j+1,k)-uri(n,i,j,k))*(1.d0/tmp1ja) 
          sm2j=(uri(n,i,j,k)-uri(n,i,j-1,k))*(1.d0/tmp1jb)

          spm2j=sp2j*sm2j
!- CF, CB  
          cfj=tmp1ja*(1.d0/tmp1jc)  
          cbj=tmp1jb*(1.d0/tmp1jd)
!
!- limiter
          if(sp2j .lt. 0.d0) then
            tmp2j=-1.d0
          else
            tmp2j=1.d0
          endif   
          
          dk2j=tmp2j*max(0.d0, min(abs(cfj*sp2j), tmp2j*cbj*sm2j, &
                         tmp2j*0.5d0*(sp2j+sm2j)))

!          urijl(n,i,j,k)=uri(n,i,j,k)+dk2j*(0.5*(x2p-x2c))
!          urijr(n,i,j-1,k)=uri(n,i,j,k)-dk2j*(0.5*(x2c-x2m))
!
          urijl(n,i,j,k)=uri(n,i,j,k)+dk2j*(x2ap-x2c)
          urijr(n,i,j-1,k)=uri(n,i,j,k)-dk2j*(x2c-x2am)
!
!!! k-th direction          
!
          if(metric .eq. 403) then
            x3p=x3b(k+1)
            x3c=x3b(k)
            x3m=x3b(k-1)
            x3ap=x3ab(k)
            x3am=x3ab(k-1)
          else
            x3p=x3(k+1)
            x3c=x3(k)
            x3m=x3(k-1)
            x3ap=x3a(k)
            x3am=x3a(k-1)
          endif
!
          tmp1ka=x3p-x3c+small
          tmp1kb=x3c-x3m+small 
          tmp1kc=x3ap-x3c+small
          tmp1kd=x3c-x3am+small
!
!          sp2k=(uri(n,i,j,k+1)-uri(n,i,j,k))/dx3(k+1)
!          sm2k=(uri(n,i,j,k)-uri(n,i,j,k-1))/dx3(k)
!                   
          sp2k=(uri(n,i,j,k+1)-uri(n,i,j,k))*(1.d0/tmp1ka)
          sm2k=(uri(n,i,j,k)-uri(n,i,j,k-1))*(1.d0/tmp1kb)   

          spm2k=sp2k*sm2k
!- CF, CB   
          cfk=tmp1ka*(1.d0/tmp1kc)  
          cbk=tmp1kb*(1.d0/tmp1kd)

!- limiter
          if(sp2k .lt. 0.d0) then
            tmp2k=-1.d0
          else
            tmp2k=1.d0
          endif   
          
          dk2k=tmp2k*max(0.d0, min(abs(cfk*sp2k), tmp2k*cbk*sm2k, &
                         tmp2k*0.5d0*(sp2k+sm2k)))      
          
!          urikl(n,i,j,k)=uri(n,i,j,k)+dk2k*(0.5*(x3p-x3c))
!          urikr(n,i,j,k-1)=uri(n,i,j,k)-dk2k*(0.5*(x3c-x3m))
!
          urikl(n,i,j,k)=uri(n,i,j,k)+dk2k*(x3ap-x3c)
          urikr(n,i,j,k-1)=uri(n,i,j,k)-dk2k*(x3c-x3am)
          
        enddo
      enddo
    enddo
  enddo

  return
end subroutine mclim2
!
!---------------------------------------------------------------------@
subroutine minlim(uri,x1,x2,x3,x1a,x2a,x3a,x1b,x2b,x3b,x1ab,x2ab,x3ab, &
                  uriir,urijr,urikr,uriil,urijl,urikl,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Reconstruction step with Minmod Limitter 
!     Calculate cell-interface variables(uriil, uriir) 
!     from cell-center variables (uri)
!
  use pram, only : imax, jmax, kmax, nv, metric
  implicit none
!
  integer :: i, j, k, n, is1, ie1, js1, je1, ks1, ke1
!
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: x1(imax), x2(jmax), x3(kmax), & !- cell center position
             x1a(imax), x2a(jmax), x3a(kmax), & !- cell boundary position
                                                !- x1a(i)=x1(i+1/2)
             x1b(imax), x2b(jmax), x3b(kmax), & !- cell center position for MKS-!
             x1ab(imax), x2ab(jmax), x3ab(kmax) !- cell boundary position for MKS-!  
!
  real(8) :: sp2i, sm2i, spm2i, dk2i, tmp2i, &
             sp2j, sm2j, spm2j, dk2j, tmp2j, &
             sp2k, sm2k, spm2k, dk2k, tmp2k, &
             x1p, x1c, x1m, x1ap, x1am, &
             x2p, x2c, x2m, x2ap, x2am, &
             x3p, x3c, x3m, x3ap, x3am
!
  real(8), parameter :: small=1.d-12 
!-----------------------------------------------------------------------

  do k=ks1+1,ke1-1
    do j=js1+1,je1-1
      do i=is1+1,ie1-1
        do n=1,nv
!
!!! i-th direction
!
          if(metric .eq. 403) then
            x1p=x1b(i+1)
            x1c=x1b(i)
            x1m=x1b(i-1)
            x1ap=x1ab(i)
            x1am=x1ab(i-1)
          else
            x1p=x1(i+1)
            x1c=x1(i)
            x1m=x1(i-1)
            x1ap=x1a(i)
            x1am=x1a(i-1)
          endif   
!           
!          sp2i=(uri(n,i+1,j,k)-uri(n,i,j,k))/dx1(i+1)
!          sm2i=(uri(n,i,j,k)-uri(n,i-1,j,k))/dx1(i)
!- S+, S-
          sp2i=(uri(n,i+1,j,k)-uri(n,i,j,k))*(1./(x1p-x1c+small))
          sm2i=(uri(n,i,j,k)-uri(n,i-1,j,k))*(1./(x1c-x1m+small))
!- check sign          
          spm2i=sp2i*sm2i

          if(sp2i .lt. 0.d0) then
            tmp2i=-1.d0
          else
            tmp2i=1.d0
          endif
!- limiter
          dk2i=tmp2i*max(0.d0, min(abs(sp2i), tmp2i*sm2i))

!          uriil(n,i,j,k)=uri(n,i,j,k)+dk2i*(0.5*(x1p-x1c))
!          uriir(n,i-1,j,k)=uri(n,i,j,k)-dk2i*(0.5*(x1c-x1m))
!          
          uriil(n,i,j,k)=uri(n,i,j,k)+dk2i*(x1ap-x1c)
          uriir(n,i-1,j,k)=uri(n,i,j,k)-dk2i*(x1c-x1am)
!
!          uriil(n,i,j,k)=uri(n,i,j,k)+dk2i*(0.5*dx1(i))
!          uriir(n,i-1,j,k)=uri(n,i,j,k)-dk2i*(0.5*dx1(i))
!
!!! j-th direction
!          
          if(metric .eq. 403) then
            x2p=x2b(j+1)
            x2c=x2b(j)
            x2m=x2b(j-1)
            x2ap=x2ab(j)
            x2am=x2ab(j-1)
          else
            x2p=x2(j+1)
            x2c=x2(j)
            x2m=x2(j-1)
            x2ap=x2a(j)
            x2am=x2a(j-1)
          endif
!
!          sp2j=(uri(n,i,j+1,k)-uri(n,i,j,k))/dx2(j+1)
!          sm2j=(uri(n,i,j,k)-uri(n,i,j-1,k))/dx2(j)
!
          sp2j=(uri(n,i,j+1,k)-uri(n,i,j,k))*(1./(x2p-x2c+small))
          sm2j=(uri(n,i,j,k)-uri(n,i,j-1,k))*(1./(x2c-x2m+small))

          spm2j=sp2j*sm2j
!- check sign
          if(sp2j .lt. 0.d0) then
            tmp2j=-1.d0
          else
            tmp2j=1.d0
          endif
!- limiter
          dk2j=tmp2j*max(0.d0, min(abs(sp2j), tmp2j*sm2j))
          
!          urijl(n,i,j,k)=uri(n,i,j,k)+dk2j*(0.5*(x2p-x2c))
!          urijr(n,i,j-1,k)=uri(n,i,j,k)-dk2j*(0.5*(x2c-x2m))
!
          urijl(n,i,j,k)=uri(n,i,j,k)+dk2j*(x2ap-x2c)
          urijr(n,i,j-1,k)=uri(n,i,j,k)-dk2j*(x2c-x2am)
!     
!          urijl(n,i,j,k)=uri(n,i,j,k)+dk2j*(0.5*dx2(j))
!          urijr(n,i,j-1,k)=uri(n,i,j,k)-dk2j*(0.5*dx2(j))
!
!!! k-th direction
!          
          if(metric .eq. 403) then
            x3p=x3b(k+1)
            x3c=x3b(k)
            x3m=x3b(k-1)
            x3ap=x3ab(k)
            x3am=x3ab(k-1)
          else
            x3p=x3(k+1)
            x3c=x3(k)
            x3m=x3(k-1)
            x3ap=x3a(k)
            x3am=x3a(k-1)
          endif
!
!          sp2k=(uri(n,i,j,k+1)-uri(n,i,j,k))/dx3(k+1)
!          sm2k=(uri(n,i,j,k)-uri(n,i,j,k-1))/dx3(k)
!          
          sp2k=(uri(n,i,j,k+1)-uri(n,i,j,k))*(1./(x3p-x3c+small))
          sm2k=(uri(n,i,j,k)-uri(n,i,j,k-1))*(1./(x3c-x3m+small))

          spm2k=sp2k*sm2k
!- check sign          
          if(sp2k .lt. 0.d0) then
            tmp2k=-1.d0
          else
            tmp2k=1.d0
          endif
!- limiter           
          dk2k=tmp2k*max(0.d0, min(abs(sp2k), tmp2k*sm2k))

!          urikl(n,i,j,k)=uri(n,i,j,k)+dk2k*(0.5*(x3p-x3c))
!          urikr(n,i,j,k-1)=uri(n,i,j,k)-dk2k*(0.5*(x3c-x3m))
!          
          urikl(n,i,j,k)=uri(n,i,j,k)+dk2k*(x3ap-x3c)
          urikr(n,i,j,k-1)=uri(n,i,j,k)-dk2k*(x3c-x3am)
          
!          urikl(n,i,j,k)=uri(n,i,j,k)+dk2k*(0.5*dx3(k))
!          urikr(n,i,j,k-1)=uri(n,i,j,k)-dk2k*(0.5*dx3(k))
!          
        enddo
      enddo
    enddo
  enddo
!
  return
end subroutine minlim
!
!---------------------------------------------------------------------@
subroutine ceno(uri,x1,x2,x3,uriir,urijr,urikr,uriil,urijl,urikl, &
                is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Reconstruction step with Convex ENO method 
!     Calculate cell-interface variables(uril, urir) 
!     from cell-center variables (uri)
!
  use pram, only : imax, jmax, kmax, nv
  implicit none
!
  integer :: i, j, k, n, is1, ie1, js1, je1, ks1, ke1
  integer :: merr
!
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: uriir1(:,:,:,:), uriil1(:,:,:,:), &
           urijr1(:,:,:,:), urijl1(:,:,:,:), &
           urikr1(:,:,:,:), urikl1(:,:,:,:)
!
  real(8), allocatable :: uriir2m(:,:,:,:), uriil2m(:,:,:,:), &
           urijr2m(:,:,:,:), urijl2m(:,:,:,:), &
           urikr2m(:,:,:,:), urikl2m(:,:,:,:)
!              
  real(8), allocatable :: uriir20(:,:,:,:), uriil20(:,:,:,:), &
           urijr20(:,:,:,:), urijl20(:,:,:,:), &
           urikr20(:,:,:,:), urikl20(:,:,:,:)
!              
  real(8), allocatable :: uriir2p(:,:,:,:), uriil2p(:,:,:,:), &
           urijr2p(:,:,:,:), urijl2p(:,:,:,:), &
           urikr2p(:,:,:,:), urikl2p(:,:,:,:)
!
  real(8), allocatable :: dirm(:,:,:,:), dilm(:,:,:,:), &
           djrm(:,:,:,:), djlm(:,:,:,:), dkrm(:,:,:,:), dklm(:,:,:,:)
!      
  real(8), allocatable :: dir0(:,:,:,:), dil0(:,:,:,:), &
           djr0(:,:,:,:), djl0(:,:,:,:),  dkr0(:,:,:,:), dkl0(:,:,:,:)
               
  real(8), allocatable :: dirp(:,:,:,:), dilp(:,:,:,:), &
           djrp(:,:,:,:), djlp(:,:,:,:), dkrp(:,:,:,:), dklp(:,:,:,:)
      
  real(8) :: x1(imax), x2(jmax), x3(kmax)
!
  real(8) :: alp1, alp2, alp3
  real(8) :: sp1i, sm1i, spm1i, dk1i, tmp1i, &
             sp1j, sm1j, spm1j, dk1j, tmp1j, &
             sp1k, sm1k, spm1k, dk1k, tmp1k, &
             d2im, d2i0, d2ip, d3im, d3i0, d3ip, tmp2im, tmp2i0, tmp2ip, &
             tmp3im, tmp3i0, tmp3ip, &
             d2jm, d2j0, d2jp, d3jm, d3j0, d3jp, tmp2jm, tmp2j0, tmp2jp, &
             tmp3jm, tmp3j0, tmp3jp, &
             d2km, d2k0, d2kp, d3km, d3k0, d3kp, tmp2km, tmp2k0, tmp2kp, &
             tmp3km, tmp3k0, tmp3kp
!
  allocate( uriir1(nv,is1:ie1,js1:je1,ks1:ke1), &
            uriil1(nv,is1:ie1,js1:je1,ks1:ke1), &
            urijr1(nv,is1:ie1,js1:je1,ks1:ke1), &
            urijl1(nv,is1:ie1,js1:je1,ks1:ke1), &
            urikr1(nv,is1:ie1,js1:je1,ks1:ke1), &
            urikl1(nv,is1:ie1,js1:je1,ks1:ke1), &
            uriir2m(nv,is1:ie1,js1:je1,ks1:ke1), & 
            uriil2m(nv,is1:ie1,js1:je1,ks1:ke1), &
            urijr2m(nv,is1:ie1,js1:je1,ks1:ke1), & 
            urijl2m(nv,is1:ie1,js1:je1,ks1:ke1), &
            urikr2m(nv,is1:ie1,js1:je1,ks1:ke1), &
            urikl2m(nv,is1:ie1,js1:je1,ks1:ke1), &
            uriir20(nv,is1:ie1,js1:je1,ks1:ke1), &
            uriil20(nv,is1:ie1,js1:je1,ks1:ke1), &
            urijr20(nv,is1:ie1,js1:je1,ks1:ke1), &
            urijl20(nv,is1:ie1,js1:je1,ks1:ke1), &
            urikr20(nv,is1:ie1,js1:je1,ks1:ke1), &
            urikl20(nv,is1:ie1,js1:je1,ks1:ke1), &   
            uriir2p(nv,is1:ie1,js1:je1,ks1:ke1), &
            uriil2p(nv,is1:ie1,js1:je1,ks1:ke1), &
            urijr2p(nv,is1:ie1,js1:je1,ks1:ke1), &
            urijl2p(nv,is1:ie1,js1:je1,ks1:ke1), &
            urikr2p(nv,is1:ie1,js1:je1,ks1:ke1), &
            urikl2p(nv,is1:ie1,js1:je1,ks1:ke1), &
            dirm(nv,is1:ie1,js1:je1,ks1:ke1), &
            dilm(nv,is1:ie1,js1:je1,ks1:ke1), &
            djrm(nv,is1:ie1,js1:je1,ks1:ke1), &
            djlm(nv,is1:ie1,js1:je1,ks1:ke1), &
            dkrm(nv,is1:ie1,js1:je1,ks1:ke1), &
            dklm(nv,is1:ie1,js1:je1,ks1:ke1), &
            dir0(nv,is1:ie1,js1:je1,ks1:ke1), &
            dil0(nv,is1:ie1,js1:je1,ks1:ke1), &
            djr0(nv,is1:ie1,js1:je1,ks1:ke1), &
            djl0(nv,is1:ie1,js1:je1,ks1:ke1), &
            dkr0(nv,is1:ie1,js1:je1,ks1:ke1), &
            dkl0(nv,is1:ie1,js1:je1,ks1:ke1), &
            dirp(nv,is1:ie1,js1:je1,ks1:ke1), &
            dilp(nv,is1:ie1,js1:je1,ks1:ke1), &
            djrp(nv,is1:ie1,js1:je1,ks1:ke1), &
            djlp(nv,is1:ie1,js1:je1,ks1:ke1), &
            dkrp(nv,is1:ie1,js1:je1,ks1:ke1), &
            dklp(nv,is1:ie1,js1:je1,ks1:ke1), stat=merr)
!
!-----------------------------------------------------------------------
!     Parameter

  alp1=1.d0
  alp2=0.7d0
  alp3=1.d0

!-----------------------------------------------------------------------

  do k=ks1+1,ke1-1
    do j=js1+1,je1-1
      do i=is1+1,ie1-1
        do n=1,nv
!
!  linear polynomial
!
          sp1i=(uri(n,i+1,j,k)-uri(n,i,j,k))/(x1(i+1)-x1(i))
          sm1i=(uri(n,i,j,k)-uri(n,i-1,j,k))/(x1(i)-x1(i-1))
          
          spm1i=sp1i*sm1i
          if(spm1i .le. 0.d0) then
            dk1i=0.d0
          else
            if(sp1i .lt. 0.d0) then
              tmp1i=-1.d0
            else
              tmp1i=1.d0
            endif
! Minmod
            dk1i=tmp1i*min(abs(sp1i), abs(sm1i))
! MC
!           dk1i=tmp1i*min(2.0*abs(sp1i), 2.0*abs(sm1i), abs(sp1i+sm1i)/2.0)
!
          endif

          uriil1(n,i,j,k)=uri(n,i,j,k)+dk1i*(0.5*(x1(i+1)-x1(i)))
          uriir1(n,i-1,j,k)=uri(n,i,j,k)-dk1i*(0.5*(x1(i)-x1(i-1)))
!          
          sp1j=(uri(n,i,j+1,k)-uri(n,i,j,k))/(x2(j+1)-x2(j))
          sm1j=(uri(n,i,j,k)-uri(n,i,j-1,k))/(x2(j)-x2(j-1))

          spm1j=sp1j*sm1j
          if(spm1j .le. 0.d0) then
            dk1j=0.d0
          else
            if(sp1j .lt. 0.d0) then
              tmp1j=-1.0
            else
              tmp1j=1.0
            endif
! Minmod
            dk1j=tmp1j*min(abs(sp1j), abs(sm1j))
! MC
!           dk1j=tmp1j*min(2.0*abs(sp1j), 2.0*abs(sm1j), abs(sp1j+sm1j)/2.0)
!
          endif
          
          urijl1(n,i,j,k)=uri(n,i,j,k)+dk1j*(0.5*(x2(j+1)-x2(j)))
          urijr1(n,i,j-1,k)=uri(n,i,j,k)-dk1j*(0.5*(x2(j)-x2(j-1)))
!    
          sp1k=(uri(n,i,j,k+1)-uri(n,i,j,k))/(x3(k+1)-x3(k))
          sm1k=(uri(n,i,j,k)-uri(n,i,j,k-1))/(x3(k)-x3(k-1))

          spm1k=sp1k*sm1k
          if(spm1k .le. 0.d0) then
            dk1k=0.d0
          else
            if(sp1k .lt. 0.d0) then
              tmp1k=-1.d0
            else
              tmp1k=1.d0
            endif
! Minmod
            dk1k=tmp1k*min(abs(sp1k), abs(sm1k))
! MC
!           dk1k=tmp1k*min(2.0*abs(sp1k), 2.0*abs(sm1k), abs(sp1k+sm1k)/2.0)
!
          endif
          
          urikl1(n,i,j,k)=uri(n,i,j,k)+dk1k*(0.5*(x3(k+1)-x3(k)))
          urikr1(n,i,j,k-1)=uri(n,i,j,k)-dk1k*(0.5*(x3(k)-x3(k-1)))
          
        enddo
      enddo
    enddo
  enddo
!
  do k=ks1+2,ke1-2
    do j=js1+2,je1-2
      do i=is1+2,ie1-2
        do n=1,nv
 
          d2im=0.5*(uri(n,i,j,k)-uri(n,i-2,j,k))
          d2i0=0.5*(uri(n,i+1,j,k)-uri(n,i-1,j,k))
          d2ip=0.5*(uri(n,i+2,j,k)-uri(n,i,j,k))
         
          d3im=uri(n,i,j,k)-2.0*uri(n,i-1,j,k)+uri(n,i-2,j,k)
          d3i0=uri(n,i+1,j,k)-2.0*uri(n,i,j,k)+uri(n,i-1,j,k)
          d3ip=uri(n,i+2,j,k)-2.0*uri(n,i+1,j,k)+uri(n,i,j,k)
         
          tmp2im=(0.5*(x1(i+1)+x1(i))-x1(i-1))/(x1(i+1)-x1(i))
          tmp2i0=(0.5*(x1(i+1)+x1(i))-x1(i))/(x1(i+1)-x1(i))
          tmp2ip=(0.5*(x1(i+1)+x1(i))-x1(i+1))/(x1(i+1)-x1(i))
         
          tmp3im=(0.5*(x1(i)+x1(i-1))-x1(i-1))/(x1(i)-x1(i-1))
          tmp3i0=(0.5*(x1(i)+x1(i-1))-x1(i))/(x1(i)-x1(i-1))
          tmp3ip=(0.5*(x1(i)+x1(i-1))-x1(i+1))/(x1(i)-x1(i-1))
         
          uriil2m(n,i,j,k)=uri(n,i,j,k)+d2im*tmp2im+0.5*d3im*tmp2im**2
          uriil20(n,i,j,k)=uri(n,i,j,k)+d2i0*tmp2i0+0.5*d3i0*tmp2i0**2
          uriil2p(n,i,j,k)=uri(n,i,j,k)+d2ip*tmp2ip+0.5*d3ip*tmp2ip**2
         
          uriir2m(n,i-1,j,k)=uri(n,i,j,k)+d2im*tmp3im+0.5*d3im*tmp3im**2
          uriir20(n,i-1,j,k)=uri(n,i,j,k)+d2i0*tmp3i0+0.5*d3i0*tmp3i0**2
          uriir2p(n,i-1,j,k)=uri(n,i,j,k)+d2ip*tmp3ip+0.5*d3ip*tmp3ip**2
         
          dilm(n,i,j,k)=alp1*(uriil2m(n,i,j,k)-uriil1(n,i,j,k))
          dil0(n,i,j,k)=alp2*(uriil20(n,i,j,k)-uriil1(n,i,j,k))
          dilp(n,i,j,k)=alp3*(uriil2p(n,i,j,k)-uriil1(n,i,j,k))
         
          dirm(n,i-1,j,k)=alp1*(uriir2m(n,i-1,j,k)-uriir1(n,i-1,j,k))
          dir0(n,i-1,j,k)=alp2*(uriir20(n,i-1,j,k)-uriir1(n,i-1,j,k))
          dirp(n,i-1,j,k)=alp3*(uriir2p(n,i-1,j,k)-uriir1(n,i-1,j,k))
!
          d2jm=0.5*(uri(n,i,j,k)-uri(n,i,j-2,k))
          d2j0=0.5*(uri(n,i,j+1,k)-uri(n,i,j-1,k))
          d2jp=0.5*(uri(n,i,j+2,k)-uri(n,i,j,k))

          d3jm=uri(n,i,j,k)-2.0*uri(n,i,j-1,k)+uri(n,i,j-2,k)
          d3j0=uri(n,i,j+1,k)-2.0*uri(n,i,j,k)+uri(n,i,j-1,k)
          d3jp=uri(n,i,j+2,k)-2.0*uri(n,i,j+1,k)+uri(n,i,j,k)
 
          tmp2jm=(0.5*(x2(j+1)+x2(j))-x2(j-1))/(x2(j+1)-x2(j))
          tmp2j0=(0.5*(x2(j+1)+x2(j))-x2(j))/(x2(j+1)-x2(j))
          tmp2jp=(0.5*(x2(j+1)+x2(j))-x2(j+1))/(x2(j+1)-x2(j))
         
          tmp3jm=(0.5*(x2(j)+x2(j-1))-x2(j-1))/(x2(j)-x2(j-1))
          tmp3j0=(0.5*(x2(j)+x2(j-1))-x2(j))/(x2(j)-x2(j-1))
          tmp3jp=(0.5*(x2(j)+x2(j-1))-x2(j+1))/(x2(j)-x2(j-1))
         
          urijl2m(n,i,j,k)=uri(n,i,j,k)+d2jm*tmp2jm+0.5*d3jm*tmp2jm**2
          urijl20(n,i,j,k)=uri(n,i,j,k)+d2j0*tmp2j0+0.5*d3j0*tmp2j0**2
          urijl2p(n,i,j,k)=uri(n,i,j,k)+d2jp*tmp2jp+0.5*d3jp*tmp2jp**2
         
          urijr2m(n,i,j-1,k)=uri(n,i,j,k)+d2jm*tmp3jm+0.5*d3jm*tmp3jm**2
          urijr20(n,i,j-1,k)=uri(n,i,j,k)+d2j0*tmp3j0+0.5*d3j0*tmp3j0**2
          urijr2p(n,i,j-1,k)=uri(n,i,j,k)+d2jp*tmp3jp+0.5*d3jp*tmp3jp**2
         
          djlm(n,i,j,k)=alp1*(urijl2m(n,i,j,k)-urijl1(n,i,j,k))
          djl0(n,i,j,k)=alp2*(urijl20(n,i,j,k)-urijl1(n,i,j,k))
          djlp(n,i,j,k)=alp3*(urijl2p(n,i,j,k)-urijl1(n,i,j,k))
         
          djrm(n,i,j-1,k)=alp1*(urijr2m(n,i,j-1,k)-urijr1(n,i,j-1,k))
          djr0(n,i,j-1,k)=alp2*(urijr20(n,i,j-1,k)-urijr1(n,i,j-1,k))
          djrp(n,i,j-1,k)=alp3*(urijr2p(n,i,j-1,k)-urijr1(n,i,j-1,k))
!
          d2km=0.5*(uri(n,i,j,k)-uri(n,i,j,k-2))
          d2k0=0.5*(uri(n,i,j,k+1)-uri(n,i,j,k-1))
          d2kp=0.5*(uri(n,i,j,k+2)-uri(n,i,j,k))
         
          d3km=uri(n,i,j,k)-2.0*uri(n,i,j,k-1)+uri(n,i,j,k-2)
          d3k0=uri(n,i,j,k+1)-2.0*uri(n,i,j,k)+uri(n,i,j,k-1)
          d3kp=uri(n,i,j,k+2)-2.0*uri(n,i,j,k+1)+uri(n,i,j,k)
         
          tmp2km=(0.5*(x3(k+1)+x3(k))-x3(k-1))/(x3(k+1)-x3(k))
          tmp2k0=(0.5*(x3(k+1)+x3(k))-x3(k))/(x3(k+1)-x3(k))
          tmp2kp=(0.5*(x3(k+1)+x3(k))-x3(k+1))/(x3(k+1)-x3(k))
         
          tmp3km=(0.5*(x3(k)+x3(k-1))-x3(k-1))/(x3(k)-x3(k-1))
          tmp3k0=(0.5*(x3(k)+x3(k-1))-x3(k))/(x3(k)-x3(k-1))
          tmp3kp=(0.5*(x3(k)+x3(k-1))-x3(k+1))/(x3(k)-x3(k-1))
         
          urikl2m(n,i,j,k)=uri(n,i,j,k)+d2km*tmp2km+0.5*d3km*tmp2km**2
          urikl20(n,i,j,k)=uri(n,i,j,k)+d2k0*tmp2k0+0.5*d3k0*tmp2k0**2
          urikl2p(n,i,j,k)=uri(n,i,j,k)+d2kp*tmp2kp+0.5*d3kp*tmp2kp**2
         
          urikr2m(n,i,j,k-1)=uri(n,i,j,k)+d2km*tmp3km+0.5*d3km*tmp3km**2
          urikr20(n,i,j,k-1)=uri(n,i,j,k)+d2k0*tmp3k0+0.5*d3k0*tmp3k0**2
          urikr2p(n,i,j,k-1)=uri(n,i,j,k)+d2kp*tmp3kp+0.5*d3kp*tmp3kp**2
         
          dklm(n,i,j,k)=alp1*(urikl2m(n,i,j,k)-urikl1(n,i,j,k))
          dkl0(n,i,j,k)=alp2*(urikl20(n,i,j,k)-urikl1(n,i,j,k))
          dklp(n,i,j,k)=alp3*(urikl2p(n,i,j,k)-urikl1(n,i,j,k))
         
          dkrm(n,i,j,k-1)=alp1*(urikr2m(n,i,j,k-1)-urikr1(n,i,j,k-1))
          dkr0(n,i,j,k-1)=alp2*(urikr20(n,i,j,k-1)-urikr1(n,i,j,k-1))
          dkrp(n,i,j,k-1)=alp3*(urikr2p(n,i,j,k-1)-urikr1(n,i,j,k-1))

        enddo
      enddo
    enddo
  enddo
!
  do k=ks1+2,ke1-2
    do j=js1+2,je1-2
      do i=is1+2,ie1-2
        do n=1,nv

          uriil(n,i,j,k)=uriil1(n,i,j,k)

          if(dilp(n,i,j,k) .gt. 0.d0 .and. dil0(n,i,j,k) .gt. 0.d0 & 
              .and. dilm(n,i,j,k) .gt. 0.d0) then

            uriil(n,i,j,k)=uriil2m(n,i,j,k)

            if(abs(dil0(n,i,j,k)) .lt. abs(dilm(n,i,j,k)) .and. &
               abs(dil0(n,i,j,k)) .lt. abs(dilp(n,i,j,k))) then

              uriil(n,i,j,k)=uriil20(n,i,j,k)

            endif

            if(abs(dilp(n,i,j,k)) .lt. abs(dilm(n,i,j,k)) .and. &
               abs(dilp(n,i,j,k)) .lt. abs(dil0(n,i,j,k))) then
      
              uriil(n,i,j,k)=uriil2p(n,i,j,k)

            endif

          endif

          if(dilm(n,i,j,k) .lt. 0.d0 .and. dil0(n,i,j,k) .lt. 0.d0 &
              .and. dilp(n,i,j,k) .lt. 0.d0) then

            uriil(n,i,j,k)=uriil2m(n,i,j,k)

            if(abs(dil0(n,i,j,k)) .lt. abs(dilm(n,i,j,k)) .and. & 
               abs(dil0(n,i,j,k)) .lt. abs(dilp(n,i,j,k))) then

              uriil(n,i,j,k)=uriil20(n,i,j,k)

            endif

            if(abs(dilp(n,i,j,k)) .lt. abs(dilm(n,i,j,k)) .and. &
               abs(dilp(n,i,j,k)) .lt. abs(dil0(n,i,j,k))) then

              uriil(n,i,j,k)=uriil2p(n,i,j,k)

            endif

          endif
!
          uriir(n,i-1,j,k)=uriir1(n,i-1,j,k)

          if(dirm(n,i-1,j,k) .gt. 0.d0 .and. dir0(n,i-1,j,k) .gt. 0.d0 &
             .and. dirp(n,i-1,j,k) .gt. 0.d0) then

            uriir(n,i-1,j,k)=uriir2m(n,i-1,j,k)

            if(abs(dir0(n,i-1,j,k)) .lt. abs(dirm(n,i-1,j,k)) .and. & 
               abs(dir0(n,i-1,j,k)) .lt. abs(dirp(n,i-1,j,k))) then

              uriir(n,i-1,j,k)=uriir20(n,i-1,j,k)

            endif

            if(abs(dirp(n,i-1,j,k)) .lt. abs(dirm(n,i-1,j,k)) .and. &
               abs(dirp(n,i-1,j,k)) .lt. abs(dir0(n,i-1,j,k))) then

              uriir(n,i-1,j,k)=uriir2p(n,i-1,j,k)

            endif

          endif

          if(dirm(n,i-1,j,k) .lt. 0.d0 .and. dir0(n,i-1,j,k) .lt. 0.d0 &
             .and. dirp(n,i-1,j,k) .lt. 0.d0) then

            uriir(n,i-1,j,k)=uriir2m(n,i-1,j,k)

            if(abs(dir0(n,i-1,j,k)) .lt. abs(dirm(n,i-1,j,k)) .and. &
               abs(dir0(n,i-1,j,k)) .lt. abs(dirp(n,i-1,j,k))) then

              uriir(n,i-1,j,k)=uriir20(n,i-1,j,k)

            endif

            if(abs(dirp(n,i-1,j,k)) .lt. abs(dirm(n,i-1,j,k)) .and. &
               abs(dirp(n,i-1,j,k)) .lt. abs(dir0(n,i-1,j,k))) then

              uriir(n,i-1,j,k)=uriir2p(n,i-1,j,k)

            endif

          endif
!
          urijl(n,i,j,k)=urijl1(n,i,j,k)
           
          if(djlm(n,i,j,k) .gt. 0.d0 .and. djl0(n,i,j,k) .gt. 0.d0 & 
               .and. djlp(n,i,j,k) .gt. 0.d0) then
            urijl(n,i,j,k)=urijl2m(n,i,j,k)
            if(abs(djl0(n,i,j,k)) .lt. abs(djlm(n,i,j,k)) .and. & 
               abs(djl0(n,i,j,k)) .lt. abs(djlp(n,i,j,k))) then
              urijl(n,i,j,k)=urijl20(n,i,j,k)
            endif
            if(abs(djlp(n,i,j,k)) .lt. abs(djlm(n,i,j,k)) .and. & 
               abs(djlp(n,i,j,k)) .lt. abs(djl0(n,i,j,k))) then
              urijl(n,i,j,k)=urijl2p(n,i,j,k)
            endif
          endif

          if(djlm(n,i,j,k) .lt. 0.d0 .and. djl0(n,i,j,k) .lt. 0.d0 &
              .and. djlp(n,i,j,k) .lt. 0.d0) then
            urijl(n,i,j,k)=urijl2m(n,i,j,k)
            if(abs(djl0(n,i,j,k)) .lt. abs(djlm(n,i,j,k)) .and. &
               abs(djl0(n,i,j,k)) .lt. abs(djlp(n,i,j,k))) then
              urijl(n,i,j,k)=urijl20(n,i,j,k)
            endif
            if(abs(djlp(n,i,j,k)) .lt. abs(djlm(n,i,j,k)) .and. &
               abs(djlp(n,i,j,k)) .lt. abs(djl0(n,i,j,k))) then
              urijl(n,i,j,k)=urijl2p(n,i,j,k)
            endif
          endif
!
          urijr(n,i,j-1,k)=urijr1(n,i,j-1,k)
           
          if(djrm(n,i,j-1,k) .gt. 0.d0 .and. djr0(n,i,j-1,k) .gt. 0.d0 &
            .and. djrp(n,i,j-1,k) .gt. 0.d0) then
            urijr(n,i,j-1,k)=urijr2m(n,i,j-1,k)
            if(abs(djr0(n,i,j-1,k)) .lt. abs(djrm(n,i,j-1,k)) .and. &
               abs(djr0(n,i,j-1,k)) .lt. abs(djrp(n,i,j-1,k))) then
              urijr(n,i,j-1,k)=urijr20(n,i,j-1,k)
            endif
            if(abs(djrp(n,i,j-1,k)) .lt. abs(djrm(n,i,j-1,k)) .and. & 
               abs(djrp(n,i,j-1,k)) .lt. abs(djr0(n,i,j-1,k))) then
              urijr(n,i,j-1,k)=urijr2p(n,i,j-1,k)
            endif
          endif

          if(djrm(n,i,j-1,k) .lt. 0.d0 .and. djr0(n,i,j-1,k) .lt. 0.d0 & 
            .and. djrp(n,i,j-1,k) .lt. 0.d0) then
            urijr(n,i,j-1,k)=urijr2m(n,i,j-1,k)
            if(abs(djr0(n,i,j-1,k)) .lt. abs(djrm(n,i,j-1,k)) .and. &
               abs(djr0(n,i,j-1,k)) .lt. abs(djrp(n,i,j-1,k))) then
              urijr(n,i,j-1,k)=urijr20(n,i,j-1,k)
            endif
            if(abs(djrp(n,i,j-1,k)) .lt. abs(djrm(n,i,j-1,k)) .and. &
               abs(djrp(n,i,j-1,k)) .lt. abs(djr0(n,i,j-1,k))) then
              urijr(n,i,j-1,k)=urijr2p(n,i,j-1,k)
            endif
          endif
!
!
          urikl(n,i,j,k)=urikl1(n,i,j,k)
           
          if(dklm(n,i,j,k) .gt. 0.d0 .and. dkl0(n,i,j,k) .gt. 0.d0 & 
             .and. dklp(n,i,j,k) .gt. 0.d0) then
            urikl(n,i,j,k)=urikl2m(n,i,j,k)
            if(abs(dkl0(n,i,j,k)) .lt. abs(dklm(n,i,j,k)) .and. &
               abs(dkl0(n,i,j,k)) .lt. abs(dklp(n,i,j,k))) then
              urikl(n,i,j,k)=urikl20(n,i,j,k)
            endif
            if(abs(dklp(n,i,j,k)) .lt. abs(dklm(n,i,j,k)) .and. & 
               abs(dklp(n,i,j,k)) .lt. abs(dkl0(n,i,j,k))) then
              urikl(n,i,j,k)=urikl2p(n,i,j,k)
            endif
          endif

          if(dklm(n,i,j,k) .lt. 0.d0 .and. dkl0(n,i,j,k) .lt. 0.d0 & 
              .and. dklp(n,i,j,k) .lt. 0.d0) then
            urikl(n,i,j,k)=urikl2m(n,i,j,k)
            if(abs(dkl0(n,i,j,k)) .lt. abs(dklp(n,i,j,k)) .and. & 
               abs(dkl0(n,i,j,k)) .lt. abs(dklm(n,i,j,k))) then
              urikl(n,i,j,k)=urikl20(n,i,j,k)
            endif
            if(abs(dklp(n,i,j,k)) .lt. abs(dklm(n,i,j,k)) .and. &
               abs(dklp(n,i,j,k)) .lt. abs(dklp(n,i,j,k))) then
              urikl(n,i,j,k)=urikl2p(n,i,j,k)
            endif
          endif
!
          urikr(n,i,j,k-1)=urikr1(n,i,j,k-1)
           
          if(dkrm(n,i,j,k-1) .gt. 0.d0 .and. dkr0(n,i,j,k-1) .gt. 0.d0 & 
            .and. dkrp(n,i,j,k-1) .gt. 0.d0) then
            urikr(n,i,j,k-1)=urikr2m(n,i,j,k-1)
            if(abs(dkr0(n,i,j,k-1)) .lt. abs(dkrm(n,i,j,k-1)) .and. & 
               abs(dkr0(n,i,j,k-1)) .lt. abs(dkrp(n,i,j,k-1))) then
              urikr(n,i,j,k-1)=urikr20(n,i,j,k-1)
            endif
            if(abs(dkrp(n,i,j,k-1)) .lt. abs(dkrm(n,i,j,k-1)) .and. & 
               abs(dkrp(n,i,j,k-1)) .lt. abs(dkr0(n,i,j,k-1))) then
              urikr(n,i,j,k-1)=urikr2p(n,i,j,k-1)
            endif
          endif

          if(dkrm(n,i,j,k-1) .lt. 0.d0 .and. dkr0(n,i,j,k-1) .lt. 0.d0 &
            .and. dkrp(n,i,j,k-1) .lt. 0.d0) then
            urikr(n,i,j,k-1)=urikr2m(n,i,j,k-1)
            if(abs(dkr0(n,i,j,k-1)) .lt. abs(dkrp(n,i,j,k-1)) .and. &
               abs(dkr0(n,i,j,k-1)) .lt. abs(dkrm(n,i,j,k-1))) then
              urikr(n,i,j,k-1)=urikr20(n,i,j,k-1)
            endif
            if(abs(dkrp(n,i,j,k-1)) .lt. abs(dkrm(n,i,j,k-1)) .and. & 
               abs(dkrp(n,i,j,k-1)) .lt. abs(dkr0(n,i,j,k-1))) then
              urikr(n,i,j,k-1)=urikr2p(n,i,j,k-1)
            endif
          endif
           
        enddo
      enddo
    enddo
  enddo
!
  deallocate( uriir1, uriil1, urijr1, urijl1, urikr1, urikl1, &
              uriir2m, uriil2m, urijr2m, urijl2m, urikr2m, urikl2m, &
              uriir20, uriil20, urijr20, urijl20, urikr20, urikl20, &   
              uriir2p, uriil2p, urijr2p, urijl2p, urikr2p, urikl2p, &
              dirm, dilm, djrm, djlm, dkrm, dklm, dir0, dil0, &
              djr0, djl0, dkr0, dkl0, dirp, dilp, djrp, djlp, &
              dkrp, dklp, stat=merr)
!
  return
end subroutine ceno 
!
!---------------------------------------------------------------------@
subroutine ppm(uri,x1,x2,x3,dx1b,dx2b,dx3b, &
               uriir,urijr,urikr,uriil,urijl,urikl, &
               is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Reconstruction step with Piecewise Parabolic method (PPM) 
!     Calculate cell-interface variables(uriil, uriir) 
!     from cell-center variables (uri)
!
  use pram, only : imax, jmax, kmax, nv
  implicit none
!
  integer :: is1, ie1, js1, je1, ks1, ke1, merr
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: dmai(:,:,:,:), dmaj(:,:,:,:), dmak(:,:,:,:)
!
  real(8) :: x1(imax), x2(jmax), x3(kmax)
  real(8) :: dx1b(2:imax-1), dx2b(2:jmax-1), dx3b(2:kmax-1)
!
  allocate( dmai(nv,is1:ie1,js1:je1,ks1:ke1), &
            dmaj(nv,is1:ie1,js1:je1,ks1:ke1), &
            dmak(nv,is1:ie1,js1:je1,ks1:ke1), stat=merr)
!
!-----------------------------------------------------------------------
      
  call interp(uri,dx1b,dx2b,dx3b,uriir,urijr,urikr,uriil,urijl,urikl, &
              dmai,dmaj,dmak,is1,ie1,js1,je1,ks1,ke1)
  call detect(uri,x1,x2,x3,dx1b,dx2b,dx3b, &
              uriir,urijr,urikr,uriil,urijl,urikl,dmai,dmaj,dmak, &
              is1,ie1,js1,je1,ks1,ke1)
  call flaten(uri,uriir,urijr,urikr,uriil,urijl,urikl, &
              is1,ie1,js1,je1,ks1,ke1)      
  call monoto(uri,uriir,urijr,urikr,uriil,urijl,urikl, &
              is1,ie1,js1,je1,ks1,ke1)
!
  deallocate( dmai, dmaj, dmak, stat=merr)
!
  return
end subroutine ppm
!
!---------------------------------------------------------------------@
subroutine interp(uri,dx1b,dx2b,dx3b,uriir,urijr,urikr,uriil,urijl,urikl, &
                  dmai,dmaj,dmak,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!      Step1 for rPPM method (Interpolation)
!        from Marti & Muller 1996,JCP, 123,1
!
  use pram, only : imax, jmax, kmax, nv
  implicit none

  integer :: i, j, k, n, is1, ie1, js1, je1, ks1, ke1
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: dmai(nv,is1:ie1,js1:je1,ks1:ke1), &
             dmaj(nv,is1:ie1,js1:je1,ks1:ke1), &
             dmak(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: dx1b(2:imax-1), dx2b(2:jmax-1), dx3b(2:kmax-1)
!
  real(8) :: tmp1i, tmp1j, tmp1k, &
             tmp1ai, tmp1bi, tmp1ci, tmp1di, sp1i, sm1i, spm1i, dai, &
             tmp1aj, tmp1bj, tmp1cj, tmp1dj, sp1j, sm1j, spm1j, daj, &
             tmp1ak, tmp1bk, tmp1ck, tmp1dk, sp1k, sm1k, spm1k, dak, &
             tmp2i, tmp3ai, tmp3bi, tmp3ci, tmp3i, tmp4i, tmp5i, &
             tmp2j, tmp3aj, tmp3bj, tmp3cj, tmp3j, tmp4j, tmp5j, &
             tmp2k, tmp3ak, tmp3bk, tmp3ck, tmp3k, tmp4k, tmp5k
!     
!-----------------------------------------------------------------------
           
  do k=ks1+1,ke1-1
    do j=js1+1,je1-1
      do i=is1+1,ie1-1
        do n=1,nv
!
!  Equation  (60) & (61) of Marti & Muller(1996)
!
          tmp1ai=dx1b(i)/(dx1b(i-1)+dx1b(i)+dx1b(i+1))
          tmp1bi=(2.0*dx1b(i-1)+dx1b(i))/(dx1b(i+1)+dx1b(i))
          tmp1ci=(dx1b(i)+2.0*dx1b(i+1))/(dx1b(i-1)+dx1b(i))
          
          sp1i=(uri(n,i+1,j,k)-uri(n,i,j,k))
          sm1i=(uri(n,i,j,k)-uri(n,i-1,j,k))
          spm1i=sp1i*sm1i
          
          dai=tmp1ai*(tmp1bi*sp1i+tmp1ci*sm1i)
          
          if(spm1i .le. 0.d0) then
            dmai(n,i,j,k)=0.d0
          else
            if(dai .lt. 0.d0) then
              tmp1i=-1.d0
            else
              tmp1i=1.d0
            endif
            tmp1di=min(abs(sm1i),abs(sp1i))
            dmai(n,i,j,k)=tmp1i*min(abs(dai),2.0*tmp1di)
          endif
          
          tmp1aj=dx2b(j)/(dx2b(j-1)+dx2b(j)+dx2b(j+1))
          tmp1bj=(2.0*dx2b(j-1)+dx2b(j))/(dx2b(j+1)+dx2b(j))
          tmp1cj=(dx2b(j)+2.0*dx2b(j+1))/(dx2b(j-1)+dx2b(j))
     
          sp1j=(uri(n,i,j+1,k)-uri(n,i,j,k))
          sm1j=(uri(n,i,j,k)-uri(n,i,j-1,k))
          spm1j=sp1j*sm1j
     
          daj=tmp1aj*(tmp1bj*sp1j+tmp1cj*sm1j)
          
          if(spm1j .le. 0.d0) then
            dmaj(n,i,j,k)=0.d0
          else
            if(daj .lt. 0.d0) then
              tmp1j=-1.d0
            else
              tmp1j=1.d0
            endif
            tmp1dj=min(abs(sm1j),abs(sp1j))
            dmaj(n,i,j,k)=tmp1j*min(abs(daj),2.0*tmp1dj)
          endif
      
          tmp1ak=dx3b(k)/(dx3b(k-1)+dx3b(k)+dx3b(k+1))
          tmp1bk=(2.0*dx3b(k-1)+dx3b(k))/(dx3b(k+1)+dx3b(k))
          tmp1ck=(dx3b(k)+2.0*dx3b(k+1))/(dx3b(k-1)+dx3b(k))
      
          sp1k=(uri(n,i,j,k+1)-uri(n,i,j,k))
          sm1k=(uri(n,i,j,k)-uri(n,i,j,k-1))
          spm1k=sp1k*sm1k
          
          dak=tmp1ak*(tmp1bk*sp1k+tmp1ck*sm1k)
          
          if(spm1k .le. 0.d0) then
            dmak(n,i,j,k)=0.d0
          else
            if(dak .lt. 0.d0) then
              tmp1k=-1.d0
            else
              tmp1k=1.d0
            endif
            tmp1dk=min(abs(sm1k),abs(sp1k))
            dmak(n,i,j,k)=tmp1k*min(abs(dak),tmp1dk)
          endif
       
        enddo
      enddo
    enddo
  enddo

  do k=ks1+1,ke1-2
    do j=js1+1,je1-2
      do i=is1+1,ie1-2
        do n=1,nv
      
!
!  Equation (51) of Marti & Muller(1996)
!
          tmp2i=(dx1b(i)/(dx1b(i)+dx1b(i+1)))*(uri(n,i+1,j,k)-uri(n,i,j,k))
          tmp3ai=1.0/(dx1b(i-1)+dx1b(i)+dx1b(i+1)+dx1b(i+2))
          tmp3bi=(2.0*dx1b(i+1)*dx1b(i))/(dx1b(i)+dx1b(i+1))
          tmp3ci=((dx1b(i-1)+dx1b(i))/(2.0*dx1b(i)+dx1b(i+1))) &
                -((dx1b(i+2)+dx1b(i+1))/(2.0*dx1b(i+1)+dx1b(i)))
          tmp3i=tmp3ai*tmp3bi*tmp3ci*(uri(n,i+1,j,k)-uri(n,i,j,k))
          tmp4i=tmp3ai*dx1b(i)*((dx1b(i-1)+dx1b(i))/(2.0*dx1b(i)+dx1b(i+1)))
          tmp5i=tmp3ai*dx1b(i+1)*((dx1b(i+1)+dx1b(i+2))/(dx1b(i)+2.0*dx1b(i+1)))
          
          uriir(n,i,j,k)=uri(n,i,j,k)+tmp2i+tmp3i-tmp4i*dmai(n,i+1,j,k) &
                        +tmp5i*dmai(n,i,j,k)
          uriil(n,i,j,k)=uriir(n,i,j,k)
          
          tmp2j=(dx2b(j)/(dx2b(j)+dx2b(j+1)))*(uri(n,i,j+1,k)-uri(n,i,j,k))
          tmp3aj=1.0/(dx2b(j-1)+dx2b(j)+dx2b(j+1)+dx2b(j+2))
          tmp3bj=(2.0*dx2b(j+1)*dx2b(j))/(dx2b(j)+dx2b(j+1))
          tmp3cj=((dx2b(j-1)+dx2b(j))/(2.0*dx2b(j)+dx2b(j+1))) &
                -((dx2b(j+2)+dx2b(j+1))/(2.0*dx2b(j+1)+dx2b(j)))
          tmp3j=tmp3aj*tmp3bj*tmp3cj*(uri(n,i,j+1,k)-uri(n,i,j,k))
          tmp4j=tmp3aj*dx2b(j)*((dx2b(j-1)+dx2b(j))/(2.0*dx2b(j)+dx2b(j+1)))
          tmp5j=tmp3aj*dx2b(j+1)*((dx2b(j+1)+dx2b(j+2))/(dx2b(j)+2.0*dx2b(j+1)))
          
          urijr(n,i,j,k)=uri(n,i,j,k)+tmp2j+tmp3j-tmp4j*dmaj(n,i,j+1,k) &
                        +tmp5j*dmaj(n,i,j,k)
          urijl(n,i,j,k)=urijr(n,i,j,k)
          
          
          tmp2k=(dx3b(k)/(dx3b(k)+dx3b(k+1)))*(uri(n,i,j,k+1)-uri(n,i,j,k))
          tmp3ak=1.0/(dx3b(k-1)+dx3b(k)+dx3b(k+1)+dx3b(k+2))
          tmp3bk=(2.0*dx3b(k+1)*dx3b(k))/(dx3b(k)+dx3b(k+1))
          tmp3ck=((dx3b(k-1)+dx3b(k))/(2.0*dx3b(k)+dx3b(k+1))) &
                -((dx3b(k+2)+dx3b(k+1))/(2.0*dx3b(k+1)+dx3b(k)))
          tmp3k=tmp3ak*tmp3bk*tmp3ck*(uri(n,i,j,k+1)-uri(n,i,j,k))
          tmp4k=tmp3ak*dx3b(k)*((dx3b(k-1)+dx3b(k))/(2.0*dx3b(k)+dx3b(k+1)))
          tmp5k=tmp3ak*dx3b(k+1)*((dx3b(k+1)+dx3b(k+2))/(dx3b(k)+2.0*dx3b(k+1)))
          
          urikr(n,i,j,k)=uri(n,i,j,k)+tmp2k+tmp3k-tmp4k*dmak(n,i,j,k+1) &
                        +tmp5k*dmak(n,i,j,k)
          urikl(n,i,j,k)=urikr(n,i,j,k)
          
        enddo
      enddo
    enddo
  enddo

  return
end subroutine interp
!
!---------------------------------------------------------------------@
subroutine detect(uri,x1,x2,x3,dx1b,dx2b,dx3b, &
                  uriir,urijr,urikr,uriil,urijl,urikl,dmai,dmaj,dmak, &
                  is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!      Step2 for rPPM method (Interpolation)
!        from Marti & Muller 1996,JCP, 123,1  
!
  use pram, only : imax, jmax, kmax, nv, gam, c0
  implicit none
!
  integer :: i, j, k, is1, ie1, js1, je1, ks1, ke1, merr

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: dmai(nv,is1:ie1,js1:je1,ks1:ke1), &
             dmaj(nv,is1:ie1,js1:je1,ks1:ke1), &
             dmak(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: ddai(:,:,:), ddaj(:,:,:), ddak(:,:,:)
                
  real(8), allocatable :: etai(:,:,:), etaj(:,:,:), etak(:,:,:)

  real(8) :: x1(imax), x2(jmax), x3(kmax)
  real(8) :: dx1b(2:imax-1), dx2b(2:jmax-1), dx3b(2:kmax-1)
!
  real(8) :: eta1, eta2, eps1, ak0
  real(8) :: tmp0ai, tmp0bi, tmp0ci, tmp1i, tmp2i, tmp2ai, tmp2bi, etatili, &
             tmp0aj, tmp0bj, tmp0cj, tmp1j, tmp2j, tmp2aj, tmp2bj, etatilj, &
             tmp0ak, tmp0bk, tmp0ck, tmp1k, tmp2k, tmp2ak, tmp2bk, etatilk
  real(8) :: tmp3i, tmp4i, tmp3j, tmp4j, tmp3k, tmp4k
!
  allocate( ddai(is1:ie1,js1:je1,ks1:ke1), &
            ddaj(is1:ie1,js1:je1,ks1:ke1), &
            ddak(is1:ie1,js1:je1,ks1:ke1), &
            etai(is1:ie1,js1:je1,ks1:ke1), &
            etaj(is1:ie1,js1:je1,ks1:ke1), &
            etak(is1:ie1,js1:je1,ks1:ke1), stat=merr)
!     
!-----------------------------------------------------------------------
!  Parameter
!-----------------------------------------------------------------------
!
  eta1=5.d0
  eta2=0.05d0
  eps1=0.1d0
  ak0=1.d0

!-----------------------------------------------------------------------

  do k=ks1+1,ke1-1
    do j=js1+1,je1-1
      do i=is1+1,ie1-1
!
!  Equation (69)
!
        tmp0ai=1.0/(dx1b(i-1)+dx1b(i)+dx1b(i+1))
        tmp0bi=(uri(1,i+1,j,k)-uri(1,i,j,k))/(dx1b(i+1)+dx1b(i))
        tmp0ci=(uri(1,i,j,k)-uri(1,i-1,j,k))/(dx1b(i)+dx1b(i-1))
          
        ddai(i,j,k)=tmp0ai*(tmp0bi-tmp0ci)
     
        tmp0aj=1.0/(dx2b(j-1)+dx2b(j)+dx2b(j+1))
        tmp0bj=(uri(1,i,j+1,k)-uri(1,i,j,k))/(dx2b(j+1)+dx2b(j))
        tmp0cj=(uri(1,i,j,k)-uri(1,i,j-1,k))/(dx2b(j)+dx2b(j-1))
          
        ddaj(i,j,k)=tmp0aj*(tmp0bj-tmp0cj)

        tmp0ak=1.0/(dx3b(k-1)+dx3b(k)+dx3b(k+1))
        tmp0bk=(uri(1,i,j,k+1)-uri(1,i,j,k))/(dx3b(k+1)+dx3b(k))
        tmp0ck=(uri(1,i,j,k)-uri(1,i,j,k-1))/(dx3b(k)+dx3b(k-1))

        ddak(i,j,k)=tmp0ak*(tmp0bk-tmp0ck)

      enddo
    enddo
  enddo

  do k=ks1+2,ke1-2
    do j=js1+2,je1-2
      do i=is1+2,ie1-2
!
!  Equation (66), (67), (68)
!
        tmp1i=ddai(i+1,j,k)*ddai(i-1,j,k)
        tmp2i=abs(uri(1,i+1,j,k)-uri(1,i-1,j,k)) &
                -eps1*min(abs(uri(1,i+1,j,k)), abs(uri(1,i-1,j,k)))

        if(tmp1i .ge. 0.d0 .and. tmp2i .le. 0.d0) then
          etatili=0.d0
        else
          tmp2ai=-(ddai(i+1,j,k)-ddai(i-1,j,k))/(x1(i+1)-x1(i-1))
          tmp2bi=((x1(i)-x1(i-1))**3+(x1(i+1)-x1(i))**3) &
                 /(uri(1,i+1,j,k)-uri(1,i-1,j,k))
          etatili=tmp2ai*tmp2bi
        endif

        etai(i,j,k)=max(0.d0, min(eta1*(etatili-eta2),1.d0))

        tmp1j=ddaj(i,j+1,k)*ddaj(i,j-1,k)
        tmp2j=abs(uri(1,i,j+1,k)-uri(1,i,j-1,k)) &
                -eps1*min(abs(uri(1,i,j+1,k)), abs(uri(1,i,j-1,k)))

        if(tmp1j .ge. 0.d0 .and. tmp2j .le. 0.d0) then
          etatilj=0.d0
        else
          tmp2aj=-(ddaj(i,j+1,k)-ddaj(i,j-1,k))/(x2(j+1)-x2(j-1))
          tmp2bj=((x2(j)-x2(j-1))**3+(x2(j+1)-x2(j))**3) &
                 /(uri(1,i,j+1,k)-uri(1,i,j-1,k))
          etatilj=tmp2aj*tmp2bj
        endif

        etaj(i,j,k)=max(0.d0, min(eta1*(etatilj-eta2),1.d0))

        tmp1k=ddak(i,j,k+1)*ddak(i,j,k-1)
        tmp2k=abs(uri(1,i,j,k+1)-uri(1,i,j,k-1)) &
                -eps1*min(abs(uri(1,i,j,k+1)), abs(uri(1,i,j,k-1)))

        if(tmp1k .ge. 0.d0 .and. tmp2k .le. 0.d0) then
          etatilk=0.d0
        else
          tmp2ak=-(ddak(i,j,k+1)-ddak(i,j,k-1))/(x3(k+1)-x3(k-1))
          tmp2bk=((x3(k)-x3(k-1))**3+(x3(k+1)-x3(k))**3) &
                 /(uri(1,i,j,k+1)-uri(1,i,j,k-1))
          etatilk=tmp2ak*tmp2bk
        endif

        etak(i,j,k)=max(0.d0, min(eta1*(etatilk-eta2),1.d0))

      enddo
    enddo
  enddo

  do k=ks1+2,ke1-2
    do j=js1+2,je1-2
      do i=is1+2,ie1-2
!
!  Equation (63)
!
        tmp3i=abs(uri(1,i+1,j,k)-uri(1,i-1,j,k)) &
              /min(uri(1,i+1,j,k),uri(1,i-1,j,k))
        tmp4i=abs(uri(5,i+1,j,k)-uri(5,i-1,j,k)) &
              /min(uri(5,i+1,j,k),uri(5,i-1,j,k))

        if(gam*ak0*tmp3i .lt. tmp4i) then
          etai(i,j,k)=0.d0
        endif
         
        tmp3j=abs(uri(1,i,j+1,k)-uri(1,i,j-1,k)) &
              /min(uri(1,i,j+1,k),uri(1,i,j-1,k))
        tmp4j=abs(uri(5,i,j+1,k)-uri(5,i,j-1,k)) &
              /min(uri(5,i,j+1,k),uri(5,i,j-1,k))

        if(gam*ak0*tmp3j .lt. tmp4j) then
          etaj(i,j,k)=0.d0
        endif

        tmp3k=abs(uri(1,i,j,k+1)-uri(1,i,j,k-1)) &
              /min(uri(1,i,j,k+1),uri(1,i,j,k-1))
        tmp4k=abs(uri(5,i,j,k+1)-uri(5,i,j,k-1)) &
              /min(uri(5,i,j,k+1),uri(5,i,j,k-1))

        if(gam*ak0*tmp3k .lt. tmp4k) then
          etak(i,j,k)=0.d0
        endif

      enddo
    enddo
  enddo

  do k=ks1+2,ke1-2
    do j=js1+2,je1-2
      do i=is1+2,ie1-2
!
!  Equation (64), (65)
!
        uriil(1,i,j,k)=uriil(1,i,j,k)*(1.0-etai(i,j,k)) &
                       +(uri(1,i+1,j,k)-0.5*dmai(1,i+1,j,k))*etai(i,j,k)
    
        urijl(1,i,j,k)=urijl(1,i,j,k)*(1.0-etaj(i,j,k)) &
                       +(uri(1,i,j+1,k)-0.5*dmaj(1,i,j+1,k))*etaj(i,j,k)
    
        urikl(1,i,j,k)=urikl(1,i,j,k)*(1.0-etak(i,j,k)) &
                       +(uri(1,i,j,k+1)-0.5*dmak(1,i,j,k+1))*etak(i,j,k)
    
        uriir(1,i-1,j,k)=uriir(1,i-1,j,k)*(1.0-etai(i,j,k)) &
                        +(uri(1,i-1,j,k)+0.5*dmai(1,i-1,j,k))*etai(i,j,k)

        urijr(1,i,j-1,k)=urijr(1,i,j-1,k)*(1.0-etaj(i,j,k)) &
                        +(uri(1,i,j-1,k)+0.5*dmaj(1,i,j-1,k))*etaj(i,j,k)
    
        urikr(1,i,j,k-1)=urikr(1,i,j,k-1)*(1.0-etak(i,j,k)) &
                        +(uri(1,i,j,k-1)+0.5*dmak(1,i,j,k-1))*etak(i,j,k)

      enddo
    enddo
  enddo
!
  deallocate( ddai, ddaj, ddak, etai, etaj, etak, stat=merr)
!
  return
end subroutine detect
!
!---------------------------------------------------------------------@
subroutine flaten(uri,uriir,urijr,urikr,uriil,urijl,urikl, &
                  is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!      Step3 for rPPM method (Flattening)
!        from Marti & Muller 1996,JCP, 123,1  
!
  use pram, only : imax, jmax, kmax, nv
  implicit none
!
  integer :: i, j, k, n, is1, ie1, js1, je1, ks1, ke1, merr

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8), allocatable :: f1ai(:,:,:), f1aj(:,:,:), f1ak(:,:,:)
     
  real(8), allocatable :: f2ai(:,:,:), f2aj(:,:,:), f2ak(:,:,:)
 
  real(8) :: eps2, omg2, omg3
  real(8) :: tmp1i, tmp2i, tmp3ia, tmp3ib, tmp4i, omg1i, &
             tmp1j, tmp2j, tmp3ja, tmp3jb, tmp4j, omg1j, &
             tmp1k, tmp2k, tmp3ka, tmp3kb, tmp4k, omg1k, &
             f1bi, f1bj, f1bk
!
  allocate( f1ai(is1:ie1,js1:je1,ks1:ke1), &
            f1aj(is1:ie1,js1:je1,ks1:ke1), &
            f1ak(is1:ie1,js1:je1,ks1:ke1), &
            f2ai(is1:ie1,js1:je1,ks1:ke1), &
            f2aj(is1:ie1,js1:je1,ks1:ke1), &
            f2ak(is1:ie1,js1:je1,ks1:ke1), stat=merr)    
!
!-----------------------------------------------------------------------
!  Parameter
!-----------------------------------------------------------------------
!
  eps2=1.d0
  omg2=0.52d0
  omg3=10.d0

!-----------------------------------------------------------------------

  do k=ks1+2,ke1-2
    do j=js1+2,je1-2
      do i=is1+2,ie1-2
!
!  Equation (71), (72)
!

!        tmp3ia=(uri(5,i+1,j,k)-uri(5,i-1,j,k))
        tmp3ia=abs(uri(5,i+1,j,k)-uri(5,i-1,j,k))
        tmp1i=eps2*min(uri(5,i+1,j,k), uri(5,i-1,j,k))-tmp3ia
        tmp2i=uri(2,i+1,j,k)-uri(2,i-1,j,k)
         
        if(tmp1i .lt. 0.d0 .and. tmp2i .lt. 0.d0) then
          omg1i=1.d0
        else
          omg1i=0.d0
        endif
          
        tmp3ib=(uri(5,i+2,j,k)-uri(5,i-2,j,k))
    
        if(tmp3ib .eq. 0.d0) then
          if(tmp3ia .eq. 0.d0) then
            tmp4i=-omg2
          else
            tmp4i=1.0-omg2
          endif
        else
          tmp4i=(tmp3ia/tmp3ib)-omg2
        endif

        f1ai(i,j,k)=min(1.d0,omg1i*max(0.d0,tmp4i*omg3))
                    
!        tmp3ja=(uri(5,i,j+1,k)-uri(5,i,j-1,k))
        tmp3ja=abs(uri(5,i,j+1,k)-uri(5,i,j-1,k))
        tmp1j=eps2*min(uri(5,i,j+1,k), uri(5,i,j-1,k))-tmp3ja
        tmp2j=uri(3,i,j+1,k)-uri(3,i,j-1,k)
          
        if(tmp1j .lt. 0.d0 .and. tmp2j .lt. 0.d0) then
          omg1j=1.d0
        else
          omg1j=0.d0
        endif
         
        tmp3jb=(uri(5,i,j+2,k)-uri(5,i,j-2,k))
   
        if(tmp3jb .eq. 0.d0) then
          if(tmp3ja .eq. 0.d0) then
            tmp4j=-omg2
          else
            tmp4j=1.0-omg2
          endif
        else
          tmp4j=(tmp3ja/tmp3jb)-omg2
        endif

        f1aj(i,j,k)=min(1.d0,omg1j*max(0.d0,tmp4j*omg3))
          
!        tmp3ka=(uri(5,i,j,k+1)-uri(5,i,j,k-1))
        tmp3ka=abs(uri(5,i,j,k+1)-uri(5,i,j,k-1))
        tmp1k=eps2*min(uri(5,i,j,k+1), uri(5,i,j,k-1))-tmp3ka
        tmp2k=uri(4,i,j,k+1)-uri(4,i,j,k-1)
          
        if(tmp1k .lt. 0.d0 .and. tmp2k .lt. 0.d0) then
          omg1k=1.d0
        else
          omg1k=0.d0
        endif
          
        tmp3kb=(uri(5,i,j,k+2)-uri(5,i,j,k-2))
    
        if(tmp3kb .eq. 0.d0) then
          if(tmp3ka .eq. 0.d0) then
            tmp4k=-omg2
          else
            tmp4k=1.0-omg2
          endif
        else
          tmp4k=(tmp3ka/tmp3kb)-omg2
        endif
    
        f1ak(i,j,k)=min(1.d0,omg1k*max(0.d0,tmp4k*omg3))
     
      enddo
    enddo
  enddo

  do k=ks1+3,ke1-3
    do j=js1+3,je1-3
      do i=is1+3,ie1-3
!
!  Equation (70)
!
        if( uri(5,i+1,j,k)-uri(5,i-1,j,k) .lt. 0.d0) then
          f1bi=f1ai(i+1,j,k)
        else
          f1bi=f1ai(i-1,j,k)
        endif
         
        f2ai(i,j,k)=max(f1ai(i,j,k), f1bi)
          
        if( uri(5,i,j+1,k)-uri(5,i,j-1,k) .lt. 0.d0) then
          f1bj=f1aj(i,j+1,k)
        else
          f1bj=f1aj(i,j-1,k)
        endif
          
        f2aj(i,j,k)=max(f1aj(i,j,k), f1bj)
        
        if( uri(5,i,j,k+1)-uri(5,i,j,k-1) .lt. 0.d0) then
          f1bk=f1ak(i,j,k+1)
        else
          f1bk=f1ak(i,j,k-1)
        endif
          
        f2ak(i,j,k)=max(f1ak(i,j,k), f1bk)
        
      enddo
    enddo
  enddo

  do k=ks1+3,ke1-3
    do j=js1+3,je1-3
      do i=is1+3,ie1-3
        do n=1,nv

          uriil(n,i,j,k)=uri(n,i,j,k)*f2ai(i,j,k) &
                        +uriil(n,i,j,k)*(1.0-f2ai(i,j,k))
          uriir(n,i-1,j,k)=uri(n,i,j,k)*f2ai(i,j,k) &
                          +uriir(n,i-1,j,k)*(1.0-f2ai(i,j,k))
          
          urijl(n,i,j,k)=uri(n,i,j,k)*f2aj(i,j,k) &
                        +urijl(n,i,j,k)*(1.0-f2aj(i,j,k))
          urijr(n,i,j-1,k)=uri(n,i,j,k)*f2aj(i,j,k) &
                          +urijr(n,i,j-1,k)*(1.0-f2aj(i,j,k))
          
          urikl(n,i,j,k)=uri(n,i,j,k)*f2ak(i,j,k) &
                        +urikl(n,i,j,k)*(1.0-f2ak(i,j,k))
          urikr(n,i,j,k-1)=uri(n,i,j,k)*f2ak(i,j,k) &
                          +urikr(n,i,j,k-1)*(1.0-f2ak(i,j,k))
  
        enddo
      enddo
    enddo
  enddo
!
  deallocate( f1ai, f1aj, f1ak, f2ai, f2aj, f2ak, stat=merr) 
!
  return
end subroutine flaten
!
!---------------------------------------------------------------------@
subroutine monoto(uri,uriir,urijr,urikr,uriil,urijl,urikl, &
                  is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!      Step4 for rPPM method (Monotonization)
!        from Marti & Muller 1996,JCP, 123,1 
!
  use pram, only : imax, jmax, nv, gam, c0
  implicit none
!
  integer :: i, j, k, n, is1, ie1, js1, je1, ks1, ke1, merr

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8), allocatable :: uriir1(:,:,:,:), uriil1(:,:,:,:), &
           urijr1(:,:,:,:), urijl1(:,:,:,:), &
           urikr1(:,:,:,:), urikl1(:,:,:,:)

  real(8) :: tmp3i, tmp4i, tmp3j, tmp4j, tmp3k, tmp4k     
!
  allocate( uriir1(nv,is1:ie1,js1:je1,ks1:ke1), &
            uriil1(nv,is1:ie1,js1:je1,ks1:ke1), &
            urijr1(nv,is1:ie1,js1:je1,ks1:ke1), &
            urijl1(nv,is1:ie1,js1:je1,ks1:ke1), &
            urikr1(nv,is1:ie1,js1:je1,ks1:ke1), &
            urikl1(nv,is1:ie1,js1:je1,ks1:ke1), &
            stat=merr )
!
!-----------------------------------------------------------------------

  do k=ks1+2,ke1-2
    do j=js1+2,je1-2
      do i=is1+2,ie1-2
        do n=1,nv
!
!  Equation (73)
!
          if( (uriil(n,i,j,k)-uri(n,i,j,k)) &
              *(uriir(n,i-1,j,k)-uri(n,i,j,k)) .ge. 0.d0 ) then
            uriil(n,i,j,k)=uri(n,i,j,k)
            uriir(n,i-1,j,k)=uri(n,i,j,k)
          endif
           
          if( (urijl(n,i,j,k)-uri(n,i,j,k)) &
              *(urijr(n,i,j-1,k)-uri(n,i,j,k)) .ge. 0.d0 ) then
            urijl(n,i,j,k)=uri(n,i,j,k)
            urijr(n,i,j-1,k)=uri(n,i,j,k)
          endif
           
          if( (urikl(n,i,j,k)-uri(n,i,j,k)) &
              *(urikr(n,i,j,k-1)-uri(n,i,j,k)) .ge. 0.d0 ) then
            urikl(n,i,j,k)=uri(n,i,j,k)
            urikr(n,i,j,k-1)=uri(n,i,j,k)
          endif
 
        enddo
      enddo
    enddo
  enddo

  do k=ks1+2,ke1-2
    do j=js1+2,je1-2
      do i=is1+2,ie1-2
        do n=1,nv
     
          if( (uriil(n,i,j,k)-uri(n,i,j,k)) &
              *(uri(n,i,j,k)-uriir(n,i-1,j,k)) .eq. 0.d0 ) then
            uriil1(n,i,j,k)=uriil(n,i,j,k)
            uriir1(n,i-1,j,k)=uriir(n,i-1,j,k)
          else
            uriil1(n,i,j,k)=3.0*uri(n,i,j,k)-2.0*uriir(n,i-1,j,k)
            uriir1(n,i-1,j,k)=3.0*uri(n,i,j,k)-2.0*uriil(n,i,j,k)
          endif
           
          if( (urijl(n,i,j,k)-uri(n,i,j,k)) &
              *(uri(n,i,j,k)-urijr(n,i,j-1,k)) .eq. 0.d0 ) then
            urijl1(n,i,j,k)=urijl(n,i,j,k)
            urijr1(n,i,j-1,k)=urijr(n,i,j-1,k)
          else
            urijl1(n,i,j,k)=3.0*uri(n,i,j,k)-2.0*urijr(n,i,j-1,k)
            urijr1(n,i,j-1,k)=3.0*uri(n,i,j,k)-2.0*urijl(n,i,j,k)
          endif
     
          if( (urikl(n,i,j,k)-uri(n,i,j,k)) &
              *(uri(n,i,j,k)-urikr(n,i,j,k-1)) .eq. 0.d0 ) then
            urikl1(n,i,j,k)=urikl(n,i,j,k)
            urikr1(n,i,j,k-1)=urikr(n,i,j,k-1)
          else
            urikl1(n,i,j,k)=3.0*uri(n,i,j,k)-2.0*urikr(n,i,j,k-1)
            urikr1(n,i,j,k-1)=3.0*uri(n,i,j,k)-2.0*urikl(n,i,j,k)
          endif
           
        enddo
      enddo
    enddo
  enddo

  do k=ks1+2,ke1-2
    do j=js1+2,je1-2
      do i=is1+2,ie1-2
        do n=1,nv
!
!  Equation (74), (75)
!
          tmp3i=6.0*(uriil(n,i,j,k)-uriir(n,i-1,j,k)) &
               *(uri(n,i,j,k)-0.5*(uriil(n,i,j,k)+uriir(n,i-1,j,k)))
          tmp4i=(uriil(n,i,j,k)-uriir(n,i-1,j,k)) &
               *(uriil(n,i,j,k)-uriir(n,i-1,j,k))

!          tmp3i=(uriil(n,i,j,k)-uriir(n,i-1,j,k)) &
!               *(uriir(n,i-1,j,k)-uriir1(n,i-1,j,k))
!          tmp4i=(uriil(n,i,j,k)-uriir(n,i-1,j,k)) &
!               *(uriil1(n,i,j,k)-uriil(n,i,j,k))
          
          if(tmp3i .gt. tmp4i) then
            uriir(n,i-1,j,k)=3.0*uri(n,i,j,k)-2.0*uriil(n,i,j,k)
          endif
          if(-tmp3i .gt. tmp4i) then
            uriil(n,i,j,k)=3.0*uri(n,i,j,k)-2.0*uriir(n,i-1,j,k)
          endif
          
          tmp3j=6.0*(urijl(n,i,j,k)-urijr(n,i,j-1,k)) &
               *(uri(n,i,j,k)-0.5*(urijl(n,i,j,k)+urijr(n,i,j-1,k)))
          tmp4j=(urijl(n,i,j,k)-urijr(n,i,j-1,k)) &
               *(urijl(n,i,j,k)-urijr(n,i,j-1,k))
          
!          tmp3j=(urijl(n,i,j,k)-urijr(n,i,j-1,k)) &
!               *(urijr(n,i,j-1,k)-urijr1(n,i,j-1,k))
!          tmp4j=(urijl(n,i,j,k)-urijr(n,i,j-1,k)) &
!               *(urijl1(n,i,j,k)-urijl(n,i,j,k))
          
          if(tmp3j .gt. tmp4j) then
            urijr(n,i,j-1,k)=3.0*uri(n,i,j,k)-2.0*urijl(n,i,j,k)
          endif
          if(-tmp3j .gt. tmp4j) then
            urijl(n,i,j,k)=3.0*uri(n,i,j,k)-2.0*urijr(n,i,j-1,k)
          endif
          
          tmp3k=6.0*(urikl(n,i,j,k)-urikr(n,i,j,k-1)) &
               *(uri(n,i,j,k)-0.5*(urikl(n,i,j,k)+urikr(n,i,j,k-1)))
          tmp4k=(urikl(n,i,j,k)-urikr(n,i,j,k-1)) &
               *(urikl(n,i,j,k)-urikr(n,i,j,k-1))
          
!          tmp3k=(urikl(n,i,j,k)-urikr(n,i,j,k-1)) &
!               *(urikr(n,i,j,k-1)-urikr1(n,i,j,k-1))
!          tmp4k=(urikl(n,i,j,k)-urikr(n,i,j,k-1)) &
!               *(urikl1(n,i,j,k)-urikl(n,i,j,k))
          
          if(tmp3k .gt. tmp4k) then
            urikr(n,i,j,k-1)=3.0*uri(n,i,j,k)-2.0*urikl(n,i,j,k)
          endif
          if(-tmp3k .gt. tmp4k) then
            urikl(n,i,j,k)=3.0*uri(n,i,j,k)-2.0*urikr(n,i,j,k-1)
          endif

        enddo
      enddo
    enddo
  enddo
!
  deallocate( uriir1, uriil1, urijr1, urijl1, urikr1, urikl1, stat=merr )
!  
  return
end subroutine monoto
!
!---------------------------------------------------------------------@
subroutine mp5(uri,uriir,urijr,urikr,uriil,urijl,urikl, &
               is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Monotonicity Preserving Scheme for Reconstruction 
!     Calculate cell-interface variables(uriil, uriir) 
!     from cell-center variables (uri)
!
  use pram, only : imax, jmax, kmax, nv
  implicit none
!
  integer :: i, j, k, n, is1, ie1, js1, je1, ks1, ke1

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)
  
  real(8) :: tmp1, tmp2, alpha, eps1
  real(8) :: uriil1, uriilmp, tmpil, dim1, di1, dip1, dm4iph, dm4imh, &
             uriilul, uriilav, uriilmd, uriillc, &
             tmpilm1, tmpilm2, tmpilm3, tmpilm4, uriilmax, uriilmin, &
             uriir1, uriirmp, tmpir, &
             uriirul, uriirav, uriirmd, uriirlc, &
             tmpirm1, tmpirm2, tmpirm3, tmpirm4, uriirmax, uriirmin, &
             urijl1, urijlmp, tmpjl, djm1, dj1, djp1, dm4jph, dm4jmh, &
             urijlul, urijlav, urijlmd, urijllc, &
             tmpjlm1, tmpjlm2, tmpjlm3, tmpjlm4, urijlmax, urijlmin, &
             urijr1, urijrmp, tmpjr, &
             urijrul, urijrav, urijrmd, urijrlc, &
             tmpjrm1, tmpjrm2, tmpjrm3, tmpjrm4, urijrmax, urijrmin, &
             urikl1, uriklmp, tmpkl, dkm1, dk1, dkp1, dm4kph, dm4kmh, &
             uriklul, uriklav, uriklmd, urikllc, &
             tmpklm1, tmpklm2, tmpklm3, tmpklm4, uriklmax, uriklmin, &
             urikr1, urikrmp, tmpkr, &
             urikrul, urikrav, urikrmd, urikrlc, &
             tmpkrm1, tmpkrm2, tmpkrm3, tmpkrm4, urikrmax, urikrmin
  real(8) :: dmm, dm4
!
!-----------------------------------------------------------------------
!    Parameter
!

  tmp1=1.d0/60.d0
  tmp2=4.d0/3.d0
  alpha=4.d0
  eps1=1.0d-10
!
  do k=ks1+2,ke1-2
    do j=js1+2,je1-2
      do i=is1+2,ie1-2
        do n=1,nv
!!
! uriil
!
          uriil1=tmp1*(2.*uri(n,i-2,j,k)-13.*uri(n,i-1,j,k) &
                 +47.*uri(n,i,j,k)+27.*uri(n,i+1,j,k)-3.*uri(n,i+2,j,k))
          uriilmp=uri(n,i,j,k)+dmm(uri(n,i+1,j,k)-uri(n,i,j,k), &
                  alpha*(uri(n,i,j,k)-uri(n,i-1,j,k)))
 
          tmpil=(uriil1-uri(n,i,j,k))*(uriil1-uriilmp)

          if(tmpil .le. eps1) then
            uriil(n,i,j,k)=uriil1
          else
            dim1=uri(n,i-2,j,k)-2.*uri(n,i-1,j,k)+uri(n,i,j,k)
            di1=uri(n,i-1,j,k)-2.*uri(n,i,j,k)+uri(n,i+1,j,k)
            dip1=uri(n,i,j,k)-2.*uri(n,i+1,j,k)+uri(n,i+2,j,k)
           
            dm4iph=dm4(4.*di1-dip1,4.*dip1-di1,di1,dip1)
            dm4imh=dm4(4.*di1-dim1,4.*dim1-di1,di1,dim1)
           
            uriilul=uri(n,i,j,k)+alpha*(uri(n,i,j,k)-uri(n,i-1,j,k))
            uriilav=0.5*(uri(n,i,j,k)+uri(n,i+1,j,k))
            uriilmd=uriilav-0.5*dm4iph
            uriillc=uri(n,i,j,k)+0.5*(uri(n,i,j,k)-uri(n,i-1,j,k))+tmp2*dm4imh
           
            tmpilm1=min(uri(n,i,j,k), uri(n,i+1,j,k), uriilmd)
            tmpilm2=min(uri(n,i,j,k), uriilul, uriillc)
            tmpilm3=max(uri(n,i,j,k), uri(n,i+1,j,k), uriilmd)
            tmpilm4=max(uri(n,i,j,k), uriilul, uriillc)
          
            uriilmax=max(tmpilm1,tmpilm2)
            uriilmin=min(tmpilm3,tmpilm4)
          
            uriil(n,i,j,k)=uriil1+dmm(uriilmin-uriil1,uriilmax-uriil1)
          endif
!!
!  uriir
!          
          uriir1=tmp1*(2.*uri(n,i+2,j,k)-13.*uri(n,i+1,j,k) &
                 +47.*uri(n,i,j,k)+27.*uri(n,i-1,j,k)-3.*uri(n,i-2,j,k))
          uriirmp=uri(n,i,j,k)+dmm(uri(n,i-1,j,k)-uri(n,i,j,k), & 
                  alpha*(uri(n,i,j,k)-uri(n,i+1,j,k)))
     
          tmpir=(uriir1-uri(n,i,j,k))*(uriir1-uriirmp)
          
          if(tmpir .le. eps1) then
            uriir(n,i-1,j,k)=uriir1
          else 
            dim1=uri(n,i-2,j,k)-2.*uri(n,i-1,j,k)+uri(n,i,j,k)
            di1=uri(n,i-1,j,k)-2.*uri(n,i,j,k)+uri(n,i+1,j,k)
            dip1=uri(n,i,j,k)-2.*uri(n,i+1,j,k)+uri(n,i+2,j,k)
           
            dm4iph=dm4(4.*di1-dip1,4.*dip1-di1,di1,dip1)
            dm4imh=dm4(4.*di1-dim1,4.*dim1-di1,di1,dim1)
           
            uriirul=uri(n,i,j,k)+alpha*(uri(n,i,j,k)-uri(n,i+1,j,k))
            uriirav=0.5*(uri(n,i,j,k)+uri(n,i-1,j,k))
           
!            uriirmd=uriirav-0.5*dm4imh
!            uriirlc=uri(n,i,j,k)+0.5*(uri(n,i,j,k)-uri(n,i+1,j,k))+tmp2*dm4iph
     
            uriirmd=uriirav+0.5*dm4imh
            uriirlc=uri(n,i,j,k)+0.5*(uri(n,i,j,k)-uri(n,i+1,j,k))-tmp2*dm4iph
     
            tmpirm1=min(uri(n,i,j,k), uri(n,i-1,j,k), uriirmd)
            tmpirm2=min(uri(n,i,j,k), uriirul, uriirlc)
            tmpirm3=max(uri(n,i,j,k), uri(n,i-1,j,k), uriirmd)
            tmpirm4=max(uri(n,i,j,k), uriirul, uriirlc)
          
            uriirmax=max(tmpirm1,tmpirm2)
            uriirmin=min(tmpirm3,tmpirm4)
          
            uriir(n,i-1,j,k)=uriir1+dmm(uriirmin-uriir1,uriirmax-uriir1)
          endif
!!
!  urijl
!
          urijl1=tmp1*(2.*uri(n,i,j-2,k)-13.*uri(n,i,j-1,k) &
                 +47.*uri(n,i,j,k)+27.*uri(n,i,j+1,k)-3.*uri(n,i,j+2,k))
          urijlmp=uri(n,i,j,k)+dmm(uri(n,i,j+1,k)-uri(n,i,j,k), & 
                  alpha*(uri(n,i,j,k)-uri(n,i,j-1,k)))
          
          tmpjl=(urijl1-uri(n,i,j,k))*(urijl1-urijlmp)
!
          if(tmpjl .le. eps1) then
            urijl(n,i,j,k)=urijl1
          else
            djm1=uri(n,i,j-2,k)-2.*uri(n,i,j-1,k)+uri(n,i,j,k)
            dj1=uri(n,i,j-1,k)-2.*uri(n,i,j,k)+uri(n,i,j+1,k)
            djp1=uri(n,i,j,k)-2.*uri(n,i,j+1,k)+uri(n,i,j+2,k)
           
            dm4jph=dm4(4.*dj1-djp1,4.*djp1-dj1,dj1,djp1)
            dm4jmh=dm4(4.*dj1-djm1,4.*djm1-dj1,dj1,djm1)
            
            urijlul=uri(n,i,j,k)+alpha*(uri(n,i,j,k)-uri(n,i,j-1,k))
            urijlav=0.5*(uri(n,i,j,k)+uri(n,i,j+1,k))
            urijlmd=urijlav-0.5*dm4jph
            urijllc=uri(n,i,j,k)+0.5*(uri(n,i,j,k)-uri(n,i,j-1,k))+tmp2*dm4jmh
           
            tmpjlm1=min(uri(n,i,j,k), uri(n,i,j+1,k), urijlmd)
            tmpjlm2=min(uri(n,i,j,k), urijlul, urijllc)
            tmpjlm3=max(uri(n,i,j,k), uri(n,i,j+1,k), urijlmd)
            tmpjlm4=max(uri(n,i,j,k), urijlul, urijllc)
          
            urijlmax=max(tmpjlm1,tmpjlm2)
            urijlmin=min(tmpjlm3,tmpjlm4)
          
            urijl(n,i,j,k)=urijl1+dmm(urijlmin-urijl1,urijlmax-urijl1)
          endif
!!
!  urijr
!
          urijr1=tmp1*(2.*uri(n,i,j+2,k)-13.*uri(n,i,j+1,k) &
                 +47.*uri(n,i,j,k)+27.*uri(n,i,j-1,k)-3.*uri(n,i,j-2,k))
          urijrmp=uri(n,i,j,k)+dmm(uri(n,i,j-1,k)-uri(n,i,j,k), &
                  alpha*(uri(n,i,j,k)-uri(n,i,j+1,k)))
          
          tmpjr=(urijr1-uri(n,i,j,k))*(urijr1-urijrmp)
          
          if(tmpjr .le. eps1) then
            urijr(n,i,j-1,k)=urijr1
          else 
            djm1=uri(n,i,j-2,k)-2.*uri(n,i,j-1,k)+uri(n,i,j,k)
            dj1=uri(n,i,j-1,k)-2.*uri(n,i,j,k)+uri(n,i,j+1,k)
            djp1=uri(n,i,j,k)-2.*uri(n,i,j+1,k)+uri(n,i,j+2,k)
            
            dm4jph=dm4(4.*dj1-djp1,4.*djp1-dj1,dj1,djp1)
            dm4jmh=dm4(4.*dj1-djm1,4.*djm1-dj1,dj1,djm1)
           
            urijrul=uri(n,i,j,k)+alpha*(uri(n,i,j,k)-uri(n,i,j+1,k))
            urijrav=0.5*(uri(n,i,j,k)+uri(n,i,j-1,k))

!           urijrmd=urijrav-0.5*dm4jmh
!           urijrlc=uri(n,i,j,k)+0.5*(uri(n,i,j,k)-uri(n,i,j+1,k))+tmp2*dm4jph
     
            urijrmd=urijrav+0.5*dm4jmh
            urijrlc=uri(n,i,j,k)+0.5*(uri(n,i,j,k)-uri(n,i,j+1,k))-tmp2*dm4jph
     
            tmpjrm1=min(uri(n,i,j,k), uri(n,i,j-1,k), urijrmd)
            tmpjrm2=min(uri(n,i,j,k), urijrul, urijrlc)
            tmpjrm3=max(uri(n,i,j,k), uri(n,i,j-1,k), urijrmd)
            tmpjrm4=max(uri(n,i,j,k), urijrul, urijrlc)
          
            urijrmax=max(tmpjrm1,tmpjrm2)
            urijrmin=min(tmpjrm3,tmpjrm4)
          
            urijr(n,i,j-1,k)=urijr1+dmm(urijrmin-urijr1,urijrmax-urijr1)
          endif
!!
!  urikl     
!     
          urikl1=tmp1*(2.*uri(n,i,j,k-2)-13.*uri(n,i,j,k-1) &
                 +47.*uri(n,i,j,k)+27.*uri(n,i,j,k+1)-3.*uri(n,i,j,k+2))
          uriklmp=uri(n,i,j,k)+dmm(uri(n,i,j,k+1)-uri(n,i,j,k), & 
                  alpha*(uri(n,i,j,k)-uri(n,i,j,k-1)))
          
          tmpkl=(urikl1-uri(n,i,j,k))*(urikl1-uriklmp)

          if(tmpkl .le. eps1) then
            urikl(n,i,j,k)=urikl1
          else
            dkm1=uri(n,i,j,k-2)-2.*uri(n,i,j,k-1)+uri(n,i,j,k)
            dk1=uri(n,i,j,k-1)-2.*uri(n,i,j,k)+uri(n,i,j,k+1)
            dkp1=uri(n,i,j,k)-2.*uri(n,i,j,k+1)+uri(n,i,j,k+2)
           
            dm4kph=dm4(4.*dk1-dkp1,4.*dkp1-dk1,dk1,dkp1)
            dm4kmh=dm4(4.*dk1-dkm1,4.*dkm1-dk1,dk1,dkm1)
           
            uriklul=uri(n,i,j,k)+alpha*(uri(n,i,j,k)-uri(n,i,j,k-1))
            uriklav=0.5*(uri(n,i,j,k)+uri(n,i,j,k+1))
            uriklmd=uriklav-0.5*dm4kph
            urikllc=uri(n,i,j,k)+0.5*(uri(n,i,j,k)-uri(n,i,j,k-1))+tmp2*dm4kmh
           
            tmpklm1=min(uri(n,i,j,k), uri(n,i,j,k+1), uriklmd)
            tmpklm2=min(uri(n,i,j,k), uriklul, urikllc)
            tmpklm3=max(uri(n,i,j,k), uri(n,i,j,k+1), uriklmd)
            tmpklm4=max(uri(n,i,j,k), uriklul, urikllc)
          
            uriklmax=max(tmpklm1,tmpklm2)
            uriklmin=min(tmpklm3,tmpklm4)
          
            urikl(n,i,j,k)=urikl1+dmm(uriklmin-urikl1,uriklmax-urikl1)
          endif
!!
!  urikr
!
          urikr1=tmp1*(2.*uri(n,i,j,k+2)-13.*uri(n,i,j,k+1) &
                 +47.*uri(n,i,j,k)+27.*uri(n,i,j,k-1)-3.*uri(n,i,j,k-2))
          urikrmp=uri(n,i,j,k)+dmm(uri(n,i,j,k-1)-uri(n,i,j,k), &
                  alpha*(uri(n,i,j,k)-uri(n,i,j,k+1)))
          
          tmpkr=(urikr1-uri(n,i,j,k))*(urikr1-urikrmp)
          
          if(tmpkr .le. eps1) then
            urikr(n,i,j,k-1)=urikr1
          else 
            dkm1=uri(n,i,j,k-2)-2.*uri(n,i,j,k-1)+uri(n,i,j,k)
            dk1=uri(n,i,j,k-1)-2.*uri(n,i,j,k)+uri(n,i,j,k+1)
            dkp1=uri(n,i,j,k)-2.*uri(n,i,j,k+1)+uri(n,i,j,k+2)
            
            dm4kph=dm4(4.*dk1-dkp1,4.*dkp1-dk1,dk1,dkp1)
            dm4kmh=dm4(4.*dk1-dkm1,4.*dkm1-dk1,dk1,dkm1)
           
            urikrul=uri(n,i,j,k)+alpha*(uri(n,i,j,k)-uri(n,i,j,k+1))
            urikrav=0.5*(uri(n,i,j,k)+uri(n,i,j,k-1))
           
!            urikrmd=urikrav-0.5*dm4kmh
!            urikrlc=uri(n,i,j,k)+0.5*(uri(n,i,j,k)-uri(n,i,j,k+1))+tmp2*dm4kph
     
            urikrmd=urikrav+0.5*dm4kmh
            urikrlc=uri(n,i,j,k)+0.5*(uri(n,i,j,k)-uri(n,i,j,k+1))-tmp2*dm4kph
     
            tmpkrm1=min(uri(n,i,j,k), uri(n,i,j,k-1), urikrmd)
            tmpkrm2=min(uri(n,i,j,k), urikrul, urikrlc)
            tmpkrm3=max(uri(n,i,j,k), uri(n,i,j,k-1), urikrmd)
            tmpkrm4=max(uri(n,i,j,k), urikrul, urikrlc)
          
            urikrmax=max(tmpkrm1,tmpkrm2)
            urikrmin=min(tmpkrm3,tmpkrm4)
           
            urikr(n,i,j,k-1)=urikr1+dmm(urikrmin-urikr1,urikrmax-urikr1)
          endif

        enddo 
      enddo
    enddo
  enddo

  return
end subroutine mp5
      
!--------------------------------------------------------------------
function dmm(xx,yy)
!--------------------------------------------------------------------
 
  implicit none
  real(8) :: dmm, xx, yy
!
  dmm=0.5*(sign(1.0d0,xx)+sign(1.0d0,yy))*min(abs(xx),abs(yy))
!
  return
end function dmm
!
!--------------------------------------------------------------------
function dm4(ww,xx,yy,zz)
!--------------------------------------------------------------------

  implicit none
  real(8) :: dm4, ww, xx, yy, zz
!
  dm4=0.125*(sign(1.0d0,ww) + sign(1.0d0,xx))* &
         abs(sign(1.0d0,ww) + sign(1.0d0,yy))* &
            (sign(1.0d0,ww) + sign(1.0d0,zz))* &
         min(abs(ww),abs(xx),abs(yy),abs(zz))
!
  return
end function dm4
!
!---------------------------------------------------------------------@
subroutine weno5(uri,uriir,urijr,urikr,uriil,urijl,urikl, &
                 is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Weighted Essentially Non-Oscillatory Reconstruction scheme
!     Calculate cell-interface variables(urir, uril) 
!     from cell-center variables (uu, uri)
!
  use pram, only : imax, jmax, kmax, nv
  implicit none
!
  integer :: i, j, k, n, is1, ie1, js1, je1, ks1, ke1

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)
!
! Jiang & Shu (1990)
!
  real(8), parameter :: cd0=1.d0/10.d0, cd1=6.d0/10.d0, cd2=3.d0/10.d0
  real(8), parameter :: c00=2.d0/6.d0, c01=-7.d0/6.d0, c02=11.d0/6.d0, &
                        c10=-1.d0/6.d0, c11=5.d0/6.d0,  c12=2.d0/6.d0, &
                        c20=2.d0/6.d0, c21=5.d0/6.d0,  c22=-1.d0/6.d0
!
! Tchekhovskoy et al. (2008)
!
!  real(8), parameter :: cd0=1.d0/16.d0, cd1=10.d0/16.d0, cd2=5.d0/16.d0
!  real(8), parameter :: c00=3.d0/8.d0, c01=-10.d0/8.d0, c02=15.d0/8.d0, &
!                        c10=-1.d0/8.d0, c11=6.d0/8.d0,  c12=3.d0/8.d0, &
!                        c20=3.d0/8.d0, c21=6.d0/8.d0,  c22=-1.d0/8.d0
!
  real(8), parameter :: epsi=1.d-6
!
  real(8) :: betail0, betail1, betail2, omegil10, omegil11, omegil12, & 
             omegitl, omegil0, omegil1, omegil2, uriil10, uriil11, uriil12, &
             betair0, betair1, betair2, omegir10, omegir11, omegir12, &
             omegitr, omegir0, omegir1, omegir2, uriir10, uriir11, uriir12, &
             betajl0, betajl1, betajl2, omegjl10, omegjl11, omegjl12, &
             omegjtl, omegjl0, omegjl1, omegjl2, urijl10, urijl11, urijl12, &
             betajr0, betajr1, betajr2, omegjr10, omegjr11, omegjr12, &
             omegjtr, omegjr0, omegjr1, omegjr2, urijr10, urijr11, urijr12, &
             betakl0, betakl1, betakl2, omegkl10, omegkl11, omegkl12, &
             omegktl, omegkl0, omegkl1, omegkl2, urikl10, urikl11, urikl12, &
             betakr0, betakr1, betakr2, omegkr10, omegkr11, omegkr12, &
             omegktr, omegkr0, omegkr1, omegkr2, urikr10, urikr11, urikr12
!
!-----------------------------------------------------------------------
!
  do k=ks1+2,ke1-2
    do j=js1+2,je1-2
      do i=is1+2,ie1-2
        do n=1,nv
!!
!! x-direction
!!
!
! smooth indicater

          betail0=(13./12.) &
                 *(uri(n,i-2,j,k)-2.*uri(n,i-1,j,k)+uri(n,i,j,k))**2 &
                 +(1./4.) &
                 *(uri(n,i-2,j,k)-4.*uri(n,i-1,j,k)+3.*uri(n,i,j,k))**2+epsi
          betail1=(13./12.) &
                 *(uri(n,i-1,j,k)-2.*uri(n,i,j,k)+uri(n,i+1,j,k))**2 &
                 +(1./4.) &
                 *(uri(n,i-1,j,k)-uri(n,i+1,j,k))**2+epsi
          betail2=(13./12.) &
                 *(uri(n,i,j,k)-2.*uri(n,i+1,j,k)+uri(n,i+2,j,k))**2 &
                 +(1./4.) &
                 *(3.*uri(n,i,j,k)-4.*uri(n,i+1,j,k)+uri(n,i+2,j,k))**2+epsi
          
          betair0=(13./12.) &
                 *(uri(n,i+2,j,k)-2.*uri(n,i+1,j,k)+uri(n,i,j,k))**2 &
                 +(1./4.) &
                 *(uri(n,i+2,j,k)-4.*uri(n,i+1,j,k)+3.*uri(n,i,j,k))**2+epsi
          betair1=(13./12.) &
                 *(uri(n,i+1,j,k)-2.*uri(n,i,j,k)+uri(n,i-1,j,k))**2 &
                 +(1./4.) &
                 *(uri(n,i+1,j,k)-uri(n,i-1,j,k))**2+epsi
          betair2=(13./12.) &
                 *(uri(n,i,j,k)-2.*uri(n,i-1,j,k)+uri(n,i-2,j,k))**2 &
                 +(1./4.) &
                 *(3.*uri(n,i,j,k)-4.*uri(n,i-1,j,k)+uri(n,i-2,j,k))**2+epsi
!     
! unnormalized weights     
!     
          omegil10=cd0/betail0**2
          omegil11=cd1/betail1**2
          omegil12=cd2/betail2**2
          omegitl=omegil10+omegil11+omegil12

          omegir10=cd0/betair0**2
          omegir11=cd1/betair1**2
          omegir12=cd2/betair2**2
          omegitr=omegir10+omegir11+omegir12
!
! normalized weights
!
          omegil0=omegil10/omegitl
          omegil1=omegil11/omegitl
          omegil2=omegil12/omegitl

          omegir0=omegir10/omegitr
          omegir1=omegir11/omegitr
          omegir2=omegir12/omegitr
!
! reconstruction function
!
          uriil10=c00*uri(n,i-2,j,k)+c01*uri(n,i-1,j,k)+c02*uri(n,i,j,k)
          uriil11=c10*uri(n,i-1,j,k)+c11*uri(n,i,j,k)+c12*uri(n,i+1,j,k)
          uriil12=c20*uri(n,i,j,k)+c21*uri(n,i+1,j,k)+c22*uri(n,i+2,j,k)
          
          uriir10=c00*uri(n,i+2,j,k)+c01*uri(n,i+1,j,k)+c02*uri(n,i,j,k)
          uriir11=c10*uri(n,i+1,j,k)+c11*uri(n,i,j,k)+c12*uri(n,i-1,j,k)
          uriir12=c20*uri(n,i,j,k)+c21*uri(n,i-1,j,k)+c22*uri(n,i-2,j,k)
!
! cell interface values
!
          uriil(n,i,j,k)=omegil0*uriil10+omegil1*uriil11 +omegil2*uriil12
          uriir(n,i-1,j,k)=omegir0*uriir10+omegir1*uriir11+omegir2*uriir12


!!
!! y-direction
!!
!
! smooth indicater
!
          betajl0=(13./12.) &
                 *(uri(n,i,j-2,k)-2.0*uri(n,i,j-1,k)+uri(n,i,j,k))**2 &
                 +(1./4.) &
                 *(uri(n,i,j-2,k)-4.*uri(n,i,j-1,k)+3.*uri(n,i,j,k))**2+epsi
          betajl1=(13./12.) &
                 *(uri(n,i,j-1,k)-2.*uri(n,i,j,k)+uri(n,i,j+1,k))**2 &
                 +(1./4.) &
                 *(uri(n,i,j-1,k)-uri(n,i,j+1,k))**2+epsi
          betajl2=(13./12.) &
                 *(uri(n,i,j,k)-2.*uri(n,i,j+1,k)+uri(n,i,j+2,k))**2 &
                 +(1./4.) &
                 *(3.*uri(n,i,j,k)-4.*uri(n,i,j+1,k)+uri(n,i,j+2,k))**2+epsi
          
          betajr0=(13./12.) &
                 *(uri(n,i,j+2,k)-2.*uri(n,i,j+1,k)+uri(n,i,j,k))**2 &
                 +(1./4.) &
                 *(uri(n,i,j+2,k)-4.*uri(n,i,j+1,k)+3.*uri(n,i,j,k))**2+epsi
          betajr1=(13./12.) &
                 *(uri(n,i,j+1,k)-2.*uri(n,i,j,k)+uri(n,i,j-1,k))**2 &
                 +(1./4.) &
                 *(uri(n,i,j+1,k)-uri(n,i,j-1,k))**2+epsi
          betajr2=(13./12.) &
                 *(uri(n,i,j,k)-2.0*uri(n,i,j-1,k)+uri(n,i,j-2,k))**2 &
                 +(1./4.) &
                 *(3.*uri(n,i,j,k)-4.*uri(n,i,j-1,k)+uri(n,i,j-2,k))**2+epsi
!     
! unnormalized weights     
!     
          omegjl10=cd0/betajl0**2
          omegjl11=cd1/betajl1**2
          omegjl12=cd2/betajl2**2
          omegjtl=omegjl10+omegjl11+omegjl12

          omegjr10=cd0/betajr0**2
          omegjr11=cd1/betajr1**2
          omegjr12=cd2/betajr2**2
          omegjtr=omegjr10+omegjr11+omegjr12
!
! normalized weights
!
          omegjl0=omegjl10/omegjtl
          omegjl1=omegjl11/omegjtl
          omegjl2=omegjl12/omegjtl

          omegjr0=omegjr10/omegjtr
          omegjr1=omegjr11/omegjtr
          omegjr2=omegjr12/omegjtr
!
! reconstruction function
!
          urijl10=c00*uri(n,i,j-2,k)+c01*uri(n,i,j-1,k)+c02*uri(n,i,j,k)
          urijl11=c10*uri(n,i,j-1,k)+c11*uri(n,i,j,k)+c12*uri(n,i,j+1,k)
          urijl12=c20*uri(n,i,j,k)+c21*uri(n,i,j+1,k)+c22*uri(n,i,j+2,k)
          
          urijr10=c00*uri(n,i,j+2,k)+c01*uri(n,i,j+1,k)+c02*uri(n,i,j,k)
          urijr11=c10*uri(n,i,j+1,k)+c11*uri(n,i,j,k)+c12*uri(n,i,j-1,k)
          urijr12=c20*uri(n,i,j,k)+c21*uri(n,i,j-1,k)+c22*uri(n,i,j-2,k)
!
! cell interface values
!
          urijl(n,i,j,k)=omegjl0*urijl10+omegjl1*urijl11+omegjl2*urijl12
          urijr(n,i,j-1,k)=omegjr0*urijr10+omegjr1*urijr11+omegjr2*urijr12
!!
!! z-direction
!!
!
! smooth indicater
!
          betakl0=(13./12.) &
                 *(uri(n,i,j,k-2)-2.0*uri(n,i,j,k-1)+uri(n,i,j,k))**2 &
                 +(1./4.) &
                 *(uri(n,i,j,k-2)-4.*uri(n,i,j,k-1)+3.*uri(n,i,j,k))**2+epsi
          betakl1=(13./12.) &
                 *(uri(n,i,j,k-1)-2.*uri(n,i,j,k)+uri(n,i,j,k+1))**2 &
                 +(1./4.) &
                 *(uri(n,i,j,k-1)-uri(n,i,j,k+1))**2+epsi
          betakl2=(13./12.) &
                 *(uri(n,i,j,k)-2.*uri(n,i,j,k+1)+uri(n,i,j,k+2))**2 &
                 +(1./4.) &
                 *(3.*uri(n,i,j,k)-4.*uri(n,i,j,k+1)+uri(n,i,j,k+2))**2+epsi
          
          betakr0=(13./12.) &
                 *(uri(n,i,j,k+2)-2.*uri(n,i,j,k+1)+uri(n,i,j,k))**2 &
                 +(1./4.) &
                 *(uri(n,i,j,k+2)-4.*uri(n,i,j,k+1)+3.*uri(n,i,j,k))**2+epsi
          betakr1=(13./12.) &
                 *(uri(n,i,j,k+1)-2.*uri(n,i,j,k)+uri(n,i,j,k-1))**2 &
                 +(1./4.) &
                 *(uri(n,i,j,k+1)-uri(n,i,j,k-1))**2+epsi
          betakr2=(13./12.) &
                 *(uri(n,i,j,k)-2.0*uri(n,i,j,k-1)+uri(n,i,j,k-2))**2 &
                 +(1./4.) &
                 *(3.*uri(n,i,j,k)-4.*uri(n,i,j,k-1)+uri(n,i,j,k-2))**2+epsi
!     
! unnormalized weights     
!     
          omegkl10=cd0/betakl0**2
          omegkl11=cd1/betakl1**2
          omegkl12=cd2/betakl2**2
          omegktl=omegkl10+omegkl11+omegkl12

          omegkr10=cd0/betakr0**2
          omegkr11=cd1/betakr1**2
          omegkr12=cd2/betakr2**2
          omegktr=omegkr10+omegkr11+omegkr12
!
! normalized weights
!
          omegkl0=omegkl10/omegktl
          omegkl1=omegkl11/omegktl
          omegkl2=omegkl12/omegktl

          omegkr0=omegkr10/omegktr
          omegkr1=omegkr11/omegktr
          omegkr2=omegkr12/omegktr
!
! reconstruction function
!
          urikl10=c00*uri(n,i,j,k-2)+c01*uri(n,i,j,k-1)+c02*uri(n,i,j,k)
          urikl11=c10*uri(n,i,j,k-1)+c11*uri(n,i,j,k)+c12*uri(n,i,j,k+1)
          urikl12=c20*uri(n,i,j,k)+c21*uri(n,i,j,k+1)+c22*uri(n,i,j,k+2)
          
          urikr10=c00*uri(n,i,j,k+2)+c01*uri(n,i,j,k+1)+c02*uri(n,i,j,k)
          urikr11=c10*uri(n,i,j,k+1)+c11*uri(n,i,j,k)+c12*uri(n,i,j,k-1)
          urikr12=c20*uri(n,i,j,k)+c21*uri(n,i,j,k-1)+c22*uri(n,i,j,k-2)
!
! cell interface values
!
          urikl(n,i,j,k)=omegkl0*urikl10+omegkl1*urikl11+omegkl2*urikl12
          urikr(n,i,j,k-1)=omegkr0*urikr10+omegkr1*urikr11+omegkr2*urikr12
    
        enddo
      enddo
    enddo
  enddo

  return
end subroutine weno5
!
!---------------------------------------------------------------------@
subroutine mpweno5(uri,uriir,urijr,urikr,uriil,urijl,urikl, &
                   is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Monotonicity-Preserving
!     Weighted Essentially Non-Oscillatory Reconstruction scheme
!     Balsara & Shu (2000)
!     Calculate cell-interface variables(urir, uril) 
!     from cell-center variables (uu, uri)
!
  use pram, only : imax, jmax, kmax, nv
  implicit none
!
  integer :: i, j, k, n, is1, ie1, js1, je1, ks1, ke1

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)
!
! Jiang & Shu (1990)
!
  real(8), parameter :: cd0=1.d0/10.d0, cd1=6.d0/10.d0, cd2=3.d0/10.d0
  real(8), parameter :: c00=2.d0/6.d0, c01=-7.d0/6.d0, c02=11.d0/6.d0, &
                        c10=-1.d0/6.d0, c11=5.d0/6.d0,  c12=2.d0/6.d0, &
                        c20=2.d0/6.d0, c21=5.d0/6.d0,  c22=-1.d0/6.d0
!
! Tchekhovskoy et al. (2008)
!
!  real(8), parameter :: cd0=1.d0/16.d0, cd1=10.d0/16.d0, cd2=5.d0/16.d0
!  real(8), parameter :: c00=3.d0/8.d0, c01=-10.d0/8.d0, c02=15.d0/8.d0, &
!                        c10=-1.d0/8.d0, c11=6.d0/8.d0,  c12=3.d0/8.d0, &
!                        c20=3.d0/8.d0, c21=6.d0/8.d0,  c22=-1.d0/8.d0
!
  real(8) :: tmp1, tmp2, alpha, eps1, epsi
!
  real(8) :: betail0, betail1, betail2, omegil10, omegil11, omegil12, &
             omegitl, omegil0, omegil1, omegil2, uriil10, uriil11, uriil12, &
             uriil1, uriilmp, tmpil, dim1, di1, dip1, dm4iph, dm4imh, &
             uriilul, uriilav, uriilmd, uriillc, &
             tmpilm1, tmpilm2, tmpilm3, tmpilm4, uriilmax, uriilmin, &
             betair0, betair1, betair2, omegir10, omegir11, omegir12, &
             omegitr, omegir0, omegir1, omegir2, uriir10, uriir11, uriir12, &
             uriir1, uriirmp, tmpir, uriirul, uriirav, uriirmd, uriirlc, &
             tmpirm1, tmpirm2, tmpirm3, tmpirm4, uriirmax, uriirmin, &
             betajl0, betajl1, betajl2, omegjl10, omegjl11, omegjl12, &
             omegjtl, omegjl0, omegjl1, omegjl2, urijl10, urijl11, urijl12, &
             urijl1, urijlmp, tmpjl, djm1, dj1, djp1, dm4jph, dm4jmh, &
             urijlul, urijlav, urijlmd, urijllc, &
             tmpjlm1, tmpjlm2, tmpjlm3, tmpjlm4, urijlmax, urijlmin, &
             betajr0, betajr1, betajr2, omegjr10, omegjr11, omegjr12, & 
             omegjtr, omegjr0, omegjr1, omegjr2, urijr10, urijr11, urijr12, &
             urijr1, urijrmp, tmpjr, urijrul, urijrav, urijrmd, urijrlc, &
             tmpjrm1, tmpjrm2, tmpjrm3, tmpjrm4, urijrmax, urijrmin, &
             betakl0, betakl1, betakl2, omegkl10, omegkl11, omegkl12, &
             omegktl, omegkl0, omegkl1, omegkl2, urikl10, urikl11, urikl12, &
             urikl1, uriklmp, tmpkl, dkm1, dk1, dkp1, dm4kph, dm4kmh, &
             uriklul, uriklav, uriklmd, urikllc, &
             tmpklm1, tmpklm2, tmpklm3, tmpklm4, uriklmax, uriklmin, &
             betakr0, betakr1, betakr2, omegkr10, omegkr11, omegkr12, &
             omegktr, omegkr0, omegkr1, omegkr2, urikr10, urikr11, urikr12, &
             urikr1, urikrmp, tmpkr, urikrul, urikrav, urikrmd, urikrlc, &
             tmpkrm1, tmpkrm2, tmpkrm3, tmpkrm4, urikrmax, urikrmin
  real(8) :: dmm, dm4
!
!-----------------------------------------------------------------------
!    Parameter
!

  tmp1=1.d0/60.d0
  tmp2=4.d0/3.d0
  alpha=4.d0
  eps1=1.0d-10
  epsi=1.0d-6

!-----------------------------------------------------------------------
!
  do k=ks1+2,ke1-2
    do j=js1+2,je1-2
      do i=is1+2,ie1-2
        do n=1,nv
!!
!! x-direction
!!
!
! smooth indicater
!
          betail0=(13./12.) &
             *(uri(n,i-2,j,k)-2.*uri(n,i-1,j,k)+uri(n,i,j,k))**2 &
             +(1./4.) &
             *(uri(n,i-2,j,k)-4.*uri(n,i-1,j,k)+3.*uri(n,i,j,k))**2 &
             +epsi
          betail1=(13./12.) &
            *(uri(n,i-1,j,k)-2.*uri(n,i,j,k)+uri(n,i+1,j,k))**2 &
            +(1./4.) &
            *(uri(n,i-1,j,k)-uri(n,i+1,j,k))**2 &
            +epsi
          betail2=(13./12.) &
            *(uri(n,i,j,k)-2.*uri(n,i+1,j,k)+uri(n,i+2,j,k))**2 &
            +(1./4.) &
            *(3.*uri(n,i,j,k)-4.*uri(n,i+1,j,k)+uri(n,i+2,j,k))**2 &
            +epsi
                    
          betair0=(13./12.) &
            *(uri(n,i+2,j,k)-2.*uri(n,i+1,j,k)+uri(n,i,j,k))**2 &
            +(1./4.) &
            *(uri(n,i+2,j,k)-4.*uri(n,i+1,j,k)+3.*uri(n,i,j,k))**2 &
            +epsi
          betair1=(13./12.) &
            *(uri(n,i+1,j,k)-2.*uri(n,i,j,k)+uri(n,i-11,j,k))**2 &
            +(1./4.) &
            *(uri(n,i+1,j,k)-uri(n,i-1,j,k))**2 &
            +epsi
          betair2=(13./12.) &
            *(uri(n,i,j,k)-2.*uri(n,i-1,j,k)+uri(n,i-2,j,k))**2 &
            +(1./4.) &
            *(3.*uri(n,i,j,k)-4.*uri(n,i-1,j,k)+uri(n,i-2,j,k))**2 &
            +epsi
!     
! unnormalized weights     
!     
          omegil10=cd0/betail0**2
          omegil11=cd1/betail1**2
          omegil12=cd2/betail2**2
          omegitl=omegil10+omegil11+omegil12

          omegir10=cd0/betair0**2
          omegir11=cd1/betair1**2
          omegir12=cd2/betair2**2
          omegitr=omegir10+omegir11+omegir12
!
! normalized weights
!
          omegil0=omegil10/omegitl
          omegil1=omegil11/omegitl
          omegil2=omegil12/omegitl

          omegir0=omegir10/omegitr
          omegir1=omegir11/omegitr
          omegir2=omegir12/omegitr
!
! reconstruction function
!
          uriil10=c00*uri(n,i-2,j,k)+c01*uri(n,i-1,j,k)+c02*uri(n,i,j,k)
          uriil11=c10*uri(n,i-1,j,k)+c11*uri(n,i,j,k)+c12*uri(n,i+1,j,k)
          uriil12=c20*uri(n,i,j,k)+c21*uri(n,i+1,j,k)+c22*uri(n,i+2,j,k)
          
          uriir10=c00*uri(n,i+2,j,k)+c01*uri(n,i+1,j,k)+c02*uri(n,i,j,k)
          uriir11=c10*uri(n,i+1,j,k)+c11*uri(n,i,j,k)+c12*uri(n,i-1,j,k)
          uriir12=c20*uri(n,i,j,k)+c21*uri(n,i-1,j,k)+c22*uri(n,i-2,j,k)
!
! cell interface values
!
          uriil1=omegil0*uriil10+omegil1*uriil11 +omegil2*uriil12
          uriilmp=uri(n,i,j,k)+dmm(uri(n,i+1,j,k)-uri(n,i,j,k), &
                  alpha*(uri(n,i,j,k)-uri(n,i-1,j,k)))
          tmpil=(uriil1-uri(n,i,j,k))*(uriil1-uriilmp)
          
          uriir1=omegir0*uriir10+omegir1*uriir11+omegir2*uriir12
          uriirmp=uri(n,i,j,k)+dmm(uri(n,i-1,j,k)-uri(n,i,j,k), &
                  alpha*(uri(n,i,j,k)-uri(n,i+1,j,k)))
          tmpir=(uriir1-uri(n,i,j,k))*(uriir1-uriirmp)
!
! Monotonicity=preserving
!
          if(tmpil .le. eps1) then
            uriil(n,i,j,k)=uriil1
          else
            dim1=uri(n,i-2,j,k)-2.*uri(n,i-1,j,k)+uri(n,i,j,k)
            di1=uri(n,i-1,j,k)-2.*uri(n,i,j,k)+uri(n,i+1,j,k)
            dip1=uri(n,i,j,k)-2.*uri(n,i+1,j,k)+uri(n,i+2,j,k)
         
            dm4iph=dm4(4.*di1-dip1,4.*dip1-di1,di1,dip1)
            dm4imh=dm4(4.*di1-dim1,4.*dim1-di1,di1,dim1)
           
            uriilul=uri(n,i,j,k)+alpha*(uri(n,i,j,k)-uri(n,i-1,j,k))
            uriilav=0.5*(uri(n,i,j,k)+uri(n,i+1,j,k))
            uriilmd=uriilav-0.5*dm4iph
            uriillc=uri(n,i,j,k)+0.5*(uri(n,i,j,k)-uri(n,i-1,j,k)) &
                   +tmp2*dm4imh
           
            tmpilm1=min(uri(n,i,j,k), uri(n,i+1,j,k), uriilmd)
            tmpilm2=min(uri(n,i,j,k), uriilul, uriillc)
            tmpilm3=max(uri(n,i,j,k), uri(n,i+1,j,k), uriilmd)
            tmpilm4=max(uri(n,i,j,k), uriilul, uriillc)
          
            uriilmax=max(tmpilm1,tmpilm2)
            uriilmin=min(tmpilm3,tmpilm4)
          
            uriil(n,i,j,k)=uriil1+dmm(uriilmin-uriil1,uriilmax-uriil1)
          endif
!
          if(tmpir .le. eps1) then
            uriir(n,i-1,j,k)=uriir1
          else 
            dim1=uri(n,i-2,j,k)-2.*uri(n,i-1,j,k)+uri(n,i,j,k)
            di1=uri(n,i-1,j,k)-2.*uri(n,i,j,k)+uri(n,i+1,j,k)
            dip1=uri(n,i,j,k)-2.*uri(n,i+1,j,k)+uri(n,i+2,j,k)
         
            dm4iph=dm4(4.*di1-dip1,4.*dip1-di1,di1,dip1)
            dm4imh=dm4(4.*di1-dim1,4.*dim1-di1,di1,dim1)
           
            uriirul=uri(n,i,j,k)+alpha*(uri(n,i,j,k)-uri(n,i+1,j,k))
            uriirav=0.5*(uri(n,i,j,k)+uri(n,i-1,j,k))
           
!           uriirmd=uriirav-0.5*dm4imh
!           uriirlc=uri(n,i,j,k)+0.5*(uri(n,i,j,k)-uri(n,i+1,j,k))
!     &             +tmp2*dm4iph
     
            uriirmd=uriirav+0.5*dm4imh
            uriirlc=uri(n,i,j,k)+0.5*(uri(n,i,j,k)-uri(n,i+1,j,k)) &
                   -tmp2*dm4iph
     
            tmpirm1=min(uri(n,i,j,k), uri(n,i-1,j,k), uriirmd)
            tmpirm2=min(uri(n,i,j,k), uriirul, uriirlc)
            tmpirm3=max(uri(n,i,j,k), uri(n,i-1,j,k), uriirmd)
            tmpirm4=max(uri(n,i,j,k), uriirul, uriirlc)
          
            uriirmax=max(tmpirm1,tmpirm2)
            uriirmin=min(tmpirm3,tmpirm4)
          
            uriir(n,i-1,j,k)=uriir1+dmm(uriirmin-uriir1,uriirmax-uriir1)
          endif
!!
!! y-direction
!!
!
! smooth indicater
!
          betajl0=(13./12.) &
            *(uri(n,i,j-2,k)-2.0*uri(n,i,j-1,k)+uri(n,i,j,k))**2 &
            +(1./4.) &
            *(uri(n,i,j-2,k)-4.*uri(n,i,j-1,k)+3.*uri(n,i,j,k))**2 &
            +epsi
          betajl1=(13./12.) &
            *(uri(n,i,j-1,k)-2.*uri(n,i,j,k)+uri(n,i,j+1,k))**2 &
            +(1./4.) &
            *(uri(n,i,j-1,k)-uri(n,i,j+1,k))**2 &
            +epsi
          betajl2=(13./12.) &
            *(uri(n,i,j,k)-2.*uri(n,i,j+1,k)+uri(n,i,j+2,k))**2 &
            +(1./4.) &
            *(3.*uri(n,i,j,k)-4.*uri(n,i,j+1,k)+uri(n,i,j+2,k))**2 &
            +epsi
          
          betajr0=(13./12.) &
            *(uri(n,i,j+2,k)-2.*uri(n,i,j+1,k)+uri(n,i,j,k))**2 &
            +(1./4.) &
            *(uri(n,i,j+2,k)-4.*uri(n,i,j+1,k)+3.*uri(n,i,j,k))**2 &
            +epsi
          betajr1=(13./12.) &
            *(uri(n,i,j+1,k)-2.*uri(n,i,j,k)+uri(n,i,j-1,k))**2 &
            +(1./4.) &
            *(uri(n,i,j+1,k)-uri(n,i,j-1,k))**2 &
            +epsi
          betajr2=(13./12.) &
            *(uri(n,i,j,k)-2.0*uri(n,i,j-1,k)+uri(n,i,j-2,k))**2 &
            +(1./4.) &
            *(3.*uri(n,i,j,k)-4.*uri(n,i,j-1,k)+uri(n,i,j-2,k))**2 &
            +epsi
!     
! unnormalized weights     
!     
          omegjl10=cd0/betajl0**2
          omegjl11=cd1/betajl1**2
          omegjl12=cd2/betajl2**2
          omegjtl=omegjl10+omegjl11+omegjl12

          omegjr10=cd0/betajr0**2
          omegjr11=cd1/betajr1**2
          omegjr12=cd2/betajr2**2
          omegjtr=omegjr10+omegjr11+omegjr12
!
! normalized weights
!
          omegjl0=omegjl10/omegjtl
          omegjl1=omegjl11/omegjtl
          omegjl2=omegjl12/omegjtl

          omegjr0=omegjr10/omegjtr
          omegjr1=omegjr11/omegjtr
          omegjr2=omegjr12/omegjtr
!
! reconstruction function
!
          urijl10=c00*uri(n,i,j-2,k)+c01*uri(n,i,j-1,k)+c02*uri(n,i,j,k)
          urijl11=c10*uri(n,i,j-1,k)+c11*uri(n,i,j,k)+c12*uri(n,i,j+1,k)
          urijl12=c20*uri(n,i,j,k)+c21*uri(n,i,j+1,k)+c22*uri(n,i,j+2,k)
          
          urijr10=c00*uri(n,i,j+2,k)+c01*uri(n,i,j+1,k)+c02*uri(n,i,j,k)
          urijr11=c10*uri(n,i,j+1,k)+c11*uri(n,i,j,k)+c12*uri(n,i,j-1,k)
          urijr12=c20*uri(n,i,j,k)+c21*uri(n,i,j-1,k)+c22*uri(n,i,j-2,k)
!
! cell interface values
!
          urijl1=omegjl0*urijl10+omegjl1*urijl11 +omegjl2*urijl12
          urijlmp=uri(n,i,j,k)+dmm(uri(n,i,j+1,k)-uri(n,i,j,k), &
                  alpha*(uri(n,i,j,k)-uri(n,i,j-1,k)))
          tmpjl=(urijl1-uri(n,i,j,k))*(urijl1-urijlmp)
          
          urijr1=omegjr0*urijr10+omegjr1*urijr11+omegjr2*urijr12
          urijrmp=uri(n,i,j,k)+dmm(uri(n,i,j-1,k)-uri(n,i,j,k), & 
                  alpha*(uri(n,i,j,k)-uri(n,i+1,j,k)))
          tmpjr=(urijr1-uri(n,i,j,k))*(urijr1-urijrmp)
!
! Monotonicity-Preserving
!
          if(tmpjl .le. eps1) then
            urijl(n,i,j,k)=urijl1
          else
            djm1=uri(n,i,j-2,k)-2.*uri(n,i,j-1,k)+uri(n,i,j,k)
            dj1=uri(n,i,j-1,k)-2.*uri(n,i,j,k)+uri(n,i,j+1,k)
            djp1=uri(n,i,j,k)-2.*uri(n,i,j+1,k)+uri(n,i,j+2,k)
         
            dm4jph=dm4(4.*dj1-djp1,4.*djp1-dj1,dj1,djp1)
            dm4jmh=dm4(4.*dj1-djm1,4.*djm1-dj1,dj1,djm1)
            
            urijlul=uri(n,i,j,k)+alpha*(uri(n,i,j,k)-uri(n,i,j-1,k))
            urijlav=0.5*(uri(n,i,j,k)+uri(n,i,j+1,k))
            urijlmd=urijlav-0.5*dm4jph
            urijllc=uri(n,i,j,k)+0.5*(uri(n,i,j,k)-uri(n,i,j-1,k))+tmp2*dm4jmh
           
            tmpjlm1=min(uri(n,i,j,k), uri(n,i,j+1,k), urijlmd)
            tmpjlm2=min(uri(n,i,j,k), urijlul, urijllc)
            tmpjlm3=max(uri(n,i,j,k), uri(n,i,j+1,k), urijlmd)
            tmpjlm4=max(uri(n,i,j,k), urijlul, urijllc)
          
            urijlmax=max(tmpjlm1,tmpjlm2)
            urijlmin=min(tmpjlm3,tmpjlm4)
          
            urijl(n,i,j,k)=urijl1+dmm(urijlmin-urijl1,urijlmax-urijl1)
          endif
!
          if(tmpjr .le. eps1) then
            urijr(n,i,j-1,k)=urijr1
          else 
            djm1=uri(n,i,j-2,k)-2.*uri(n,i,j-1,k)+uri(n,i,j,k)
            dj1=uri(n,i,j-1,k)-2.*uri(n,i,j,k)+uri(n,i,j+1,k)
            djp1=uri(n,i,j,k)-2.*uri(n,i,j+1,k)+uri(n,i,j+2,k)
         
            dm4jph=dm4(4.*dj1-djp1,4.*djp1-dj1,dj1,djp1)
            dm4jmh=dm4(4.*dj1-djm1,4.*djm1-dj1,dj1,djm1)
           
            urijrul=uri(n,i,j,k)+alpha*(uri(n,i,j,k)-uri(n,i,j+1,k))
            urijrav=0.5*(uri(n,i,j,k)+uri(n,i,j-1,k))

!           urijrmd=urijrav-0.5*dm4jmh
!           urijrlc=uri(n,i,j,k)+0.5*(uri(n,i,j,k)-uri(n,i,j+1,k))
!     &             +tmp2*dm4jph
     
            urijrmd=urijrav+0.5*dm4jmh
            urijrlc=uri(n,i,j,k)+0.5*(uri(n,i,j,k)-uri(n,i,j+1,k))-tmp2*dm4jph
     
            tmpjrm1=min(uri(n,i,j,k), uri(n,i,j-1,k), urijrmd)
            tmpjrm2=min(uri(n,i,j,k), urijrul, urijrlc)
            tmpjrm3=max(uri(n,i,j,k), uri(n,i,j-1,k), urijrmd)
            tmpjrm4=max(uri(n,i,j,k), urijrul, urijrlc)
          
            urijrmax=max(tmpjrm1,tmpjrm2)
            urijrmin=min(tmpjrm3,tmpjrm4)
          
            urijr(n,i,j-1,k)=urijr1+dmm(urijrmin-urijr1,urijrmax-urijr1)
          endif
!!
!! z-direction
!!
!
! smooth indicater
!
          betakl0=(13./12.) &
            *(uri(n,i,j,k-2)-2.0*uri(n,i,j,k-1)+uri(n,i,j,k))**2 &
            +(1./4.) &
            *(uri(n,i,j,k-2)-4.*uri(n,i,j,k-1)+3.*uri(n,i,j,k))**2 &
            +epsi
          betakl1=(13./12.) &
            *(uri(n,i,j,k-1)-2.*uri(n,i,j,k)+uri(n,i,j,k+1))**2 &
            +(1./4.) &
            *(uri(n,i,j,k-1)-uri(n,i,j,k+1))**2 &
            +epsi
          betakl2=(13./12.) &
            *(uri(n,i,j,k)-2.*uri(n,i,j,k+1)+uri(n,i,j,k+2))**2 &
            +(1./4.) &
            *(3.*uri(n,i,j,k)-4.*uri(n,i,j,k+1)+uri(n,i,j,k+2))**2 &
            +epsi
          
          betakr0=(13./12.) &
            *(uri(n,i,j,k+2)-2.*uri(n,i,j,k+1)+uri(n,i,j,k))**2 &
            +(1./4.) &
            *(uri(n,i,j,k+2)-4.*uri(n,i,j,k+1)+3.*uri(n,i,j,k))**2 &
            +epsi
          betakr1=(13./12.) &
            *(uri(n,i,j,k+1)-2.*uri(n,i,j,k)+uri(n,i,j,k-1))**2 &
            +(1./4.) &
            *(uri(n,i,j,k+1)-uri(n,i,j,k-1))**2 &
            +epsi
          betakr2=(13./12.) &
            *(uri(n,i,j,k)-2.0*uri(n,i,j,k-1)+uri(n,i,j,k-2))**2 &
            +(1./4.) &
            *(3.*uri(n,i,j,k)-4.*uri(n,i,j,k-1)+uri(n,i,j,k-2))**2 &
            +epsi
!     
! unnormalized weights     
!     
          omegkl10=cd0/betakl0**2
          omegkl11=cd1/betakl1**2
          omegkl12=cd2/betakl2**2
          omegktl=omegkl10+omegkl11+omegkl12

          omegkr10=cd0/betakr0**2
          omegkr11=cd1/betakr1**2
          omegkr12=cd2/betakr2**2
          omegktr=omegkr10+omegkr11+omegkr12
!
! normalized weights
!
          omegkl0=omegkl10/omegktl
          omegkl1=omegkl11/omegktl
          omegkl2=omegkl12/omegktl

          omegkr0=omegkr10/omegktr
          omegkr1=omegkr11/omegktr
          omegkr2=omegkr12/omegktr
!
! reconstruction function
!
          urikl10=c00*uri(n,i,j,k-2)+c01*uri(n,i,j,k-1)+c02*uri(n,i,j,k)
          urikl11=c10*uri(n,i,j,k-1)+c11*uri(n,i,j,k)+c12*uri(n,i,j,k+1)
          urikl12=c20*uri(n,i,j,k)+c21*uri(n,i,j,k+1)+c22*uri(n,i,j,k+2)
          
          urikr10=c00*uri(n,i,j,k+2)+c01*uri(n,i,j,k+1)+c02*uri(n,i,j,k)
          urikr11=c10*uri(n,i,j,k+1)+c11*uri(n,i,j,k)+c12*uri(n,i,j,k-1)
          urikr12=c20*uri(n,i,j,k)+c21*uri(n,i,j,k-1)+c22*uri(n,i,j,k-2)
!
! cell interface values
!
          urikl1=omegkl0*urikl10+omegkl1*urikl11 +omegkl2*urikl12
          uriklmp=uri(n,i,j,k)+dmm(uri(n,i,j,k+1)-uri(n,i,j,k), &
                  alpha*(uri(n,i,j,k)-uri(n,i,j,k-1)))
          tmpkl=(urikl1-uri(n,i,j,k))*(urikl1-uriklmp)
          
          urikr1=omegkr0*urikr10+omegkr1*urikr11+omegkr2*urikr12
          urikrmp=uri(n,i,j,k)+dmm(uri(n,i,j,k-1)-uri(n,i,j,k), & 
                  alpha*(uri(n,i,j,k)-uri(n,i,j,k+1)))
          tmpkr=(urikr1-uri(n,i,j,k))*(urikr1-urikrmp)
!
! Monotonicity-Preserving
!
          if(tmpkl .le. eps1) then
            urikl(n,i,j,k)=urikl1
          else
            dkm1=uri(n,i,j,k-2)-2.*uri(n,i,j,k-1)+uri(n,i,j,k)
            dk1=uri(n,i,j,k-1)-2.*uri(n,i,j,k)+uri(n,i,j,k+1)
            dkp1=uri(n,i,j,k)-2.*uri(n,i,j,k+1)+uri(n,i,j,k+2)
         
            dm4kph=dm4(4.*dk1-dkp1,4.*dkp1-dk1,dk1,dkp1)
            dm4kmh=dm4(4.*dk1-dkm1,4.*dkm1-dk1,dk1,dkm1)
           
            uriklul=uri(n,i,j,k)+alpha*(uri(n,i,j,k)-uri(n,i,j,k-1))
            uriklav=0.5*(uri(n,i,j,k)+uri(n,i,j,k+1))
            uriklmd=uriklav-0.5*dm4kph
            urikllc=uri(n,i,j,k)+0.5*(uri(n,i,j,k)-uri(n,i,j,k-1))+tmp2*dm4kmh
           
            tmpklm1=min(uri(n,i,j,k), uri(n,i,j,k+1), uriklmd)
            tmpklm2=min(uri(n,i,j,k), uriklul, urikllc)
            tmpklm3=max(uri(n,i,j,k), uri(n,i,j,k+1), uriklmd)
            tmpklm4=max(uri(n,i,j,k), uriklul, urikllc)
          
            uriklmax=max(tmpklm1,tmpklm2)
            uriklmin=min(tmpklm3,tmpklm4)
          
            urikl(n,i,j,k)=urikl1+dmm(uriklmin-urikl1,uriklmax-urikl1)
          endif
!
          if(tmpkr .le. eps1) then
            urikr(n,i,j,k-1)=urikr1
          else 
            dkm1=uri(n,i,j,k-2)-2.*uri(n,i,j,k-1)+uri(n,i,j,k)
            dk1=uri(n,i,j,k-1)-2.*uri(n,i,j,k)+uri(n,i,j,k+1)
            dkp1=uri(n,i,j,k)-2.*uri(n,i,j,k+1)+uri(n,i,j,k+2)
         
            dm4kph=dm4(4.*dk1-dkp1,4.*dkp1-dk1,dk1,dkp1)
            dm4kmh=dm4(4.*dk1-dkm1,4.*dkm1-dk1,dk1,dkm1)
           
            urikrul=uri(n,i,j,k)+alpha*(uri(n,i,j,k)-uri(n,i,j,k+1))
            urikrav=0.5*(uri(n,i,j,k)+uri(n,i,j,k-1))
           
!           urikrmd=urikrav-0.5*dm4kmh
!           urikrlc=uri(n,i,j,k)+0.5*(uri(n,i,j,k)-uri(n,i,j,k+1))
!     &             +tmp2*dm4kph
     
            urikrmd=urikrav+0.5*dm4kmh
            urikrlc=uri(n,i,j,k)+0.5*(uri(n,i,j,k)-uri(n,i,j,k+1))-tmp2*dm4kph
     
            tmpkrm1=min(uri(n,i,j,k), uri(n,i,j,k-1), urikrmd)
            tmpkrm2=min(uri(n,i,j,k), urikrul, urikrlc)
            tmpkrm3=max(uri(n,i,j,k), uri(n,i,j,k-1), urikrmd)
            tmpkrm4=max(uri(n,i,j,k), urikrul, urikrlc)
          
            urikrmax=max(tmpkrm1,tmpkrm2)
            urikrmin=min(tmpkrm3,tmpkrm4)
          
            urikr(n,i,j,k-1)=urikr1+dmm(urikrmin-urikr1,urikrmax-urikr1)
          endif
          
        enddo
      enddo
    enddo
  enddo

  return
end subroutine mpweno5
!
!---------------------------------------------------------------------@
subroutine weno5z(uri,uriir,urijr,urikr,uriil,urijl,urikl, &
                  is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Weighted Essentially Non-Oscillatory Reconstruction scheme
!     5th order + modified weight (Borges et al. 2008)
!     Calculate cell-interface variables(urir, uril) 
!     from cell-center variables (uu, uri)
!
  use pram, only : imax, jmax, kmax, nv
  implicit none
!
  integer :: i, j, k, n, is1, ie1, js1, je1, ks1, ke1

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)
!
! Jiang & Shu (1990)
!
  real(8), parameter :: cd0=1.d0/10.d0, cd1=6.d0/10.d0, cd2=3.d0/10.d0
  real(8), parameter :: c00=2.d0/6.d0, c01=-7.d0/6.d0, c02=11.d0/6.d0, &
                        c10=-1.d0/6.d0, c11=5.d0/6.d0,  c12=2.d0/6.d0, &
                        c20=2.d0/6.d0, c21=5.d0/6.d0,  c22=-1.d0/6.d0
!
! Tchekhovskoy et al. (2008)
!
!  real(8), parameter :: cd0=1.d0/16.d0, cd1=10.d0/16.d0, cd2=5.d0/16.d0
!  real(8), parameter :: c00=3.d0/8.d0, c01=-10.d0/8.d0, c02=15.d0/8.d0, &
!                        c10=-1.d0/8.d0, c11=6.d0/8.d0,  c12=3.d0/8.d0, &
!                        c20=3.d0/8.d0, c21=6.d0/8.d0,  c22=-1.d0/8.d0
!
  real(8), parameter :: epsi=1.d-6, p0=1.d0
!
  real(8) :: betail0, betail1, betail2, omegil10, omegil11, omegil12, tauil,& 
             omegitl, omegil0, omegil1, omegil2, uriil10, uriil11, uriil12, &
             betair0, betair1, betair2, omegir10, omegir11, omegir12, tauir,&
             omegitr, omegir0, omegir1, omegir2, uriir10, uriir11, uriir12, &
             betajl0, betajl1, betajl2, omegjl10, omegjl11, omegjl12, taujl,&
             omegjtl, omegjl0, omegjl1, omegjl2, urijl10, urijl11, urijl12, &
             betajr0, betajr1, betajr2, omegjr10, omegjr11, omegjr12, taujr,&
             omegjtr, omegjr0, omegjr1, omegjr2, urijr10, urijr11, urijr12, &
             betakl0, betakl1, betakl2, omegkl10, omegkl11, omegkl12, taukl,&
             omegktl, omegkl0, omegkl1, omegkl2, urikl10, urikl11, urikl12, &
             betakr0, betakr1, betakr2, omegkr10, omegkr11, omegkr12, taukr,&
             omegktr, omegkr0, omegkr1, omegkr2, urikr10, urikr11, urikr12     
!
!-----------------------------------------------------------------------
!
  do k=ks1+2,ke1-2
    do j=js1+2,je1-2
      do i=is1+2,ie1-2
        do n=1,nv
!!
!! x-direction
!!
!
! smooth indicater
!
          betail0=(13./12.) &
            *(uri(n,i-2,j,k)-2.*uri(n,i-1,j,k)+uri(n,i,j,k))**2 &
            +(1./4.) &
            *(uri(n,i-2,j,k)-4.*uri(n,i-1,j,k)+3.*uri(n,i,j,k))**2 
          betail1=(13./12.) &
            *(uri(n,i-1,j,k)-2.*uri(n,i,j,k)+uri(n,i+1,j,k))**2 &
            +(1./4.) &
            *(uri(n,i-1,j,k)-uri(n,i+1,j,k))**2 
          betail2=(13./12.) &
            *(uri(n,i,j,k)-2.*uri(n,i+1,j,k)+uri(n,i+2,j,k))**2 &
            +(1./4.) &
            *(3.*uri(n,i,j,k)-4.*uri(n,i+1,j,k)+uri(n,i+2,j,k))**2
          
          betair0=(13./12.) &
            *(uri(n,i+2,j,k)-2.*uri(n,i+1,j,k)+uri(n,i,j,k))**2 &
            +(1./4.) &
            *(uri(n,i+2,j,k)-4.*uri(n,i+1,j,k)+3.*uri(n,i,j,k))**2 
          betair1=(13./12.) &
            *(uri(n,i+1,j,k)-2.*uri(n,i,j,k)+uri(n,i-1,j,k))**2 &
            +(1./4.) &
            *(uri(n,i+1,j,k)-uri(n,i-1,j,k))**2
          betair2=(13./12.) &
            *(uri(n,i,j,k)-2.*uri(n,i-1,j,k)+uri(n,i-2,j,k))**2 &
            +(1./4.) &
            *(3.*uri(n,i,j,k)-4.*uri(n,i-1,j,k)+uri(n,i-2,j,k))**2
!     
! unnormalized weights     
!     
          tauil=abs(betail0-betail2)
          omegil10=cd0*(1.0+(tauil/(betail0+epsi)**p0))
          omegil11=cd1*(1.0+(tauil/(betail1+epsi)**p0))
          omegil12=cd2*(1.0+(tauil/(betail2+epsi)**p0))
          omegitl=omegil10+omegil11+omegil12
!
          tauir=abs(betair0-betair2)
          omegir10=cd0*(1.0+(tauir/(betair0+epsi)**p0))
          omegir11=cd1*(1.0+(tauir/(betair1+epsi)**p0))
          omegir12=cd2*(1.0+(tauir/(betair2+epsi)**p0))
          omegitr=omegir10+omegir11+omegir12
!
! normalized weights
!
          omegil0=omegil10/omegitl
          omegil1=omegil11/omegitl
          omegil2=omegil12/omegitl

          omegir0=omegir10/omegitr
          omegir1=omegir11/omegitr
          omegir2=omegir12/omegitr
!
! reconstruction function
!
          uriil10=c00*uri(n,i-2,j,k)+c01*uri(n,i-1,j,k)+c02*uri(n,i,j,k)
          uriil11=c10*uri(n,i-1,j,k)+c11*uri(n,i,j,k)+c12*uri(n,i+1,j,k)
          uriil12=c20*uri(n,i,j,k)+c21*uri(n,i+1,j,k)+c22*uri(n,i+2,j,k)
          
          uriir10=c00*uri(n,i+2,j,k)+c01*uri(n,i+1,j,k)+c02*uri(n,i,j,k)
          uriir11=c10*uri(n,i+1,j,k)+c11*uri(n,i,j,k)+c12*uri(n,i-1,j,k)
          uriir12=c20*uri(n,i,j,k)+c21*uri(n,i-1,j,k)+c22*uri(n,i-2,j,k)
!
! cell interface values
!
          uriil(n,i,j,k)=omegil0*uriil10+omegil1*uriil11+omegil2*uriil12
          uriir(n,i-1,j,k)=omegir0*uriir10+omegir1*uriir11+omegir2*uriir12

!!
!! y-direction
!!
!
! smooth indicater
!
          betajl0=(13./12.) &
            *(uri(n,i,j-2,k)-2.0*uri(n,i,j-1,k)+uri(n,i,j,k))**2 &
            +(1./4.) &
            *(uri(n,i,j-2,k)-4.*uri(n,i,j-1,k)+3.*uri(n,i,j,k))**2
          betajl1=(13./12.) &
            *(uri(n,i,j-1,k)-2.*uri(n,i,j,k)+uri(n,i,j+1,k))**2 &
            +(1./4.) &
            *(uri(n,i,j-1,k)-uri(n,i,j+1,k))**2
          betajl2=(13./12.) &
            *(uri(n,i,j,k)-2.*uri(n,i,j+1,k)+uri(n,i,j+2,k))**2 &
            +(1./4.) &
            *(3.*uri(n,i,j,k)-4.*uri(n,i,j+1,k)+uri(n,i,j+2,k))**2
          
          betajr0=(13./12.) &
            *(uri(n,i,j+2,k)-2.*uri(n,i,j+1,k)+uri(n,i,j,k))**2 &
            +(1./4.) &
            *(uri(n,i,j+2,k)-4.*uri(n,i,j+1,k)+3.*uri(n,i,j,k))**2
          betajr1=(13./12.) &
            *(uri(n,i,j+1,k)-2.*uri(n,i,j,k)+uri(n,i,j-1,k))**2 &
            +(1./4.) &
            *(uri(n,i,j+1,k)-uri(n,i,j-1,k))**2
          betajr2=(13./12.) &
            *(uri(n,i,j,k)-2.0*uri(n,i,j-1,k)+uri(n,i,j-2,k))**2 &
            +(1./4.) &
            *(3.*uri(n,i,j,k)-4.*uri(n,i,j-1,k)+uri(n,i,j-2,k))**2
!     
! unnormalized weights     
!     
          taujl=abs(betajl0-betajl2)
          omegjl10=cd0*(1.0+(taujl/(betajl0+epsi)**p0))
          omegjl11=cd1*(1.0+(taujl/(betajl1+epsi)**p0))
          omegjl12=cd2*(1.0+(taujl/(betajl2+epsi)**p0))
          omegjtl=omegjl10+omegjl11+omegjl12

          taujr=abs(betajr0-betajr2)
          omegjr10=cd0*(1.0+(taujr/(betajr0+epsi)**p0))
          omegjr11=cd1*(1.0+(taujr/(betajr1+epsi)**p0))
          omegjr12=cd2*(1.0+(taujr/(betajr2+epsi)**p0))
          omegjtr=omegjr10+omegjr11+omegjr12
!
! normalized weights
!
          omegjl0=omegjl10/omegjtl
          omegjl1=omegjl11/omegjtl
          omegjl2=omegjl12/omegjtl

          omegjr0=omegjr10/omegjtr
          omegjr1=omegjr11/omegjtr
          omegjr2=omegjr12/omegjtr
!
! reconstruction function
!
          urijl10=c00*uri(n,i,j-2,k)+c01*uri(n,i,j-1,k)+c02*uri(n,i,j,k)
          urijl11=c10*uri(n,i,j-1,k)+c11*uri(n,i,j,k)+c12*uri(n,i,j+1,k)
          urijl12=c20*uri(n,i,j,k)+c21*uri(n,i,j+1,k)+c22*uri(n,i,j+2,k)
          
          urijr10=c00*uri(n,i,j+2,k)+c01*uri(n,i,j+1,k)+c02*uri(n,i,j,k)
          urijr11=c10*uri(n,i,j+1,k)+c11*uri(n,i,j,k)+c12*uri(n,i,j-1,k)
          urijr12=c20*uri(n,i,j,k)+c21*uri(n,i,j-1,k)+c22*uri(n,i,j-2,k)
!
! cell interface values
!
          urijl(n,i,j,k)=omegjl0*urijl10+omegjl1*urijl11 +omegjl2*urijl12
          urijr(n,i,j-1,k)=omegjr0*urijr10+omegjr1*urijr11+omegjr2*urijr12
!!
!! z-direction
!!
!
! smooth indicater
!
          betakl0=(13./12.) &
            *(uri(n,i,j,k-2)-2.0*uri(n,i,j,k-1)+uri(n,i,j,k))**2 &
            +(1./4.) &
            *(uri(n,i,j,k-2)-4.*uri(n,i,j,k-1)+3.*uri(n,i,j,k))**2
          betakl1=(13./12.) &
            *(uri(n,i,j,k-1)-2.*uri(n,i,j,k)+uri(n,i,j,k+1))**2 &
            +(1./4.) &
            *(uri(n,i,j,k-1)-uri(n,i,j,k+1))**2
          betakl2=(13./12.) &
            *(uri(n,i,j,k)-2.*uri(n,i,j,k+1)+uri(n,i,j,k+2))**2 &
            +(1./4.) &
            *(3.*uri(n,i,j,k)-4.*uri(n,i,j,k+1)+uri(n,i,j,k+2))**2
          
          betakr0=(13./12.) &
            *(uri(n,i,j,k+2)-2.*uri(n,i,j,k+1)+uri(n,i,j,k))**2 &
            +(1./4.) &
            *(uri(n,i,j,k+2)-4.*uri(n,i,j,k+1)+3.*uri(n,i,j,k))**2
          betakr1=(13./12.) &
            *(uri(n,i,j,k+1)-2.*uri(n,i,j,k)+uri(n,i,j,k-1))**2 &
            +(1./4.) &
            *(uri(n,i,j,k+1)-uri(n,i,j,k-1))**2
          betakr2=(13./12.) &
            *(uri(n,i,j,k)-2.0*uri(n,i,j,k-1)+uri(n,i,j,k-2))**2 &
            +(1./4.) &
            *(3.*uri(n,i,j,k)-4.*uri(n,i,j,k-1)+uri(n,i,j,k-2))**2
!     
! unnormalized weights     
!         
          taukl=abs(betakl0-betakl2)
          omegkl10=cd0*(1.0+(taukl/(betakl0+epsi)**p0))
          omegkl11=cd1*(1.0+(taukl/(betakl1+epsi)**p0))
          omegkl12=cd2*(1.0+(taukl/(betakl2+epsi)**p0))
          omegktl=omegkl10+omegkl11+omegkl12

          taukr=abs(betakr0-betakr2)
          omegkr10=cd0*(1.0+(taukr/(betakr0+epsi)**p0))
          omegkr11=cd1*(1.0+(taukr/(betakr1+epsi)**p0))
          omegkr12=cd2*(1.0+(taukr/(betakr2+epsi)**p0))
          omegktr=omegkr10+omegkr11+omegkr12
!
! normalized weights
!
          omegkl0=omegkl10/omegktl
          omegkl1=omegkl11/omegktl
          omegkl2=omegkl12/omegktl

          omegkr0=omegkr10/omegktr
          omegkr1=omegkr11/omegktr
          omegkr2=omegkr12/omegktr
!
! reconstruction function
!
          urikl10=c00*uri(n,i,j,k-2)+c01*uri(n,i,j,k-1)+c02*uri(n,i,j,k)
          urikl11=c10*uri(n,i,j,k-1)+c11*uri(n,i,j,k)+c12*uri(n,i,j,k+1)
          urikl12=c20*uri(n,i,j,k)+c21*uri(n,i,j,k+1)+c22*uri(n,i,j,k+2)
          
          urikr10=c00*uri(n,i,j,k+2)+c01*uri(n,i,j,k+1)+c02*uri(n,i,j,k)
          urikr11=c10*uri(n,i,j,k+1)+c11*uri(n,i,j,k)+c12*uri(n,i,j,k-1)
          urikr12=c20*uri(n,i,j,k)+c21*uri(n,i,j,k-1)+c22*uri(n,i,j,k-2)
!
! cell interface values
!
          urikl(n,i,j,k)=omegkl0*urikl10+omegkl1*urikl11+omegkl2*urikl12
          urikr(n,i,j,k-1)=omegkr0*urikr10+omegkr1*urikr11+omegkr2*urikr12
    
        enddo
      enddo
    enddo
  enddo
  
  return
end subroutine weno5z
!
!---------------------------------------------------------------------@
subroutine weno5m(uri,uriir,urijr,urikr,uriil,urijl,urikl, &
                  is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Mapped Weighted Essentially Non-Oscillatory scheme 
!     5th order + modified weight
!     (Henrick et al. 2005, JCP, 207, 542)
!     Calculate cell-interface variables(urir, uril) 
!     from cell-center variables (uu, uri)
!
  use pram, only : imax, jmax, kmax, nv
  implicit none
!
  integer :: i, j, k, n, is1, ie1, js1, je1, ks1, ke1

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)
!
! Jiang & Shu (1990)
!
  real(8), parameter :: cd0=1.d0/10.d0, cd1=6.d0/10.d0, cd2=3.d0/10.d0
  real(8), parameter :: c00=2.d0/6.d0, c01=-7.d0/6.d0, c02=11.d0/6.d0, &
                        c10=-1.d0/6.d0, c11=5.d0/6.d0,  c12=2.d0/6.d0, &
                        c20=2.d0/6.d0, c21=5.d0/6.d0,  c22=-1.d0/6.d0
!
  real(8), parameter :: epsi=1.d-40
!
  real(8) :: betail0, betail1, betail2, omegail10, omegail11, omegail12, &
             omegaitl, omegail0, omegail1, omegail2, &
             omegil10, omegil11, omegil12, omegitl, omegil0, omegil1, omegil2,& 
             uriil10, uriil11, uriil12, &
             betair0, betair1, betair2, omegair10, omegair11, omegair12, & 
             omegaitr, omegair0, omegair1, omegair2, &
             omegir10, omegir11, omegir12, omegitr, omegir0, omegir1, omegir2,& 
             uriir10, uriir11, uriir12, &
             betajl0, betajl1, betajl2, omegajl10, omegajl11, omegajl12, & 
             omegajtl, omegajl0, omegajl1, omegajl2, &
             omegjl10, omegjl11, omegjl12, omegjtl, omegjl0, omegjl1, omegjl2,&
             urijl10, urijl11, urijl12, &
             betajr0, betajr1, betajr2, omegajr10, omegajr11, omegajr12, & 
             omegajtr, omegajr0, omegajr1, omegajr2, &
             omegjr10, omegjr11, omegjr12, omegjtr, omegjr0, omegjr1, omegjr2,& 
             urijr10, urijr11, urijr12, &
             betakl0, betakl1, betakl2, omegakl10, omegakl11, omegakl12, & 
             omegaktl, omegakl0, omegakl1, omegakl2, &
             omegkl10, omegkl11, omegkl12, omegktl, omegkl0, omegkl1, omegkl2,&
             urikl10, urikl11, urikl12, &
             betakr0, betakr1, betakr2, omegakr10, omegakr11, omegakr12, & 
             omegaktr, omegakr0, omegakr1, omegakr2, &
             omegkr10, omegkr11, omegkr12, omegktr, omegkr0, omegkr1, omegkr2,&
             urikr10, urikr11, urikr12
      
!
!-----------------------------------------------------------------------
!
  do k=ks1+2,ke1-2
    do j=js1+2,je1-2
      do i=is1+2,ie1-2
        do n=1,nv
!!
!! x-direction
!!
!
! smooth indicater
!
          betail0=(13./12.) &
            *(uri(n,i-2,j,k)-2.*uri(n,i-1,j,k)+uri(n,i,j,k))**2 &
            +(1./4.) &
            *(uri(n,i-2,j,k)-4.*uri(n,i-1,j,k)+3.*uri(n,i,j,k))**2
          betail1=(13./12.) &
            *(uri(n,i-1,j,k)-2.*uri(n,i,j,k)+uri(n,i+1,j,k))**2 &
            +(1./4.) &
            *(uri(n,i-1,j,k)-uri(n,i+1,j,k))**2
          betail2=(13./12.) &
            *(uri(n,i,j,k)-2.*uri(n,i+1,j,k)+uri(n,i+2,j,k))**2 &
            +(1./4.) &
            *(3.*uri(n,i,j,k)-4.*uri(n,i+1,j,k)+uri(n,i+2,j,k))**2
            
          betair0=(13./12.) &
            *(uri(n,i+2,j,k)-2.*uri(n,i+1,j,k)+uri(n,i,j,k))**2 &
            +(1./4.) &
            *(uri(n,i+2,j,k)-4.*uri(n,i+1,j,k)+3.*uri(n,i,j,k))**2
          betair1=(13./12.) &
            *(uri(n,i+1,j,k)-2.*uri(n,i,j,k)+uri(n,i-1,j,k))**2 &
            +(1./4.) &
            *(uri(n,i+1,j,k)-uri(n,i-1,j,k))**2
          betair2=(13./12.) &
            *(uri(n,i,j,k)-2.*uri(n,i-1,j,k)+uri(n,i-2,j,k))**2 &
            +(1./4.) &
            *(3.*uri(n,i,j,k)-4.*uri(n,i-1,j,k)+uri(n,i-2,j,k))**2
!     
! unnormalized weights     
!     
          omegail10=cd0/(betail0+epsi)**2
          omegail11=cd1/(betail1+epsi)**2
          omegail12=cd2/(betail2+epsi)**2
          omegaitl=omegail10+omegail11+omegail12

          omegair10=cd0/(betair0+epsi)**2
          omegair11=cd1/(betair1+epsi)**2
          omegair12=cd2/(betair2+epsi)**2
          omegaitr=omegair10+omegair11+omegair12
!
! normalized weights
!
          omegail0=omegail10/omegaitl
          omegail1=omegail11/omegaitl
          omegail2=omegail12/omegaitl

          omegair0=omegair10/omegaitr
          omegair1=omegair11/omegaitr
          omegair2=omegair12/omegaitr
!
! Mapped weights
!
          omegil10=omegail0*(cd0+cd0**2-3.*cd0*omegail0+omegail0**2) &
                   /(cd0**2+omegail0*(1.-2.*cd0)) 
          omegil11=omegail1*(cd1+cd1**2-3.*cd1*omegail1+omegail1**2) &
                   /(cd1**2+omegail1*(1.-2.*cd1))
          omegil12=omegail2*(cd2+cd2**2-3.*cd2*omegail2+omegail2**2) &
                   /(cd2**2+omegail2*(1.-2.*cd2))
          omegitl=omegil10+omegil11+omegil12

          omegir10=omegair0*(cd0+cd0**2-3.*cd0*omegair0+omegair0**2) &
                   /(cd0**2+omegair0*(1.-2.*cd0)) 
          omegir11=omegair1*(cd1+cd1**2-3.*cd1*omegair1+omegair1**2) &
                   /(cd1**2+omegair1*(1.-2.*cd1))
          omegir12=omegair2*(cd2+cd2**2-3.*cd2*omegair2+omegair2**2) &
                   /(cd2**2+omegair2*(1.-2.*cd2))
          omegitr=omegir10+omegir11+omegir12
 
          omegil0=omegil10/omegitl
          omegil1=omegil11/omegitl
          omegil2=omegil12/omegitl

          omegir0=omegir10/omegitr
          omegir1=omegir11/omegitr
          omegir2=omegir12/omegitr 
!
! reconstruction function
!
          uriil10=c00*uri(n,i-2,j,k)+c01*uri(n,i-1,j,k)+c02*uri(n,i,j,k)
          uriil11=c10*uri(n,i-1,j,k)+c11*uri(n,i,j,k)+c12*uri(n,i+1,j,k)
          uriil12=c20*uri(n,i,j,k)+c21*uri(n,i+1,j,k)+c22*uri(n,i+2,j,k)
          
          uriir10=c00*uri(n,i+2,j,k)+c01*uri(n,i+1,j,k)+c02*uri(n,i,j,k)
          uriir11=c10*uri(n,i+1,j,k)+c11*uri(n,i,j,k)+c12*uri(n,i-1,j,k)
          uriir12=c20*uri(n,i,j,k)+c21*uri(n,i-1,j,k)+c22*uri(n,i-2,j,k)
!
! cell interface values
!
          uriil(n,i,j,k)=omegil0*uriil10+omegil1*uriil11+omegil2*uriil12
          uriir(n,i-1,j,k)=omegir0*uriir10+omegir1*uriir11+omegir2*uriir12

!!
!! y-direction
!!
!
! smooth indicater
!
          betajl0=(13./12.) &
            *(uri(n,i,j-2,k)-2.0*uri(n,i,j-1,k)+uri(n,i,j,k))**2 &
            +(1./4.) &
            *(uri(n,i,j-2,k)-4.*uri(n,i,j-1,k)+3.*uri(n,i,j,k))**2
          betajl1=(13./12.) &
            *(uri(n,i,j-1,k)-2.*uri(n,i,j,k)+uri(n,i,j+1,k))**2 &
            +(1./4.) &
            *(uri(n,i,j-1,k)-uri(n,i,j+1,k))**2
          betajl2=(13./12.) &
            *(uri(n,i,j,k)-2.*uri(n,i,j+1,k)+uri(n,i,j+2,k))**2 &
            +(1./4.) &
            *(3.*uri(n,i,j,k)-4.*uri(n,i,j+1,k)+uri(n,i,j+2,k))**2
          
          betajr0=(13./12.) &
            *(uri(n,i,j+2,k)-2.*uri(n,i,j+1,k)+uri(n,i,j,k))**2 &
            +(1./4.) &
            *(uri(n,i,j+2,k)-4.*uri(n,i,j+1,k)+3.*uri(n,i,j,k))**2 
          betajr1=(13./12.) &
            *(uri(n,i,j+1,k)-2.*uri(n,i,j,k)+uri(n,i,j-1,k))**2 &
            +(1./4.) &
            *(uri(n,i,j+1,k)-uri(n,i,j-1,k))**2
          betajr2=(13./12.) &
            *(uri(n,i,j,k)-2.0*uri(n,i,j-1,k)+uri(n,i,j-2,k))**2 &
            +(1./4.) &
            *(3.*uri(n,i,j,k)-4.*uri(n,i,j-1,k)+uri(n,i,j-2,k))**2
!     
! unnormalized weights     
!     
          omegajl10=cd0/(betajl0+epsi)**2
          omegajl11=cd1/(betajl1+epsi)**2
          omegajl12=cd2/(betajl2+epsi)**2
          omegajtl=omegajl10+omegajl11+omegajl12

          omegajr10=cd0/(betajr0+epsi)**2
          omegajr11=cd1/(betajr1+epsi)**2
          omegajr12=cd2/(betajr2+epsi)**2
          omegajtr=omegajr10+omegajr11+omegajr12
!
! normalized weights
!
          omegajl0=omegajl10/omegajtl
          omegajl1=omegajl11/omegajtl
          omegajl2=omegajl12/omegajtl

          omegajr0=omegajr10/omegajtr
          omegajr1=omegajr11/omegajtr
          omegajr2=omegajr12/omegajtr
!
! Mapped weights
!
          omegjl10=omegajl0*(cd0+cd0**2-3.*cd0*omegajl0+omegajl0**2) &
                   /(cd0**2+omegajl0*(1.-2.*cd0)) 
          omegjl11=omegajl1*(cd1+cd1**2-3.*cd1*omegajl1+omegajl1**2) &
                   /(cd1**2+omegajl1*(1.-2.*cd1))
          omegjl12=omegajl2*(cd2+cd2**2-3.*cd2*omegajl2+omegajl2**2) &
                   /(cd2**2+omegajl2*(1.-2.*cd2))
          omegjtl=omegjl10+omegjl11+omegjl12

          omegjr10=omegajr0*(cd0+cd0**2-3.*cd0*omegajr0+omegajr0**2) &
                   /(cd0**2+omegajr0*(1.-2.*cd0)) 
          omegjr11=omegajr1*(cd1+cd1**2-3.*cd1*omegajr1+omegajr1**2) &
                   /(cd1**2+omegajr1*(1.-2.*cd1))
          omegjr12=omegajr2*(cd2+cd2**2-3.*cd2*omegajr2+omegajr2**2) &
                   /(cd2**2+omegajr2*(1.-2.*cd2))
          omegjtr=omegjr10+omegjr11+omegjr12
 
          omegjl0=omegjl10/omegjtl
          omegjl1=omegjl11/omegjtl
          omegjl2=omegjl12/omegjtl

          omegjr0=omegjr10/omegjtr
          omegjr1=omegjr11/omegjtr
          omegjr2=omegjr12/omegjtr 
!
! reconstruction function
!
          urijl10=c00*uri(n,i,j-2,k)+c01*uri(n,i,j-1,k)+c02*uri(n,i,j,k)
          urijl11=c10*uri(n,i,j-1,k)+c11*uri(n,i,j,k)+c12*uri(n,i,j+1,k)
          urijl12=c20*uri(n,i,j,k)+c21*uri(n,i,j+1,k)+c22*uri(n,i,j+2,k)
          
          urijr10=c00*uri(n,i,j+2,k)+c01*uri(n,i,j+1,k)+c02*uri(n,i,j,k)
          urijr11=c10*uri(n,i,j+1,k)+c11*uri(n,i,j,k)+c12*uri(n,i,j-1,k)
          urijr12=c20*uri(n,i,j,k)+c21*uri(n,i,j-1,k)+c22*uri(n,i,j-2,k)
!
! cell interface values
!
          urijl(n,i,j,k)=omegjl0*urijl10+omegjl1*urijl11+omegjl2*urijl12
          urijr(n,i,j-1,k)=omegjr0*urijr10+omegjr1*urijr11+omegjr2*urijr12
!!
!! z-direction
!!
!
! smooth indicater
!
          betakl0=(13./12.) &
            *(uri(n,i,j,k-2)-2.0*uri(n,i,j,k-1)+uri(n,i,j,k))**2 &
            +(1./4.) &
            *(uri(n,i,j,k-2)-4.*uri(n,i,j,k-1)+3.*uri(n,i,j,k))**2
          betakl1=(13./12.) &
            *(uri(n,i,j,k-1)-2.*uri(n,i,j,k)+uri(n,i,j,k+1))**2 &
            +(1./4.) &
            *(uri(n,i,j,k-1)-uri(n,i,j,k+1))**2
          betakl2=(13./12.) &
            *(uri(n,i,j,k)-2.*uri(n,i,j,k+1)+uri(n,i,j,k+2))**2 &
            +(1./4.) &
            *(3.*uri(n,i,j,k)-4.*uri(n,i,j,k+1)+uri(n,i,j,k+2))**2
          
          betakr0=(13./12.) &
            *(uri(n,i,j,k+2)-2.*uri(n,i,j,k+1)+uri(n,i,j,k))**2 &
            +(1./4.) &
            *(uri(n,i,j,k+2)-4.*uri(n,i,j,k+1)+3.*uri(n,i,j,k))**2 
          betakr1=(13./12.) &
            *(uri(n,i,j,k+1)-2.*uri(n,i,j,k)+uri(n,i,j,k-1))**2 &
            +(1./4.) &
            *(uri(n,i,j,k+1)-uri(n,i,j,k-1))**2
          betakr2=(13./12.) &
            *(uri(n,i,j,k)-2.0*uri(n,i,j,k-1)+uri(n,i,j,k-2))**2 &
            +(1./4.) &
            *(3.*uri(n,i,j,k)-4.*uri(n,i,j,k-1)+uri(n,i,j,k-2))**2
!     
! unnormalized weights     
!     
          omegakl10=cd0/(betakl0+epsi)**2
          omegakl11=cd1/(betakl1+epsi)**2
          omegakl12=cd2/(betakl2+epsi)**2
          omegaktl=omegakl10+omegakl11+omegakl12

          omegakr10=cd0/(betakr0+epsi)**2
          omegakr11=cd1/(betakr1+epsi)**2
          omegakr12=cd2/(betakr2+epsi)**2
          omegaktr=omegakr10+omegakr11+omegakr12
!
! normalized weights
!
          omegakl0=omegakl10/omegaktl
          omegakl1=omegakl11/omegaktl
          omegakl2=omegakl12/omegaktl

          omegakr0=omegakr10/omegaktr
          omegakr1=omegakr11/omegaktr
          omegakr2=omegakr12/omegaktr
!
! Mapped weights
!
          omegkl10=omegakl0*(cd0+cd0**2-3.*cd0*omegakl0+omegakl0**2) &
                   /(cd0**2+omegakl0*(1.-2.*cd0)) 
          omegkl11=omegakl1*(cd1+cd1**2-3.*cd1*omegakl1+omegakl1**2) &
                   /(cd1**2+omegakl1*(1.-2.*cd1))
          omegkl12=omegakl2*(cd2+cd2**2-3.*cd2*omegakl2+omegakl2**2) &
                   /(cd2**2+omegakl2*(1.-2.*cd2))
          omegktl=omegkl10+omegkl11+omegkl12

          omegkr10=omegakr0*(cd0+cd0**2-3.*cd0*omegakr0+omegakr0**2) &
                   /(cd0**2+omegakr0*(1.-2.*cd0)) 
          omegkr11=omegakr1*(cd1+cd1**2-3.*cd1*omegakr1+omegakr1**2) &
                   /(cd1**2+omegakr1*(1.-2.*cd1))
          omegkr12=omegakr2*(cd2+cd2**2-3.*cd2*omegakr2+omegakr2**2) &
                   /(cd2**2+omegakr2*(1.-2.*cd2))
          omegktr=omegkr10+omegkr11+omegkr12
 
          omegkl0=omegkl10/omegktl
          omegkl1=omegkl11/omegktl
          omegkl2=omegkl12/omegktl

          omegkr0=omegkr10/omegktr
          omegkr1=omegkr11/omegktr
          omegkr2=omegkr12/omegktr 
!
! reconstruction function
!
          urikl10=c00*uri(n,i,j,k-2)+c01*uri(n,i,j,k-1)+c02*uri(n,i,j,k)
          urikl11=c10*uri(n,i,j,k-1)+c11*uri(n,i,j,k)+c12*uri(n,i,j,k+1)
          urikl12=c20*uri(n,i,j,k)+c21*uri(n,i,j,k+1)+c22*uri(n,i,j,k+2)
          
          urikr10=c00*uri(n,i,j,k+2)+c01*uri(n,i,j,k+1)+c02*uri(n,i,j,k)
          urikr11=c10*uri(n,i,j,k+1)+c11*uri(n,i,j,k)+c12*uri(n,i,j,k-1)
          urikr12=c20*uri(n,i,j,k)+c21*uri(n,i,j,k-1)+c22*uri(n,i,j,k-2)
!
! cell interface values
!
          urikl(n,i,j,k)=omegkl0*urikl10+omegkl1*urikl11+omegkl2*urikl12
          urikr(n,i,j,k-1)=omegkr0*urikr10+omegkr1*urikr11+omegkr2*urikr12
    
        enddo
      enddo
    enddo
  enddo

  return
end subroutine weno5m
!
!---------------------------------------------------------------------@
subroutine lim03(uri,x1,x2,x3,x1b,x2b,x3b,uriir,urijr,urikr,uriil,urijl,urikl, &
                 is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     3rd-order Limiter function (Cada & Torrilhon 2009)
!     Calculate cell-interface variables(urir, uril) 
!     from cell-center variables (uu, uri)
!
  use pram, only : imax, jmax, kmax, nv, metric
  implicit none
!
  integer :: i, j, k, n, is1, ie1, js1, je1, ks1, ke1

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: x1(imax), x2(jmax), x3(kmax), &
             x1b(imax), x2b(jmax), x3b(kmax)
  real(8), parameter :: epsi=1.d-12, r0=1.0d0

  real(8) :: thetail, thetair, p3il, p3ir, etail, etair, & 
             thetajl, thetajr, p3jl, p3jr, etajl, etajr, &
             thetakl, thetakr, p3kl, p3kr, etakl, etakr, &
             tmp1ail, tmp1air, tmp1bil, tmp1bir, tmp3il, tmp3ir, &
             tmp1ajl, tmp1ajr, tmp1bjl, tmp1bjr, tmp3jl, tmp3jr, &
             tmp1akl, tmp1akr, tmp1bkl, tmp1bkr, tmp3kl, tmp3kr, &
             phiil, phiir, alphail, alphair, &
             phijl, phijr, alphajl, alphajr, &
             phikl, phikr, alphakl, alphakr
  real(8) :: x1p, x1c, x1m, x2p, x2c, x2m, x3p, x3c, x3m
  real(8) :: uriip, uriim, urijp, urijm, urikp, urikm
!
!-----------------------------------------------------------------------
!
  do k=ks1+1,ke1-1
    do j=js1+1,je1-1
      do i=is1+1,ie1-1
        do n=1,nv
!!
!! x-direction
!!
          if(metric .eq. 403) then
            x1p=x1b(i+1)
            x1c=x1b(i)
            x1m=x1b(i-1)
          else
            x1p=x1(i+1)
            x1c=x1(i)
            x1m=x1(i-1)
          endif  
!
          uriip=uri(n,i+1,j,k)-uri(n,i,j,k)
          uriim=uri(n,i,j,k)-uri(n,i-1,j,k)
!
! slope ratio
!
          if(uriip .eq. 0.d0) then
            thetail=0.d0
          else
            thetail=uriim*(1./uriip)
          endif
          if(uriim .eq. 0.d0) then
            thetair=0.d0
          else
            thetair=(-uriip)*(1./(-uriim))
          endif
!
! building block
!          
          p3il=(2.0+thetail)/3.0
          p3ir=(2.0+thetair)/3.0
!
! 3rd-order limitter
!
          if(thetail .ge. 0.d0) then
            tmp1ail=min(p3il, 2.d0*thetail, 1.6d0)
            phiil=max(0.d0,tmp1ail)
          else
            tmp1bil=min(p3il,-0.5d0*thetail)
            phiil=max(0.d0,tmp1bil)
          endif
          
          if(thetair .ge. 0.d0) then
            tmp1air=min(p3ir, 2.d0*thetair, 1.6d0)
            phiir=max(0.d0,tmp1air)
          else
            tmp1bir=min(p3ir,-0.5d0*thetair)
            phiir=max(0.d0,tmp1bir)
          endif
!
! smooth switch
!
          etail=((uriim*uriim)+(uriip*uriip)) &
               *(1./(r0*(x1p-x1c))**2)
          tmp3il=min(1.d0, 0.5d0+(etail-1.d0)*(1.d0/(2.d0*epsi)))
          alphail=max(0.d0, tmp3il)
          
          etair=(((-uriip)*(-uriip))+((-uriim)*(-uriim))) &
               *(1./(r0*(x1c-x1m))**2)
          tmp3ir=min(1.d0, 0.5d0+(etair-1.d0)*(1.d0/(2.d0*epsi)))
          alphair=max(0.d0, tmp3ir)
          
          uriil(n,i,j,k)=uri(n,i,j,k) &
                         +0.5*uriip*(p3il+alphail*(phiil-p3il))
          uriir(n,i-1,j,k)=uri(n,i,j,k) &
                         +0.5*(-uriim)*(p3ir+alphair*(phiir-p3ir))
!!
!! y-direction
!!
          if(metric .eq. 403) then
            x2p=x2b(j+1)
            x2c=x2b(j)
            x2m=x2b(j-1)
          else
            x2p=x2(j+1)
            x2c=x2(j)
            x2m=x2(j-1)
          endif  
!
          urijp=uri(n,i,j+1,k)-uri(n,i,j,k)
          urijm=uri(n,i,j,k)-uri(n,i,j-1,k)
!
          if(urijp .eq. 0.d0) then
            thetajl=0.d0
          else
            thetajl=urijm*(1./urijp)
          endif
          if(urijm .eq. 0.d0) then
            thetajr=0.d0
          else
            thetajr=(-urijp)*(1./(-urijm))
          endif
          
          p3jl=(2.0+thetajl)/3.0
          p3jr=(2.0+thetajr)/3.0
          
          if(thetajl .ge. 0.d0) then
            tmp1ajl=min(p3jl, 2.d0*thetajl, 1.6d0)
            phijl=max(0.d0, tmp1ajl)
          else
            tmp1bjl=min(p3jl,-0.d5*thetajl)
            phijl=max(0.d0,tmp1bjl)
          endif
          
          if(thetajr .ge. 0.d0) then
            tmp1ajr=min(p3jr, 2.d0*thetajr, 1.6d0)
            phijr=max(0.d0,tmp1ajr)
          else
            tmp1bjr=min(p3jr,-0.d5*thetajr)
            phijr=max(0.d0,tmp1bjr)
          endif
          
          etajl=((urijm*urijm)+(urijp*urijp)) &
               *(1./(r0*(x2p-x2c))**2)
          tmp3jl=min(1.d0, 0.5d0+(etajl-1.d0)/(2.d0*epsi))
          alphajl=max(0.d0, tmp3jl)
          
          etajr=(((-urijp)*(-urijp))+((-urijm)*(-urijm))) &
                *(1./(r0*(x2c-x2m))**2)
          tmp3jr=min(1.d0, 0.5d0+(etajr-1.d0)/(2.d0*epsi))
          alphajr=max(0.d0, tmp3jr)
          
          urijl(n,i,j,k)=uri(n,i,j,k) &
                        +0.5*urijp*(p3jl+alphajl*(phijl-p3jl))
          urijr(n,i,j-1,k)=uri(n,i,j,k) &
                        +0.5*(-urijm)*(p3jr+alphajr*(phijr-p3jr))
!!
!! z-direction
!!
          if(metric .eq. 403) then
            x3p=x3b(k+1)
            x3c=x3b(k)
            x3m=x3b(k-1)
          else
            x3p=x3(k+1)
            x3c=x3(k)
            x3m=x3(k-1)
          endif  
!
          urikp=uri(n,i,j,k+1)-uri(n,i,j,k)
          urikm=uri(n,i,j,k)-uri(n,i,j,k-1)
          
          if(urikp .eq. 0.d0) then
            thetakl=0.d0
          else
            thetakl=urikm*(1./urikp)
          endif
          if(urikm .eq. 0.d0) then
            thetakr=0.d0
          else
            thetakr=(-urikp)*(1./(-urikm))
          endif

          p3kl=(2.0+thetakl)/3.0
          p3kr=(2.0+thetakr)/3.0
          
          if(thetakl .ge. 0.d0) then
            tmp1akl=min(p3kl, 2.d0*thetakl, 1.6d0)
            phikl=max(0.d0,tmp1akl)
          else
            tmp1bkl=min(p3kl,-0.5d0*thetakl)
            phikl=max(0.d0,tmp1bkl)
          endif
          
          if(thetakr .ge. 0.d0) then
            tmp1akr=min(p3kr, 2.d0*thetakr, 1.6d0)
            phikr=max(0.d0,tmp1akr)
          else
            tmp1bkr=min(p3kr,-0.5d0*thetakr)
            phikr=max(0.d0,tmp1bkr)
          endif
          
          etakl=((urikm*urikm)+(urikp*urikp)) &
                *(1./(r0*(x3p-x3c))**2)
          tmp3kl=min(1.d0, 0.5d0+(etakl-1.d0)/(2.d0*epsi))
          alphakl=max(0.d0, tmp3kl)
          
          etakr=(((-urikp)*(-urikp))+((-urikm)*(-urikm))) &
                *(1./(r0*(x3c-x3m))**2)
          tmp3kr=min(1.d0, 0.5d0+(etakr-1.d0)/(2.d0*epsi))
          alphakr=max(0.d0, tmp3kr)
          
          urikl(n,i,j,k)=uri(n,i,j,k) &
                        +0.5*urikp*(p3kl+alphakl*(phikl-p3kl))
          urikr(n,i,j,k-1)=uri(n,i,j,k) &
                        +0.5*(-urikm)*(p3kr+alphakr*(phikr-p3kr))

        enddo
      enddo
    enddo
  enddo

  return
end subroutine lim03
!
!---------------------------------------------------------------------@
subroutine lim03a(uri,x1,x2,x3,x1a,x2a,x3a,x1b,x2b,x3b,x1ab,x2ab,x3ab,&
                 uriir,urijr,urikr,uriil,urijl,urikl,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Reconstruction step with Lim03 limiter 
!     Calculate cell-interface variables(uril, urir, wwl, wwr) 
!     from cell-center variables (uri)
!
  use pram, only : imax, jmax, kmax, nv, metric
  implicit none
!
  integer :: i, j, k, n, is1, ie1, js1, je1, ks1, ke1
!
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: x1(imax), x2(jmax), x3(kmax), & !- cell center position
             x1a(imax), x2a(jmax), x3a(kmax), & !- cell boundary position
                                               !- x1a(i)=x1(i+1/2)
             x1b(imax), x2b(jmax), x3b(kmax), & !- cell center position for MKS-!
             x1ab(imax), x2ab(jmax), x3ab(kmax) !- cell boundary position for MKS-! 
!
  real(8) :: sp2i, sm2i, spm2i, dk2fi, dk2bi, tmp2i, &
             sp2j, sm2j, spm2j, dk2fj, dk2bj, tmp2j, &
             sp2k, sm2k, spm2k, dk2fk, dk2bk, tmp2k, &
             cfi, cbi, cfj, cbj, cfk, cbk, &
             x1p, x1c, x1m, x1ap, x1am, &
             x2p, x2c, x2m, x2ap, x2am, &
             x3p, x3c, x3m, x3ap, x3am, &
             tmp1ia, tmp1ib, tmp1ic, tmp1id, tmp1i, &
             tmp1ja, tmp1jb, tmp1jc, tmp1jd, tmp1j, &
             tmp1ka, tmp1kb, tmp1kc, tmp1kd, tmp1k, &
             cadalfa, cadbeta, cadgamma
!
  real(8), parameter :: small=1.d-12
 
!
!-----------------------------------------------------------------------
!- Parameter

  cadalfa=0.5d0
  cadbeta=2.d0
  cadgamma=1.6d0
  
! 
  do k=ks1+1,ke1-1
    do j=js1+1,je1-1
      do i=is1+1,ie1-1
        do n=1,nv
!
!!! i-th direction
!           
          if(metric .eq. 403) then
            x1p=x1b(i+1)
            x1c=x1b(i)
            x1m=x1b(i-1)
            x1ap=x1ab(i)
            x1am=x1ab(i-1)
          else
            x1p=x1(i+1)
            x1c=x1(i)
            x1m=x1(i-1)
            x1ap=x1a(i)
            x1am=x1a(i-1)
          endif
!         
          tmp1ia=x1p-x1c+small
          tmp1ib=x1c-x1m+small
          tmp1ic=x1ap-x1c+small
          tmp1id=x1c-x1am+small 
!
!- S+, S-          
          sp2i=(uri(n,i+1,j,k)-uri(n,i,j,k))*(1.d0/tmp1ia)   
          sm2i=(uri(n,i,j,k)-uri(n,i-1,j,k))*(1.d0/tmp1ib)  

          spm2i=sp2i*sm2i          
!- CF, CB  
          cfi=tmp1ia*(1.d0/tmp1ic)
          cbi=tmp1ib*(1.d0/tmp1id)
!
!- limiter
          if(sp2i .lt. 0.d0) then
            tmp2i=-1.d0
          else
            tmp2i=1.d0
          endif   
          
          dk2bi=max(0.d0, min( (2.d0*abs(sp2i)+tmp2i*sm2i)*(1.d0/3.d0), &
                           max(-cadalfa*sm2i*tmp2i, &
                            min(cadbeta*sm2i*tmp2i, &
                                (2.d0*abs(sp2i)+tmp2i*sm2i)*(1.d0/3.d0), &
                                cadgamma*abs(sp2i)))))    

          dk2fi=max(0.d0, min( (abs(sp2i)+tmp2i*2.d0*sm2i)*(1.d0/3.d0), &
                           max(-cadalfa*abs(sp2i), &
                            min(cadbeta*abs(sp2i), &
                                (abs(sp2i)+tmp2i*2.d0*sm2i)*(1.d0/3.d0), &
                                cadgamma*tmp2i*sm2i)))) 
          
!           uriil(n,i,j,k)=uri(n,i,j,k)+dk2i*(0.5*(x1p-x1c))
!           uriir(n,i-1,j,k)=uri(n,i,j,k)-dk2i*(0.5*(x1c-x1p))
!
          uriil(n,i,j,k)=uri(n,i,j,k)+dk2fi*(x1ap-x1c)
          uriir(n,i-1,j,k)=uri(n,i,j,k)-dk2bi*(x1c-x1am)
!
!!! j-th direction
!
          if(metric .eq. 403) then
            x2p=x2b(j+1)
            x2c=x2b(j)
            x2m=x2b(j-1)
            x2ap=x2ab(j)
            x2am=x2ab(j-1)
          else
            x2p=x2(j+1)
            x2c=x2(j)
            x2m=x2(j-1)
            x2ap=x2a(j)
            x2am=x2a(j-1)
          endif  

          tmp1ja=x2p-x2c+small
          tmp1jb=x2c-x2m+small
          tmp1jc=x2ap-x2c+small
          tmp1jd=x2c-x2am+small          
!
!          sp2j=(uri(n,i,j+1,k)-uri(n,i,j,k))/dx2(j+1)
!          sm2j=(uri(n,i,j,k)-uri(n,i,j-1,k))/dx2(j)
!  
          sp2j=(uri(n,i,j+1,k)-uri(n,i,j,k))*(1.d0/tmp1ja) 
          sm2j=(uri(n,i,j,k)-uri(n,i,j-1,k))*(1.d0/tmp1jb)

          spm2j=sp2j*sm2j
!- CF, CB  
          cfj=tmp1ja*(1.d0/tmp1jc)  
          cbj=tmp1jb*(1.d0/tmp1jd)
!
!- limiter
          if(sp2j .lt. 0.d0) then
            tmp2j=-1.d0
          else
            tmp2j=1.d0
          endif   

          dk2bj=max(0.d0, min( (2.d0*abs(sp2j)+tmp2j*sm2j)*(1.d0/3.d0), &
                           max(-cadalfa*sm2j*tmp2j, &
                            min(cadbeta*sm2j*tmp2j, &
                                (2.d0*abs(sp2j)+tmp2j*sm2j)*(1.d0/3.d0), &
                                cadgamma*abs(sp2j)))))    

          dk2fj=max(0.d0, min( (abs(sp2j)+tmp2j*2.d0*sm2j)*(1.d0/3.d0), &
                           max(-cadalfa*abs(sp2j), &
                            min(cadbeta*abs(sp2j), &
                                (abs(sp2j)+tmp2j*2.d0*sm2j)*(1.d0/3.d0), &
                                cadgamma*tmp2j*sm2j)))) 

!          urijl(n,i,j,k)=uri(n,i,j,k)+dk2j*(0.5*(x2p-x2c))
!          urijr(n,i,j-1,k)=uri(n,i,j,k)-dk2j*(0.5*(x2c-x2m))
!
          urijl(n,i,j,k)=uri(n,i,j,k)+dk2fj*(x2ap-x2c)
          urijr(n,i,j-1,k)=uri(n,i,j,k)-dk2bj*(x2c-x2am)
!
!!! k-th direction          
!
          if(metric .eq. 403) then
            x3p=x3b(k+1)
            x3c=x3b(k)
            x3m=x3b(k-1)
            x3ap=x3ab(k)
            x3am=x3ab(k-1)
          else
            x3p=x3(k+1)
            x3c=x3(k)
            x3m=x3(k-1)
            x3ap=x3a(k)
            x3am=x3a(k-1)
          endif
!
          tmp1ka=x3p-x3c+small
          tmp1kb=x3c-x3m+small 
          tmp1kc=x3ap-x3c+small
          tmp1kd=x3c-x3am+small
!
!          sp2k=(uri(n,i,j,k+1)-uri(n,i,j,k))/dx3(k+1)
!          sm2k=(uri(n,i,j,k)-uri(n,i,j,k-1))/dx3(k)
!                   
          sp2k=(uri(n,i,j,k+1)-uri(n,i,j,k))*(1.d0/tmp1ka)
          sm2k=(uri(n,i,j,k)-uri(n,i,j,k-1))*(1.d0/tmp1kb)   

          spm2k=sp2k*sm2k
!- CF, CB   
          cfk=tmp1ka*(1.d0/tmp1kc)  
          cbk=tmp1kb*(1.d0/tmp1kd)

!- limiter
          if(sp2k .lt. 0.d0) then
            tmp2k=-1.d0
          else
            tmp2k=1.d0
          endif        

          dk2bk=max(0.d0, min( (2.d0*abs(sp2k)+tmp2k*sm2k)*(1.d0/3.d0), &
                           max(-cadalfa*sm2k*tmp2k, &
                            min(cadbeta*sm2k*tmp2k, &
                                (2.d0*abs(sp2k)+tmp2k*sm2k)*(1.d0/3.d0), &
                                cadgamma*abs(sp2k)))))    

          dk2fk=max(0.d0, min( (abs(sp2k)+tmp2k*2.d0*sm2k)*(1.d0/3.d0), &
                           max(-cadalfa*abs(sp2k), &
                            min(cadbeta*abs(sp2k), &
                                (abs(sp2k)+tmp2k*2.d0*sm2k)*(1.d0/3.d0), &
                                cadgamma*tmp2k*sm2k)))) 
          
!          urikl(n,i,j,k)=uri(n,i,j,k)+dk2k*(0.5*(x3p-x3c))
!          urikr(n,i,j,k-1)=uri(n,i,j,k)-dk2k*(0.5*(x3c-x3m))
!
          urikl(n,i,j,k)=uri(n,i,j,k)+dk2fk*(x3ap-x3c)
          urikr(n,i,j,k-1)=uri(n,i,j,k)-dk2bk*(x3c-x3am)
          
        enddo
      enddo
    enddo
  enddo

  return
end subroutine lim03a
!---------------------------------------------------------------------@
subroutine mvllim(uri,x1,x2,x3,x1a,x2a,x3a,x1b,x2b,x3b,x1ab,x2ab,x3ab,&
                  uriir,urijr,urikr,uriil,urijl,urikl,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Reconstruction step with modified van-Leer flux Limitter 
!     Calculate cell-interface variables(uriil, uriir) 
!     from cell-center variables (uri)
!
  use pram, only : imax, jmax, kmax, nv, metric
  implicit none
!
  integer :: i, j, k, n, is1, ie1, js1, je1, ks1, ke1
!
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: x1(imax), x2(jmax), x3(kmax), & !- cell center position
             x1a(imax), x2a(jmax), x3a(kmax), & !- cell boundary position
                                             !- x1a(i)=x1(i+1/2) 
             x1b(imax), x2b(jmax), x3b(kmax), & !- cell center position for MKS-!
             x1ab(imax), x2ab(jmax), x3ab(kmax) !- cell boundary position for MKS-! 
!
  real(8) :: sp2i, sm2i, spm2i, dk2i, tmp2i, &
             sp2j, sm2j, spm2j, dk2j, tmp2j, &
             sp2k, sm2k, spm2k, dk2k, tmp2k, &
             cfi, cbi, cfj, cbj, cfk, cbk, svi, svj, svk, &
             x1p, x1c, x1m, x1ap, x1am, &
             x2p, x2c, x2m, x2ap, x2am, &
             x3p, x3c, x3m, x3ap, x3am, &
             tmp1ia, tmp1ib, tmp1ic, tmp1id, tmp1i, &
             tmp1ja, tmp1jb, tmp1jc, tmp1jd, tmp1j, &
             tmp1ka, tmp1kb, tmp1kc, tmp1kd, tmp1k
!
  real(8), parameter :: small=1.d-12 
!-----------------------------------------------------------------------

  do k=ks1+1,ke1-1
    do j=js1+1,je1-1
      do i=is1+1,ie1-1
        do n=1,nv
!
!!! i-th direction
!           
          if(metric .eq. 403) then
            x1p=x1b(i+1)
            x1c=x1b(i)
            x1m=x1b(i-1)
            x1ap=x1ab(i)
            x1am=x1ab(i-1)
          else
            x1p=x1(i+1)
            x1c=x1(i)
            x1m=x1(i-1)
            x1ap=x1a(i)
            x1am=x1a(i-1)
          endif 
!
!          sp2i=(uri(n,i+1,j,k)-uri(n,i,j,k))/dx1(i+1)
!          sm2i=(uri(n,i,j,k)-uri(n,i-1,j,k))/dx1(i)
!- S+, S-
          tmp1ia=x1p-x1c+small
          tmp1ib=x1c-x1m+small
  
          sp2i=(uri(n,i+1,j,k)-uri(n,i,j,k))/tmp1ia   
          sm2i=(uri(n,i,j,k)-uri(n,i-1,j,k))/tmp1ib  

          spm2i=sp2i*sm2i
!- CF, CB
          tmp1ic=x1ap-x1c+small
          tmp1id=x1c-x1am+small 
   
          cfi=tmp1ia/tmp1ic   
          cbi=tmp1ib/tmp1id

!- check sign          
          if(spm2i .le. 0.d0) then
            dk2i=0.0
          else
!- limiter
!             dk2i=2.*max(spm2i,0.d0)/(sp2i+sm2i+small)
             
             tmp1i=sm2i*sm2i+(cfi+cbi-2.d0)*sm2i*sp2i+sp2i*sp2i
!  
             dk2i=(sm2i*sp2i*(cfi*sm2i+cbi*sp2i))/(tmp1i+small)
!            
          endif

!          uriil(n,i,j,k)=uri(n,i,j,k)+dk2i*(0.5*(x1p-x1c))
!          uriir(n,i-1,j,k)=uri(n,i,j,k)-dk2i*(0.5*(x1c-x1m))
!          
          uriil(n,i,j,k)=uri(n,i,j,k)+dk2i*tmp1ic
          uriir(n,i-1,j,k)=uri(n,i,j,k)-dk2i*tmp1id
!
!!! j-th direction
!          
          if(metric .eq. 403) then
            x2p=x2b(j+1)
            x2c=x2b(j)
            x2m=x2b(j-1)
            x2ap=x2ab(j)
            x2am=x2ab(j-1)
          else
            x2p=x2(j+1)
            x2c=x2(j)
            x2m=x2(j-1)
            x2ap=x2a(j)
            x2am=x2a(j-1)
          endif 
!
!          sp2j=(uri(n,i,j+1,k)-uri(n,i,j,k))/dx2(j+1)
!          sm2j=(uri(n,i,j,k)-uri(n,i,j-1,k))/dx2(j)
!
          tmp1ja=x2p-x2c+small
          tmp1jb=x2c-x2m+small  

          sp2j=(uri(n,i,j+1,k)-uri(n,i,j,k))/tmp1ja 
          sm2j=(uri(n,i,j,k)-uri(n,i,j-1,k))/tmp1jb

          spm2j=sp2j*sm2j
!- CF, CB
          tmp1jc=x2ap-x2c+small
          tmp1jd=x2c-x2am+small
  
          cfj=tmp1ja/tmp1jc
          cbj=tmp1jb/tmp1jd
   
!- check sign
          if(spm2j .le. 0.d0) then
            dk2j=0.d0
          else
!- limiter
!             dk2j=2.*max(spm2j,0.d0)/(sp2j+sm2j+small)
             
             tmp1j=sm2j*sm2j+(cfj+cbj-2.d0)*sm2j*sp2j+sp2j*sp2j
             
             dk2j=(sm2j*sp2j*(cfj*sm2j+cbj*sp2j))/(tmp1j+small)

          endif
          
!          urijl(n,i,j,k)=uri(n,i,j,k)+dk2j*(0.5*(x2p-x2c))
!          urijr(n,i,j-1,k)=uri(n,i,j,k)-dk2j*(0.5*(x2c-x2m))
!
          urijl(n,i,j,k)=uri(n,i,j,k)+dk2j*tmp1jc
          urijr(n,i,j-1,k)=uri(n,i,j,k)-dk2j*tmp1jd
!
!!! k-th direction
!          
          if(metric .eq. 403) then
            x3p=x3b(k+1)
            x3c=x3b(k)
            x3m=x3b(k-1)
            x3ap=x3ab(k)
            x3am=x3ab(k-1)
          else
            x3p=x3(k+1)
            x3c=x3(k)
            x3m=x3(k-1)
            x3ap=x3a(k)
            x3am=x3a(k-1)
          endif
!
!          sp2k=(uri(n,i,j,k+1)-uri(n,i,j,k))/dx3(k+1)
!          sm2k=(uri(n,i,j,k)-uri(n,i,j,k-1))/dx3(k)
!          
          tmp1ka=x3p-x3c+small
          tmp1kb=x3c-x3m+small       
 
          sp2k=(uri(n,i,j,k+1)-uri(n,i,j,k))/tmp1ka
          sm2k=(uri(n,i,j,k)-uri(n,i,j,k-1))/tmp1kb

          spm2k=sp2k*sm2k
!- CF, CB
          tmp1kc=x3ap-x3c+small
          tmp1kd=x3c-x3am+small
  
          cfk=tmp1ka/tmp1kc  
          cbk=tmp1kb/tmp1kd
   
!- check sign          
          if(spm2k .le. 0.d0) then
            dk2k=0.d0
          else
!- limiter           
!             dk2k=2.*max(spm2k,0.d0)/(sp2k+sm2k+small)
!
             tmp1k=sm2k*sm2k+(cfk+cbk-2.d0)*sm2k*sp2k+sp2k*sp2k
!  
             dk2k=(sm2k*sp2k*(cfk*sm2k+cbk*sp2k))/(tmp1k+small)

          endif

!          urikl(n,i,j,k)=uri(n,i,j,k)+dk2k*(0.5*(x3p-x3c))
!          urikr(n,i,j,k-1)=uri(n,i,j,k)-dk2k*(0.5*(x3c-x3m))
!          
          urikl(n,i,j,k)=uri(n,i,j,k)+dk2k*tmp1kc
          urikr(n,i,j,k-1)=uri(n,i,j,k)-dk2k*tmp1kd
!          
        enddo
      enddo
    enddo
  enddo
!
  return
end subroutine mvllim
!
!---------------------------------------------------------------------@
subroutine koren(uri,x1,x2,x3,x1a,x2a,x3a,x1b,x2b,x3b,x1ab,x2ab,x3ab,&
                 uriir,urijr,urikr,uriil,urijl,urikl,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!
!     Reconstruction step with Koren limiter 
!     Calculate cell-interface variables(uril, urir, wwl, wwr) 
!     from cell-center variables (uri)
!
  use pram, only : imax, jmax, kmax, nv, metric
  implicit none
!
  integer :: i, j, k, n, is1, ie1, js1, je1, ks1, ke1
!
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: x1(imax), x2(jmax), x3(kmax), & !- cell center position
             x1a(imax), x2a(jmax), x3a(kmax), & !- cell boundary position
                                               !- x1a(i)=x1(i+1/2)
             x1b(imax), x2b(jmax), x3b(kmax), & !- cell center position for MKS-!
             x1ab(imax), x2ab(jmax), x3ab(kmax) !- cell boundary position for MKS-! 
!
  real(8) :: sp2i, sm2i, spm2i, dk2fi, dk2bi, tmp2i, &
             sp2j, sm2j, spm2j, dk2fj, dk2bj, tmp2j, &
             sp2k, sm2k, spm2k, dk2fk, dk2bk, tmp2k, &
             cfi, cbi, cfj, cbj, cfk, cbk, &
             x1p, x1c, x1m, x1ap, x1am, &
             x2p, x2c, x2m, x2ap, x2am, &
             x3p, x3c, x3m, x3ap, x3am, &
             tmp1ia, tmp1ib, tmp1ic, tmp1id, tmp1i, &
             tmp1ja, tmp1jb, tmp1jc, tmp1jd, tmp1j, &
             tmp1ka, tmp1kb, tmp1kc, tmp1kd, tmp1k
!
  real(8), parameter :: small=1.d-12 
!
!-----------------------------------------------------------------------

  do k=ks1+1,ke1-1
    do j=js1+1,je1-1
      do i=is1+1,ie1-1
        do n=1,nv
!
!!! i-th direction
!           
          if(metric .eq. 403) then
            x1p=x1b(i+1)
            x1c=x1b(i)
            x1m=x1b(i-1)
            x1ap=x1ab(i)
            x1am=x1ab(i-1)
          else
            x1p=x1(i+1)
            x1c=x1(i)
            x1m=x1(i-1)
            x1ap=x1a(i)
            x1am=x1a(i-1)
          endif
!         
          tmp1ia=x1p-x1c+small
          tmp1ib=x1c-x1m+small
          tmp1ic=x1ap-x1c+small
          tmp1id=x1c-x1am+small 
!
!- S+, S-          
          sp2i=(uri(n,i+1,j,k)-uri(n,i,j,k))*(1.d0/tmp1ia)   
          sm2i=(uri(n,i,j,k)-uri(n,i-1,j,k))*(1.d0/tmp1ib)  

          spm2i=sp2i*sm2i          
!- CF, CB  
          cfi=tmp1ia*(1.d0/tmp1ic)
          cbi=tmp1ib*(1.d0/tmp1id)
!
!- limiter
          if(sp2i .lt. 0.d0) then
            tmp2i=-1.d0
          else
            tmp2i=1.d0
          endif   
          
!          dk2fi=tmp2i*max(0.d0, min(2.d0*abs(sp2i), tmp2i*2.d0*sm2i, &
!                          (2.d0*abs(sp2i)+tmp2i*sm2i)*(1.d0/3.d0) ))
!          dk2bi=tmp2i*max(0.d0, min(2.0*abs(sp2i), tmp2i*2.d0*sm2i, &
!                          (abs(sp2i)+2.d0*tmp2i*sm2i)*(1.d0/3.d0)  ))

          dk2fi=tmp2i*max(0.d0, min(abs(cfi*sp2i), tmp2i*cbi*sm2i, &
                          (2.d0*abs(sp2i)+tmp2i*sm2i)*(1.d0/3.d0) ))
          dk2bi=tmp2i*max(0.d0, min(abs(cfi*sp2i), tmp2i*cbi*sm2i, &
                          (abs(sp2i)+2.d0*tmp2i*sm2i)*(1.d0/3.d0)  ))
          
          
!           uriil(n,i,j,k)=uri(n,i,j,k)+dk2i*(0.5*(x1p-x1c))
!           uriir(n,i-1,j,k)=uri(n,i,j,k)-dk2i*(0.5*(x1c-x1p))
!
          uriil(n,i,j,k)=uri(n,i,j,k)+dk2fi*(x1ap-x1c)
          uriir(n,i-1,j,k)=uri(n,i,j,k)-dk2bi*(x1c-x1am)
!
!!! j-th direction
!
          if(metric .eq. 403) then
            x2p=x2b(j+1)
            x2c=x2b(j)
            x2m=x2b(j-1)
            x2ap=x2ab(j)
            x2am=x2ab(j-1)
          else
            x2p=x2(j+1)
            x2c=x2(j)
            x2m=x2(j-1)
            x2ap=x2a(j)
            x2am=x2a(j-1)
          endif  

          tmp1ja=x2p-x2c+small
          tmp1jb=x2c-x2m+small
          tmp1jc=x2ap-x2c+small
          tmp1jd=x2c-x2am+small          
!
!          sp2j=(uri(n,i,j+1,k)-uri(n,i,j,k))/dx2(j+1)
!          sm2j=(uri(n,i,j,k)-uri(n,i,j-1,k))/dx2(j)
!  
          sp2j=(uri(n,i,j+1,k)-uri(n,i,j,k))*(1.d0/tmp1ja) 
          sm2j=(uri(n,i,j,k)-uri(n,i,j-1,k))*(1.d0/tmp1jb)

          spm2j=sp2j*sm2j
!- CF, CB  
          cfj=tmp1ja*(1.d0/tmp1jc)  
          cbj=tmp1jb*(1.d0/tmp1jd)
!
!- limiter
          if(sp2j .lt. 0.d0) then
            tmp2j=-1.d0
          else
            tmp2j=1.d0
          endif   
          
!          dk2fj=tmp2j*max(0.d0, min(2.d0*abs(sp2j), tmp2j*2.d0*sm2j, &
!                          (2.d0*abs(sp2j)+tmp2j*sm2j)*(1.d0/3.d0) ))
!          dk2bj=tmp2j*max(0.d0, min(2.d0*abs(sp2j), tmp2j*2.d0*sm2j, &
!                          (abs(sp2j)+2.d0*tmp2i*sm2j)*(1.d0/3.d0)  ))

          dk2fj=tmp2j*max(0.d0, min(abs(cfj*sp2j), tmp2j*cbj*sm2j, &
                          (2.d0*abs(sp2j)+tmp2j*sm2j)*(1.d0/3.d0) ))
          dk2bj=tmp2j*max(0.d0, min(abs(cfj*sp2j), tmp2j*cbj*sm2j, &
                          (abs(sp2j)+2.d0*tmp2i*sm2j)*(1.d0/3.d0)  ))
          
!          urijl(n,i,j,k)=uri(n,i,j,k)+dk2j*(0.5*(x2p-x2c))
!          urijr(n,i,j-1,k)=uri(n,i,j,k)-dk2j*(0.5*(x2c-x2m))
!
          urijl(n,i,j,k)=uri(n,i,j,k)+dk2fj*(x2ap-x2c)
          urijr(n,i,j-1,k)=uri(n,i,j,k)-dk2bj*(x2c-x2am)
!
!!! k-th direction          
!
          if(metric .eq. 403) then
            x3p=x3b(k+1)
            x3c=x3b(k)
            x3m=x3b(k-1)
            x3ap=x3ab(k)
            x3am=x3ab(k-1)
          else
            x3p=x3(k+1)
            x3c=x3(k)
            x3m=x3(k-1)
            x3ap=x3a(k)
            x3am=x3a(k-1)
          endif
!
          tmp1ka=x3p-x3c+small
          tmp1kb=x3c-x3m+small 
          tmp1kc=x3ap-x3c+small
          tmp1kd=x3c-x3am+small
!
!          sp2k=(uri(n,i,j,k+1)-uri(n,i,j,k))/dx3(k+1)
!          sm2k=(uri(n,i,j,k)-uri(n,i,j,k-1))/dx3(k)
!                   
          sp2k=(uri(n,i,j,k+1)-uri(n,i,j,k))*(1.d0/tmp1ka)
          sm2k=(uri(n,i,j,k)-uri(n,i,j,k-1))*(1.d0/tmp1kb)   

          spm2k=sp2k*sm2k
!- CF, CB   
          cfk=tmp1ka*(1.d0/tmp1kc)  
          cbk=tmp1kb*(1.d0/tmp1kd)

!- limiter
          if(sp2k .lt. 0.d0) then
            tmp2k=-1.d0
          else
            tmp2k=1.d0
          endif        

!          dk2fk=tmp2k*max(0.d0, min(2.d0*abs(sp2k), tmp2k*2.d0*sm2k, &
!                          (2.d0*abs(sp2k)+tmp2k*sm2k)*(1.d0/3.d0) ))
!          dk2bk=tmp2k*max(0.d0, min(2.d0*abs(sp2k), tmp2k*2.d0*sm2k, &
!                          (abs(sp2k)+2.d0*tmp2k*sm2k)*(1.d0/3.d0)  ))

          dk2fk=tmp2k*max(0.d0, min(abs(cfk*sp2k), tmp2k*cbk*sm2k, &
                          (2.d0*abs(sp2k)+tmp2k*sm2k)*(1.d0/3.d0) ))
          dk2bk=tmp2k*max(0.d0, min(abs(cfk*sp2k), tmp2k*cbk*sm2k, &
                          (abs(sp2k)+2.d0*tmp2k*sm2k)*(1.d0/3.d0)  ))
          
!          urikl(n,i,j,k)=uri(n,i,j,k)+dk2k*(0.5*(x3p-x3c))
!          urikr(n,i,j,k-1)=uri(n,i,j,k)-dk2k*(0.5*(x3c-x3m))
!
          urikl(n,i,j,k)=uri(n,i,j,k)+dk2fk*(x3ap-x3c)
          urikr(n,i,j,k-1)=uri(n,i,j,k)-dk2bk*(x3c-x3am)
          
        enddo
      enddo
    enddo
  enddo

  return
end subroutine koren
!
