!--------------------------------------------------------------------
subroutine calflx(uriir,urijr,urikr,uriil,urijl,urikl,&
                  uuir,uujr,uukr,uuil,uujl,uukl,&
                  wwir,wwjr,wwkr,wwil,wwjl,wwkl, &
                  gcovi,gcovj,gcovk,gconi,gconj,gconk, &
                  detgi,detgj,detgk,nm0,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos  
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1
  integer :: nm0
  integer :: merr
       
  real(8) :: uriir(nv,is1:ie1,js1:je1,ks1:ke1), &
             uriil(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urijl(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikr(nv,is1:ie1,js1:je1,ks1:ke1), &
             urikl(nv,is1:ie1,js1:je1,ks1:ke1)
!
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

  real(8) :: gconi(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gconj(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gconk(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detgi(is1:ie1,js1:je1,ks1:ke1), & 
             detgj(is1:ie1,js1:je1,ks1:ke1), & 
             detgk(is1:ie1,js1:je1,ks1:ke1)

  real(8), allocatable :: gcov1i(:,:), gcov1j(:,:), gcov1k(:,:)
  real(8), allocatable :: gcon1i(:,:), gcon1j(:,:), gcon1k(:,:)
  real(8), allocatable :: deltai(:,:), deltaj(:,:), deltak(:,:)
  real(8), allocatable :: utilir(:), utiljr(:), utilkr(:), &
                          utilil(:), utiljl(:), utilkl(:)
  real(8), allocatable :: vconir(:), vconjr(:), vconkr(:), &
                          vconil(:), vconjl(:), vconkl(:)
  real(8), allocatable :: bconir(:), bconjr(:), bconkr(:), &
                          bconil(:), bconjl(:), bconkl(:)
  real(8), allocatable :: ucovir(:), ucovjr(:), ucovkr(:), &
                          ucovil(:), ucovjl(:), ucovkl(:), &
                          uconir(:), uconjr(:), uconkr(:), &
                          uconil(:), uconjl(:), uconkl(:)
  real(8), allocatable :: bbcovir(:), bbcovjr(:), bbcovkr(:), &
                          bbcovil(:), bbcovjl(:), bbcovkl(:), &
                          bbconir(:), bbconjr(:), bbconkr(:), &
                          bbconil(:), bbconjl(:), bbconkl(:) 

  real(8) :: deir, dejr, dekr, deil, dejl, dekl
  real(8) :: prir, prjr, prkr, pril, prjl, prkl
  real(8) :: rohir, rohjr, rohkr, rohil, rohjl, rohkl 
  real(8) :: roeir, roejr, roekr, roeil, roejl, roekl  
  real(8) :: gflir, gfljr, gflkr, gflil, gfljl, gflkl
  real(8) :: bbsqir, bbsqjr, bbsqkr, bbsqil, bbsqjl, bbsqkl
!
!- allocate variables -!
  allocate( gcov1i(0:3,0:3), gcov1j(0:3,0:3), gcov1k(0:3,0:3), &
            gcon1i(0:3,0:3), gcon1j(0:3,0:3), gcon1k(0:3,0:3), &
            deltai(0:3,0:3), deltaj(0:3,0:3), deltak(0:3,0:3), stat=merr)
  allocate( utilir(1:3), utiljr(1:3), utilkr(1:3), &
            utilil(1:3), utiljl(1:3), utilkl(1:3), stat=merr)
  allocate( vconir(1:3), vconjr(1:3), vconkr(1:3), &
            vconil(1:3), vconjl(1:3), vconkl(1:3), stat=merr)   
  allocate( bconir(1:3), bconjr(1:3), bconkr(1:3), &
            bconil(1:3), bconjl(1:3), bconkl(1:3), stat=merr)
  allocate( ucovir(0:3), ucovjr(0:3), ucovkr(0:3), &
            ucovil(0:3), ucovjl(0:3), ucovkl(0:3), &
            uconir(0:3), uconjr(0:3), uconkr(0:3), &
            uconil(0:3), uconjl(0:3), uconkl(0:3), stat=merr)
  allocate( bbcovir(0:3), bbcovjr(0:3), bbcovkr(0:3), &
            bbcovil(0:3), bbcovjl(0:3), bbcovkl(0:3), &
            bbconir(0:3), bbconjr(0:3), bbconkr(0:3), &
            bbconil(0:3), bbconjl(0:3), bbconkl(0:3), stat=merr)   
!
!=====================================================================

  do k=ks1+nm0-1,ke1-nm0
    do j=js1+nm0-1,je1-nm0
      do i=is1+nm0-1,ie1-nm0
!
!- copy of primitive variables -!
        deir=uriir(1,i,j,k)
        dejr=urijr(1,i,j,k)
        dekr=urikr(1,i,j,k)
        deil=uriil(1,i,j,k)
        dejl=urijl(1,i,j,k)
        dekl=urikl(1,i,j,k)
!
        do m=1,3
          utilir(m)=uriir(m+1,i,j,k)
          utiljr(m)=urijr(m+1,i,j,k)
          utilkr(m)=urikr(m+1,i,j,k)
          utilil(m)=uriil(m+1,i,j,k)
          utiljl(m)=urijl(m+1,i,j,k)
          utilkl(m)=urikl(m+1,i,j,k)
!
!          vconir(m)=uriir(m+1,i,j,k)
!          vconjr(m)=urijr(m+1,i,j,k)
!          vconkr(m)=urikr(m+1,i,j,k)
!          vconil(m)=uriil(m+1,i,j,k)
!          vconjl(m)=urijl(m+1,i,j,k)
!          vconkl(m)=urikl(m+1,i,j,k)
        enddo
!
        prir=uriir(5,i,j,k)
        prjr=urijr(5,i,j,k)
        prkr=urikr(5,i,j,k)
        pril=uriil(5,i,j,k)
        prjl=urijl(5,i,j,k)
        prkl=urikl(5,i,j,k)
!
        do m=1,3
          bconir(m)=uriir(m+6,i,j,k)
          bconjr(m)=urijr(m+6,i,j,k)
          bconkr(m)=urikr(m+6,i,j,k)
          bconil(m)=uriil(m+6,i,j,k)
          bconjl(m)=urijl(m+6,i,j,k)
          bconkl(m)=urikl(m+6,i,j,k)
        enddo

        if(ieos .eq. 0 .or. ieos .eq. 3) then
          rohir=deir+(gam/(gam-1.0))*prir
          rohjr=dejr+(gam/(gam-1.0))*prjr
          rohkr=dekr+(gam/(gam-1.0))*prkr
          rohil=deil+(gam/(gam-1.0))*pril
          rohjl=dejl+(gam/(gam-1.0))*prjl
          rohkl=dekl+(gam/(gam-1.0))*prkl
        elseif(ieos .eq. 1) then
          rohir=(5./2.)*prir+sqrt((9./4.)*prir**2+deir**2)
          rohjr=(5./2.)*prjr+sqrt((9./4.)*prjr**2+dejr**2)
          rohkr=(5./2.)*prkr+sqrt((9./4.)*prkr**2+dekr**2)
          rohil=(5./2.)*pril+sqrt((9./4.)*pril**2+deil**2)
          rohjl=(5./2.)*prjl+sqrt((9./4.)*prjl**2+dejl**2)
          rohkl=(5./2.)*prkl+sqrt((9./4.)*prkl**2+dekl**2)
        elseif(ieos .eq. 2) then
          roeir=(3./2.)*(prir+((3.*prir**2)/(2.0*deir+ &
                sqrt(2.*prir**2+4.*deir**2)) ))
          roejr=(3./2.)*(prjr+((3.*prjr**2)/(2.0*dejr+ &
                sqrt(2.*prjr**2+4.*dejr**2)) ))
          roekr=(3./2.)*(prkr+((3.*prkr**2)/(2.0*dekr+ &
                sqrt(2.*prkr**2+4.*dekr**2)) ))
          roeil=(3./2.)*(pril+((3.*pril**2)/(2.0*deil+ &
                sqrt(2.*pril**2+4.*deil**2)) ))
          roejl=(3./2.)*(prjl+((3.*prjl**2)/(2.0*dejl+ &
                sqrt(2.*prjl**2+4.*dejl**2)) ))
          roekl=(3./2.)*(prkl+((3.*prkl**2)/(2.0*dekl+ &
                sqrt(2.*prkl**2+4.*dekl**2)) ))

          rohir=deir+roeir+prir
          rohjr=dejr+roejr+prjr
          rohkr=dekr+roekr+prkr
          rohir=deil+roeil+pril
          rohjr=deil+roejl+prjl
          rohkr=dekl+roekl+prkl
        endif
 
!- copy of metric terms -!
        do m=0,3
          do n=0,3
            gcov1i(m,n)=gcovi(m,n,i,j,k)
            gcov1j(m,n)=gcovj(m,n,i,j,k)
            gcov1k(m,n)=gcovk(m,n,i,j,k)
            gcon1i(m,n)=gconi(m,n,i,j,k)
            gcon1j(m,n)=gconj(m,n,i,j,k)
            gcon1k(m,n)=gconk(m,n,i,j,k)
          enddo
        enddo
!
!- cal of lorentz factor -!
        call calgfl(utilir,gflir,gcov1i) 
        call calgfl(utiljr,gfljr,gcov1j) 
        call calgfl(utilkr,gflkr,gcov1k) 
        call calgfl(utilil,gflil,gcov1i) 
        call calgfl(utiljl,gfljl,gcov1j) 
        call calgfl(utilkl,gflkl,gcov1k) 
!- cal of contravariant 4-velocity -!
        call calucon(utilir,uconir,gflir,gcon1i) 
        call calucon(utiljr,uconjr,gfljr,gcon1j)
        call calucon(utilkr,uconkr,gflkr,gcon1k)
        call calucon(utilil,uconil,gflil,gcon1i) 
        call calucon(utiljl,uconjl,gfljl,gcon1j)
        call calucon(utilkl,uconkl,gflkl,gcon1k)
!- cal of covariant 4-velocity -!
        call lower(uconir,ucovir,gcov1i)
        call lower(uconjr,ucovjr,gcov1j)
        call lower(uconkr,ucovkr,gcov1k)
        call lower(uconil,ucovil,gcov1i)
        call lower(uconjl,ucovjl,gcov1j)
        call lower(uconkl,ucovkl,gcov1k)
!- cal of contravariant 4-magnetic field -!
        call calbbcon(bconir,bbconir,ucovir,uconir) 
        call calbbcon(bconjr,bbconjr,ucovjr,uconjr) 
        call calbbcon(bconkr,bbconkr,ucovkr,uconkr) 
        call calbbcon(bconil,bbconil,ucovil,uconil) 
        call calbbcon(bconjl,bbconjl,ucovjl,uconjl) 
        call calbbcon(bconkl,bbconkl,ucovkl,uconkl) 
!- cal of covariant 4-magnetic field -!
        call lower(bbconir,bbcovir,gcov1i) 
        call lower(bbconjr,bbcovjr,gcov1j) 
        call lower(bbconkr,bbcovkr,gcov1k) 
        call lower(bbconil,bbcovil,gcov1i) 
        call lower(bbconjl,bbcovjl,gcov1j) 
        call lower(bbconkl,bbcovkl,gcov1k) 
!- cal of 4-magnetic field square -! 
        call calbbsq(bbconir,bbcovir,bbsqir)
        call calbbsq(bbconjr,bbcovjr,bbsqjr)
        call calbbsq(bbconkr,bbcovkr,bbsqkr)
        call calbbsq(bbconil,bbcovil,bbsqil)
        call calbbsq(bbconjl,bbcovjl,bbsqjl)
        call calbbsq(bbconkl,bbcovkl,bbsqkl)
!- cal of delta function -!
        call caldelta(deltai,gcov1i,gcon1i)
        call caldelta(deltaj,gcov1j,gcon1j)
        call caldelta(deltak,gcov1k,gcon1k)
!
!- Calculation of numerical flux at cell-boundary -!
!
        wwir(1,i,j,k)=deir*uconir(1)
        wwir(2,i,j,k)=(rohir+bbsqir)*uconir(1)*ucovir(1) &
                     +(prir+0.5*bbsqir)*deltai(1,1)-bbconir(1)*bbcovir(1)
        wwir(3,i,j,k)=(rohir+bbsqir)*uconir(1)*ucovir(2) &
                     +(prir+0.5*bbsqir)*deltai(1,2)-bbconir(1)*bbcovir(2)
        wwir(4,i,j,k)=(rohir+bbsqir)*uconir(1)*ucovir(3) &
                     +(prir+0.5*bbsqir)*deltai(1,3)-bbconir(1)*bbcovir(3)
!        wwir(5,i,j,k)=(rohir+bbsqir)*uconir(1)*ucovir(0) &
!                     +(prir+0.5*bbsqir)*deltai(1,0)-bbconir(1)*bbcovir(0) &
!                     -deir*uconir(1)   
        wwir(5,i,j,k)=(rohir+bbsqir)*uconir(1)*ucovir(0) &
                     +(prir+0.5*bbsqir)*deltai(1,0)-bbconir(1)*bbcovir(0) &
                     +deir*uconir(1) 
!        wwir(5,i,j,k)=(rohir+bbsqir)*uconir(1)*ucovir(0)-bbconir(1)*bbcovir(0)
        wwir(6,i,j,k)=deir*(prir/(deir**(gam)))*uconir(1) 
        wwir(7,i,j,k)=0.d0
        wwir(8,i,j,k)=(uconir(1)*(1./uconkr(0)))*bconir(2) &
                     -(uconir(2)*(1./uconkr(0)))*bconir(1)
        wwir(9,i,j,k)=(uconir(1)*(1./uconkr(0)))*bconir(3) &
                     -(uconir(3)*(1./uconkr(0)))*bconir(1)

        wwjr(1,i,j,k)=dejr*uconjr(2)
        wwjr(2,i,j,k)=(rohjr+bbsqjr)*uconjr(2)*ucovjr(1) &
                     +(prjr+0.5*bbsqjr)*deltaj(2,1)-bbconjr(2)*bbcovjr(1)
        wwjr(3,i,j,k)=(rohjr+bbsqjr)*uconjr(2)*ucovjr(2) &
                     +(prjr+0.5*bbsqjr)*deltaj(2,2)-bbconjr(2)*bbcovjr(2)
        wwjr(4,i,j,k)=(rohjr+bbsqjr)*uconjr(2)*ucovjr(3) &
                     +(prjr+0.5*bbsqjr)*deltaj(2,3)-bbconjr(2)*bbcovjr(3)
!        wwjr(5,i,j,k)=(rohjr+bbsqjr)*uconjr(2)*ucovjr(0) &
!                     +(prjr+0.5*bbsqjr)*deltaj(2,0)-bbconjr(2)*bbcovjr(0) &
!                     -dejr*uconjr(2)
        wwjr(5,i,j,k)=(rohjr+bbsqjr)*uconjr(2)*ucovjr(0) &
                     +(prjr+0.5*bbsqjr)*deltaj(2,0)-bbconjr(2)*bbcovjr(0) &
                     +dejr*uconjr(2)
!        wwjr(5,i,j,k)=(rohjr+bbsqjr)*uconjr(2)*ucovjr(0)-bbconjr(2)*bbcovjr(0)
        wwjr(6,i,j,k)=dejr*(prjr/(dejr**(gam)))*uconjr(2)  
        wwjr(7,i,j,k)=(uconjr(2)*(1./uconkr(0)))*bconjr(1) &
                     -(uconjr(1)*(1./uconkr(0)))*bconjr(2)
        wwjr(8,i,j,k)=0.d0
        wwjr(9,i,j,k)=(uconjr(2)*(1./uconkr(0)))*bconjr(3) &
                     -(uconjr(3)*(1./uconkr(0)))*bconjr(2)

        wwkr(1,i,j,k)=dekr*uconkr(3)
        wwkr(2,i,j,k)=(rohkr+bbsqkr)*uconkr(3)*ucovkr(1) &
                     +(prkr+0.5*bbsqkr)*deltak(3,1)-bbconkr(3)*bbcovkr(1)
        wwkr(3,i,j,k)=(rohkr+bbsqkr)*uconkr(3)*ucovkr(2) &
                     +(prkr+0.5*bbsqkr)*deltak(3,2)-bbconkr(3)*bbcovkr(2)
        wwkr(4,i,j,k)=(rohkr+bbsqkr)*uconkr(3)*ucovkr(3) &
                     +(prkr+0.5*bbsqkr)*deltak(3,3)-bbconkr(3)*bbcovkr(3)
!        wwkr(5,i,j,k)=(rohkr+bbsqkr)*uconkr(3)*ucovkr(0) &
!                     +(prkr+0.5*bbsqkr)*deltak(3,0)-bbconkr(3)*bbcovkr(0) &
!                     -dekr*uconkr(3)
        wwkr(5,i,j,k)=(rohkr+bbsqkr)*uconkr(3)*ucovkr(0) &
                     +(prkr+0.5*bbsqkr)*deltak(3,0)-bbconkr(3)*bbcovkr(0) &
                     +dekr*uconkr(3)
!        wwkr(5,i,j,k)=(rohkr+bbsqkr)*uconkr(3)*ucovkr(0)-bbconkr(3)*bbcovkr(0)
        wwkr(6,i,j,k)=dekr*(prkr/(dekr**(gam)))*uconkr(3)  
        wwkr(7,i,j,k)=(uconkr(3)*(1./uconkr(0)))*bconkr(1) &
                     -(uconkr(1)*(1./uconkr(0)))*bconkr(3)
        wwkr(8,i,j,k)=(uconkr(3)*(1./uconkr(0)))*bconkr(2) &
                     -(uconkr(2)*(1./uconkr(0)))*bconkr(3)
        wwkr(9,i,j,k)=0.d0
!
        wwil(1,i,j,k)=deil*uconil(1)
        wwil(2,i,j,k)=(rohil+bbsqil)*uconil(1)*ucovil(1) &
                     +(pril+0.5*bbsqil)*deltai(1,1)-bbconil(1)*bbcovil(1)
        wwil(3,i,j,k)=(rohil+bbsqil)*uconil(1)*ucovil(2) &
                     +(pril+0.5*bbsqil)*deltai(1,2)-bbconil(1)*bbcovil(2)
        wwil(4,i,j,k)=(rohil+bbsqil)*uconil(1)*ucovil(3) &
                     +(pril+0.5*bbsqil)*deltai(1,3)-bbconil(1)*bbcovil(3)
!        wwil(5,i,j,k)=(rohil+bbsqil)*uconil(1)*ucovil(0) &
!                     +(pril+0.5*bbsqil)*deltai(1,0)-bbconil(1)*bbcovil(0) &
!                     -deil*uconil(1) 
        wwil(5,i,j,k)=(rohil+bbsqil)*uconil(1)*ucovil(0) &
                     +(pril+0.5*bbsqil)*deltai(1,0)-bbconil(1)*bbcovil(0) &
                     +deil*uconil(1) 
!        wwil(5,i,j,k)=(rohil+bbsqil)*uconil(1)*ucovil(0)-bbconil(1)*bbcovil(0)
        wwil(6,i,j,k)=deil*(pril/(deil**(gam)))*uconil(1) 
        wwil(7,i,j,k)=0.d0
        wwil(8,i,j,k)=(uconil(1)*(1./uconil(0)))*bconil(2) &
                     -(uconil(2)*(1./uconil(0)))*bconil(1)
        wwil(9,i,j,k)=(uconil(1)*(1./uconil(0)))*bconil(3) &
                     -(uconil(3)*(1./uconil(0)))*bconil(1)

        wwjl(1,i,j,k)=dejl*uconjl(2)
        wwjl(2,i,j,k)=(rohjl+bbsqjl)*uconjl(2)*ucovjl(1) &
                     +(prjl+0.5*bbsqjl)*deltaj(2,1)-bbconjl(2)*bbcovjl(1)
        wwjl(3,i,j,k)=(rohjl+bbsqjl)*uconjl(2)*ucovjl(2) &
                     +(prjl+0.5*bbsqjl)*deltaj(2,2)-bbconjl(2)*bbcovjl(2)
        wwjl(4,i,j,k)=(rohjl+bbsqjl)*uconjl(2)*ucovjl(3) &
                     +(prjl+0.5*bbsqjl)*deltaj(2,3)-bbconjl(2)*bbcovjl(3)
!        wwjl(5,i,j,k)=(rohjl+bbsqjl)*uconjl(2)*ucovjl(0) &
!                     +(prjl+0.5*bbsqjl)*deltaj(2,0)-bbconjl(2)*bbcovjl(0) &
!                     -dejl*uconjl(2)
        wwjl(5,i,j,k)=(rohjl+bbsqjl)*uconjl(2)*ucovjl(0) &
                     +(prjl+0.5*bbsqjl)*deltaj(2,0)-bbconjl(2)*bbcovjl(0) &
                     +dejl*uconjl(2)
!        wwjl(5,i,j,k)=(rohjl+bbsqjl)*uconjl(2)*ucovjl(0)-bbconjl(2)*bbcovjl(0)
        wwjl(6,i,j,k)=dejl*(prjl/(dejl**(gam)))*uconjl(2)  
        wwjl(7,i,j,k)=(uconjl(2)*(1./uconjl(0)))*bconjl(1) &
                     -(uconjl(1)*(1./uconjl(0)))*bconjl(2)
        wwjl(8,i,j,k)=0.d0
        wwjl(9,i,j,k)=(uconjl(2)*(1./uconjl(0)))*bconjl(3) &
                     -(uconjl(3)*(1./uconjl(0)))*bconjl(2)

        wwkl(1,i,j,k)=dekl*uconkl(3)
        wwkl(2,i,j,k)=(rohkl+bbsqkl)*uconkl(3)*ucovkl(1) &
                     +(prkl+0.5*bbsqkl)*deltak(3,1)-bbconkl(3)*bbcovkl(1)
        wwkl(3,i,j,k)=(rohkl+bbsqkl)*uconkl(3)*ucovkl(2) &
                     +(prkl+0.5*bbsqkl)*deltak(3,2)-bbconkl(3)*bbcovkl(2)
        wwkl(4,i,j,k)=(rohkl+bbsqkl)*uconkl(3)*ucovkl(3) &
                     +(prkl+0.5*bbsqkl)*deltak(3,3)-bbconkl(3)*bbcovkl(3)
!        wwkl(5,i,j,k)=(rohkl+bbsqkl)*uconkl(3)*ucovkl(0) &
!                     -(prkl+0.5*bbsqkl)*deltak(3,0)-bbconkl(3)*bbcovkl(0) &
!                     -dekl*uconkl(3)
        wwkl(5,i,j,k)=(rohkl+bbsqkl)*uconkl(3)*ucovkl(0) &
                     -(prkl+0.5*bbsqkl)*deltak(3,0)-bbconkl(3)*bbcovkl(0) &
                     +dekl*uconkl(3)
!        wwkl(5,i,j,k)=(rohkl+bbsqkl)*uconkl(3)*ucovkl(0)-bbconkl(3)*bbcovkl(0)
        wwkl(6,i,j,k)=dekl*(prkl/(dekl**(gam)))*uconkl(3) 
        wwkl(7,i,j,k)=(uconkl(3)*(1./uconkl(0)))*bconkl(1) &
                     -(uconkl(1)*(1./uconkl(0)))*bconkl(3)
        wwkl(8,i,j,k)=(uconkl(3)*(1./uconkl(0)))*bconkl(2) &
                     -(uconkl(2)*(1./uconkl(0)))*bconkl(3)
        wwkl(9,i,j,k)=0.d0
!
        do n=1,9
!
          wwir(n,i,j,k)=detgi(i,j,k)*wwir(n,i,j,k)
          wwjr(n,i,j,k)=detgj(i,j,k)*wwjr(n,i,j,k)
          wwkr(n,i,j,k)=detgk(i,j,k)*wwkr(n,i,j,k)
!
          wwil(n,i,j,k)=detgi(i,j,k)*wwil(n,i,j,k)
          wwjl(n,i,j,k)=detgj(i,j,k)*wwjl(n,i,j,k)
          wwkl(n,i,j,k)=detgk(i,j,k)*wwkl(n,i,j,k)
!
        enddo   
!
      enddo
    enddo
  enddo
!
  deallocate( gcov1i, gcov1j, gcov1k, gcon1i, gcon1j, gcon1k, &
              deltai, deltaj, deltak, &
              utilir, utiljr, utilkr, utilil, utiljl, utilkl, &
              vconir, vconjr, vconkr, vconil, vconjl, vconkl, &   
              bconir, bconjr, bconkr, bconil, bconjl, bconkl, &
              ucovir, ucovjr, ucovkr, ucovil, ucovjl, ucovkl, &
              uconir, uconjr, uconkr, uconil, uconjl, uconkl, &
              bbcovir, bbcovjr, bbcovkr, bbcovil, bbcovjl, bbcovkl, &
              bbconir, bbconjr, bbconkr, bbconil, bbconjl, bbconkl,stat=merr)
!      
  return
end subroutine calflx
!
!--------------------------------------------------------------------
subroutine calflx3i(urii,uui,wwi,gcovi,gconi,detgi,nm0, &
                    is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos  
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1
  integer :: nm0
  integer :: merr
       
  real(8) :: urii(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uui(nv,is1:ie1,js1:je1,ks1:ke1)  
  real(8) :: wwi(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: gcovi(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gconi(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detgi(is1:ie1,js1:je1,ks1:ke1)

  real(8), allocatable :: gcov1i(:,:), gcon1i(:,:), deltai(:,:)
  real(8), allocatable :: utili(:), vconi(:), bconi(:)
  real(8), allocatable :: ucovi(:), uconi(:), bbcovi(:), bbconi(:)   

  real(8) :: dei, pri, rohi, roei, gfli, bbsqi, tmp1, tmp2
  real(8) :: alpha1i, dde
!      
!- allocate variables -!
  allocate( gcov1i(0:3,0:3),  gcon1i(0:3,0:3), deltai(0:3,0:3), stat=merr)
  allocate( utili(1:3), vconi(1:3), bconi(1:3), stat=merr)
  allocate( ucovi(0:3), uconi(0:3), bbcovi(0:3), bbconi(0:3), stat=merr)
!
!=====================================================================

  do k=ks1+nm0-1,ke1-nm0
    do j=js1+nm0-1,je1-nm0
      do i=is1+nm0-1,ie1-nm0
!
!- copy of primitive variables -!
        dei=urii(1,i,j,k)
        do m=1,3
          utili(m)=urii(m+1,i,j,k)
!          vconi(m)=urii(m+1,i,j,k)
        enddo
        pri=urii(5,i,j,k)
        do m=1,3
          bconi(m)=urii(m+6,i,j,k)
        enddo
       
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          rohi=dei+(gam/(gam-1.0))*pri
        elseif(ieos .eq. 1) then
          rohi=(5./2.)*pri+sqrt((9./4.)*pri**2+dei**2)
        elseif(ieos .eq. 2) then
          roei=(3./2.)*(pri+((3.*pri**2)/(2.0*dei+sqrt(2.*pri**2+4.*dei**2)) ))
          rohi=dei+roei+pri
        endif
!
!- copy of metric terms -!
        do m=0,3
          do n=0,3
            gcov1i(m,n)=gcovi(m,n,i,j,k)
            gcon1i(m,n)=gconi(m,n,i,j,k)
          enddo
        enddo
!
!- cal of lorentz factor -!
        call calgfl(utili,gfli,gcov1i) 
!- cal of contravariant 4-velocity -!
        call calucon(utili,uconi,gfli,gcon1i) 
!- cal of covariant 4-velocity -!
        call lower(uconi,ucovi,gcov1i)
!- cal of contravariant 4-magnetic field -!
        call calbbcon(bconi,bbconi,ucovi,uconi) 
!- cal of covariant 4-magnetic field -!
        call lower(bbconi,bbcovi,gcov1i) 
!- cal of 4-magnetic field square -! 
        call calbbsq(bbconi,bbcovi,bbsqi)
!- cal of delta function -!
        call caldelta(deltai,gcov1i,gcon1i)
!
!- Calculation of numerical flux at cell-boundary -!
!
        wwi(1,i,j,k)=dei*uconi(1)
        wwi(2,i,j,k)=(rohi+bbsqi)*uconi(1)*ucovi(1) &
                    +(pri+0.5*bbsqi)*deltai(1,1)-bbconi(1)*bbcovi(1)
        wwi(3,i,j,k)=(rohi+bbsqi)*uconi(1)*ucovi(2) &
                    +(pri+0.5*bbsqi)*deltai(1,2)-bbconi(1)*bbcovi(2)
        wwi(4,i,j,k)=(rohi+bbsqi)*uconi(1)*ucovi(3) &
                    +(pri+0.5*bbsqi)*deltai(1,3)-bbconi(1)*bbcovi(3)
!        wwi(5,i,j,k)=(rohi+bbsqi)*uconi(1)*ucovi(0) &
!                    +(pri+0.5*bbsqi)*deltai(1,0)-bbconi(1)*bbcovi(0) &
!                    -dei*uconi(1)
        wwi(5,i,j,k)=(rohi+bbsqi)*uconi(1)*ucovi(0) &
                    +(pri+0.5*bbsqi)*deltai(1,0)-bbconi(1)*bbcovi(0) &
                    +dei*uconi(1)
!        wwi(5,i,j,k)=(rohi+bbsqi)*uconi(1)*ucovi(0) &
!                    +(pri+0.5*bbsqi)*deltai(1,0)-bbconi(1)*bbcovi(0) 
        wwi(6,i,j,k)=dei*(pri/(dei**(gam)))*uconi(1) 

!        wwi(7,i,j,k)=0.d0                                  !wrong?
!        wwi(8,i,j,k)=uconi(1)*bbconi(2)-uconi(2)*bbconi(1)
!        wwi(9,i,j,k)=uconi(1)*bbconi(3)-uconi(3)*bbconi(1)
!
!        wwi(7,i,j,k)=0.d0
!        wwi(8,i,j,k)=bbconi(1)*uconi(2)-bbconi(2)*uconi(1)
!        wwi(9,i,j,k)=bbconi(1)*uconi(3)-bbconi(3)*uconi(1)
!
        wwi(7,i,j,k)=0.d0                                  !wrong?
        wwi(8,i,j,k)=(uconi(1)*(1./uconi(0)))*bconi(2) &
                    -(uconi(2)*(1./uconi(0)))*bconi(1)
        wwi(9,i,j,k)=(uconi(1)*(1./uconi(0)))*bconi(3) &
                    -(uconi(3)*(1./uconi(0)))*bconi(1)
!
!        wwi(7,i,j,k)=0.d0
!        wwi(8,i,j,k)=(uconi(2)*(1./uconi(0)))*bconi(1) &
!                    -(uconi(1)*(1./uconi(0)))*bconi(2)
!        wwi(9,i,j,k)=(uconi(3)*(1./uconi(0)))*bconi(1) &
!                    -(uconi(1)*(1./uconi(0)))*bconi(3)        
!
        tmp1=rohi*uconi(0)*gfli
        tmp2=tmp1*ucovi(1)
!
        do n=1,9
          wwi(n,i,j,k)=detgi(i,j,k)*wwi(n,i,j,k)
        enddo
!
        
        alpha1i=1./sqrt(-gcon1i(0,0))
        dde=dei*gfli*uconi(1)*(1./gfli)
!
      enddo
    enddo
  enddo
!      
  deallocate( gcov1i, gcon1i, deltai, utili, vconi, bconi, &
              ucovi, uconi, bbcovi, bbconi, stat=merr)
  return
end subroutine calflx3i
!
!--------------------------------------------------------------------
subroutine calflx3j(urij,uuj,wwj,gcovj,gconj,detgj,nm0, &
                    is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos  
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1
  integer :: nm0
  integer :: merr

  real(8) :: urij(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uuj(nv,is1:ie1,js1:je1,ks1:ke1)  
  real(8) :: wwj(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: gcovj(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gconj(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: detgj(is1:ie1,js1:je1,ks1:ke1)

  real(8), allocatable :: gcov1j(:,:), gcon1j(:,:), deltaj(:,:)
  real(8), allocatable :: utilj(:), vconj(:), bconj(:)
  real(8), allocatable :: ucovj(:), uconj(:), bbcovj(:), bbconj(:)  

  real(8) :: dej, prj, rohj, roej, gflj,  bbsqj
!
!- allocate variables -!
  allocate( gcov1j(0:3,0:3),  gcon1j(0:3,0:3), deltaj(0:3,0:3), stat=merr)
  allocate( utilj(1:3), vconj(1:3), bconj(1:3), stat=merr)
  allocate( ucovj(0:3), uconj(0:3), bbcovj(0:3), bbconj(0:3), stat=merr)
!      
!=====================================================================

  do k=ks1+nm0-1,ke1-nm0
    do j=js1+nm0-1,je1-nm0
      do i=is1+nm0-1,ie1-nm0
!
!- copy of primitive variables -!
        dej=urij(1,i,j,k)
        do m=1,3
          utilj(m)=urij(m+1,i,j,k)
!          vconi(m)=urii(m+1,i,j,k)
        enddo
        prj=urij(5,i,j,k)
        do m=1,3
          bconj(m)=urij(m+6,i,j,k)
        enddo

        if(ieos .eq. 0 .or. ieos .eq. 3) then
          rohj=dej+(gam/(gam-1.0))*prj
        elseif(ieos .eq. 1) then
          rohj=(5./2.)*prj+sqrt((9./4.)*prj**2+dej**2)
        elseif(ieos .eq. 2) then
          roej=(3./2.)*(prj+((3.*prj**2)/(2.0*dej+sqrt(2.*prj**2+4.*dej**2)) ))
          rohj=dej+roej+prj
        endif
!
!- copy of metric terms -!
        do m=0,3
          do n=0,3
            gcov1j(m,n)=gcovj(m,n,i,j,k)
            gcon1j(m,n)=gconj(m,n,i,j,k)
          enddo
        enddo
!
!- cal of lorentz factor -!
        call calgfl(utilj,gflj,gcov1j) 
!- cal of contravariant 4-velocity -!
        call calucon(utilj,uconj,gflj,gcon1j) 
!- cal of covariant 4-velocity -!
        call lower(uconj,ucovj,gcov1j)
!- cal of contravariant 4-magnetic field -!
        call calbbcon(bconj,bbconj,ucovj,uconj) 
!- cal of covariant 4-magnetic field -!
        call lower(bbconj,bbcovj,gcov1j) 
!- cal of 4-magnetic field square -! 
        call calbbsq(bbconj,bbcovj,bbsqj)
!- cal of delta function -!
        call caldelta(deltaj,gcov1j,gcon1j) 
!
!- Calculation of numerical flux at cell-boundary -!
!
        wwj(1,i,j,k)=dej*uconj(2)
        wwj(2,i,j,k)=(rohj+bbsqj)*uconj(2)*ucovj(1) &
                    +(prj+0.5*bbsqj)*deltaj(2,1)-bbconj(2)*bbcovj(1)
        wwj(3,i,j,k)=(rohj+bbsqj)*uconj(2)*ucovj(2) &
                    +(prj+0.5*bbsqj)*deltaj(2,2)-bbconj(2)*bbcovj(2)
        wwj(4,i,j,k)=(rohj+bbsqj)*uconj(2)*ucovj(3) &
                    +(prj+0.5*bbsqj)*deltaj(2,3)-bbconj(2)*bbcovj(3)
!        wwj(5,i,j,k)=(rohj+bbsqj)*uconj(2)*ucovj(0) &
!                    +(prj+0.5*bbsqj)*deltaj(2,0)-bbconj(2)*bbcovj(0) &
!                    -dej*uconj(2)
        wwj(5,i,j,k)=(rohj+bbsqj)*uconj(2)*ucovj(0) &
                    +(prj+0.5*bbsqj)*deltaj(2,0)-bbconj(2)*bbcovj(0) &
                    +dej*uconj(2)
!        wwj(5,i,j,k)=(rohj+bbsqj)*uconj(2)*ucovj(0) &
!                    +(prj+0.5*bbsqj)*deltaj(2,0)-bbconj(2)*bbcovj(0)
        wwj(6,i,j,k)=dej*(prj/(dej**(gam)))*uconj(2) 
        
!        wwj(7,i,j,k)=uconj(2)*bbconj(1)-uconj(1)*bbconj(2) !wrong?
!        wwj(8,i,j,k)=0.d0
!        wwj(9,i,j,k)=uconj(2)*bbconj(3)-uconj(3)*bbconj(2)
!
!        wwj(7,i,j,k)=bbconj(2)*uconj(1)-bbconj(1)*uconj(2)
!        wwj(8,i,j,k)=0.d0
!        wwj(9,i,j,k)=bbconj(2)*uconj(3)-bbconj(3)*uconj(2)
!
        wwj(7,i,j,k)=(uconj(2)*(1./uconj(0)))*bconj(1) &   !wrong?
                    -(uconj(1)*(1./uconj(0)))*bconj(2)
        wwj(8,i,j,k)=0.d0
        wwj(9,i,j,k)=(uconj(2)*(1./uconj(0)))*bconj(3) &
                    -(uconj(3)*(1./uconj(0)))*bconj(2)
!
!        wwj(7,i,j,k)=(uconj(1)*(1./uconj(0)))*bconj(2) &
!                    -(uconj(2)*(1./uconj(0)))*bconj(1)
!        wwj(8,i,j,k)=0.d0
!        wwj(9,i,j,k)=(uconj(3)*(1./uconj(0)))*bconj(2) &
!                    -(uconj(2)*(1./uconj(0)))*bconj(3)
!        
        do n=1,9
          wwj(n,i,j,k)=detgj(i,j,k)*wwj(n,i,j,k)
        enddo
!
      enddo
    enddo
  enddo
!
  deallocate( gcov1j, gcon1j, deltaj, utilj, vconj, bconj, &
              ucovj, uconj, bbcovj, bbconj, stat=merr)
!      
  return
end subroutine calflx3j
!
!--------------------------------------------------------------------
subroutine calflx3k(urik,uuk,wwk,gcovk,gconk,detgk,nm0, &
                    is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos  
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1
  integer :: nm0
  integer :: merr

  real(8) :: urik(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uuk(nv,is1:ie1,js1:je1,ks1:ke1)  
  real(8) :: wwk(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: gcovk(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gconk(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: detgk(is1:ie1,js1:je1,ks1:ke1)

  real(8), allocatable :: gcov1k(:,:), gcon1k(:,:), deltak(:,:)
  real(8), allocatable :: utilk(:), vconk(:), bconk(:)
  real(8), allocatable :: ucovk(:), uconk(:), bbcovk(:), bbconk(:)   

  real(8) :: dek, prk, rohk, roek, gflk, bbsqk
!      
!- allocate variables -!
  allocate( gcov1k(0:3,0:3),  gcon1k(0:3,0:3), deltak(0:3,0:3), stat=merr)
  allocate( utilk(1:3), vconk(1:3), bconk(1:3), stat=merr)
  allocate( ucovk(0:3), uconk(0:3), bbcovk(0:3), bbconk(0:3), stat=merr)
!
!=====================================================================

  do k=ks1+nm0-1,ke1-nm0
    do j=js1+nm0-1,je1-nm0
      do i=is1+nm0-1,ie1-nm0
!
!- copy of primitive variables -!
        dek=urik(1,i,j,k)
        do m=1,3
          utilk(m)=urik(m+1,i,j,k)
!          vconk(m)=urik(m+1,i,j,k)
        enddo
        prk=urik(5,i,j,k)
        do m=1,3
          bconk(m)=urik(m+6,i,j,k)
        enddo

        if(ieos .eq. 0 .or. ieos .eq. 3) then
          rohk=dek+(gam/(gam-1.0))*prk
        elseif(ieos .eq. 1) then
          rohk=(5./2.)*prk+sqrt((9./4.)*prk**2+dek**2)
        elseif(ieos .eq. 2) then
          roek=(3./2.)*(prk+((3.*prk**2)/(2.0*dek+sqrt(2.*prk**2+4.*dek**2)) ))
          rohk=dek+roek+prk
        endif

!
!- copy of metric terms -!
        do m=0,3
          do n=0,3
            gcov1k(m,n)=gcovk(m,n,i,j,k)
            gcon1k(m,n)=gconk(m,n,i,j,k)
          enddo
        enddo
!
!- cal of lorentz factor -!
        call calgfl(utilk,gflk,gcov1k) 
!- cal of contravariant 4-velocity -!
        call calucon(utilk,uconk,gflk,gcon1k) 
!- cal of covariant 4-velocity -!
        call lower(uconk,ucovk,gcov1k)
!- cal of contravariant 4-magnetic field -!
        call calbbcon(bconk,bbconk,ucovk,uconk) 
!- cal of covariant 4-magnetic field -!
        call lower(bbconk,bbcovk,gcov1k) 
!- cal of 4-magnetic field square -! 
        call calbbsq(bbconk,bbcovk,bbsqk)
!- cal of delta function -!
        call caldelta(deltak,gcov1k,gcon1k)
!
!- Calculation of numerical flux at cell-boundary -!
!
        wwk(1,i,j,k)=dek*uconk(3)
        wwk(2,i,j,k)=(rohk+bbsqk)*uconk(3)*ucovk(1) &
                    +(prk+0.5*bbsqk)*deltak(3,1)-bbconk(3)*bbcovk(1)
        wwk(3,i,j,k)=(rohk+bbsqk)*uconk(3)*ucovk(2) &
                    +(prk+0.5*bbsqk)*deltak(3,2)-bbconk(3)*bbcovk(2)
        wwk(4,i,j,k)=(rohk+bbsqk)*uconk(3)*ucovk(3) &
                    +(prk+0.5*bbsqk)*deltak(3,3)-bbconk(3)*bbcovk(3)
!        wwk(5,i,j,k)=(rohk+bbsqk)*uconk(3)*ucovk(0) &
!                    +(prk+0.5*bbsqk)*deltak(3,0)-bbconk(3)*bbcovk(0) &
!                    -dek*uconk(3)
        wwk(5,i,j,k)=(rohk+bbsqk)*uconk(3)*ucovk(0) &
                    +(prk+0.5*bbsqk)*deltak(3,0)-bbconk(3)*bbcovk(0) &
                    +dek*uconk(3)
!        wwk(5,i,j,k)=(rohk+bbsqk)*uconk(3)*ucovk(0) &
!                    +(prk+0.5*bbsqk)*deltak(3,0)-bbconk(3)*bbcovk(0) 
        wwk(6,i,j,k)=dek*(prk/(dek**(gam)))*uconk(3) 
        
!        wwk(7,i,j,k)=uconk(3)*bbconk(1)-uconk(1)*bbconk(3)  !wrong?
!        wwk(8,i,j,k)=uconk(3)*bbconk(2)-uconk(2)*bbconk(3)
!        wwk(9,i,j,k)=0.d0
!
!        wwk(7,i,j,k)=bbconk(3)*uconk(1)-bbconk(1)*uconk(3)
!        wwk(8,i,j,k)=bbconk(3)*uconk(2)-bbconk(2)*uconk(3)
!        wwk(9,i,j,k)=0.d0
!
        wwk(7,i,j,k)=(uconk(3)*(1./uconk(0)))*bconk(1) &    !wrong?
                    -(uconk(1)*(1./uconk(0)))*bconk(3)
        wwk(8,i,j,k)=(uconk(3)*(1./uconk(0)))*bconk(2) &
                    -(uconk(2)*(1./uconk(0)))*bconk(3)
        wwk(9,i,j,k)=0.d0
!
!        wwk(7,i,j,k)=(uconk(1)*(1./uconk(0)))*bconk(3) &
!                    -(uconk(3)*(1./uconk(0)))*bconk(1)
!        wwk(8,i,j,k)=(uconk(2)*(1./uconk(0)))*bconk(3) &
!                    -(uconk(3)*(1./uconk(0)))*bconk(2)
!        wwk(9,i,j,k)=0.d0
!        
        do n=1,9
          wwk(n,i,j,k)=detgk(i,j,k)*wwk(n,i,j,k)
        enddo
!
      enddo
    enddo
  enddo
!      
  deallocate( gcov1k, gcon1k, deltak, utilk, vconk, bconk, &
              ucovk, uconk, bbcovk, bbconk, stat=merr)
!
  return
end subroutine calflx3k
!
!--------------------------------------------------------------------
subroutine calwwo(uri,uu,wwo,gcov,gcon,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos  
  implicit none

  integer :: i, j, k, m, n, is1, ie1, js1, je1, ks1, ke1
  integer :: merr
       
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: wwo(3,nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8), allocatable :: gcov1(:,:), gcon1(:,:), delta(:,:)
  real(8), allocatable :: util(:), vcon(:), bcon(:)
  real(8), allocatable :: ucov(:), ucon(:), bbcov(:), bbcon(:)  

  real(8) :: de, pr, roh, roe, gfl, bbsq
!      
!- allocate variables -!
  allocate( gcov1(0:3,0:3),  gcon1(0:3,0:3), delta(0:3,0:3), stat=merr)
  allocate( util(1:3), vcon(1:3), bcon(1:3), stat=merr)
  allocate( ucov(0:3), ucon(0:3), bbcov(0:3), bbcon(0:3), stat=merr)
!
!=====================================================================

  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
!
!- copy of primitive variables -!
        de=uri(1,i,j,k)
        do m=1,3
          util(m)=uri(m+1,i,j,k)
!          vcon(m)=uri(m+1,i,j,k)
        enddo        
        pr=uri(5,i,j,k)
        do m=1,3
          bcon(m)=uri(m+6,i,j,k)
        enddo
        
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          roh=de+(gam/(gam-1.0))*pr
        elseif(ieos .eq. 1) then
          roh=(5./2.)*pr+sqrt((9./4.)*pr**2+de**2)
        elseif(ieos .eq. 2) then
          roe=(3./2.)*(pr+((3.*pr**2)/(2.0*de+sqrt(2.*pr**2+4.*de**2)) ))
          roh=de+roe+pr
        endif
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
!- cal of delta function -!
        call caldelta(delta,gcov1,gcon1)
!
!- Calculation of numerical flux at cell-center -!
!
        wwo(1,1,i,j,k)=de*ucon(1)
        wwo(2,1,i,j,k)=de*ucon(2)
        wwo(3,1,i,j,k)=de*ucon(3)

        wwo(1,2,i,j,k)=(roh+bbsq)*ucon(1)*ucov(1)+(pr+0.5*bbsq)*delta(1,1) &
                      -bbcon(1)*bbcov(1)
        wwo(1,3,i,j,k)=(roh+bbsq)*ucon(1)*ucov(2)+(pr+0.5*bbsq)*delta(1,2) &
                      -bbcon(1)*bbcov(2)
        wwo(1,4,i,j,k)=(roh+bbsq)*ucon(1)*ucov(3)+(pr+0.5*bbsq)*delta(1,3) &
                      -bbcon(1)*bbcov(3)

        wwo(2,2,i,j,k)=(roh+bbsq)*ucon(2)*ucov(1)+(pr+0.5*bbsq)*delta(2,1) &
                      -bbcon(2)*bbcov(1)
        wwo(2,3,i,j,k)=(roh+bbsq)*ucon(2)*ucov(2)+(pr+0.5*bbsq)*delta(2,2) &
                      -bbcon(2)*bbcov(2)
        wwo(2,4,i,j,k)=(roh+bbsq)*ucon(2)*ucov(3)+(pr+0.5*bbsq)*delta(2,3) &
                      -bbcon(2)*bbcov(3)

        wwo(3,2,i,j,k)=(roh+bbsq)*ucon(3)*ucov(1)+(pr+0.5*bbsq)*delta(3,1) &
                      -bbcon(3)*bbcov(1)
        wwo(3,3,i,j,k)=(roh+bbsq)*ucon(3)*ucov(2)+(pr+0.5*bbsq)*delta(3,2) &
                      -bbcon(3)*bbcov(2)
        wwo(3,4,i,j,k)=(roh+bbsq)*ucon(3)*ucov(3)+(pr+0.5*bbsq)*delta(3,3) &
                      -bbcon(3)*bbcov(3)

!        wwo(1,5,i,j,k)=(roh+bbsq)*ucon(1)*ucov(0)-bbcon(1)*bbcov(0) &
!                      -de*ucon(1)      
!        wwo(2,5,i,j,k)=(roh+bbsq)*ucon(2)*ucov(0)-bbcon(2)*bbcov(0) &
!                      -de*ucon(2) 
!        wwo(3,5,i,j,k)=(roh+bbsq)*ucon(3)*ucov(0)-bbcon(3)*bbcov(0) &
!                      -de*ucon(3)

        wwo(1,5,i,j,k)=(roh+bbsq)*ucon(1)*ucov(0)-bbcon(1)*bbcov(0) &
                      +de*ucon(1)      
        wwo(2,5,i,j,k)=(roh+bbsq)*ucon(2)*ucov(0)-bbcon(2)*bbcov(0) &
                      +de*ucon(2) 
        wwo(3,5,i,j,k)=(roh+bbsq)*ucon(3)*ucov(0)-bbcon(3)*bbcov(0) &
                      +de*ucon(3)

!        wwo(1,5,i,j,k)=(roh+bbsq)*ucon(1)*ucov(0)-bbcon(1)*bbcov(0)       
!        wwo(2,5,i,j,k)=(roh+bbsq)*ucon(2)*ucov(0)-bbcon(2)*bbcov(0) 
!        wwo(3,5,i,j,k)=(roh+bbsq)*ucon(3)*ucov(0)-bbcon(3)*bbcov(0)

        wwo(1,6,i,j,k)=de*(pr/(de**(gam)))*ucon(1) 
        wwo(2,6,i,j,k)=de*(pr/(de**(gam)))*ucon(2)  
        wwo(3,6,i,j,k)=de*(pr/(de**(gam)))*ucon(3)  
           
        wwo(1,7,i,j,k)=0.d0
        wwo(1,8,i,j,k)=(ucon(1)*(1./ucon(0)))*bcon(2) &
                      -(ucon(2)*(1./ucon(0)))*bcon(1)
        wwo(1,9,i,j,k)=(ucon(1)*(1./ucon(0)))*bcon(3) &
                      -(ucon(3)*(1./ucon(0)))*bcon(1)

        wwo(2,7,i,j,k)=(ucon(2)*(1./ucon(0)))*bcon(1) &
                      -(ucon(1)*(1./ucon(0)))*bcon(2)
        wwo(2,8,i,j,k)=0.d0
        wwo(2,9,i,j,k)=(ucon(2)*(1./ucon(0)))*bcon(3) &
                      -(ucon(3)*(1./ucon(0)))*bcon(2)

        wwo(3,7,i,j,k)=(ucon(3)*(1./ucon(0)))*bcon(1) &
                      -(ucon(1)*(1./ucon(0)))*bcon(3)
        wwo(3,8,i,j,k)=(ucon(3)*(1./ucon(0)))*bcon(2) &
                      -(ucon(2)*(1./ucon(0)))*bcon(3)
        wwo(3,9,i,j,k)=0.d0
!
      enddo
    enddo
  enddo
!      
  deallocate( gcov1, gcon1, delta, util, vcon, bcon, &
              ucov, ucon, bbcov, bbcon, stat=merr)
!      
  return
end subroutine calwwo
!
