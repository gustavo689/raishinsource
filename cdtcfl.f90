!-----------------------------------------------------------------------
subroutine cdtcfl4(uri,gcov,gcon,dtcfl,dx1,dx2,dx3,dtx1,dtx2,dtx3,nm1, &
                   is1,ie1,js1,je1,ks1,ke1)
!-----------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv, gam, c0, icha  
  implicit none

  integer :: i, j, k, nm1, nmax, is1, ie1, js1, je1, ks1, ke1
  integer :: merr

  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)
!
  real(8) :: gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)
!      
  real(8) :: dx1(imax), dx2(jmax), dx3(kmax)
!
  real(8), allocatable :: dtx1a(:,:,:), dtx2a(:,:,:), dtx3a(:,:,:)
!
  real(8), allocatable :: vfip(:,:,:), vfim(:,:,:), &
           vfjp(:,:,:), vfjm(:,:,:), vfkp(:,:,:), vfkm(:,:,:)
!     
  real(8) :: pmin, dmax, vmin1, tmp1, vmin, dtx1, dtx2, dtx3, dtcfl
  real(8) :: vfi, vfj, vfk 
!
  real(8), parameter :: small=1.d-12 
!
!-----------------------------------------------------------------------
  allocate( dtx1a(is1:ie1,js1:je1,ks1:ke1), dtx2a(is1:ie1,js1:je1,ks1:ke1), &
            dtx3a(is1:ie1,js1:je1,ks1:ke1), &
            vfip(is1:ie1,js1:je1,ks1:ke1), vfim(is1:ie1,js1:je1,ks1:ke1), &
            vfjp(is1:ie1,js1:je1,ks1:ke1), vfjm(is1:ie1,js1:je1,ks1:ke1), &
            vfkp(is1:ie1,js1:je1,ks1:ke1), vfkm(is1:ie1,js1:je1,ks1:ke1), &
            stat=merr )

!-----------------------------------------------------------------------
!  
  nmax=nm1*2+(nm1-1)
!  nmax=nm1*2+2
      
  vmin1=uri(5,is1,js1,ks1)/uri(1,is1,js1,ks1)
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
        tmp1=uri(5,i,j,k)/uri(1,i,j,k)
        vmin1=max(vmin1,tmp1)
      enddo
    enddo
  enddo
!
  vmin=sqrt(gam*vmin1)
!
  dtcfl=1.d0
  dtx1=1.d0
  dtx2=1.d0
  dtx3=1.d0

  if(icha .eq. 0 .or. icha .eq. 1 .or. icha .eq. 2) then
    call calcha4(uri,gcov,gcon,&
                 vfip,vfim,vfjp,vfjm,vfkp,vfkm,is1,ie1,js1,je1,ks1,ke1)
  elseif(icha .eq. 3) then
    call calcha5(uri,gcov,gcon,&
                 vfip,vfim,vfjp,vfjm,vfkp,vfkm,is1,ie1,js1,je1,ks1,ke1)    
  endif
     
  do k=ks1+1,ke1-1
    do j=js1+1,je1-1
      do i=is1+1,ie1-1

        vfi=max(0.d0,abs(vfip(i,j,k)),abs(vfim(i,j,k)))
        vfj=max(0.d0,abs(vfjp(i,j,k)),abs(vfjm(i,j,k)))
        vfk=max(0.d0,abs(vfkp(i,j,k)),abs(vfkm(i,j,k)))

        if( imax.gt.nmax .and. vfi.gt.0.d0 ) then
          if(vfi .gt. 1.d0) then
            vfi=1.d0
          endif   
          dtx1a(i,j,k)=abs(dx1(i)/(vfi+small))
        else
!          dtx1a(i,j,k)=dtx1
          dtx1a(i,j,k)=1.d10
        endif
!
        if( jmax.gt.nmax .and. vfj.gt.0.d0 ) then
          if(vfj .gt. 1.d0) then
            vfj=1.d0
          endif   
          dtx2a(i,j,k)=abs(dx2(j)/(vfj+small))
        else          
!           dtx2a(i,j,k)=dtx2
!           dtx2a(i,j,k)=dtcfl
          dtx2a(i,j,k)=1.d10
        endif
        if ( kmax.gt.nmax .and. vfk.gt.0.d0 ) then
          if(vfk .gt. 1.d0) then
            vfk=1.d0
          endif   
          dtx3a(i,j,k)=abs(dx3(k)/(vfk+small))
        else
!          dtx3a(i,j,k)=dtx3
!          dtx3a(i,j,k)=dtcfl
          dtx3a(i,j,k)=1.d10 
        endif

      enddo
    enddo
  enddo
          
  do k=ks1+1,ke1-1
    do j=js1+1,je1-1
      do i=is1+1,ie1-1

        dtcfl=min(dtcfl,dtx1a(i,j,k))
        dtcfl=min(dtcfl,dtx2a(i,j,k))
        dtcfl=min(dtcfl,dtx3a(i,j,k))
!
        dtx1=min(dtx1,dtx1a(i,j,k))
        dtx2=min(dtx2,dtx2a(i,j,k))
        dtx3=min(dtx3,dtx3a(i,j,k))
!
      enddo
    enddo
  enddo

  deallocate( dtx1a,dtx2a,dtx3a,vfip,vfim,vfjp,vfjm,vfkp,vfkm,stat=merr )
!
  return
end  subroutine cdtcfl4
!
