!--------------------------------------------------------------------
subroutine ident(uu,un,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv
  implicit none

  integer :: i, j, k, m, is1, ie1, js1, je1, ks1, ke1
  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1), un(nv,is1:ie1,js1:je1,ks1:ke1)
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
        do m=1,nv
          un(m,i,j,k) = uu(m,i,j,k)
        enddo
      enddo
    enddo
  enddo
!
  return
end subroutine ident
!
!--------------------------------------------------------------------
subroutine ident2(uu,un,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv
  implicit none

  integer :: i, j, k, m, is1, ie1, js1, je1, ks1, ke1

  real(8) :: uu(nv,is1-1:ie1+2,js1-1:je1+2,ks1-1:ke1+2)
  real(8) :: un(nv,is1-1:ie1+2,js1-1:je1+2,ks1-1:ke1+2)

  do k=ks1-2,ke1+2
    do j=js1-2,je1+2
      do i=is1-2,ie1+2
        do m=1,nv
          un(m,i,j,k) = uu(m,i,j,k)
        enddo
      enddo
    enddo
  enddo
!
  return
end subroutine ident2
!
!--------------------------------------------------------------------
subroutine ident3(uu,un,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv
  implicit none

  integer :: i, j, k, m, is1, ie1, js1, je1, ks1, ke1

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: un(nv,is1-1:ie1+2,js1-1:je1+2,ks1-1:ke1+2)
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
        do m=1,nv
          un(m,i,j,k) = uu(m,i,j,k)
        enddo
      enddo
    enddo
  enddo
!
  return
end subroutine ident3
!
!--------------------------------------------------------------------
subroutine uextend(uu,is1,ie1,js1,je1,ks1,ke1,myranki,myrankj,myrankk)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv, iprocs, jprocs, kprocs
  implicit none

  integer :: i, j, k, m, is1, ie1, js1, je1, ks1, ke1
  integer :: myranki, myrankj, myrankk

  real(8) :: uu(nv,is1-1:ie1+2,js1-1:je1+2,ks1-1:ke1+2)
!
  do m=1,nv
!
    do k=ks1,ke1
      do i=is1,ie1
        if(myrankj .eq. 0) then
          uu(m,i,-1,k)=uu(m,i,1,k)
          uu(m,i,0,k)=uu(m,i,1,k)
        endif
        if(myrankj .eq. jprocs-1) then
          uu(m,i,jmax+1,k)=uu(m,i,jmax,k)
          uu(m,i,jmax+2,k)=uu(m,i,jmax,k)
        endif
      enddo
    enddo
!
    do k=ks1,ke1
      do j=js1,je1
        if(myranki .eq. 0) then
          uu(m,-1,j,k)=uu(m,1,j,k)
          uu(m,0,j,k)=uu(m,1,j,k)
        endif
        if(myranki .eq. iprocs-1) then
          uu(m,imax+1,j,k)=uu(m,imax,j,k)
          uu(m,imax+2,j,k)=uu(m,imax,j,k)
        endif
      enddo
    enddo
!
    do j=js1,je1
      do i=is1,ie1
        if(myrankk .eq. 0) then
          uu(m,i,j,-1)=uu(m,i,j,1)
          uu(m,i,j,0)=uu(m,i,j,1)
        endif
        if(myrankk .eq. kprocs-1) then
          uu(m,i,j,kmax+1)=uu(m,i,j,kmax)
          uu(m,i,j,kmax+2)=uu(m,i,j,kmax)
        endif
      enddo
    enddo
!
  enddo
!
  return
end subroutine uextend
!
!***********************************************************************
!           GENERAL FUNCTIONS
!***********************************************************************
!--------------------------------------------------------------------
function fjump(xx,bi,ee,aa,dd,bb)
!--------------------------------------------------------------------
  implicit none

  real(8) :: xx, bi, ee, aa, dd, bb, fjump
!
  if( xx .le. aa-bb ) then
    fjump=bi
  elseif( xx .ge. aa+bb ) then
    fjump=ee
  else
    fjump=1.0+tanh((xx-aa)/dd)/tanh(bb/dd)
    fjump=bi+0.5*fjump*(ee-bi)
  endif
!
  return
end function fjump
!--------------------------------------------------------------------
function ran0(idum)
!--------------------------------------------------------------------
  implicit none

  integer :: idum, k
  integer, parameter :: ia=16807
  integer, parameter :: im=2147483647
  integer, parameter :: iq=127773
  integer, parameter :: ir=2836
  integer, parameter :: mask=123459876
  real(8) :: ran0, am

  am=1.d0/float(im)
  idum=ieor(idum, mask)
  k=idum/iq
  idum=ia*(idum-k*iq)-ir*k
  if(idum .lt. 0) idum=idum+im
  ran0=am*float(idum)
  idum=ieor(idum, mask)
! write(*,*) 'ran0=',ran0    
  return  
end function ran0
!--------------------------------------------------------------------
function ran1(idum)
!--------------------------------------------------------------------
  implicit none

  integer :: idum,idum2,j,k
  real(8) :: ran1

  integer, parameter :: im1=2147483563
  integer, parameter :: im2=2147483399
  integer, parameter :: imm1=im1-1  
  integer, parameter :: ia1=40014
  integer, parameter :: ia2=40692
  integer, parameter :: iq1=53668
  integer, parameter :: iq2=52774 
  integer, parameter :: ir1=12211 
  integer, parameter :: ir2=3791
  integer, parameter :: ntab=32
  integer, parameter :: ndiv=1+imm1/ntab

  real(8), parameter :: am=1.d0/float(im1)  
  real(8), parameter :: eps=1.2d-7

  real(8), parameter :: rnmx=1.0-eps
  integer :: iy, iv(ntab)  
! 
  idum2=123456789
  iy=0
  do j=1,ntab
    iv(j)=ntab*0
  enddo
    
  if (idum.le.0) then 
    idum=max(-idum,1) 
    idum2=idum 
    do j=ntab+8,1,-1 
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1 
      if (idum.lt.0) idum=idum+im1 
      if (j.le.ntab) iv(j)=idum 
    enddo  
    iy=iv(1) 
  endif 
     
  k=idum/iq1
  idum=ia1*(idum-k*iq1)-k*ir1
  if (idum.lt.0) idum=idum+im1 
  k=idum2/iq2 
  idum2=ia2*(idum2-k*iq2)-k*ir2 
  if (idum2.lt.0) idum2=idum2+im2 
  j=1+iy/ndiv 
  iy=iv(j)-idum2 
  iv(j)=idum 
  if(iy.lt.1) iy=iy+imm1 

  ran1=min(am*float(iy),rnmx) 
! write(*,*) 'ran1=',ran1      
  return 
end function ran1 
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function ran2(IDUM)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  implicit none
  
  integer :: IFF,IY,IR(97),IDUM,J
  real(8) :: RAN2
  integer, parameter :: M=714025,IA=1366,IC=150889
  real(8), parameter :: RM=1.4005112D-6
      
  IFF=0
  IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
    IFF=1
    IDUM=MOD(IC-IDUM,M)
    DO J=1,97
      IDUM=MOD(IA*IDUM+IC,M)
      IR(J)=IDUM
    ENDDO
    IDUM=MOD(IA*IDUM+IC,M)
    IY=IDUM
  ENDIF
  J=1+(97*IY)/M
  IF(J.GT.97.OR.J.LT.1) THEN
   stop
  ELSE
    IY=IR(J)
    RAN2=float(IY)*RM
    IDUM=MOD(IA*IDUM+IC,M)
    IR(J)=IDUM
  ENDIF

! write(*,*) 'RAN2=',RAN2   
  RETURN
end function ran2
!
!--------------------------------------------------------------------
subroutine identw(ww,wn,is1,ie1,js1,je1,ks1,ke1)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, nv
  implicit none

  integer :: i, j, k, m, is1, ie1, js1, je1, ks1, ke1

  real(8) :: ww(3,nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: wn(3,nv,is1:ie1,js1:je1,ks1:ke1)
!
  do k=ks1,ke1
    do j=js1,je1
      do i=is1,ie1
        do m=1,nv
          wn(1,m,i,j,k) = ww(1,m,i,j,k)
          wn(2,m,i,j,k) = ww(2,m,i,j,k)
          wn(3,m,i,j,k) = ww(3,m,i,j,k)
        enddo
      enddo
    enddo
  enddo
!
  return
end subroutine identw
!
!***********************************************************************
!           SPECIAL FUNCTIONS
!***********************************************************************
!***********************************************************************
function fsg(rr,zz,r0,hh)
!***********************************************************************
  implicit none
 
  real(8) :: fsg, rr, zz, r0, hh, zet
!
  zet=zz/r0*tanh(zz/hh)
  fsg=sqrt((rr/r0)**2+(1.0+zet)**2)-(1.0+zet)
!
  return
end function fsg
!
!----------------------------------------------------------------------
subroutine kclock(ig)
!----------------------------------------------------------------------
  implicit none
!      real*8 cpu8
  real(4) :: cpu4
!      real*4 dim(2)
  integer :: ig
!
!      ig (sec)    :output
!
  call rclock(cpu4)
  ig=int(cpu4)
!
  return
end subroutine kclock
!----------------------------------------------------------------------
subroutine rclock(cpu4)
!----------------------------------------------------------------------
!     subroutine for the calcuration of CPU time
!     should chenge the appropriate function for your using computer

  implicit none
!  real*8 cpu8
  real(4) :: cpu4, dim(2), etime
!
!      cpu4 (sec)    :output
!
!
! default
!def       cpu4=0.0
!
! sx4
! sx5
! Altix
  cpu4=etime(dim)
           
! ibm in NCSA
!      cpu4=etime_(dim)
! sx3
!sx3       call clock(cpu8)
!sx3       cpu4=cpu8
!
! uxp
!uxp       call clock(ig,0,0)
!uxp       cpu4=ig
!
  return
end subroutine rclock
