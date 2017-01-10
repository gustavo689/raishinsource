!----------------------------------------------------------------------
subroutine grid(x1,x2,x3,dx1,dx2,dx3,dx1a,dx2a,dx3a, &
                dx1b,dx2b,dx3b,x1a,x2a,x3a,x1b,x2b,x3b,x1ab,x2ab,x3ab)
!----------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, xmin, xmax, ymin, ymax, zmin, zmax, &
                   ix1, ix2, ix3, R0, hslope, pi, metric
  implicit none
!
  integer :: i, j, k
  integer :: merr
!
  real(8) ::  x1(imax), x2(jmax), x3(kmax)
  real(8) :: x1a(imax), x2a(jmax), x3a(kmax)
  real(8) :: x1b(imax), x2b(jmax), x3b(kmax)
  real(8) :: x1ab(imax), x2ab(jmax), x3ab(kmax)
  real(8) :: dx1(imax), dx2(jmax), dx3(kmax)
!
  real(8) :: dx1b(imax-1), dx2b(jmax-1), dx3b(kmax-1)
  real(8) :: dx1a, dx2a, dx3a
  real(8) :: xx0, x1stat, x1end, dx1t, x1t, x3stat, x3end, dx3t, x3t
  real(8) :: aa0, aa1, aa2, R01, x1tt, p1, hth, p2
!
!- x1-direction
!
  if(ix1 .eq. 1) then !- uniform grid -!
    do  i=1,imax
      dx1(i)=(xmax-xmin)/float(imax-1)
      x1(i)=xmin+dx1(i)*float(i-1)
      x1b(i)=x1(i)
    enddo 
    dx1a=dx1(2)
  elseif(ix1 .eq. 2) then !- logarithmic spacing grid -!
!    R01=0.d0
    R01=R0
     
    xx0=R01
    x1stat=log(xmin-xx0)
    x1end=log(xmax-xx0)
    dx1t=(x1end-x1stat)*(1./float(imax-1))
    do i=1,imax
      x1t=x1stat+float(i-1)*dx1t
      x1(i)=exp(x1t)+xx0  
      x1b(i)=x1t
!      write(*,*) "i, x1=", i, x1(i)  
    enddo
!
    do i=2,imax
      if(metric .eq. 403) then !- for Modified KS metric -!
        dx1(i)=x1b(i)-x1b(i-1)
      else
        dx1(i)=x1(i)-x1(i-1)
      endif
    enddo
    dx1(1)=dx1(2)    
!    
  elseif(ix1 .eq. 3) then !- non-uniform spacing grid -!
    do i=1,imax
      dx1t=(xmax-xmin)/float(imax-1)          
      x1t=xmin+dx1t*float(i-1) 
      aa0=5.d0
      aa1=(xmax-xmin)/aa0
      aa2=atan(aa0)
      x1tt=(x1t-xmin)/(xmax-xmin)      
      x1(i)=xmin+aa1*tan(aa2*x1tt)
    enddo
!   
    do i=2,imax
      dx1(i)=x1(i)-x1(i-1)
    enddo
    dx1(1)=dx1(2)
!
  elseif(ix1 .eq. 4) then !- logarithmic spacing grid ver.2-!
    R01=0.3d0
!     
    xx0=R01
    x1stat=log10(xmin-xx0)
    x1end=log10(xmax-xx0)
!    dx1t=log10((xmax-xx0)/(xmin-xx0))/float(imax-1)
    dx1t=(x1end-x1stat)*(1./float(imax-1))
    do i=1,imax
      x1t=x1stat+float(i-1)*dx1t
      x1(i)=10.0**(x1t)+xx0  
      x1b(i)=x1t
!      write(*,*) "i, x1=", i, x1(i)  
    enddo
!
    do i=2,imax
      if(metric .eq. 403) then !- for Modified KS metric -!
        dx1(i)=x1b(i)-x1b(i-1)
      else
        dx1(i)=x1(i)-x1(i-1)
      endif
    enddo
    dx1(1)=dx1(2) 

  elseif(ix1 .eq. 5) then !- logarithmic spacing grid ver.3-!
    R01=0.3d0
    p1=2.d0 

    xx0=R01
    x1stat=log(xmin-xx0)**(1./p1)
    x1end=log(xmax-xx0)**(1./p1)
    dx1t=(x1end-x1stat)*(1./float(imax-1))
    write(*,*) "x1stat, x1end, dx1t=", x1stat, x1end, dx1t    
    do i=1,imax
      x1t=x1stat+float(i-1)*dx1t
      x1(i)=exp(x1t**p1)+xx0 
      x1b(i)=x1t
!      write(*,*) "x1stat, x1end=", x1stat, x1end
!      write(*,*) "i, x1t, x1=", i, x1t, x1(i)  
    enddo
!
    do i=2,imax
      if(metric .eq. 403) then !- for Modified KS metric -!
        dx1(i)=x1b(i)-x1b(i-1)
      else
        dx1(i)=x1(i)-x1(i-1)
      endif
    enddo
    dx1(1)=dx1(2) 

  endif
!!!!
!- Cell-boundary position 
  do i=1,imax-1
    x1a(i)=0.5d0*(x1(i)+x1(i+1))
  enddo
  x1a(imax)=x1a(imax-1)+dx1(imax)
!
!- Cell boundary position for MKS  
  if(metric .eq. 403) then !- for Modified KS metric -!
    do i=1,imax-1
      x1ab(i)=0.5d0*(x1b(i)+x1b(i+1))
    enddo
    x1ab(imax)=x1ab(imax-1)+dx1(imax)
  endif
!
!- distance between cell-boundary position   
  do i=2,imax
    if(metric .eq. 403) then !- for Modified KS metric -!
      dx1b(i)=x1ab(i)-x1ab(i-1)
    else
      dx1b(i)=x1a(i)-x1a(i-1)
    endif
  enddo
  dx1b(1)=dx1b(2)
!
!- x2-direction
!
  if(ix2 .eq. 1) then !- uniform grid -!
    do j=1,jmax
      dx2(j)=(ymax-ymin)/float(jmax-1)
      x2(j)=ymin+dx2(j)*float(j-1)
      x2b(j)=x2(j)
!      write(*,*) "j, phi=", j, x2(j)*(180./pi)  
    enddo 
    dx2a=dx2(2)
  elseif(ix2 .eq. 2) then

  endif
!
!- Cell-boundary position  
  do j=1,jmax-1
     x2a(j)=0.5d0*(x2(j)+x2(j+1))
  enddo
    x2a(jmax)=x2a(jmax-1)+dx2(jmax)
!
!- Cell boudanry position for MKS
  do j=1,jmax
    x2ab(j)=x2a(j)
  enddo
!  
!- distance between cell-boundary position
  do j=2,jmax
    dx2b(j)=x2a(j)-x2a(j-1)
  enddo
  dx2b(1)=dx2b(2)
!
!- x3-direction
!
  if(ix3 .eq. 1) then !- uniform grid -!
    do k=1,kmax
      dx3(k)=(zmax-zmin)/float(kmax-1)
      x3(k)=zmin+dx3(k)*float(k-1)
      x3b(k)=x3(k)
    enddo
    dx3a=dx3(2)
!
  elseif(ix3 .eq. 2) then !- non-uniform grid -!
! 
!    hslope=0.3d0    
!    x3stat=(zmin)/pi
    x3stat=zmin/pi
    x3end=zmax/pi
!    dx3t=1./float(kmax)
!    dx3t=((zmax/pi)-(zmin/pi))/float(kmax-1)
    dx3t=(x3end-x3stat)/float(kmax-1)
    do k=1, kmax
      x3t=pi*(x3stat+float(k-1)*dx3t)
!      x3(k)=x3t+0.5*(1.-hslope)*sin(2.*x3t)
      x3(k)=x3t+0.5*(1.-hslope)*sin(2.*x3t)
      x3b(k)=x3t
!         write(*,*) "k, x3b, x3=", k, x3b(k)/pi, x3(k)/pi

    enddo
    do k=2,kmax
      if(metric .eq. 403) then !- for Modified KS metric -!
        dx3(k)=x3b(k)-x3b(k-1)
      else
        dx3(k)=x3(k)-x3(k-1)
      endif
    enddo
   dx3(1)=dx3(2)  
!
  elseif(ix3 .eq. 3) then !- non-uniform grid -!
! 
    hth=0.15
    p2=5.0
    
    x3stat=zmin/pi
    x3end=zmax/pi 
    dx3t=(x3end-x3stat)/float(kmax-1)
    do k=1, kmax
      x3t=pi*(x3stat+float(k-1)*dx3t)
      x3(k)=(0.5*pi)*(hth*(2.*x3t-1.)+(1.-hth)*((2.*x3t-1.)**p2)+1.)
      x3b(k)=x3t
         write(*,*) "k, x3b, x3=", k, x3(k)
    enddo    
    do k=2,kmax
      if(metric .eq. 403) then !- for Modified KS metric -!
        dx3(k)=x3b(k)-x3b(k-1)
      else
        dx3(k)=x3(k)-x3(k-1)
      endif
    enddo
    dx3(1)=dx3(2)  

  endif
!!!
!- Cell-boundary position
  do k=1,kmax-1
    x3a(k)=0.5d0*(x3(k)+x3(k+1))
  enddo
  x3a(kmax)=x3a(kmax-1)+dx3(kmax)
!
!- Cell boundary position for MKS
  if(metric .eq. 403) then !- for Modified KS metric -!    
    do k=1,kmax-1
      x3ab(k)=0.5d0*(x3b(k)+x3b(k+1))
    enddo
    x3ab(kmax)=x3ab(kmax-1)+dx3(kmax)     
  endif
!
!- distance between cell-bounary position
  do k=2,kmax
    if(metric .eq. 403) then !- for Modified KS metric -!  
      dx3b(k)=x3ab(k)-x3ab(k-1)
    else
      dx3b(k)=x3a(k)-x3a(k-1)
    endif
  enddo
  dx3b(1)=dx3b(2)

!
  return
end subroutine grid
