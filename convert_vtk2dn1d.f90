!***********************************************************************
!      RAISHIN code: 3D General Relativistic MHD (3D MPI) version 
!      program for output to VTK format
!      written by Y. Mizuno
!      ver. 150715
!      new version
!      using "uria" in main program (output1)
!***********************************************************************
!
program convert_vtk2dn1d
  use pram, only: imax, jmax, kmax, iprocs, jprocs, kprocs, akm, metric, &
                  gam, ieos

  implicit none

!======================================================================@
!    Difinition for variables
!======================================================================@
!
  integer, parameter :: nv=9
  integer, parameter :: npe=iprocs*jprocs*kprocs ! number of cpus 
  integer, parameter :: ns=0, ne=200 ! start and end data file number
  integer, parameter :: dataformat=0 !- 0= ascii, 1=binary -!
  integer, parameter :: idirec=1 ! direction (1=x, 2=y, 3=z)

!  integer, parameter :: nq=13
!
  integer :: i, j, k, l, m, n
  integer :: is, ie, js, je, ks, ke
  integer :: myrank, nd, npe1, myrank1, merr
  integer :: jh, kh, jh1, kh1, mmax, nmax
!  
  character*256 :: filename, filename1 
!
  real(8) :: uri(nv,imax,jmax,kmax)
  real(8) :: detg(imax,jmax,kmax)
  real(8) :: x1(imax), x2(jmax), x3(kmax)
  real(8) :: rr, th, ph, sth, cth, del, sig, aa, de, pr, roh
  real(8) :: alpha1, gfl1a
  real(8) :: util(1:3), gcov1(0:3,0:3), gcon1(0:3,0:3), ucon(0:3), ucov(0:3), &
             bcon(1:3), bbcon(0:3), bbcov(0:3), utiln(1:3)  
!
  real(8), allocatable :: uri1(:,:,:,:), detg1(:,:,:), x1a(:), x2a(:), x3a(:)
  real(8) :: time
  real(8), parameter :: bmin1=1.d-12 
  
!  real(8) :: qq(nq,imax,jmax,kmax)
!
  real(4) :: xx1(imax,ne+1), xx2(imax,ne+1), xx3(imax,ne+1), &
             dd(imax,ne+1), pp(imax,ne+1), ink(imax,ne+1), &
             v1(imax,ne+1), v2(imax,ne+1), v3(imax,ne+1), &
             b1(imax,ne+1), b2(imax,ne+1), b3(imax,ne+1)
!
!
!======================================================================@
!     File Open
!======================================================================@
!
!  open(unit=6,file='convert.outdat',status='unknown')
!  open(unit=8,file='structr.outdat',status='unknown',form='unformatted')
!
    write(filename1,991) 1
    open( unit=9, file=filename1,status='replace',form="formatted")
991 format('ok',i3.3,'.vtk')   
!
!======================================================================@
!
! Data initialize
!
!  do k=1,kmax
!    do j=1,jmax
!      do i=1,imax
!        do n=1,nv
!          uri(n,i,j,k)=1.d0
!        enddo
!        detg(i,j,k)=0.d0
!      enddo
!    enddo
!  enddo
!
!- Parameter
!  
!- set half grid position -!
  jh1=(jmax+1)/2
  kh1=(kmax+1)/2
!- set j-th direction =1 -!
  mmax=1
!
!======================================================================@
!
! Start main loop
!
  do nd=ns,ne
    write(*,*) 'working on ',nd,'th data'
! 
! read each process data
!
    do myrank=0,npe-1
      write(*,*) 'reading ', myrank, 'th process data'

      write(filename,990) myrank, nd
990   format('structr',i3.3,'-',i3.3,'.outdat')
      open( unit=7,file=filename,form='unformatted',status='old')

      read(7) time, myrank1, npe1
!
      if(myrank .ne. myrank1) then
        write(*,*) 'reading process data is wrong'
         stop
      endif
      if(npe .ne. npe1) then
        write(*,*) 'total process number is wrong'
        stop
      endif   
!
      read(7) is,ie,js,je,ks,ke
      allocate(uri1(nv,0:ie-is,0:je-js,0:ke-ks), &
               detg1(0:ie-is,0:je-js,0:ke-ks), &
               x1a(0:ie-is), x2a(0:je-js), x3a(0:ke-ks),stat=merr)

      do i=0,ie-is
        do j=0,je-js
          do k=0,ke-ks
            read(7) x1a(i),x2a(j),x3a(k),uri1(1,i,j,k),uri1(2,i,j,k), &
                   uri1(3,i,j,k),uri1(4,i,j,k),uri1(5,i,j,k),uri1(6,i,j,k), &
                   uri1(7,i,j,k),uri1(8,i,j,k),uri1(9,i,j,k), detg1(i,j,k) 
          enddo
        enddo
      enddo
!
!- calculation of Lorentz facor -!
      do i=0,ie-is
        do j=0,je-js
          do k=0,ke-ks            
!
            rr=x1a(i)
            ph=x2a(j)
            th=x3a(k)

            do n=0,3
              do m=0,3
                gcov1(m,n)=0.d0
                gcon1(m,n)=0.d0
              enddo
            enddo  
!
!            sth=sin(th)
            sth=1.d0
!            cth=cos(th)
            cth=0.d0

!-for Kerr BH -!
            if(metric .eq. 3) then
              del=rr*rr
              sig=rr*rr
              aa=(rr*rr)**2

              gcov1(0,0)=-1.0
              gcov1(1,1)=1.0
              gcov1(2,2)=(aa*(1./sig))*sth**2
              gcov1(3,3)=sig
!
              gcon1(0,0)=-1.0
              gcon1(1,1)=1.0
              gcon1(2,2)=1.0/(del*sth**2)
              gcon1(3,3)=1./sig   

            elseif(metric .eq. 203) then
              del=rr*rr-(2.*rr)+akm*akm
              sig=rr*rr+akm*cth*akm*cth
              aa=((rr*rr+akm*akm)**2)-akm*akm*del*sth*sth


              gcov1(0,0)=-(sig-2.*rr)*(1./sig)     
              gcov1(1,1)=sig*(1./del)               
              gcov1(2,2)=(aa*(1./sig))*sth**2 
              gcov1(0,2)=-(2.*akm*rr*(sth**2))*(1./sig) 
              gcov1(2,0)=gcov1(0,2)              
              gcov1(3,3)=sig                  
!
              gcon1(0,0)=-aa*(1./(sig*del))
              gcon1(1,1)=del*(1./sig)                          
              gcon1(2,2)=(sig-2.*rr)*(1./(del*sig*sth*sth)) 
              gcon1(0,2)=-(2.*akm*rr)*(1./(sig*del))           
              gcon1(2,0)=gcon1(0,2)   
              gcon1(3,3)=1./sig
            endif  
!
            de=uri1(1,i,j,k)            
            util(1)=uri1(2,i,j,k)
            util(2)=uri1(3,i,j,k)
            util(3)=uri1(4,i,j,k)
            pr=uri1(5,i,j,k)
!
!- trans 4-velocity to 3-velocity -!
!            
!        
            call calgfl(util,gfl1a,gcov1)
            uri1(2,i,j,k)=util(1)*(1./gfl1a)
            uri1(3,i,j,k)=util(2)*(1./gfl1a)
            uri1(4,i,j,k)=util(3)*(1./gfl1a)
            
         enddo
        enddo            
     enddo
!     
!- convert data on full array -!       
      do i=0,ie-is
        x1(is+i)=x1a(i)
      enddo
!
      do j=0,je-js
        x2(js+j)=x2a(j)
      enddo
!
      do k=0,ke-ks
        x3(ks+k)=x3a(k)
      enddo
!
      do k=0,ke-ks
        do j=0,je-js
          do i=0,ie-is
            uri(1,is+i,js+j,ks+k)=uri1(1,i,j,k)
            uri(2,is+i,js+j,ks+k)=uri1(2,i,j,k)
            uri(3,is+i,js+j,ks+k)=uri1(3,i,j,k)
            uri(4,is+i,js+j,ks+k)=uri1(4,i,j,k)
            uri(5,is+i,js+j,ks+k)=uri1(5,i,j,k)
            uri(6,is+i,js+j,ks+k)=uri1(6,i,j,k)
            uri(7,is+i,js+j,ks+k)=uri1(7,i,j,k)
            uri(8,is+i,js+j,ks+k)=uri1(8,i,j,k)
            uri(9,is+i,js+j,ks+k)=uri1(9,i,j,k)
          enddo
        enddo
      enddo    
!
      close(7)
      deallocate(uri1,detg1,x1a,x2a,x3a,stat=merr)
    enddo
!
!    do i=1,imax
!      xx1(i,nd+1)=x1(i)
!      xx2(i,nd+1)=time
!      xx3(i,nd+1)=0.d0
!      dd(i,nd+1)=uri(1,i,jh1,kh1)
!    enddo   
!
!  enddo



!
!----------------------------------------------------------------------@
!- Output data for analysis -!
!
!- set output file -!
!    write(filename1,991) nd
!    open( unit=9, file=filename1,status='replace',form="formatted")
!991 format('ok',i3.3,'.vtk')  
!
    if(idirec .eq. 1) then !- xz plane -!
!
!- convert the data array to 1D -!
!    
      do i=1, imax
!    
        xx1(i,nd+1)=x1(i)
        xx2(i,nd+1)=time
        xx3(i,nd+1)=0.d0
!            
        dd(i,nd+1)=uri(1,i,jh1,kh1) !- density -!            
        v1(i,nd+1)=uri(2,i,jh1,kh1)
                 !- 1st component of velocity -!
        v2(i,nd+1)=uri(3,i,jh1,kh1)
                 !- 2nd component of velocity -!
        v3(i,nd+1)=uri(4,i,jh1,kh1)
        pp(i,nd+1)=uri(5,i,jh1,kh1) !- gas pressure -!
        b1(i,nd+1)=uri(7,i,jh1,kh1)
                 !- 1st component of B-field -!
        b2(i,nd+1)=uri(8,i,jh1,kh1)
                 !- 2nd component of B-field -!
        b3(i,nd+1)=uri(9,i,jh1,kh1)
                 !- 3rd component of B-field -!
      enddo 
!       
    endif
!    
  enddo    
!
!  do i=1,imax
!    write(*,*) 'xx1,xx2,xx3=',xx1(i,200),xx2(i,200),xx3(i,200)
!  enddo
 
  call write_vtk2d(imax,ne+1,mmax,xx1,xx2,xx3,dd,v1,v2,v3,pp,&
                   b1,b2,b3,dataformat,filename)    
! 
  stop    
end program convert_vtk2dn1d

!------------------------------------------------------------
subroutine write_vtk2d(xmax,ymax,zmax,xx,yy,zz,dd,vx,vy,vz,pp,&
                       bx,by,bz,dataformat,filename)
!------------------------------------------------------------
!   
!- write data for VTK format -!
!  The gird points does not need to be regular.!
!
  implicit none
!
  integer :: i, j, k, xmax, ymax, zmax
  integer :: merr
  integer :: dataformat
  character*256 :: filename
!  
  real(4) :: xx(xmax,ymax), yy(xmax,ymax), zz(xmax,ymax)
  real(4) :: dd(xmax,ymax), pp(xmax,ymax)  
  real(4) :: vx(xmax,ymax), vy(xmax,ymax), vz(xmax,ymax) 
  real(4) :: bx(xmax,ymax), by(xmax,ymax), bz(xmax,ymax)
!  
  real(4) :: dum_r
  character(LEN=1), parameter :: newline = achar(10) ! newline symbol!
  
!  There are five basic parts to the VTK file format.  
!  1. Write file version and identifier (Header)
!     
    write(9, "(a)") "# vtk DataFile Version 2.0"    
!
!  2. Tilte
!    
    write(9, "(a)") "RMHD Simulation Result"    
!
!  3. Data type (ASCII or BINARY)
!    
    if(dataformat .eq. 0) then
      write(9, "(a)") "ASCII" 
    else
      write(9, "(a)") "BINARY"
    endif
!
!  4.  Dataset structure 
!
10  format(a11, 3(1x, i6))
11  format(a, 2x, i9, 2x, a5)
    write(9, "(a)") "DATASET STRUCTURED_GRID"
    write(9, 10)    "DIMENSIONS ", xmax, ymax, zmax
    write(9, 11)    "POINTS", xmax*ymax*zmax, "float"
!
    write(6,*) 'after header'
!
! -- now dumps coordnates
    if(dataformat .ne. 0) then 
      close(9)
! reopen it as unformatted and append
      open(unit=9, file=filename, status="old", form="unformatted", &
           position="append")
!           position="append", access="stream")
    endif
! -- now dumps data
    do j = 1, ymax
      do i = 1, xmax
        if(dataformat .eq. 0) then                 
          write(9,*) xx(i,j), yy(i,j), zz(i,j)
        else
          write(9) xx(i,j), yy(i,j), zz(i,j)
        endif
      enddo
    enddo
!
    write(6,*) 'after coordinates'
!    
!  5 . Data
!    
    if(dataformat .ne. 0) then
      close(9)
! reopen it as formatted and append
      open(unit=9, file=filename, status="old", form="formatted", &
           position="append")          
    endif
!
12  format(a, a11, 1x, i12)
! You need to add newline character for safty after writing data in binary
!   write(9, "(a)") 
    write(9, 12) newline, "POINT_DATA ", xmax*ymax*zmax
!
    write(6,*) 'after new line'   
!
! Write density
!
    write(9, "(a)") "SCALARS density float"
    write(9, "(a)") "LOOKUP_TABLE default"
!
! --- now dumps data
    if(dataformat .ne. 0) then
      close(9) 
! reopen it as unformatted and append
      open(unit=9, file=filename, status="old", form="unformatted", &
           position="append")   
!           position="append", access="stream")  
    endif
    if(dataformat .eq. 0) then
      do j = 1, ymax
        do i = 1, xmax             
          dum_r = max(dd(i,j), 1e-30) ! for safety
          write(9,*) dum_r
        enddo
      enddo
    else
      do j = 1, ymax
        do i = 1, xmax
          dum_r = max(dd(i,j), 1e-30) ! for safety
          write(9) dum_r
        enddo
      enddo
    endif
!
    write(6,*) 'after density'  
!
! Write velocity vector
!
!    if(dataformat .ne. 0) then
!      close(9) 
! reopen it as formatted and append
!      open(unit=9, file=filename, status="old", form="formatted", &
!           position="append")   
!    endif
!    
!    write(9, "(a, a)") newline, "VECTORS velocity float"
!    
!    if(dataformat .ne. 0) then
!      close(9) 
! reopen it as unformatted and append
!      open(unit=9, file=filename, status="old", form="unformatted", &
!           position="append")   
!           position="append", access="stream")   
!    endif
! -- now dumps data
!    if(dataformat .eq. 0) then
!      do j = 1, ymax
!        do i = 1, xmax
!          write(9,*) vx(i,j), vy(i,j), vz(i,j)
!        enddo           
!      enddo
!    else
!      do j = 1, ymax
!        do i = 1, xmax
!          write(9) vx(i,j), vy(i,j), vz(i,j)
!        enddo
!      enddo           
!    endif
!
!    write(6,*) 'after velocity vector'  
!
! Write B-field scalars (each components, Bx)
!
    write(9, "(a)") "SCALARS Vx float"
    write(9, "(a)") "LOOKUP_TABLE default"
!
! --- now dumps data
    if(dataformat .ne. 0) then
      close(9) 
! reopen it as unformatted and append
      open(unit=9, file=filename, status="old", form="unformatted", &
           position="append")   
!           position="append", access="stream")  
    endif
    if(dataformat .eq. 0) then
      do j = 1, ymax
        do i = 1, xmax             
          dum_r = vx(i,j) ! for safety
          write(9,*) dum_r
        enddo
      enddo
    else
      do j = 1, ymax
        do i = 1, xmax
          dum_r = vx(i,j) ! for safety
          write(9) dum_r
        enddo
      enddo
    endif
!
    write(6,*) 'after V-field scalar (Vx)'  
!
! Write pressure
!
    write(9, "(a)") "SCALARS pressure float"
    write(9, "(a)") "LOOKUP_TABLE default"
!
! --- now dumps data
    if(dataformat .ne. 0) then
      close(9) 
! reopen it as unformatted and append
      open(unit=9, file=filename, status="old", form="unformatted", &
           position="append")   
!           position="append", access="stream")  
    endif
    if(dataformat .eq. 0) then
      do j = 1, ymax
        do i = 1, xmax             
          dum_r = max(pp(i,j), 1e-30) ! for safety
          write(9,*) dum_r
        enddo
      enddo
    else
      do j = 1, ymax
        do i = 1, xmax
          dum_r = max(pp(i,j), 1e-30) ! for safety
          write(9) dum_r
        enddo
      enddo
    endif
!
    write(6,*) 'after pressure'  
!
! Write Lorentz factor
!
!    write(9, "(a)") "SCALARS LorentzW float"
!    write(9, "(a)") "LOOKUP_TABLE default"
!
! --- now dumps data
!    if(dataformat .ne. 0) then
!      close(9) 
! reopen it as unformatted and append
!      open(unit=9, file=filename, status="old", form="unformatted", &
!           position="append")   
!!           position="append", access="stream")  
!    endif
!    if(dataformat .eq. 0) then
!      do j = 1, ymax
!        do i = 1, xmax             
!          dum_r = max(gf(i,j), 1e-30) ! for safety
!          write(9,*) dum_r
!        enddo
!      enddo
!    else
!      do j = 1, ymax
!        do i = 1, xmax
!          dum_r = max(gf(i,j), 1e-30) ! for safety
!          write(9) dum_r
!        enddo
!      enddo
!    endif
!
!    write(6,*) 'after Lorentz factor'  
!
!
! Write B-field vector
!
!    if(dataformat .ne. 0) then
!      close(9) 
! reopen it as formatted and append
!      open(unit=9, file=filename, status="old", form="formatted", &
!           position="append")   
!    endif
!    
!    write(9, "(a, a)") newline, "VECTORS B-field float"
!    
!    if(dataformat .ne. 0) then
!      close(9) 
! reopen it as unformatted and append
!      open(unit=9, file=filename, status="old", form="unformatted", &
!           position="append")   
! !          position="append", access="stream") 
!    endif
! -- now dumps data
!    if(dataformat .eq. 0) then
!      do j = 1, ymax
!        do i = 1, xmax
!          write(9,*) bx(i,j), by(i,j), bz(i,j)
!        enddo           
!      enddo
!    else
!      do j = 1, ymax
!        do i = 1, xmax
!          write(9) bx(i,j), by(i,j), bz(i,j)
!        enddo
!      enddo           
!    endif
!
!    write(6,*) 'after B-field vector'  
!
! Write B-field scalars (each components, Bx)
!
!    write(9, "(a)") "SCALARS Bx float"
!    write(9, "(a)") "LOOKUP_TABLE default"
!
! --- now dumps data
!    if(dataformat .ne. 0) then
!      close(9) 
! reopen it as unformatted and append
!      open(unit=9, file=filename, status="old", form="unformatted", &
!           position="append")   
!!           position="append", access="stream")  
!    endif
!    if(dataformat .eq. 0) then
!      do j = 1, ymax
!        do i = 1, xmax             
!          dum_r = bx(i,j) ! for safety
!          write(9,*) dum_r
!        enddo
!      enddo
!    else
!      do j = 1, ymax
!        do i = 1, xmax
!          dum_r = bx(i,j) ! for safety
!          write(9) dum_r
!        enddo
!      enddo
!    endif
!
!    write(6,*) 'after B-field scalar (Bx): mh=', mh  
!
! Write B-field scalars (each components, By)
!
!    write(9, "(a)") "SCALARS By float"
!    write(9, "(a)") "LOOKUP_TABLE default"
!
! --- now dumps data
!    if(dataformat .ne. 0) then
!      close(9) 
! reopen it as unformatted and append
!      open(unit=9, file=filename, status="old", form="unformatted", &
!           position="append")   
!!           position="append", access="stream")  
!    endif
!    if(dataformat .eq. 0) then
!      do j = 1, ymax
!        do i = 1, xmax             
!          dum_r = by(i,j) ! for safety
!          write(9,*) dum_r
!        enddo
!      enddo
!    else
!      do j = 1, ymax
!        do i = 1, xmax
!          dum_r = by(i,j) ! for safety
!          write(9) dum_r
!        enddo
!      enddo
!    endif   
!
!    write(6,*) 'after B-field scalar (By)'  
!
! Write B-field scalars (each components, Bz)
!
!    write(9, "(a)") "SCALARS Bz float"
!    write(9, "(a)") "LOOKUP_TABLE default"
!
! --- now dumps data
!    if(dataformat .ne. 0) then
!      close(9) 
! reopen it as unformatted and append
!      open(unit=9, file=filename, status="old", form="unformatted", &
!           position="append")   
!!           position="append", access="stream")  
!    endif
!    if(dataformat .eq. 0) then
!      do j = 1, ymax
!        do i = 1, xmax             
!          dum_r = bz(i,j) ! for safety
!          write(9,*) dum_r
!        enddo
!      enddo
!    else
!      do j = 1, ymax
!        do i = 1, xmax
!          dum_r = bz(i,j) ! for safety
!          write(9) dum_r
!        enddo
!      enddo
!    endif
!
!    write(6,*) 'after B-field scalar (Bz)'  
!
! Write b^2
!
!    write(9, "(a)") "SCALARS b^2 float"
!    write(9, "(a)") "LOOKUP_TABLE default"
!
! --- now dumps data
!    if(dataformat .ne. 0) then
!      close(9) 
! reopen it as unformatted and append
!      open(unit=9, file=filename, status="old", form="unformatted", &
!           position="append")   
!!           position="append", access="stream")  
!    endif
!    if(dataformat .eq. 0) then
!      do j = 1, ymax
!        do i = 1, xmax             
!          dum_r = max(bbsq(i,j), 1e-30) ! for safety
!          write(9,*) dum_r
!        enddo
!      enddo
!    else
!      do j = 1, ymax
!        do i = 1, xmax
!          dum_r = max(bbsq(i,j), 1e-30) ! for safety
!          write(9) dum_r
!        enddo
!      enddo
!    endif
!
!    write(6,*) 'after b^2 scalar'  
!
    close(9)
!
  return
end subroutine write_vtk2d
!
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
  if(abs(utsq) .gt. 1.0d10) then
    utsq=1.0d9
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
subroutine transks2bl1(util,utiln,x1aa,x2aa,x3aa)
!--------------------------------------------------------------------
!- convert 4-velocity in KS to 4-velocity in BL coordinates -!
!-  reading whole primitive variables -!
!
!- util: \tilde{u}_KS (input)-!
!- uria: \tilde{u}_BL (output)-! 
!
!  
  use pram, only : metric, akm, ix1, ix2, ix3, R0, hslope, pi  
  implicit none

  integer :: m, n, merr
       
  real(8) :: util(1:3), utiln(1:3)

  real(8), allocatable :: ucon(:), ucon1(:)
  real(8), allocatable :: gcov1(:,:), gcon1(:,:), gcon1a(:,:), beta1(:)
  real(8), allocatable :: tmp1(:), trans(:,:)

  real(8) :: x1aa, x2aa, x3aa, rr, ph, th, rbh1, alpha1, gfl, tmp2, &
             sth, cth, del, sig, aa
!
!- allocate variables -!
  allocate(ucon(0:3), ucon1(0:3), gcov1(0:3,0:3), gcon1(0:3,0:3), &
           gcon1a(0:3,0:3), beta1(1:3), tmp1(0:3), trans(0:3,0:3), stat=merr)
!
!=====================================================================
!
   rbh1=1.+sqrt(1.-akm*akm)

   rr=x1aa
   ph=x2aa
   th=x3aa   
!
   if(rr .gt. rbh1) then
!
!- get KS metric -!

!           
     do n=0,3
       do m=0,3
         gcov1(m,n)=0.d0
         gcon1(m,n)=0.d0
       enddo
     enddo  
!      
     sth=sin(th)
     cth=cos(th)
     del=rr*rr-(2.*rr)+akm*akm
     sig=rr*rr+akm*cth*akm*cth
     aa=((rr*rr+akm*akm)**2)-akm*akm*del*sth*sth

     gcov1(0,0)=-1.+2.*rr/sig     
     gcov1(0,1)=2.*rr/sig          
     gcov1(0,2)=-2.*akm*rr*sth*sth/sig 
     gcov1(1,0)=gcov1(0,1)        
     gcov1(1,1)=1.+2.*rr/sig      
     gcov1(1,2)=-1.*akm*sth*sth*(1.+2.*rr/sig) 
     gcov1(2,0)=gcov1(0,2)        
     gcov1(2,1)=gcov1(1,2)        
     gcov1(2,2)=aa*sth*sth/sig 

     gcon1(0,0)=-1.-2.*rr/sig
     gcon1(0,1)=2.*rr/sig
     gcon1(1,0)=gcon1(0,1)
     gcon1(1,1)=del/sig
     gcon1(1,2)=akm/sig
     gcon1(2,1)=gcon1(1,2)
     gcon1(2,2)=1./(sig*sth*sth)
     gcon1(3,3)=1./sig

!- calculation of contravariant 4-velocity in KS coordinates -!        
     call calgfl(util,gfl,gcov1)
     call calucon(util,ucon,gfl,gcon1)
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
!
     do n=0,3
       do m=0,3
         gcon1a(m,n)=0.d0
       enddo
     enddo 

     gcon1a(0,0)=-aa/(sig*del)
     gcon1a(1,1)=del/sig                          
     gcon1a(2,2)=(sig-2.*rr)/(del*sig*sin(th)*sin(th)) 
     gcon1a(0,2)=-(2.*akm*rr)/(sig*del)           
     gcon1a(2,0)=gcon1(0,2)   
     gcon1a(3,3)=1./sig 
!
     alpha1=1./sqrt(-gcon1a(0,0))
     gfl=ucon1(0)*alpha1
!
     beta1(1)=alpha1*alpha1*gcon1a(0,1)
     beta1(2)=alpha1*alpha1*gcon1a(0,2)
     beta1(3)=alpha1*alpha1*gcon1a(0,3)
!
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
  deallocate(  ucon, ucon1, gcov1, gcon1, gcon1a, beta1, &
               tmp1, trans, stat=merr)
      
  return
end subroutine transks2bl1
