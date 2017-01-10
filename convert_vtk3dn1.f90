!***********************************************************************
!      RAISHIN code: 3D general Relativistic MHD (3D MPI) version 
!      program for output to VTK format (output in Cartesian cooddinates)
!      written by Y. Mizuno
!      ver. 150517
!      for older output case (output "uri"or "uria" in main program,
!      output1)
!***********************************************************************
!
program convert_vtk3dn1
  use pram, only: imax, jmax, kmax, iprocs, jprocs, kprocs, metric

  implicit none

!======================================================================@
!    Difinition for variables
!======================================================================@
!
  integer, parameter :: nv=9
  integer, parameter :: npe=iprocs*jprocs*kprocs ! number of cpus
  integer, parameter :: ns=0, ne=1 ! start and end data file number
  integer, parameter :: dataformat=0 !- 0= ascii, 1=binary -!
  integer, parameter :: idirec=123 ! direction (1=x, 2=y, 3=z)

!  integer, parameter :: nq=13
!
  integer :: i, j, k, l, m, n
  integer :: is, ie, js, je, ks, ke
  integer :: myrank, nd, npe1, myrank1, merr
  integer :: jh1, mmax, nmax
!  
  character*256 :: filename, filename1 
!
  real(8) :: uri(nv,imax,jmax,kmax)
  real(8) :: detg(imax,jmax,kmax)
  real(8) :: x1(imax), x2(jmax), x3(kmax)
  real(8), allocatable :: uri1(:,:,:,:), detg1(:,:,:), x1a(:), x2a(:), x3a(:)
  real(8) :: time
  
!  real(8) :: qq(nq,imax,jmax,kmax)
! 
  real(4), allocatable :: xx(:,:,:), yy(:,:,:), zz(:,:,:)
  real(4), allocatable :: dd(:,:,:), pp(:,:,:), ink(:,:,:)  
  real(4), allocatable :: ux(:,:,:), uy(:,:,:), uz(:,:,:) 
  real(4), allocatable :: bx(:,:,:), by(:,:,:), bz(:,:,:) 
!

!
!
!======================================================================@
!     File Open
!======================================================================@
!
  open(unit=6,file='convert.outdat',status='unknown')
!  open(unit=8,file='structr.outdat',status='unknown',form='unformatted')
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
               detg1(0:ie-is,0:je-js,0:ke-ks),&
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
            detg(is+i,js+j,ks+k)=detg1(i,j,k)  
          enddo
        enddo
      enddo    
!
      close(7)
      deallocate(uri1,detg1,x1a,x2a,x3a,stat=merr)
    enddo
!
!----------------------------------------------------------------------@
!- Output data for analysis -!
!
!- set output file -!
    write(filename1,991) nd
    open( unit=9, file=filename1,status='replace',form="formatted")
991 format('ok',i3.3,'.vtk')  
!
    if(idirec .eq. 123) then
!
!- set j-th direction =jmax - overlap region
!      mmax=jmax-6
      mmax=jmax-5
!- data allocate -!
      allocate(xx(imax,mmax,kmax), yy(imax,mmax,kmax), zz(imax,mmax,kmax), &
               dd(imax,mmax,kmax), pp(imax,mmax,kmax), ink(imax,mmax,kmax), &
               ux(imax,mmax,kmax), uy(imax,mmax,kmax), uz(imax,mmax,kmax), &
               bx(imax,mmax,kmax), by(imax,mmax,kmax), bz(imax,mmax,kmax), &
               stat=merr)
!
!- convert the data array -!
!
      do k=1, kmax 
        do j=1, mmax
          do i=1, imax
!            jh1=j+3
            jh1=j+2
    
            xx(i,j,k)=real(x1(i)*sin(x3(k))*cos(x2(jh1)))
            yy(i,j,k)=real(x1(i)*sin(x3(k))*sin(x2(jh1)))
            zz(i,j,k)=real(x1(i)*cos(x3(k)))
         
            dd(i,j,k)=real(uri(1,i,jh1,k)) !- density -!
          ux(i,j,k)=real(uri(2,i,jh1,k)*sin(x3(k))*cos(x2(jh1)) &
                        -uri(3,i,jh1,k)*cos(x3(k)) &
                        +uri(4,i,jh1,k)*cos(x3(k))*cos(x2(jh1)) &
                        ) !- 1st component of velocity -!
          uy(i,j,k)=real(uri(2,i,jh1,k)*sin(x3(k))*sin(x2(jh1)) &
                        +uri(3,i,jh1,k)*cos(x3(k)) &
                        +uri(4,i,jh1,k)*cos(x3(k))*sin(x2(jh1)) &
                        ) !- 2nd component of velocity -!
          uz(i,j,k)=real(uri(2,i,jh1,k)*cos(x3(k)) &
                        -uri(4,i,jh1,k)*sin(x3(k)) &
                        ) !- 3rd component of velocity -!
            pp(i,j,k)=real(uri(5,i,jh1,k)) !- gas pressure -!
            ink(i,j,k)=real(uri(6,i,jh1,k)) !- tracer scalar -!
          bx(i,j,k)=real(uri(7,i,jh1,k)*sin(x3(k))*cos(x2(jh1)) &
                        -uri(8,i,jh1,k)*cos(x3(k)) &
                        +uri(9,i,jh1,k)*cos(x3(k))*cos(x2(jh1)) &
                        ) !- 1st component of B-field -!
          by(i,j,k)=real(uri(7,i,jh1,k)*sin(x3(k))*sin(x2(jh1)) &
                        +uri(8,i,jh1,k)*cos(x3(k)) &
                        +uri(9,i,jh1,k)*cos(x3(k))*sin(x2(jh1)) &
                        ) !- 2nd component of B-field -!
          bz(i,j,k)=real(uri(7,i,jh1,k)*cos(x3(k)) &
                        -uri(9,i,jh1,k)*sin(x3(k)) &
                        ) !- 3rd component of B-field -!
          enddo
        enddo
      enddo
!     
      call write_vtk3d(imax,mmax,kmax,xx,yy,zz,dd,ux,uy,uz,pp,&
                       bx,by,bz,nd,dataformat,filename1)
!
      deallocate(xx, yy, zz, dd, pp, ink, ux, uy, uz, bx, by, bz, stat=merr)
!            
    endif
!      
  enddo    

  stop    
end program convert_vtk3dn1

!------------------------------------------------------------
subroutine write_vtk3d(xmax,ymax,zmax,xx,yy,zz,dd,ux,uy,uz,pp,&
                       bx,by,bz,nd,dataformat,filename)
!------------------------------------------------------------
!   
!- write data for VTK format -!
!  Now the data must be cartesian coordinates!
!  The gird points does not need to be regular.!
!
  implicit none
!
  integer :: i, j, k, xmax, ymax, zmax, nd
  integer :: merr
  integer :: dataformat
  character*256 :: filename
!  
  real(4) :: xx(xmax,ymax,zmax), yy(xmax,ymax,zmax), zz(xmax,ymax,zmax)
  real(4) :: dd(xmax,ymax,zmax), pp(xmax,ymax,zmax)  
  real(4) :: ux(xmax,ymax,zmax), uy(xmax,ymax,zmax), uz(xmax,ymax,zmax) 
  real(4) :: bx(xmax,ymax,zmax), by(xmax,ymax,zmax), bz(xmax,ymax,zmax)
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
    write(9, "(a)") "GRMHD Simulation Result"    
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
    write(6,*) 'after header: nd=', nd
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
    do k = 1, zmax
      do j = 1, ymax
        do i = 1, xmax
          if(dataformat .eq. 0) then                 
            write(9,*) xx(i,j,k), yy(i,j,k), zz(i,j,k)
          else
            write(9) xx(i,j,k), yy(i,j,k), zz(i,j,k)
          endif
        enddo  
      enddo
    enddo
!
    write(6,*) 'after coordinates: nd=', nd
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
    write(6,*) 'after new line: nd=', nd   
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
      do k = 1, zmax 
        do j = 1, ymax
          do i = 1, xmax             
            dum_r = max(dd(i,j,k), 1e-30) ! for safety
            write(9,*) dum_r
          enddo
        enddo 
      enddo
    else
      do k = 1, zmax
        do j = 1, ymax
          do i = 1, xmax
            dum_r = max(dd(i,j,k), 1e-30) ! for safety
            write(9) dum_r
          enddo  
        enddo
      enddo
    endif
!
    write(6,*) 'after density: nd=', nd  
!
! Write velocity
!
    if(dataformat .ne. 0) then
      close(9) 
! reopen it as formatted and append
      open(unit=9, file=filename, status="old", form="formatted", &
           position="append")   
    endif
    
    write(9, "(a, a)") newline, "VECTORS velocity float"
    
    if(dataformat .ne. 0) then
      close(9) 
! reopen it as unformatted and append
      open(unit=9, file=filename, status="old", form="unformatted", &
           position="append")   
!           position="append", access="stream")   
    endif
! -- now dumps data
    if(dataformat .eq. 0) then
      do k = 1, zmax 
        do j = 1, ymax
          do i = 1, xmax
            write(9,*) ux(i,j,k), uy(i,j,k), uz(i,j,k)
          enddo           
        enddo
      enddo 
    else
      do k = 1, zmax
        do j = 1, ymax
          do i = 1, xmax
            write(9) ux(i,j,k), uy(i,j,k), uz(i,j,k)
          enddo
        enddo
      enddo 
    endif
!
    write(6,*) 'after velocity vector: nd=', nd  
!
! Write B-field vector
!
    if(dataformat .ne. 0) then
      close(9) 
! reopen it as formatted and append
      open(unit=9, file=filename, status="old", form="formatted", &
           position="append")   
    endif
    
    write(9, "(a, a)") newline, "VECTORS B-field float"
    
    if(dataformat .ne. 0) then
      close(9) 
! reopen it as unformatted and append
      open(unit=9, file=filename, status="old", form="unformatted", &
           position="append")   
 !          position="append", access="stream") 
    endif
! -- now dumps data
    if(dataformat .eq. 0) then
      do k = 1, zmax 
        do j = 1, ymax
          do i = 1, xmax
            write(9,*) bx(i,j,k), by(i,j,k), bz(i,j,k)
          enddo
        enddo 
      enddo
    else
      do k = 1, zmax
        do j = 1, ymax
          do i = 1, xmax
            write(9) bx(i,j,k), by(i,j,k), bz(i,j,k)
          enddo
        enddo 
      enddo           
    endif
!
    write(6,*) 'end of the data: nd=', nd  
!
    close(9)
!
  return
end subroutine write_vtk3d   
