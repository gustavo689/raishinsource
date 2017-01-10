!***********************************************************************
!      RAISHIN code: 3D General Relativistic MHD (3D MPI) version 
!      program for output to VTK format
!      written by Y. Mizuno
!      ver. 150715
!      new version
!      using "uria" in main program (output1)
!***********************************************************************
!
program convert_idl
  use pram, only: imax, jmax, kmax, iprocs, jprocs, kprocs, metric

  implicit none

!======================================================================@
!    Difinition for variables
!======================================================================@
!
  integer, parameter :: nv=9
  integer, parameter :: npe=iprocs*jprocs*kprocs ! number of cpus
  integer, parameter :: ns=0, ne=1 ! start and end data file number
  integer, parameter :: ndim=1 ! output data dimension
  integer, parameter :: idirec=1 ! direction (1=x, 2=y, 3=z)

!  integer, parameter :: nq=13
!
  integer :: i, j, k, l, m, n
  integer :: is, ie, js, je, ks, ke
  integer :: myrank, nd, npe1, myrank1, merr
  integer :: ih, jh, kh, jh1, kh1
!  
  character*256 :: filename, filename1 
!
!  real(8) :: uri(nv,imax,jmax,kmax)
  real(8), allocatable :: uri1(:,:,:,:), detg1(:,:,:), x1a(:), x2a(:), x3a(:)
  real(8) :: time
!
  real(8) :: dd(imax,jmax,kmax), pp(imax,jmax,kmax), ink(imax,jmax,kmax)  
  real(8) :: ux(imax,jmax,kmax), uy(imax,jmax,kmax), uz(imax,jmax,kmax) 
  real(8) :: bx(imax,jmax,kmax), by(imax,jmax,kmax), bz(imax,jmax,kmax) 
  real(8) :: detg(imax,jmax,kmax)
  real(8) :: x1(imax), x2(jmax), x3(kmax)
 
  
!  real(8) :: qq(nq,imax,jmax,kmax)
!

  real(8), allocatable :: dir1dx(:), dir1dy(:), dir1dz(:), &
                          de1d(:), pp1d(:), ink1d(:), &  
                          ux1d(:), uy1d(:), uz1d(:), & 
                          bx1d(:), by1d(:), bz1d(:) 
  
  real(8), allocatable :: dir2dx(:,:), dir2dy(:,:), dir2dz(:,:), &
                          de2d(:,:), pp2d(:,:), ink2d(:,:), &  
                          ux2d(:,:), uy2d(:,:), uz2d(:,:), & 
                          bx2d(:,:), by2d(:,:), bz2d(:,:) 
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
990   format('structr',i4.4,'-',i4.4,'.outdat')
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
            dd(is+i,js+j,ks+k)=uri1(1,i,j,k)
            ux(is+i,js+j,ks+k)=uri1(2,i,j,k)
            uy(is+i,js+j,ks+k)=uri1(3,i,j,k)
            uz(is+i,js+j,ks+k)=uri1(4,i,j,k)
            pp(is+i,js+j,ks+k)=uri1(5,i,j,k)
            ink(is+i,js+j,ks+k)=uri1(6,i,j,k)
            bx(is+i,js+j,ks+k)=uri1(7,i,j,k)
            by(is+i,js+j,ks+k)=uri1(8,i,j,k)
            bz(is+i,js+j,ks+k)=uri1(9,i,j,k)
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
!- Output data for analysis in IDL -!
!
!- set output file -!
    write(filename1,991) nd
    open( unit=9, file=filename1,status='unknown')
991 format('ok',i3.3)  
!
!- convert the data array -!
!
    if (ndim .eq. 1) then !- 1D output -!
!
      if(idirec .eq. 1) then !- x direction -! 
!
!- allocate of 1D array -!
!
        allocate(dir1dx(imax), de1d(imax), ux1d(imax), uy1d(imax), uz1d(imax), &
                 pp1d(imax), ink1d(imax), bx1d(imax), by1d(imax), bz1d(imax), &
                 stat=merr)
!
!- set half grid position -!
!
        jh=(jmax+1)/2
        kh=(kmax+1)/2
!
!- convert -!
!        
        if(metric .eq. 1 .or. metric .eq. 2) then        
!- Cartesian, clyndrical => Cartesian             
          do i=1,imax
            dir1dx(i)=x1(i)
!
            de1d(i)=dd(i,jh,kh)
            ux1d(i)=ux(i,jh,kh)
            uy1d(i)=uy(i,jh,kh)
            uz1d(i)=uz(i,jh,kh)
            pp1d(i)=pp(i,jh,kh)
            ink1d(i)=ink(i,jh,kh)
            bx1d(i)=bx(i,jh,kh)
            by1d(i)=by(i,jh,kh)
            bz1d(i)=bz(i,jh,kh)
          enddo  
        elseif(metric .eq. 103 .or. metric .eq. 203 .or. metric .eq. 303) then
!- Spherical => Cartesian
          do i=1,imax
            dir1dx(i)=x1(i)
!
            de1d(i)=dd(i,jh,kh)
            ux1d(i)=ux(i,jh,kh)
            uy1d(i)=uy(i,jh,kh)
            uz1d(i)=-uz(i,jh,kh)
            pp1d(i)=pp(i,jh,kh)
            ink1d(i)=ink(i,jh,kh)
            bx1d(i)=bx(i,jh,kh)
            by1d(i)=by(i,jh,kh)
            bz1d(i)=-bz(i,jh,kh)

            write(*,*) 'i,pp=', i, pp(i,jh,kh)
          enddo
        endif
!
!- output -!
!
        do i=1,imax
          write(9,800) dir1dx(i),de1d(i),ux1d(i),uy1d(i),uz1d(i), &
                       pp1d(i),bx1d(i),by1d(i),bz1d(i)
        enddo
        write(9,*) time,'  time'
!
        deallocate(dir1dx, de1d, ux1d, uy1d, uz1d, &
                   pp1d, ink1d, bx1d, by1d, bz1d, stat=merr)
!!!
      elseif(idirec .eq. 2) then !- y-direction -!
!
!- allocate of 1D array -!
!
        allocate(dir1dy(jmax), de1d(jmax), ux1d(jmax), uy1d(jmax), uz1d(jmax), &
                 pp1d(jmax), ink1d(jmax), bx1d(jmax), by1d(jmax), bz1d(jmax), &
                 stat=merr)
!
!- set half grid position -!
!
        ih=(imax+1)/2
        kh=(kmax+1)/2
!
!- convert -!
!
        if(metric .eq. 1) then
!- Cartesian, => Cartesian
          do j=1,jmax
            dir1dy(j)=x2(j)
!
            de1d(j)=dd(ih,j,kh)
            ux1d(j)=ux(ih,j,kh)
            uy1d(j)=uy(ih,j,kh)
            uz1d(j)=uz(ih,j,kh)
            pp1d(j)=pp(ih,j,kh)
            ink1d(j)=ink(ih,j,kh)
            bx1d(j)=bx(ih,j,kh)
            by1d(j)=by(ih,j,kh)
            bz1d(j)=bz(ih,j,kh)                   
          enddo
        endif
!       
!- output -!
!
        do j=1,jmax
          write(9,800) dir1dy(j),de1d(j),ux1d(j),uy1d(j),uz1d(j), &
                       pp1d(j),bx1d(j),by1d(j),bz1d(j)
        enddo
        write(9,*) time,'  time'
!
        deallocate(dir1dy, de1d, ux1d, uy1d, uz1d, &
                   pp1d, ink1d, bx1d, by1d, bz1d, stat=merr)
!!!
      elseif(idirec .eq. 3) then !- z-direction -!
!
!- allocate of 1D array -!
!
        allocate(dir1dz(kmax), de1d(kmax), ux1d(kmax), uy1d(kmax), uz1d(kmax), &
                 pp1d(kmax), ink1d(kmax), bx1d(kmax), by1d(kmax), bz1d(kmax), &
                 stat=merr)
!        
!- set half grid position -!
!
        ih=(imax+1)/2
        jh=(jmax+1)/2
!
!- convert -!
!
        if(metric .eq. 1) then
!- Cartesian, => Cartesian 
          do k=1,kmax  
            dir1dz(k)=x3(k)
!
            de1d(k)=dd(ih,jh,k)
            ux1d(k)=ux(ih,jh,k)
            uy1d(k)=uy(ih,jh,k)
            uz1d(k)=uz(ih,jh,k)
            pp1d(k)=pp(ih,jh,k)
            ink1d(k)=ink(ih,jh,k)
            bx1d(k)=bx(ih,jh,k)
            by1d(k)=by(ih,jh,k)
            bz1d(k)=bz(ih,jh,k)  
          enddo
        endif
!       
!- output -!
!
        do k=1,kmax
          write(9,800) dir1dz(k),de1d(k),ux1d(k),uy1d(k),uz1d(k), &
                       pp1d(k),bx1d(k),by1d(k),bz1d(k)
        enddo
        write(9,*) time,'  time'
!
        deallocate(dir1dz, de1d, ux1d, uy1d, uz1d, &
                   pp1d, ink1d, bx1d, by1d, bz1d, stat=merr)
!
      endif
!
!=========
!
    elseif(ndim .eq. 2) then !- 2D output -!
!!!
      if(idirec .eq. 13) then !- xz plane -!
!
!- set half grid position -!
        jh1=(jmax+1)/2
!
!- convert -!
!
        if(metric .eq. 1 .or. metric .eq. 2) then         
!
!- data allocate -!
          allocate(dir2dx(imax,kmax), dir2dz(imax,kmax), de2d(imax,kmax), &
                   ux2d(imax,kmax), uy2d(imax,kmax), uz2d(imax,kmax), &
                   pp2d(imax,kmax), ink2d(imax,kmax), &
                   bx2d(imax,kmax), by2d(imax,kmax), bz2d(imax,kmax), &
                   stat=merr)
!
          do k=1,kmax
            do i=1,imax
              dir2dx(i,k)=x1(i)
              dir2dz(i,k)=x3(k)
!
              de2d(i,k)=dd(i,jh,k)
              ux2d(i,k)=ux(i,jh,k)
              uy2d(i,k)=uy(i,jh,k)
              uz2d(i,k)=uz(i,jh,k)
              pp2d(i,k)=pp(i,jh,k)
              ink2d(i,k)=ink(i,jh,k)
              bx2d(i,k)=bx(i,jh,k)
              by2d(i,k)=by(i,jh,k)
              bz2d(i,k)=bz(i,jh,k) 
            enddo
          enddo
!         
!- Output -!
!
          do k=1,kmax
            do i=1,imax

              write(9,800) dir2dx(i,k),dir2dz(i,k),de2d(i,k), &
                           ux2d(i,k),uy2d(i,k),uz2d(i,k), &
                           pp2d(i,k),bx2d(i,k),by2d(i,k),bz2d(i,k)
            enddo
          enddo
          write(9,*) time,'  time'
!
          deallocate(dir2dx, dir2dz, de2d, ux2d, uy2d, uz2d, &
                     pp2d, ink2d, bx2d, by2d, bz2d, stat=merr)
!      
        endif
!!!
      elseif(idirec .eq. 23) then !- yz-direction -!
!
!- set half grid position -!
!
         ih=(imax+1)/2
!
!- convert -!
!         
        if(metric .eq. 1) then                 
!- Cartesian => Cartesian
!- allocate of 2D array -!
          allocate(dir2dy(jmax,kmax), dir2dz(jmax,kmax), de2d(jmax,kmax), &
                   ux2d(jmax,kmax), uy2d(jmax,kmax), uz2d(jmax,kmax), &
                   pp2d(jmax,kmax), ink2d(jmax,kmax), &
                   bx2d(jmax,kmax), by2d(jmax,kmax), bz2d(jmax,kmax), &
                   stat=merr)

          do k=1,kmax
            do j=1,jmax
              dir2dy(j,k)=x2(j)
              dir2dz(j,k)=x3(k)
!
              de2d(j,k)=dd(ih,j,k)
              ux2d(j,k)=ux(ih,j,k)
              uy2d(j,k)=uy(ih,j,k)
              uz2d(j,k)=uz(ih,j,k)
              pp2d(j,k)=pp(ih,j,k)
              ink2d(j,k)=ink(ih,j,k)
              bx2d(j,k)=bx(ih,j,k)
              by2d(j,k)=by(ih,j,k)
              bz2d(j,k)=bz(ih,j,k) 
!                     
            enddo
          enddo
!         
!- Output -!
!
          do j=1,jmax
            do k=1,kmax
              write(9,800) dir2dy(j,k),dir2dz(j,k),de2d(j,k), &
                           ux2d(j,k),uy2d(j,k),uz2d(j,k), &
                           pp2d(j,k),bx2d(j,k),by2d(j,k),bz2d(j,k)
            enddo
          enddo
          write(9,*) time,'  time'
!
          deallocate(dir2dy, dir2dz, de2d, ux2d, uy2d, uz2d, &
                     pp2d, ink2d, bx2d, by2d, bz2d, stat=merr)
!
        endif
!!!
      elseif(idirec .eq. 12) then !- xy-direction -!
!
!- set half grid position -!
!
        kh=(kmax+1)/2
!
!- convert -!
!        
        if(metric .eq. 1) then
!- Cartesian => Cartesian
!- allocate of 2D array -!
          allocate(dir2dx(imax,jmax), dir2dy(imax,jmax), de2d(imax,jmax), &
                   ux2d(imax,jmax), uy2d(imax,jmax), uz2d(imax,jmax), &
                   pp2d(imax,jmax), ink2d(imax,jmax), &
                   bx2d(imax,jmax), by2d(imax,jmax), bz2d(imax,jmax), &
                   stat=merr)

          do j=1,jmax
            do i=1,imax
              dir2dx(i,j)=x1(i)
              dir2dy(i,j)=x2(j)
!
              de2d(i,j)=dd(i,j,kh)
              ux2d(i,j)=ux(i,j,kh)
              uy2d(i,j)=uy(i,j,kh)
              uz2d(i,j)=uz(i,j,kh)
              pp2d(i,j)=pp(i,j,kh)
              ink2d(i,j)=ink(i,j,kh)
              bx2d(i,j)=bx(i,j,kh)
              by2d(i,j)=by(i,j,kh)
              bz2d(i,j)=bz(i,j,kh)                      
            enddo
          enddo
!         
!- Output -!
!
          do i=1,imax
            do j=1,jmax
              write(9,800) dir2dx(i,j),dir2dy(i,j),de2d(i,j), &
                           ux2d(i,j),uy2d(i,j),uz2d(i,j), &
                           pp2d(i,j),bx2d(i,j),by2d(i,j),bz2d(i,j)
            enddo
          enddo
          write(9,*) time,'  time'
!
          deallocate(dir2dx, dir2dy, de2d, ux2d, uy2d, uz2d, &
                     pp2d, ink2d, bx2d, by2d, bz2d, stat=merr)
!
        endif
!
      endif     
       
    endif   
!   
    close(9)
!          
 enddo
!
 800 format(1h ,1pe12.4,21(1pe12.4))
! 
  stop
!
end program convert_idl

 

