!***********************************************************************
!      RAISHIN code: 3D General Relativistic MHD (3D MPI) version 
!      program for output to VTK format
!      written by Y. Mizuno
!      ver. 150715
!      new version
!      using "uria" in main program (output1)
!***********************************************************************
!
program convert_idl2
  use pram, only: imax, jmax, kmax, iprocs, jprocs, kprocs, akm, metric, &
                  gam, ieos, ix1, ix2, ix3, hslope, R0, pi

  implicit none

!======================================================================@
!    Difinition for variables
!======================================================================@
!
  integer, parameter :: nv=9
  integer, parameter :: npe=iprocs*jprocs*kprocs ! number of cpus 
  integer, parameter :: ns=0, ne=1 ! start and end data file number
  integer, parameter :: dataformat=0 !- 0= ascii, 1=binary -!
  integer, parameter :: idirec=1 ! direction (1=x, 2=y, 3=z)

!  integer, parameter :: nq=13
!
  integer :: i, j, k, l, m, n
  integer :: is, ie, js, je, ks, ke
  integer :: myrank, nd, npe1, myrank1, merr
  integer :: jh, kh, jh1, kh1, mmax, nmax
  integer :: kint, nn
!  
  character*256 :: filename, filename1 
!
  real(8) :: uri(nv,imax,jmax,kmax)
  real(8) :: detg(imax,jmax,kmax), gfla(imax,jmax,kmax), gflb(imax,jmax,kmax), &
             bsq(imax,jmax,kmax), hut(imax,jmax,kmax)
  real(8) :: x1(imax), x2(jmax), x3(kmax)
  real(8) :: rr, th, ph, sth, cth, del, sig, aa, gfl2a, gfl2b, &
             bbsq1a, de, pr, roh
  real(8) :: alpha1, hut1a, bsq1a, &
             gcon1a00, gcon1a01, gcon1a02, gcon1a03, gcov1a00
  real(8) :: util(1:3), gcov1(0:3,0:3), gcon1(0:3,0:3), ucon(0:3), ucov(0:3), &
             bcon(1:3), bcov(1:3), bbcon(0:3), bbcov(0:3), uconBL(0:3), &
             beta1(1:3), utilKS(1:3)
  real(8) :: gcov3a(1:3,1:3), gcon3a(1:3,1:3)
!
  real(8), allocatable :: uri1(:,:,:,:), detg1(:,:,:), x1a(:), x2a(:), x3a(:)
  real(8), allocatable :: gfl1a(:,:,:), gfl1b(:,:,:), bsq1(:,:,:), hut1(:,:,:)
  real(8) :: time, rbh1
  real(8), parameter :: bmin1=1.d-20 
  real(8) :: aph1
  real(8) :: x1t, x2t, x3t, tfac, rfac, pfac, hfac
  
!  real(8) :: qq(nq,imax,jmax,kmax)
! 
  real(8), allocatable :: dir1dx(:), dir1dy(:), dir1dz(:)
  real(8), allocatable :: de1d(:), pp1d(:), ink1d(:)  
  real(8), allocatable :: u11d(:), u21d(:), u31d(:), gfa1d(:), gfb1d(:) 
  real(8), allocatable :: b11d(:), b21d(:), b31d(:), bbsq1d(:), huut1d(:), &
                          te1d(:), pbeta1d(:)
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
               detg1(0:ie-is,0:je-js,0:ke-ks), &
               gfl1a(0:ie-is,0:je-js,0:ke-ks), gfl1b(0:ie-is,0:je-js,0:ke-ks), &
               bsq1(0:ie-is,0:je-js,0:ke-ks), &
               hut1(0:ie-is,0:je-js,0:ke-ks), &
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
            rbh1=1.+sqrt(1.-akm*akm)
!
             if(metric .eq. 403) then
              if(ix1 .eq. 1) then
                rr=x1a(i)
                x1t=x1a(i) 
              elseif(ix1 .eq. 2) then
                rr=R0+exp(x1a(i))
                x1t=x1a(i)                 
              endif  
!
              ph=x2a(j)
              x2t=x2a(j)
!               
              if(ix3 .eq. 1) then
                th=x3a(k)
                x3t=x3a(k) 
              elseif(ix3 .eq. 2) then
                th=x3a(k)+0.5*(1.-hslope)*sin(2.*x3a(k))
                x3t=x3a(k) 
              endif   
            else   
              rr=x1a(i)
              ph=x2a(j)
              th=x3a(k)
            endif  

!- set metric -!            
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
              
            if(metric .eq. 203) then

              gcov1(0,0)=-(sig-2.*rr)/sig     
              gcov1(1,1)=sig/del               
              gcov1(2,2)=(aa/sig)*sth**2 
              gcov1(0,2)=-(2.*akm*rr*(sth**2))/sig 
              gcov1(2,0)=gcov1(0,2)              
              gcov1(3,3)=sig                  
!
              gcon1(0,0)=-aa/(sig*del)
              gcon1(1,1)=del/sig                          
              gcon1(2,2)=(sig-2.*rr)/(del*sig*sin(th)*sin(th)) 
              gcon1(0,2)=-(2.*akm*rr)/(sig*del)           
              gcon1(2,0)=gcon1(0,2)   
              gcon1(3,3)=1./sig 
!               
            elseif(metric .eq. 303) then

              gcov1(0,0)=-1.+2.*rr/sig     
              gcov1(0,1)=2.*rr/sig          
              gcov1(0,2)=-2.*akm*rr*sth*sth/sig 
              gcov1(1,0)=gcov1(0,1)        
              gcov1(1,1)=1.+2.*rr/sig      
              gcov1(1,2)=-1.*akm*sth*sth*(1.+2.*rr/sig) 
              gcov1(2,0)=gcov1(0,2)        
              gcov1(2,1)=gcov1(1,2)        
              gcov1(2,2)=aa*sth*sth/sig 
              gcov1(3,3)=sig
              
              gcon1(0,0)=-1.-2.*rr/sig
              gcon1(0,1)=2.*rr/sig
              gcon1(1,0)=gcon1(0,1)
              gcon1(1,1)=del/sig
              gcon1(1,2)=akm/sig
              gcon1(2,1)=gcon1(1,2)
              gcon1(2,2)=1./(sig*sth*sth)
              gcon1(3,3)=1./sig

            elseif(metric .eq. 403) then
!           
!              call get_mks_position(x1t,x2t,x3t,i,j,k)
!              
!- get modification parameter for Modified Kerr-Schild coordinates -!
               
              tfac=1.d0

              if(ix1 .eq. 1) then
                rfac=1.d0
              elseif(ix1 .eq. 2) then
                rfac=rr
              endif  

              pfac=1.d0

              if(ix3 .eq. 1 ) then
                hfac=1.d0
              elseif(ix3 .eq. 2) then
                hfac=1.+(1.-hslope)*cos(2.*x3a(k)) 
              endif              
              
              gcov1(0,0)=(-1.+2.*rr/sig)*tfac*tfac     
              gcov1(0,1)=(2.*rr/sig)*tfac*rfac          
              gcov1(0,2)=(-2.*akm*rr*sth*sth*(1./sig))*tfac*pfac 
              gcov1(1,0)=gcov1(0,1)        
              gcov1(1,1)=(1.+2.*rr/sig)*rfac*rfac      
              gcov1(1,2)=(-1.*akm*sth*sth*(1.+2.*rr/sig))*rfac*pfac 
              gcov1(2,0)=gcov1(0,2)        
              gcov1(2,1)=gcov1(1,2)        
              gcov1(2,2)=(aa*sth*sth*(1./sig))*pfac*pfac 
              gcov1(3,3)=(sig)*hfac*hfac
              call invert_matrix2(gcov1,gcon1,4)              
               
            endif   
!
            call cal3met(gcov1,gcon1,gcov3a,gcon3a)
!
            de=uri1(1,i,j,k)            
            util(1)=uri1(2,i,j,k)
            util(2)=uri1(3,i,j,k)
            util(3)=uri1(4,i,j,k)
            pr=uri1(5,i,j,k)
!
!- calculation of relativistic enthalpy -!        
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          roh=de+(gam/(gam-1.0))*pr
        elseif(ieos .eq. 1) then
          roh=(5./2.)*pr+sqrt((9./4.)*pr**2+de**2)
        endif            
!
!- set minimum magnetic field            
!
            if(abs(uri1(7,i,j,k)) .lt. bmin1) then
              uri1(7,i,j,k)=0.d0               
            endif   
            if(abs(uri1(8,i,j,k)) .lt. bmin1) then
              uri1(8,i,j,k)=0.d0
            endif
            if(abs(uri1(9,i,j,k)) .lt. bmin1) then
              uri1(9,i,j,k)=0.d0
            endif   

            bcon(1)=uri1(7,i,j,k)
            bcon(2)=uri1(8,i,j,k)
            bcon(3)=uri1(9,i,j,k)
!            
!- calculation of Lorentz factor -!
            call calgfl(util,gfl2a,gcov1)
!- cal of contravariant 4-velocity -!
            call calucon(util,ucon,gfl2a,gcon1)
!- cal of covariant 4-velocity -! 
            call lower(ucon,ucov,gcov1)
!- cal of contravariant 4-magnetic field -!
            call calbbcon(bcon,bbcon,ucov,ucon)
!- cal of covariant 4-magnetic field -! 
            call lower(bbcon,bbcov,gcov1) 
!- cal of 4-magnetic field square -! 
            call calbbsq(bbcon,bbcov,bbsq1a)           

            if(bbsq1a .lt. 0.d0) then
              write(*,*) "bbsq < 0", bbsq1a, rr, ph, th
            endif
            
            bsq1(i,j,k)=bbsq1a
!-
!- cal of covariant 3-magnetic field -!
            call lower1(bcon,bcov,gcov3a)
!
!- cal of 3-magentic field square -!
            call calbsq(bcon,bcov,bsq1a)
            if(bsq1a .lt. 0.d0) then
              write(*,*) "bsq < 0", bsq1a, rr, ph, th
            endif
!            if(j .eq. 4) then
!              write(*,*) "bbsq, bsq=", bbsq1a, bsq1a, rr, th 
!            endif   
!
!
            if(metric .eq. 303 .or. metric .eq. 403) then
!- trans KS velocity to BL velocity & calculation of Lorentz factor-!
!               
              if(metric .eq. 403) then
                call transmks2ks1(util,utilKS,x1t,x2t,x3t)
                call transks2bl1(utilKS,uconBL,rr,ph,th) 
              elseif(metric .eq. 303) then
                call transks2bl1(util,uconBL,rr,ph,th)                
              endif  
!             
              if(rr .gt. rbh1) then
                sth=sin(th)
                cth=cos(th)
                del=rr*rr-(2.*rr)+akm*akm
                sig=rr*rr+akm*cth*akm*cth
                aa=((rr*rr+akm*akm)**2)-akm*akm*del*sth*sth
!            
                gcon1a00=-aa/(sig*del) !- g^tt in BL coordinate
                gcon1a01=0.d0
                gcon1a02=-(2.*akm*rr)/(sig*del) ! g^t\phi in BL coordinate
                gcon1a03=0.d0
                gcov1a00=-(sig-2.*rr)/sig !- g_tt in BL coordinate
              
                alpha1=1./sqrt(-gcon1a00)
                
                beta1(1)=alpha1*alpha1*gcon1a01
                beta1(2)=alpha1*alpha1*gcon1a02
                beta1(3)=alpha1*alpha1*gcon1a03
!
                gfl2a=alpha1*uconBL(0)  
                gfl2b=uconBL(0)*sqrt(-gcov1a00) !- Lorentz factor at static obs at infinity -!
                uri1(2,i,j,k)=uconBL(1)+gfl2a*beta1(1)*(1./alpha1)
                uri1(3,i,j,k)=uconBL(2)+gfl2a*beta1(2)*(1./alpha1)
                uri1(4,i,j,k)=uconBL(3)+gfl2a*beta1(3)*(1./alpha1)
!                
              else
                gfl2a=1.d0
                gfl2b=1.d0 
                uri1(2,i,j,k)=0.d0
                uri1(3,i,j,k)=0.d0
                uri1(4,i,j,k)=0.d0   
              endif

            else
!- Kerr-BL coordinates
              alpha1=1./sqrt(-gcon1(0,0))
              gfl2a=alpha1*ucon(0)            !- Lorentz factor in ZAMO obs -!
              gfl2b=ucon(0)*sqrt(-gcov1(0,0)) !- Lorentz factor at static obs at infinity -!
            endif
          
            hut1a=-(roh/de)*ucov(0)
            
            gfl1a(i,j,k)=gfl2a
            gfl1b(i,j,k)=gfl2b
            
            hut1(i,j,k)=hut1a
            
          enddo
        enddo            
     enddo
!     
!- convert data on full array -!       
      do i=0,ie-is
        if(metric .eq. 403) then
          if(ix1 .eq. 1) then
            x1(is+i)=x1a(i)
          elseif(ix1 .eq. 2) then
            x1(is+i)=R0+exp(x1a(i))
          endif
        else
          x1(is+i)=x1a(i)
        endif
      enddo
!
      do j=0,je-js
        x2(js+j)=x2a(j)
      enddo
!
      do k=0,ke-ks
        if(metric .eq. 403) then
          if(ix3 .eq. 1) then
            x3(ks+k)=x3a(k)  
          elseif(ix3 .eq. 2) then
            x3(ks+k)=x3a(k)+0.5*(1.-hslope)*sin(2.*x3a(k))  
          endif
        else   
          x3(ks+k)=x3a(k)
        endif   
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
            gfla(is+i,js+j,ks+k)=gfl1a(i,j,k)
            gflb(is+i,js+j,ks+k)=gfl1b(i,j,k)
            bsq(is+i,js+j,ks+k)=bsq1(i,j,k)
            hut(is+i,js+j,ks+k)=hut1(i,j,k)
          enddo
        enddo
      enddo    
!
      close(7)
      deallocate(uri1,detg1,gfl1a,gfl1b,bsq1,hut1,x1a,x2a,x3a,stat=merr)
    enddo
!
!----------------------------------------------------------------------@
!- Output data for analysis -!
!
!- set output file -!
    write(filename1,991) nd
    open( unit=9, file=filename1,status='unknown')
991 format('ok',i3.3)  
!
    if(idirec .eq. 1) then !- x plane -!
!
!- set half grid position -!
      jh1=(jmax+1)/2
      kh1=(kmax+1)/2       
!- set j-th direction =1 -!
      mmax=1
!
!- data allocate -!
      allocate(dir1dx(imax), &
               de1d(imax), pp1d(imax), ink1d(imax), &
               u11d(imax), u21d(imax), u31d(imax), &
               b11d(imax), b21d(imax), b31d(imax), &
               gfa1d(imax), gfb1d(imax), bbsq1d(imax), &
               huut1d(imax), te1d(imax), pbeta1d(imax), stat=merr)
!
!- convert the data array -!
!
        do i=1, imax
!    
          dir1dx(i)=x1(i)
            
          de1d(i)=uri(1,i,jh1,kh1) !- density -!            
          u11d(i)=uri(2,i,jh1,kh1)
          u21d(i)=uri(4,i,jh1,kh1)
          u31d(i)=uri(3,i,jh1,kh1)
          pp1d(i)=uri(5,i,jh1,kh1) !- gas pressure -!
          b11d(i)=uri(7,i,jh1,kh1)
          b21d(i)=uri(9,i,jh1,kh1)
          b31d(i)=uri(8,i,jh1,kh1)
                
          gfa1d(i)=gfla(i,jh1,kh1)
          gfb1d(i)=gflb(i,jh1,kh1)
          bbsq1d(i)=bsq(i,jh1,kh1)
          huut1d(i)=hut(i,jh1,kh1)

          te1d(i)=pp1d(i)/de1d(i)
          if(bbsq1d(i) .ne. 0.d0) then
            pbeta1d(i)=2.*pp1d(i)/bbsq1d(i)
          else
            pbeta1d(i)=1.d10
          endif   
        enddo
!
!- output -!
!
        write(9,*)'# r, \rho, p_g, temperature, sqrt(b^2), p_mag, beta' 
        do i=1,imax
!          write(9,800) dir1dx(i),de1d(i),u11d(i),u21d(i),u31d(i), &
!                       pp1d(i),b11d(i),b21d(i),b31d(i),bbsq1d(i)
           write(9,800) dir1dx(i),de1d(i),pp1d(i),te1d(i), &
                        sqrt(bbsq1d(i)), 0.5d0*bbsq1d(i), pbeta1d(i)
        enddo
!        write(9,*) time,'  time'     
!
      deallocate(dir1dx, de1d, pp1d, ink1d, u11d, u21d, u31d, &
                 b11d, b21d, b31d, gfa1d, gfb1d, bbsq1d, huut1d, stat=merr)
!
    elseif(idirec .eq. 2) then !- xy plane -!

    endif
!      
  enddo    
!
 800 format(1h ,1pe12.4,21(1pe12.4))
!  
  stop    
end program convert_idl2

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
subroutine calbsq(bcon,bcov,bsq)
!--------------------------------------------------------------------
!- Calculation of square of 3-magnetic field -!
!
  implicit none

!- bcon: contravariant 3-magnetic field -!
!- bcov: covariant 3-magnetic field -!

  real(8) :: bcon(1:3), bcov(1:3)
  real(8) :: bsq

  bsq=bcon(1)*bcov(1)+bcon(2)*bcov(2)+bcon(3)*bcov(3)

  return
end subroutine calbsq
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
subroutine lower1(vcon,vcov,gcov3a)
!--------------------------------------------------------------------
!- contravariant 4-vector => covariant 4-vector -!
!
  implicit none
!
!- vcon: contravariant 3-vector (input) -!
!- vcov: covariant 3-vector (output) -!
!- gcov3a: covariant 3-metric (input) -! 
!
  real(8) :: vcon(1:3), vcov(1:3) 
  real(8) :: gcov3a(1:3,1:3)


  vcov(1)=gcov3a(1,1)*vcon(1) &
         +gcov3a(1,2)*vcon(2) &
         +gcov3a(1,3)*vcon(3)

  vcov(2)=gcov3a(2,1)*vcon(1) &
         +gcov3a(2,2)*vcon(2) &
         +gcov3a(2,3)*vcon(3)

  vcov(3)=gcov3a(3,1)*vcon(1) &
         +gcov3a(3,2)*vcon(2) &
         +gcov3a(3,3)*vcon(3)

  return
end subroutine lower1
!
!--------------------------------------------------------------------
subroutine transks2bl1(util,uconBL,x1aa,x2aa,x3aa)
!--------------------------------------------------------------------
!- convert 4-velocity in KS to 4-velocity in BL coordinates -!
!-  reading whole primitive variables -!
!
!- util: \tilde{u}_KS (input)-!
!- uconBL: \u_BL (output)-! 
!
!  
  use pram, only : metric, akm, ix1, ix2, ix3, R0, hslope, pi  
  implicit none

  integer :: m, n, merr
       
  real(8) :: util(1:3), uconBL(0:3)

  real(8), allocatable :: ucon(:)
  real(8), allocatable :: gcov1(:,:), gcon1(:,:), gcon1a(:,:)
  real(8), allocatable :: tmp1(:), trans(:,:)

  real(8) :: x1aa, x2aa, x3aa, rr, ph, th, rbh1, alpha1, gfl, tmp2, &
             sth, cth, del, sig, aa
!
!- allocate variables -!
  allocate(ucon(0:3), gcov1(0:3,0:3), gcon1(0:3,0:3), &
           gcon1a(0:3,0:3), tmp1(0:3), trans(0:3,0:3), stat=merr)
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
     gcov1(3,3)=sig

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
       uconBL(m)=0.d0
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
       uconBL(m)=tmp1(m)
     enddo
!
   else
     uconBL(0)=ucon(0)     
     uconBL(1)=0.d0
     uconBL(2)=0.d0
     uconBL(3)=0.d0
    endif
!
  deallocate(  ucon, gcov1, gcon1, gcon1a, tmp1, trans, stat=merr)
      
  return
end subroutine transks2bl1
!
!--------------------------------------------------------------------
subroutine transmks2ks1(util,utilKS,x1t,x2t,x3t)
!--------------------------------------------------------------------
!- Transformation from 3-velocity in MKS to 3-velocity in KS coordinates -!
!
  use pram, only: ix1, ix2, ix3, hslope, pi
  implicit none
!
  real(8) :: util(1:3), utilKS(1:3)
  real(8) :: x1t, x2t, x3t, rfac, pfac, hfac
!!!!!
  
  if(ix1 .eq. 1) then
    rfac=1.d0
  elseif(ix1 .eq. 2) then
    rfac=exp(x1t) 
  endif
!                
  pfac=1.d0
!
  if(ix3 .eq. 1) then
    hfac=1.d0
  elseif(ix3 .eq. 2) then
    hfac=1.+(1.-hslope)*cos(2.*x3t)                                       
  endif

  utilKS(1)=rfac*util(1)
  utilKS(2)=pfac*util(2)
  utilKS(3)=hfac*util(3)
  
  return
end subroutine transmks2ks1  
!
!--------------------------------------------------------------------
subroutine cal3met(gcov1,gcon1,gcov3a,gcon3a)
!--------------------------------------------------------------------
!- Calculation of spatial 3-metric -!
!
  implicit none

  integer :: m, n, merr

  real(8) ::  gcov1(0:3,0:3) !- convariant 4-metric -!
  real(8) ::  gcon1(0:3,0:3) !- contravariant 4-metric -!
  real(8) ::  gcov3a(1:3,1:3) !- convariant 3-metric -!
  real(8) ::  gcon3a(1:3,1:3) !- contravariant 3-metric -!

  real(8), allocatable :: oncov(:), oncon(:)
  real(8), allocatable :: beta1(:)
  real(8) :: alpha1 
!
  allocate(oncov(0:3), oncon(0:3), beta1(1:3), stat=merr)
!
  alpha1=1./sqrt(-gcon1(0,0))
  beta1(1)=alpha1*alpha1*gcon1(0,1)
  beta1(2)=alpha1*alpha1*gcon1(0,2)
  beta1(3)=alpha1*alpha1*gcon1(0,3)   
!
  oncon(0)=1./alpha1
  oncon(1)=-beta1(1)*(1./alpha1)
  oncon(2)=-beta1(2)*(1./alpha1)
  oncon(3)=-beta1(3)*(1./alpha1)
!
  oncov(0)=-alpha1
  oncov(1)=0.d0
  oncov(2)=0.d0
  oncov(3)=0.d0
!
  do m=1,3
    do n=1,3
      gcov3a(m,n)=gcov1(m,n)+(oncov(m)*oncov(n))
      gcon3a(m,n)=gcon1(m,n)+(oncon(m)*oncon(n))
    enddo
  enddo
!
  deallocate(oncov, oncon, beta1, stat=merr)
  return
end subroutine cal3met
!
!--------------------------------------------------------------------
subroutine get_mks_position(x1t,x2t,x3t,i,j,k)
!--------------------------------------------------------------------
  use pram, only : imax, jmax, kmax, xmin, xmax, ymin, ymax, zmin, zmax, &
                   ix1, ix2, ix3, R0, hslope, pi
  
  implicit none
!  
  integer :: i, j, k
!  
  real(8) :: x1t, x2t, x3t, dx1, dx2, dx3, xx0, x1stat, x1end, p1, &
             x3stat, x3end, hth, p2
  
  if(ix1 .eq. 1) then
     
    dx1=(xmax-xmin)/float(imax-1)
    x1t=xmin+dx1*float(i-1)

  elseif(ix1 .eq. 2) then

    xx0=R0
    x1stat=log(xmin-xx0)
    x1end=log(xmax-xx0)
    dx1=(x1end-x1stat)*(1./float(imax-1))
    x1t=x1stat+float(i-1)*dx1

  elseif(ix1 .eq. 3) then
    
  elseif(ix1 .eq. 4) then

    xx0=R0 
    x1stat=log10(xmin-xx0)
    x1end=log10(xmax-xx0)
!    dx1=log10((xmax-xx0)/(xmin-xx0))/float(imax-1)
    dx1=(x1end-x1stat)*(1./float(imax-1))
    x1t=x1stat+float(i-1)*dx1
    
  elseif(ix1 .eq. 5) then
    
    p1=2.d0  
    xx0=R0
    x1stat=log(xmin-xx0)**(1./p1)
    x1end=log(xmax-xx0)**(1./p1)
    dx1=(x1end-x1stat)*(1./float(imax-1))
    x1t=x1stat+float(i-1)*dx1   

  endif
!
  if(ix2 .eq. 1) then

    dx2=(ymax-ymin)/float(jmax-1)
    x2t=ymin+dx2*float(j-1)  

  elseif(ix2 .eq. 2) then

  endif   
!
  if(ix3 .eq. 1) then

    dx3=(zmax-zmin)/float(kmax-1)
    x3t=zmin+dx3*float(k-1)

  elseif(ix3 .eq. 2) then

    x3stat=zmin/pi
    x3end=zmax/pi
    dx3=(x3end-x3stat)/float(kmax-1)
    x3t=pi*(x3stat+float(k-1)*dx3)

  elseif(ix3 .eq. 3) then

    hth=0.15
    p2=5.0
    
    x3stat=zmin/pi
    x3end=zmax/pi 
    dx3=(x3end-x3stat)/float(kmax-1)
    x3t=pi*(x3stat+float(k-1)*dx3)
    
  endif
!
  return
end subroutine get_mks_position
  
!--------------------------------------------------------------------
subroutine invert_matrix2(a,axn,n)
!--------------------------------------------------------------------
!
  implicit none
!
  integer :: n, indx(1:n)
  integer :: i, j
  real(8) :: an(1:n,1:n), a(0:n-1,0:n-1), ax(1:n,1:n), axn(0:n-1,0:n-1)
  real(8) :: y(1:n,1:n), b(1:n)           
  real(8) :: d
!  
!- input matrix A copy to matrix An

  do i=1,n
    do j=1,n 
      an(i,j)=a(i-1,j-1)
    enddo
  enddo  
!
!- set up identity matrix
  do i=1,n
    do j=1,n
      if(i .eq. j) then  
        y(i,j)=1.d0
      else
        y(i,j)=0.d0
      endif   
    enddo
  enddo  
!
  call ludcmp(an,n,indx,d)
!
  do j=1,n
    do i=1,n 
      b(i)=y(i,j)
    enddo 
!   
    call lubksb(an,n,indx,b)
!
    do i=1,n
      ax(i,j)=b(i)
    enddo         
  enddo
!
  do j=1,n
    do i=1,n
      axn(i-1,j-1)=ax(i,j)
    enddo
  enddo
!        
  return
end subroutine invert_matrix2
!
!--------------------------------------------------------------------
subroutine ludcmp(a,n,indx,d)
!--------------------------------------------------------------------
!- input :: a(n,n), n
!- output :: a(n,n), indx, d  
  
  implicit none
!
  integer :: n, indx(1:n)
  integer, parameter :: nmax=500
  real(8) :: d, a(1:n,1:n)
  real(8), parameter :: tiny=1.e-20
!  
  integer :: i, imax, j, k
  real(8) :: aamax, dum, sum, vv(nmax)
!
!- parameter
  d=1.d0
!
  do i=1,n
!
    aamax=0.d0
!
    do j=1,n 
      if (abs(a(i,j)) .gt. aamax) then
        aamax=abs(a(i,j))
      endif
    enddo
!
    if (aamax .eq. 0.d0) then
      write(*,*) "sigular matrix in ludcmp"
    endif   
    vv(i)=1./aamax
!
  enddo   
!  
  do j=1,n
!
     do i=0,j-1
      sum=a(i,j)
      do k=1,i-1
        sum=sum-a(i,k)*a(k,j) 
      enddo
      a(i,j)=sum
    enddo
!
    aamax=0.d0
!
    do i=j,n
      sum=a(i,j)
      do k=1,j-1
         sum=sum-a(i,k)*a(k,j)
      enddo
      a(i,j)=sum
      dum=vv(i)*abs(sum)
      if(dum.ge.aamax) then
        imax=i
        aamax=dum
      endif  
    enddo
!
    if(j .ne. imax) then
      do k=1,n
        dum=a(imax,k)
        a(imax,k)=a(j,k)
        a(j,k)=dum
      enddo
      d=-d
      vv(imax)=vv(j)
    endif
!
    indx(j)=imax
!
    if(a(j,j) .eq. 0.d0) then
      a(j,j)=tiny  
    endif
!
    if(j .ne. n) then
      dum=1./a(j,j)
      do i=j+1,n
        a(i,j)=a(i,j)*dum
      enddo   
    endif   
!    
  enddo
!
  return
end subroutine ludcmp  
!
!--------------------------------------------------------------------
subroutine lubksb(a,n,indx,b)
!--------------------------------------------------------------------
!
  implicit none
!
  integer :: n, indx(n)
  real(8) :: a(1:n,1:n), b(1:n)
  integer :: i,ii,j,ll
  real(8) :: sum
!
  ii=0
!  
  do i=1,n
    ll=indx(i)
    sum=b(ll)
    b(ll)=b(i)
!
    if(ii .ne. 0) then
       
      do j=ii, i-1
        sum=sum-a(i,j)*b(j)
      enddo

    else if (sum .ne. 0) then
        ii=i  
    endif
!
    b(i)=sum
!
  enddo
!
  do i=n,1,-1
    sum=b(i)
    do j=i+1, n
      sum=sum-a(i,j)*b(j) 
    enddo   
!
    b(i)=sum/a(i,i)
!    
  enddo
!
  return
end subroutine lubksb
