!-----------------------------------------------------------------------
subroutine recov(uu,uri,gcov,gcon,detg,x1,x2,x3,nm1, &
                 is1,ie1,js1,je1,ks1,ke1)
!-----------------------------------------------------------------------
!     Calculation of primitive variables from conserved variables
!     
  use pram, only : imax, jmax, kmax, nv, iwvec, iter, ieos
  implicit none
!
  integer :: nm1, is1, ie1, js1, je1, ks1, ke1
  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)

  real(8) :: x1(imax), x2(jmax), x3(kmax)
!  
!- 2 variables (W, v^2) following Noble et al. (2006)  
  if( iwvec.eq.1 ) then
!     call recov1(uu,uri,gcov,gcon,detg,x1,x2,x3,nm1,is1,ie1,js1,je1,ks1,ke1)
     call recov1a(uu,uri,gcov,gcon,detg,x1,x2,x3,nm1,is1,ie1,js1,je1,ks1,ke1)
!     call recov1b(uu,uri,gcov,gcon,detg,x1,x2,x3,nm1,is1,ie1,js1,je1,ks1,ke1)
!- 1 variable (W) following Mignone & McKinney (2007)
  elseif( iwvec.eq.2 ) then
     call recov2(uu,uri,gcov,gcon,detg,x1,x2,x3,nm1,is1,ie1,js1,je1,ks1,ke1)
!- 1 variables (W) following Noble et al. (2006)  
  elseif( iwvec.eq.3 ) then
     call recov3(uu,uri,gcov,gcon,detg,x1,x2,x3,nm1,is1,ie1,js1,je1,ks1,ke1)
!- 1 variables (v^2) for polytropic EoS  
  elseif( iwvec.eq.4 .and. ieos .eq. 3) then
     call recov4(uu,uri,gcov,gcon,detg,x1,x2,x3,nm1,is1,ie1,js1,je1,ks1,ke1)
!- 1 variables (v^2) using entropy eqs.  
  elseif( iwvec.eq.5 ) then
     call recov5(uu,uri,gcov,gcon,detg,x1,x2,x3,nm1,is1,ie1,js1,je1,ks1,ke1)   
!- 1 variables (W) using entropy eqs. 
  elseif( iwvec.eq.6 ) then
     call recov6(uu,uri,gcov,gcon,detg,x1,x2,x3,nm1,is1,ie1,js1,je1,ks1,ke1)
  else
    write(*,*) "choose wrong inversion procedure, iwvec, ieos=", iwvec, ieos
    stop 
  endif
!
  return
end subroutine recov
!
!---------------------------------------------------------------------@
subroutine recov1(uu,uri,gcov,gcon,detg,x1,x2,x3,nm1,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!- Following Noble et al. (2006), 2-variables inversion procedure -!
!- (only ideal EoS)
  
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos, iter, dmin, pmin, &
                   kpol, akm
  implicit none
!
  integer :: i, j, k, m, n, nnn, nm1, is1, ie1, js1, je1, ks1, ke1
  integer :: merr, igue, irecov

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)

  real(8) :: x1(imax), x2(jmax), x3(kmax)

  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: qcon(:), qcov(:), bcon(:), bcov(:), &
                          oncov(:), oncon(:) 
  real(8), allocatable :: util(:), util1(:), qtcon(:)
      
  real(8) :: safty, saftybig, gam1, gam2, dv, alpha1, detg1 
  real(8) :: dro, ro, pp, ro1, pp1, roh, roe, bsq, qdotb, qdotbsq, qdotn, &
             qsq, qtsq, utsq, gfsq, gf, vsq, wg, wsq, xsq, qsq1, &
             delta, deltaend, df, f, xxn, yyn, xo, yo, w2, w3, gtmp, ptmp, &
             f1, g1, dfx, dfy, dgx, dgy, det, dx, dy, dpdw, dpdvsq 
  real(8) :: tmp1a, tmp1b, tmp1c, tmp1d, wgp, wgm, vsq1, tmp3a, tmp3b, &
             tmp3c, tmp4a, romin, prmin
  real(8) :: rin, rbh1
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), stat=merr)
  allocate(qcon(0:3), qcov(0:3), bcon(0:3), bcov(0:3), &
           oncov(0:3), oncon(0:3), util(1:3), util1(1:3), qtcon(1:3), stat=merr)
!
!--- define variable of safe number ----
      safty = 1.0d-4
      saftybig = 1.0d5
      dv=1.0d-2
      
      gam1=gam/(gam-1.0)
      gam2=1.0/gam1
      wg=0.d0

      igue=1
!
      rin=15.d0
      rbh1=1.+sqrt(1.-akm*akm)       
!
!=====================================================================@
!
! ----  Calculation of variables ---
!
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
       do i=is1+nm1,ie1-nm1

!- set parameter for Newton Raphson -!

        delta=1.d0
        deltaend = 1.d-8
        df=1.d0
        f=1.d0
!        irecov=0
!        
!- check the position within the BH radius -!
!
!        if(x1(i) .le. rbh1) then
!          irecov=2
!        else
          irecov=0
!        endif  
          
!- copy of metric terms -!
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo
!
        alpha1=1./sqrt(-gcon1(0,0))
!
        detg1=detg(i,j,k)
      
        dro=alpha1*uu(1,i,j,k)*(1./detg1)

!        qcov(0)=alpha1*(uu(5,i,j,k)+uu(1,i,j,k))*(1./detg1)
        qcov(0)=alpha1*(uu(5,i,j,k)-uu(1,i,j,k))*(1./detg1)
!        qcov(0)=alpha1*uu(5,i,j,k)*(1./detg1)
        qcov(1)=alpha1*uu(2,i,j,k)*(1./detg1)
        qcov(2)=alpha1*uu(3,i,j,k)*(1./detg1)
        qcov(3)=alpha1*uu(4,i,j,k)*(1./detg1)

        tmp4a=alpha1*uu(5,i,j,k)*(1./detg1)
        
        bcon(0)=0.d0
        bcon(1)=alpha1*uu(7,i,j,k)*(1./detg1)
        bcon(2)=alpha1*uu(8,i,j,k)*(1./detg1)
        bcon(3)=alpha1*uu(9,i,j,k)*(1./detg1)
         
        ro=uri(1,i,j,k)
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
!        vcon(1)=uri(2,i,j,k)
!        vcon(2)=uri(3,i,j,k)
!        vcon(3)=uri(4,i,j,k)
        pp=uri(5,i,j,k)
        
        ro1=uri(1,i,j,k)
        util1(1)=uri(2,i,j,k)
        util1(2)=uri(3,i,j,k)
        util1(3)=uri(4,i,j,k)
!        vcon1(1)(i,j,k)=uri(2,i,j,k)
!        vcon1(2)(i,j,k)=uri(3,i,j,k)
!        vcon1(3)(i,j,k)=uri(4,i,j,k)
        pp1=uri(5,i,j,k)
!
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          roh=ro+(gam/(gam-1.0))*pp
        elseif(ieos .eq. 1) then
          roh=(5./2.)*pp &
              +sqrt((9./4.)*pp**2+ro**2)
        elseif(ieos .eq. 2) then
          roe=(3./2.)*(pp+((3.*pp**2) &
              /(2.0*ro+sqrt(2.*pp**2+4.*ro**2)) ))
          roh=ro+roe+pp
        endif

!- cal of covariant 3-magnetic field -! 
        call lower(bcon,bcov,gcov1)
!- cal of contravariant 4-q -!
        call upper(qcov,qcon,gcon1)
!- cal of B-field square (3-vector) -!
        bsq=bcon(1)*bcov(1)+bcon(2)*bcov(2)+bcon(3)*bcov(3)
!- cal of q*B -!
        qdotb=qcov(0)*bcon(0)+qcov(1)*bcon(1)+qcov(2)*bcon(2)+qcov(3)*bcon(3)
!- cal of (Q*B)^2 -!
        qdotbsq=qdotb*qdotb
!- set of covariant 4-vector n -!
        oncov(0)=-1.d0*alpha1
        oncov(1)=0.d0
        oncov(2)=0.d0
        oncov(3)=0.d0
!- cal of contravariant 4-vector n -!
        call upper(oncov,oncon,gcon1)
!- cal of Q*n -!
        qdotn=qcon(0)*oncov(0)
!
!- cal of square of Q -!
        qsq=0.d0
        qsq=qcov(0)*qcon(0)+qcov(1)*qcon(1) &
           +qcov(2)*qcon(2)+qcov(3)*qcon(3)
!- cal of Q^2+Q*n -!
        qtsq=qsq+(qdotn*qdotn)
!
!        qsq1=0.d0
!        qsq1=qcov(1)*qcon(1)+qcov(2)*qcon(2)+qcov(3)*qcon(3)
!
!
!- Check inequaility of updated conserved variables (qdotn) -!
!- set atmosphere value of density and pressure
        romin=dmin*(x1(i)**(-3./2.))
        prmin=pmin*(gam-1.)*(x1(i)**(-5./2.))
!
        if(dro .le. 0.d0) then
          irecov=2
        endif
!
        if(irecov .eq. 0) then
          call inequality(dro,bsq,qtsq,qdotb,qdotn,romin,prmin)
        endif        
!- Calculate W from last timestep and use for inital guess -!      
!
!- calculation of 4-velocity square 
        utsq=0.d0
!        vsq=0.d0
        do m=1,3
          do n=1,3
            utsq=utsq+gcov1(m,n)*util(m)*util(n)
!            vsq=vsq+gcov1(m,n)*util(m)*util(n)
          enddo
        enddo
!
        if(utsq .lt. 0.d0 .and. abs(utsq) .le. safty) then
          utsq=abs(utsq)
        endif
        if(abs(utsq) .gt. 1.d4) then
          utsq=1.d4
        endif 
!- calculation of Lorentz factor -!
        gfsq=1.+utsq
        gf=sqrt(gfsq)
!        gf=1./sqrt(1.-vsq)
!        gfsq=gf**2
        vsq=(gfsq-1.)*(1./gfsq)
!- calculation of W -!      
!        if(vsq .lt. 0.9d0 .and. bsq/roh .lt. 1.d0) then 
! 
!          ro=dro*(1./gf)
!          roh=ro+gam1*pp
!          wg=abs(gfsq*roh)
!
!        else
          
!          if(igue .eq. 1) then
!            vsq1=1.d0
!          else
            vsq1=vsq
!          endif
!            
          tmp1a=4.-vsq1
          tmp1b=4.*(bsq+qdotn)
          tmp1c=qtsq+bsq*bsq+2.*bsq*qdotn
          tmp1d=(tmp1b*tmp1b)-4.*tmp1a*tmp1c
          if(tmp1d .le. 0.d0) then
            tmp1d=abs(tmp1d)
          endif
          wgp=(-tmp1b+sqrt(tmp1d))*(1./(2.*tmp1a))
          wgm=(-tmp1b-sqrt(tmp1d))*(1./(2.*tmp1a))
!
!         if(tmp1d .le. 0.d0) then
!           write(*,*) "inside sqrt < 0 for cal wg", tmp1b, tmp1c
!          endif
!
          if(wgp .gt. 0.d0) then
            wg=wgp
          elseif(wgm .ge. 0.d0) then
            wg=wgm
          endif
          
!        endif
!
! Make sure that W is large enough so that v^2 < 1
!       
95      continue

        if(ieos .eq. 0) then !- gamma-law EoS -!
          pp=gam2*(wg*gf-dro*sqrt(gf))
        elseif(ieos .eq. 3) then !- polytropic EoS -!   
          pp=kpol*((dro*sqrt(gf))**gam)
        endif
       
        tmp3a=(wg*wg*wg)*(wg+2.*bsq)-qdotbsq*(2.*wg+bsq)
        tmp3b=(wg*wg)*(qtsq-bsq*bsq)
!
        tmp3c=(wg*wg*wg)+(wg*wg)*(0.5*bsq+qdotn-pp)-0.5*qdotbsq
!        
!        do nnn=1,iter
          if(tmp3a .le. tmp3b) then
!             write(*,*) "added wg at i,j,k=",i,j,k
             wg= abs(1.5*wg)
             go to 95
!             tmp3a=(wg*wg*wg)*(wg+2.*bsq)-qdotbsq*(2.*wg+bsq)
!             tmp3b=(wg*wg)*(qtsq-bsq*bsq)
          endif
!
!        enddo
!          
        if(tmp3c .gt. 0.d0) then
          wg=abs(0.9*wg)
          go to 95 
        endif   
!
! Calculate W & v^2 for initial guess
!
        wsq=wg*wg
        xsq=(bsq+wg)*(bsq+wg)
        vsq=abs((wsq*qtsq+qdotbsq*(bsq+2.*wg))*(1./(wsq*xsq)))
!       
        if(vsq .gt. 1.d0) then
          write(6,*) 'vsq > 1.0 at initial guess', i, j, k 
          vsq=1.d0-dv
        endif
!
!=====================================================================@
!- Newton Rapson calculation -!

        xxn=wg
        yyn=vsq
        xo=wg
        yo=vsq

!- Newton-Raphson routine start -!
         
        do nnn=1,iter

          if(abs(delta) .ge. deltaend .and. irecov .eq. 0) then
           
            if(nnn .eq. iter) then
              irecov=1

            else
              wg=xxn
              vsq=yyn

              wsq=wg*wg
              w3=wsq*wg
              gtmp=1.-vsq
              if(gtmp .lt. 0.d0 .or. gtmp .lt. safty) then
                gtmp=safty
                vsq=1.-gtmp 
              endif   

              if(ieos .eq. 0) then !- gamma-law EoS -!
                ptmp=gam2*(wg*gtmp-dro*sqrt(gtmp))
              elseif(ieos .eq. 3) then !- polytropic EoS -!   
                ptmp=kpol*((dro*sqrt(gtmp))**gam)
              endif
              
              f1=qtsq-vsq*(wg+bsq)*(wg+bsq)+(qdotbsq*(bsq+2.*wg))*(1./wsq)
              g1=-qdotn-0.5*bsq*(1.+vsq)+0.5*(qdotbsq*(1./wsq))-wg+ptmp

              if(ieos .eq. 0) then !- gamma-law EoS -!
                dpdw=gam2*gtmp
                dpdvsq=-gam2*wg+0.5*gam2*dro*(1./(sqrt(gtmp)))
              elseif(ieos .eq. 3) then !- polytropic EoS -!
                dpdw=0.d0
                dpdvsq=-0.5*kpol*((dro/(sqrt(gtmp)))**gam)  
              endif
              
              dfx=-2.*(wg+bsq)*(vsq+qdotbsq*(1./w3))
              dfy=-((wg+bsq)*(wg+bsq))
              dgx=-1.-(qdotbsq*(1./w3))+dpdw
              dgy=-0.5*bsq+dpdvsq
              
              det=dfx*dgy-dfy*dgx

              if(det .eq. 0.d0) then
                write(*,*) "det=0 in recov1 at", i, j, k
                det=1.d0 
              endif   
              dx=(-dgy*f1+dfy*g1)*(1./det)
              dy=(dgx*f1-dfx*g1)*(1./det)

              df=-f1*f1-g1*g1
              f=-0.5*df
            
              delta=0.5*(abs(dx)+abs(dy))
!              delta=abs(dx)
 
              xo=wg
              yo=vsq
            
              xxn=wg+dx
              yyn=vsq+dy
!
              if(xxn .lt. 0.d0) then
!                irecov=2
                 xxn=abs(xxn)
              endif
!         
              if(xxn .gt. saftybig) then
                xxn=saftybig
              endif
!                         
              if(yyn .lt. 0.d0) then
                yyn=abs(yyn)
              endif  
!
              if(yyn .gt. 1.d0) then
                dv=1.0d-2
                yyn=1.0-dv
              endif
             
           endif
          endif

        enddo
!
!- Newton-Raphson routine end -!
!
        if(irecov .eq. 0) then
          wg=xxn
          vsq=yyn
!
          gtmp=(1.-vsq)
!
          if(gtmp .lt. 0.d0 .or. gtmp .lt. safty) then
            gtmp=safty
          endif 
!
          gf=1./sqrt(gtmp)
          ro=dro/gf
!
          w2=wg*gtmp
!          if(ieos .eq. 3) then !-polytropic EoS
!            pp=kpol*(ro**gam)
!          elseif(ieos .eq. 0) then !- gamma-law EoS 
            pp=gam2*(w2-ro)
!          endif
            
          if(ro .lt. 0.d0) then
            write(*,*) 'negative density', ro, 'at i, j, k', i, j ,k
!            ro=dmin*(x1(i)**(-3./2.))
!            ro=dmin
            ro=romin
            pp=prmin
          endif

!          if(irecov .ne. 3) then
            qtcon(1)=qcon(1)+oncon(1)*qdotn
            qtcon(2)=qcon(2)+oncon(2)*qdotn
            qtcon(3)=qcon(3)+oncon(3)*qdotn

            util(1)=(gf*(1./(wg+bsq)))*(qtcon(1)+qdotb*bcon(1)*(1./wg))
            util(2)=(gf*(1./(wg+bsq)))*(qtcon(2)+qdotb*bcon(2)*(1./wg))
            util(3)=(gf*(1./(wg+bsq)))*(qtcon(3)+qdotb*bcon(3)*(1./wg))

!            vcon(1)=(1./(xxn+bsq))*(qtcon(1)+qdotb*bcon(1)/wg)
!            vcon(2)=(1./(xxn+bsq))*(qtcon(2)+qdotb*bcon(2)/wg)
!            vcon(3)=(1./(xxn+bsq))*(qtcon(3)+qdotb*bcon(3)/wg)
!           endif

          if (pp .lt. 0.d0) then
!            irecov=3
            write(*,*) 'negative pressure', pp, 'at i, j, k', i, j ,k
            if(ieos .eq. 0) then
!              pp=pmin*(gam-1.)*(x1(i)**(-5./2.))
!              pp=pmin 
              pp=prmin
              ro=romin
              util(1)=0.d0
              util(2)=0.d0
              util(3)=0.d0
            elseif(ieos .eq. 3) then
              pp=kpol*(dmin**gam)
            endif   
!            pp=pmin
          endif

        endif
                 
        wg=xxn
        vsq=yyn   

!---------------------------------------------------------------------
!- Error message -!                
        if(irecov .eq. 1) then
          write(6,*) ' >> Not convergence in recov1' &
                    ,' at i, j, k:',i, j, k
          write(6,*) ' >> delta =', delta &
               ,', at x1, x2, x3:',x1(i),x2(j),x3(k)
          if(x1(i) .le. 2.*rbh1) then
            irecov=2
          endif   
         
        elseif(irecov .eq. 2) then
          write(6,*) ' >> wg  < 0 in recov1' &
                    ,' at i, j, k:',i, j, k
          write(6,*) ' >> wg', wg &
                    ,', at x1, x2, x3:',x1(i),x2(j),x3(k)

!         elseif(irecov .eq. 3) then
!           write(6,*) ' >> ro or pr < 0 in recov1d' &
!                    ,' at i, j, k:',i, j, k
!           write(6,*) ' >> ro, pr', ro, pp &
!                     ,' at x1, x2, x3: ',x1(i),x2(j),x3(k)
        endif
!
!---------------------------------------------------------------------@
!   Set Primitive variables
!---------------------------------------------------------------------@
!
        if(irecov .eq. 0) then
         
          uri(1,i,j,k)=ro
          uri(2,i,j,k)=util(1)
          uri(3,i,j,k)=util(2)
          uri(4,i,j,k)=util(3)
!          uri(2,i,j,k)=vcon(1)
!          uri(3,i,j,k)=vcon(2)
!          uri(4,i,j,k)=vcon(3)
          uri(5,i,j,k)=pp
          uri(6,i,j,k)=pp/(ro**(gam-1.))
          uri(7,i,j,k)=uu(7,i,j,k)*(1./detg1)
          uri(8,i,j,k)=uu(8,i,j,k)*(1./detg1)
          uri(9,i,j,k)=uu(9,i,j,k)*(1./detg1)
!
        elseif(irecov .eq. 2) then !- dro < 0 region
!          
          uri(1,i,j,k)=romin
          uri(2,i,j,k)=0.d0
          uri(3,i,j,k)=0.d0
          uri(4,i,j,k)=0.d0
          uri(5,i,j,k)=prmin
          uri(6,i,j,k)=prmin/(romin**(gam-1.))
          uri(7,i,j,k)=uu(7,i,j,k)*(1./detg1)
          uri(8,i,j,k)=uu(8,i,j,k)*(1./detg1)
          uri(9,i,j,k)=uu(9,i,j,k)*(1./detg1)
          
        else
         
          uri(1,i,j,k)=ro1
          uri(2,i,j,k)=util1(1)
          uri(3,i,j,k)=util1(2)
          uri(4,i,j,k)=util1(3)
!          uri(2,i,j,k)=vcon1(1)
!          uri(3,i,j,k)=vcon1(2)
!          uri(4,i,j,k)=vcon1(3)
          uri(5,i,j,k)=pp1
          uri(6,i,j,k)=pp1/(ro1**(gam-1.)) 
          uri(7,i,j,k)=uu(7,i,j,k)*(1./detg1)
          uri(8,i,j,k)=uu(8,i,j,k)*(1./detg1)
          uri(9,i,j,k)=uu(9,i,j,k)*(1./detg1)
         
        endif
!
      enddo
    enddo
  enddo
!
!=====================================================================@
!
  deallocate( gcov1, gcon1, qcon, qcov, bcon, bcov, &
              oncov, oncon, util, util1, qtcon, stat=merr )
!
  return
end subroutine recov1
!
!---------------------------------------------------------------------@
subroutine recov1a(uu,uri,gcov,gcon,detg,x1,x2,x3,nm1,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!- Following Noble et al. (2006), 2-variables inversion procedure -!
!- Switch to entropy inversion procedure (v^2) when it is failed -!  
!- (only ideal EoS)
  
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos, iter, dmin, pmin, &
                   kpol, akm, pi
  implicit none
!
  integer :: i, j, k, m, n, nnn, nm1, is1, ie1, js1, je1, ks1, ke1
  integer :: merr, igue, irecov

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)

  real(8) :: x1(imax), x2(jmax), x3(kmax)

  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: qcon(:), qcov(:), bcon(:), bcov(:), &
                          oncov(:), oncon(:) 
  real(8), allocatable :: util(:), util1(:), qtcon(:), ucon(:)
      
  real(8) :: safty, saftybig, gam1, gam2, dv, alpha1, detg1 
  real(8) :: dro, ro, pp, ro1, pp1, roh, roe, bsq, qdotb, qdotbsq, qdotn, &
             qsq, qtsq, utsq, gfsq, gf, vsq, wg, wsq, xsq, qsq1, &
             delta, deltaend, df, f, xxn, yyn, xo, yo, w2, w3, gtmp, ptmp, &
             f1, g1, dfx, dfy, dgx, dgy, det, dx, dy, dpdw, dpdvsq, &
             dgfdvsq, dwdvsq, gss, ss, ss1, wgo, vsqo
  real(8) :: tmp1a, tmp1b, tmp1c, tmp1d, wgp, wgm, vsq1, tmp3a, tmp3b,  &
             tmp3c, tmp4a
  real(8) :: rin, rbh1, rr, th
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), stat=merr)
  allocate(qcon(0:3), qcov(0:3), bcon(0:3), bcov(0:3), &
           oncov(0:3), oncon(0:3), util(1:3), util1(1:3), qtcon(1:3), &
           ucon(0:3), stat=merr)
!
!--- define variable of safe number ----
      safty = 1.0d-4
      saftybig = 1.0d5
      dv=1.0d-2
      
      gam1=gam/(gam-1.0)
      gam2=1.0/gam1
      wg=0.d0

      igue=1
!
      rin=15.d0
      rbh1=1.+sqrt(1.-akm*akm)      
!
!      
!=====================================================================@
!
! ----  Calculation of variables ---
!
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
      do i=is1+nm1,ie1-nm1

!- set parameter for Newton Raphson -!

        delta=1.d0
        deltaend = 1.d-7
        df=1.d0
        f=1.d0
!        irecov=0

!- check the position within the BH radius -!
!
!        if(x1(i) .le. rbh1) then
!          irecov=2
!        else
          irecov=0
!        endif  
!      
!- copy of metric terms -!
!
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo
!
        alpha1=1./sqrt(-gcon1(0,0))
!
        detg1=detg(i,j,k)
!
!- rearrange of conserved variables -!       
!
        dro=alpha1*uu(1,i,j,k)*(1./detg1)

!        qcov(0)=alpha1*(uu(5,i,j,k)+uu(1,i,j,k))*(1./detg1)
        qcov(0)=alpha1*(uu(5,i,j,k)-uu(1,i,j,k))*(1./detg1)
!        qcov(0)=alpha1*uu(5,i,j,k)*(1./detg1)
        qcov(1)=alpha1*uu(2,i,j,k)*(1./detg1)
        qcov(2)=alpha1*uu(3,i,j,k)*(1./detg1)
        qcov(3)=alpha1*uu(4,i,j,k)*(1./detg1)

        gss=alpha1*uu(6,i,j,k)*(1./detg1)

        tmp4a=alpha1*uu(5,i,j,k)*(1./detg1)
        
        bcon(0)=0.d0
        bcon(1)=alpha1*uu(7,i,j,k)*(1./detg1)
        bcon(2)=alpha1*uu(8,i,j,k)*(1./detg1)
        bcon(3)=alpha1*uu(9,i,j,k)*(1./detg1)
!        
!- copy of primitive variables (previous time)         
!
        ro=uri(1,i,j,k)
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
!        vcon(1)=uri(2,i,j,k)
!        vcon(2)=uri(3,i,j,k)
!        vcon(3)=uri(4,i,j,k)
        pp=uri(5,i,j,k)
        ss=uri(6,i,j,k)
!
!- store for previous time        
!
        ro1=uri(1,i,j,k)
        util1(1)=uri(2,i,j,k)
        util1(2)=uri(3,i,j,k)
        util1(3)=uri(4,i,j,k)
!        vcon1(1)(i,j,k)=uri(2,i,j,k)
!        vcon1(2)(i,j,k)=uri(3,i,j,k)
!        vcon1(3)(i,j,k)=uri(4,i,j,k)
        pp1=uri(5,i,j,k)
        ss1=uri(6,i,j,k)
!
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          roh=ro+(gam/(gam-1.0))*pp
        elseif(ieos .eq. 1) then
          roh=(5./2.)*pp &
              +sqrt((9./4.)*pp**2+ro**2)
        elseif(ieos .eq. 2) then
          roe=(3./2.)*(pp+((3.*pp**2) &
              /(2.0*ro+sqrt(2.*pp**2+4.*ro**2)) ))
          roh=ro+roe+pp
        endif

!- cal of covariant 3-magnetic field -! 
        call lower(bcon,bcov,gcov1)
!- cal of contravariant 4-q -!
        call upper(qcov,qcon,gcon1)
!- cal of B-field square (3-vector) -!
        bsq=bcon(1)*bcov(1)+bcon(2)*bcov(2)+bcon(3)*bcov(3)
!- cal of q*B -!
        qdotb=qcov(0)*bcon(0)+qcov(1)*bcon(1)+qcov(2)*bcon(2)+qcov(3)*bcon(3)
!- cal of (Q*B)^2 -!
        qdotbsq=qdotb*qdotb
!- set of covariant 4-vector n -!
        oncov(0)=-1.d0*alpha1
        oncov(1)=0.d0
        oncov(2)=0.d0
        oncov(3)=0.d0
!- cal of contravariant 4-vector n -!
        call upper(oncov,oncon,gcon1)
!- cal of Q*n -!
        qdotn=qcon(0)*oncov(0)
!
!- cal of square of Q -!
        qsq=0.d0
        qsq=qcov(0)*qcon(0)+qcov(1)*qcon(1) &
           +qcov(2)*qcon(2)+qcov(3)*qcon(3)
!- cal of Q^2+Q*n -!
        qtsq=qsq+(qdotn*qdotn)
!
!        qsq1=0.d0
!        qsq1=qcov(1)*qcon(1)+qcov(2)*qcon(2)+qcov(3)*qcon(3)
!
!
!- Calculate W from last timestep and use for inital guess -!      
!
!- calculation of 4-velocity square 
        utsq=0.d0
!        vsq=0.d0
        do m=1,3
          do n=1,3
            utsq=utsq+gcov1(m,n)*util(m)*util(n)
!            vsq=vsq+gcov1(m,n)*util(m)*util(n)
          enddo
        enddo
!
        if(utsq .lt. 0.d0 .and. abs(utsq) .le. safty) then
          utsq=abs(utsq)
        endif
        if(abs(utsq) .gt. 1.d4) then
          utsq=1.d4
        endif 
!- calculation of Lorentz factor -!
        gfsq=1.+utsq
        gf=sqrt(gfsq)
!        gf=1./sqrt(1.-vsq)
!        gfsq=gf**2
        vsq=(gfsq-1.)*(1./gfsq)
!- calculation of W -!      
!        if(vsq .lt. 0.9d0 .and. bsq/roh .lt. 1.d0) then 
! 
!          ro=dro*(1./gf)
!          roh=ro+gam1*pp
!          wg=abs(gfsq*roh)
!
!        else
          
!          if(igue .eq. 1) then
!            vsq1=1.d0
!          else
            vsq1=vsq
!          endif
!            
          tmp1a=4.-vsq1
          tmp1b=4.*(bsq+qdotn)
          tmp1c=qtsq+bsq*bsq+2.*bsq*qdotn
          tmp1d=(tmp1b*tmp1b)-4.*tmp1a*tmp1c
          if(tmp1d .le. 0.d0) then
            tmp1d=abs(tmp1d)
          endif
          wgp=(-tmp1b+sqrt(tmp1d))*(1./(2.*tmp1a))
          wgm=(-tmp1b-sqrt(tmp1d))*(1./(2.*tmp1a))
!
         if(tmp1d .le. 0.d0) then
           write(*,*) "inside sqrt < 0 for cal wg", tmp1b, tmp1c
          endif
!
          if(wgp .gt. 0.d0) then
            wg=wgp
          elseif(wgm .ge. 0.d0) then
            wg=wgm
          endif
          
!        endif
!
! Make sure that W is large enough so that v^2 < 1
!       
95      continue

        if(ieos .eq. 0) then !- gamma-law EoS -!
          pp=gam2*(wg*gf-dro*sqrt(gf))
        elseif(ieos .eq. 3) then !- polytropic EoS -!   
          pp=kpol*((dro*sqrt(gf))**gam)
        endif
       
        tmp3a=(wg*wg*wg)*(wg+2.*bsq)-qdotbsq*(2.*wg+bsq)
        tmp3b=(wg*wg)*(qtsq-bsq*bsq)
!
        tmp3c=(wg*wg*wg)+(wg*wg)*(0.5*bsq+qdotn-pp)-0.5*qdotbsq
!        
!        do nnn=1,iter
          if(tmp3a .le. tmp3b) then
!             write(*,*) "added wg at i,j,k=",i,j,k
             wg= abs(1.5*wg)
             go to 95
!             tmp3a=(wg*wg*wg)*(wg+2.*bsq)-qdotbsq*(2.*wg+bsq)
!             tmp3b=(wg*wg)*(qtsq-bsq*bsq)
          endif

!        enddo
 
        if(tmp3c .gt. 0.d0) then
          wg=abs(0.9*wg)
          go to 95 
        endif            
!
! Calculate W & v^2 for initial guess
!
        wsq=wg*wg
        xsq=(bsq+wg)*(bsq+wg)
        vsq=abs((wsq*qtsq+qdotbsq*(bsq+2.*wg))*(1./(wsq*xsq)))
!       
        if(vsq .gt. 1.d0) then
          write(6,*) 'vsq > 1.0 at initial guess', i, j, k 
          vsq=1.d0-dv
        endif
!
!- Obtained initial guess -!
        wgo=wg
        vsqo=vsq
!         
!=====================================================================@
!- Newton Rapson calculation -!

        xxn=wg
        yyn=vsq
        xo=wg
        yo=vsq

!- Newton-Raphson routine start -!
         
        do nnn=1,iter

          if(abs(delta) .ge. deltaend .and. irecov .eq. 0) then
           
            if(nnn .eq. iter) then
              irecov=1

            else
              wg=xxn
              vsq=yyn

              wsq=wg*wg
              w3=wsq*wg
              gtmp=1.-vsq
              if(gtmp .lt. 0.d0 .or. gtmp .lt. safty) then
                gtmp=safty
                vsq=1.-gtmp 
              endif      
              
              if(ieos .eq. 0) then !- gamma-law EoS -!
                ptmp=gam2*(wg*gtmp-dro*sqrt(gtmp))
              elseif(ieos .eq. 3) then !- polytropic EoS -!   
                ptmp=kpol*(dro**gam)*sqrt(gtmp)
              endif
              
              f1=qtsq-vsq*(wg+bsq)*(wg+bsq)+(qdotbsq*(bsq+2.*wg))*(1./wsq)
              g1=-qdotn-0.5*bsq*(1.+vsq)+0.5*(qdotbsq*(1./wsq))-wg+ptmp

              if(ieos .eq. 0) then !- gamma-law EoS -!
                dpdw=gam2*gtmp
                dpdvsq=-gam2*wg+0.5*gam2*dro*(1./(sqrt(gtmp)))
              elseif(ieos .eq. 3) then !- polytropic EoS -!
                dpdw=0.d0
                dpdvsq=-0.5*kpol*(dro**gam)*(1./(sqrt(gtmp)))  
              endif
              
              dfx=-2.*(wg+bsq)*(vsq+qdotbsq*(1./w3))
              dfy=-((wg+bsq)*(wg+bsq))
              dgx=-1.-(qdotbsq*(1./w3))+dpdw
              dgy=-0.5*bsq+dpdvsq
              
              det=dfx*dgy-dfy*dgx

              if(det .eq. 0.d0) then
                write(*,*) "det=0 in recov1 at", i, j, k
                det=1.d0 
              endif   
              dx=(-dgy*f1+dfy*g1)*(1./det)
              dy=(dgx*f1-dfx*g1)*(1./det)

              df=-f1*f1-g1*g1
              f=-0.5*df
            
              delta=0.5*(abs(dx)+abs(dy))
!              delta=abs(dx)
 
              xo=wg
              yo=vsq
            
              xxn=wg+dx
              yyn=vsq+dy
!
              if(xxn .lt. 0.d0) then
!                irecov=2
                 xxn=abs(xxn)
              endif
!         
              if(xxn .gt. saftybig) then
                xxn=saftybig
              endif
!            
              if(yyn .lt. 0.d0) then
                yyn=abs(yyn)
              endif
!             
              if(abs(yyn) .gt. 1.d0) then
                dv=1.0d-2
                yyn=1.0-dv
              endif
 !            
            endif
          endif

        enddo
!
!- Newton-Raphson routine end -!
!
        if(irecov .eq. 0) then
          wg=xxn
          vsq=yyn
!
          gtmp=1.-vsq

          if(gtmp .lt. 0.d0 .or. gtmp .lt. safty) then
            gtmp=safty
            vsq=1.-gtmp 
          endif
           
          gf=1./sqrt(gtmp)
          ro=dro*(1./gf)
!
          w2=wg*(gtmp)
!          if(ieos .eq. 3) then !-polytropic EoS
!            pp=kpol*(ro**gam)
!          elseif(ieos .eq. 0) then !- gamma-law EoS 
            pp=gam2*(w2-ro)
!          endif

          if(ro .lt. 0.d0 .or. pp .lt. 0.d0) then  !- check pressure & density -!           
!          if(ro .lt. 0.d0) then
!            write(*,*) 'negative density', ro, 'at i, j, k', i, j ,k
!            ro=dmin*(x1(i)**(-3./2.))
!!            ro=dmin
!          endif
            irecov=1
          endif  
!          if (pp .lt. 0.d0) then
!!            irecov=3
!            write(*,*) 'negative pressure', pp, 'at i, j, k', i, j ,k
!            pp=(gam-1.)*pmin*(x1(i)**(-5./2.))
!!            pp=pmin
!          endif

          if(irecov .eq. 0) then !- velocity, pressure, & density all clear -! 
            qtcon(1)=qcon(1)+oncon(1)*qdotn
            qtcon(2)=qcon(2)+oncon(2)*qdotn
            qtcon(3)=qcon(3)+oncon(3)*qdotn

            util(1)=(gf*(1./(wg+bsq)))*(qtcon(1)+qdotb*bcon(1)*(1./wg))
            util(2)=(gf*(1./(wg+bsq)))*(qtcon(2)+qdotb*bcon(2)*(1./wg))
            util(3)=(gf*(1./(wg+bsq)))*(qtcon(3)+qdotb*bcon(3)*(1./wg))

!            vcon(1)=(1./(xxn+bsq))*(qtcon(1)+qdotb*bcon(1)/wg)
!            vcon(2)=(1./(xxn+bsq))*(qtcon(2)+qdotb*bcon(2)/wg)
!            vcon(3)=(1./(xxn+bsq))*(qtcon(3)+qdotb*bcon(3)/wg)
          endif
                 
          wg=xxn
          vsq=yyn   
!
        endif       
!
!=====================================================================@        
!- Another Newton Raphson Calculation -!

        if(irecov .eq. 1) then
           
!- set parameter  -!

          delta=1.d0
          deltaend = 1.d-6
          df=1.d0
          f=1.d0
          irecov=0

!          xxn=wgo
!          xo=wgo
          xxn=vsqo
          xo=vsqo
          
!- Newton-Raphson routine start -!
         
          do nnn=1,iter

            if(abs(delta) .ge. deltaend .and. irecov .eq. 0) then
           
              if(nnn .eq. iter) then
                irecov=1

              else
                vsq=xxn
!               
!- calculation of W -!
!
                gtmp=1.-vsq
                if(gtmp .lt. 0.d0 .or. gtmp .lt. safty) then
                  gtmp=safty
                  vsq=1.-gtmp 
                endif
            
                gf=1./sqrt(gtmp)
                gfsq=gf*gf
              
                ro=dro*(1./gf)
                pp=(gss/dro)*(ro**(gam))
                roh=ro+gam1*pp
                wg=abs(gfsq*roh)
              
                wsq=wg*wg
                w3=wsq*wg
!
!-calculation of dW/dv^2 -!
!
                dgfdvsq=0.5*(gtmp)**(-3./2.)
!                dwdvsq=2.*gf*roh*dgfdvsq
                dwdvsq=dro*dgfdvsq &
                     +gam1*gss*(dro**(gam-1.)) &
                     *(2.-gam)*(gf**(1.-gam))*dgfdvsq
!              
!- calculation of f and df -!
!              
                f=qtsq-vsq*(wg+bsq)*(wg+bsq)+(qdotbsq*(bsq+2.*wg))/wsq
                tmp1a=-(wg+bsq)*(wg+bsq)
                tmp1b=-2.*vsq*(wg+bsq)
                tmp1c=-2.*qdotbsq*(bsq+2.*wg)*(1./w3)
                tmp1d=2.*qdotbsq*(1./wsq)
                df=tmp1a+(tmp1b+tmp1c+tmp1d)*dwdvsq
!
!- calculation of delta -!
!              
                delta=f/df
! 
                xo=vsq
                xxn=vsq-delta
!
                if(xxn .lt. 0.d0) then
                  xxn=abs(xxn)
                endif
!               
                if(abs(xxn) .ge. 1.d0) then
!                  irecov=2
                  xxn=1.-dv
                endif 
!              
              endif
            endif

          enddo
!
!- Newton-Raphson routine end -!
!
          if(irecov .eq. 0) then
            vsq=xxn
            gtmp=1.-vsq

            if(gtmp.lt. 0.d0 .or. gtmp .lt. safty) then
              gtmp=safty 
            endif
         
            gf=1./sqrt(gtmp)
            gfsq=gf*gf
              
            ro=dro*(1./gf)
            pp=(gss*(1./gf))*(ro**(gam-1.))
            roh=ro+gam1*pp

            if(ro .lt. 0.d0) then
              write(*,*) 'negative density', ro, 'at i, j, k', i, j ,k
              ro=dmin*(x1(i)**(-3./2.))
!!            ro=dmin
            endif
            if (pp .lt. 0.d0) then
!!            irecov=3
              write(*,*) 'negative pressure', pp, 'at i, j, k', i, j ,k
              pp=(gam-1.)*pmin*(x1(i)**(-5./2.))
!!            pp=pmin
            endif

            
            wg=abs(gfsq*roh)

            qtcon(1)=qcon(1)+oncon(1)*qdotn
            qtcon(2)=qcon(2)+oncon(2)*qdotn
            qtcon(3)=qcon(3)+oncon(3)*qdotn

            util(1)=(gf/(wg+bsq))*(qtcon(1)+qdotb*bcon(1)*(1./wg))
            util(2)=(gf/(wg+bsq))*(qtcon(2)+qdotb*bcon(2)*(1./wg))
            util(3)=(gf/(wg+bsq))*(qtcon(3)+qdotb*bcon(3)*(1./wg))
          endif
         
        endif
        
!        wg=xxn  
!
!---------------------------------------------------------------------
!- Error message -!                
        if(irecov .eq. 1) then
          write(6,*) ' >> Not convergence in recov1-3' &
                    ,' at i, j, k:',i, j, k
          write(6,*) ' >> delta =', delta &
                    ,', at x1, x2, x3:',x1(i),x2(j),x3(k)
         
        elseif(irecov .eq. 2) then
          write(6,*) ' >> wg  < 0 in recov1-3' &
                    ,' at i, j, k:',i, j, k
          write(6,*) ' >> wg', wg &
                    ,', at x1, x2, x3:',x1(i),x2(j),x3(k)

!         elseif(irecov .eq. 3) then
!           write(6,*) ' >> ro or pr < 0 in recov1d' &
!                    ,' at i, j, k:',i, j, k
!           write(6,*) ' >> ro, pr', ro, pp &
!                     ,' at x1, x2, x3: ',x1(i),x2(j),x3(k)
        endif
!
!---------------------------------------------------------------------@
!   Set Primitive variables
!---------------------------------------------------------------------@
!
        if(irecov .eq. 0) then
         
          uri(1,i,j,k)=ro
          uri(2,i,j,k)=util(1)
          uri(3,i,j,k)=util(2)
          uri(4,i,j,k)=util(3)
!          uri(2,i,j,k)=vcon(1)
!          uri(3,i,j,k)=vcon(2)
!          uri(4,i,j,k)=vcon(3)
          uri(5,i,j,k)=pp
          uri(6,i,j,k)=ro*pp/(ro**(gam-1.))
          uri(7,i,j,k)=uu(7,i,j,k)*(1./detg1)
          uri(8,i,j,k)=uu(8,i,j,k)*(1./detg1)
          uri(9,i,j,k)=uu(9,i,j,k)*(1./detg1)
!        
!- cal of contravariant 4-velocity -!
          call calucon(util,ucon,gf,gcon1)
!- update conserved enthalpy -!
          uu(6,i,j,k)=detg1*uri(6,i,j,k)*ucon(0)
!          
       else
         
          uri(1,i,j,k)=ro1
          uri(2,i,j,k)=util1(1)
          uri(3,i,j,k)=util1(2)
          uri(4,i,j,k)=util1(3)
!          uri(2,i,j,k)=vcon1(1)
!          uri(3,i,j,k)=vcon1(2)
!          uri(4,i,j,k)=vcon1(3)
          uri(5,i,j,k)=pp1
          uri(6,i,j,k)=ro1*pp1/(ro1**(gam)) 
          uri(7,i,j,k)=uu(7,i,j,k)*(1./detg1)
          uri(8,i,j,k)=uu(8,i,j,k)*(1./detg1)
          uri(9,i,j,k)=uu(9,i,j,k)*(1./detg1)
!       
!- cal of contravariant 4-velocity -!
          call calucon(util,ucon,gf,gcon1)
!- update conserved enthalpy -!
          uu(6,i,j,k)=detg1*ro1*uri(6,i,j,k)*ucon(0)
!          
        endif
!
      enddo
    enddo
  enddo
!
!=====================================================================@
!
  deallocate( gcov1, gcon1, qcon, qcov, bcon, bcov, &
              oncov, oncon, util, util1, qtcon, ucon, stat=merr )
!
  return
end subroutine recov1a
!
!---------------------------------------------------------------------@
subroutine recov1b(uu,uri,gcov,gcon,detg,x1,x2,x3,nm1,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!- Following Noble et al. (2006), 2-variables inversion procedure -!
!- Switch to entropy inversion procedure (W) when it is failed -!  
!- (only ideal EoS)
  
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos, iter, dmin, pmin, &
                   kpol, akm, pi
  implicit none
!
  integer :: i, j, k, m, n, nnn, nm1, is1, ie1, js1, je1, ks1, ke1
  integer :: merr, igue, irecov

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)

  real(8) :: x1(imax), x2(jmax), x3(kmax)

  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: qcon(:), qcov(:), bcon(:), bcov(:), &
                          oncov(:), oncon(:) 
  real(8), allocatable :: util(:), util1(:), qtcon(:), ucon(:)
      
  real(8) :: safty, saftybig, gam1, gam2, dv, alpha1, detg1 
  real(8) :: dro, ro, pp, ro1, pp1, roh, roe, bsq, qdotb, qdotbsq, qdotn, &
             qsq, qtsq, utsq, gfsq, gf, vsq, wg, wsq, xsq, qsq1, &
             delta, deltaend, df, f, xxn, yyn, xo, yo, w2, w3, gtmp, ptmp, &
             f1, g1, dfx, dfy, dgx, dgy, det, dx, dy, dpdw, dpdvsq, &
             gss, ss, ss1, wgo, vsqo, dpdw2, dvsqdw, drodw
  real(8) :: tmp1a, tmp1b, tmp1c, tmp1d, wgp, wgm, vsq1, tmp3a, tmp3b,  &
             tmp3c, tmp4a
  real(8) :: rin, rbh1, rr, th
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), stat=merr)
  allocate(qcon(0:3), qcov(0:3), bcon(0:3), bcov(0:3), &
           oncov(0:3), oncon(0:3), util(1:3), util1(1:3), qtcon(1:3), &
           ucon(0:3), stat=merr)
!
!--- define variable of safe number ----
      safty = 1.0d-4
      saftybig = 1.0d5
      dv=1.0d-2
      
      gam1=gam/(gam-1.0)
      gam2=1.0/gam1
      wg=0.d0

      igue=1
!
      rin=15.d0
      rbh1=1.+sqrt(1.-akm*akm)      
!
!      
!=====================================================================@
!
! ----  Calculation of variables ---
!
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
      do i=is1+nm1,ie1-nm1

!- set parameter for Newton Raphson -!

        delta=1.d0
        deltaend = 1.d-7
        df=1.d0
        f=1.d0
!        irecov=0

!- check the position within the BH radius -!
!
!        if(x1(i) .le. rbh1) then
!          irecov=2
!        else
          irecov=0
!        endif  
!      
!- copy of metric terms -!
!
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo
!
        alpha1=1./sqrt(-gcon1(0,0))
!
        detg1=detg(i,j,k)
!
!- rearrange of conserved variables -!       
!
        dro=alpha1*uu(1,i,j,k)*(1./detg1)

!        qcov(0)=alpha1*(uu(5,i,j,k)+uu(1,i,j,k))*(1./detg1)
        qcov(0)=alpha1*(uu(5,i,j,k)-uu(1,i,j,k))*(1./detg1)
!        qcov(0)=alpha1*uu(5,i,j,k)*(1./detg1)
        qcov(1)=alpha1*uu(2,i,j,k)*(1./detg1)
        qcov(2)=alpha1*uu(3,i,j,k)*(1./detg1)
        qcov(3)=alpha1*uu(4,i,j,k)*(1./detg1)

        gss=alpha1*uu(6,i,j,k)*(1./detg1)

        tmp4a=alpha1*uu(5,i,j,k)*(1./detg1)
        
        bcon(0)=0.d0
        bcon(1)=alpha1*uu(7,i,j,k)*(1./detg1)
        bcon(2)=alpha1*uu(8,i,j,k)*(1./detg1)
        bcon(3)=alpha1*uu(9,i,j,k)*(1./detg1)
!        
!- copy of primitive variables (previous time)         
!
        ro=uri(1,i,j,k)
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
!        vcon(1)=uri(2,i,j,k)
!        vcon(2)=uri(3,i,j,k)
!        vcon(3)=uri(4,i,j,k)
        pp=uri(5,i,j,k)
        ss=uri(6,i,j,k)
!
!- store for previous time        
!
        ro1=uri(1,i,j,k)
        util1(1)=uri(2,i,j,k)
        util1(2)=uri(3,i,j,k)
        util1(3)=uri(4,i,j,k)
!        vcon1(1)(i,j,k)=uri(2,i,j,k)
!        vcon1(2)(i,j,k)=uri(3,i,j,k)
!        vcon1(3)(i,j,k)=uri(4,i,j,k)
        pp1=uri(5,i,j,k)
        ss1=uri(6,i,j,k)
!
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          roh=ro+(gam/(gam-1.0))*pp
        elseif(ieos .eq. 1) then
          roh=(5./2.)*pp &
              +sqrt((9./4.)*pp**2+ro**2)
        elseif(ieos .eq. 2) then
          roe=(3./2.)*(pp+((3.*pp**2) &
              /(2.0*ro+sqrt(2.*pp**2+4.*ro**2)) ))
          roh=ro+roe+pp
        endif

!- cal of covariant 3-magnetic field -! 
        call lower(bcon,bcov,gcov1)
!- cal of contravariant 4-q -!
        call upper(qcov,qcon,gcon1)
!- cal of B-field square (3-vector) -!
        bsq=bcon(1)*bcov(1)+bcon(2)*bcov(2)+bcon(3)*bcov(3)
!- cal of q*B -!
        qdotb=qcov(0)*bcon(0)+qcov(1)*bcon(1)+qcov(2)*bcon(2)+qcov(3)*bcon(3)
!- cal of (Q*B)^2 -!
        qdotbsq=qdotb*qdotb
!- set of covariant 4-vector n -!
        oncov(0)=-1.d0*alpha1
        oncov(1)=0.d0
        oncov(2)=0.d0
        oncov(3)=0.d0
!- cal of contravariant 4-vector n -!
        call upper(oncov,oncon,gcon1)
!- cal of Q*n -!
        qdotn=qcon(0)*oncov(0)
!
!- cal of square of Q -!
        qsq=0.d0
        qsq=qcov(0)*qcon(0)+qcov(1)*qcon(1) &
           +qcov(2)*qcon(2)+qcov(3)*qcon(3)
!- cal of Q^2+Q*n -!
        qtsq=qsq+(qdotn*qdotn)
!
!        qsq1=0.d0
!        qsq1=qcov(1)*qcon(1)+qcov(2)*qcon(2)+qcov(3)*qcon(3)
!
!
!- Calculate W from last timestep and use for inital guess -!      
!
!- calculation of 4-velocity square 
        utsq=0.d0
!        vsq=0.d0
        do m=1,3
          do n=1,3
            utsq=utsq+gcov1(m,n)*util(m)*util(n)
!            vsq=vsq+gcov1(m,n)*util(m)*util(n)
          enddo
        enddo
!
        if(utsq .lt. 0.d0 .and. abs(utsq) .le. safty) then
          utsq=abs(utsq)
        endif
        if(abs(utsq) .gt. 1.d4) then
          utsq=1.d4
        endif 
!- calculation of Lorentz factor -!
        gfsq=1.+utsq
        gf=sqrt(gfsq)
!        gf=1./sqrt(1.-vsq)
!        gfsq=gf**2
        vsq=(gfsq-1.)*(1./gfsq)
!- calculation of W -!      
!        if(vsq .lt. 0.9d0 .and. bsq/roh .lt. 1.d0) then 
! 
!          ro=dro*(1./gf)
!          roh=ro+gam1*pp
!          wg=abs(gfsq*roh)
!
!        else
          
!          if(igue .eq. 1) then
!            vsq1=1.d0
!          else
            vsq1=vsq
!          endif
!            
          tmp1a=4.-vsq1
          tmp1b=4.*(bsq+qdotn)
          tmp1c=qtsq+bsq*bsq+2.*bsq*qdotn
          tmp1d=(tmp1b*tmp1b)-4.*tmp1a*tmp1c
          if(tmp1d .le. 0.d0) then
            tmp1d=abs(tmp1d)
          endif
          wgp=(-tmp1b+sqrt(tmp1d))*(1./(2.*tmp1a))
          wgm=(-tmp1b-sqrt(tmp1d))*(1./(2.*tmp1a))
!
         if(tmp1d .le. 0.d0) then
           write(*,*) "inside sqrt < 0 for cal wg", tmp1b, tmp1c
          endif
!
          if(wgp .gt. 0.d0) then
            wg=wgp
          elseif(wgm .ge. 0.d0) then
            wg=wgm
          endif
          
!        endif
!
! Make sure that W is large enough so that v^2 < 1
!       
95      continue

        if(ieos .eq. 0) then !- gamma-law EoS -!
          pp=gam2*(wg*gf-dro*sqrt(gf))
        elseif(ieos .eq. 3) then !- polytropic EoS -!   
          pp=kpol*((dro*sqrt(gf))**gam)
        endif
       
        tmp3a=(wg*wg*wg)*(wg+2.*bsq)-qdotbsq*(2.*wg+bsq)
        tmp3b=(wg*wg)*(qtsq-bsq*bsq)
!
        tmp3c=(wg*wg*wg)+(wg*wg)*(0.5*bsq+qdotn-pp)-0.5*qdotbsq
!        
!        do nnn=1,iter
          if(tmp3a .le. tmp3b) then
!             write(*,*) "added wg at i,j,k=",i,j,k
             wg= abs(1.5*wg)
             go to 95
!             tmp3a=(wg*wg*wg)*(wg+2.*bsq)-qdotbsq*(2.*wg+bsq)
!             tmp3b=(wg*wg)*(qtsq-bsq*bsq)
          endif

!        enddo
 
        if(tmp3c .gt. 0.d0) then
          wg=abs(0.9*wg)
          go to 95 
        endif            
!
! Calculate W & v^2 for initial guess
!
        wsq=wg*wg
        xsq=(bsq+wg)*(bsq+wg)
        vsq=abs((wsq*qtsq+qdotbsq*(bsq+2.*wg))*(1./(wsq*xsq)))
!       
        if(vsq .gt. 1.d0) then
          write(6,*) 'vsq > 1.0 at initial guess', i, j, k 
          vsq=1.d0-dv
        endif
!
!- Obtained initial guess -!
        wgo=wg
        vsqo=vsq
!         
!=====================================================================@
!- Newton Rapson calculation -!

        xxn=wg
        yyn=vsq
        xo=wg
        yo=vsq

!- Newton-Raphson routine start -!
         
        do nnn=1,iter

          if(abs(delta) .ge. deltaend .and. irecov .eq. 0) then
           
            if(nnn .eq. iter) then
              irecov=1

            else
              wg=xxn
              vsq=yyn

              wsq=wg*wg
              w3=wsq*wg
              gtmp=1.-vsq
              if(gtmp .lt. 0.d0 .or. gtmp .lt. safty) then
                gtmp=safty
                vsq=1.-gtmp 
              endif      
              
              if(ieos .eq. 0) then !- gamma-law EoS -!
                ptmp=gam2*(wg*gtmp-dro*sqrt(gtmp))
              elseif(ieos .eq. 3) then !- polytropic EoS -!   
                ptmp=kpol*(dro**gam)*sqrt(gtmp)
              endif
              
              f1=qtsq-vsq*(wg+bsq)*(wg+bsq)+(qdotbsq*(bsq+2.*wg))*(1./wsq)
              g1=-qdotn-0.5*bsq*(1.+vsq)+0.5*(qdotbsq*(1./wsq))-wg+ptmp

              if(ieos .eq. 0) then !- gamma-law EoS -!
                dpdw=gam2*gtmp
                dpdvsq=-gam2*wg+0.5*gam2*dro*(1./(sqrt(gtmp)))
              elseif(ieos .eq. 3) then !- polytropic EoS -!
                dpdw=0.d0
                dpdvsq=-0.5*kpol*(dro**gam)*(1./(sqrt(gtmp)))  
              endif
              
              dfx=-2.*(wg+bsq)*(vsq+qdotbsq*(1./w3))
              dfy=-((wg+bsq)*(wg+bsq))
              dgx=-1.-(qdotbsq*(1./w3))+dpdw
              dgy=-0.5*bsq+dpdvsq
              
              det=dfx*dgy-dfy*dgx

              if(det .eq. 0.d0) then
                write(*,*) "det=0 in recov1 at", i, j, k
                det=1.d0 
              endif   
              dx=(-dgy*f1+dfy*g1)*(1./det)
              dy=(dgx*f1-dfx*g1)*(1./det)

              df=-f1*f1-g1*g1
              f=-0.5*df
            
              delta=0.5*(abs(dx)+abs(dy))
!              delta=abs(dx)
 
              xo=wg
              yo=vsq
            
              xxn=wg+dx
              yyn=vsq+dy
!
              if(xxn .lt. 0.d0) then
!                irecov=2
                 xxn=abs(xxn)
              endif
!         
              if(xxn .gt. saftybig) then
                xxn=saftybig
              endif
!            
              if(yyn .lt. 0.d0) then
                yyn=abs(yyn)
              endif
!             
              if(abs(yyn) .gt. 1.d0) then
                dv=1.0d-2
                yyn=1.0-dv
              endif
 !            
            endif
          endif

        enddo
!
!- Newton-Raphson routine end -!
!
        if(irecov .eq. 0) then
          wg=xxn
          vsq=yyn
!
          gtmp=1.-vsq

          if(gtmp .lt. 0.d0 .or. gtmp .lt. safty) then
            gtmp=safty
            vsq=1.-gtmp 
          endif
           
          gf=1./sqrt(gtmp)
          ro=dro*(1./gf)
!
          w2=wg*(gtmp)
!          if(ieos .eq. 3) then !-polytropic EoS
!            pp=kpol*(ro**gam)
!          elseif(ieos .eq. 0) then !- gamma-law EoS 
            pp=gam2*(w2-ro)
!          endif

          if(ro .lt. 0.d0 .or. pp .lt. 0.d0) then  !- check pressure & density -!           
!          if(ro .lt. 0.d0) then
!            write(*,*) 'negative density', ro, 'at i, j, k', i, j ,k
!            ro=dmin*(x1(i)**(-3./2.))
!!            ro=dmin
!          endif
            irecov=1
          endif  
!          if (pp .lt. 0.d0) then
!!            irecov=3
!            write(*,*) 'negative pressure', pp, 'at i, j, k', i, j ,k
!            pp=(gam-1.)*pmin*(x1(i)**(-5./2.))
!!            pp=pmin
!          endif

          if(irecov .eq. 0) then !- velocity, pressure, & density all clear -! 
            qtcon(1)=qcon(1)+oncon(1)*qdotn
            qtcon(2)=qcon(2)+oncon(2)*qdotn
            qtcon(3)=qcon(3)+oncon(3)*qdotn

            util(1)=(gf*(1./(wg+bsq)))*(qtcon(1)+qdotb*bcon(1)*(1./wg))
            util(2)=(gf*(1./(wg+bsq)))*(qtcon(2)+qdotb*bcon(2)*(1./wg))
            util(3)=(gf*(1./(wg+bsq)))*(qtcon(3)+qdotb*bcon(3)*(1./wg))

!            vcon(1)=(1./(xxn+bsq))*(qtcon(1)+qdotb*bcon(1)/wg)
!            vcon(2)=(1./(xxn+bsq))*(qtcon(2)+qdotb*bcon(2)/wg)
!            vcon(3)=(1./(xxn+bsq))*(qtcon(3)+qdotb*bcon(3)/wg)
          endif
                 
          wg=xxn
          vsq=yyn   
!
        endif       
!
!=====================================================================@        
!- Another Newton Raphson Calculation -!

        if(irecov .eq. 1) then
           
!- set parameter  -!

          delta=1.d0
          deltaend = 1.d-6
          df=1.d0
          f=1.d0
          irecov=0

          xxn=wgo
          xo=wgo
!          xxn=vsqo
!          xo=vsqo
          
!- Newton-Raphson routine start -!
         
          do nnn=1,iter

            if(abs(delta) .ge. deltaend .and. irecov .eq. 0) then
           
              if(nnn .eq. iter) then
                irecov=1

              else
                wg=xxn

                wsq=wg*wg
                w3=wsq*wg
!
!- calculation of v^2 -!
!
                xsq=(bsq+wg)*(bsq+wg)
                vsq=abs((wsq*qtsq+qdotbsq*(bsq+2.*wg))*(1./(wsq*xsq)))
!       
                if(vsq .gt. 1.d0) then
                  write(6,*) 'vsq > 1.0 at', i, j, k 
                  vsq=1.d0-dv
                endif
              
                gtmp=1.-vsq

                if(gtmp .lt. safty) then
                  gtmp=safty
                  vsq=1.-gtmp 
                endif
!             
                gf=1./sqrt(gtmp)
                gfsq=gf*gf
!
                ro=dro*(1./gf)              
!
!- calculation of pressure -!
!              
                if(ieos .eq. 0) then !- gamma-law EoS -!
                  ptmp=gam2*(wg*gtmp-dro*sqrt(gtmp))
!                  ptmp=(gss/dro))*(ro**(gam))                 
                elseif(ieos .eq. 3) then !- polytropic EoS -!   
                  ptmp=kpol*((dro*sqrt(gtmp))**gam)
!                  ptmp=(gss/dro))*(ro**(gam)) 
                endif
!
!-calculation of dv^2/dw -!
!
                tmp1a=(bsq+wg)
                tmp1b=(bsq+wg)*(bsq+wg)*(bsq+wg)
                tmp1c=-2.*qtsq*(1./tmp1b)
                tmp1d=-2.*qdotbsq*(3.*wg*tmp1a+bsq*bsq)*(1./(tmp1b*w3))
                dvsqdw=tmp1c+tmp1d
!
!- calculation of dp/dw, dp/dv^2 -!
!              
                if(ieos .eq. 0) then !- gamma-law EoS -!
                  dpdw=gam2*gtmp
                  dpdvsq=-gam2*wg+0.5*gam2*dro*(1./(sqrt(gtmp)))
                  dpdw2=dpdw+dpdvsq*dvsqdw                
                endif
!
!- calculation of drho/dw -!
!              
                drodw=-0.5*dro*gf*dvsqdw
!             
!- calculation of f and df -!
!              
                f=dro*ptmp/(ro**(gam))-gss
                df=(dro/(ro**(gam)))*dpdw2-gam*ptmp*(ro**(gam+1.))*drodw
!
!- calculation of delta -!
!              
                delta=f/df
! 
                xo=wg
                xxn=wg-delta
!
                if(xxn .lt. 0.d0) then
!                 irecov=2
                  xxn=abs(xxn)
                endif
!         
                if(xxn .gt. 1.d6) then
                  xxn=1.d6
                endif
!              
              endif
            endif

          enddo
!
!- Newton-Raphson routine end -!
!
          if(irecov .eq. 0) then
            wg=xxn
            wsq=wg*wg
            xsq=(bsq+wg)*(bsq+wg)
            vsq=abs((wsq*qtsq+qdotbsq*(bsq+2.*wg))*(1./(wsq*xsq)))
!       
            if(vsq .gt. 1.d0) then
              write(6,*) 'vsq > 1.0 at', i, j, k 
              vsq=1.d0-dv
            endif           
!
            gtmp=1.-vsq

            if(gtmp .lt. safty) then
              gtmp=safty 
            endif
           
            gf=1./sqrt(gtmp)
            gfsq=gf*gf
          
            ro=dro/gf
          
            w2=wg*(gtmp)
            if(ieos .eq. 3) then !-polytropic EoS
               pp=kpol*(ro**gam)
!               pp=(gen*(1./gf))*(ro**(gam-1.))
            elseif(ieos .eq. 0) then !- gamma-law EoS 
               pp=gam2*(w2-ro)
!               pp=(gen*(1./gf))*(ro**(gam-1.))
            endif
            
            if(ro .lt. 0.d0) then
              write(*,*) 'negative density', ro, 'at', i, j ,k
              ro=dmin*(x1(i)**(-3./2.))
!              ro=dmin
            endif
            if (pp .lt. 0.d0) then
!              irecov=3
              write(*,*) 'negative pressure', pp, 'at', i, j ,k
              pp=pmin*(gam-1.)*(x1(i)**(-5./2.))
!              pp=pmin
            endif
          
            qtcon(1)=qcon(1)+oncon(1)*qdotn
            qtcon(2)=qcon(2)+oncon(2)*qdotn
            qtcon(3)=qcon(3)+oncon(3)*qdotn

            util(1)=(gf*(1./(wg+bsq)))*(qtcon(1)+qdotb*bcon(1)*(1./wg))
            util(2)=(gf*(1./(wg+bsq)))*(qtcon(2)+qdotb*bcon(2)*(1./wg))
            util(3)=(gf*(1./(wg+bsq)))*(qtcon(3)+qdotb*bcon(3)*(1./wg))
          endif
!
        endif          
!
!---------------------------------------------------------------------
!- Error message -!                
        if(irecov .eq. 1) then
          write(6,*) ' >> Not convergence in recov1-3' &
                    ,' at i, j, k:',i, j, k
          write(6,*) ' >> delta =', delta &
                    ,', at x1, x2, x3:',x1(i),x2(j),x3(k)
         
        elseif(irecov .eq. 2) then
          write(6,*) ' >> wg  < 0 in recov1-3' &
                    ,' at i, j, k:',i, j, k
          write(6,*) ' >> wg', wg &
                    ,', at x1, x2, x3:',x1(i),x2(j),x3(k)

!         elseif(irecov .eq. 3) then
!           write(6,*) ' >> ro or pr < 0 in recov1d' &
!                    ,' at i, j, k:',i, j, k
!           write(6,*) ' >> ro, pr', ro, pp &
!                     ,' at x1, x2, x3: ',x1(i),x2(j),x3(k)
        endif
!
!---------------------------------------------------------------------@
!   Set Primitive variables
!---------------------------------------------------------------------@
!
        if(irecov .eq. 0) then
         
          uri(1,i,j,k)=ro
          uri(2,i,j,k)=util(1)
          uri(3,i,j,k)=util(2)
          uri(4,i,j,k)=util(3)
!          uri(2,i,j,k)=vcon(1)
!          uri(3,i,j,k)=vcon(2)
!          uri(4,i,j,k)=vcon(3)
          uri(5,i,j,k)=pp
          uri(6,i,j,k)=ro*pp/(ro**(gam-1.))
          uri(7,i,j,k)=uu(7,i,j,k)*(1./detg1)
          uri(8,i,j,k)=uu(8,i,j,k)*(1./detg1)
          uri(9,i,j,k)=uu(9,i,j,k)*(1./detg1)
!        
!- cal of contravariant 4-velocity -!
          call calucon(util,ucon,gf,gcon1)
!- update conserved enthalpy -!
          uu(6,i,j,k)=detg1*uri(6,i,j,k)*ucon(0)
!          
       else
         
          uri(1,i,j,k)=ro1
          uri(2,i,j,k)=util1(1)
          uri(3,i,j,k)=util1(2)
          uri(4,i,j,k)=util1(3)
!          uri(2,i,j,k)=vcon1(1)
!          uri(3,i,j,k)=vcon1(2)
!          uri(4,i,j,k)=vcon1(3)
          uri(5,i,j,k)=pp1
          uri(6,i,j,k)=ro1*pp1/(ro1**(gam)) 
          uri(7,i,j,k)=uu(7,i,j,k)*(1./detg1)
          uri(8,i,j,k)=uu(8,i,j,k)*(1./detg1)
          uri(9,i,j,k)=uu(9,i,j,k)*(1./detg1)
!       
!- cal of contravariant 4-velocity -!
          call calucon(util,ucon,gf,gcon1)
!- update conserved enthalpy -!
          uu(6,i,j,k)=detg1*ro1*uri(6,i,j,k)*ucon(0)
!          
        endif
!
      enddo
    enddo
  enddo
!
!=====================================================================@
!
  deallocate( gcov1, gcon1, qcon, qcov, bcon, bcov, &
              oncov, oncon, util, util1, qtcon, ucon, stat=merr )
!
  return
end subroutine recov1b
!
!---------------------------------------------------------------------@
subroutine recov2(uu,uri,gcov,gcon,detg,x1,x2,x3,nm1,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!- Following Mignone & McKinney (2007) 1 variables inversion procedure -!
!- with variable EoS  

  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos, iter, dmin, pmin, &
                   kpol, akm
  implicit none
!
  integer :: i, j, k, m, n, nnn, nm1, is1, ie1, js1, je1, ks1, ke1
  integer :: merr, igue, irecov

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)

  real(8) :: x1(imax), x2(jmax), x3(kmax)

  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: qcon(:), qcov(:), bcon(:), bcov(:), &
                          oncov(:), oncon(:) 
  real(8), allocatable :: util(:), util1(:), qtcon(:)
      
  real(8) :: safty, saftybig, gam1, gam2, dv, alpha1, detg1 
  real(8) :: dro, ro, pp, ro1, pp1, roh, roe, bsq, qdotb, qdotbsq, qdotn, &
             qsq, qtsq, utsq, gfsq, gf, vsq, usq, wg, wsq, xsq, &
             delta, deltaend, df, f, xxn, xo, wgd, w3, gtmp, &
             xi, dpdxi, dpdro, dvdwd, dxidwd, drodwd, dpdwd
  real(8) :: tmp1a, tmp1b, tmp1c, tmp1d, wgp, wgm, vsq1, tmp3a, tmp3b, tmp3c, &
             tmp4a, tmp4b, tmp4c, tmp4d, tmp5a, tmp5b, tmp5c, tmp5d, &
             vsqo, rbh1
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), stat=merr)
  allocate(qcon(0:3), qcov(0:3), bcon(0:3), bcov(0:3), &
           oncov(0:3), oncon(0:3), util(1:3), util1(1:3), qtcon(1:3), stat=merr)
!
!--- define variable of safe number ----
      safty = 1.0d-5
      saftybig = 1.0d5
      dv=1.0d-4
      
      gam1=gam/(gam-1.0)
      gam2=1.0/gam1
      wg=0.d0

      igue=1
!
      rbh1=1.+sqrt(1.-akm*akm)
!      
!=====================================================================@
!
! ----  Calculation of variables ---
!
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
      do i=is1+nm1,ie1-nm1
!
!- set parameter for Newton Raphson -!
!
        delta=1.d0
        deltaend =1.d-8
        df=1.d0
        f=1.d0
!        irecov=0
!         
!- check the position within the BH radius -!
!
!        if(x1(i) .le. rbh1) then
!          irecov=2
!        else
          irecov=0
!        endif           
!         
!- copy of metric terms -!
!
         do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo
!
        alpha1=1./sqrt(-gcon1(0,0))
!
        detg1=detg(i,j,k)
!
!- rearrange of conserved variables -!       
!      
        dro=alpha1*uu(1,i,j,k)/detg1

!        qcov(0)=alpha1*(uu(5,i,j,k)+uu(1,i,j,k))*(1./detg1)
        qcov(0)=alpha1*(uu(5,i,j,k)-uu(1,i,j,k))*(1./detg1)
!        qcov(0)=alpha1*uu(5,i,j,k)*(1./detg1)
        qcov(1)=alpha1*uu(2,i,j,k)*(1./detg1)
        qcov(2)=alpha1*uu(3,i,j,k)*(1./detg1)
        qcov(3)=alpha1*uu(4,i,j,k)*(1./detg1)

        tmp4a=alpha1*uu(5,i,j,k)*(1./detg1)
        
        bcon(0)=0.d0
        bcon(1)=alpha1*uu(7,i,j,k)*(1./detg1)
        bcon(2)=alpha1*uu(8,i,j,k)*(1./detg1)
        bcon(3)=alpha1*uu(9,i,j,k)*(1./detg1)
!        
!- copy of primitive variables (previous time)         
!         
        ro=uri(1,i,j,k)
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
!        vcon(1)=uri(2,i,j,k)
!        vcon(2)=uri(3,i,j,k)
!        vcon(3)=uri(4,i,j,k)
        pp=uri(5,i,j,k)
        
        ro1=uri(1,i,j,k)
        util1(1)=uri(2,i,j,k)
        util1(2)=uri(3,i,j,k)
        util1(3)=uri(4,i,j,k)
!        vcon1(1)(i,j,k)=uri(2,i,j,k)
!        vcon1(2)(i,j,k)=uri(3,i,j,k)
!        vcon1(3)(i,j,k)=uri(4,i,j,k)
        pp1=uri(5,i,j,k)
!
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          roh=ro+(gam/(gam-1.0))*pp
        elseif(ieos .eq. 1) then
          roh=(5./2.)*pp &
              +sqrt((9./4.)*pp**2+ro**2)
        elseif(ieos .eq. 2) then
          roe=(3./2.)*(pp+((3.*pp**2) &
              /(2.0*ro+sqrt(2.*pp**2+4.*ro**2)) ))
          roh=ro+roe+pp
        endif

!- cal of covariant 3-magnetic field -! 
        call lower(bcon,bcov,gcov1)
!- cal of contravariant 4-q -!
        call upper(qcov,qcon,gcon1)
!- cal of B-field square (3-vector) -!
        bsq=bcon(1)*bcov(1)+bcon(2)*bcov(2)+bcon(3)*bcov(3)
!- cal of q*B -!
        qdotb=qcov(0)*bcon(0)+qcov(1)*bcon(1)+qcov(2)*bcon(2)+qcov(3)*bcon(3)
!- cal of (Q*B)^2 -!
        qdotbsq=qdotb*qdotb
!- set of covariant 4-vector n -!
        oncov(0)=-1.d0*alpha1
        oncov(1)=0.d0
        oncov(2)=0.d0
        oncov(3)=0.d0
!- cal of contravariant 4-vector n -!
        call upper(oncov,oncon,gcon1)
!- cal of Q*n -!
        qdotn=qcon(0)*oncov(0)
!
!- cal of square of Q -!
        qsq=0.d0
        qsq=qcov(0)*qcon(0)+qcov(1)*qcon(1) &
           +qcov(2)*qcon(2)+qcov(3)*qcon(3)
!- cal of Q^2+Q*n -!
        qtsq=qsq+(qdotn*qdotn)
!
!- Calculate W from last timestep and use for inital guess -!      
!
!- calculation of 4-velocity square 
        utsq=0.d0
!        vsq=0.d0
        do m=1,3
          do n=1,3
            utsq=utsq+gcov1(m,n)*util(m)*util(n)
!            vsq=vsq+gcov1(m,n)*util(m)*util(n)
          enddo
        enddo
!
        if(utsq .lt. 0.d0 .and. abs(utsq) .le. safty) then
          utsq=abs(utsq)
        endif
        if(utsq .lt. 0.d0 .or. abs(utsq) .gt. 1.d4) then
          utsq=1.d4
        endif 
!- calculation of Lorentz factor -!
        gfsq=1.+utsq
        gf=sqrt(gfsq)
!        gf=1./sqrt(1.-vsq)
!        gfsq=gf**2
        vsq=(gfsq-1.)*(1./gfsq)
!- calculation of W -!      
!        if(vsq .lt. 0.9d0 .and. bsq/roh .lt. 1.d0) then 
! 
!          ro=dro*(1./gf)
!          roh=ro+gam1*pp
!          wg=abs(gfsq*roh)
!
!        else
          
!          if(igue .eq. 1) then
!            vsq1=1.d0
!          else
            vsq1=vsq
!          endif
!            
          tmp1a=4.-vsq1
          tmp1b=4.*(bsq+qdotn)
          tmp1c=qtsq+bsq*bsq+2.*bsq*qdotn
          tmp1d=(tmp1b*tmp1b)-4.*tmp1a*tmp1c
          if(tmp1d .le. 0.d0) then
            tmp1d=abs(tmp1d)
          endif
          wgp=(-tmp1b+sqrt(tmp1d))*(1./(2.*tmp1a))
          wgm=(-tmp1b-sqrt(tmp1d))*(1./(2.*tmp1a))
!
         if(tmp1d .le. 0.d0) then
           write(*,*) "inside sqrt < 0 for cal wg", tmp1b, tmp1c
         endif
!
          if(wgp .gt. 0.d0) then
            wg=wgp
          elseif(wgm .ge. 0.d0) then
            wg=wgm
          endif
!         
!        endif
!
! Make sure that W is large enough so that v^2 < 1
!       
95      continue

        if(ieos .eq. 0) then !- gamma-law EoS -!
          pp=gam2*(wg*gf-dro*sqrt(gf))
        elseif(ieos .eq. 3) then !- polytropic EoS -!   
          pp=kpol*((dro*sqrt(gf))**gam)
        endif
       
        tmp3a=(wg*wg*wg)*(wg+2.*bsq)-qdotbsq*(2.*wg+bsq)
        tmp3b=(wg*wg)*(qtsq-bsq*bsq)
!
        tmp3c=(wg*wg*wg)+(wg*wg)*(0.5*bsq+qdotn-pp)-0.5*qdotbsq
!        
!        do nnn=1,iter
          if(tmp3a .le. tmp3b) then
!             write(*,*) "added wg at i,j,k=",i,j,k
             wg= abs(1.5*wg)
             go to 95
!             tmp3a=(wg*wg*wg)*(wg+2.*bsq)-qdotbsq*(2.*wg+bsq)
!             tmp3b=(wg*wg)*(qtsq-bsq*bsq)
          endif

!        enddo
 
        if(tmp3c .gt. 0.d0) then
          wg=abs(0.9*wg)
          go to 95 
        endif 
!
!=====================================================================@
!- Newton Rapson calculation -!

        wgd=wg-dro
        xxn=wgd
        xo=wgd
!        xxn=wg
!        xo=wg
        
!- Newton-Raphson routine start -!

        do nnn=1,iter

          if(abs(delta) .ge. deltaend .and. irecov .eq. 0) then
           
            if(nnn .eq. iter) then
              irecov=1

           else
              vsqo=vsq
              wgd=xxn
              wg=wgd+dro
!              wg=xxn
!              
!- Calculation of velocity and Lorentz factor -!
              wsq=wg*wg
              w3=wsq*wg
              xsq=(bsq+wg)*(bsq+wg)
              vsq=abs((wsq*qtsq+qdotbsq*(bsq+2.*wg))*(1./(wsq*xsq)))
!       
              if(abs(vsq) .gt. 1.d0) then
                write(6,*) 'vsq > 1.0 at', i, j, k 
                vsq=1.d0-dv
!                vsq=vsqo
              endif
              if(vsq .lt. 0.d0) then
                vsq=abs(vsq)
              endif   
!         
              gtmp=1.-vsq
              if(gtmp .lt. 0.d0 .or. gtmp .lt. safty) then
                 gtmp=safty
                 vsq=1.-gtmp 
              endif
!
!              gf=1./sqrt(gtmp)
              usq=(1.-gtmp)/(gtmp)
              gf=sqrt(1.+usq)
              gfsq=gf*gf
!- Calculation of density, xi -!
              ro=dro/gf
!
!              xi=(wgd*(1./gfsq))-(dro*usq)*(1./((1.0+gf)*gfsq))
              xi=(wg-dro*gf)*(1./gfsq)
!- Calculation of pressure from EoS -!              
              if(ieos .eq. 0) then
                pp=gam2*xi
              elseif(ieos .eq. 1 .or. ieos .eq. 2) then
                tmp4a=2.*xi*(xi+2.*ro)
                tmp4b=5.*(xi+ro)+sqrt(9.*(xi**2)+18.*ro*xi+25.*(ro**2))
                pp=tmp4a*(1./tmp4b)
              endif
!
!- calculation of dp/dxi and dp/dro -!              
              if(ieos .eq. 0) then !- gamma-law EoS -!
                dpdxi=gam2
                dpdro=0.d0
              elseif(ieos .eq. 1 .or. ieos .eq. 2) then !- TM EoS -!
                dpdxi=(2.*xi+2.*ro-5.*pp)*(1./(5.*ro+5.*xi-8.*pp))
                dpdro=(2.*xi-5.*pp)*(1./(5.*ro+5.*xi-8.*pp))
              endif
!- calculation of dv^2/dw -!
              tmp1a=(bsq+wg)
              tmp1b=(bsq+wg)*(bsq+wg)*(bsq+wg)
              if(qtsq .eq. 0.d0) then
                tmp1c=0.d0
              else   
                tmp1c=-2.*qtsq*(1./tmp1b)
              endif   
              if(qdotbsq .eq. 0.d0) then
                tmp1d=0.d0
              else   
                tmp1d=-2.*qdotbsq*(3.*wg*tmp1a+bsq*bsq)*(1./(tmp1b*w3))
              endif   
              
              dvdwd=tmp1c+tmp1d
!- calculation of dxi/dw and dro/dw -!          
              dxidwd=(1./gfsq)-0.5*gf*(dro+2.*gf*xi)*dvdwd
              drodwd=-0.5*dro*gf*dvdwd
!              drodwd=0.5*dro*gf*dvdwd
!- calculation of dp/dw -!              
              dpdwd=dpdxi*dxidwd+dpdro*drodwd
          
! --- calculation of f(W') and df(W')/dW' ----
!
              tmp5a=bsq*qtsq-qdotbsq
              tmp5b=bsq+wgd+dro
!              tmp5b=bsq+wg
!            
              if(tmp5a .eq. 0.d0) then
                tmp5c=0.d0
                tmp5d=0.d0
              else
                tmp5c=0.5*tmp5a*(1./(tmp5b*tmp5b))
                tmp5d=tmp5a*(1./(tmp5b*tmp5b*tmp5b))
              endif
!
!- Mignone & McKinney
              f=-(qdotn+dro)-wgd+pp-0.5*bsq-tmp5c
!              f=-qdotn-wg+pp-0.5*bsq-tmp5c
              df=-1.+dpdwd+tmp5d
!
!- Noble et al (2006)
!              f=-wg-0.5*bsq*(1.+vsq)+(0.5*qdotbsq/wsq)-qdotn+pp
!              df=-1.+dpdwd-(qdotbsq/w3)-0.5*bsq*dvdwd
!             
              delta=f/df                            
!             
              xo=wgd
!              xo=wg
              xxn=wgd-delta
!              xxn=wg-delta
!  
              if(xxn .lt. 0.d0) then
!                irecov=2
                 xxn=abs(xxn)
              endif
!        
              if(xxn .gt. 1.d8) then
                xxn=1.d8
              endif             
!
            endif
          endif

        enddo
!
!- Newton-Raphson routine end -!
!
        if(irecov .eq. 0) then
          wgd=xxn
          wg=wgd+dro
!          wg=xxn
          wsq=wg*wg
          xsq=(bsq+wg)*(bsq+wg)
          vsq=abs((wsq*qtsq+qdotbsq*(bsq+2.*wg))*(1./(wsq*xsq)))
!       
          if(vsq .lt. 0.d0) then
            vsq=abs(vsq)
          endif   
          
          if(vsq .gt. 1.d0) then
            write(6,*) 'vsq > 1.0 at', i, j, k 
            vsq=1.d0-dv
!            vsq=vsqo
          endif
!         
          gtmp=1.-vsq

          if(gtmp .lt. safty) then
            gtmp=safty 
          endif
         
!          gf=1./sqrt(gtmp)
          usq=(1.-gtmp)/gtmp
          gf=sqrt(1.+usq)
          gfsq=gf*gf
!- Calculation of density, xi -!
          ro=dro/gf
          xi=(wgd*(1./gfsq))-(dro*usq)*(1./((1.+gf)*gfsq))
!          xi=(wg-dro*gf)*(1./gfsq)
          xi=wg*(1./gfsq)-ro
!- Calculation of pressure from EoS -!              
          if(ieos .eq. 0) then
            pp=gam2*xi
          elseif(ieos .eq. 1 .or. ieos .eq. 2) then
            tmp4a=2.*xi*(xi+2.*ro)
            tmp4b=5.*(xi+ro)+sqrt(9.*(xi*xi)+18.*ro*xi+25.*(ro*ro))
            pp=tmp4a*(1./tmp4b)
          endif

          if(ro .lt. 0.d0) then
            write(*,*) 'negative density', ro, 'at i, j, k', i, j ,k
            ro=dmin*(x1(i)**(-3./2.))
!            ro=dmin
          endif
          if (pp .lt. 0.d0) then
!            irecov=3
            write(*,*) 'negative pressure', pp, 'at i, j, k', i, j ,k
            write(*,*) 'dro, gf, wg, wgd, xi=', dro, gf, wg, wgd, xi 
            pp=pmin*(gam-1.)*(x1(i)**(-5./2.))
!            pp=pmin
          endif
!        
          if(irecov .ne. 3) then
            qtcon(1)=qcon(1)+oncon(1)*qdotn
            qtcon(2)=qcon(2)+oncon(2)*qdotn
            qtcon(3)=qcon(3)+oncon(3)*qdotn

            util(1)=(gf*(1./(wg+bsq)))*(qtcon(1)+qdotb*bcon(1)*(1./wg))
            util(2)=(gf*(1./(wg+bsq)))*(qtcon(2)+qdotb*bcon(2)*(1./wg))
            util(3)=(gf*(1./(wg+bsq)))*(qtcon(3)+qdotb*bcon(3)*(1./wg))

!            vcon(1)=(1./(xxn+bsq))*(qtcon(1)+qdotb*bcon(1)/wg)
!            vcon(2)=(1./(xxn+bsq))*(qtcon(2)+qdotb*bcon(2)/wg)
!            vcon(3)=(1./(xxn+bsq))*(qtcon(3)+qdotb*bcon(3)/wg)
           endif
        endif   
!
!---------------------------------------------------------------------
!- Error message -!                
        if(irecov .eq. 1) then  
          write(6,*) ' >> Not convergence in recov2' &
                    ,' at i, j, k:',i, j, k
          write(6,*) ' >> delta =', delta &
                    ,', at x1, x2, x3:',x1(i),x2(j),x3(k)
         
        elseif(irecov .eq. 2) then
          write(6,*) ' >> wg  < 0 in recov2' &
                    ,' at i, j, k:',i, j, k
          write(6,*) ' >> wg', wg &
                    ,', at x1, x2, x3:',x1(i),x2(j),x3(k)

!         elseif(irecov .eq. 3) then
!           write(6,*) ' >> ro or pr < 0 in recov2d' &
!                     ,' at i, j, k:',i, j, k
!           write(6,*) ' >> ro, pr', ro, pp &
!                     ,' at x1, x2, x3: ',x1(i),x2(j),x3(k)
        endif
!
!---------------------------------------------------------------------@
!   Set Primitive variables
!---------------------------------------------------------------------@
!
        if(irecov .eq. 0) then
         
          uri(1,i,j,k)=ro
          uri(2,i,j,k)=util(1)
          uri(3,i,j,k)=util(2)
          uri(4,i,j,k)=util(3)
!          uri(2,i,j,k)=vcon(1)
!          uri(3,i,j,k)=vcon(2)
!          uri(4,i,j,k)=vcon(3)
          uri(5,i,j,k)=pp
          uri(6,i,j,k)=pp/(ro**(gam-1.))
          uri(7,i,j,k)=uu(7,i,j,k)*(1./detg1)
          uri(8,i,j,k)=uu(8,i,j,k)*(1./detg1)
          uri(9,i,j,k)=uu(9,i,j,k)*(1./detg1)
        
        else
         
          uri(1,i,j,k)=ro1
          uri(2,i,j,k)=util1(1)
          uri(3,i,j,k)=util1(2)
          uri(4,i,j,k)=util1(3)
!          uri(2,i,j,k)=vcon1(1)
!          uri(3,i,j,k)=vcon1(2)
!          uri(4,i,j,k)=vcon1(3)
          uri(5,i,j,k)=pp1
          uri(6,i,j,k)=pp1/(ro1**(gam-1.))
          uri(7,i,j,k)=uu(7,i,j,k)*(1./detg1)
          uri(8,i,j,k)=uu(8,i,j,k)*(1./detg1)
          uri(9,i,j,k)=uu(9,i,j,k)*(1./detg1)
         
        endif
!
      enddo
    enddo
  enddo
!
!=====================================================================@
!
  deallocate( gcov1, gcon1, qcon, qcov, bcon, bcov, &
              oncov, oncon, util, util1, qtcon, stat=merr )
!
  return
end subroutine recov2
!
!---------------------------------------------------------------------@
subroutine recov3(uu,uri,gcov,gcon,detg,x1,x2,x3,nm1,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!- Noble et al. (2006) 1 variable inversion procedure -!
!- (ideal EoS only)
  
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos, iter, dmin, pmin, &
                   kpol, akm, pi
  implicit none
!
  integer :: i, j, k, m, n, nnn, nm1, is1, ie1, js1, je1, ks1, ke1
  integer :: merr, igue, irecov

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)

  real(8) :: x1(imax), x2(jmax), x3(kmax)

  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: qcon(:), qcov(:), bcon(:), bcov(:), &
                          oncov(:), oncon(:) 
  real(8), allocatable :: util(:), util1(:), qtcon(:)
      
  real(8) :: safty, saftybig, gam1, gam2, dv, alpha1, detg1 
  real(8) :: dro, ro, pp, ro1, pp1, roh, roe, bsq, qdotb, qdotbsq, qdotn, &
             qsq, qtsq, utsq, gfsq, gf, vsq, wg, wsq, xsq, qsq1, &
             delta, deltaend, df, f, xxn, xo, w2, w3, gtmp, ptmp, &
             dpdw, dpdw2, dpdvsq, dvsqdw 
  real(8) :: tmp1a, tmp1b, tmp1c, tmp1d, wgp, wgm, vsq1, tmp3a, tmp3b, &
             tmp3c, tmp4a
  real(8) :: rin, rbh1, rr, th
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), stat=merr)
  allocate(qcon(0:3), qcov(0:3), bcon(0:3), bcov(0:3), &
           oncov(0:3), oncon(0:3), util(1:3), util1(1:3), qtcon(1:3), stat=merr)
!
!--- define variable of safe number ----
      safty = 1.0d-6
      saftybig = 1.0d6
      dv=1.0d-2
!      
!--- Parameters ---!
!      
      gam1=gam/(gam-1.0)
      gam2=1.0/gam1
      wg=0.d0

      igue=1
!
      rin=15.d0
      rbh1=1.+sqrt(1.-akm*akm)
!
!
!=====================================================================@
!
! ----  Calculation of variables ---
!
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
       do i=is1+nm1,ie1-nm1
!
!- set parameter for Newton Raphson -!
!
      delta=1.d0
      deltaend = 1.d-6
      df=1.d0
      f=1.d0
!      irecov=0
! 
!- check the position within the BH radius -!
!
        if(x1(i) .le. rbh1) then
          irecov=1
        else
          irecov=0
        endif   
!          
!- copy of metric terms -!
!
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo
!
        alpha1=1./sqrt(-gcon1(0,0))
!
        detg1=detg(i,j,k)
!
!- rearrange of conserved variables -!       
!      
        dro=alpha1*uu(1,i,j,k)*(1./detg1)

!        qcov(0)=alpha1*(uu(5,i,j,k)+uu(1,i,j,k))*(1./detg1)
        qcov(0)=alpha1*(uu(5,i,j,k)-uu(1,i,j,k))*(1./detg1)
!        qcov(0)=alpha1*uu(5,i,j,k)*(1./detg1)
        qcov(1)=alpha1*uu(2,i,j,k)*(1./detg1)
        qcov(2)=alpha1*uu(3,i,j,k)*(1./detg1)
        qcov(3)=alpha1*uu(4,i,j,k)*(1./detg1)

        tmp4a=alpha1*uu(5,i,j,k)*(1./detg1)
        
        bcon(0)=0.d0
        bcon(1)=alpha1*uu(7,i,j,k)*(1./detg1)
        bcon(2)=alpha1*uu(8,i,j,k)*(1./detg1)
        bcon(3)=alpha1*uu(9,i,j,k)*(1./detg1)
!        
!- copy of primitive variables (previous time)         
!         
        ro=uri(1,i,j,k)
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
!        vcon(1)=uri(2,i,j,k)
!        vcon(2)=uri(3,i,j,k)
!        vcon(3)=uri(4,i,j,k)
        pp=uri(5,i,j,k)
!
!- store for previous time        
!        
        ro1=uri(1,i,j,k)
        util1(1)=uri(2,i,j,k)
        util1(2)=uri(3,i,j,k)
        util1(3)=uri(4,i,j,k)
!        vcon1(1)(i,j,k)=uri(2,i,j,k)
!        vcon1(2)(i,j,k)=uri(3,i,j,k)
!        vcon1(3)(i,j,k)=uri(4,i,j,k)
        pp1=uri(5,i,j,k)
!
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          roh=ro+(gam/(gam-1.0))*pp
        elseif(ieos .eq. 1) then
          roh=(5./2.)*pp &
              +sqrt((9./4.)*pp**2+ro**2)
        elseif(ieos .eq. 2) then
          roe=(3./2.)*(pp+((3.*pp**2) &
              /(2.0*ro+sqrt(2.*pp**2+4.*ro**2)) ))
          roh=ro+roe+pp
        endif
!
!- check B-field minimum -!
!
!        if(abs(bcon(1)) .lt. bmin) then
!          bcon(1)=0.d0
!        endif
!        if(abs(bcon(2)) .lt. bmin) then
!          bcon(2)=0.d0
!        endif       
!        if(abs(bcon(3)) .lt. bmin) then
!          bcon(3)=0.d0
!        endif   
        
!- cal of covariant 3-magnetic field -! 
        call lower(bcon,bcov,gcov1)
!- cal of contravariant 4-q -!
        call upper(qcov,qcon,gcon1)
!- cal of B-field square (3-vector) -!
        bsq=bcon(1)*bcov(1)+bcon(2)*bcov(2)+bcon(3)*bcov(3)
!- cal of q*B -!
        qdotb=qcov(0)*bcon(0)+qcov(1)*bcon(1)+qcov(2)*bcon(2)+qcov(3)*bcon(3)
!- cal of (Q*B)^2 -!
        qdotbsq=qdotb*qdotb
!- set of covariant 4-vector n -!
        oncov(0)=-1.d0*alpha1
        oncov(1)=0.d0
        oncov(2)=0.d0
        oncov(3)=0.d0
!- cal of contravariant 4-vector n -!
        call upper(oncov,oncon,gcon1)
!- cal of Q*n -!
        qdotn=qcon(0)*oncov(0)
!
!- cal of square of Q -!
        qsq=0.d0
        qsq=qcov(0)*qcon(0)+qcov(1)*qcon(1) &
           +qcov(2)*qcon(2)+qcov(3)*qcon(3)
!- cal of Q^2+Q*n -!
        qtsq=qsq+(qdotn*qdotn)
!
!        qsq1=0.d0
!        qsq1=qcov(1)*qcon(1)+qcov(2)*qcon(2)+qcov(3)*qcon(3)
!
!
!- Calculate W from last timestep and use for inital guess -!      
!
!- calculation of 4-velocity square 
        utsq=0.d0
!        vsq=0.d0
        do m=1,3
          do n=1,3
            utsq=utsq+gcov1(m,n)*util(m)*util(n)
!            vsq=vsq+gcov1(m,n)*util(m)*util(n)
          enddo
        enddo
!
        if(utsq .lt. 0.d0 .and. abs(utsq) .le. safty) then
          utsq=abs(utsq)
        endif
        if(utsq .lt. 0.d0 .or. abs(utsq) .gt. saftybig) then
          utsq=saftybig
        endif 
!- calculation of Lorentz factor -!
        gfsq=1.+utsq
        gf=sqrt(gfsq)
!        gf=1./sqrt(1.-vsq)
!        gfsq=gf**2
        vsq=(gfsq-1.)*(1./gfsq)
!- calculation of W -!      
!        if(vsq .lt. 0.9d0 .and. bsq/roh .lt. 1.d0) then 
! 
!          ro=dro/gf
!          roh=ro+gam1*pp
!          wg=abs(gfsq*roh)
!
!        else
          
!          if(igue .eq. 1) then
!            vsq1=1.d0
!          else
            vsq1=vsq
!          endif
!            
          tmp1a=4.-vsq1
          tmp1b=4.*(bsq+qdotn)
          tmp1c=qtsq+bsq*bsq+2.*bsq*qdotn
          tmp1d=(tmp1b*tmp1b)-4.*tmp1a*tmp1c
          if(tmp1d .le. 0.d0) then
            tmp1d=abs(tmp1d)
          endif
          wgp=(-tmp1b+sqrt(tmp1d))*(1./(2.*tmp1a))
          wgm=(-tmp1b-sqrt(tmp1d))*(1./(2.*tmp1a))
!
         if(tmp1d .le. 0.d0) then
           write(*,*) "inside sqrt < 0 for cal wg", tmp1b, tmp1c
          endif
!
          if(wgp .gt. 0.d0) then
            wg=wgp
          elseif(wgm .ge. 0.d0) then
            wg=wgm
          endif
          
!        endif
!
! Make sure that W is large enough so that v^2 < 1
!       
  95    continue

        if(ieos .eq. 0) then !- gamma-law EoS -!
          pp=gam2*(wg*gf-dro*sqrt(gf))
        elseif(ieos .eq. 3) then !- polytropic EoS -!   
          pp=kpol*((dro*sqrt(gf))**gam)
        endif
!       
        tmp3a=(wg*wg*wg)*(wg+2.*bsq)-qdotbsq*(2.*wg+bsq)
        tmp3b=(wg*wg)*(qtsq-bsq*bsq)
!        
        tmp3c=(wg*wg*wg)+(wg*wg)*(0.5*bsq+qdotn-pp)-0.5*qdotbsq
!        
!        do nnn=1,iter
          if(tmp3a .le. tmp3b) then
!             write(*,*) "added wg at i,j,k=",i,j,k
             wg= abs(1.5*wg)
             go to 95
!             tmp3a=(wg*wg*wg)*(wg+2.*bsq)-qdotbsq*(2.*wg+bsq)
!             tmp3b=(wg*wg)*(qtsq-bsq*bsq)
          endif
!        enddo
!          
         if(tmp3c .gt. 0.d0) then
           wg=abs(0.9*wg)
           go to 95 
         endif 
!
!=====================================================================@
!- Newton Rapson calculation -!

        xxn=wg
        xo=wg

!- Newton-Raphson routine start -!
         
        do nnn=1,iter

          if(abs(delta) .ge. deltaend .and. irecov .eq. 0) then
           
            if(nnn .eq. iter) then
              irecov=1

            else
              wg=xxn

              wsq=wg*wg
              w3=wsq*wg
!
!- calculation of v^2 -!
!
              xsq=(bsq+wg)*(bsq+wg)
              vsq=abs((wsq*qtsq+qdotbsq*(bsq+2.*wg))*(1./(wsq*xsq)))
!       
              if(vsq .gt. 1.d0) then
                write(6,*) 'vsq > 1.0 at', i, j, k 
                vsq=1.d0-dv
              endif
              
              gtmp=1.-vsq

              if(gtmp .lt. safty) then
                gtmp=safty
                vsq=1.-gtmp 
              endif  
!
!-calculation of dv^2/dw -!
!
              tmp1a=(bsq+wg)
              tmp1b=(bsq+wg)*(bsq+wg)*(bsq+wg)
              tmp1c=-2.*qtsq*(1./tmp1b)
              tmp1d=-2.*qdotbsq*(3.*wg*tmp1a+bsq*bsq)*(1./(tmp1b*w3))
              dvsqdw=tmp1c+tmp1d
!
!- calculation of pressure -!
!              
              if(ieos .eq. 0) then !- gamma-law EoS -!
                ptmp=gam2*(wg*gtmp-dro*sqrt(gtmp))
              elseif(ieos .eq. 3) then !- polytropic EoS -!   
                ptmp=kpol*((dro*sqrt(gtmp))**gam)
              endif
!
!- calculation of dp/dw, dp/dv^2 -!
!              
              if(ieos .eq. 0) then !- gamma-law EoS -!
                dpdw=gam2*gtmp
                dpdvsq=-gam2*wg+0.5*gam2*dro*(1./(sqrt(gtmp)))
                dpdw2=dpdw+dpdvsq*dvsqdw                
              endif
!             
!- calculation of f and df -!
!              
              f=-wg-0.5*bsq*(1.+vsq)+(0.5*qdotbsq*(1./wsq))-qdotn+ptmp
              df=-1.+dpdw2-(qdotbsq*(1./w3))-0.5*bsq*dvsqdw
!
!- calculation of delta -!
!              
              delta=f*(1./df)
! 
              xo=wg
              xxn=wg-delta
!
              if(xxn .lt. 0.d0) then
!                irecov=2
                 xxn=abs(xxn)
              endif
!         
              if(xxn .gt. 1.d6) then
                xxn=1.d6
              endif
!              
            endif
          endif

        enddo
!
!- Newton-Raphson routine end -!
!
        if(irecov .eq. 0) then
          wg=xxn
          wsq=wg*wg
          xsq=(bsq+wg)*(bsq+wg)
          vsq=abs((wsq*qtsq+qdotbsq*(bsq+2.*wg))*(1./(wsq*xsq)))
!       
          if(vsq .gt. 1.d0) then
            write(6,*) 'vsq > 1.0 at', i, j, k 
            vsq=1.d0-dv
          endif           
!
          gtmp=1.-vsq

          if(gtmp .lt. safty) then
            gtmp=safty 
          endif
           
          gf=1./sqrt(gtmp)
          ro=dro/gf
!
          w2=wg*(gtmp)
          if(ieos .eq. 3) then !-polytropic EoS
            pp=kpol*(ro**gam)
          elseif(ieos .eq. 0) then !- gamma-law EoS 
            pp=gam2*(w2-ro)
          endif
            
          if(ro .lt. 0.d0) then
            write(*,*) 'negative density', ro, 'at', i, j ,k
            ro=dmin*(x1(i)**(-3./2.))
!            ro=dmin
          endif
          if (pp .lt. 0.d0) then
!            irecov=3
            write(*,*) 'negative pressure', pp, 'at', i, j ,k
            pp=pmin*(gam-1.)*(x1(i)**(-5./2.))
!            pp=pmin
          endif

!          if(irecov .ne. 3) then
            qtcon(1)=qcon(1)+oncon(1)*qdotn
            qtcon(2)=qcon(2)+oncon(2)*qdotn
            qtcon(3)=qcon(3)+oncon(3)*qdotn

            util(1)=(gf*(1./(wg+bsq)))*(qtcon(1)+qdotb*bcon(1)*(1./wg))
            util(2)=(gf*(1./(wg+bsq)))*(qtcon(2)+qdotb*bcon(2)*(1./wg))
            util(3)=(gf*(1./(wg+bsq)))*(qtcon(3)+qdotb*bcon(3)*(1./wg))

!            vcon(1)=(1./(xxn+bsq))*(qtcon(1)+qdotb*bcon(1)/wg)
!            vcon(2)=(1./(xxn+bsq))*(qtcon(2)+qdotb*bcon(2)/wg)
!            vcon(3)=(1./(xxn+bsq))*(qtcon(3)+qdotb*bcon(3)/wg)
!           endif
        endif
                 
        wg=xxn  

!---------------------------------------------------------------------
!- Error message -!                
        if(irecov .eq. 1) then
          write(6,*) ' >> Not convergence in recov3' &
                    ,' at:',i, j, k
          write(6,*) ' >> delta =', delta &
                    ,', at x1, x2, x3:',x1(i),x2(j),x3(k)
         
        elseif(irecov .eq. 2) then
          write(6,*) ' >> wg  < 0 in recov3' &
                    ,' at:',i, j, k
          write(6,*) ' >> wg', wg &
                    ,', at x1, x2, x3:',x1(i),x2(j),x3(k)

!         elseif(irecov .eq. 3) then
!           write(6,*) ' >> ro or pr < 0 in recov1d' &
!                    ,' at i, j, k:',i, j, k
!           write(6,*) ' >> ro, pr', ro, pp &
!                     ,' at x1, x2, x3: ',x1(i),x2(j),x3(k)
        endif
!
!---------------------------------------------------------------------@
!   Set Primitive variables
!---------------------------------------------------------------------@
!
        if(irecov .eq. 0) then
         
          uri(1,i,j,k)=ro
          uri(2,i,j,k)=util(1)
          uri(3,i,j,k)=util(2)
          uri(4,i,j,k)=util(3)
!          uri(2,i,j,k)=vcon(1)
!          uri(3,i,j,k)=vcon(2)
!          uri(4,i,j,k)=vcon(3)
          uri(5,i,j,k)=pp
          uri(6,i,j,k)=pp/(ro**(gam-1.))
!          uri(7,i,j,k)=bcon(1)
!          uri(8,i,j,k)=bcon(2)
!          uri(9,i,j,k)=bcon(3)
!          
          uri(7,i,j,k)=uu(7,i,j,k)*(1./detg1)
          uri(8,i,j,k)=uu(8,i,j,k)*(1./detg1)
          uri(9,i,j,k)=uu(9,i,j,k)*(1./detg1)
!
        rr=x1(i)
        th=x3(k)*180.d0/pi 
!        if(x1(i) .ge. 35.d0 .and. x1(i) .le. 35.5d0) then
!          if((uri(7,i,j,k)) .gt. 0.d0 .or. abs(uri(9,i,j,k)) .gt. 0.d0) then
!            write(*,*) 'rr, th, br, bph, bth in recov3a=', rr, th, uri(7,i,j,k), uri(8,i,j,k), uri(9,i,j,k)
!          endif           
!        endif
       
        else
         
          uri(1,i,j,k)=ro1
          uri(2,i,j,k)=util1(1)
          uri(3,i,j,k)=util1(2)
          uri(4,i,j,k)=util1(3)
!          uri(2,i,j,k)=vcon1(1)
!          uri(3,i,j,k)=vcon1(2)
!          uri(4,i,j,k)=vcon1(3)
          uri(5,i,j,k)=pp1
          uri(6,i,j,k)=pp1/(ro1**(gam-1.))
          uri(7,i,j,k)=uu(7,i,j,k)*(1./detg1)
          uri(8,i,j,k)=uu(8,i,j,k)*(1./detg1)
          uri(9,i,j,k)=uu(9,i,j,k)*(1./detg1)
         
        endif
!
      enddo
    enddo
  enddo
!
!=====================================================================@
!
  deallocate( gcov1, gcon1, qcon, qcov, bcon, bcov, &
              oncov, oncon, util, util1, qtcon, stat=merr )
!
  return
end subroutine recov3
!
!---------------------------------------------------------------------@
subroutine recov4(uu,uri,gcov,gcon,detg,x1,x2,x3,nm1,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!- 1 variable inversion procedure -!
!- (Polytropic EoS only)
  
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos, iter, dmin, pmin, &
                   kpol, akm
  implicit none
!
  integer :: i, j, k, m, n, nnn, nm1, is1, ie1, js1, je1, ks1, ke1
  integer :: merr, igue, irecov

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)

  real(8) :: x1(imax), x2(jmax), x3(kmax)

  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: qcon(:), qcov(:), bcon(:), bcov(:), &
                          oncov(:), oncon(:) 
  real(8), allocatable :: util(:), util1(:), qtcon(:)
      
  real(8) :: safty, saftybig, gam1, gam2, dv, alpha1, detg1 
  real(8) :: dro, ro, pp, ro1, pp1, roh, roe, bsq, qdotb, qdotbsq, qdotn, &
             qsq, qtsq, utsq, gfsq, gf, vsq, wg, wsq, xsq, qsq1, &
             delta, deltaend, df, f, xxn, xo, w3, gtmp, &
             dgfdvsq, dwdvsq 
  real(8) :: tmp1a, tmp1b, tmp1c, tmp1d, wgp, wgm, vsq1, tmp3a, tmp3b, &
             tmp3c, tmp4a
  real(8) :: rin, rbh1
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), stat=merr)
  allocate(qcon(0:3), qcov(0:3), bcon(0:3), bcov(0:3), &
           oncov(0:3), oncon(0:3), util(1:3), util1(1:3), qtcon(1:3), stat=merr)
!
!--- define variable of safe number ----
      safty = 1.0d-6
      saftybig = 1.0d6
      dv=1.0d-5
      
      gam1=gam/(gam-1.0)
      gam2=1.0/gam1
      wg=0.d0

      igue=1
!
      rin=15.d0
      rbh1=1.+sqrt(1.-akm*akm)      

!=====================================================================@
!
! ----  Calculation of variables ---
!
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
      do i=is1+nm1,ie1-nm1
!
!- set parameter for Newton Raphson  -!
!
        delta=1.d0
        deltaend = 1.d-8
        df=1.d0
        f=1.d0
!        irecov=0
! 
!- check the position within the BH radius -!
!
!        if(x1(i) .le. rbh1) then
!          irecov=1
!        else
          irecov=0
!        endif
!
!- copy of metric terms -!
!
         do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo
!
        alpha1=1./sqrt(-gcon1(0,0))
!
        detg1=detg(i,j,k)
      
        dro=alpha1*uu(1,i,j,k)*(1./detg1)

!        qcov(0)=alpha1*(uu(5,i,j,k)+uu(1,i,j,k))*(1./detg1)
        qcov(0)=alpha1*(uu(5,i,j,k)-uu(1,i,j,k))*(1./detg1)
!        qcov(0)=alpha1*uu(5,i,j,k)*(1./detg1)
        qcov(1)=alpha1*uu(2,i,j,k)*(1./detg1)
        qcov(2)=alpha1*uu(3,i,j,k)*(1./detg1)
        qcov(3)=alpha1*uu(4,i,j,k)*(1./detg1)

        tmp4a=alpha1*uu(5,i,j,k)*(1./detg1)
        
        bcon(0)=0.d0
        bcon(1)=alpha1*uu(7,i,j,k)*(1./detg1)
        bcon(2)=alpha1*uu(8,i,j,k)*(1./detg1)
        bcon(3)=alpha1*uu(9,i,j,k)*(1./detg1)
         
        ro=uri(1,i,j,k)
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
!        vcon(1)=uri(2,i,j,k)
!        vcon(2)=uri(3,i,j,k)
!        vcon(3)=uri(4,i,j,k)
        pp=uri(5,i,j,k)
        
        ro1=uri(1,i,j,k)
        util1(1)=uri(2,i,j,k)
        util1(2)=uri(3,i,j,k)
        util1(3)=uri(4,i,j,k)
!        vcon1(1)(i,j,k)=uri(2,i,j,k)
!        vcon1(2)(i,j,k)=uri(3,i,j,k)
!        vcon1(3)(i,j,k)=uri(4,i,j,k)
        pp1=uri(5,i,j,k)
!
        roh=ro+gam1*pp

!- cal of covariant 3-magnetic field -! 
        call lower(bcon,bcov,gcov1)
!- cal of contravariant 4-q -!
        call upper(qcov,qcon,gcon1)
!- cal of B-field square (3-vector) -!
        bsq=bcon(1)*bcov(1)+bcon(2)*bcov(2)+bcon(3)*bcov(3)
!- cal of q*B -!
        qdotb=qcov(0)*bcon(0)+qcov(1)*bcon(1)+qcov(2)*bcon(2)+qcov(3)*bcon(3)
!- cal of (Q*B)^2 -!
        qdotbsq=qdotb*qdotb
!- set of covariant 4-vector n -!
        oncov(0)=-1.d0*alpha1
        oncov(1)=0.d0
        oncov(2)=0.d0
        oncov(3)=0.d0
!- cal of contravariant 4-vector n -!
        call upper(oncov,oncon,gcon1)
!- cal of Q*n -!
        qdotn=qcon(0)*oncov(0)
!
!- cal of square of Q -!
        qsq=0.d0
        qsq=qcov(0)*qcon(0)+qcov(1)*qcon(1) &
           +qcov(2)*qcon(2)+qcov(3)*qcon(3)
!- cal of Q^2+Q*n -!
        qtsq=qsq+(qdotn*qdotn)
!
!        qsq1=0.d0
!        qsq1=qcov(1)*qcon(1)+qcov(2)*qcon(2)+qcov(3)*qcon(3)
!
!
!- Calculate W from last timestep and use for inital guess -!      
!
!- calculation of 4-velocity square 
        utsq=0.d0
!        vsq=0.d0
        do m=1,3
          do n=1,3
            utsq=utsq+gcov1(m,n)*util(m)*util(n)
!            vsq=vsq+gcov1(m,n)*util(m)*util(n)
          enddo
        enddo
!
        if(utsq .lt. 0.d0 .and. abs(utsq) .le. safty) then
          utsq=abs(utsq)
        endif
        if(utsq .lt. 0.d0 .or. abs(utsq) .gt. saftybig) then
          utsq=saftybig
        endif 
!- calculation of Lorentz factor -!
        gfsq=1.+utsq
        gf=sqrt(gfsq)
!        gf=1./sqrt(1.-vsq)
!        gfsq=gf**2
        vsq=(gfsq-1.)*(1./gfsq)
!
!- calculation of W -!      
!
!        if(vsq .lt. 0.9d0 .and. bsq/roh .lt. 1.d0) then 
! 
!          ro=dro/gf
!          pp=kpol*(ro**gam)
!          roh=ro+gam1*pp
!          wg=abs(gfsq*roh)
!
!        else
          
!          if(igue .eq. 1) then
!            vsq1=1.d0
!          else
            vsq1=vsq
!          endif
!            
          tmp1a=4.-vsq1
          tmp1b=4.*(bsq+qdotn)
          tmp1c=qtsq+bsq*bsq+2.*bsq*qdotn
          tmp1d=(tmp1b*tmp1b)-4.*tmp1a*tmp1c
          if(tmp1d .le. 0.d0) then
            tmp1d=abs(tmp1d)
          endif
          wgp=(-tmp1b+sqrt(tmp1d))*(1./(2.*tmp1a))
          wgm=(-tmp1b-sqrt(tmp1d))*(1./(2.*tmp1a))
!
         if(tmp1d .le. 0.d0) then
!           write(*,*) "inside sqrt < 0 for cal wg", tmp1b, tmp1c
          endif
!
          if(wgp .gt. 0.d0) then
            wg=wgp
          elseif(wgm .ge. 0.d0) then
            wg=wgm
          endif
          
!        endif
!
! Make sure that W is large enough so that v^2 < 1
!       
95      continue
  
        pp=kpol*((dro*sqrt(gf))**gam)
       
        tmp3a=(wg*wg*wg)*(wg+2.*bsq)-qdotbsq*(2.*wg+bsq)
        tmp3b=(wg*wg)*(qtsq-bsq*bsq)
!
        tmp3c=(wg*wg*wg)+(wg*wg)*(0.5*bsq+qdotn-pp)-0.5*qdotbsq
!        
!        do nnn=1,iter
          if(tmp3a .le. tmp3b) then
!             write(*,*) "added wg at i,j,k=",i,j,k
             wg= abs(1.5*wg)
             go to 95
!             tmp3a=(wg*wg*wg)*(wg+2.*bsq)-qdotbsq*(2.*wg+bsq)
!             tmp3b=(wg*wg)*(qtsq-bsq*bsq)
          endif

!        enddo
 
        if(tmp3c .gt. 0.d0) then
          wg=abs(0.9*wg)
          go to 95 
        endif 
!
! Calculate W & v^2 for initial guess
!
        wsq=wg*wg
        xsq=(bsq+wg)*(bsq+wg)
        vsq=abs((wsq*qtsq+qdotbsq*(bsq+2.*wg))*(1./(wsq*xsq)))
!       
        if(vsq .gt. 1.d0) then
          write(6,*) 'vsq > 1.0 at initial guess', i, j, k 
          vsq=1.d0-dv
        endif
!
!=====================================================================@
!- Newton Rapson calculation -!

        xxn=vsq
        xo=vsq

!- Newton-Raphson routine start -!
         
        do nnn=1,iter

          if(abs(delta) .ge. deltaend .and. irecov .eq. 0) then
           
            if(nnn .eq. iter) then
              irecov=1

            else
              vsq=xxn
!               
!- calculation of W -!
!
              gtmp=1.-vsq

              if(gtmp .lt. safty) then
                gtmp=safty
                vsq=1.-gtmp 
              endif
              
              gf=1./sqrt(gtmp)
              gfsq=gf*gf
              
              ro=dro*(1./gf)
              pp=kpol*(ro**gam)
              roh=ro+gam1*pp
              wg=abs(gfsq*roh)
              
              wsq=wg*wg
              w3=wsq*wg
!
!-calculation of dW/dv^2 -!
!
              dgfdvsq=0.5*(gtmp)**(-3./2.)
!              dwdvsq=2.*gf*roh*dgfdvsq
              dwdvsq=dro*dgfdvsq &
                    +gam1*kpol*(dro**gam)*(2.-gam)*(gf**(1.-gam))*dgfdvsq
!              
!- calculation of f and df -!
!              
              f=qtsq-vsq*(wg+bsq)*(wg+bsq)+(qdotbsq*(bsq+2.*wg))*(1./wsq)
              tmp1a=-(wg+bsq)*(wg+bsq)
              tmp1b=-2.*vsq*(wg+bsq)
              tmp1c=-2.*qdotbsq*(bsq+2.*wg)*(1./w3)
              tmp1d=2.*qdotbsq*(1./wsq)
              df=tmp1a+(tmp1b+tmp1c+tmp1d)*dwdvsq
!
!- calculation of delta -!
!              
              delta=f/df
! 
              xo=vsq
              xxn=vsq-delta
!
              if(xxn .ge. 1.d0) then
!                irecov=2
                xxn=1.-dv
              endif
!
              if(xxn .lt. 0.d0) then
                xxn=abs(xxn)
              endif
              
            endif
          endif

        enddo
!
!- Newton-Raphson routine end -!
!
        if(irecov .eq. 0) then
          vsq=xxn
          gtmp=1.-vsq
          
          if(gtmp .lt. safty) then
            gtmp=safty
            vsq=1.-gtmp 
          endif
          
          gf=1./sqrt(gtmp)
          gfsq=gf*gf
              
          ro=dro*(1./gf)
          pp=kpol*(ro**gam)
          roh=ro+gam1*pp

          wg=abs(gfsq*roh) 
            
          if(ro .lt. 0.d0) then
            write(*,*) 'negative density', ro, 'at i, j, k', i, j ,k
            ro=dmin*(x1(i)**(-3./2.))
!            ro=dmin
          endif
          if (pp .lt. 0.d0) then
!            irecov=3
            write(*,*) 'negative pressure', pp, 'at i, j, k', i, j ,k
!            pp=(gam-1.)*pmin*(x1(i)**(-5./2.))
            pp=kpol*(dmin**gam)          
          endif

!          if(irecov .ne. 3) then
            qtcon(1)=qcon(1)+oncon(1)*qdotn
            qtcon(2)=qcon(2)+oncon(2)*qdotn
            qtcon(3)=qcon(3)+oncon(3)*qdotn

            util(1)=(gf*(1./(wg+bsq)))*(qtcon(1)+qdotb*bcon(1)*(1./wg))
            util(2)=(gf*(1./(wg+bsq)))*(qtcon(2)+qdotb*bcon(2)*(1./wg))
            util(3)=(gf*(1./(wg+bsq)))*(qtcon(3)+qdotb*bcon(3)*(1./wg))

!            vcon(1)=(1./(xxn+bsq))*(qtcon(1)+qdotb*bcon(1)/wg)
!            vcon(2)=(1./(xxn+bsq))*(qtcon(2)+qdotb*bcon(2)/wg)
!            vcon(3)=(1./(xxn+bsq))*(qtcon(3)+qdotb*bcon(3)/wg)
!           endif
        endif
                 
!        wg=xxn  

!---------------------------------------------------------------------
!- Error message -!                
        if(irecov .eq. 1) then
          write(6,*) ' >> Not convergence in recov4' &
                    ,' at i, j, k:',i, j, k
          write(6,*) ' >> delta =', delta &
                    ,', at x1, x2, x3:',x1(i),x2(j),x3(k)
         
!        elseif(irecov .eq. 2) then
!          write(6,*) ' >> wg  < 0 in recov4' &
!                    ,' at i, j, k:',i, j, k
!          write(6,*) ' >> wg', wg &
!                    ,', at x1, x2, x3:',x1(i),x2(j),x3(k)

!         elseif(irecov .eq. 3) then
!           write(6,*) ' >> ro or pr < 0 in recov4' &
!                    ,' at i, j, k:',i, j, k
!           write(6,*) ' >> ro, pr', ro, pp &
!                     ,' at x1, x2, x3: ',x1(i),x2(j),x3(k)
        endif
!
!---------------------------------------------------------------------@
!   Set Primitive variables
!---------------------------------------------------------------------@
!
        if(irecov .eq. 0) then
         
          uri(1,i,j,k)=ro
          uri(2,i,j,k)=util(1)
          uri(3,i,j,k)=util(2)
          uri(4,i,j,k)=util(3)
!          uri(2,i,j,k)=vcon(1)
!          uri(3,i,j,k)=vcon(2)
!          uri(4,i,j,k)=vcon(3)
          uri(5,i,j,k)=pp
          uri(6,i,j,k)=pp*(1./(ro**(gam-1.)))
          uri(7,i,j,k)=uu(7,i,j,k)*(1./detg1)
          uri(8,i,j,k)=uu(8,i,j,k)*(1./detg1)
          uri(9,i,j,k)=uu(9,i,j,k)*(1./detg1)
        
        else
         
          uri(1,i,j,k)=ro1
          uri(2,i,j,k)=util1(1)
          uri(3,i,j,k)=util1(2)
          uri(4,i,j,k)=util1(3)
!          uri(2,i,j,k)=vcon1(1)
!          uri(3,i,j,k)=vcon1(2)
!          uri(4,i,j,k)=vcon1(3)
          uri(5,i,j,k)=pp1
          uri(6,i,j,k)=pp1*(1./(ro1**(gam-1.)))
          uri(7,i,j,k)=uu(7,i,j,k)*(1./detg1)
          uri(8,i,j,k)=uu(8,i,j,k)*(1./detg1)
          uri(9,i,j,k)=uu(9,i,j,k)*(1./detg1)
         
        endif
!
      enddo
    enddo
  enddo
!
!=====================================================================@
!
  deallocate( gcov1, gcon1, qcon, qcov, bcon, bcov, &
              oncov, oncon, util, util1, qtcon, stat=merr )
!
  return
end subroutine recov4
!
!---------------------------------------------------------------------@
subroutine recov5(uu,uri,gcov,gcon,detg,x1,x2,x3,nm1,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!- 1 variable inversion procedure using S-!
  
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos, iter, dmin, pmin, &
                   kpol, akm
  implicit none
!
  integer :: i, j, k, m, n, nnn, nm1, is1, ie1, js1, je1, ks1, ke1
  integer :: merr, igue, irecov

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)

  real(8) :: x1(imax), x2(jmax), x3(kmax)

  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: qcon(:), qcov(:), bcon(:), bcov(:), &
                          oncov(:), oncon(:) 
  real(8), allocatable :: util(:), util1(:), qtcon(:), ucon(:)
      
  real(8) :: safty, saftybig, gam1, gam2, dv, alpha1, detg1 
  real(8) :: dro, ro, pp, ro1, pp1, roh, roe, bsq, qdotb, qdotbsq, qdotn, &
             qsq, qtsq, utsq, gfsq, gf, vsq, wg, wsq, qsq1, xsq, &
             delta, deltaend, df, f, xxn, xo, w3, gtmp, &
             dgfdvsq, dwdvsq, gen, ent, ent1 
  real(8) :: tmp1a, tmp1b, tmp1c, tmp1d, wgp, wgm, vsq1, tmp3a, tmp3b, tmp3c, &
             tmp4a
  real(8) :: rin, rbh1
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), stat=merr)
  allocate(qcon(0:3), qcov(0:3), bcon(0:3), bcov(0:3), &
           oncov(0:3), oncon(0:3), util(1:3), util1(1:3), qtcon(1:3), &
           ucon(0:3), stat=merr)
!
!--- define variable of safe number ----
      safty = 1.0d-6
      saftybig = 1.0d6
      dv=1.0d-5
      
      gam1=gam/(gam-1.0)
      gam2=1.0/gam1
      wg=0.d0

      igue=1
!
      rin=15.d0
      rbh1=1.+sqrt(1.-akm*akm)
!      
!=====================================================================@
!
! ----  Calculation of variables ---
!
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
      do i=is1+nm1,ie1-nm1
!
!- set parameter for Newton Raphson -!
!
        delta=1.d0
        deltaend = 1.d-8
        df=1.d0
        f=1.d0
!        irecov=0
! 
!- check the position within the BH radius -!
!
!        if(x1(i) .le. rbh1) then
!          irecov=1
!        else
          irecov=0
!        endif
!        
!- copy of metric terms -!
!
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo
!
        alpha1=1./sqrt(-gcon1(0,0))
!
        detg1=detg(i,j,k)
!
!- rearrange of conserved variables -!       
!       
        dro=alpha1*uu(1,i,j,k)*(1./detg1)

!        qcov(0)=alpha1*(uu(5,i,j,k)+uu(1,i,j,k))*(1./detg1)
        qcov(0)=alpha1*(uu(5,i,j,k)-uu(1,i,j,k))*(1./detg1)
!        qcov(0)=alpha1*uu(5,i,j,k)*(1./detg1)
        qcov(1)=alpha1*uu(2,i,j,k)*(1./detg1)
        qcov(2)=alpha1*uu(3,i,j,k)*(1./detg1)
        qcov(3)=alpha1*uu(4,i,j,k)*(1./detg1)

        gen=alpha1*uu(6,i,j,k)*(1./detg1)

        tmp4a=alpha1*uu(5,i,j,k)*(1./detg1)
        
        bcon(0)=0.d0
        bcon(1)=alpha1*uu(7,i,j,k)*(1./detg1)
        bcon(2)=alpha1*uu(8,i,j,k)*(1./detg1)
        bcon(3)=alpha1*uu(9,i,j,k)*(1./detg1)
!        
!- copy of primitive variables (previous time)         
!         
        ro=uri(1,i,j,k)
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
!        vcon(1)=uri(2,i,j,k)
!        vcon(2)=uri(3,i,j,k)
!        vcon(3)=uri(4,i,j,k)
        pp=uri(5,i,j,k)
        ent=uri(6,i,j,k)
        
        ro1=uri(1,i,j,k)
        util1(1)=uri(2,i,j,k)
        util1(2)=uri(3,i,j,k)
        util1(3)=uri(4,i,j,k)
!        vcon1(1)(i,j,k)=uri(2,i,j,k)
!        vcon1(2)(i,j,k)=uri(3,i,j,k)
!        vcon1(3)(i,j,k)=uri(4,i,j,k)
        pp1=uri(5,i,j,k)
        ent1=uri(6,i,j,k)
!
!        roh=ro+gam1*pp

!- cal of covariant 3-magnetic field -! 
        call lower(bcon,bcov,gcov1)
!- cal of contravariant 4-q -!
        call upper(qcov,qcon,gcon1)
!- cal of B-field square (3-vector) -!
        bsq=bcon(1)*bcov(1)+bcon(2)*bcov(2)+bcon(3)*bcov(3)
!- cal of q*B -!
        qdotb=qcov(0)*bcon(0)+qcov(1)*bcon(1)+qcov(2)*bcon(2)+qcov(3)*bcon(3)
!- cal of (Q*B)^2 -!
        qdotbsq=qdotb*qdotb
!- set of covariant 4-vector n -!
        oncov(0)=-1.d0*alpha1
        oncov(1)=0.d0
        oncov(2)=0.d0
        oncov(3)=0.d0
!- cal of contravariant 4-vector n -!
        call upper(oncov,oncon,gcon1)
!- cal of Q*n -!
        qdotn=qcon(0)*oncov(0)
!
!- cal of square of Q -!
        qsq=0.d0
        qsq=qcov(0)*qcon(0)+qcov(1)*qcon(1) &
           +qcov(2)*qcon(2)+qcov(3)*qcon(3)
!- cal of Q^2+Q*n -!
        qtsq=qsq+(qdotn*qdotn)
!
!        qsq1=0.d0
!        qsq1=qcov(1)*qcon(1)+qcov(2)*qcon(2)+qcov(3)*qcon(3)
!
!
!- Calculate W from last timestep and use for inital guess -!      
!
!- calculation of 4-velocity square 
        utsq=0.d0
!        vsq=0.d0
        do m=1,3
          do n=1,3
            utsq=utsq+gcov1(m,n)*util(m)*util(n)
!            vsq=vsq+gcov1(m,n)*util(m)*util(n)
          enddo
        enddo
!
        if(utsq .lt. 0.d0 .and. abs(utsq) .le. safty) then
          utsq=abs(utsq)
        endif
        if(utsq .lt. 0.d0 .or. abs(utsq) .gt. saftybig) then
          utsq=saftybig
        endif 
!- calculation of Lorentz factor -!
        gfsq=1.+utsq
        gf=sqrt(gfsq)
!        gf=1./sqrt(1.-vsq)
!        gfsq=gf**2
        vsq=(gfsq-1.)*(1./gfsq)
!
!- calculation of W -!
!        
!        if(vsq .lt. 0.9d0 .and. bsq/roh .lt. 1.d0) then 
! 
!          ro=dro/gf
!          pp=(gen/gf)*ro**(gam-1.)
!          roh=ro+gam1*pp
!          wg=abs(gfsq*roh)
!
!        else
          
!          if(igue .eq. 1) then
!            vsq1=1.d0
!          else
            vsq1=vsq
!          endif
!            
          tmp1a=4.-vsq1
          tmp1b=4.*(bsq+qdotn)
          tmp1c=qtsq+bsq*bsq+2.*bsq*qdotn
          tmp1d=(tmp1b*tmp1b)-4.*tmp1a*tmp1c
          if(tmp1d .le. 0.d0) then
            tmp1d=abs(tmp1d)
          endif
          wgp=(-tmp1b+sqrt(tmp1d))/(2.*tmp1a)
          wgm=(-tmp1b-sqrt(tmp1d))/(2.*tmp1a)
!
         if(tmp1d .le. 0.d0) then
!           write(*,*) "inside sqrt < 0 for cal wg", tmp1b, tmp1c
          endif
!
          if(wgp .gt. 0.d0) then
            wg=wgp
          elseif(wgm .ge. 0.d0) then
            wg=wgm
          endif
          
!        endif
!
! Make sure that W is large enough so that v^2 < 1
!       
95      continue

        if(ieos .eq. 0) then !- gamma-law EoS -!
          pp=gam2*(wg*gf-dro*sqrt(gf))
        elseif(ieos .eq. 3) then !- polytropic EoS -!   
          pp=kpol*((dro*sqrt(gf))**gam)
        endif
       
        tmp3a=(wg*wg*wg)*(wg+2.*bsq)-qdotbsq*(2.*wg+bsq)
        tmp3b=(wg*wg)*(qtsq-bsq*bsq)
!
        tmp3c=(wg*wg*wg)+(wg*wg)*(0.5*bsq+qdotn-pp)-0.5*qdotbsq
!        
!        do nnn=1,iter
          if(tmp3a .le. tmp3b) then
!             write(*,*) "added wg at i,j,k=",i,j,k
             wg= abs(1.5*wg)
             go to 95
!             tmp3a=(wg*wg*wg)*(wg+2.*bsq)-qdotbsq*(2.*wg+bsq)
!             tmp3b=(wg*wg)*(qtsq-bsq*bsq)
          endif

!        enddo
 
        if(tmp3c .gt. 0.d0) then
          wg=abs(0.9*wg)
          go to 95 
        endif 
!
! Calculate W & v^2 for initial guess
!
        wsq=wg*wg
        xsq=(bsq+wg)*(bsq+wg)
        vsq=abs((wsq*qtsq+qdotbsq*(bsq+2.*wg))*(1./(wsq*xsq)))
!       
        if(vsq .gt. 1.d0) then
          write(6,*) 'vsq > 1.0 at initial guess', i, j, k 
          vsq=1.d0-dv
        endif
!
!=====================================================================@
!- Newton Rapson calculation -!

        xxn=vsq
        xo=vsq

!- Newton-Raphson routine start -!
         
        do nnn=1,iter

          if(abs(delta) .ge. deltaend .and. irecov .eq. 0) then
           
            if(nnn .eq. iter) then
              irecov=1

            else
              vsq=xxn
!               
!- calculation of W -!
!              
              gtmp=1.-vsq

              if(gtmp .lt. safty) then
                gtmp=safty
                vsq=1.-gtmp 
              endif
           
              gf=1./sqrt(gtmp)
              gfsq=gf*gf
              
              ro=dro*(1./gf)
              pp=(gen*(1./gf))*(ro**(gam-1.))
              roh=ro+gam1*pp
              wg=abs(gfsq*roh)
              
              wsq=wg*wg
              w3=wsq*wg
!
!-calculation of dW/dv^2 -!
!
              dgfdvsq=0.5*(gtmp)**(-3./2.)
!              dwdvsq=2.*gf*roh*dgfdvsq
              dwdvsq=dro*dgfdvsq &
                   +gam1*gen*(dro**(gam-1.)) &
                   *(2.-gam)*(gf**(1.-gam))*dgfdvsq
!              
!- calculation of f and df -!
!              
              f=qtsq-vsq*(wg+bsq)*(wg+bsq)+(qdotbsq*(bsq+2.*wg))*(1./wsq)
              tmp1a=-(wg+bsq)*(wg+bsq)
              tmp1b=-2.*vsq*(wg+bsq)
              tmp1c=-2.*qdotbsq*(bsq+2.*wg)*(1./w3)
              tmp1d=2.*qdotbsq*(1./wsq)
              df=tmp1a+(tmp1b+tmp1c+tmp1d)*dwdvsq
!
!- calculation of delta -!
!              
              delta=f/df
! 
              xo=vsq
              xxn=vsq-delta
!
              if(xxn .ge. 1.d0) then
!                irecov=2
                xxn=1.-dv
              endif
!
              if(xxn .le. 0.d0) then
                xxn=abs(xxn)
              endif
!             
            endif
          endif

        enddo
!
!- Newton-Raphson routine end -!
!
        if(irecov .eq. 0) then
          vsq=xxn
          gtmp=1.-vsq

          if(gtmp .lt. safty) then
            gtmp=safty 
          endif
           
          gf=1./sqrt(gtmp)
          gfsq=gf*gf
              
          ro=dro*(1./gf)
          pp=(gen*(1./gf))*(ro**(gam-1.))
          roh=ro+gam1*pp

          wg=abs(gfsq*roh) 
            
          if(ro .lt. 0.d0) then
            write(*,*) 'negative density', ro, 'at i, j, k', i, j ,k
            ro=dmin*(x1(i)**(-3./2.))
!            ro=dmin
          endif
          if (pp .lt. 0.d0) then
!            irecov=3
            write(*,*) 'negative pressure', pp, 'at i, j, k', i, j ,k
            if(ieos .eq. 0) then
              pp=pmin*(gam-1.)*(x1(i)**(-5./2.))
!              pp=pmin 
            elseif(ieos .eq. 3) then
              pp=kpol*(dmin**gam)
            endif
!            pp=(gam-1.)*pmin
          endif

!          if(irecov .ne. 3) then
            qtcon(1)=qcon(1)+oncon(1)*qdotn
            qtcon(2)=qcon(2)+oncon(2)*qdotn
            qtcon(3)=qcon(3)+oncon(3)*qdotn

            util(1)=(gf/(wg+bsq))*(qtcon(1)+qdotb*bcon(1)*(1./wg))
            util(2)=(gf/(wg+bsq))*(qtcon(2)+qdotb*bcon(2)*(1./wg))
            util(3)=(gf/(wg+bsq))*(qtcon(3)+qdotb*bcon(3)*(1./wg))

!            vcon(1)=(1./(xxn+bsq))*(qtcon(1)+qdotb*bcon(1)/wg)
!            vcon(2)=(1./(xxn+bsq))*(qtcon(2)+qdotb*bcon(2)/wg)
!            vcon(3)=(1./(xxn+bsq))*(qtcon(3)+qdotb*bcon(3)/wg)
!           endif
        endif
                 
        wg=xxn  

!---------------------------------------------------------------------
!- Error message -!                
        if(irecov .eq. 1) then
          write(6,*) ' >> Not convergence in recov5' &
                    ,' at i, j, k:',i, j, k
          write(6,*) ' >> delta =', delta &
                    ,', at x1, x2, x3:',x1(i),x2(j),x3(k)
         
!        elseif(irecov .eq. 2) then
!          write(6,*) ' >> wg  < 0 in recov3' &
!                    ,' at i, j, k:',i, j, k
!          write(6,*) ' >> wg', wg &
!                    ,', at x1, x2, x3:',x1(i),x2(j),x3(k)

!         elseif(irecov .eq. 3) then
!           write(6,*) ' >> ro or pr < 0 in recov1d' &
!                    ,' at i, j, k:',i, j, k
!           write(6,*) ' >> ro, pr', ro, pp &
!                     ,' at x1, x2, x3: ',x1(i),x2(j),x3(k)
        endif
!
!---------------------------------------------------------------------@
!   Set Primitive variables
!---------------------------------------------------------------------@
!
        if(irecov .eq. 0) then
         
          uri(1,i,j,k)=ro
          uri(2,i,j,k)=util(1)
          uri(3,i,j,k)=util(2)
          uri(4,i,j,k)=util(3)
!          uri(2,i,j,k)=vcon(1)
!          uri(3,i,j,k)=vcon(2)
!          uri(4,i,j,k)=vcon(3)
          uri(5,i,j,k)=pp
          uri(6,i,j,k)=pp*(1./(ro**(gam-1.)))
          uri(7,i,j,k)=uu(7,i,j,k)*(1./detg1)
          uri(8,i,j,k)=uu(8,i,j,k)*(1./detg1)
          uri(9,i,j,k)=uu(9,i,j,k)*(1./detg1)
!        
!- cal of contravariant 4-velocity -!
          call calucon(util,ucon,gf,gcon1)
!- update conserved enthalpy -!
          uu(6,i,j,k)=detg1*((pp*(1./(ro**(gam-1.))))*ucon(0))
!          
        else
         
          uri(1,i,j,k)=ro1
          uri(2,i,j,k)=util1(1)
          uri(3,i,j,k)=util1(2)
          uri(4,i,j,k)=util1(3)
!          uri(2,i,j,k)=vcon1(1)
!          uri(3,i,j,k)=vcon1(2)
!          uri(4,i,j,k)=vcon1(3)
          uri(5,i,j,k)=pp1
          uri(6,i,j,k)=pp1/(ro1**(gam-1.))
          uri(7,i,j,k)=uu(7,i,j,k)/detg1
          uri(8,i,j,k)=uu(8,i,j,k)/detg1
          uri(9,i,j,k)=uu(9,i,j,k)/detg1
!        
!- cal of contravariant 4-velocity -!
          call calucon(util,ucon,gf,gcon1)
!- update conserved enthalpy -!
          uu(6,i,j,k)=detg1*uri(6,i,j,k)*ucon(0)
!         
        endif
!
      enddo
    enddo
  enddo
!
!=====================================================================@
!
  deallocate( gcov1, gcon1, qcon, qcov, bcon, bcov, &
              oncov, oncon, util, util1, qtcon, ucon, stat=merr )
!
  return
end subroutine recov5
!---------------------------------------------------------------------@
subroutine recov6(uu,uri,gcov,gcon,detg,x1,x2,x3,nm1,is1,ie1,js1,je1,ks1,ke1)
!---------------------------------------------------------------------@
!- 1 variable (W) inversion procedure using enthalpy conservation -!
!- (ideal EoS only)
  
  use pram, only : imax, jmax, kmax, nv, gam, c0, ieos, iter, dmin, pmin, &
                   kpol, akm, pi
  implicit none
!
  integer :: i, j, k, m, n, nnn, nm1, is1, ie1, js1, je1, ks1, ke1
  integer :: merr, igue, irecov

  real(8) :: uu(nv,is1:ie1,js1:je1,ks1:ke1)
  real(8) :: uri(nv,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: gcon(0:3,0:3,is1:ie1,js1:je1,ks1:ke1), &
             gcov(0:3,0:3,is1:ie1,js1:je1,ks1:ke1)

  real(8) :: detg(is1:ie1,js1:je1,ks1:ke1)

  real(8) :: x1(imax), x2(jmax), x3(kmax)

  real(8), allocatable :: gcov1(:,:), gcon1(:,:)
  real(8), allocatable :: qcon(:), qcov(:), bcon(:), bcov(:), &
                          oncov(:), oncon(:) 
  real(8), allocatable :: util(:), util1(:), qtcon(:)
      
  real(8) :: safty, saftybig, gam1, gam2, dv, alpha1, detg1 
  real(8) :: dro, ro, pp, ro1, pp1, roh, roe, bsq, qdotb, qdotbsq, qdotn, &
             qsq, qtsq, utsq, gfsq, gf, vsq, wg, wsq, xsq, qsq1, &
             delta, deltaend, df, f, xxn, xo, w2, w3, gtmp, ptmp, &
             dpdw, dpdw2, dpdvsq, dvsqdw, gss, ss, ss1, drodw 
  real(8) :: tmp1a, tmp1b, tmp1c, tmp1d, wgp, wgm, vsq1, tmp3a, tmp3b, &
             tmp3c, tmp4a
  real(8) :: rin, rbh1, rr, th
!
!- allocate variables -!
  allocate(gcov1(0:3,0:3), gcon1(0:3,0:3), stat=merr)
  allocate(qcon(0:3), qcov(0:3), bcon(0:3), bcov(0:3), &
           oncov(0:3), oncon(0:3), util(1:3), util1(1:3), qtcon(1:3), stat=merr)
!
!--- define variable of safe number ----
      safty = 1.0d-6
      saftybig = 1.0d6
      dv=1.0d-2
!      
!--- Parameters ---!
!      
      gam1=gam/(gam-1.0)
      gam2=1.0/gam1
      wg=0.d0

      igue=1
!
      rin=15.d0
      rbh1=1.+sqrt(1.-akm*akm)
!
!
!=====================================================================@
!
! ----  Calculation of variables ---
!
  do k=ks1+nm1,ke1-nm1
    do j=js1+nm1,je1-nm1
       do i=is1+nm1,ie1-nm1
!
!- set parameter for Newton Raphson -!
!
      delta=1.d0
      deltaend = 1.d-6
      df=1.d0
      f=1.d0
!      irecov=0
! 
!- check the position within the BH radius -!
!
!        if(x1(i) .le. rbh1) then
!          irecov=1
!        else
          irecov=0
!        endif   
!          
!- copy of metric terms -!
!
        do m=0,3
          do n=0,3
            gcov1(m,n)=gcov(m,n,i,j,k)
            gcon1(m,n)=gcon(m,n,i,j,k)
          enddo
        enddo
!
        alpha1=1./sqrt(-gcon1(0,0))
!
        detg1=detg(i,j,k)
!
!- rearrange of conserved variables -!       
!      
        dro=alpha1*uu(1,i,j,k)*(1./detg1)

!        qcov(0)=alpha1*(uu(5,i,j,k)+uu(1,i,j,k))*(1./detg1)
        qcov(0)=alpha1*(uu(5,i,j,k)-uu(1,i,j,k))*(1./detg1)
!        qcov(0)=alpha1*uu(5,i,j,k)*(1./detg1)
        qcov(1)=alpha1*uu(2,i,j,k)*(1./detg1)
        qcov(2)=alpha1*uu(3,i,j,k)*(1./detg1)
        qcov(3)=alpha1*uu(4,i,j,k)*(1./detg1)

        gss=alpha1*uu(6,i,j,k)*(1./detg1)

        tmp4a=alpha1*uu(5,i,j,k)*(1./detg1)
        
        bcon(0)=0.d0
        bcon(1)=alpha1*uu(7,i,j,k)*(1./detg1)
        bcon(2)=alpha1*uu(8,i,j,k)*(1./detg1)
        bcon(3)=alpha1*uu(9,i,j,k)*(1./detg1)
!        
!- copy of primitive variables (previous time)         
!         
        ro=uri(1,i,j,k)
        util(1)=uri(2,i,j,k)
        util(2)=uri(3,i,j,k)
        util(3)=uri(4,i,j,k)
!        vcon(1)=uri(2,i,j,k)
!        vcon(2)=uri(3,i,j,k)
!        vcon(3)=uri(4,i,j,k)
        pp=uri(5,i,j,k)
        ss=uri(6,i,j,k)
!
!- store for previous time        
!        
        ro1=uri(1,i,j,k)
        util1(1)=uri(2,i,j,k)
        util1(2)=uri(3,i,j,k)
        util1(3)=uri(4,i,j,k)
!        vcon1(1)(i,j,k)=uri(2,i,j,k)
!        vcon1(2)(i,j,k)=uri(3,i,j,k)
!        vcon1(3)(i,j,k)=uri(4,i,j,k)
        pp1=uri(5,i,j,k)
        ss1=uri(6,i,j,k)
!
        if(ieos .eq. 0 .or. ieos .eq. 3) then
          roh=ro+(gam/(gam-1.0))*pp
        elseif(ieos .eq. 1) then
          roh=(5./2.)*pp &
              +sqrt((9./4.)*pp**2+ro**2)
        elseif(ieos .eq. 2) then
          roe=(3./2.)*(pp+((3.*pp**2) &
              /(2.0*ro+sqrt(2.*pp**2+4.*ro**2)) ))
          roh=ro+roe+pp
        endif
!
!- check B-field minimum -!
!
!        if(abs(bcon(1)) .lt. bmin) then
!          bcon(1)=0.d0
!        endif
!        if(abs(bcon(2)) .lt. bmin) then
!          bcon(2)=0.d0
!        endif       
!        if(abs(bcon(3)) .lt. bmin) then
!          bcon(3)=0.d0
!        endif   
        
!- cal of covariant 3-magnetic field -! 
        call lower(bcon,bcov,gcov1)
!- cal of contravariant 4-q -!
        call upper(qcov,qcon,gcon1)
!- cal of B-field square (3-vector) -!
        bsq=bcon(1)*bcov(1)+bcon(2)*bcov(2)+bcon(3)*bcov(3)
!- cal of q*B -!
        qdotb=qcov(0)*bcon(0)+qcov(1)*bcon(1)+qcov(2)*bcon(2)+qcov(3)*bcon(3)
!- cal of (Q*B)^2 -!
        qdotbsq=qdotb*qdotb
!- set of covariant 4-vector n -!
        oncov(0)=-1.d0*alpha1
        oncov(1)=0.d0
        oncov(2)=0.d0
        oncov(3)=0.d0
!- cal of contravariant 4-vector n -!
        call upper(oncov,oncon,gcon1)
!- cal of Q*n -!
        qdotn=qcon(0)*oncov(0)
!
!- cal of square of Q -!
        qsq=0.d0
        qsq=qcov(0)*qcon(0)+qcov(1)*qcon(1) &
           +qcov(2)*qcon(2)+qcov(3)*qcon(3)
!- cal of Q^2+Q*n -!
        qtsq=qsq+(qdotn*qdotn)
!
!        qsq1=0.d0
!        qsq1=qcov(1)*qcon(1)+qcov(2)*qcon(2)+qcov(3)*qcon(3)
!
!
!- Calculate W from last timestep and use for inital guess -!      
!
!- calculation of 4-velocity square 
        utsq=0.d0
!        vsq=0.d0
        do m=1,3
          do n=1,3
            utsq=utsq+gcov1(m,n)*util(m)*util(n)
!            vsq=vsq+gcov1(m,n)*util(m)*util(n)
          enddo
        enddo
!
        if(utsq .lt. 0.d0 .and. abs(utsq) .le. safty) then
          utsq=abs(utsq)
        endif
        if(utsq .lt. 0.d0 .or. abs(utsq) .gt. saftybig) then
          utsq=saftybig
        endif 
!- calculation of Lorentz factor -!
        gfsq=1.+utsq
        gf=sqrt(gfsq)
!        gf=1./sqrt(1.-vsq)
!        gfsq=gf**2
        vsq=(gfsq-1.)*(1./gfsq)
!- calculation of W -!      
!        if(vsq .lt. 0.9d0 .and. bsq/roh .lt. 1.d0) then 
! 
!          ro=dro/gf
!          roh=ro+gam1*pp
!          wg=abs(gfsq*roh)
!
!        else
          
!          if(igue .eq. 1) then
!            vsq1=1.d0
!          else
            vsq1=vsq
!          endif
!            
          tmp1a=4.-vsq1
          tmp1b=4.*(bsq+qdotn)
          tmp1c=qtsq+bsq*bsq+2.*bsq*qdotn
          tmp1d=(tmp1b*tmp1b)-4.*tmp1a*tmp1c
          if(tmp1d .le. 0.d0) then
            tmp1d=abs(tmp1d)
          endif
          wgp=(-tmp1b+sqrt(tmp1d))*(1./(2.*tmp1a))
          wgm=(-tmp1b-sqrt(tmp1d))*(1./(2.*tmp1a))
!
         if(tmp1d .le. 0.d0) then
           write(*,*) "inside sqrt < 0 for cal wg", tmp1b, tmp1c
          endif
!
          if(wgp .gt. 0.d0) then
            wg=wgp
          elseif(wgm .ge. 0.d0) then
            wg=wgm
          endif
          
!        endif
!
! Make sure that W is large enough so that v^2 < 1
!       
  95    continue

        if(ieos .eq. 0) then !- gamma-law EoS -!
          pp=gam2*(wg*gf-dro*sqrt(gf))
        elseif(ieos .eq. 3) then !- polytropic EoS -!   
          pp=kpol*((dro*sqrt(gf))**gam)
        endif
!       
        tmp3a=(wg*wg*wg)*(wg+2.*bsq)-qdotbsq*(2.*wg+bsq)
        tmp3b=(wg*wg)*(qtsq-bsq*bsq)
!        
        tmp3c=(wg*wg*wg)+(wg*wg)*(0.5*bsq+qdotn-pp)-0.5*qdotbsq
!        
!        do nnn=1,iter
          if(tmp3a .le. tmp3b) then
!             write(*,*) "added wg at i,j,k=",i,j,k
             wg= abs(1.5*wg)
             go to 95
!             tmp3a=(wg*wg*wg)*(wg+2.*bsq)-qdotbsq*(2.*wg+bsq)
!             tmp3b=(wg*wg)*(qtsq-bsq*bsq)
          endif
!        enddo
!          
         if(tmp3c .gt. 0.d0) then
           wg=abs(0.9*wg)
           go to 95 
         endif 
!
!=====================================================================@
!- Newton Rapson calculation -!

        xxn=wg
        xo=wg

!- Newton-Raphson routine start -!
         
        do nnn=1,iter

          if(abs(delta) .ge. deltaend .and. irecov .eq. 0) then
           
            if(nnn .eq. iter) then
              irecov=1

            else
              wg=xxn

              wsq=wg*wg
              w3=wsq*wg
!
!- calculation of v^2 -!
!
              xsq=(bsq+wg)*(bsq+wg)
              vsq=abs((wsq*qtsq+qdotbsq*(bsq+2.*wg))*(1./(wsq*xsq)))
!       
              if(vsq .gt. 1.d0) then
                write(6,*) 'vsq > 1.0 at', i, j, k 
                vsq=1.d0-dv
              endif
              
              gtmp=1.-vsq

              if(gtmp .lt. safty) then
                gtmp=safty
                vsq=1.-gtmp 
              endif
!             
              gf=1./sqrt(gtmp)
              gfsq=gf*gf
!
              ro=dro*(1./gf)              
!
!- calculation of pressure -!
!              
              if(ieos .eq. 0) then !- gamma-law EoS -!
                ptmp=gam2*(wg*gtmp-dro*sqrt(gtmp))
!                ptmp=(gss/dro))*(ro**(gam))                 
              elseif(ieos .eq. 3) then !- polytropic EoS -!   
                ptmp=kpol*((dro*sqrt(gtmp))**gam)
!                ptmp=(gss/dro))*(ro**(gam)) 
              endif
!
!-calculation of dv^2/dw -!
!
              tmp1a=(bsq+wg)
              tmp1b=(bsq+wg)*(bsq+wg)*(bsq+wg)
              tmp1c=-2.*qtsq*(1./tmp1b)
              tmp1d=-2.*qdotbsq*(3.*wg*tmp1a+bsq*bsq)*(1./(tmp1b*w3))
              dvsqdw=tmp1c+tmp1d

!
!- calculation of dp/dw, dp/dv^2 -!
!              
              if(ieos .eq. 0) then !- gamma-law EoS -!
                dpdw=gam2*gtmp
                dpdvsq=-gam2*wg+0.5*gam2*dro*(1./(sqrt(gtmp)))
                dpdw2=dpdw+dpdvsq*dvsqdw                
              endif
!
!- calculation of drho/dw -!
!              
              drodw=-0.5*dro*gf*dvsqdw
              
!             
!- calculation of f and df -!
!              
              f=dro*ptmp/(ro**(gam))-gss
              df=(dro/(ro**(gam)))*dpdw2-gam*ptmp*(ro**(gam+1.))*drodw
!
!- calculation of delta -!
!              
              delta=f*(1./df)
! 
              xo=wg
              xxn=wg-delta
!
              if(xxn .lt. 0.d0) then
!                irecov=2
                 xxn=abs(xxn)
              endif
!         
              if(xxn .gt. 1.d6) then
                xxn=1.d6
              endif
!              
            endif
          endif

        enddo
!
!- Newton-Raphson routine end -!
!
        if(irecov .eq. 0) then
          wg=xxn
          wsq=wg*wg
          xsq=(bsq+wg)*(bsq+wg)
          vsq=abs((wsq*qtsq+qdotbsq*(bsq+2.*wg))*(1./(wsq*xsq)))
!       
          if(vsq .gt. 1.d0) then
            write(6,*) 'vsq > 1.0 at', i, j, k 
            vsq=1.d0-dv
          endif           
!
          gtmp=1.-vsq

          if(gtmp .lt. safty) then
            gtmp=safty 
          endif
           
          gf=1./sqrt(gtmp)
          gfsq=gf*gf
          
          ro=dro/gf
          
          w2=wg*(gtmp)
          if(ieos .eq. 3) then !-polytropic EoS
             pp=kpol*(ro**gam)
!             pp=(gen*(1./gf))*(ro**(gam-1.))
          elseif(ieos .eq. 0) then !- gamma-law EoS 
             pp=gam2*(w2-ro)
!             pp=(gen*(1./gf))*(ro**(gam-1.))
          endif
            
          if(ro .lt. 0.d0) then
            write(*,*) 'negative density', ro, 'at', i, j ,k
            ro=dmin*(x1(i)**(-3./2.))
!            ro=dmin
          endif
          if (pp .lt. 0.d0) then
!            irecov=3
            write(*,*) 'negative pressure', pp, 'at', i, j ,k
            pp=pmin*(gam-1.)*(x1(i)**(-5./2.))
!            pp=pmin
          endif
          
!          if(irecov .ne. 3) then
            qtcon(1)=qcon(1)+oncon(1)*qdotn
            qtcon(2)=qcon(2)+oncon(2)*qdotn
            qtcon(3)=qcon(3)+oncon(3)*qdotn

            util(1)=(gf*(1./(wg+bsq)))*(qtcon(1)+qdotb*bcon(1)*(1./wg))
            util(2)=(gf*(1./(wg+bsq)))*(qtcon(2)+qdotb*bcon(2)*(1./wg))
            util(3)=(gf*(1./(wg+bsq)))*(qtcon(3)+qdotb*bcon(3)*(1./wg))

!            vcon(1)=(1./(xxn+bsq))*(qtcon(1)+qdotb*bcon(1)/wg)
!            vcon(2)=(1./(xxn+bsq))*(qtcon(2)+qdotb*bcon(2)/wg)
!            vcon(3)=(1./(xxn+bsq))*(qtcon(3)+qdotb*bcon(3)/wg)
!           endif
        endif
                 
        wg=xxn  

!---------------------------------------------------------------------
!- Error message -!                
        if(irecov .eq. 1) then
          write(6,*) ' >> Not convergence in recov3' &
                    ,' at:',i, j, k
          write(6,*) ' >> delta =', delta &
                    ,', at x1, x2, x3:',x1(i),x2(j),x3(k)
         
        elseif(irecov .eq. 2) then
          write(6,*) ' >> wg  < 0 in recov3' &
                    ,' at:',i, j, k
          write(6,*) ' >> wg', wg &
                    ,', at x1, x2, x3:',x1(i),x2(j),x3(k)

!         elseif(irecov .eq. 3) then
!           write(6,*) ' >> ro or pr < 0 in recov1d' &
!                    ,' at i, j, k:',i, j, k
!           write(6,*) ' >> ro, pr', ro, pp &
!                     ,' at x1, x2, x3: ',x1(i),x2(j),x3(k)
        endif
!
!---------------------------------------------------------------------@
!   Set Primitive variables
!---------------------------------------------------------------------@
!
        if(irecov .eq. 0) then
         
          uri(1,i,j,k)=ro
          uri(2,i,j,k)=util(1)
          uri(3,i,j,k)=util(2)
          uri(4,i,j,k)=util(3)
!          uri(2,i,j,k)=vcon(1)
!          uri(3,i,j,k)=vcon(2)
!          uri(4,i,j,k)=vcon(3)
          uri(5,i,j,k)=pp
          uri(6,i,j,k)=ro*pp/(ro**(gam))
!          uri(7,i,j,k)=bcon(1)
!          uri(8,i,j,k)=bcon(2)
!          uri(9,i,j,k)=bcon(3)
!          
          uri(7,i,j,k)=uu(7,i,j,k)*(1./detg1)
          uri(8,i,j,k)=uu(8,i,j,k)*(1./detg1)
          uri(9,i,j,k)=uu(9,i,j,k)*(1./detg1)
!
        rr=x1(i)
        th=x3(k)*180.d0/pi 
!        if(x1(i) .ge. 35.d0 .and. x1(i) .le. 35.5d0) then
!          if((uri(7,i,j,k)) .gt. 0.d0 .or. abs(uri(9,i,j,k)) .gt. 0.d0) then
!            write(*,*) 'rr, th, br, bph, bth in recov3a=', rr, th, uri(7,i,j,k), uri(8,i,j,k), uri(9,i,j,k)
!          endif           
!        endif
       
        else
         
          uri(1,i,j,k)=ro1
          uri(2,i,j,k)=util1(1)
          uri(3,i,j,k)=util1(2)
          uri(4,i,j,k)=util1(3)
!          uri(2,i,j,k)=vcon1(1)
!          uri(3,i,j,k)=vcon1(2)
!          uri(4,i,j,k)=vcon1(3)
          uri(5,i,j,k)=pp1
          uri(6,i,j,k)=ro1*pp1/(ro1**(gam))
          uri(7,i,j,k)=uu(7,i,j,k)*(1./detg1)
          uri(8,i,j,k)=uu(8,i,j,k)*(1./detg1)
          uri(9,i,j,k)=uu(9,i,j,k)*(1./detg1)
         
        endif
!
      enddo
    enddo
  enddo
!
!=====================================================================@
!
  deallocate( gcov1, gcon1, qcon, qcov, bcon, bcov, &
              oncov, oncon, util, util1, qtcon, stat=merr )
!
  return
end subroutine recov6
!
!---------------------------------------------------------------------@
subroutine inequality(dro,bsq,qtsq,qdotb,qdotn,romin,prmin)  
!---------------------------------------------------------------------@
!- Modification of updated conserved variables before starting inversion
!- procedure by checking inequality
!  
  use pram, only : gam, iter
!    
  implicit none
!  
  integer :: nnn
  real(8) :: tau, taum, taumin, tauatm, taun
  real(8) :: qdotn, dro, bsq, qtsq, qdotb, romin, prmin !- input
  real(8) :: drosq, qdotbsq, rohatm, wo, wm, wmsq, xmsq, qtsqm, &
             tmp1, tmp2, small

!=====================================================================@  
!  
  small=1.d-10
!  
  tau=-qdotn
  drosq=dro*dro
  qdotbsq=qdotb*qdotb
  
  rohatm=romin+(gam/(gam-1.))*prmin
  tauatm=rohatm-prmin
!
!- check inequality for tau-!
!  
  if(tau .lt. 0.5d0*bsq) then
    tau=tauatm+0.5*bsq
    qdotn=-tau 
  endif   

  taumin=tau-0.5*bsq -(bsq*qtsq)/(2.*(dro+bsq)*(dro+bsq))
  tmp1=sqrt(qtsq+drosq) - dro
  
  if(taumin .lt. tmp1) then
!- calculation of Wm     

!- initial guess for wm
    wo=rohatm
    wm=sqrt(qdotbsq -drosq)
!
    do nnn=1, iter 
      if(abs((wo-wm)/wo) .gt. 1.e-10) then
        wo=wm  
        wmsq=wm*wm
        xmsq=(bsq+wm)*(bsq+wm)
        qtsqm=(wmsq*qtsq+qdotbsq*(bsq+2.*wm))/xmsq
        wm=sqrt(qtsqm+drosq)
      endif    
    enddo
   
    taum=tau-0.5*bsq -((bsq*qtsq)-qdotbsq)/(2.*(wm+bsq)*(wm+bsq))

    if(taum .le. tmp1) then
!- set new tau -!
!       taun=tau-taum+tmp1+small
       taun=tau-taum+tmp1+tauatm
       qdotn=-taun
    endif   
!     
  endif 
!  
  return
end subroutine  
