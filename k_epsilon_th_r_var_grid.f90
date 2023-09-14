!!!************************************************************
!!!*						                	                *
!!!*                     K Epsilon 1D      	         	    *
!!!*							        	                    *
!!!*                Author: Dr. Andrea Boghi  		        *
!!!*							        	                    *
!!!*								                            *
!!!************************************************************
    
Program main_K_epsilon
    implicit none
    real*8,allocatable :: r(:),U(:),kt(:),eps(:),detadr(:),d2etadr2(:)
    real*8,allocatable :: U0(:),kt0(:),eps0(:)
    real*8,allocatable :: nut(:),f2(:),dnutdr(:),dUdr(:),dnutdeta(:),dUdeta(:)
    real*8,allocatable :: tau_mu(:),tau_R(:),Pk(:),Tk(:),Dk(:)
    real*8,allocatable :: Peps(:),Teps(:),Deps(:),epseps(:)
    real*8,allocatable :: T(:),Th2(:),T0(:),Th20(:),lambda(:),dlambdadT(:),d2lambdadT2(:)
    real*8,allocatable :: dTdeta(:),dTdr(:),d2Tdeta2(:),d2Tdr2(:),dTh2deta(:),dTh2dr(:),d2Th2deta2(:),d2Th2dr2(:)
    real*8,allocatable :: q_lam(:),q_R(:),q_new(:),P_Th2(:),eps_Th2(:),T_Th2(:),D_Th2(:),H_Th2(:)
    integer j,nr,iter,niter
    real*8 Re_tau,Pr,Bp,Cp,Dp,Ep,dr_min,sigmak,sigmae,Cmu,Ce1,Ce2,f1,alphaU,alphaKt,alphaeps
    real*8 resU,resK,resE,resT,resTh2,deta,aU_w,aU_e,sU,aK_w,aK_e,sK,aE_w,aE_e,sE, conv_fac
    real*8 resU_old,resK_old,resE_old,resT_old,resTh2_old
    real*8 aT_e,aT_w,sT,aTh2_e,aTh2_w,sTh2
    real*8 sigmaT,sigmaTh2,alphaT,alphaTh2,U_max
    CHARACTER(len=80)::fname_ke
    CHARACTER(len=80)::fname_th
    CHARACTER(len=80)::fname_res
    logical flag

    open(1,file='imp_ke_th_var.dat')
    read(1,*) flag
    read(1,*) nr
    read(1,*) niter
    read(1,*) Re_tau
    read(1,*) dr_min
    read(1,*) Pr
    read(1,*) Bp
    read(1,*) Cp
    read(1,*) Dp
    read(1,*) Ep
    read(1,*) alphaU
    read(1,*) alphaKt
    read(1,*) alphaeps 
    read(1,*) alphaT
    read(1,*) alphaTh2
    close(1)

    print*, ' niter =', niter, ' nr =', nr 
    print*, ' Re_tau =', Re_tau, ' flag =', flag 

    allocate(r(1:nr),U(1:nr),kt(1:nr),eps(1:nr),nut(1:nr),dnutdr(1:nr),dUdr(1:nr),f2(1:nr))
    allocate(U0(1:nr),kt0(1:nr),eps0(1:nr),detadr(1:nr),d2etadr2(1:nr),dnutdeta(1:nr),dUdeta(1:nr))
    allocate(T(1:nr),Th2(1:nr),T0(1:nr),Th20(1:nr),lambda(1:nr),dlambdadT(1:nr),d2lambdadT2(1:nr))
    allocate(dTdeta(1:nr),dTdr(1:nr),d2Tdeta2(1:nr),d2Tdr2(1:nr),dTh2deta(1:nr),dTh2dr(1:nr),d2Th2deta2(1:nr),d2Th2dr2(1:nr))
    allocate(q_lam(1:nr),q_R(1:nr),q_new(1:nr),P_Th2(1:nr),eps_Th2(1:nr),T_Th2(1:nr),D_Th2(1:nr),H_Th2(1:nr))

    call initialization(flag,nr,Re_tau,Pr,Bp,Cp,dr_min,r,detadr,d2etadr2,U,Kt,eps,T,Th2,deta)

    call nagano_takawa_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu,f1)

    sigmaT = 1.d0
    sigmaTh2 = 1.d0

    conv_fac = 1.d0

    write(fname_res,110)'residualsNT_var_Re_tau=',int(Re_tau),'_log10(Pr)=',int(log10(Pr)),'.csv'
    110 format(a23,i3,a11,i2,a4)

    open(11,file=fname_res)
    write(11,*) '"iter","resU","resK","resE","resT","resTh2"'

    do iter=1,niter

        resU_old = resU
        resK_old = resK
        resE_old = resE
        resT_old = resT
        resTh2_old = resTh2

        U0 = U
        Kt0 = kt
        eps0 = eps
        T0 = T
        Th20 = Th2

        do j=1,nr
            call thernal_diffusivity(lambda(j),dlambdadT(j),d2lambdadT2(j),T(j),Bp,Cp,Dp,Ep)
            if(lambda(j) < 0) then
                print*, 'Warning: at j= ',j,' lambda<0'
                lambda(j) = -lambda(j)
            endif
        enddo

        call nagano_takawa_k_epsilon_functions(nut,f2,nr,r,kt,eps,Cmu)
        call ddeta(nr,nut,dnutdeta,deta)
        dnutdr = dnutdeta*detadr
        call ddeta(nr,U,dUdeta,deta)
        dUdr = dUdeta*detadr
        call ddeta(nr,T,dTdeta,deta)
        dTdr = dTdeta*detadr
        call ddeta(nr,Th2,dTh2deta,deta)
        dTh2dr = dTh2deta*detadr
        call d2deta2(nr,T,d2Tdeta2,deta)
        d2Tdr2 = d2Tdeta2*detadr**2.d0 + dTdeta*d2etadr2
        call d2deta2(nr,Th2,d2Th2deta2,deta)
        d2Th2dr2 = d2Th2deta2*detadr**2.d0 + dTh2deta*d2etadr2

        U_max = maxval(U, dim=1, mask=(U>0))

        call solve_u(U,nut,dnutdr,r,detadr,d2etadr2,deta,Re_tau,nr) 
        call solve_Kt(Kt,eps,dUdr,nut,dnutdr,r,detadr,d2etadr2,deta,sigmak,nr) 
        call solve_eps(Kt,eps,dUdr,nut,dnutdr,r,detadr,d2etadr2,deta,sigmae,ce1,ce2,f1,f2,nr)
        call solve_T(T,nut,dnutdr,lambda,dlambdadT,d2lambdadT2,dTh2dr,d2Th2dr2,dTdr,r,Re_tau,Pr,sigmaT,deta,d2etadr2,detadr,nr)
        call solve_Th2(Th2,nut,r,dnutdr,lambda,dlambdadT,d2lambdadT2,dTdr,d2Tdr2,Pr,sigmaTh2,deta,d2etadr2,detadr,eps,Kt,nr)

        call residuals(nr,U,U0,resU)
        call residuals(nr,Kt,Kt0,resK)
        call residuals(nr,eps,eps0,resE)
        call residuals(nr,T,T0,resT)
        call residuals(nr,Th2,Th20,resTh2)
        write(11,102) conv_fac*iter,',',resU,',',resK,',',resE,',',resT,',',resTh2

        U = alphaU*U +(1.d0-alphaU)*U0
        Kt = alphaKt*Kt +(1.d0-alphaKt)*Kt0
        eps = alphaeps*eps +(1.d0-alphaeps)*eps0
        !!! T can be negative
        T = alphaT*T +(1.d0-alphaT)*T0
        Th2 = alphaTh2*Th2 +(1.d0-alphaTh2)*Th20
        print*, ' completed =', 100*real(iter)/real(niter), ' resU = ', resU, ' resK = ', resK, ' resE = ', resE, &
        ' resT = ', resT, ' resTh2 = ', resTh2
        
    enddo
    close(11)

    open(512,file='point_ke_var.dat',form='unformatted')
    write(512) r,detadr,d2etadr2,U,Kt,eps,nut,f2
    close(512)

    open(512,file='point_T_var.dat',form='unformatted')
    write(512) r,T,Th2
    close(512)

    allocate(tau_mu(1:nr),tau_R(1:nr),Pk(1:nr),Tk(1:nr),Dk(1:nr))
    allocate(Peps(1:nr),Teps(1:nr),Deps(1:nr),epseps(1:nr))

    call output_fields_k_eps(nr,r,U,Kt,eps,nut,f1,f2,deta,sigmaK,sigmaE,Ce1,Ce2,tau_mu,tau_R,Pk,Tk,Dk,Peps,Teps,Deps,epseps, &
        detadr,d2etadr2)

    write(fname_ke,111)'momentumNT_var_Re_tau=',int(Re_tau),'_log10(Pr)=',int(log10(Pr)),'.csv'
    111 format(a22,i3,a11,i2,a4)

    open(14,file=fname_ke,form='formatted')
    write(14,*) '"r","U","k","epsilon","nut","tau_mu","tau_R","Pk","Tk","Dk","Peps","Teps","Deps","epsEps"'
    do j=1,nr
       write(14,101) r(j),',',U(j),',',Kt(j),',',eps(j),',',nut(j),',',tau_mu(j),',',tau_R(j),',',Pk(j),',',Tk(j),',',Dk(j),',', &
       Peps(j),',',Teps(j),',',Deps(j),',',epsEps(j)
    enddo
    close(14)

    call output_fields_thermal(nr,r,T,Th2,lambda,dlambdadT,d2lambdadT2,Kt,eps,nut,detadr,d2etadr2,deta,sigmaT,sigmaTh2, &
        Pr, q_lam,q_R,q_new,P_Th2,eps_Th2,T_Th2,D_Th2,H_Th2)

    write(fname_th,112)'thermalNT_var_Re_tau=',int(Re_tau),'_log10(Pr)=',int(log10(Pr)),'.csv'
    112 format(a21,i3,a11,i2,a4)

    open(15,file=fname_th,form='formatted')
    write(15,*) '"r","T","Th2","q_lam","q_R","q_new","P_Th2","eps_Th2","T_Th2","D_Th2","H_Th2","lambda"'
    do j=1,nr
       write(15,103) r(j),',',T(j),',',Th2(j),',',q_lam(j),',',q_R(j),',',q_new(j),',',P_Th2(j),',',eps_Th2(j),',',T_Th2(j), & 
       ',',D_Th2(j),',',H_Th2(j),',',lambda(j)
    enddo
    close(15)

    101 format(e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10, &
    A,e18.10)

    102 format(e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10)

    103 format(e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10)

    end

!!!*************************************************
!!!*						         	             *
!!!*            initialization                        *
!!!*								                 *
!!!*************************************************

subroutine initialization(flag,nr,Re_tau,Pr,Bp,Cp,dr_min,r,detadr,d2etadr2,U,Kt,eps,T,Th2,deta)
    implicit none
    logical, intent(in) :: flag
    integer, intent(in) :: nr
    real*8, intent(in) :: Re_tau,dr_min,Pr,Bp,Cp
    real*8, intent(out) :: r(1:nr),detadr(1:nr),d2etadr2(1:nr),U(1:nr),kt(1:nr),T(1:nr),Th2(1:nr),eps(1:nr),deta
    integer j
    real*8 Kappa, Cmu,Ce1,Ce2,nut(1:nr),f2(1:nr),r_mid,dUdeta(1:nr),dUdr(1:nr), r_wall, a_co,b_co

    Kappa = 4.d-1
    Cmu = 9.d-2
    Ce1=1.45d0
    Ce2=1.9d0
    r_mid = 11.635

    !initial conditions
    if (flag) then
        
        print *,'continuation'
  
        open(511,file='point_ke_var.dat',form='unformatted')
        read(511) r(1:nr),detadr(1:nr),d2etadr2(1:nr),U(1:nr),Kt(1:nr),eps(1:nr),nut(1:nr),f2(1:nr)
        close(511)        
        
        open(511,file='point_T_var.dat',form='unformatted')
        read(511) r(1:nr),T(1:nr),Th2(1:nr)
        close(511)

        deta = Re_tau/(nr-1)

        !T = T - T(1)
  
    else
  
        print *,'new run'

        call grid(nr,dr_min,Re_tau,r,detadr,d2etadr2,deta)

        a_co = 3.d0 * ( 3.d0 - Re_tau / (Re_tau - r_mid) ) - r_mid / Kappa
        b_co = 1.d0 / Kappa - (2.d0 / r_mid) * ( 2.d0 - Re_tau / (Re_tau - r_mid) )

        a_co = ((Re_tau - r_mid)**2.d0) / ( Re_tau * (4.d0 * Re_tau - 5.d0 * r_mid) * r_mid**2.d0) * a_co
        b_co = ((Re_tau - r_mid)**2.d0) / ( Re_tau * (4.d0 * Re_tau - 5.d0 * r_mid) * r_mid**2.d0) * b_co

        do j=nr,1,-1
            r_wall = r(nr) - r(j)
            if(r_wall<=r_mid) then
                U(j) = r_wall
                !eps(j) = (r_wall*Kappa -r_wall*r_wall*Kappa/Re_tau)/(Kappa*r_mid)**2.d0 ! 0.1335d0 / ( 1.d0 + ( ( r(j) - 15.515d0 )**2.d0 ) / 166.7634d0 ) !
                eps(j) = 2.d0*( (Re_tau - r_mid) / (Re_tau * r_mid**2.d0) ) + &
                    r_wall * ( (Re_tau - r_mid) / (Re_tau * r_mid**2.d0) ) * (1.d0 / Kappa - 2.d0 / r_mid) !2.d0 * a_co * ( 2.d0 - Re_tau / r(j) ) + 3.d0 * b_co * r_wall  * ( 3.d0 - Re_tau / r(j) )  !4.d-04 * dexp( -( (r_wall - r_mid)/(dsqrt(8.d0*Re_tau)) )**2.d0 )
                Kt(j) = ( (Re_tau - r_mid) / (Re_tau * r_mid**2.d0) ) * r_wall**2.d0 !+ b_co * r_wall**3.d0 ! ((r_wall/r_mid)**2.d0) * dsqrt( (Re_tau-r_mid)*r_mid*(Kappa/(Cmu*Re_tau))*eps(j) )!dsqrt( r_wall*r_wall*Kappa*eps(j)/Cmu ) ! 0.019678d0 * r(j) * r(j) / ( 1.d0 + ( ( r(j) - 7.28d0 )**2.d0 ) / 88.263d0 ) !
            else
                U(j) = (1.d0/Kappa)*dlog(r_wall) +5.5d0 
                eps(j) = ( r(j) / r_wall )/(Kappa*Re_tau)
                Kt(j) = r(j) / Re_tau  ! 0.019678d0 * r(j) * r(j) / ( 1.d0 + ( ( r(j) - 7.28d0 )**2.d0 ) / 88.263d0 ) !
            endif
            !eps(j) = 4.d-04
            !kt(j) = 1.d-04 * r_wall**2.d0
            Kt(j) = Kt(j)/dsqrt(Cmu)
        enddo
        
        call ddeta(nr,U,dUdeta,deta)
        dUdr = dUdeta*detadr

        do j=1,nr
            Th2(j) = Kt(j)
        enddo

        !!!********************************************************
        !!!
        !!! dr^+ = -K^+(T)dT^+
        !!!
        !!!********************************************************
  
        do j=1,nr
            T(j) = - Pr*( U(j) - U(1) )
        enddo

        !T = T - T(1)

    endif

    end

!!!*************************************************
!!!*						         	             *
!!!*                    ddeta                        *
!!!*								                 *
!!!*************************************************
    
subroutine  ddeta(nr,A,DA,deta)
    implicit none
    integer, intent(in) :: nr
    real*8, intent(in) :: A(1:nr),deta
    real*8, intent(out) :: DA(1:nr)
    integer j

    DA(1) =  -(3.d0*A(1) -4.d0*A(2) +A(3))/(2.d0*deta)

    do j=2,nr-1
        DA(j) = (A(j+1) -A(j-1))/(2.d0*deta)
    enddo

    DA(nr) = (3.d0*A(nr) -4.d0*A(nr-1) +A(nr-2))/(2.d0*deta)

    end

!!!*************************************************
!!!*						         	             *
!!!*                    d2deta2                        *
!!!*								                 *
!!!*************************************************

subroutine  d2deta2(nr,A,D2A,deta)
    implicit none
    integer, intent(in) :: nr
    real*8, intent(in) :: A(1:nr),deta
    real*8, intent(out) :: D2A(1:nr)
    real*8 deta2
    integer j

    deta2 = deta*deta
    
    D2A(1) =  (12.d0*a(1) -30.d0*a(2) +24.d0*a(3) -6.d0*a(4))/(6.d0*deta2)
    
    do j=2,nr-1
        D2A(j) = (a(j+1) - 2.d0*a(j) + a(j-1))/deta2
    enddo
    
    D2A(nr) = (12.d0*a(nr) -30.d0*a(nr-1) +24.d0*a(nr-2) -6.d0*a(nr-3))/(6.d0*deta2)
    
    end

!!!***************************************************
!!!*						         	               *
!!!*       Nagano Takawa K - Epsilon Constants 	       	   *
!!!*								                   *
!!!***************************************************
    
subroutine  nagano_takawa_k_epsilon_constants(sigmak,sigmae,Ce1,Ce2,Cmu,f1)
    implicit none
    real*8, intent(out) :: sigmak,sigmae,Ce1,Ce2,Cmu,f1

    sigmaK= 1.0d0
    sigmae= 1.3d0
    Ce1=1.45d0
    Ce2=1.9d0
    Cmu=0.09d0
    f1=1.d0
    
    end

!!!***************************************************
!!!*						         	               *
!!!*       Nagano Takawa K - Epsilon Functions 	       	   *
!!!*								                   *
!!!***************************************************
    
subroutine  nagano_takawa_k_epsilon_functions(nut,f2,nr,r,kt,eps,Cmu)
    implicit none
    integer, intent(in) :: nr
    real*8, intent(in) :: r(1:nr),Kt(1:nr),eps(1:nr),Cmu
    real*8, intent(out) :: nut(1:nr),f2(1:nr)
    real*8 Ret(1:nr), fmu(1:nr), Ret_min, eps_min, r_wall
    integer j

    eps_min = 1.d-60

    do j=1,nr
        if (eps(j) <= eps_min) then
            Ret(j)= dabs(Kt(j)*Kt(j)/eps_min)
        else
            Ret(j)= dabs(Kt(j)*Kt(j)/eps(j))
        endif
    enddo

    Ret_min = 1.d-60
    do j=1,nr
        r_wall = r(nr) - r(j)
        if (Ret(j) <= Ret_min) then
            fmu(j)= (1.d0 +4.1d0/Ret_min**0.75d0)*(1.d0 -dexp(-r_wall/26.d0))**2.d0
        else
            fmu(j)= (1.d0 +4.1d0/Ret(j)**0.75d0)*(1.d0 -dexp(-r_wall/26.d0))**2.d0
        endif
    enddo

    do j=1,nr
        nuT(j)= Cmu*fmu(j)*Ret(j)
    enddo

    do j=1,nr
        r_wall = r(nr) - r(j)
        f2(j)= (1.d0 -0.3d0*dexp(-(Ret(j)/6.5d0)**2.d0))*(1.d0 -dexp(-r_wall/6.d0))**2.d0
    enddo
    
    end

!!!***************************************************
!!!*						         	               *
!!!*                U coefficients	       	   *
!!!*								                   *
!!!***************************************************
subroutine  u_coefficients(aU_w,aU_e,sU,nut,dnutdr,deta,Re_tau,d2etadr2,detadr,r)
    implicit none
    real*8, intent(in) :: nut,dnutdr,deta,Re_tau,d2etadr2,detadr,r
    real*8, intent(out) :: aU_w,aU_e,sU
    real*8 dev

    if(r>0.d0) then
        dev = deta*( d2etadr2 + (1.d0/r + dnutdr/(1.d0+nut))*detadr )/(4.d0*detadr**2.d0)
        sU = (deta*deta)/(Re_tau*(1.d0+nut)*(detadr)**2.d0)
    else 
        dev = (deta / 4.d0) * d2etadr2 / (detadr)**2.d0
        sU = (deta*deta)/(2.d0*Re_tau*(1.d0+nut)*(detadr)**2.d0)
    endif

    aU_w = 5.d-1 - dev
    aU_e = 5.d-1 + dev
    
    end

!!!***************************************************
!!!*						         	               *
!!!*                K coefficients	       	   *
!!!*								                   *
!!!***************************************************
subroutine  K_coefficients(aK_w,aK_e,sK,eps,nut,dnutdr,dUdr,deta,sigmak,d2etadr2,detadr,r)
    implicit none
    real*8, intent(in) :: eps,nut,dnutdr,dUdr,deta,sigmak,d2etadr2,detadr,r
    real*8, intent(out) :: aK_w,aK_e,sK
    real*8 dev

    if(r>0.d0) then
        dev = deta*( d2etadr2 + (1.d0/r + dnutdr/(sigmak+nut))*detadr )/(4.d0*detadr**2.d0)
        sK = (nut*dUdr*dUdr - eps)*(deta*deta)/(2.d0*(1.d0+nut/sigmak)*detadr**2.d0)
    else
        dev = (deta / 4.d0) * d2etadr2 / (detadr)**2.d0
        sK = - eps*(deta*deta)/(4.d0*(1.d0+nut/sigmak)*detadr**2.d0)
    endif

    aK_w = 5.d-1 - dev
    aK_e = 5.d-1 + dev
    

    end

!!!***************************************************
!!!*						         	               *
!!!*                E coefficients	       	   *
!!!*								                   *
!!!***************************************************
subroutine  E_coefficients(aE_w,aE_e,sE,eps,Kt,nut,dnutdr,dUdr,deta,sigmae,Ce1,f1,Ce2,f2,d2etadr2,detadr,r)
    implicit none
    real*8, intent(in) :: eps,Kt,nut,dnutdr,dUdr,deta,sigmae,Ce1,f1,Ce2,f2,d2etadr2,detadr,r
    real*8, intent(out) :: aE_w,aE_e,sE
    real*8 K_min, Kb, dev
    logical method1

    K_min = 1.d-60

    if(r>0.d0) then
        dev = deta*( d2etadr2 + (1.d0/r + dnutdr/(sigmae+nut))*detadr )/(4.d0*detadr**2.d0)
        Kb = (Ce1*f1*nut*dUdr*dUdr -Ce2*f2*eps)*(deta*deta/(2.d0*(1.d0+nut/sigmae)*detadr**2.d0))
    else
        dev = (deta / 4.d0) * d2etadr2 / (detadr)**2.d0
        Kb = -Ce2*f2*eps*(deta*deta/(4.d0*(1.d0+nut/sigmae)*detadr**2.d0))
    endif

    method1 = .false.
    if (method1) then
        aE_w = 5.d-1 - dev
        aE_e = 5.d-1 + dev
        if (Kt<=K_min) then
            sE = Kb*eps/K_min
        else
            sE = Kb*eps/Kt
        endif
    else
        aE_w = (5.d-1 - dev)/(1.d0 - Kb/Kt)
        aE_e = (5.d-1 + dev)/(1.d0 - Kb/Kt)
        sE = 0.d0 
    endif
    
    end

!!!***************************************************
!!!*						         	             *
!!!*                T coefficients	       	         *
!!!*								                 *
!!!***************************************************

subroutine  T_coefficients(aT_w,aT_e,sT,r,Re_tau,nut,dnutdr,lambda,dlambdadT,d2lambdadT2,dTh2dr,d2Th2dr2,dTdr, &
    Pr,sigmaT,deta,d2etadr2,detadr)
    implicit none
    real*8, intent(in) :: r,Re_tau
    real*8, intent(in) :: nut,dnutdr,lambda,dlambdadT,d2lambdadT2,dTh2dr,d2Th2dr2,dTdr,Pr,sigmaT,deta,d2etadr2,detadr
    real*8, intent(out) :: aT_w,aT_e,sT
    real*8 A1,A2,A3, dev

    if(r>0.d0) then
        A1 = lambda / Pr + nut / sigmaT
        A2 = A1 / r + ( dlambdadT * dTdr + d2lambdadT2 * dTh2dr ) / Pr + dnutdr / sigmaT
        A3 =  2.d0 / Re_tau + ( dlambdadT / Pr ) * ( d2Th2dr2 + dTh2dr / r )
    else
        A1 = 2.d0 * ( lambda / Pr + nut / sigmaT )
        A2 = 0.d0
        A3 =  2.d0 / Re_tau + 2.d0 * ( dlambdadT / Pr ) * d2Th2dr2
    endif

    !print*, ' A1 = ',A1,' A2 = ',A2,' A3 = ',A3

    dev = deta * ( d2etadr2 / detadr + A2 / A1 ) / ( 4.d0 * detadr )

    aT_w = ( 5.d-1 - dev )
    aT_e = ( 5.d-1 + dev )
    sT = ( ( deta * deta ) / ( 2.d0 * detadr**2.d0 ) ) * A3 / A1

    end

!!!***************************************************
!!!*						         	             *
!!!*                Th2 coefficients	       	     *
!!!*								                 *
!!!***************************************************

subroutine  Th2_coefficients(aTh2_w,aTh2_e,sTh2,r,nut,dnutdr,lambda,dlambdadT,d2lambdadT2,dTdr,d2Tdr2, &
    Pr,sigmaTh2,deta,d2etadr2,detadr,eps,Kt)
    implicit none
    real*8, intent(in) :: nut,dnutdr,lambda,dlambdadT,d2lambdadT2,dTdr,d2Tdr2,Pr,sigmaTh2,deta,d2etadr2,detadr
    real*8, intent(in) :: r,eps,Kt
    real*8, intent(out) :: aTh2_w,aTh2_e,sTh2
    real*8 A1,A2,A3,A4,dev, den

    if(r>0.d0) then
        A1 = lambda / Pr + nut / sigmaTh2
        A2 = 2.d0 * dlambdadT * dTdr / Pr + dnutdr / sigmaTh2 + A1 / r
        A3 =  ( nut / sigmaTh2 ) * dTdr**2.d0
        A4 = ( ( eps / Kt ) * lambda - 2.d0 * dlambdadT * d2Tdr2 - 2.d0 *d2lambdadT2 * dTdr**2.d0 &
            - 2.d0 *dlambdadT * dTdr / r ) / Pr
    else 
        A1 = 2.d0 * ( lambda / Pr + nut / sigmaTh2 )
        A2 = 0.d0
        A3 =  0.d0
        A4 = ( ( eps / Kt ) * lambda - 4.d0 * dlambdadT * d2Tdr2 ) / Pr
    endif

    dev = deta * ( d2etadr2 / detadr + A2 / A1 ) / ( 4.d0 * detadr )

    den = ( 2.d0 * detadr**2.d0 + deta * deta * ( A4 / A1 ) )

    aTh2_w = ( 5.d-1  - dev ) * ( 2.d0 * detadr**2.d0 ) / den
    aTh2_e = ( 5.d-1  + dev ) * ( 2.d0 * detadr**2.d0 ) / den
    sTh2 = deta * deta * ( A3 / A1 ) / den

    end

!!!*************************************************
!!!*						         	           *
!!!*          output fields k-eps                  *
!!!*								               *
!!!*************************************************

subroutine  output_fields_k_eps(nr,r,U,Kt,eps,nut,f1,f2,deta,sigmaK,sigmaE,Ce1,Ce2,tau_mu,tau_R,Pk,Tk,Dk,Peps,Teps,Deps,epseps, &
    detadr, d2etadr2)
    implicit none
    integer, intent(in) :: nr
    real*8, intent(in) :: deta,sigmaK,sigmaE,Ce1,Ce2,f1,f2(1:nr),r(1:nr),Kt(1:nr),eps(1:nr),nut(1:nr),U(1:nr)
    real*8, intent(in) :: detadr(1:nr),d2etadr2(1:nr)
    real*8, INTENT(OUT) :: tau_mu(1:nr),tau_R(1:nr),Pk(1:nr),Tk(1:nr),Dk(1:nr)
    real*8, INTENT(OUT) :: Peps(1:nr),Teps(1:nr),Deps(1:nr),epseps(1:nr)
    real*8 dUdr(1:nr),d2Ktdeta2(1:nr),d2epsdeta2(1:nr),dKtdeta(1:nr),depsdeta(1:nr),dnutdr(1:nr),dUdeta(1:nr),dnutdeta(1:nr)

    call ddeta(nr,U,dUdeta,deta)
    dUdr = dUdeta*detadr
    call ddeta(nr,nut,dnutdeta,deta)
    dnutdr = dnutdeta*detadr

    tau_mu = dUdr
    tau_R = nut*dUdr
    Pk = tau_R*dUdr

    call d2deta2(nr,Kt,d2Ktdeta2,deta)
    call ddeta(nr,Kt,dKtdeta,deta)

    Dk = d2Ktdeta2*detadr**2.d0 + dKtdeta*d2etadr2
    Dk(1) = Dk(1) + d2Ktdeta2(1)*detadr(1)**2.d0 + dKtdeta(1)*d2etadr2(1)
    Dk(2:nr) = Dk(2:nr) + dKtdeta(2:nr) * detadr(2:nr) / r(2:nr)
    Tk = (nut/sigmaK)*Dk + (dKtdeta*detadr/sigmaK)*dnutdr  

    Peps(1:nr-1) = f1*Ce1*(eps(1:nr-1)/Kt(1:nr-1))*Pk(1:nr-1)
    Peps(nr) = Peps(nr-1)

    call d2deta2(nr,eps,D2epsdeta2,deta)
    call ddeta(nr,eps,depsdeta,deta)

    Deps = d2epsdeta2*detadr**2.d0 + depsdeta*d2etadr2
    Deps(1) = Deps(1) + d2epsdeta2(1)*detadr(1)**2.d0 + depsdeta(1)*d2etadr2(1)
    Deps(2:nr) = Deps(2:nr) + depsdeta(2:nr) * detadr(2:nr) / r(2:nr)
    Teps = (nut/sigmaE)*Deps + (depsdeta*detadr/sigmaE)*dnutdr 
    epsEps(1:nr-1) = -f2(1:nr-1)*Ce2*(eps(1:nr-1)/Kt(1:nr-1))*eps(1:nr-1)
    epsEps(nr) = epsEps(nr-1)
    
    end

!!!*************************************************
!!!*						         	           *
!!!*          output fields thermal                  *
!!!*								               *
!!!*************************************************

subroutine  output_fields_thermal(nr,r,T,Th2,lambda,dlambdadT,d2lambdadT2,Kt,eps,nut,detadr,d2etadr2,deta,sigmaT,sigmaTh2, &
    Pr, q_lam,q_R,q_new,P_Th2,eps_Th2,T_Th2,D_Th2,H_Th2)
    implicit none
    integer, intent(in) :: nr
    real*8, intent(in) :: r(1:nr),deta,Kt(1:nr),eps(1:nr),nut(1:nr),detadr(1:nr),d2etadr2(1:nr)
    real*8, intent(in) :: T(1:nr),Th2(1:nr),lambda(1:nr),dlambdadT(1:nr),d2lambdadT2(1:nr),sigmaT,sigmaTh2,Pr
    real*8, INTENT(OUT) :: q_lam(1:nr),q_R(1:nr),q_new(1:nr),P_Th2(1:nr),eps_Th2(1:nr),T_Th2(1:nr),D_Th2(1:nr),H_Th2(1:nr)
    real*8 dTdr(1:nr),dTdeta(1:nr),d2Tdr2(1:nr),d2Tdeta2(1:nr),dTh2dr(1:nr),dTh2deta(1:nr),d2Th2dr2(1:nr),d2Th2deta2(1:nr)
    real*8 dnutdr(1:nr),dnutdeta(1:nr)

    call ddeta(nr,T,dTdeta,deta)
    dTdr = dTdeta*detadr
    call ddeta(nr,Th2,dTh2deta,deta)
    dTh2dr = dTh2deta*detadr
    call ddeta(nr,nut,dnutdeta,deta)
    dnutdr = dnutdeta*detadr
    call d2deta2(nr,T,d2Tdeta2,deta)
    d2Tdr2 = d2Tdeta2*detadr**2.d0 + dTdeta*d2etadr2
    call d2deta2(nr,Th2,d2Th2deta2,deta)
    d2Th2dr2 = d2Th2deta2*detadr**2.d0 + dTh2deta*d2etadr2

    q_lam = ( lambda / Pr ) * dTdr
    q_R = ( nut / sigmaT ) * dTdr
    q_new = ( dlambdadT / Pr ) * dTh2dr

    P_Th2 = ( nut / sigmaT ) * dTdr**2.d0
    eps_Th2(1:nr-1) = ( eps(1:nr-1) / Kt(1:nr-1) ) * ( lambda(1:nr-1) / Pr ) * Th2(1:nr-1)
    eps_Th2(nr) = eps_Th2(nr-1)

    T_Th2 = ( nut / sigmaTh2 )*d2Th2dr2 + ( dTh2dr / sigmaTh2 ) * dnutdr
    T_Th2(1) = T_Th2(1) + ( nut(1) / sigmaTh2 )*d2Th2dr2(1) + ( dTh2dr(1) / sigmaTh2 ) * dnutdr(1)
    T_Th2(2:nr) = T_Th2(2:nr) + ( nut / sigmaTh2 ) * dTh2dr(2:nr) / r(2:nr)
    D_Th2 = ( lambda / Pr )*d2Th2dr2 + ( dTh2dr / Pr ) * dlambdadT * dTdr
    D_Th2(1) = D_Th2(1) + ( lambda(1) / Pr )*d2Th2dr2(1) + ( dTh2dr(1) / Pr ) * dlambdadT(1) * dTdr(1)
    D_Th2(2:nr) = D_Th2(2:nr) + ( lambda(2:nr) / Pr ) * dTh2dr(2:nr) / r(2:nr)

    H_Th2 = ( 2.d0 / Pr ) * d2lambdadT2 * Th2 * dTdr**2.d0 + ( 2.d0 / Pr ) * dlambdadT * Th2 * d2Tdr2 &
    + ( 1.d0 / Pr ) * dlambdadT * dTh2dr * dTdr
    H_Th2(1) = H_Th2(1) + ( 2.d0 / Pr ) * dlambdadT(1) * Th2(1) * d2Tdr2(1)
    H_Th2(2:nr) =  H_Th2(2:nr) + ( 2.d0 / Pr ) * dlambdadT(2:nr) * Th2(2:nr) * dTdr(2:nr) / r(2:nr)
    
    end


!!!*************************************************
!!!*						         	           *
!!!*               residuals                   *
!!!*								               *
!!!*************************************************

subroutine  residuals(nr,A,A0,resA)
    implicit none
    integer, intent(in) :: nr
    real*8, intent(in) :: A(1:nr),A0(1:nr)
    real*8, INTENT(OUT) :: resA
    real*8 sumN, sumD
    integer j

    sumN = 0
    do j=1,nr
        sumN = sumN +dabs(A(j)- A0(j))
    enddo

    sumD = 0
    do j=1,nr
        sumD = sumD + dabs(A0(j))
    enddo

    resA = sumN/sumD

    end

!!!*************************************************
!!!*						         	             *
!!!*            grid                       *
!!!*								                 *
!!!*************************************************

subroutine grid(nr,dr_min,Re_tau,r,detadr,d2etadr2,deta)
    implicit none
    integer, intent(in) :: nr
    real*8, intent(in) :: Re_tau,dr_min
    real*8, intent(out) :: r(1:nr),detadr(1:nr),d2etadr2(1:nr),deta
    integer j
    real*8 a, b, c, d, e, eta, r_max

    r_max = Re_tau

    deta = r_max/(nr-1)

    a = (dr_min - deta)/(Re_tau*deta - deta*deta)
    b = (2.d0 * Re_tau * deta - dr_min*Re_tau - deta*deta)/(Re_tau*deta - deta*deta)

    do j=1,nr
        eta = deta*(j-1)
        r(j) = a*eta**2.d0 + b*eta
        detadr(j) = 1.d0/(2.d0*a*eta + b)
        d2etadr2(j) = -2.d0*a/(2.d0*a*eta + b)**3.d0
    enddo

    print*, ' dr_min =', r(nr)-r(nr-1), ' dr_max =', r(2)-r(1), ' ratio =', (r(2)-r(1))/(r(nr)-r(nr-1))

    end

!!!*************************************************
!!!*						         	           *
!!!*            Thermal Diffusivity                *
!!!*								               *
!!!*************************************************
    
subroutine  thernal_diffusivity(lambda,dlambdadT,d2lambdadT2,T,Bp,Cp,Dp,Ep)
    implicit none
    real*8, intent(in) :: T,Bp,Cp,Dp,Ep
    real*8, intent(out) :: lambda,dlambdadT,d2lambdadT2
    integer j

    lambda = 1.d0 + Bp * T + Cp * T**2.d0 + Dp * T**3.d0 + Ep * T**4.d0
    dlambdadT = Bp + 2.d0 * Cp * T + 3.d0 * Dp * T**2.d0 + 4.d0 * Ep * T**3.d0
    d2lambdadT2 = 2.d0 * Cp + 6.d0 * Dp * T + 12.d0 * Ep * T**2.d0

    end

!!!*************************************************
!!!*						         	           *
!!!*            Thomas Algorithm                *
!!!*								               *
!!!*************************************************
    
subroutine  thomas_algorithm(var,nr,a_e,a_w,a_p,S)
    implicit none
    integer, intent(in) :: nr
    real*8, intent(inout) :: var(1:nr),a_e,a_w,a_p,S
    real*8 A(1:nr),C_apex(1:nr)
    integer j

    A(1) = 0.d0
    A(nr) = 0.d0
    C_apex(1) = var(1)
    C_apex(nr) = var(nr)

    do j=2,nr-1
        A(j) = a_e / ( a_p - a_w * A(j-1) )
        C_apex(j) = ( a_w * C_apex(j-1) + S ) / ( a_p - a_w * A(j-1) )
        var(j) = A(j) * var(j+1) + C_apex(j)
    enddo

    end

!!!*************************************************
!!!*						         	           *
!!!*            Correct Residuals               *
!!!*								               *
!!!*************************************************
    
subroutine  correct_residuals(alpha,res,res_old)
    implicit none
    real*8, intent(in) :: res,res_old
    real*8, intent(inout) :: alpha
    real*8 increment

    increment = 1.d-06

    if(res < res_old) then
        alpha = alpha * ( 1.d0 + increment )
    elseif(res > res_old) then
        alpha = alpha * ( 1.d0 + increment )
    else
        alpha = alpha
    endif

    alpha = min(1.d0,max(alpha,increment))

    end

!!!*************************************************
!!!*						         	           *
!!!*                 Solve U                        *
!!!*								               *
!!!*************************************************
    
subroutine  solve_u(U,nut,dnutdr,r,detadr,d2etadr2,deta,Re_tau,nr)
    implicit none
    integer, intent(in) :: nr
    real*8, intent(in) :: nut(1:nr),dnutdr(1:nr),detadr(1:nr),d2etadr2(1:nr),r(1:nr),deta,Re_tau
    real*8, intent(inout) :: U(1:nr)
    real*8 aU_w,aU_e,sU
    real*8 A(1:nr),C_apex(1:nr),denominator
    integer j

    call u_coefficients(aU_w,aU_e,sU,nut(1),dnutdr(1),deta,Re_tau,d2etadr2(1),detadr(1),r(1))
                       
    A(1) = aU_e + aU_w
    C_apex(1) = sU
    do j =2,nr
        call u_coefficients(aU_w,aU_e,sU,nut(j),dnutdr(j),deta,Re_tau,d2etadr2(j),detadr(j),r(j))
        denominator = ( 1.d0 - aU_w * A(j-1) )
        A(j) = aU_e / denominator
        C_apex(j) = ( aU_w * C_apex(j-1) + sU ) / denominator
    enddo
    U(nr) = 0.d0

    !U(1) = C_apex(1) + A(1)*U(2)
    !do j =2,nr-1
    !    call u_coefficients(aU_w,aU_e,sU,nut(j),dnutdr(j),deta,Re_tau,d2etadr2(j),detadr(j),r(j))
    !    U(j) =  sU + aU_e*U(j+1) + aU_w*U(j-1)
    !    !!!U(j) = A(j) * U(j+1) + C_apex(j)
    !enddo
    U(nr) = 0.d0

    do j =nr-1,1,-1 ! 1,nr-1 !
        U(j) = A(j) * U(j+1) + C_apex(j)
    enddo

    end

!!!*************************************************
!!!*						         	           *
!!!*                 Solve Kt                        *
!!!*								               *
!!!*************************************************
    
subroutine  solve_Kt(Kt,eps,dUdr,nut,dnutdr,r,detadr,d2etadr2,deta,sigmak,nr)
    implicit none
    integer, intent(in) :: nr
    real*8, intent(in) :: eps(1:nr),dUdr(1:nr),nut(1:nr),dnutdr(1:nr),detadr(1:nr),d2etadr2(1:nr),r(1:nr),deta,sigmak
    real*8, intent(inout) :: Kt(1:nr)
    real*8 aK_w,aK_e,sK, temp
    real*8 A(1:nr),C_apex(1:nr),denominator
    integer j
    
    call K_coefficients(aK_w,aK_e,sK,eps(1),nut(1),dnutdr(1),dUdr(1),deta,sigmak,d2etadr2(1),detadr(1),r(1))
    A(1) = aK_e + aK_w
    C_apex(1) = sK
    do j =2,nr
        call K_coefficients(aK_w,aK_e,sK,eps(j),nut(j),dnutdr(j),dUdr(j),deta,sigmak,d2etadr2(j),detadr(j),r(j))
        denominator = ( 1.d0 - aK_w * A(j-1) )
        A(j) = aK_e / denominator
        C_apex(j) = ( aK_w * C_apex(j-1) + sK ) / denominator
    enddo
    Kt(nr) = 0.d0
    Kt(1) = C_apex(1) + A(1)*Kt(2)
    do j =2,nr-1
        call K_coefficients(aK_w,aK_e,sK,eps(j),nut(j),dnutdr(j),dUdr(j),deta,sigmak,d2etadr2(j),detadr(j),r(j))
        temp = sK + aK_e*Kt(j+1) + aK_w*Kt(j-1)
        if(temp < 0.d0) then 
            print*, ' j =',j,'; sK = ',sK,'; aK_e = ',aK_e,'; aK_w = ',aK_w
        else
            Kt(j) =  temp
        endif
    enddo
    
    !do j =nr-1,1,-1 ! 1,nr-1 !
    !    temp = A(j) * Kt(j+1) + C_apex(j)
    !    if(temp < 0.d0) then 
    !        print*, ' j =',j,'; A(j) = ',A(j),'; C_apex(j) = ',C_apex(j)
    !    else
    !        Kt(j) =  temp
    !    endif
    !enddo

    end

!!!*************************************************
!!!*						         	           *
!!!*                 Solve Eps                        *
!!!*								               *
!!!*************************************************
    
subroutine  solve_eps(Kt,eps,dUdr,nut,dnutdr,r,detadr,d2etadr2,deta,sigmae,ce1,ce2,f1,f2,nr)
    implicit none
    integer, intent(in) :: nr
    real*8, intent(in) :: Kt(1:nr),dUdr(1:nr),nut(1:nr),dnutdr(1:nr),detadr(1:nr),d2etadr2(1:nr),r(1:nr),deta,sigmae
    real*8, intent(in) :: ce1,ce2,f1,f2(1:nr)
    real*8, intent(inout) :: eps(1:nr)
    real*8 aE_w,aE_e,sE
    real*8 A(1:nr),C_apex(1:nr),denominator,temp
    integer j

    call E_coefficients(aE_w,aE_e,sE,eps(1),Kt(1),nut(1),dnutdr(1),dUdr(1),deta,sigmae,Ce1,f1,Ce2,f2(1),d2etadr2(1), &
        detadr(1),r(1))
    A(1) = aE_e + aE_w
    C_apex(1) = sE
    do j =2,nr
        call E_coefficients(aE_w,aE_e,sE,eps(j),Kt(j),nut(j),dnutdr(j),dUdr(j),deta,sigmae,Ce1,f1,Ce2,f2(j),d2etadr2(j), &
        detadr(j),r(j))
        denominator = ( 1.d0 - aE_w * A(j-1) )
        A(j) = aE_e / denominator
        C_apex(j) = ( aE_w * C_apex(j-1) + sE ) / denominator
    enddo

    !print*, ' Kt(nr) = ',Kt(nr), '; Kt(nr-1) = ',Kt(nr-1), '; Kt(nr-2) = ',Kt(nr-2)
    eps(nr) = 2.d0*( ( (-3.d0*dsqrt(dabs(Kt(nr)))+4.d0*dsqrt(dabs(Kt(nr-1)))-dsqrt(dabs(Kt(nr-2))))/(2.d0*deta) ) &
        *detadr(nr) )**2.d0
    eps(1) = C_apex(1) + A(1)*eps(2)
    do j =2,nr-1
        call E_coefficients(aE_w,aE_e,sE,eps(j),Kt(j),nut(j),dnutdr(j),dUdr(j),deta,sigmae,Ce1,f1,Ce2,f2(j),d2etadr2(j), &
        detadr(j),r(j))
        eps(j) = sE + aE_e*eps(j+1) + aE_w*eps(j-1)
    enddo
    !print*, ' eps(nr) = ',eps(nr)
    !do j =nr-1,1,-1 ! 1,nr-1 !
    !    temp = A(j) * eps(j+1) + C_apex(j)
    !    if(temp < 0.d0) then 
    !        print*, ' j =',j,'; A(j) = ',A(j),'; C_apex(j) = ',C_apex(j)
    !    else
    !        eps(j) =  temp
    !    endif
    !enddo    

    end

!!!*************************************************
!!!*						         	           *
!!!*                 Solve T                       *
!!!*								               *
!!!*************************************************
    
subroutine  solve_T(T,nut,dnutdr,lambda,dlambdadT,d2lambdadT2,dTh2dr,d2Th2dr2,dTdr,r,Re_tau,Pr,sigmaT,deta,d2etadr2,detadr,nr)
    implicit none
    integer, intent(in) :: nr
    real*8, intent(in) :: r(1:nr),dTdr(1:nr),nut(1:nr),dnutdr(1:nr),detadr(1:nr),d2etadr2(1:nr),deta,sigmaT
    real*8, intent(in) :: Re_tau,Pr,lambda(1:nr),dlambdadT(1:nr),d2lambdadT2(1:nr),dTh2dr(1:nr),d2Th2dr2(1:nr)
    real*8, intent(inout) :: T(1:nr)
    real*8 aT_w,aT_e,sT
    real*8 A(1:nr),C_apex(1:nr),denominator
    integer j

    call T_coefficients(aT_w,aT_e,sT,r(1),Re_tau,nut(1),dnutdr(1),lambda(1),dlambdadT(1),d2lambdadT2(1),dTh2dr(1),d2Th2dr2(1), & 
        dTdr(1),Pr,sigmaT,deta,d2etadr2(1),detadr(1))
    A(1) = aT_e + aT_w
    C_apex(1) = sT
    do j =2,nr
        call T_coefficients(aT_w,aT_e,sT,r(j),Re_tau,nut(j),dnutdr(j),lambda(j),dlambdadT(j),d2lambdadT2(j),dTh2dr(j),d2Th2dr2(j), & 
            dTdr(j),Pr,sigmaT,deta,d2etadr2(j),detadr(j))
        denominator = ( 1.d0 - aT_w * A(j-1) )
        A(j) = aT_e / denominator
        C_apex(j) = ( aT_w * C_apex(j-1) + sT ) / denominator
    enddo
    !call T_coefficients(aT_w,aT_e,sT,r(nr),Re_tau,nut(nr),dnutdr(nr),lambda(nr),dlambdadT(nr),d2lambdadT2(nr),dTh2dr(nr), & 
    !   d2Th2dr2(nr),dTdr(nr), Pr,sigmaT,deta,d2etadr2(nr),detadr(nr))
    !denominator = ( 1.d0 - aT_w * A(nr-1) )
    A(nr) = 0.d0
    !C_apex(nr) = ( aT_w * C_apex(nr-1) + sT ) / denominator

    !do j =2,nr-1
    !    call T_coefficients(aT_w,aT_e,sT,r(j),Re_tau,nut(j),dnutdr(j),lambda(j),dlambdadT(j),d2lambdadT2(j),dTh2dr(j),d2Th2dr2(j), &
    !       dTdr(j),Pr,sigmaT,deta,d2etadr2(j),detadr(j))
    !    T(j) = sT + aT_e*T(j+1) + aT_w*T(j-1)
    !enddo
    !call T_coefficients(aT_w,aT_e,sT,r(nr),Re_tau,nut(nr),dnutdr(nr),lambda(nr),dlambdadT(nr),d2lambdadT2(nr),dTh2dr(nr), &
    !    d2Th2dr2(nr), dTdr(nr), Pr,sigmaT,deta,d2etadr2(nr),detadr(nr))
    !T(nr) = sT + aT_w*T(nr-1) + aT_e*( T(nr-1) + 2.d0 * deta * Pr / ( lambda(nr) * detadr(nr) ) ) 
    T(nr) = 0.d0
    
    do j =nr-1,1,-1
        T(j) = A(j) * T(j+1) + C_apex(j)
    enddo

    !if (mod(nr,2).eq.0) then
    !    T = T - ( T(nr-1) + T(nr) ) / 2.d0
    !else
    !    T = T - T(nr/2)
    !endif
    !T = T - T(1)

    end

!!!*************************************************
!!!*						         	           *
!!!*                 Solve Th2                       *
!!!*								               *
!!!*************************************************
    
subroutine  solve_Th2(Th2,nut,r,dnutdr,lambda,dlambdadT,d2lambdadT2,dTdr,d2Tdr2,Pr,sigmaTh2,deta,d2etadr2,detadr,eps,Kt,nr)
    implicit none
    integer, intent(in) :: nr
    real*8, intent(in) :: r(1:nr),dTdr(1:nr),nut(1:nr),dnutdr(1:nr),detadr(1:nr),d2etadr2(1:nr),deta,sigmaTh2
    real*8, intent(in) :: Pr,lambda(1:nr),dlambdadT(1:nr),d2lambdadT2(1:nr),d2Tdr2(1:nr),Kt(1:nr),eps(1:nr)
    real*8, intent(inout) :: Th2(1:nr)
    real*8 aTh2_w,aTh2_e,sTh2
    real*8 A(1:nr),C_apex(1:nr),denominator
    integer j

    call Th2_coefficients(aTh2_w,aTh2_e,sTh2,r(1),nut(1),dnutdr(1),lambda(1),dlambdadT(1),d2lambdadT2(1),dTdr(1),d2Tdr2(1), &
        Pr,sigmaTh2,deta,d2etadr2(1),detadr(1),eps(1),Kt(1))
   
    A(1) = aTh2_e + aTh2_w
    C_apex(1) = sTh2
    do j =2,nr
        call Th2_coefficients(aTh2_w,aTh2_e,sTh2,r(j),nut(j),dnutdr(j),lambda(j),dlambdadT(j),d2lambdadT2(j),dTdr(j),d2Tdr2(j), &
        Pr,sigmaTh2,deta,d2etadr2(j),detadr(j),eps(j),Kt(j))
        denominator = ( 1.d0 - aTh2_w * A(j-1) )
        A(j) = aTh2_e / denominator
        C_apex(j) = ( aTh2_w * C_apex(j-1) + sTh2 ) / denominator
    enddo
    !call Th2_coefficients(aTh2_w,aTh2_e,sTh2,r(nr),nut(nr),dnutdr(nr),lambda(nr),dlambdadT(nr),d2lambdadT2(nr),dTdr(nr),d2Tdr2(nr), &
    !    Pr,sigmaTh2,deta,d2etadr2(nr),detadr(nr),eps(nr),Kt(nr))
    !denominator = ( 1.d0 - aTh2_w * A(nr-1) )
    A(nr-1) = 0.d0
    !C_apex(nr) = ( aTh2_w * C_apex(nr-1) + sTh2 ) / denominator
    Th2(nr) = 0.d0
    
    !do j =2,nr-1
    !    call Th2_coefficients(aTh2_w,aTh2_e,sTh2,r(j),Snut(j),dnutdr(j),lambda(j),dlambdadT(j),d2lambdadT2(j),dTdr(j),d2Tdr2(j), &
    !    Pr,sigmaTh2,deta,d2etadr2(j),detadr(j),eps(j),Kt(j))
    !    Th2(j) = sTh2 + aTh2_e*Th2(j+1) + aTh2_w*Th2(j-1)
    !enddo

    do j =nr-1,1,-1
        Th2(j) = A(j) * Th2(j+1) + C_apex(j)
    enddo
    

    end
