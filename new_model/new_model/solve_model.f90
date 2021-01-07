subroutine solve_model(a_policy,g_policy,lfc_x,u_x,beq100_policy)
    use dimensions;use nrtype; use structural_p2; use grids; use structural_p1; use MD_reform
    implicit none
    real(SP),dimension(nkk,clusters,nzz,nzz2,L_gender,L_PI2,f_t,generations),intent(out)::g_policy,beq100_policy
    integer,dimension(nkk,clusters,nzz,nzz2,L_gender,L_PI2,f_t,generations),intent(out)::a_policy  !Policy for assets
    real(SP),dimension(nkk,clusters,nzz2,f_t,L_PI2),intent(out)::u_x
    real(SP),dimension(nkk,clusters,nzz2,f_t,L_PI2),intent(out)::lfc_x
    real(SP),dimension(nkk,clusters+1,nzz,nzz2,2)::V
    real(SP),dimension(nkk,clusters,nzz2,f_t,L_PI2)::c_x
    integer::x_l,f_l,h_l,k_l2,t_l,i_l,ps_l,ge_l,k2_wo_MD,i_l2,k2_l_min,z_l
    real(SP)::V_wo_MD,V_MD,ECV_k2_no_beq,V_k2
    !Timer
    integer::calc
    real::calctime
    interface
        real function tock(t)
            integer, intent(in) :: t
        end function tock
    end interface
    
    !Compute medicaid provision
    !print*,'computing medicaid provision'
    call medicaid_provision()
    
    
    !Solve for the intratemporal formal care/consumption decision given a cash to spend coh_grid(x_l)
    !print*,'Solving intratemporal'
    do i_l=1,L_PI2;do f_l=1,f_t;do h_l=1,clusters;do z_l=1,nzz2; do x_l=1,nkk;
        if (x_l==1) then
            u_x(x_l,h_l,z_l,f_l,i_l)=-1.0_sp/0.0_sp
            lfc_x(x_l,h_l,z_l,f_l,i_l)=0.0_sp
            c_x(x_l,h_l,z_l,f_l,i_l)=0.0_sp
        else
            call solve_intratemporal(p_fc,coh_grid(x_l),h_l,z_l,l_ic(f_l,h_l),u_x(x_l,h_l,z_l,f_l,i_l),lfc_x(x_l,h_l,z_l,f_l,i_l),c_x(x_l,h_l,z_l,f_l,i_l))
        end if
    end do; end do; end do; end do;end do
    
    print*,'vfi'
    !Solve intertemporal problem by standard VFI
    a_policy=-9
    g_policy=-9
    beq100_policy=-9
    call tick(calc)
    
    !$OMP PARALLEL default(none) private(i_l,t_l,ps_l,h_l,x_l,z_l,k2_l_min,V_k2,ECV_k2_no_beq,V_MD,V_wo_MD,k2_wo_MD,k_l2,ge_l,f_l,V) shared(g_policy,a_policy,u_x,u_bar,coh_grid,delta,lambda,sigma,sigma_beq,omega,beta,V_70,sigma_varep)
    !$OMP  DO collapse(3)
    do ge_l=1,L_gender;
    do f_l=1,f_t 
    do i_l=1,L_PI2
        !Solve for the value of leaving bequests
        V=-9.0_sp
        do x_l=1,nkk
            V(x_l,clusters+1,:,:,1:2)=lambda(f_l)*(coh_grid(x_l)+delta(f_l))**(1.0_sp-sigma_beq)/(1.0_sp-sigma_beq)
        end do
    do t_l=generations,1,-1 
    do z_l=1,nzz2
    do ps_l=1,nzz
    do h_l=1,clusters
    do x_l=1,nkk
        !Std VFI
        if (x_l>2) then
            k2_l_min=a_policy(x_l-1,h_l,ps_l,z_l,ge_l,i_l,f_l,t_l) 
        else
            k2_l_min=1
        end if
        call std_vfi(u_x,k2_l_min,x_l,ge_l,i_l,f_l,t_l,h_l,ps_l,z_l,V,V_wo_MD,k2_wo_MD)
        !Value of Medicaid
        k_l2=1
        call ECV_V_k_l2(u_x(x_l-k_l2+1,h_l,z_l,f_l,i_l),V,ge_l,i_l,f_l,t_l,h_l,ps_l,x_l,k_l2,V_k2,ECV_k2_no_beq)
        V_MD=u_bar(h_l,f_l,z_l)+beta*ECV_k2_no_beq
        !if (isnan(V_MD)) then
        !    print*,''
        !    call ECV_V_k_l2(u_x(x_l-k_l2+1,h_l,z_l,f_l,i_l),V,ge_l,i_l,f_l,t_l,h_l,ps_l,x_l,k_l2,V_k2,ECV_k2_no_beq)
        !end if
        !Discrete choice
        V(x_l,h_l,ps_l,z_l,1)=sigma_varep(h_l)*log(exp(V_MD/sigma_varep(h_l))+exp(V_wo_MD/sigma_varep(h_l)))
        if (V(x_l,h_l,ps_l,z_l,1)==-1.0d0/0.0d0 .or. V(x_l,h_l,ps_l,z_l,1)==1.0d0/0.0d0 .or. isnan(V(x_l,h_l,ps_l,z_l,1))) then
            if (V_MD>=V_wo_MD) then
                V(x_l,h_l,ps_l,z_l,1)=V_MD
                a_policy(x_l,h_l,ps_l,z_l,ge_l,i_l,f_l,t_l)=1
                g_policy(x_l,h_l,ps_l,z_l,ge_l,i_l,f_l,t_l)=1.0_sp
            else
                V(x_l,h_l,ps_l,z_l,1)=V_wo_MD
                a_policy(x_l,h_l,ps_l,z_l,ge_l,i_l,f_l,t_l)=k2_wo_MD
                g_policy(x_l,h_l,ps_l,z_l,ge_l,i_l,f_l,t_l)=0.0_sp
            end if
        else
            g_policy(x_l,h_l,ps_l,z_l,ge_l,i_l,f_l,t_l)=exp(V_MD/sigma_varep(h_l))/(exp(V_MD/sigma_varep(h_l))+exp(V_wo_MD/sigma_varep(h_l)))
            a_policy(x_l,h_l,ps_l,z_l,ge_l,i_l,f_l,t_l)=k2_wo_MD
        end if
        !if (V(x_l,h_l,ps_l,z_l,1)==-1.0d0/0.0d0 ) then
        !    print*,''
        !end if   
    end do !x_l
    end do !h_l
    end do !ps_l
    end do !z_l
    V(:,:,:,:,2)=V(:,:,:,:,1)
    end do !t_l
    V_70(:,:,:,:,ge_l,i_l,f_l)=V(:,:,:,:,1)
    end do !i_l
    end do !f_l
    end do !ge_l
    !$OMP END DO  
    !$OMP END PARALLEL
    
    calctime = tock(calc)
    print *,'Timing summary'
    print *,'Calc: ', calctime
end subroutine
    
    
subroutine tick(t)
    integer, intent(OUT) :: t
    call system_clock(t)
end subroutine tick

! returns time in seconds from now to time described by t
real function tock(t)
    integer, intent(in) :: t
    integer :: now, clock_rate
    call system_clock(now,clock_rate)
    tock = real(now - t)/real(clock_rate)
end function tock