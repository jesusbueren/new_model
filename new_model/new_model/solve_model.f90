subroutine solve_model(a_policy,g_policy,lfc_x,u_x,beq100_policy)
    use dimensions;use nrtype; use structural_p2; use grids; use structural_p1; use MD_reform
    implicit none
    real(SP),dimension(nkk,clusters,nzz,L_gender,L_PI2,f_t,generations),intent(out)::beq100_policy
    integer,dimension(nkk,clusters,nzz,L_gender,L_PI2,f_t,generations),intent(out)::a_policy,g_policy  !Policy for assets
    real(SP),dimension(nkk,clusters,f_t,L_PI2),intent(out)::u_x
    real(SP),dimension(nkk,clusters,f_t,L_PI2),intent(out)::lfc_x
    real(SP),dimension(nkk,clusters+1,nzz,2)::V,beq_aux
    real(SP),dimension(nkk,clusters,f_t,L_PI2)::c_x
    integer::x_l,f_l,h_l,k_l2,t_l,i_l,ps_l,ge_l,k2_wo_MD,i_l2,k2_l_min
    real(SP)::V_wo_MD,V_MD,ECV_k2_no_beq,V_k2,beq_wo_md,beq_MD
    character::pause_k
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
    do i_l=1,L_PI2;do f_l=1,f_t;do h_l=1,clusters;
        do x_l=1,nkk;
        if (x_l==1) then
            u_x(x_l,h_l,f_l,i_l)=-1.0_sp/0.0_sp
            lfc_x(x_l,h_l,f_l,i_l)=0.0_sp
            c_x(x_l,h_l,f_l,i_l)=0.0_sp
        else
            call solve_intratemporal(p_fc,coh_grid(x_l),h_l,l_ic(f_l,h_l),u_x(x_l,h_l,f_l,i_l),lfc_x(x_l,h_l,f_l,i_l),c_x(x_l,h_l,f_l,i_l))
        end if
    end do; end do; end do; end do

    print*,'vfi'
    !Solve intertemporal problem by standard VFI
    a_policy=-9
    g_policy=-9
    beq100_policy=-9
    call tick(calc)
    
    !$OMP PARALLEL default(none) private(i_l,t_l,ps_l,h_l,x_l,k2_l_min,V_k2,ECV_k2_no_beq,V_MD,V_wo_MD,k2_wo_MD,k_l2,ge_l,f_l,V,beq_aux,beq_wo_md,beq_md,pause_k) shared(g_policy,beq100_policy,a_policy,u_x,u_bar,coh_grid,delta,lambda,sigma,sigma_beq,omega,beta,V_70)
    !$OMP  DO collapse(3)
    do ge_l=1,L_gender;
    do f_l=1,f_t 
    do i_l=1,L_PI2
        !Solve for the value of leaving bequests
        beq_aux=-9.0_sp
        V=-9.0_sp
        do x_l=1,nkk
            V(x_l,clusters+1,:,1:2)=lambda(f_l)*(coh_grid(x_l)+delta(f_l))**(1.0_sp-sigma_beq)/(1.0_sp-sigma_beq)
            if (coh_grid(x_l)>0.0d0) then
                beq_aux(x_l,clusters+1,:,1:2)=1.0_sp
            else
                beq_aux(x_l,clusters+1,:,1:2)=0.0_sp
            end if
        end do
    do t_l=generations,1,-1 
    do ps_l=1,nzz
    do h_l=1,clusters
    do x_l=1,nkk
        !Std VFI
        if (x_l>2) then
            k2_l_min=a_policy(x_l-1,h_l,ps_l,ge_l,i_l,f_l,t_l) 
        else
            k2_l_min=1
        end if

        call std_vfi(u_x,k2_l_min,x_l,ge_l,i_l,f_l,t_l,h_l,ps_l,V,V_wo_MD,k2_wo_MD,beq_aux,beq_wo_md)
        !Value of Medicaid
        k_l2=1
        call ECV_V_k_l2(u_x(x_l-k_l2+1,h_l,f_l,i_l),V,ge_l,i_l,f_l,t_l,h_l,ps_l,x_l,k_l2,V_k2,ECV_k2_no_beq,beq_aux,beq_MD)
        V_MD=u_bar(h_l,f_l)+beta*ECV_k2_no_beq
        !if (isnan(V_MD)) then
        !    print*,''
        !    call ECV_V_k_l2(u_x(x_l-k_l2+1,h_l,f_l,i_l),V,ge_l,i_l,f_l,t_l,h_l,ps_l,x_l,k_l2,V_k2,ECV_k2_no_beq)
        !end if
        !Discrete choice
        if (V_MD>=V_wo_MD) then
            V(x_l,h_l,ps_l,1)=V_MD
            a_policy(x_l,h_l,ps_l,ge_l,i_l,f_l,t_l)=1
            g_policy(x_l,h_l,ps_l,ge_l,i_l,f_l,t_l)=1
            beq100_policy(x_l,h_l,ps_l,ge_l,i_l,f_l,t_l)=beq_MD
        else
            V(x_l,h_l,ps_l,1)=V_wo_MD
            a_policy(x_l,h_l,ps_l,ge_l,i_l,f_l,t_l)=k2_wo_MD
            g_policy(x_l,h_l,ps_l,ge_l,i_l,f_l,t_l)=0
            beq100_policy(x_l,h_l,ps_l,ge_l,i_l,f_l,t_l)=beq_wo_md 
        end if

        !if (V(x_l,h_l,ps_l,1)==-1.0d0/0.0d0 ) then
        !    print*,''
        !end if   
    end do !x_l
    end do !h_l
    end do !ps_l
    V(:,:,:,2)=V(:,:,:,1)
    beq_aux(:,1:clusters,:,2)=beq100_policy(:,:,:,ge_l,i_l,f_l,t_l)
    end do !t_l
    V_70(:,:,:,ge_l,i_l,f_l)=V(:,:,:,1)
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