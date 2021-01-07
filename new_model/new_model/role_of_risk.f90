!subroutine role_of_risk(parameters)
!use dimensions; use nrtype; use structural_p1; use structural_p2; use grids; use simulation_input
!implicit none
!real(SP),dimension(parameters_to_est),intent(in)::parameters
!integer,dimension(1)::seed=254
!real(SP)::c,l_fc,mu_av,u
!integer,dimension(nkk,clusters,nzz,nzz2,L_gender,L_PI2,f_t,generations)::a_policy,a_policy_new
!real(SP),dimension(nkk,clusters,nzz,nzz2,L_gender,L_PI2,f_t,generations)::g_policy,g_policy_new
!real(SP),dimension(nkk,clusters,nzz2,f_t,L_PI2)::u_x,u_x_new
!real(SP),dimension(nkk,clusters,nzz2,f_t,L_PI2)::lfc_x,lfc_x_new
!integer,parameter::samples1=100,samples2=50
!integer::h_l,z_l,s_l,f_l,ind,t_l,i_l
!real(SP),dimension(generations)::mean_fc_age,mean_fc_age_no_ltc,mean_assets2,mean_c2
!real(SP),dimension(generations,L_gender,L_PI2,clusters,nzz,nzz)::m_exp_all_or
!real(SP),dimension(generations,indv,samples2)::mean_assets_it,mean_c_it
!real(SP),dimension(f_t,generations,L_gender,L_PI,clusters,nzz,nzz)::m_exp_all_f
!real(SP)::x_ini
!integer::h_ini,PI_q_ii,PI_q_ii2,gender_ii
!character::pause_k
!
!!simulate the model 1000 times for a particular indv in the sample and compute the expected
!!expenses in LTC
!    !Set the seed
!    call random_seed(PUT=seed) 
!    
!    m_exp_all_or=m_exp_all
!    H_av(clusters+1,clusters+1,:,:,:)=1.0_sp
!    
!    sigma=parameters(1)
!    nu=parameters(2)
!    delta(1:f_t)=parameters(3)
!    x_bar(1:clusters)=parameters(4:7)
!    lambda(1:3)=exp(parameters(8:10))
!    alpha_mu(1:3)=parameters(11:13)
!    beta=parameters(14)
!    omega=parameters(15)
!    sigma_varep=parameters(16)
!    sigma_beq=sigma
!    sigma2_varep(1:3)=0.0_sp
!    kappa_h=0.0_sp
!    delta_h=1.0_sp
!    
!    do h_l=1,clusters
!        if (h_l==1) then
!            mu_av=0.0_sp
!        elseif (h_l==2) then
!            mu_av=exp(alpha_mu(1)+sigma2_varep(1)/2.0_sp)
!        elseif (h_l==3) then
!            mu_av=exp(alpha_mu(2)+sigma2_varep(2)/2.0_sp)
!        elseif (h_l==4) then
!            mu_av=exp(alpha_mu(3)+sigma2_varep(3)/2.0_sp)
!        end if
!        call solve_intratemporal_av(p_fc,x_bar(h_l),h_l,mu_av,0.0_sp,u_bar_no_f(h_l),l_fc,c)
!    end do
!    
!    !Discretize var_ep and compute preference shifter for care
!    if (nzz2>1) then
!        do h_l=2,clusters
!            call discretize_var_ep(sigma2_varep(h_l-1),nzz2,varep_grid(:,h_l-1),pr_varep)
!        end do
!    else
!        varep_grid=0.0_sp
!        pr_varep=1.0_sp
!    end if
!    mu=0.0_sp
!    do z_l=1,nzz2; do h_l=2,clusters
!        mu(h_l,z_l)=exp(alpha_mu(h_l-1)+varep_grid(z_l,h_l-1))
!    end do;end do
!    
!    !Solve the model given a set of parameters
!    call solve_model(a_policy,g_policy,lfc_x,u_x)
!    !
!    !do i_l=1,indv
!    !    print*,'individual',i_l,'out of ',indv
!    !    if (group_i(i_l)==1 .and. PI_q_i(i_l)==4) then
!    !        !Solve the model for each indv and each family type
!    !        do f_l=1,f_t
!    !            !Simulate the model
!    !            h_ini=maxloc(s_h_i(i_l,1,:),1) 
!    !            x_ini=min(x_i(i_l),coh_grid(nkk))
!    !            PI_q_ii=PI_q_i(i_l)
!    !            PI_q_ii2=PI_q_i2(i_l)
!    !            gender_ii=gender_i(i_l)
!    !            m_exp_all=m_exp_all_or
!    !            call simulate_i(a_policy,g_policy,lfc_x,samples1,h_ini,x_ini,f_l,PI_q_ii,PI_q_ii2,gender_ii,mean_fc_age,mean_assets2,mean_c2)
!    !            !Take out average formal care expenses and set the m.u. of care to zero
!    !            do t_l=1,generations
!    !                m_exp_all(t_l,:,:,:,:,:)=m_exp_all_or(t_l,:,:,:,:,:)+mean_fc_age(t_l)*p_fc
!    !            end do
!    !            m_exp_all_f(f_l,:,:,:,:,:,:)=m_exp_all
!    !            mu=0.0_sp
!    !            call solve_model_no_risk(a_policy_new,g_policy_new,lfc_x_new,u_x_new,f_l)
!    !        end do
!    !        !Simulate the model sample2 times for each indv taking a random draw of its family type
!    !        !print*,'simulation wo risk'
!    !        do s_l=1,samples2
!    !            !Sample the family type 
!    !            call RANDOM_NUMBER(u)
!    !            if (u<IC_pr_i(i_l,1)) then
!    !                f_l=1
!    !            elseif (u<sum(IC_pr_i(i_l,1:2))) then
!    !                f_l=2
!    !            else
!    !                f_l=3
!    !            end if
!    !            m_exp_all=m_exp_all_f(f_l,:,:,:,:,:,:)
!    !            call simulate_i(a_policy_new,g_policy_new,lfc_x_new,1,h_ini,x_ini,f_l,PI_q_ii,PI_q_ii2,gender_ii,mean_fc_age,mean_assets_it(:,i_l,s_l),mean_c_it(:,i_l,s_l))
!    !        end do
!    !    end if
!    !end do
!    !
!    !open(unit=9,file='mean_assets_it.txt')
!    !    write(9,*) mean_assets_it
!    !close(9)
!    !open(unit=9,file='mean_c_it.txt')
!    !    write(9,*) mean_c_it
!    !close(9)
!    
!    !Compute median assets and mean consumption
!    open(unit=9,file='mean_assets_it.txt')
!        read(9,*) mean_assets_it
!    close(9)
!    open(unit=9,file='mean_c_it.txt')
!        read(9,*) mean_c_it
!    close(9)
!    
!    call compute_distribution(mean_assets_it,mean_c_it,samples2)
!    
!    !Average Formal Care expenses paid by female top pi quartile
!    h_ini=1
!    x_ini=280
!    PI_q_ii=4
!    PI_q_ii2=7
!    gender_ii=2
!    f_l=2
!    m_exp_all=m_exp_all_or
!    call simulate_i(a_policy,g_policy,lfc_x,10000,h_ini,x_ini,f_l,PI_q_ii,PI_q_ii2,gender_ii,mean_fc_age,mean_assets2,mean_c2)
!    open(unit=9,file='fc_exp_f.txt')
!    do t_l=1,generations
!        write(9,*),mean_fc_age(t_l)*p_fc
!    end do
!    close(9)
!    
!    
!end subroutine
!    
!subroutine solve_model_no_risk(a_policy,g_policy,lfc_x,u_x,f_i)
!    use dimensions;use nrtype; use structural_p2; use grids; use structural_p1; use MD_reform
!    implicit none
!    integer,intent(in)::f_i
!    integer,dimension(nkk,clusters,nzz,nzz2,L_gender,L_PI2,f_t,generations),intent(inout)::a_policy  !Policy for Medicaid & for assets
!    real(SP),dimension(nkk,clusters,nzz,nzz2,L_gender,L_PI2,f_t,generations),intent(inout)::g_policy
!    real(SP),dimension(nkk,clusters,nzz2,f_t,L_PI2),intent(inout)::u_x
!    real(SP),dimension(nkk,clusters,nzz2,f_t,L_PI2),intent(inout)::lfc_x
!    real(SP),dimension(nkk,clusters+1,nzz,nzz2,2)::V
!    real(SP),dimension(nkk,clusters,nzz2,f_t,L_PI2)::c_x
!    integer::x_l,f_l,h_l,k_l2,t_l,i_l,ps_l,ge_l,k2_wo_MD,i_l2,k2_l_min,z_l
!    real(SP)::V_wo_MD,V_MD,ECV_k2_no_beq,V_k2
!    !Timer
!    integer::calc
!    real::calctime
!    interface
!        real function tock(t)
!            integer, intent(in) :: t
!        end function tock
!    end interface
!    
!    
!    !Solve for the intratemporal formal care/consumption decision given a cash to spend coh_grid(x_l)
!    !print*,'Solving intratemporal'
!    do i_l=1,L_PI2;do f_l=f_i,f_i;do h_l=1,clusters;do z_l=1,nzz2; do x_l=1,nkk;
!        if (x_l==1) then
!            u_x(x_l,h_l,z_l,f_l,i_l)=-1.0_sp/0.0_sp
!            lfc_x(x_l,h_l,z_l,f_l,i_l)=0.0_sp
!            c_x(x_l,h_l,z_l,f_l,i_l)=0.0_sp
!        else
!            call solve_intratemporal(p_fc,coh_grid(x_l),h_l,z_l,l_ic(f_l,h_l),u_x(x_l,h_l,z_l,f_l,i_l),lfc_x(x_l,h_l,z_l,f_l,i_l),c_x(x_l,h_l,z_l,f_l,i_l))
!        end if
!    end do; end do; end do; end do;end do
!    
!    print*,'vfi'
!    !Solve intertemporal problem by standard VFI
!    call tick(calc)
!    
!    f_l=f_i
!    !$OMP PARALLEL default(none) private(i_l,t_l,ps_l,h_l,x_l,z_l,k2_l_min,V_k2,ECV_k2_no_beq,V_MD,V_wo_MD,k2_wo_MD,k_l2,ge_l,V) shared(g_policy,a_policy,u_x,u_bar,coh_grid,delta,lambda,beta,sigma,V_70,f_l)
!    !$OMP  DO collapse(2)
!    do ge_l=1,L_gender;
!    do i_l=7,L_PI2
!        !Solve for the value of leaving bequests
!        V=-9.0_sp
!        do x_l=1,nkk
!            V(x_l,clusters+1,:,:,1:2)=lambda(f_l)*(coh_grid(x_l)+delta(f_l))**(1.0_sp-sigma)/(1.0_sp-sigma)
!        end do
!    do t_l=generations,1,-1 
!    do z_l=1,nzz2
!    do ps_l=1,nzz
!    do h_l=1,clusters
!    do x_l=1,nkk
!        !Std VFI
!        if (x_l>2) then
!            k2_l_min=a_policy(x_l-1,h_l,ps_l,z_l,ge_l,i_l,f_l,t_l) 
!        else
!            k2_l_min=1
!        end if
!        call std_vfi(u_x,k2_l_min,x_l,ge_l,i_l,f_l,t_l,h_l,ps_l,z_l,V,V_wo_MD,k2_wo_MD)
!        !Value of Medicaid
!        k_l2=1
!        call ECV_V_k_l2(u_x(x_l-k_l2+1,h_l,z_l,f_l,i_l),V,ge_l,i_l,f_l,t_l,h_l,ps_l,x_l,k_l2,V_k2,ECV_k2_no_beq)
!        V_MD=u_bar(h_l,f_l,z_l)+beta*ECV_k2_no_beq
!        !Agent chooses the maximum
!        if (V_MD>=V_wo_MD)then
!            V(x_l,h_l,ps_l,z_l,1)=V_MD
!            a_policy(x_l,h_l,ps_l,z_l,ge_l,i_l,f_l,t_l)=1
!            g_policy(x_l,h_l,ps_l,z_l,ge_l,i_l,f_l,t_l)=1
!        else
!            V(x_l,h_l,ps_l,z_l,1)=V_wo_MD
!            a_policy(x_l,h_l,ps_l,z_l,ge_l,i_l,f_l,t_l)=k2_wo_MD
!            g_policy(x_l,h_l,ps_l,z_l,ge_l,i_l,f_l,t_l)=0
!        end if
!    end do !x_l
!    end do !h_l
!    end do !ps_l
!    end do !z_l
!    V(:,:,:,:,2)=V(:,:,:,:,1)
!    end do !t_l
!    V_70(:,:,:,:,ge_l,i_l,f_l)=V(:,:,:,:,1)
!    end do !i_l
!    end do !ge_l
!    !$OMP END DO  
!    !$OMP END PARALLEL
!    
!    calctime = tock(calc)
!    print *,'Timing summary'
!    print *,'Calc: ', calctime
!    end subroutine
!    
!subroutine compute_distribution(mean_assets_it,mean_c_it,samples2)
!use simulation_input; use nrtype; use structural_p1
!implicit none
!integer,intent(in)::samples2
!real(SP),dimension(generations,indv,samples2),intent(in)::mean_assets_it,mean_c_it
!integer,parameter::indv_c=3000
!real(SP),dimension(generations,indv_c)::assets_all_age,c_all_age
!integer::i_l2,s_l,i_l,t_l,ind
!integer,dimension(generations)::counter_all_age,h_s
!real(sp)::u
!real(SP),dimension(generations-1,2)::p50_75_assets_all_age
!real(SP),dimension(generations-1)::p50_c_all_age
!
!counter_all_age=0
!do i_l2=1,indv_c
!    h_s=-9
!1   call RANDOM_NUMBER(u)
!    u=u*(real(indv-1)+0.5_sp)+1.0_sp
!    i_l=int(u)
!    if (group_i(i_l)/=1 .or. PI_q_i(i_l)/=4) then
!        go to 1
!    end if
!    call RANDOM_NUMBER(u)
!    u=u*(real(samples2-1)+0.5_sp)+1.0_sp
!    s_l=int(u)
!    do t_l=1,generations-1
!        if (mean_assets_it(t_l,i_l,s_l)/=-9.0_sp) then
!            if (isnan(mean_assets_it(t_l,i_l,s_l))) then
!            else
!            counter_all_age(t_l)=counter_all_age(t_l)+1
!            assets_all_age(t_l,counter_all_age(t_l))=mean_assets_it(t_l,i_l,s_l)
!            c_all_age(t_l,counter_all_age(t_l))=mean_c_it(t_l,i_l,s_l)
!            end if
!        end if
!    end do
!end do
!
!    do t_l=1,generations
!        print*,t_l
!        if (counter_all_age(t_l)>1) then
!            !p50_75_assets_all_age(t_l,1)=sum(assets_all_age(t_l,1:counter_all_age(t_l)))/counter_all_age(t_l)
!            call compute_percentile(assets_all_age(t_l,1:counter_all_age(t_l)), &
!                                    counter_all_age(t_l),&
!                                    50, & 
!                                    p50_75_assets_all_age(t_l,1))
!            call compute_percentile(assets_all_age(t_l,1:counter_all_age(t_l)), &
!                                    counter_all_age(t_l),&
!                                    75, & 
!                                    p50_75_assets_all_age(t_l,2))
!            !p50_c_all_age(t_l)=sum(c_all_age(t_l,1:counter_all_age(t_l)))/counter_all_age(t_l)
!            call compute_percentile(c_all_age(t_l,1:counter_all_age(t_l)), &
!                                    counter_all_age(t_l),&
!                                    50, & 
!                                    p50_c_all_age(t_l))
!        end if
!    end do
!    open(unit=9,file='no_risk.txt')
!    do t_l=1,18
!        write(9,*) p50_75_assets_all_age(t_l,1),p50_75_assets_all_age(t_l,2),p50_c_all_age(t_l)
!    end do
!    close(9)
!    
!end subroutine
!
!
!    
!    
!    
!    
!    
