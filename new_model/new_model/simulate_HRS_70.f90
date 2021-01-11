subroutine simulate_HRS_70(parameters,p50_75_assets_ic_age,p50_75_assets_all_age,med_c_pi_age,EDP,CV,ind_h)
    use dimensions;use grids; use nrtype; use simulation_input; use structural_p1; use pdfs; use structural_p2; use MD_reform
    use targets; use optimization
    implicit none
    real(SP),dimension(parameters_to_est),intent(in)::parameters
    integer,intent(in)::ind_h
    real(SP),dimension(L_PI+1),intent(out)::EDP,CV
    integer,parameter::indv_c=150000
    real(SP),dimension(nkk,clusters,nzz,L_gender,L_PI2,f_t,generations)::beq100_policy
    integer,dimension(nkk,clusters,nzz,L_gender,L_PI2,f_t,generations)::a_policy,g_policy
    real(SP),dimension(nkk,clusters,f_t,L_PI2)::u_x
    real(SP),dimension(nkk,clusters,f_t,L_PI2)::lfc_x
    integer::i_l,xi_l,i_l2,t_l,ind,k2_l,pos_x,xi_l2,f_l,ts_l2,pi_l,nwq_l,h_l,pos_x2,k2_l2
    real(SP)::u,x_pr_md,sum_t,sum_t_md
    real(SP),dimension(indv_c)::tr_i
    integer,dimension(indv_c)::pi_i2
    integer,dimension(generations)::h_s,g_it
    real(SP),dimension(generations)::x_it,tr_it,a_it
    integer,dimension(L_PI,generations)::counter_pi_age
    integer,dimension(f_t,generations)::counter_ic_age
    integer,dimension(generations)::counter_all_age
    real(SP),dimension(L_PI,generations,indv_c/2)::assets_pi_age,g_tr,c_pi_age
    real(SP),dimension(f_t,generations,indv_c/2)::assets_ic_age
    real(SP),dimension(L_PI,indv_c)::lambda_pi,ltci_pi
    real(SP),dimension(generations,indv_c)::assets_all_age
    real(SP),dimension(L_PI,generations-1)::med_assets_pi_age,med_c_pi_age
    real(SP),dimension(f_t,generations-1,2),intent(out)::p50_75_assets_ic_age
    real(SP),dimension(generations-1,2),intent(out)::p50_75_assets_all_age
    real(SP)::c,l_fc,mu_av
    integer,dimension(1)::seed=254
    character::pause_k
    real(SP)::sub_per_h
    real(SP),dimension(moment_conditions,1)::model_moments1,model_moments, &
                                             model_moments_new,model_moments1_new
    real(SP),dimension(L_PI,obs,groups)::moments_PI_ut
    real(SP),dimension(L_PI,generations)::m_av,obs_m_av
    
    !Set the seed
    call random_seed(PUT=seed) 
    
    H_av(clusters+1,clusters+1,:,:,:)=1.0_sp
    
    sigma=parameters(1)
    x_bar(1)=parameters(2)
    alpha_mu(1:3)=parameters(3:5)
    delta(1:f_t)=parameters(6)
    nu=parameters(7)
    lambda(1:2)=exp(parameters(8:9))   
    sigma_beq=sigma
    beta=0.95_sp
    omega=1.0_sp 
    kappa_h=0.0_sp
    delta_h=1.0_sp
    
    !Utility at the floor when healthy
    do h_l=1,1
        if (h_l==1) then
            mu_av=0.0_sp
        elseif (h_l==2) then
            mu_av=exp(alpha_mu(1))
        elseif (h_l==3) then
            mu_av=exp(alpha_mu(2))
        elseif (h_l==4) then
            mu_av=exp(alpha_mu(3))
        end if
        call solve_intratemporal_av(p_fc,x_bar(h_l),h_l,mu_av,0.0_sp,u_bar_no_f(h_l),l_fc,c)
    end do
    u_bar_no_f=u_bar_no_f(1)
    
    print*,beta,sigma,nu,delta(1),x_bar(1), &
           lambda,alpha_mu,omega
    
    mu=0.0_sp
    do h_l=2,clusters
        mu(h_l)=exp(alpha_mu(h_l-1))
    end do
    
    !Change price of formal care
    sub_per_h=p_fc-p_fc*(1.0_sp-p_sub)
    p_fc=p_fc*(1.0_sp-p_sub)
    
    !Solve the model given a set of parameters
    call solve_model(a_policy,g_policy,lfc_x,u_x,beq100_policy)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print*,'Just for simulation of moments'
    call simulate_model(a_policy,g_policy,lfc_x,beq100_policy, &
                        model_moments1,model_moments)
    print*,'got here2'
    do i_l=1,moment_conditions
        if (data_moments(i_l,1)==-9.0_sp) then
            model_moments(i_l,1)=-9.0_sp
            model_moments1(i_l,1)=-9.0_sp
        end if
    end do
    call empty_missing(model_moments,model_moments_new ,int(moment_conditions),real_moments)
    
    print*,'obj fct',real(indv)/(1.0_dp+1.0_dp/real(samples_per_i))*matmul(matmul(transpose(model_moments_new(1:real_moments,1:1)), & 
                                    W_opt(1:real_moments,1:real_moments)), &
                                    model_moments_new(1:real_moments,1:1))
    print*,'got here3'
    open(unit=9,file='model_moments.txt')
        write(9,*) model_moments1
    close(9)
    print*,'close window, just for simulation of moments'
    read*,pause_k
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (ind_or==1) then
        V_70_or=V_70
    else
        V_70_new=V_70
    end if

    print*,'simulating'
    !Simulate decisions
    counter_pi_age=0
    counter_ic_age=0
    counter_all_age=0
    tr_i=0.0_sp
    lambda_pi=0.0_sp
    ltci_pi=0.0_sp
    m_av=0.0_sp
    obs_m_av=0.0_sp
    do i_l2=1,indv_c
        x_it=-9.0_sp
        a_it=-9.0_sp
        h_s=-9
1       call RANDOM_NUMBER(u)
        u=u*(real(indv-1)+0.5_sp)+1
        i_l=int(u)
        if (group_i(i_l)/=1) then
            go to 1
        end if
        !Sample the family type 
        call RANDOM_NUMBER(u)
        if (u<IC_pr_i(i_l,1)) then
            f_l=1
        elseif (u<sum(IC_pr_i(i_l,1:2))) then
            f_l=2
        else
            f_l=3
        end if
        h_s=-9
        g_it=-9
        a_it(1)=min(x_i(i_l),coh_grid(nkk))
        do t_l=1,generations-1
            !Sample health type
            if (t_l==1) then
                h_s(1)=maxloc(s_h_i(i_l,t_l,:),1)
            else 
                call RANDOM_NUMBER(u)
                ind=1
                do while (h_s(t_l)==-9)
                    if (u<sum(H_av(h_s(t_l-1),1:ind,t_l,PI_q_i2(i_l),gender_i(i_l)))) then 
                            h_s(t_l)=ind
                    else
                        ind=ind+1
                    end if
                end do
            end if
        end do
        !Mortality bias
        if (h_s(16)==clusters+1) then
            go to 1
        end if
        tr_it=0.0_sp
        do t_l=1,generations-1
            if (h_s(t_l)<clusters+1 ) then 
                !Sample initial medical persistent shock to chracterize the state variable
                if (t_l==1) then
                    call RANDOM_NUMBER(u)
                    xi_l=-9
                    ind=1
                    do while (xi_l==-9)
                        if (u<sum(pr0_xi(1:ind,1)).or. ind==nzz) then
                            xi_l=ind
                        else
                            ind=ind+1
                        end if
                    end do                      
                end if
                !Sample persistent & transitory medical shock to simulate expenditures
                xi_l2=-9
                ind=1
                call RANDOM_NUMBER(u)
                do while (xi_l2==-9)
                    if (u<sum(pr_pxi(xi_l,1:ind)) .or. ind==nzz) then
                        xi_l2=ind
                        xi_l=xi_l2
                    else
                        ind=ind+1
                    end if
                end do
                ts_l2=-9
                ind=1
                call RANDOM_NUMBER(u)
                do while (ts_l2==-9 )
                    if (u<sum(pr_ep(1:ind,1)) .or. ind==nzz) then
                        ts_l2=ind
                    else
                        ind=ind+1
                    end if
                end do
                
                !Cash on hand
                x_it(t_l)=(1+r)*a_it(t_l)-m_exp_all(t_l,gender_i(i_l),PI_q_i2(i_l),h_s(t_l),xi_l2,ts_l2)+b(PI_q_i2(i_l),gender_i(i_l))
                obs_m_av(PI_q_i(i_l),t_l)=obs_m_av(PI_q_i(i_l),t_l)+1.0_sp
                m_av(PI_q_i(i_l),t_l)=(obs_m_av(PI_q_i(i_l),t_l)-1.0_sp)/obs_m_av(PI_q_i(i_l),t_l)*m_av(PI_q_i(i_l),t_l) + &
                                        1.0_sp/obs_m_av(PI_q_i(i_l),t_l)*m_exp_all(t_l,gender_i(i_l),PI_q_i2(i_l),h_s(t_l),xi_l2,ts_l2)
                pos_x=int(max(min(x_it(t_l),coh_grid(nkk)),0.0_sp)/(coh_grid(2)-coh_grid(1))+1.00000001_sp)
                !Draw on the discrete choice
                call RANDOM_NUMBER(u)
                if (g_policy(pos_x,h_s(t_l),xi_l,gender_i(i_l),PI_q_i2(i_l),f_l,t_l)==1) then
                    g_it(t_l)=1
                    k2_l=1
                else
                    g_it(t_l)=0
                    k2_l=a_policy(pos_x,h_s(t_l),xi_l,gender_i(i_l),PI_q_i2(i_l),f_l,t_l)
                end if 
                a_it(t_l+1)=coh_grid(k2_l)
                if (g_it(t_l)==1) then
                    tr_it(t_l)=max(c_bar(h_s(t_l))+p_or*l_bar(h_s(t_l))- &
                                x_it(t_l),0.0_sp)
                    x_it(t_l)=0.0_sp
                else
                    tr_it(t_l)=lfc_x(pos_x-k2_l+1,h_s(t_l),f_l,PI_q_i2(i_l))*(p_or-p_fc)
                end if
                tr_it(t_l)=tr_it(t_l)/(1.0_sp+r)**(t_l-1)
                
                !Store assets
                counter_pi_age(PI_q_i(i_l),t_l)=counter_pi_age(PI_q_i(i_l),t_l)+1
                assets_pi_age(PI_q_i(i_l),t_l,counter_pi_age(PI_q_i(i_l),t_l))=a_it(t_l)
                if (g_it(t_l)==0) then
                    c_pi_age(PI_q_i(i_l),t_l,counter_pi_age(PI_q_i(i_l),t_l))=coh_grid(pos_x)-coh_grid(k2_l)-lfc_x(pos_x-k2_l+1,h_s(t_l),f_l,PI_q_i2(i_l))*p_fc
                else
                    c_pi_age(PI_q_i(i_l),t_l,counter_pi_age(PI_q_i(i_l),t_l))=c_bar(h_s(t_l))
                end if
                if (PI_q_i(i_l)==4) then
                    counter_ic_age(f_l,t_l)=counter_ic_age(f_l,t_l)+1
                    assets_ic_age(f_l,t_l,counter_ic_age(f_l,t_l))=a_it(t_l)
                end if
                counter_all_age(t_l)=counter_all_age(t_l)+1
                assets_all_age(t_l,counter_all_age(t_l))=a_it(t_l)
                
                !Compute Compensating Variation
                if (ind_or==0 .and. t_l==1) then
                    if (V_70_or(pos_x,h_s(t_l),xi_l,gender_i(i_l),PI_q_i2(i_l),f_l)==V_70_new(pos_x,h_s(t_l),xi_l,gender_i(i_l),PI_q_i2(i_l),f_l)) then
                        lambda_pi(PI_q_i(i_l),counter_pi_age(PI_q_i(i_l),t_l))=0.0_sp
                        ltci_pi(PI_q_i(i_l),counter_pi_age(PI_q_i(i_l),t_l))=0.0_sp
                    elseif (V_70_or(pos_x,h_s(t_l),xi_l,gender_i(i_l),PI_q_i2(i_l),f_l)<V_70_new(1,h_s(t_l),xi_l,gender_i(i_l),PI_q_i2(i_l),f_l)) then
                        lambda_pi(PI_q_i(i_l),counter_pi_age(PI_q_i(i_l),t_l))=-coh_grid(pos_x) 
                    else
                        lambda_pi(PI_q_i(i_l),counter_pi_age(PI_q_i(i_l),t_l))=&
                        coh_grid(minloc(abs(V_70_new(:,h_s(t_l),xi_l,gender_i(i_l),PI_q_i2(i_l),f_l)-V_70_or(pos_x,h_s(t_l),xi_l,gender_i(i_l),PI_q_i2(i_l),f_l)),1))- &
                        coh_grid(pos_x)
                    end if
                end if
            end if
        end do
        !Sum transfers across periods
        pi_i2(i_l2)=PI_q_i(i_l)
        tr_i(i_l2)=sum(tr_it)
    end do

    do t_l=1,generations-1;
        do pi_l=1,L_PI; 
            if (counter_pi_age(pi_l,t_l)>1) then
                call compute_percentile(assets_pi_age(pi_l,t_l,1:counter_pi_age(pi_l,t_l)), &
                                        counter_pi_age(pi_l,t_l),&
                                        50, & 
                                        med_assets_pi_age(pi_l,t_l))
                call compute_percentile(c_pi_age(pi_l,t_l,1:counter_pi_age(pi_l,t_l)), &
                                        counter_pi_age(pi_l,t_l),&
                                        50, & 
                                        med_c_pi_age(pi_l,t_l)) 
            end if
        end do

        do f_l=1,f_t; 
            if (counter_ic_age(f_l,t_l)>1) then
                call compute_percentile(assets_ic_age(f_l,t_l,1:counter_ic_age(f_l,t_l)), &
                                        counter_ic_age(f_l,t_l),&
                                        50, & 
                                        p50_75_assets_ic_age(f_l,t_l,1))
                call compute_percentile(assets_ic_age(f_l,t_l,1:counter_ic_age(f_l,t_l)), &
                                        counter_ic_age(f_l,t_l),&
                                        75, & 
                                        p50_75_assets_ic_age(f_l,t_l,2))
            end if
        end do
        if (counter_all_age(t_l)>1) then
            call compute_percentile(assets_all_age(t_l,1:counter_all_age(t_l)), &
                                    counter_all_age(t_l),&
                                    50, & 
                                    p50_75_assets_all_age(t_l,1))
            call compute_percentile(assets_all_age(t_l,1:counter_all_age(t_l)), &
                                    counter_all_age(t_l),&
                                    75, & 
                                    p50_75_assets_all_age(t_l,2))
        end if
    end do
 

    do pi_l=1,L_PI+1
        ind=0
        sum_t=0.0_sp
        print*,'pi_l',pi_l
        do i_l=1,indv_c
            if (pi_l<L_PI+1) then
                if (pi_i2(i_l)==pi_l)then
                    ind=ind+1
                    sum_t=sum_t+tr_i(i_l)
                end if
            else
                ind=ind+1
                sum_t=sum_t+tr_i(i_l)
            end if
        end do
        EDP(pi_l)=sum_t/real(ind)
        print*,'EDPV of transfers to income q ',pi_l,': ',EDP(pi_l)
    end do; 
    CV=0.0_sp
    if (ind_or==0)then
        do pi_l=1,L_PI+1
            if (pi_l<=L_PI) then
                CV(pi_l)=sum(lambda_pi(pi_l,:))/real(counter_pi_age(pi_l,1))  
                print*,'LTCI demand', pi_l,sum(ltci_pi(pi_l,:))/real(counter_pi_age(pi_l,1))  
            else
                CV(pi_l)=sum(sum(lambda_pi,2),1)/real(sum(counter_pi_age(:,1)))  
            end if
            print*,'Compensated Variation to income q ',pi_l,': ',CV(pi_l)
        end do
    end if

    
    open(unit=9,file='m_av.txt')
    do t_l=1,generations-1
        write(9,'(F10.2,F10.2,F10.2,F10.2,F10.2)') m_av(1,t_l),m_av(2,t_l),m_av(3,t_l),m_av(4,t_l),m_av(5,t_l)
    end do
    close(9)
        
    !read*,pause_k
    
end subroutine