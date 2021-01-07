function SMM(parameters)
    use dimensions;use nrtype; use simulation_input; use structural_p2; use optimization; use targets
    implicit none
    real(SP),dimension(parameters_to_est),intent(in)::parameters
    integer,dimension(nkk,clusters,nzz,nzz2,L_gender,L_PI2,f_t,generations)::a_policy
    real(SP),dimension(nkk,clusters,nzz,nzz2,L_gender,L_PI2,f_t,generations)::g_policy
    integer::i_l,h_l,z_l
    real(SP),dimension(nkk,clusters,nzz2,f_t,L_PI2)::lfc_x
    real(SP),dimension(moment_conditions,1)::model_moments1,model_moments, &
                                             model_moments_new,model_moments1_new,r2
    real(SP),dimension(1,1)::obj_fct
    real(SP)::SMM,c,l_fc,mu_av
    real(SP),dimension(parameters_to_est)::parameters_n
    real(SP),dimension(f_t,2)::moments_IC_ut
    real(SP),dimension(nkk,clusters,nzz2,f_t,L_PI2)::u_x
    character::pause_key
    
    !Transform real line into parameters
    call R2p(parameters,parameters_n)

    !print*,'new parameter values'
    print*,beta,sigma,nu,delta(1),x_bar(1), &
           lambda,alpha_mu,omega,sigma_varep

    !Discretize var_ep and compute preference shifter for care
    if (nzz2>1) then
        do h_l=2,clusters
            call discretize_var_ep(sigma2_varep(h_l-1),nzz2,varep_grid(:,h_l-1),pr_varep)
            !print*,exp(alpha_mu(h_l-1)+sigma2_varep(h_l-1)/2)
        end do
    else
        varep_grid=0.0_sp
        pr_varep=1.0_sp
    end if
    mu=0.0_sp
    do z_l=1,nzz2; do h_l=2,clusters
        mu(h_l,z_l)=exp(alpha_mu(h_l-1)+varep_grid(z_l,h_l-1))
        !print*,z_l,h_l,alpha_mu(h_l-1),varep_grid(z_l,h_l-1),alpha_mu(h_l-1)+varep_grid(z_l,h_l-1)
    end do;end do
    
    !Solve the model given a set of parameters
    call solve_model(a_policy,g_policy,lfc_x,u_x)

    !Simulate model given the policy function
    !print*,'Simulating Model'
    call simulate_model(a_policy,g_policy,lfc_x, & !policy fcts
                        model_moments1,model_moments) !moments
        
    !Missing moments in the data to missing in the model
    do i_l=1,moment_conditions
        if (data_moments(i_l,1)==-9.0_sp) then
            model_moments(i_l,1)=-9.0_sp
            !model_moments1(i_l,1)=-9.0_sp
        end if
    end do
    
    call empty_missing(model_moments,model_moments_new ,int(moment_conditions),real_moments)
    !call empty_missing(model_moments1,model_moments1_new,int(moment_conditions),real_moments)
        
    !Compute objective function
    obj_fct=real(indv)/(1.0_dp+1.0_dp/real(samples_per_i))*matmul(matmul(transpose(model_moments_new(1:real_moments,1:1)), & 
                                    W_opt(1:real_moments,1:real_moments)), &
                                    model_moments_new(1:real_moments,1:1))
    
    
    if (isnan(obj_fct(1,1))) then
        print*,'problem evaluating objective fct'
        obj_fct(1,1)=200000.0_sp
        print*,'press any key to continue'
        read*,pause_key
    end if
    SMM=obj_fct(1,1)
    r1=0.0_sp
    r1(1:real_moments,1)=0.0!model_moments1_new(1:real_moments,1)
    r1(1:real_moments,1)=0.0!(model_moments_new(1:real_moments,1)-data_moments_new(1:real_moments,1))/data_moments_new(1:real_moments,1)
    
    !Store results
    if (SMM<obj_fct_min) then
        open(unit=9,file='model_moments.txt')
            write(9,*) model_moments1
        close(9)
        open(unit=9,file='parameters_prelim.txt')
            write(9,*) parameters_n
        close(9)
        obj_fct_min=SMM
    end if
        
    print*,'objective fct: ',SMM,' minimum',obj_fct_min

end function