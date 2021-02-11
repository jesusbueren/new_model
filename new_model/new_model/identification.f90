subroutine identification(parameters_original)
    use dimensions; use nrtype
    implicit none
    real(SP),dimension(parameters_to_est),intent(in)::parameters_original
    real(SP),dimension(parameters_to_est)::parameters
    real(SP),dimension(moment_conditions,1)::model_moments1
    
    !!Benchmark
    !parameters=parameters_original
    !call p2R(parameters)
    !call get_moments(parameters,model_moments1)
    !open(unit=9,file='model_moments.txt')
    !    write(9,*) model_moments1
    !close(9)
    !
    !!Stronger bequest motives
    !parameters=parameters_original
    !parameters(8:9)=parameters(8:9)*1.20d0
    !call p2R(parameters)
    !call get_moments(parameters,model_moments1)
    !open(unit=9,file='moments_high_beq.txt')
    !    write(9,*) model_moments1
    !close(9)
    
    !Larger LTC needs
    parameters=parameters_original
    parameters(5)=parameters(5)*1.05d0
    call p2R(parameters)
    call get_moments(parameters,model_moments1)
    open(unit=9,file='moments_high_LTC.txt')
        write(9,*) model_moments1
    close(9)
    
end subroutine
    

subroutine get_moments(parameters,model_moments1)  
    use dimensions;use nrtype; use simulation_input; use structural_p2; use optimization; use targets
    implicit none
    real(SP),dimension(parameters_to_est),intent(in)::parameters
    integer,dimension(nkk,clusters,nzz,L_gender,L_PI2,f_t,generations)::a_policy,g_policy
    real(SP),dimension(nkk,clusters,nzz,L_gender,L_PI2,f_t,generations)::beq100_policy
    integer::i_l,h_l,z_l
    real(SP),dimension(nkk,clusters,f_t,L_PI2)::lfc_x
    real(SP),dimension(moment_conditions,1),intent(out)::model_moments1
    real(SP),dimension(moment_conditions,1)::model_moments, &
                                             model_moments_new,model_moments1_new,r2
    real(SP),dimension(1,1)::obj_fct
    real(SP)::SMM,c,l_fc,mu_av
    real(SP),dimension(parameters_to_est)::parameters_n
    real(SP),dimension(f_t,2)::moments_IC_ut
    real(SP),dimension(nkk,clusters,f_t,L_PI2)::u_x
    character::pause_key
    
    !Transform real line into parameters
    call R2p(parameters,parameters_n)

    !print*,'new parameter values'
    print*,beta,sigma,nu,delta(1),x_bar(1), &
           lambda,alpha_mu,share_p,subs_p

    mu=0.0_sp
    do h_l=2,clusters
        mu(h_l)=exp(alpha_mu(h_l-1))
    end do
    
    !Solve the model given a set of parameters
    call solve_model(a_policy,g_policy,lfc_x,u_x,beq100_policy)

    !Simulate model given the policy function
    !print*,'Simulating Model'
    call simulate_model(a_policy,g_policy,lfc_x,beq100_policy, & !policy fcts
                        model_moments1,model_moments) !moments

end subroutine