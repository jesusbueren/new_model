subroutine SMM_se(parameters,moments,densities)
    use dimensions; use structural_p2;use structural_p1
    implicit none
    real(SP),dimension(moment_conditions,1),intent(out)::moments,densities
    real(SP),dimension(parameters_to_est),intent(in)::parameters
    integer::h_l
    integer,dimension(nkk,clusters,nzz,L_gender,L_PI2,f_t,generations)::a_policy,g_policy
    real(SP),dimension(nkk,clusters,nzz,L_gender,L_PI2,f_t,generations)::beq100_policy
    real(SP),dimension(nkk,clusters,f_t,L_PI2)::lfc_x
    character::pause_k
    real(SP)::c,l_fc,mu_av
    real(SP),dimension(nkk,clusters,f_t,L_PI2)::u_x
    
    H_av(clusters+1,clusters+1,:,:,:)=1.0_sp
    
    sigma=parameters(1)
    x_bar(1)=parameters(2)
    alpha_mu(1:3)=parameters(3:5)
    delta(1:f_t)=parameters(6)
    nu=parameters(7)
    lambda(1:f_t)=exp(parameters(8:9)) 
    share_p=parameters(10)
    subs_p=parameters(11)
        
    sigma_beq=sigma
    beta=0.95_sp
    omega=1.0_sp 
    kappa_h=0.0_sp
    delta_h=1.0_sp
    
    !Utility at the floor when healthy
    do h_l=1,1
        call solve_intratemporal(p_fc,x_bar(h_l),h_l,0.0_sp,u_bar_no_f(h_l),l_fc,c)
    end do
    u_bar_no_f=u_bar_no_f(1)
    
    print*,beta,sigma,nu,delta(1),x_bar(1), &
           lambda,alpha_mu,share_p,subs_p
    
    mu=0.0_sp
    do h_l=2,clusters
        mu(h_l)=exp(alpha_mu(h_l-1))
    end do
    
    !Solve the model given a set of parameters
    call solve_model(a_policy,g_policy,lfc_x,u_x,beq100_policy)
    
    !Simulate model given the policy function to obtain moments and densities
    print*,'Simulating Model'
    call simulate_model_se(a_policy,g_policy,lfc_x,beq100_policy, &
                            moments,densities)

end subroutine
    