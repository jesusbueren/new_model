subroutine SMM_se(parameters,moments,densities)
    use dimensions; use structural_p2;use structural_p1
    implicit none
    real(SP),dimension(moment_conditions,1),intent(out)::moments,densities
    real(SP),dimension(parameters_to_est),intent(in)::parameters
    integer::h_l,z_l
    integer,dimension(nkk,clusters,nzz,nzz2,L_gender,L_PI2,f_t,generations)::a_policy
    real(SP),dimension(nkk,clusters,nzz,nzz2,L_gender,L_PI2,f_t,generations)::g_policy
    real(SP),dimension(nkk,clusters,nzz2,f_t,L_PI2)::lfc_x
    character::pause_k
    real(SP)::c,l_fc,mu_av
    real(SP),dimension(nkk,clusters,nzz2,f_t,L_PI2)::u_x
    
    H_av(clusters+1,clusters+1,:,:,:)=1.0_sp
    
    sigma=parameters(1)
    x_bar(1)=parameters(2)
    alpha_mu(1:3)=parameters(3:5)
    delta(1:f_t)=parameters(6)
    nu=parameters(7)
    lambda(1:3)=exp(parameters(8)) 
        
    sigma_beq=sigma
    beta=0.95_sp
    omega=1.0_sp 
    sigma_varep(1:4)=0.0000001_sp
    sigma2_varep(1:3)=0.0_sp
    kappa_h=0.0_sp
    delta_h=1.0_sp
    
    !Utility at the floor when healthy
    do h_l=1,1
        if (h_l==1) then
            mu_av=0.0_sp
        elseif (h_l==2) then
            mu_av=exp(alpha_mu(1)+sigma2_varep(1)/2.0_sp)
        elseif (h_l==3) then
            mu_av=exp(alpha_mu(2)+sigma2_varep(2)/2.0_sp)
        elseif (h_l==4) then
            mu_av=exp(alpha_mu(3)+sigma2_varep(3)/2.0_sp)
        end if
        call solve_intratemporal_av(p_fc,x_bar(h_l),h_l,mu_av,0.0_sp,u_bar_no_f(h_l),l_fc,c)
    end do
    u_bar_no_f=u_bar_no_f(1)
    
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
    end do;end do
    
    !Solve the model given a set of parameters
    call solve_model(a_policy,g_policy,lfc_x,u_x)
    
    !Simulate model given the policy function to obtain moments and densities
    print*,'Simulating Model'
    call simulate_model_se(a_policy,g_policy,lfc_x, &
                            moments,densities)

end subroutine
    