subroutine solve_intratemporal(p2_fc,cts,h,z,l_i,u,l_fc,c)
    use structural_p1; use structural_p2
    implicit none
    real(SP),intent(in)::cts,l_i,p2_fc
    integer,intent(in)::h,z
    real(SP),intent(out)::u,l_fc,c
    real(SP)::l_fc_min,l_fc_max,l_fc_g,mu_c,mu_l_fc
    real(SP)::TINY=1.0e-3
    integer::it
    character::end_k
    
    l_fc_min=max(0.0_sp,-omega*l_i-kappa_h(h))
    !Check corner at l_fc_min
    l_fc_g=l_fc_min
    c=cts-p2_fc*l_fc_g
    mu_c=c
    mu_l_fc=p2_fc**(1.0_sp/sigma)*(mu(h,z)/delta_h(h))**(-1.0_sp/sigma)*(l_fc_g+omega*l_i+kappa_h(h))**(nu/sigma)
    if (mu_c<mu_l_fc .or. mu(h,z)==0.0_sp) then
        go to 1
    end if
    
    !Check corner at l_fc_max
    l_fc_g=cts/p2_fc
    if (cts/p2_fc+omega*l_i+kappa_h(h)<0.0_sp)then
        go to 1
    end if
    
    !Look for the interior solution   
    it=0
    l_fc_max=cts/p2_fc
2   l_fc_g=(l_fc_min+l_fc_max)/2
    c=cts-p2_fc*l_fc_g
    mu_c=c
    mu_l_fc=p2_fc**(1.0_sp/sigma)*(mu(h,z)/delta_h(h))**(-1.0_sp/sigma)*(l_fc_g+omega*l_i+kappa_h(h))**(nu/sigma)
    if (abs(mu_c-mu_l_fc)/mu_c>TINY .and. it<2000) then
        if (mu_c<mu_l_fc) then
            l_fc_max=l_fc_g
        else
            l_fc_min=l_fc_g
        end if
        it=it+1
        go to 2 
    elseif (it==2000) then
        print*, 'intratemporal not converged'
        !read*,end_k
    end if
1   l_fc=l_fc_g
    
    if (mu(h,z)<=0.0_sp)then
        u=delta_h(h)*c**(1.0_sp-sigma)/(1.0_sp-sigma)
    elseif (cts/p2_fc+omega*l_i+kappa_h(h)<0.0_sp) then
        u=-1.0_sp/0.0_sp
    else
        u=delta_h(h)*c**(1.0_sp-sigma)/(1.0_sp-sigma)+mu(h,z)*(l_fc+omega*l_i+kappa_h(h))**(1.0_sp-nu)/(1.0_sp-nu)
    end if
    
    if (isnan(u))then
        print*,'problem in intratemporal'
        print*,delta_h(h),c,sigma,mu(h,z),l_fc,omega*l_i,kappa_h(h),nu
        !read*,end_k
    end if
end subroutine
    
subroutine solve_intratemporal_av(p2_fc,cts,h,mu_av,l_i,u,l_fc,c)
    use structural_p1; use structural_p2
    implicit none
    real(SP),intent(in)::cts,l_i,p2_fc,mu_av
    integer,intent(in)::h
    real(SP),intent(out)::u,l_fc,c
    real(SP)::l_fc_min,l_fc_max,l_fc_g,mu_c,mu_l_fc
    real(SP)::TINY=1.0e-3
    integer::it
    
    l_fc_min=max(0.0_sp,-omega*l_i-kappa_h(h))
    !Check corner at l_fc_min
    l_fc_g=l_fc_min
    c=cts-p2_fc*l_fc_g
    mu_c=c
    mu_l_fc=p2_fc**(1.0_sp/sigma)*(mu_av/delta_h(h))**(-1.0_sp/sigma)*(l_fc_g+omega*l_i+kappa_h(h))**(nu/sigma)
    if (mu_c<mu_l_fc .or. mu_av==0.0_sp) then
        go to 1
    end if
    
    !Check corner at l_fc_max
    l_fc_g=cts/p2_fc
    if (cts/p2_fc+omega*l_i+kappa_h(h)<0.0_sp)then
        go to 1
    end if
    
    !Look for the interior solution   
    it=0
    l_fc_max=cts/p2_fc
2   l_fc_g=(l_fc_min+l_fc_max)/2
    c=cts-p2_fc*l_fc_g
    mu_c=c
    mu_l_fc=p2_fc**(1.0_sp/sigma)*(mu_av/delta_h(h))**(-1.0_sp/sigma)*(l_fc_g+omega*l_i+kappa_h(h))**(nu/sigma)
    if (abs(mu_c-mu_l_fc)/mu_c>TINY .and. it<2000) then
        if (mu_c<mu_l_fc) then
            l_fc_max=l_fc_g
        else
            l_fc_min=l_fc_g
        end if
        it=it+1
        go to 2 
    elseif (it==2000) then
        print*, 'intratemporal not converged'
    end if
1   l_fc=l_fc_g
    
    if (mu_av<=0.0_sp)then
        u=delta_h(h)*c**(1.0_sp-sigma)/(1.0_sp-sigma)
    elseif (cts/p2_fc+omega*l_i+kappa_h(h)<0.0_sp) then
        u=-1.0_sp/0.0_sp
    else
        u=delta_h(h)*c**(1.0_sp-sigma)/(1.0_sp-sigma)+mu_av*(l_fc+omega*l_i+kappa_h(h))**(1.0_sp-nu)/(1.0_sp-nu)
    end if
    
    if (isnan(u))then
        print*,'problem in intratemporal'
    end if
end subroutine    