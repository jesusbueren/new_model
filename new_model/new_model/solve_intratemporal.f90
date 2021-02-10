subroutine solve_intratemporal(p2_fc,cts,h,l_i,u,l_fc,c)
    use structural_p1; use structural_p2
    implicit none
    real(SP),intent(in)::cts,l_i,p2_fc
    integer,intent(in)::h
    real(SP),intent(out)::u,l_fc,c
    real(SP)::l_fc_min,l_fc_max,l_fc_g,mu_c,mu_l_fc
    real(SP)::TINY=1.0e-3
    integer::it
    character::end_k
    
    l_fc_min=0.001_sp
    !Check corner at l_fc_min
    l_fc_g=l_fc_min
    c=cts
    mu_c=p2_fc*c**(-sigma)
    if (subs_p==0.0_sp) then
        mu_l_fc=mu(h)*share_p*l_fc_g**(share_p-1.0_sp)*l_i**(1.0_sp-share_p)* &
                (l_fc_g**share_p*l_i**(1.0_sp-share_p))**(-nu)
    else
        mu_l_fc=mu(h)*share_p*l_fc_g**(subs_p-1.0_sp)* &
                (share_p*l_fc_g**subs_p+(1.0_sp-share_p)*l_i**subs_p)**((1.0_sp-nu-subs_p)/subs_p)
    end if
    if (mu_c>mu_l_fc .or. mu(h)==0.0_sp) then
        go to 1
    end if
    
    
    !Look for the interior solution   
    it=0
    l_fc_max=cts/p2_fc
2   l_fc_g=(l_fc_min+l_fc_max)/2.0_sp
    c=cts-p2_fc*l_fc_g
    mu_c=p2_fc*c**(-sigma)
    if (subs_p==0.0_sp) then
        mu_l_fc=mu(h)*share_p*l_fc_g**(share_p-1.0_sp)*l_i**(1.0_sp-share_p)* &
                (l_fc_g**share_p*l_i**(1.0_sp-share_p))**(-nu)
    else
        mu_l_fc=mu(h)*share_p*l_fc_g**(subs_p-1.0_sp)* &
                (share_p*l_fc_g**subs_p+(1.0_sp-share_p)*l_i**subs_p)**((1.0_sp-nu-subs_p)/subs_p)
    end if
    if (abs(mu_c-mu_l_fc)/mu_c>TINY .and. it<2000) then
        if (mu_c>mu_l_fc) then
            l_fc_max=l_fc_g
        else
            l_fc_min=l_fc_g
        end if
        it=it+1
        go to 2 
    elseif (it==5000) then
        print*, 'intratemporal not converged',mu_c,mu_l_fc,l_fc_g,cts,l_fc_max,l_fc_min
        print*, 'intratemporal not converged',mu(h)*share_p*l_fc_g**(subs_p-1.0_sp)* &
            (share_p*l_fc_g**subs_p+(1.0_sp-share_p)*l_i**subs_p)**((1.0_sp-nu-subs_p)/subs_p)
        !read*,end_k
    end if
1   l_fc=l_fc_g
    
    if (mu(h)<=0.0_sp)then
        u=c**(1.0_sp-sigma)/(1.0_sp-sigma)
    else
        if (subs_p==0.0_sp) then
            u=c**(1.0_sp-sigma)/(1.0_sp-sigma)+mu(h)*(l_fc_g**share_p*l_i**(1.0_sp-share_p))**(1.0_sp-nu)/(1.0_sp-nu)
        else
            u=c**(1.0_sp-sigma)/(1.0_sp-sigma)+mu(h)*((share_p*l_fc_g**subs_p+(1.0_sp-share_p)*l_i**subs_p)**(1.0_sp/subs_p))**(1.0_sp-nu)/(1.0_sp-nu)
        end if
    end if
    
    if (isnan(u))then
        print*,'problem in intratemporal'
        print*,delta_h(h),c,sigma,mu(h),l_fc,omega*l_i,kappa_h(h),nu
        !read*,end_k
    end if
end subroutine
    