subroutine medicaid_provision()
    use nrtype; use dimensions; use structural_p2; use structural_p1
    implicit none
    integer::h_l,f_l,z_l,it
    real(SP)::g,R_max,R_min,u_g,R_g
        
        !Find hours of care and consumption floor provided by Medicaid 
        !if no care from relatives to get to the floor for the different shocks
        do h_l=1,clusters; do z_l=1,nzz2
            !guess of max budget
            it=0
            g=1.0_sp
    1       R_max=(30.0_sp+p_fc*2.0_sp*365.0_sp*2.0_sp)*g
            call solve_intratemporal(p_fc,R_max,h_l,z_l,0.0_sp,u_g,l_bar(h_l,z_l),c_bar(h_l,z_l))
            if (u_g<u_bar_no_f(h_l) .and. it<2000) then
                g=g*2.0_sp
                it=it+1
                go to 1
            elseif (it==2000) then
                print*, 'Problem in medicaid provision'
            end if
            R_min=0.0_sp
            !Bisection
            it=0
    2       R_g=(R_max+R_min)/2.0_sp
            call solve_intratemporal(p_fc,R_g,h_l,z_l,0.0_sp,u_g,l_bar(h_l,z_l),c_bar(h_l,z_l))
            if (abs((u_g-u_bar_no_f(h_l))/u_bar_no_f(h_l))>0.0001_sp .and. it<2000) then
                if (u_g>u_bar_no_f(h_l)) then
                    R_max=R_g
                else
                    R_min=R_g
                end if
                it=it+1
                go to 2
            elseif (it==2000) then
                print*, 'medicaid provision not converged'
            end if
        end do; end do
        

        !Utility of medicaid for individuals
        do h_l=1,clusters; do f_l=1,f_t;do z_l=1,nzz2
            if (h_l==1) then
                u_bar(h_l,f_l,z_l)=c_bar(h_l,z_l)**(1.0_sp-sigma)/(1.0_sp-sigma)
            else
                u_bar(h_l,f_l,z_l)=c_bar(h_l,z_l)**(1.0_sp-sigma)/(1.0_sp-sigma)+mu(h_l,z_l)*(l_bar(h_l,z_l)+omega*l_ic(f_l,h_l)+kappa_h(h_l))**(1.0_sp-nu)/(1.0_sp-nu)
            end if
        end do;end do;end do

end subroutine