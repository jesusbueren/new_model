subroutine simulate_i(a_policy,g_policy,lfc_x,samples,h_ini,x_ini,f_l,PI_q_ii,PI_q_ii2,gender_ii,mean_fc_age,mean_assets,mean_c)
    use dimensions; use nrtype;use structural_p1;use grids; use structural_p2; use simulation_input
    implicit none
    integer,dimension(nkk,clusters,nzz,L_gender,L_PI2,f_t,generations),intent(in)::a_policy,g_policy
    real(SP),dimension(nkk,clusters,f_t,L_PI2),intent(in)::lfc_x
    integer,intent(in)::samples,h_ini,f_l,PI_q_ii,PI_q_ii2,gender_ii
    real(SP),intent(in)::x_ini
    real(SP),dimension(generations),intent(out)::mean_assets,mean_c
    real(SP),dimension(samples,generations)::assets,cons
    real(SP),dimension(generations)::fc_age,x_it,median_assets,counter_impaired,a_it
    integer,dimension(generations)::counter_age,h_s,g_it,counter_fc_age
    integer::h_l,s_l,ind,t_l,xi_l,xi_l2,ts_l2,pos_x,k2_l
    real(SP)::u
    real(SP),dimension(generations),intent(out)::mean_fc_age
    character::pause_k
    
    counter_age=0
    counter_fc_age=0
    counter_impaired=0.0_sp
    fc_age=0.0_sp
    assets=-9.0_sp
    cons=-9.0_sp
    
    do s_l=1,samples
        !Set the individual you want to simulate
        a_it(1)=x_ini
        !Sample health type
        h_s=-9
        do t_l=1,generations-1
            if (t_l==1) then
                h_s(1)=h_ini
            else 
                call RANDOM_NUMBER(u)
                ind=1
                do while (h_s(t_l)==-9)
                    if ((u<sum(H_av(h_s(t_l-1),1:ind,t_l,PI_q_ii,gender_ii))) .or. (ind==clusters+1)) then
                            h_s(t_l)=ind
                    else
                        ind=ind+1
                    end if
                end do
            end if
        end do
        !Simulate decision
        do t_l=1,generations-1
            if (h_s(t_l)<clusters+1) then
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
                !Store assets
                counter_age(t_l)=counter_age(t_l)+1
                assets(counter_age(t_l),t_l)=a_it(t_l)
                !Cash on hand
                x_it(t_l)=(1+r)*a_it(t_l)-m_exp_all(t_l,gender_ii,PI_q_ii,h_s(t_l),xi_l,ts_l2)+b(PI_q_ii2,gender_ii)
                pos_x=int(max(min(x_it(t_l),coh_grid(nkk)),0.0_sp)/(coh_grid(2)-coh_grid(1))+1.00000001_sp)
                !Draw on the discrete choice
                call RANDOM_NUMBER(u)
                if (g_policy(pos_x,h_s(t_l),xi_l,gender_ii,PI_q_ii2,f_l,t_l)==1) then
                    g_it(t_l)=1
                    k2_l=1
                    cons(counter_age(t_l),t_l)=c_bar(h_s(t_l),1)
                else
                    g_it(t_l)=0
                    k2_l=a_policy(pos_x,h_s(t_l),xi_l,gender_ii,PI_q_ii2,f_l,t_l)
                    cons(counter_age(t_l),t_l)=coh_grid(pos_x)-coh_grid(k2_l)-lfc_x(pos_x-k2_l+1,h_s(t_l),f_l,PI_q_ii2)*p_fc(1)
                end if
                a_it(t_l+1)=coh_grid(k2_l)
                !Formal care decision
                if (g_it(t_l)==0) then
                    counter_fc_age(t_l)=counter_fc_age(t_l)+1
                    fc_age(t_l)=fc_age(t_l)+lfc_x(pos_x-k2_l+1,h_s(t_l),f_l,PI_q_ii2)
                end if
                !counter impaired
                if (h_s(t_l)==4)then
                    counter_impaired(t_l)=counter_impaired(t_l)+1.0_sp
                end if
            end if
        end do
    end do
    
    !print*,'average hours of fc per day over age'
    
    do t_l=1,generations
        mean_fc_age(t_l)=fc_age(t_l)/real(counter_fc_age(t_l))
        !print*,t_l,counter_age(t_l),mean_fc_age(t_l)/2.0_sp/365.0_sp
    end do
    close(9)
    
    !print*,'mean and median assets'
    mean_assets=-9.0_sp
    mean_c=-9.0_sp
    do t_l=1,generations
        if (counter_age(t_l)>=1) then
            mean_assets(t_l)=sum(assets(1:counter_age(t_l),t_l))/counter_age(t_l)
            mean_c(t_l)=sum(cons(1:counter_age(t_l),t_l))/counter_age(t_l)
            counter_impaired(t_l)=counter_impaired(t_l)/real(counter_age(t_l))
        end if
    end do
    
    open(unit=9,file='sim_i.txt')
    do t_l=1,generations
        write(9,*) mean_fc_age(t_l),mean_assets(t_l),counter_impaired(t_l)
    end do
    close(9)
    !
    !print*,'paused in simulate_i'
    !read*,pause_k
    
end subroutine