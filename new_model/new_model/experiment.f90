subroutine experiment(parameters_original)
    use dimensions;use nrtype; use simulation_input;use structural_p1;use MD_reform; use grids; use files_savings
    implicit none
    real(SP),dimension(parameters_to_est),intent(in)::parameters_original
    real(SP),dimension(parameters_to_est)::parameters
    real(SP)::r_max,r_min
    real(SP),dimension(f_t,generations-1,2)::p50_75_assets_ic_age
    real(SP),dimension(generations-1,2)::p50_75_assets_all_age
    real(SP),dimension(L_PI,generations-1)::med_assets_pi_age,med_c_pi_age
    character::pause_k
    integer::t_l,p_l, ge_l,i_l,h_l,ps_l,ts_l,ind_h
    real(SP),dimension(L_PI+1)::EDP_bnk,CV_bnk

    !Select only healthy individuals
    ind_h=1
    
    !Benchmark model
    print*,'++++++++++++++++++++++++++++++++++++++'
    print*,'++++++++++++++++++++++++++++++++++++++'
    print*,'Simulating benchmark'
    parameters=parameters_original
    call simulate_HRS_70(parameters,p50_75_assets_ic_age,p50_75_assets_all_age,med_c_pi_age,EDP_bnk,CV_bnk,ind_h)
    open(unit=9,file='benchmark_h.txt')
    do t_l=1,generations-1
        write(9,'(F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)') p50_75_assets_all_age(t_l,1), &
                    p50_75_assets_all_age(t_l,2), &
                    p50_75_assets_ic_age(1,t_l,2),&
                    p50_75_assets_ic_age(2,t_l,2),&
                    med_c_pi_age(L_PI,t_l)
    end do
    close(9)
   
    !Break correlation health and survival pr
    print*,'++++++++++++++++++++++++++++++++++++++'
    print*,'++++++++++++++++++++++++++++++++++++++'
    print*,'Simulating Break correlation health and survival pr'
    !Upload health transition probabilities
    open(unit=9,file=path_h//'H_av2.txt')
        read(9,*)  H_av
    close(9)
    parameters=parameters_original
    call simulate_HRS_70(parameters,p50_75_assets_ic_age,p50_75_assets_all_age,med_c_pi_age,EDP_bnk,CV_bnk,ind_h)
    open(unit=9,file='health_exp.txt')
    do t_l=1,generations-1
        write(9,'(F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)') p50_75_assets_all_age(t_l,1), &
                    p50_75_assets_all_age(t_l,2), &
                    p50_75_assets_ic_age(1,t_l,2),&
                    p50_75_assets_ic_age(2,t_l,2),&
                    med_c_pi_age(L_PI,t_l)
    end do
    close(9)
    open(unit=9,file=path_h//'H_av.txt')
        read(9,*)  H_av
    close(9)
    
end subroutine