subroutine counterfactuals(parameters_original)
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
    
    !Select all individuals
    ind_h=clusters
    
    !Benchmark model
    print*,'++++++++++++++++++++++++++++++++++++++'
    print*,'++++++++++++++++++++++++++++++++++++++'
    print*,'Simulating benchmark'
    ind_or=1
    parameters=parameters_original
    call simulate_HRS_70(parameters,p50_75_assets_ic_age,p50_75_assets_all_age,med_c_pi_age,EDP_bnk,CV_bnk,ind_h)
    open(unit=9,file='benchmark.txt')
    do t_l=1,generations-1
        write(9,'(F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)') p50_75_assets_all_age(t_l,1), &
                    p50_75_assets_all_age(t_l,2), &
                    p50_75_assets_ic_age(1,t_l,1),&
                    p50_75_assets_ic_age(2,t_l,1),&
                    p50_75_assets_ic_age(1,t_l,2),&
                    p50_75_assets_ic_age(2,t_l,2),&
                    med_c_pi_age(L_PI,t_l)
    end do
    close(9)
    
    !print*,'++++++++++++++++++++++++++++++++++++++'
    !print*,'++++++++++++++++++++++++++++++++++++++'
    !print*,'Simulating no ltc'
    !ind_or=0
    !parameters=parameters_original
    !!LTC need shifter equal to zero
    !parameters(3:5)=-1.0_sp/0.0_sp
    !p_nh=0.0d0
    !call simulate_HRS_70(parameters,p50_75_assets_ic_age,p50_75_assets_all_age,med_c_pi_age,EDP_bnk,CV_bnk,ind_h)
    !open(unit=9,file='noLTC.txt')
    !do t_l=1,generations-1
    !    write(9,'(F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)') p50_75_assets_all_age(t_l,1), &
    !                p50_75_assets_all_age(t_l,2), &
    !                p50_75_assets_ic_age(1,t_l,1),&
    !                p50_75_assets_ic_age(2,t_l,1),&
    !                p50_75_assets_ic_age(1,t_l,2),&
    !                p50_75_assets_ic_age(2,t_l,2),&
    !                med_c_pi_age(L_PI,t_l)
    !end do
    !close(9)
    !p_nh=(/8.0_sp/365.0_sp*443.0_sp,8.0_sp/365.0_sp*437.0_sp,8.0/365.0_sp*463.0_sp,8.0/365.0_sp*527.0_sp,8.0/365.0_sp*527.0_sp/2.0_sp/)
    !
    !!Counterfactual for w/o bequest motives
    !print*,'++++++++++++++++++++++++++++++++++++++'
    !print*,'++++++++++++++++++++++++++++++++++++++'
    !print*,'Simulating no bequests'
    !ind_or=0
    !parameters=parameters_original
    !parameters(8:9)=-1.0_sp/0.0_sp
    !p_nh=(/8.0_sp/365.0_sp*443.0_sp,8.0_sp/365.0_sp*437.0_sp,8.0/365.0_sp*463.0_sp,8.0/365.0_sp*527.0_sp,8.0/365.0_sp*527.0_sp/2.0_sp/)
    !call simulate_HRS_70(parameters,p50_75_assets_ic_age,p50_75_assets_all_age,med_c_pi_age,EDP_bnk,CV_bnk,ind_h)
    !open(unit=9,file='noBeq.txt')
    !do t_l=1,generations-1
    !    write(9,'(F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)') p50_75_assets_all_age(t_l,1), &
    !                p50_75_assets_all_age(t_l,2), &
    !                p50_75_assets_ic_age(1,t_l,1),&
    !                p50_75_assets_ic_age(2,t_l,1),&
    !                p50_75_assets_ic_age(1,t_l,2),&
    !                p50_75_assets_ic_age(2,t_l,2),&
    !                med_c_pi_age(L_PI,t_l)
    !end do
    !close(9)
    !
    !
    !!Counterfactual for no med
    !print*,'++++++++++++++++++++++++++++++++++++++'
    !print*,'++++++++++++++++++++++++++++++++++++++'
    !print*,'Simulating no Med exp'
    !ind_or=0
    !parameters=parameters_original
    !m_exp_all_or=m_exp_all
    !m_exp_all=0.0_sp
    !call simulate_HRS_70(parameters,p50_75_assets_ic_age,p50_75_assets_all_age,med_c_pi_age,EDP_bnk,CV_bnk,ind_h)
    !open(unit=9,file='noMed.txt')
    !do t_l=1,generations-1
    !    write(9,'(F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)') p50_75_assets_all_age(t_l,1), &
    !                p50_75_assets_all_age(t_l,2), &
    !                p50_75_assets_ic_age(1,t_l,1),&
    !                p50_75_assets_ic_age(2,t_l,1),&
    !                p50_75_assets_ic_age(1,t_l,2),&
    !                p50_75_assets_ic_age(2,t_l,2),&
    !                med_c_pi_age(L_PI,t_l)
    !end do
    !close(9)   
    !m_exp_all=m_exp_all_or
    !
    !20% expansion Consumption floor
    print*,'++++++++++++++++++++++++++++++++++++++'
    print*,'++++++++++++++++++++++++++++++++++++++'
    print*,'Simulating policy reform: change medicaid provision'
    ind_or=0
    parameters=parameters_original
    call simulate_HRS_70(parameters,p50_75_assets_ic_age,p50_75_assets_all_age,med_c_pi_age,EDP_bnk,CV_bnk,ind_h)
    
    !!No close families
    !print*,'++++++++++++++++++++++++++++++++++++++'
    !print*,'++++++++++++++++++++++++++++++++++++++'
    !print*,'Simulating no female'
    !ind_or=0
    !parameters=parameters_original
    !ind_no_f=1
    !call load_hrs_data()
    !call simulate_HRS_70(parameters,p50_75_assets_ic_age,p50_75_assets_all_age,med_c_pi_age,EDP_bnk,CV_bnk,ind_h)
    !open(unit=9,file='noFemale.txt')
    !do t_l=1,generations-1
    !    write(9,'(F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)') p50_75_assets_all_age(t_l,1), &
    !                p50_75_assets_all_age(t_l,2), &
    !                p50_75_assets_ic_age(1,t_l,1),&
    !                p50_75_assets_ic_age(2,t_l,1),&
    !                p50_75_assets_ic_age(1,t_l,2),&
    !                p50_75_assets_ic_age(2,t_l,2),&
    !                med_c_pi_age(L_PI,t_l)
    !end do
    !close(9) 
    !
    !print*,'++++++++++++++++++++++++++++++++++++++'
    !print*,'++++++++++++++++++++++++++++++++++++++'
    !print*,'Simulating no ltc no close families'
    !ind_or=0
    !parameters=parameters_original
    !!LTC need shifter equal to zero
    !parameters(3:5)=-1.0_sp/0.0_sp
    !p_nh=0.0d0
    !call simulate_HRS_70(parameters,p50_75_assets_ic_age,p50_75_assets_all_age,med_c_pi_age,EDP_bnk,CV_bnk,ind_h)
    !open(unit=9,file='noFemale_noLTC.txt')
    !do t_l=1,generations-1
    !    write(9,'(F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)') p50_75_assets_all_age(t_l,1), &
    !                p50_75_assets_all_age(t_l,2), &
    !                p50_75_assets_ic_age(1,t_l,1),&
    !                p50_75_assets_ic_age(2,t_l,1),&
    !                p50_75_assets_ic_age(1,t_l,2),&
    !                p50_75_assets_ic_age(2,t_l,2),&
    !                med_c_pi_age(L_PI,t_l)
    !end do
    !close(9)
    !
    !ind_no_f=0
    !call load_hrs_data()
    !p_nh=(/8.0_sp/365.0_sp*443.0_sp,8.0_sp/365.0_sp*437.0_sp,8.0/365.0_sp*463.0_sp,8.0/365.0_sp*527.0_sp,8.0/365.0_sp*527.0_sp/2.0_sp/)
    
    
end subroutine