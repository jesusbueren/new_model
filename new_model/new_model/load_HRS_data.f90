    subroutine load_hrs_data()
    use simulation_input; use files_savings
    implicit none
    integer:: i_l
    real(SP),dimension(8,indv)::data_all
    real(SP),dimension(clusters*2+1,obs,indv)::data_h
    real(SP),dimension(indv,obs,clusters*2+1):: pr_h_i
    
    !Load HRS data from asset moments.do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    open(unit=9,file=fdir_inputs//'input_sim.txt')
        read(9,*) data_all
    close(9)
    gender_i=int(data_all(1,:))
    generation_i=(data_all(2,:)-70)/2+1
    PI_q_i=int(data_all(3,:))
    IC_pr_i=reshape(data_all(4:5,:), (/indv, f_t/), order = (/ 2, 1 /))
    x_i=data_all(6,:)/1000.0_sp
    iwendy_i=int(data_all(7,:))
    PI_q_i2=int(data_all(8,:))
    
    open(unit=9,file=fdir_inputs//'input_h.txt')
        read(9,*) data_h
    close(9)
    
    pr_h_i=reshape(data_h, (/indv, obs,clusters*2+1/), order = (/ 3,2, 1 /))
    s_h_i=pr_h_i(:,:,1:clusters)
    f_h_i=pr_h_i(:,:,clusters+1:2*clusters)
    dead_i=pr_h_i(:,:,2*clusters+1)

    do i_l=1,indv
        if (data_all(2,i_l)<=74) then
            group_i(i_l)=1
        elseif (data_all(2,i_l)<=79) then
            group_i(i_l)=2
        elseif (data_all(2,i_l)<=84) then
            group_i(i_l)=3
        else
            group_i(i_l)=4
        end if
    end do
    
    
    
end subroutine