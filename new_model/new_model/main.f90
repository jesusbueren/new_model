program jmp
    use dimensions;use nrtype; use structural_p2; use optimization
    implicit none
    integer,dimension(1)::seed=254
    real(SP),dimension(parameters_to_est)::parameters
    character::end_key
    
    !Set the seed
    call random_seed(PUT=seed) 
    
    !Build grid for cash on hand and shocks
    call build_grids()
    
    !Load 1st step parameters
    call first_step_p()
    
    !Load simulation inputs, data moments and  optimal weigthing matrix
    call data_moments_W()
    
    !!Optimization with Nelder-Mead algorithm
    !Don't forget changing the number of parameters in amoeba.f90 and amebsa.f90
    call optimization_sa(parameters)
    open(unit=9,file='parameters_p18.txt')
        write(9,*) parameters
    close(9)
    
    !Compute standard errors of estimated parameters
    !open(unit=9,file='parameters_p18.txt')
    !    read(9,*) parameters
    !close(9)
    !call compute_se(parameters)
    
    !Create paper's tables in Latex
    !open(unit=9,file='parameters_p18.txt')
    !    read(9,*) parameters
    !close(9)
    !call create_tables_draft(parameters)
    
    
    !Compute the counterfactuals to identify key mechanisms
    open(unit=9,file='parameters_p18.txt')
        read(9,*) parameters
    close(9)
    call counterfactuals(parameters)
    
    !Break correlation between health and death to analyze the effect of the correlation in precautionary savings
    !open(unit=9,file='parameters_p18.txt')
    !    read(9,*) parameters
    !close(9)
    !call experiment(parameters)
    
    
    !Role of LTC risk
    !open(unit=9,file='parameters_p18Wopt.txt')
    !    read(9,*) parameters
    !close(9)
    !call role_of_risk(parameters)
    
    
    print*,'end program, press key to close window'
    read*,end_key
    
end program