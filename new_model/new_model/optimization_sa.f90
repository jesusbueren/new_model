subroutine optimization_sa(parameters_end)
    use dimensions; use structural_p2; use nrtype; use optimization
    implicit none
    real(SP),dimension(parameters_to_est),intent(out)::parameters_end
    real(SP),dimension(parameters_to_est)::p_g
    real(SP),dimension(parameters_to_est,2)::limits
    real(SP),dimension(parameters_to_est,parameters_to_est)::xi
    real(SP),dimension(parameters_to_est+1,parameters_to_est)::parameters
    real(SP),dimension((parameters_to_est+1)*parameters_to_est)::parameters_v
    real(SP),dimension(parameters_to_est+1)::obj_fct
    real(SP)::ftol=0.0005_sp,temptr,obj_fct_end
    real(DP),dimension(parameters_to_est)::u
    interface
        function SMM(parameters)
        use nrtype;use dimensions
        implicit none
        real(SP),dimension(parameters_to_est),intent(in)::parameters
        real(SP)::SMM
        end function
    end interface
    integer::p_l,iter,p_l2,it
    integer,parameter:: INIT=100
    real(SP),dimension(INIT,parameters_to_est):: final_p
    real(SP),dimension(INIT):: value_SMM
    character::end_key
    !integer,allocatable::seed(:)
    integer,dimension(parameters_to_est)::seed=254
    integer::M
    
    !Set seed
    !call random_seed(size=M)
    !ALLOCATE (seed(M))
    call random_seed(PUT=seed)
    
    !Set initial simplex
    !
    !!Solve for INIT different initial conditions
    !!sigma
    !limits(1,1)=1.2_sp
    !limits(1,2)=4.0_sp
    !!transfer in medicaid
    !limits(2,1)=1.0_sp
    !limits(2,2)=20.0_sp
    !!alpha mu
    !limits(3,1)=1.0_sp
    !limits(3,2)=1.2_sp
    !limits(4,1)=7.0_sp
    !limits(4,2)=8.0_sp
    !limits(5,1)=10.0_sp
    !limits(5,2)=11.0_sp
    !!delta
    !limits(6,1)=50.0_sp
    !limits(6,2)=1000.0_sp
    !!nu
    !limits(7,1)=1.2_sp
    !limits(7,2)=7.3_sp
    !!lambda
    !limits(8,1)=1.0_sp
    !limits(8,2)=7.0_sp
    !
    !!!sigma ep
    !!limits(14:17,1)=0.001_sp
    !!limits(14:17,2)=0.01_sp
    !!Beta
    !!limits(14,1)=0.90_sp
    !!limits(14,2)=0.98_sp
    !
    !
    !open(unit=11,file='param_op.txt')
    !open(unit=12,file='value_op.txt')
    !do it=1,INIT
    !    print*,'optimization',it,'out of',INIT
    !    do p_l=1,parameters_to_est+1
    !        call random_seed(PUT=seed) 
    !        call RANDOM_NUMBER(u)
    !        call random_seed(GET=seed)
    !        do p_l2=1,parameters_to_est
    !            parameters(p_l,p_l2)=limits(p_l2,1)+(limits(p_l2,2)-limits(p_l2,1))*real(u(p_l2))
    !        end do
    !        call p2R(parameters(p_l,:))
    !        obj_fct(p_l)=SMM(parameters(p_l,:))
    !        !read*,end_key
    !    end do
    !    print*,'dont forget to change parameter number in amoeba if changed the number of parameters'
    !    call amoeba(parameters,obj_fct,ftol,SMM,iter)
    !    call R2p(parameters(1,:),parameters_end)
    !    final_p(it,:)=parameters_end
    !    value_SMM(it)=obj_fct(1)
    !    print*,final_p(it,:)
    !    write(11,*) final_p(it,:)
    !    write(12,*) value_SMM(it)
    !end do
    !parameters_end=final_p(minloc(value_SMM,1),:)
    print*,'Optimization started'
    open(unit=11,file='initial_simplex.txt')
        read(11,*) parameters_v
    close(11)
    parameters=reshape(parameters_v,(/parameters_to_est+1,parameters_to_est/),order=(/2,1/))
    do p_l=1,parameters_to_est+1
        print*,'p_l',p_l
        call p2R(parameters(p_l,:))
        obj_fct(p_l)=SMM(parameters(p_l,:))        
        !print*,'pause'
        !read*,end_key
    end do
    print*,'amoeba started'
    call amoeba(parameters,obj_fct,ftol,SMM,iter)
    call R2p(parameters(1,:),parameters_end)

end subroutine 
    
    