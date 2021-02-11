subroutine data_moments_W()
    use nrtype; use dimensions; use simulation_input; use pdfs; use structural_p1; use pdfs;use targets;use structural_p2
    implicit none   
    integer::i_l,ns,j_l,real_moments,nrot
    
    !Initialize value of moments
    data_NW_PI1=-9.0_sp
    data_NW_PI1b=-9.0_sp
    data_NW_PI=-9.0_sp
    data_NW_PIb=-9.0_sp
    data_beq100_IC=-9.0_sp
    data_NW_IC1=-9.0_sp
    data_NW_IC=-9.0_sp
    data_lfc_PI=-9.0_sp
    data_lfc_IC=-9.0_sp
    pdf_all_NW_PI=-9.0_sp
    data_NW_h_ut=-9.0_sp
    
    call load_hrs_data()
    !Load input for simulation and compute vector of data moments
    print*,'main sample'
    call charge_simulation_input_moments(data_NW_PI1,data_NW_PI, & !networth by PI  for cohorts obs in 1998
                                         data_NW_PI1b,data_NW_PIb,& !networth by PI for cohorts after 1998
                                         data_beq100_IC,& ! pr of large beq by IC
                                         data_NW_IC1,data_NW_IC,& ! Networth by IC
                                         data_lfc_PI,& !Formal care by PI
                                         data_lfc_IC,& !Formal care by informal groups
                                         data_NW_h_ut) !NW by health

    data_moments1(:,1)=(/reshape(data_NW_PI1,(/L_PI*obs*groups,1/)),& 
                         reshape(data_NW_PI1b,(/L_PI*obs*groups,1/)),& 
                         reshape(data_beq100_IC,(/L_PI*f_t,1/)), &
                         reshape(data_NW_IC1,(/f_t*obs*groups,1/)), &
                         reshape(data_lfc_PI,(/L_PI*clusters,1/)), &
                         reshape(data_lfc_IC,(/f_t*clusters,1/))/)

    open(unit=9,file='data_moments1.txt')
        write(9,*) data_moments1
    close(9)

    open(unit=9,file='data_moments_ut.txt')
        write(9,*) data_NW_IC1
    close(9)
    
    !select moments to untarget
    data_NW_PIb=-9.0_sp
    data_NW_IC=-9.0_sp
    
    data_moments(:,1)=(/reshape(data_NW_PI,(/L_PI*obs*groups,1/)),&
                        reshape(data_NW_PIb,(/L_PI*obs*groups,1/)),&
                        reshape(data_beq100_IC,(/L_PI*f_t,1/)),&
                        reshape(data_NW_IC,(/f_t*obs*groups,1/)),&
                        reshape(data_lfc_PI,(/L_PI*clusters,1/)),&
                        reshape(data_lfc_IC,(/f_t*clusters,1/))/)
    

    do i_l=1,moment_conditions
        if (data_moments1(i_l,1)==-9.0_sp) then
            data_moments(i_l,1)=-9.0_sp
        end if
    end do
    call empty_missing(data_moments,data_moments_new,int(moment_conditions),real_moments)
        
    !Compute optimal weighting matrix using bootstrap   
    call optimal_W()
    open(unit=9,file='W_opt.txt')
        write(9,*) W_opt
    close(9)
    open(unit=9,file='W_opt.txt')
        read(9,*) W_opt
    close(9)
    
    !!Or identity matrix
    !W_opt=0.0_sp
    !do i_l=1,moment_conditions
    !    W_opt(i_l,i_l)=1.0_sp
    !end do
    
end subroutine
    
subroutine empty_missing(arr,arr_new,moment_conditions,real_moments)
    use nrtype
        implicit none
        integer,intent(in)::moment_conditions
        integer,intent(out)::real_moments
        real(SP),intent(in)::arr(moment_conditions,1)
        real(SP),intent(out)::arr_new(moment_conditions,1)
        integer::ind,i_l
        ind=0
        arr_new=-9.0_sp
        do i_l=1,moment_conditions
            if (arr(i_l,1)/=-9.0_sp) then
                ind=ind+1
                arr_new(ind,1)=arr(i_l,1)
            end if
        end do   
        real_moments=ind
    end subroutine
    
    
subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
real(kind=8) :: a(n,n), c(n,n)
real(kind=8) :: L(n,n), U(n,n), b(n), d(n), x(n)
real(kind=8) :: coeff
integer :: i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse
