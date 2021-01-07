subroutine compute_percentile(arr_unsorted,length,percentile,p_percentile)
use nrtype;use dimensions;use grids
implicit none
integer,intent(in)::length,percentile
real(SP),intent(in)::arr_unsorted(length)
real(SP),intent(out)::p_percentile
real(SP)::arr(length)
interface
    function kernel_d(u)
    use nrtype
    implicit none
    real(SP)::kernel_d
    real(SP)::u
    end function
end interface

arr=arr_unsorted
call sort(arr,length)
p_percentile=arr(int(real(length)*real(percentile)/100.0_sp))

end subroutine

subroutine compute_pdf(arr_unsorted,length,point,density_p)
    use nrtype;use dimensions;use grids
    implicit none
    integer,intent(in)::length
    real(SP),intent(in)::arr_unsorted(length),point
    real(SP),intent(out)::density_p
    real(SP)::arr(length)
    real(SP)::h_i(length),pdf_i(length)
    real(SP),dimension(nkk)::x
    real(SP)::h,std
    integer::x_l,x_l2
    character::pause_k
    interface
        function kernel_d(u)
        use nrtype
        implicit none
        real(SP)::kernel_d
        real(SP)::u
        end function
    end interface

    arr=arr_unsorted
    call sort(arr,length)

    !Optimal bandwidth
    std=sum((arr-sum(arr)/dble(length))**2)/(dble(length)-1)
    if (arr(int(real(length)/4.0_sp*3.0_sp))-arr(int(real(length)/4.0_sp))>0.0_sp) then
        h=0.9_sp*real(length)**(-1.0_sp/5.0_sp)*min(std,1.34_sp*(arr(int(real(length)/4.0_sp*3.0_sp))-arr(int(real(length)/4.0_sp))))
    else
        h=0.9_sp*real(length)**(-1.0_sp/5.0_sp)*std
    end if
    do x_l=1,length
        pdf_i(x_l)=0.0_sp
        do x_l2=1,length
            pdf_i(x_l)=pdf_i(x_l)+kernel_d((arr(x_l)-arr(x_l2))/h)
        end do
        pdf_i(x_l)=pdf_i(x_l)/real(length)/h
    end do
    !Variable bandwidth following Silverman 
    h_i=log(pdf_i)-1.0/dble(length)*sum(log(pdf_i))
    h_i=exp(-0.5_sp*h_i)*h !sensitivity parameter
    density_p=0.0_sp
    do x_l2=1,length
        density_p=density_p+1.0_sp/h_i(x_l2)*kernel_d((point-arr(x_l2))/h_i(x_l2))  !sum(exp(-((x(x_l)-arr)/h_i)**2.0_sp/2.0_sp)/sqrt(2.0_sp*PI)/h_i)/real(length)
    end do
    density_p=density_p/real(length)
    
    if (isnan(density_p)) then
        print*,'error in compute_pdf'
        print*,point
        print*,arr
        read*,pause_k
    end if

    
end subroutine
    
function kernel_d(u)
    use nrtype
    implicit none
    real(SP)::kernel_d
    real(SP)::u
    
    !Epanechnikov
    if (abs(u)<=1.0_sp) then
        kernel_d=3.0_sp/4.0_sp*(1.0_sp-u**2.0_sp)
    else
        kernel_d=0.0_sp
    end if
    !Gaussian Kernel
    !kernel_d=1.0_sp/sqrt(2*pi)*exp(-0.5_sp*u**2)
    
end function

 
    