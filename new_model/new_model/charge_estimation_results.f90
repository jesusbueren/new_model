subroutine charge_estimation_results(H_av,p_av,dist_init_av)
    use dimensions; use nrtype
    implicit none
    REAL(SP),dimension(variables,clusters),intent(out)::p_av
    REAL(SP),dimension(clusters+1,clusters+1,generations,L_PI,L_gender),intent(out)::H_av
    REAL(SP),dimension(clusters,generations,L_PI,L_gender),intent(out)::dist_init_av
    integer,parameter::iterations=149712,end_it=65000,burn=50000
    REAL(SP),dimension(variables,clusters,iterations)::p
    REAL(SP),dimension(clusters*(clusters+1)*generations*L_PI*L_gender)::health_tr_vec
    REAL(SP),dimension(clusters,clusters,generations,L_PI,L_gender)::H_new
    integer::c_l,g_l,e_l,ge_l,v_l,c_l2,it
    character(LEN=55)::path="C:\Users\jbueren\Google Drive\JMP\Code\Health Dynamics\"
    !character(LEN=53)::path="C:\Users\Jesus\Google Drive\JMP\Code\Health Dynamics\"
    character(LEN=8)::path_s="Results\"
    character(LEN=1)::s_c
    
    Write( s_c, '(I1)' )  clusters
    open(unit=9,file=path//path_s//'p'//s_c//'.txt')
        read(9,'(F20.5)') p
    close(9)
    
    !Mean pr of i-adl by cluster
    do v_l=1,variables; do c_l=1,clusters
        p_av(v_l,c_l)=sum(p(v_l,c_l,burn:end_it))/dble(end_it-burn)
    end do;end do

    !Mean transition probabilitities (this file comes from convergence_check.m)
    open(unit=10,file=path//path_s//"transitions.txt")
        read(10,*) health_tr_vec
    close(10)
    H_av=0.0_sp
    H_av(1:clusters,:,:,:,:)=reshape(health_tr_vec,(/clusters,clusters+1,generations,L_PI,2/),order=(/2,3,4,5,1/))
    !initial distribution for mean probabilities
    do c_l=1,clusters ; do g_l=1,generations;  do e_l=1,L_PI; do ge_l=1,2
        do c_l2=1,clusters
            H_new(c_l,c_l2,g_l,e_l,ge_l)=H_av(c_l,c_l2,g_l,e_l,ge_l)/sum(H_av(c_l,1:clusters,g_l,e_l,ge_l))
        end do
        H_av(c_l,:,g_l,e_l,ge_l)=H_av(c_l,:,g_l,e_l,ge_l)/sum(H_av(c_l,:,g_l,e_l,ge_l))
    end do; end do; end do; end do 
      
    dist_init_av=1/dble(clusters)
    do g_l=1,generations; do e_l=1,L_PI; do ge_l=1,2
        if (g_l==1) then
            do it=1,40
                dist_init_av(:,g_l,e_l,ge_l)=matmul(dist_init_av(:,g_l,e_l,ge_l),H_new(:,:,g_l,e_l,ge_l))     
            end do
        else
            dist_init_av(:,g_l,e_l,ge_l)=matmul(dist_init_av(:,g_l-1,e_l,ge_l),H_new(:,:,g_l,e_l,ge_l))
        end if 
        dist_init_av(:,g_l,e_l,ge_l)=dist_init_av(:,g_l,e_l,ge_l)/sum(dist_init_av(:,g_l,e_l,ge_l))
    end do; end do; end do
    H_av(:,1:clusters,generations,:,:)=0.0_sp
    H_av(:,clusters+1,generations,:,:)=1.0_sp
end subroutine