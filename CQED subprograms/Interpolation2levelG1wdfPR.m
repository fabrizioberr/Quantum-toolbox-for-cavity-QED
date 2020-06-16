%% Memory preallocation (allows gaining in calculation time)
interpolated_G1_reflected_vs_time_tau = 0*time_list'*tau_list;
interpolated_reflected_coherent_vs_time_tau = 0*time_list'*tau_list;
interpolated_G1_transmitted_vs_time_tau = 0*time_list'*tau_list;
interpolated_transmitted_coherent_vs_time_tau = 0*time_list'*tau_list;
interpolated_G1_emitted_vs_time_tau = 0*time_list'*tau_list;
interpolated_emitted_coherent_vs_time_tau = 0*time_list'*tau_list;
interpolated_G1_diffracted_vs_time_tau = 0*time_list'*tau_list;
interpolated_diffracted_coherent_vs_time_tau = 0*time_list'*tau_list;

%From the definition of time = (t1+t2)/2 and tau=t2-t1, it is found that for a
%matrix tau_index = t2_index-t1_index + N and time_index =
%t1_index+t2_index-1. Therefore, it is obtained that t1_index =
%(time_index-tau_index+1+1)/2 and t2_index = (time_index+tau_index+1-N)/2. 

%NB: The following interpolation is based on a simple arithmetic mean between
% consecutives values at half-integer values. "Fancier" interpolations
% would surely increase the overall accuracy.

for time_index = 1:length(time_list)
    for tau_index = 1:length(tau_list)
        t1_index = (time_index-((tau_index+1)/2)+1+nb_points_time)/2; %obs: tau_index from the previous formula is acqually "(tau_index+1)/2" since we have set tau_step = t_step/2
        t2_index = (time_index+((tau_index+1)/2)+1-nb_points_time)/2; %same as above
        
        if t1_index <= nb_points_time && t1_index >= 1 && t2_index <= nb_points_time && t2_index >= 1
            
            if floor(t1_index)==t1_index && floor(t2_index)==t2_index
                interpolated_G1_reflected_vs_time_tau(time_index, tau_index) = G1_reflected_vs_t1_t2(t1_index, t2_index);
                interpolated_reflected_coherent_vs_time_tau(time_index, tau_index) = expect_b_out_dag_t1_times_expect_b_out_t2(t1_index, t2_index);
                interpolated_G1_transmitted_vs_time_tau(time_index, tau_index) = G1_transmitted_vs_t1_t2(t1_index, t2_index);
                interpolated_transmitted_coherent_vs_time_tau(time_index, tau_index) = expect_c_out_dag_t1_times_expect_c_out_t2(t1_index, t2_index);
                interpolated_G1_emitted_vs_time_tau(time_index, tau_index) = G1_emitted_vs_t1_t2(t1_index, t2_index);
                interpolated_emitted_coherent_vs_time_tau(time_index, tau_index) = expect_e_out_dag_t1_times_expect_e_out_t2(t1_index, t2_index);
                interpolated_G1_diffracted_vs_time_tau(time_index, tau_index) = G1_diffracted_vs_t1_t2(t1_index, t2_index);
                interpolated_diffracted_coherent_vs_time_tau(time_index, tau_index) = expect_d_out_dag_t1_times_expect_d_out_t2(t1_index, t2_index);
            
            elseif floor(t1_index)==t1_index && floor(t2_index)~= t2_index
                interpolated_G1_reflected_vs_time_tau(time_index, tau_index) = (G1_reflected_vs_t1_t2(t1_index, floor(t2_index)) + G1_reflected_vs_t1_t2(t1_index, ceil(t2_index)))/2;
                interpolated_reflected_coherent_vs_time_tau(time_index, tau_index) = (expect_b_out_dag_t1_times_expect_b_out_t2(t1_index, floor(t2_index)) +expect_b_out_dag_t1_times_expect_b_out_t2(t1_index, ceil(t2_index)))/2;
                interpolated_G1_transmitted_vs_time_tau(time_index, tau_index) = (G1_transmitted_vs_t1_t2(t1_index, floor(t2_index)) + G1_transmitted_vs_t1_t2(t1_index, ceil(t2_index)))/2;
                interpolated_transmitted_coherent_vs_time_tau(time_index, tau_index) = (expect_c_out_dag_t1_times_expect_c_out_t2(t1_index, floor(t2_index)) +expect_c_out_dag_t1_times_expect_c_out_t2(t1_index, ceil(t2_index)))/2;
                interpolated_G1_emitted_vs_time_tau(time_index, tau_index) = (G1_emitted_vs_t1_t2(t1_index, floor(t2_index)) + G1_emitted_vs_t1_t2(t1_index, ceil(t2_index)))/2;
                interpolated_emitted_coherent_vs_time_tau(time_index, tau_index) = (expect_e_out_dag_t1_times_expect_e_out_t2(t1_index, floor(t2_index)) +expect_e_out_dag_t1_times_expect_e_out_t2(t1_index, ceil(t2_index)))/2;
                interpolated_G1_diffracted_vs_time_tau(time_index, tau_index) = (G1_diffracted_vs_t1_t2(t1_index, floor(t2_index)) + G1_diffracted_vs_t1_t2(t1_index, ceil(t2_index)))/2;
                interpolated_diffracted_coherent_vs_time_tau(time_index, tau_index) = (expect_d_out_dag_t1_times_expect_d_out_t2(t1_index, floor(t2_index)) +expect_d_out_dag_t1_times_expect_d_out_t2(t1_index, ceil(t2_index)))/2;            
            
            elseif floor(t2_index)==t2_index && floor(t1_index)~= t1_index
                interpolated_G1_reflected_vs_time_tau(time_index, tau_index) = (G1_reflected_vs_t1_t2(floor(t1_index), t2_index)+G1_reflected_vs_t1_t2(ceil(t1_index), t2_index))/2;
                interpolated_reflected_coherent_vs_time_tau(time_index, tau_index) = (expect_b_out_dag_t1_times_expect_b_out_t2(floor(t1_index), t2_index)+expect_b_out_dag_t1_times_expect_b_out_t2(ceil(t1_index), t2_index))/2;
                interpolated_G1_transmitted_vs_time_tau(time_index, tau_index) = (G1_transmitted_vs_t1_t2(floor(t1_index), t2_index)+G1_transmitted_vs_t1_t2(ceil(t1_index), t2_index))/2;
                interpolated_transmitted_coherent_vs_time_tau(time_index, tau_index) = (expect_c_out_dag_t1_times_expect_c_out_t2(floor(t1_index), t2_index)+expect_c_out_dag_t1_times_expect_c_out_t2(ceil(t1_index), t2_index))/2;
                interpolated_G1_emitted_vs_time_tau(time_index, tau_index) = (G1_emitted_vs_t1_t2(floor(t1_index), t2_index)+G1_emitted_vs_t1_t2(ceil(t1_index), t2_index))/2;
                interpolated_emitted_coherent_vs_time_tau(time_index, tau_index) = (expect_e_out_dag_t1_times_expect_e_out_t2(floor(t1_index), t2_index)+expect_e_out_dag_t1_times_expect_e_out_t2(ceil(t1_index), t2_index))/2;
                interpolated_G1_diffracted_vs_time_tau(time_index, tau_index) = (G1_diffracted_vs_t1_t2(floor(t1_index), t2_index)+G1_diffracted_vs_t1_t2(ceil(t1_index), t2_index))/2;
                interpolated_diffracted_coherent_vs_time_tau(time_index, tau_index) = (expect_b_out_dag_t1_times_expect_d_out_t2(floor(t1_index), t2_index)+expect_b_out_dag_t1_times_expect_d_out_t2(ceil(t1_index), t2_index))/2;           
            else
                interpolated_G1_reflected_vs_time_tau(time_index, tau_index) = (G1_reflected_vs_t1_t2(floor(t1_index), floor(t2_index))+G1_reflected_vs_t1_t2(ceil(t1_index), floor(t2_index))+G1_reflected_vs_t1_t2(floor(t1_index), ceil(t2_index))+G1_reflected_vs_t1_t2(ceil(t1_index), ceil(t2_index)))/4;
                interpolated_reflected_coherent_vs_time_tau(time_index, tau_index) = (expect_b_out_dag_t1_times_expect_b_out_t2(floor(t1_index), floor(t2_index))+expect_b_out_dag_t1_times_expect_b_out_t2(ceil(t1_index), floor(t2_index)) + expect_b_out_dag_t1_times_expect_b_out_t2(floor(t1_index), ceil(t2_index))+expect_b_out_dag_t1_times_expect_b_out_t2(ceil(t1_index), ceil(t2_index)))/4;
                interpolated_G1_transmitted_vs_time_tau(time_index, tau_index) = (G1_transmitted_vs_t1_t2(floor(t1_index), floor(t2_index))+G1_transmitted_vs_t1_t2(ceil(t1_index), floor(t2_index))+G1_transmitted_vs_t1_t2(floor(t1_index), ceil(t2_index))+G1_transmitted_vs_t1_t2(ceil(t1_index), ceil(t2_index)))/4;
                interpolated_transmitted_coherent_vs_time_tau(time_index, tau_index) = (expect_c_out_dag_t1_times_expect_c_out_t2(floor(t1_index), floor(t2_index))+expect_c_out_dag_t1_times_expect_c_out_t2(ceil(t1_index), floor(t2_index)) + expect_c_out_dag_t1_times_expect_c_out_t2(floor(t1_index), ceil(t2_index))+expect_c_out_dag_t1_times_expect_c_out_t2(ceil(t1_index), ceil(t2_index)))/4;
                interpolated_G1_emitted_vs_time_tau(time_index, tau_index) = (G1_emitted_vs_t1_t2(floor(t1_index), floor(t2_index))+G1_emitted_vs_t1_t2(ceil(t1_index), floor(t2_index))+G1_emitted_vs_t1_t2(floor(t1_index), ceil(t2_index))+G1_emitted_vs_t1_t2(ceil(t1_index), ceil(t2_index)))/4;
                interpolated_emitted_coherent_vs_time_tau(time_index, tau_index) = (expect_e_out_dag_t1_times_expect_e_out_t2(floor(t1_index), floor(t2_index))+expect_e_out_dag_t1_times_expect_e_out_t2(ceil(t1_index), floor(t2_index)) + expect_e_out_dag_t1_times_expect_e_out_t2(floor(t1_index), ceil(t2_index))+expect_e_out_dag_t1_times_expect_e_out_t2(ceil(t1_index), ceil(t2_index)))/4;
                interpolated_G1_diffracted_vs_time_tau(time_index, tau_index) = (G1_diffracted_vs_t1_t2(floor(t1_index), floor(t2_index))+G1_diffracted_vs_t1_t2(ceil(t1_index), floor(t2_index))+G1_diffracted_vs_t1_t2(floor(t1_index), ceil(t2_index))+G1_diffracted_vs_t1_t2(ceil(t1_index), ceil(t2_index)))/4;
                interpolated_diffracted_coherent_vs_time_tau(time_index, tau_index) = (expect_d_out_dag_t1_times_expect_d_out_t2(floor(t1_index), floor(t2_index))+ expect_d_out_dag_t1_times_expect_d_out_t2(ceil(t1_index), floor(t2_index)) + expect_d_out_dag_t1_times_expect_d_out_t2(floor(t1_index), ceil(t2_index))+expect_d_out_dag_t1_times_expect_d_out_t2(ceil(t1_index), ceil(t2_index)))/4;           
            end
        end
    end
end