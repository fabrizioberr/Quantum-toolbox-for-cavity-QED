% Parameters of the incoming gaussian pulse
omega_pulse_ev=omega_d_ev+detuning_pulse_QD_muev*1e-6; %center energy of the incoming pulse, in eV
omega_pulse=omega_pulse_ev*ev/hbar*1e-12; %in rad/ps

% Parameters for the computation of time evolutions
t_step=(t_max_ps-t_min)/(nb_points_time-1); % Duration of a time step
t_list=linspace(t_min,t_max_ps,nb_points_time); % list of all the times considered in the computation and plots

% Parameters for the computation of preliminary time evolution,
% between 0 (long before the pulse) and t_min (time at which we want
% to start plotting and integrating the physical quantities
nb_points_time_before_t_min = 5; % (Low) time resolution for first evolution of the system for initialization
t_list_before_t_min = linspace(0,t_min,nb_points_time_before_t_min); % time array for first evolution of the sistem



% Initialization of qo array of identity operator
Id_vs_time = qo;
for time_index = 1:nb_points_time
    Id_vs_time{time_index} = Id;
end


full_tau_list = [-flip(t_list(2:end)) t_list(1:end)]; %t_list with both negative and positive delays for plotting

%%% Preallocation of the memory to save computing time

% Initialization of g2 vs (t1,t2)
g2_reflected_vs_t1_t2 = zeros(nb_points_time,nb_points_time);
g2_transmitted_vs_t1_t2 = zeros(nb_points_time,nb_points_time);
g2_emitted_vs_t1_t2 = zeros(nb_points_time,nb_points_time);

% Initialization of conditional occupation probabilities vs (t1,t2)
occupation_ground_vs_t1_vs_t2_after_click_b_out_at_t1 = zeros(nb_points_time,nb_points_time);
occupation_excited_vs_t1_vs_t2_after_click_b_out_at_t1 = zeros(nb_points_time,nb_points_time);
occupation_ground_vs_t1_vs_t2_after_click_c_out_at_t1 = zeros(nb_points_time,nb_points_time);
occupation_excited_vs_t1_vs_t2_after_click_c_out_at_t1 = zeros(nb_points_time,nb_points_time);
occupation_ground_vs_t1_vs_t2_after_click_e_out_at_t1 = zeros(nb_points_time,nb_points_time);
occupation_excited_vs_t1_vs_t2_after_click_e_out_at_t1 = zeros(nb_points_time,nb_points_time);

% Initialization of the "zero" density matrix which will be used to fill the 
% conditional density matrices, using 0 values for t2 < t1
zero_density_matrix_vs_t2_before_t1 = qo;

% Initialization of correlated normalized g2(tau)
normalized_g2_vs_delay_emitted = zeros(1,nb_points_time);
normalized_g2_vs_delay_reflected = zeros(1,nb_points_time);
normalized_g2_vs_delay_transmitted = zeros(1,nb_points_time);

% Initialization of uncorrelated normalized g2(tau)
normalized_g2_vs_delay_uncorrelated_emitted = zeros(1,nb_points_time);
normalized_g2_vs_delay_uncorrelated_reflected = zeros(1,nb_points_time);
normalized_g2_vs_delay_uncorrelated_transmitted = zeros(1,nb_points_time);
