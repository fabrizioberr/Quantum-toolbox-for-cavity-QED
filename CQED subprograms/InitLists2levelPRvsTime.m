% Parameters of the incoming gaussian pulse
omega_pulse_ev=omega_d_ev+detuning_pulse_QD_muev*1e-6; %center energy of the incoming pulse, in eV
omega_pulse=omega_pulse_ev*ev/hbar*1e-12; %in rad/ps

% Parameters for the computation of time evolutions
t_delay = 2*FWHM; % Time at which the pulse is maximally intense, so that the computation starts when the pulse has not arrived yet
t_max_ps=t_delay + 4*FWHM +0.5/gamma_sp; % Final time where we stop the computation and plots of time evolutions
nb_points_time = 1000; % Time resolution/Number of iterations / <100000 otherwise the integrating the master equation gets difficult (odesolve))
t_min = 0*FWHM; % Initial time considered for the computations and plots of time evolutions
t_step=(t_max_ps-t_min)/(nb_points_time-1); % Duration of a time step
t_list=linspace(t_min,t_max_ps,nb_points_time); % list of all the times considered in the computation and plots

% Initialization of qo array of identity operator
Id_vs_time = qo;
for time_index = 1:nb_points_time
    Id_vs_time{time_index} = Id;
end