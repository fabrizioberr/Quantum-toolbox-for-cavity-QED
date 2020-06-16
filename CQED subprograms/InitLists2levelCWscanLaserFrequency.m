%Parameters for the calculation of spectra
omega_laser_min_ev=omega_d_ev+min_detuning_muev*1e-6;
omega_laser_max_ev=omega_d_ev+max_detuning_muev*1e-6;

%%%% Initialization of lists to calculate and plot the spectra as a function of omega_laser
omega_laser_list_ev = linspace(omega_laser_min_ev,omega_laser_max_ev,nb_points_spectrum);% list of laser photon energies for the plots in eV
omega_laser_list = omega_laser_list_ev*ev/hbar*1e-12;% list of laser angular frequencies in rad/ps, the unit used for calculations
omega_step = (max(omega_laser_list) - min(omega_laser_list)) / (nb_points_spectrum-1); %in rad/ps, step for the calculation of integrals

% List to plot the spectra as a function of the detuning omega_laser-omega_d, in  mueV
detuning_list_muev = (omega_laser_list_ev - omega_d_ev)*1e6; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Memory preallocation (allows gaining in calculation time)
% The  "_vs_omega" indicates here that this is a list of values related to the defferent_values of omega_laser

total_flux_reflected_photons_vs_omega = zeros(1,length(omega_laser_list)); %flux in ps^(-1)
total_flux_transmitted_photons_vs_omega = zeros(1,length(omega_laser_list)); %flux in ps^(-1)
total_flux_diffracted_photons_vs_omega = zeros(1,length(omega_laser_list)); %flux in ps^(-1)
total_flux_emitted_photons_vs_omega = zeros(1,length(omega_laser_list)); %flux in ps^(-1)

occupation_excited_state_vs_omega = zeros(1,length(omega_laser_list)); % for the occupation of the excited state 
occupation_ground_state_vs_omega = zeros(1,length(omega_laser_list)); % for the occupation of the ground state


%%%% OPTIONAL : Memory preallocation for the coherent part of the output fields
% Here the term "laser_coherent" means that it corresponds to the part of the flux that is 
% coherent with the incoming excitation laser, and not the "total" flux. Obviously the incoherent part
% is just given by substrating the "laser_coherent" part from the "total" flux

flux_reflected_photons_laser_coherent_vs_omega  = zeros(1,length(omega_laser_list)); %flux in ps^(-1)
flux_transmitted_photons_laser_coherent_vs_omega = zeros(1,length(omega_laser_list)); %flux en s^(-1)
flux_diffracted_photons_laser_coherent_vs_omega = zeros(1,length(omega_laser_list)); %flux en s^(-1)
flux_emitted_photons_laser_coherent_vs_omega  = zeros(1,length(omega_laser_list)); %flux en s^(-1)
