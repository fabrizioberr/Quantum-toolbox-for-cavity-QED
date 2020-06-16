clear
clc
close all

%% Important note: these paths must be modified if needed
addpath(genpath('..\QotoolboxV015'));
addpath(genpath('..\CQED subprograms'));
addpath(genpath('..\CQED device parameters'))
savepath

% In addition, for the mesolve function to operate the executable files
% (.exe) and batch files (.bat) contained in  '[...] \QotoolboxV015\bin'
% have to be copied to  a folder that is on the Windows system path, in the
% main hard  drive where Windows is installed. This can be for example in:
% 'C:\Program Files\Matlab\R2014a\bin'.

% Warning: for the adiabatic version to converge, the tolerance in
% mesolve.m function must be reduced compared to the defaut values. For
% example:
% ode2file('ode_input.dat',L,rho0,t_list,struct('reltol',7e-8,'abstol',8e-8));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  g1SDPR : g1(t1,t2) and spectral densities in PR  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%This section labeled "g1SDPR" studies the resonant resonance fluorescence
%spectra in pulsed regime, by prior evaluation of the normalized
%correlation function g1(t1,t2). By setting t=(t1-t2)/2 and tau = t2-t1, an
%interpolated (unnormalized) G1(t,tau) is computed from G1(t1,t2). The
%Fourier Transform is evaluated over tau to find the associated
%Wigner-Ville-Distribution WDF(t,omega). To check its normalization,
%WDF(t,omega) is integrated over omega and compared with the photon fluxes.
%WDF is integrated over t for determining the spectral energy density
%ESD(omega), which is what is measured by the spectrometer. It is verified
%that the area of this spectrum gives the total number of photons. At last,
%the spectra of the total and the coherent component of the fluxes are
%compared for the reflected, transmitted, diffracted and emitted fields.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

%% Choice of full model 'F' or adiabatic model 'A'
model = 'A'; 

%% Experimental conditions (to be edited)
detuning_QD_C_muev = 0; %Detuning between the QD and cavity frequencies, in mueV
detuning_pulse_QD_muev = 0; %Detuning between the pulse central frequency and the QD frequency, in mueV
eta_in = 1; % Injection efficiency for the incoming photons (depends on experimentally-achieved spatial coupling)
Nb_photons_pulse = 1;  % Average number of incoming photons in a pulse. This quantity should be multiplied by eta_in to know the number of incoming photons actually coupled to the optical mode
FWHM_pulse = 15; %in ps, full width at half-maximum of the incoming Gaussian pulse intensity (unit: ps since angular frequencies are in rad/ps)
% Parameters for the evaluation of the temporal evolution
% 
% NB1: the maximum delay "t_max_ps" will also dictate the frequency
% resolution of the spectra, given by the angular frequency step
% "omega_step". The corresponding angular frequency lists are defined in
% the "Init_lists_..." subprogram (see also below details on the
% calculation and Fast Fourier Transform (FFT) algorithm)
%
% NB2: one should be careful that "t_max_ps" is large enough to include a
% good approximation of "infinite delays" (check that the g1(tau) function
% has had enough time to truly converge), while keeping a number of points
% large enough to ensure a good temporal resolution. This is especially
% important for high input powers where artifacts can appear.
t_delay = 2*FWHM_pulse; % Time at which the pulse is maximally intense, so that the computation starts when the pulse has not arrived yet
t_max_ps = (t_delay + 4*FWHM_pulse)*5; % Final time where we stop the computation and plots of time evolutions
nb_points_time = 2^7; % Time resolution/Number of iterations / <100000 otherwise the integrating the master equation gets difficult (odesolve))
% --> This should be such that 4*nb_points_time-3 is less than a power of 2
% if one wants to reduce the amount of zero padding, the latter being
% required for FFT optimization.
t_min = 0.1*FWHM_pulse; % Initial time considered for the computations and plots of time evolutions
width_spectral_window_muev = 300; % width of the spectral window to be displayed, centered on omega_pulse
%% Initialization of parameters, operators, arrays, etc...
Init2levelDeviceParametersOngoingTest;
Init2levelHilbertSpaceAndOperators;
InitMaps2levelG1wdfPR;
%%
% Definition of the input field in ps^(-1/2), in the form of a fseries (necessary for integrating the master equation)
Standard_deviation_b_in_PR = FWHM_pulse/(2*sqrt(log(2))); %Deduced from the properties of a Gaussian function
b_in_fn = fn('gauss',t_delay,Standard_deviation_b_in_PR) * sqrt( eta_in*Nb_photons_pulse / ( sqrt(pi) * Standard_deviation_b_in_PR ) ); % square root of the incoming photon number per time unit, in ps-1/2
b_in_vs_time = fsval(b_in_fn,t_list); %scalar array representing b_in vs time

% Initial density matrix before the pulse has started
switch model
    case 'F' %Full model
        psi0 = tensor(Vacuum_state,g_ket); % Initial state: tensorial product of photonic vacuum and QD ground state
        rho0 = psi0*psi0'; % Density matrix corresponding to the initial pure state
    case 'A' % Adiabatic model
        rho0 = g_ket*g_ket'; % Density matrix corresponding to the initial pure state
end

%% %%%%%%%% System Hamiltonian and time-dependent operators  %%%%%%%%%%%%
%
% The system Hamiltonian is time-dependent due to the function b_in_fn
% describing the input field b_in(t). 

% In addition, in the case of adiabatic elimination of the cavity mode an effective
% operator "a" is defined, acting on the QD subspace, based on the formula for
% adiabatic elimination. Since this formula depends
% on b_in(t), we define a time-dependent quantity "a_vs_time", which is an
% array containing, for each time of t_list, the corresponding operator
% "a". 
% In the full model case, to simplify the following calculations, we define the 
% same quantity a_vs_time, yet this time this array contains the same operator 
% (annihilation operator "a" acting on the cavity subspace), replicated for all
% times of t_list.

switch model 
        case 'A' % Adiabatic model
           
            Delta = 2*(omega_pulse-omega_c)/kappa; %normalized pulse detuning            
            
            %%%%%%%%% Definition of the adiabatic-model Hamiltonian (in the frame rotating at angular frequency omega_pulse)
            H_PR = (omega_eff-omega_pulse)*sigma_dag*sigma - 1i*sqrt(Gamma_0*eta_top)*...
                    ((1-1i*Delta)^(-1)*b_in_fn*sigma_dag-(1+1i*Delta)^(-1)*b_in_fn'*sigma); % Hamiltonian 
            
            %%% Definition of the time-dependent operator "a_vs_time" acting in the QD subspace 
            a_vs_time = -2*(kappa*(1-1i*Delta_QDC))^(-1)*g*sigma*Id_vs_time-2*sqrt(kappa_top)*(kappa*(1-1i*Delta))^(-1)*fsval(b_in_fn,t_list).*Id_vs_time; %annihilation operator a in adiabatic approximation, as fseries
            
        case 'F'  % Full model
            
            %%%%%%%%% Definition of the full-model Hamiltonian (in the frame rotating at angular frequency omega_pulse)
            H_PR = (omega_d-omega_pulse)*sigma_dag*sigma...
                 + (omega_c-omega_pulse)*a_dag*a...
                 + 1i*g*(sigma_dag*a-a_dag*sigma)...
                 - 1i*sqrt(kappa_top)*b_in_fn*(a_dag-a); 
             
            % Definition of the operator a_vs_time, even though "a" is constant in the full model, to reduce the number of "switch" in the following code
             a_vs_time = a*Id_vs_time;
             
end % end of the "switch model"

%Superoperator associated to the coherent processes (Hamiltonian)
L_coh = -1i * (spre(H_PR) - spost(H_PR)); 

%%%%%%%%% Calculation of the Liouvillian superoperator
Liouvillian  = L_coh + L_incoh; % Total Liouvillian superoperator including both coherent processes (Hamiltonian) and incoherent processes (dissipative jumps)

%%%%%%%%%%% Numerical Integration of the Master Equation %%%%%%%%

% Computation of preliminary time evolution, between 0 (long before the pulse)
% and t_min (time at which we want to start plotting and integrating the physical
% quantities. Such a computation can be performed with very low time resolution, 
% i.e. the corresponding t_list_before_t_min has a very low number of points,
% since the density matrix almost doesn't evolve between 0 and t_min. 
%
% NB: such a preliminary time evolution is mandatory to avoid having strictly
% zero expectation values for some quantities (such as the input or output fields),
% during the following time evolution between t_min and t_max). Indeed, zero
% expectation values lead to NaN errors when used in normalizing physical
% quantities, such as conditional density matrices (see below), or
% Stokes/Bloch coordinates in the Poincare'/Bloch sphere.

density_matrix_vs_time_before_t_min = mesolve(Liouvillian,rho0,t_list_before_t_min); %master equation solver based on odesolve: first evolution of the system
density_matrix_at_t_min = density_matrix_vs_time_before_t_min{nb_points_time_before_t_min};

% Computation of the density matrix vs time with t_list, i.e. between t_min
% and t_max, requiring a large enough time resolution.
density_matrix_vs_time = mesolve(Liouvillian,density_matrix_at_t_min,t_list); % second evolution of the system

%%%%%%%%%%%%%%%%%% Definition of the output operators %%%%%%%%%%%%%%%%%%%%%
% These are general formulas for both the adiabatic and full model, depending on "a_vs_time"

b_out_vs_time = b_in_vs_time*Id_vs_time + sqrt(kappa_top)*a_vs_time; % definition of the operator b_out, i.e. the output operator for the reflected light, in ps^(-1/2)
c_out_vs_time = sqrt(kappa_bottom)*a_vs_time; % definition of the operator c_out_vs_time, i.e. the output operator for the transmitted light, in ps^(-1/2)
d_out_vs_time = sqrt(kappa_loss)*a_vs_time; % UNUSED HERE: definition of the operator d_out_vs_time, i.e. the output operator for the diffracted/lost light, in ps^(-1/2)
% e_out = sqrt(gamma_sp)*sigma %  output operator for the light spontaneously emitted outside the cavity mode, in ps^(-1/2), already defined in Init_2level_Hilbert_space_and_operators.m

% Calculation of the total photon flux as a function of time, for the various optical fields 
flux_injected_photons_vs_time = b_in_vs_time.^2; % total flux of injected photons taking into account  eta_in (so only the photons coupled to the cavity mode), in ps(-1)
flux_reflected_photons_vs_time = real(expect(b_out_vs_time'*b_out_vs_time,density_matrix_vs_time)); %flux in ps^(-1)
flux_transmitted_photons_vs_time = real(expect(c_out_vs_time'*c_out_vs_time,density_matrix_vs_time)); %flux in ps^(-1)
flux_diffracted_photons_vs_time = real(expect(d_out_vs_time'*d_out_vs_time,density_matrix_vs_time)); %flux in ps^(-1)
flux_emitted_photons_vs_time = real(expect(e_out'*e_out,density_matrix_vs_time)); %flux in ps^(-1)

% Calculation of the total coherent photon flux as a function of time, for the various optical fields 
flux_reflected_photons_pulse_coherent_vs_time= abs(expect(b_out_vs_time,density_matrix_vs_time)).^2; % coherent flux = <b_out'(t)><b_out(t)> in ps^{-1}
flux_transmitted_photons_pulse_coherent_vs_time=abs(expect(c_out_vs_time,density_matrix_vs_time)).^2; % coherent flux = <c_out'(t1)><c_out(t2)> in ps^{-1}
flux_diffracted_photons_pulse_coherent_vs_time=abs(expect(d_out_vs_time,density_matrix_vs_time)).^2; % coherent flux = <d_out'(t1)><d_out(t2)> in ps^{-1}
flux_emitted_photons_pulse_coherent_vs_time=abs(expect(e_out,density_matrix_vs_time)).^2; % coherent flux = <e_out'(t1)><e_out(t2)> in ps^{-1}

%Computing the number of reflected/transmitted/emitted photons, by integrating over all times in t_list
Nb_reflected_photons = trapz(t_list,flux_reflected_photons_vs_time);
Nb_transmitted_photons = trapz(t_list,flux_transmitted_photons_vs_time);
Nb_diffracted_photons = trapz(t_list,flux_diffracted_photons_vs_time);
Nb_emitted_photons = trapz(t_list,flux_emitted_photons_vs_time);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Evaluation of g1(t1,t2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For a given output field operator b, the two-time first-order auto-correlation function, 
% g1(t1,t2), also called the "degree of first-order coherence", is defined by:
%      g1(t1,t2)= < b'(t1) b(t2) > / sqrt[ <b'(t1)> <b(t2)> ]
% % 
% To calculate such a quantity, one can use the Heisenberg representation of a two-time
% correlation function < A(t1) B(t2) >, where operators A and B are time-dependent, and
% the Schroedinger approach that we have to use here  (where the density matrix varies).
% The "recipe" to deduce < A(t1)B(t2) > consists in:
%     - Letting the density matrix evolve from time 0 to t1
%     - Replacing rho(t1) by a fictitious density matrix rho(t1)*A
%     - Computing the evolution of this fictitious density matrix between times t1 and t2
%     - Calculating the expectation value of operator B using this fictitious density matrix
%
% Note that even if rho(t1)*A is not a real/valid density matrix (it's not even Hermitian),
% we can at least make it of the order of unity, to ensure an optimal numerical convergence
% (especially important in the pulsed regime where for example the operator b_out is
% extremely small at the beginning of the pulse). Looking at the definition of g1(t1,t2),
% we see  that this is readily obtained by taking the operator B as b/sqrt( <b'b>), with
% the consequence that operator A has to be taken equal to b'/sqrt( <b'b> ). 

% NB1: The approach used in the "g2CW_vs_delay" and "g2PR_vs_t1_t2" programs, to calculate 
% second-order autocorrelation functions, is more focused on the physical/experimental
% definition of these quantities. It also makes use of real/normalized/valid density matrices,
% contrary to the fictitious density matrices used below. But the theory behind is also 
% entirely linked to the use of the Quantum Regression Theorem.

%%%%%%%%%%%  Calculation of fictitious density matrix rho(t1)*A  %%%%%%%%%%
% In the following the "density_matrix_times_OPERATOR_dag_vs_t1"-s are
% defined for each value of time t1 in t_list.

% NB: notice the normalization  b'/sqrt( <b'b> ) to optimize numerical
% convergence
density_matrix_times_b_out_dag_vs_t1 = density_matrix_vs_time*b_out_vs_time'./sqrt(flux_reflected_photons_vs_time); 
density_matrix_times_c_out_dag_vs_t1 = density_matrix_vs_time*c_out_vs_time'./sqrt(flux_transmitted_photons_vs_time);
density_matrix_times_d_out_dag_vs_t1 = density_matrix_vs_time*d_out_vs_time'./sqrt(flux_diffracted_photons_vs_time); 
density_matrix_times_e_out_dag_vs_t1 = density_matrix_vs_time*e_out'./sqrt(flux_emitted_photons_vs_time);

% NB: "tic" is used as a "start" time for the measurement of the computing time between "tic" and "toc"
%tic

%Cycle over all times t1 in t_list, corresponding to the moment where a first click occurred
for t1_index = 1:nb_points_time 
    
    % Both t1 and t2 are values in t_list. However, to compute the
    % "fictitious density matrix B*rho(t1)" vs t2, we consider
    % only t2 >= t1, and we need to  deal with the fact that there are less
    % and less remaining values of t2 in t_list, when t1 increases. To keep
    % all quantities defined in the full t_list, for each time t1 < t2 we
    % have a "zero" density matrix, i.e. a fictitious density matrix with
    % only zero elements-

    % Incrementing the array of zero density matrices, each time t1_index
    % is increased, to fill the density matrix for times t1 < t2
    
    if (t1_index >=2) % No need to include a zero density matrix at the first value of t1
        zero_density_matrix_vs_t2_before_t1{t1_index-1} = 0*Id; %null density matrix with the same dimensions of the involved Hilbert space
    end
    
    % Array of normalized density matrices, conditioned on the detection of a
    % click at time t1, as a function of t2 >= t1 (and zero otherwise)
    
    % NB: notice the normalization  b'/sqrt( <b'b> ) to optimize numerical
    % convergence
    density_matrix_times_b_out_dag_vs_t2 = [zero_density_matrix_vs_t2_before_t1  mesolve(Liouvillian,density_matrix_times_b_out_dag_vs_t1{t1_index},t_list(t1_index:end))]./sqrt(flux_reflected_photons_vs_time);
    density_matrix_times_c_out_dag_vs_t2 = [zero_density_matrix_vs_t2_before_t1  mesolve(Liouvillian,density_matrix_times_c_out_dag_vs_t1{t1_index},t_list(t1_index:end))]./sqrt(flux_transmitted_photons_vs_time);
    density_matrix_times_d_out_dag_vs_t2 = [zero_density_matrix_vs_t2_before_t1  mesolve(Liouvillian,density_matrix_times_d_out_dag_vs_t1{t1_index},t_list(t1_index:end))]./sqrt(flux_diffracted_photons_vs_time);
    density_matrix_times_e_out_dag_vs_t2 = [zero_density_matrix_vs_t2_before_t1  mesolve(Liouvillian,density_matrix_times_e_out_dag_vs_t1{t1_index},t_list(t1_index:end))]./sqrt(flux_emitted_photons_vs_time); 
    
    
    % Evaluation of g1(t1,t2) as described in the general notes above, for t2 >= t1 (and zero otherwise)
    % NB: the normalization by the photon flux is already considered in density_matrix_times_OPERATOR_dag_vs_t2
    g1_reflected_vs_t1_t2(t1_index,:) = expect(b_out_vs_time,density_matrix_times_b_out_dag_vs_t2); 
    g1_transmitted_vs_t1_t2(t1_index,:) = expect(c_out_vs_time,density_matrix_times_c_out_dag_vs_t2); 
    g1_diffracted_vs_t1_t2(t1_index,:) = expect(d_out_vs_time,density_matrix_times_d_out_dag_vs_t2); 
    g1_emitted_vs_t1_t2(t1_index,:) = expect(e_out,density_matrix_times_e_out_dag_vs_t2); 
end
%%
%%%%% Completion of previous partially-calculated maps vs t1,t2, to include the case where t2 < t1  

% Due to the symmetry between t1 and t2, we have to use properties like g1(t1,t2)=g1(t2,t1)* to fill
% the voids in the quantities that we have only partially calculated yet
% (since we systematically considered a zero value when t2 < t1. For each
% map function of t1 and t2, this completion is obtained by adding it to
% its transpose conjugate (to replace the zeros at t2 < t1) and dividing by
% 2 the elements along the diagonal (to avoid counting twice the case where
% t2=t1). This is done via an ad-hoc matrix idx below:

idx = ones(nb_points_time)-0.5*diag(ones(nb_points_time,1)); % to divide by 2 the elements along the diagonal

% Completed g1(t1,t2) for the various optical fields

g1_reflected_vs_t1_t2 = (g1_reflected_vs_t1_t2+g1_reflected_vs_t1_t2').*idx;
g1_transmitted_vs_t1_t2 = (g1_transmitted_vs_t1_t2+g1_transmitted_vs_t1_t2').*idx;
g1_diffracted_vs_t1_t2 = (g1_diffracted_vs_t1_t2+g1_diffracted_vs_t1_t2').*idx;
g1_emitted_vs_t1_t2 = (g1_emitted_vs_t1_t2+g1_emitted_vs_t1_t2').*idx; 

% G1 correlation vs (t1,t2), e.g. <b_out'(t1) b_out(t2)>
G1_reflected_vs_t1_t2 = g1_reflected_vs_t1_t2.*(flux_reflected_photons_vs_time'*flux_reflected_photons_vs_time).^0.5;
G1_transmitted_vs_t1_t2 = g1_transmitted_vs_t1_t2.*(flux_transmitted_photons_vs_time'*flux_transmitted_photons_vs_time).^0.5;
G1_diffracted_vs_t1_t2 = g1_diffracted_vs_t1_t2.*(flux_diffracted_photons_vs_time'*flux_diffracted_photons_vs_time).^0.5;
G1_emitted_vs_t1_t2 = g1_emitted_vs_t1_t2.*(flux_emitted_photons_vs_time'*flux_emitted_photons_vs_time).^0.5;

% coherent component cross product, e.g <b_out'(t1)> <b_out(t2)>
expect_b_out_dag_t1_times_expect_b_out_t2 =  expect(b_out_vs_time',density_matrix_vs_time).'*expect(b_out_vs_time,density_matrix_vs_time);
expect_c_out_dag_t1_times_expect_c_out_t2 =  expect(c_out_vs_time',density_matrix_vs_time).'*expect(c_out_vs_time,density_matrix_vs_time);
expect_d_out_dag_t1_times_expect_d_out_t2 =  expect(d_out_vs_time',density_matrix_vs_time).'*expect(d_out_vs_time,density_matrix_vs_time);
expect_e_out_dag_t1_times_expect_e_out_t2 =  expect(e_out',density_matrix_vs_time).'*expect(e_out,density_matrix_vs_time);

%% %%%%%%%%%%%%%%% Wigner distribution function%%%%%%%%%%%%%%%%%%%%%%%%%
%The Wigner distribution function (WDF) is used in signal processing as a
%transform in time-frequency analysis. the Wigner distribution function
%provides the highest possible temporal vs frequency resolution which is
%mathematically possible within the limitations of uncertainty in quantum
%wave theory. Given a non-stationary autocorrelation function C(t1,t2), by
%defining t=(t1+t2)/2 and tau= t2-t1, Fourier transforming the lag is
%obtained :
%WDF(t,omega) = \int C(t+tau/2,t-tau/2)*exp(-i*tau*omega)) d\tau. 
%More info can be found at https://en.wikipedia.org/wiki/Wigner_distribution_function
%Notice that tau is defined based on the definition of C(t1,t2) in this
%script, with opposite sign with respect to the wiki page one.
%% Mapping C(t1,t2) into C(t,tau)
%The Fourier transform of the WDF will be performed over tau/2 and not
%simply tau. This implies that in addition to mapping C(t1,t2) into
%C(t,tau), C(t1,t2) must be interpolated also not only over time and tau, but also
%over the tau/2. If this were not done, aliasing would occur due to
%under sampling.  For more information, section "More about" at
% https://www.mathworks.com/help/signal/ref/wvd.html
interpolation_2level_g1SDPR_vs_t1_t2;

%Side note:  MATLAB has implemented its own WDF since R2018b, however it is
%defined for a single time series signal x(t) and not C(t1,t2). Still, it
%can be used for comparison to evaluate the WDF of the coherent component
%of the flux. Below it is shown how it is done for the reflected coherent
%field.

% expect_b_out_vs_time_dag = expect(b_out_vs_time',density_matrix_vs_time);
% sampling_frequency = 1/t_step; % Sampling frequency ps^(-1)
% [WVD_b_out_vs_time_coherent,frequency_WVD,t_list_WVD] = wvd(expect_b_out_vs_time_dag,sampling_frequency); %returns the smoothed pseudo Wigner-Ville distribution
% nb_points_WVD_spectrum = nb_points_time; 
% omega_step_WVD = 2*pi*(sampling_frequency/2)/nb_points_WVD_spectrum; 
% omega_step_WVD_muev = omega_step_WVD/ev*hbar/1e-18; %in mueV
% omega_spectrum_WVD = omega_pulse + frequency_WVD*2*pi; % rad/ps, arrays of angular frequencies. Obs: the zero value from fft corresponds to omega_laser
% omega_spectrum_WVD_muev = omega_spectrum_WVD/ev*hbar/1e-18; % mueV, arrays of angular frequencies, zero-centered
% figure,
% surf(t_list_WVD,omega_spectrum_WVD_muev - omega_pulse_ev*1e6,WVD_b_out_vs_time_coherent/(2*pi/ev*hbar/1e-18),'EdgeColor','none')
% xlabel('time [ps]')
% ylabel('\omega-\omega_p [\mueV]')
% title('Wigner-Ville distribution reflected coherent photons - MATLAB')
% view(2)
% colorbar

%% Fourier Transforming over the delay tau to obtain the WDF
% The technique used here for the Fourier Transform is the Fast Fourier Transform (FFT) algorithm, which 
% can provide accurate spectra in a very fast way provided we use a number of points that is 
% a power of 2 (same number of points in the time and frequency domain). The spectral width of
% the FFT spectrum (here denoted as "FFT_sampling_frequency") and its resolution (here related 
% to the step in angular frequency, "omega_step") are fixed by the width and resolution used in
% in the time domain (related to "t_max_ps" and "t_delay" characterizing the "full_tau_list").
% In practice, we are interested only in a spectral window of width "width_spectral_window_muev",
% and the list of angular frequencies omega in "omega_list" (that we use to plot quantities 
% "versus omega" and calculate integrals) is just a subset of the full FFT spectrum, centered on
% omega_pulse. Also note that the raw FFT spectra have to be adequately transformed into real 
% physical quantities:
%    a) Proper normalization of the FFT result has to be ensured, so that the integral of the
%       spectral density indeed corresponds to the optical flux in ps^-1 (or more precisely 
%       the incoherent component of the optical flux, as discussed above).
%    b) One needs to ensure that the optical spectra are centered on the laser frequency/photon
%       energy, to take into account the  fact that we worked in the rotating frame at this laser 
%       frequency. Using the "fftshift" function just ensures that the first half of the FFT 
%       spectrum corresponds to negative frequencies, and the second half to positive frequencies.
%       In addition, we ensure that the lists of angular frequencies ("omega_list_full_spectrum"
%       for the complete list in rad/ps, "omega_list" for the selected spectral window in rad/ps,
%       and "omega_list_eV" for the selected energy window in eV) are defined to be centered.
%       on the laser frequency/photon energy.
%    c) The list of delays in "full_tau_list" include negative and positive delays, with the 
%       zero delay in the middle of the list. But the FFT algorithm considers that this is the 
%       first point of the list which is the time zero, which leads to a frequency-dependent 
%       phase shift along the full FFT spectrum. This phase shift has to be compensated (the 
%       spectral density of flux is a real quantity), which is done through a multiplication
%       by "phase_shift_compensation_full_FFT_spectrum".
% --> All the angular frequency lists and quantities related to the calculation of the FFT spectra
%     are defined in the "Init_lists_..." subprogram, for clarity.
%WDF-s for reflected field
WDF_interpolated_G1_reflected_vs_time_tau_full_spectrum = 1/(2*pi)*fftshift(fft(interpolated_G1_reflected_vs_time_tau,nb_points_full_spectrum,2),2)./(FFT_sampling_frequency)./phase_shift_compensation_vs_omega_full_spectrum*(ev/hbar*1e-18); % 1/(mueV*ps)
WDF_interpolated_reflected_coherent_vs_time_tau_full_spectrum = 1/(2*pi)*fftshift(fft(interpolated_reflected_coherent_vs_time_tau,nb_points_full_spectrum,2),2)./(FFT_sampling_frequency)./phase_shift_compensation_vs_omega_full_spectrum*(ev/hbar*1e-18); % 1/(mueV*ps)
WDF_interpolated_G1_reflected_vs_time_tau = WDF_interpolated_G1_reflected_vs_time_tau_full_spectrum(:,index_min_zoomed_spectrum:index_max_zoomed_spectrum); % 1/(mueV*ps)
WDF_interpolated_reflected_coherent_vs_time_tau = WDF_interpolated_reflected_coherent_vs_time_tau_full_spectrum(:,index_min_zoomed_spectrum:index_max_zoomed_spectrum); % 1/(mueV*ps)

%WDF-s for transmitted field
WDF_interpolated_G1_transmitted_vs_time_tau_full_spectrum = 1/(2*pi)*fftshift(fft(interpolated_G1_transmitted_vs_time_tau,nb_points_full_spectrum,2),2)./(FFT_sampling_frequency)./phase_shift_compensation_vs_omega_full_spectrum*(ev/hbar*1e-18); % 1/(mueV*ps)
WDF_interpolated_transmitted_coherent_vs_time_tau_full_spectrum = 1/(2*pi)*fftshift(fft(interpolated_transmitted_coherent_vs_time_tau,nb_points_full_spectrum,2),2)./(FFT_sampling_frequency)./phase_shift_compensation_vs_omega_full_spectrum*(ev/hbar*1e-18); % 1/(mueV*ps)
WDF_interpolated_G1_transmitted_vs_time_tau = WDF_interpolated_G1_transmitted_vs_time_tau_full_spectrum(:,index_min_zoomed_spectrum:index_max_zoomed_spectrum);% 1/(mueV*ps)
WDF_interpolated_transmitted_coherent_vs_time_tau = WDF_interpolated_transmitted_coherent_vs_time_tau_full_spectrum(:,index_min_zoomed_spectrum:index_max_zoomed_spectrum);% 1/(mueV*ps)

%WDF-s for emitted field
WDF_interpolated_G1_emitted_vs_time_tau_full_spectrum = 1/(2*pi)*fftshift(fft(interpolated_G1_emitted_vs_time_tau,nb_points_full_spectrum,2),2)./(FFT_sampling_frequency)./phase_shift_compensation_vs_omega_full_spectrum*(ev/hbar*1e-18); % 1/(mueV*ps)
WDF_interpolated_emitted_coherent_vs_time_tau_full_spectrum = 1/(2*pi)*fftshift(fft(interpolated_emitted_coherent_vs_time_tau,nb_points_full_spectrum,2),2)./(FFT_sampling_frequency)./phase_shift_compensation_vs_omega_full_spectrum*(ev/hbar*1e-18); % 1/(mueV*ps)
WDF_interpolated_G1_emitted_vs_time_tau = WDF_interpolated_G1_emitted_vs_time_tau_full_spectrum(:,index_min_zoomed_spectrum:index_max_zoomed_spectrum);% 1/(mueV*ps)
WDF_interpolated_emitted_coherent_vs_time_tau = WDF_interpolated_emitted_coherent_vs_time_tau_full_spectrum(:,index_min_zoomed_spectrum:index_max_zoomed_spectrum);% 1/(mueV*ps)

%WDF-s for diffracted field
WDF_interpolated_G1_diffracted_vs_time_tau_full_spectrum = 1/(2*pi)*fftshift(fft(interpolated_G1_diffracted_vs_time_tau,nb_points_full_spectrum,2),2)./(FFT_sampling_frequency)./phase_shift_compensation_vs_omega_full_spectrum*(ev/hbar*1e-18); % 1/(mueV*ps)
WDF_interpolated_diffracted_coherent_vs_time_tau_full_spectrum = 1/(2*pi)*fftshift(fft(interpolated_diffracted_coherent_vs_time_tau,nb_points_full_spectrum,2),2)./(FFT_sampling_frequency)./phase_shift_compensation_vs_omega_full_spectrum*(ev/hbar*1e-18); % 1/(mueV*ps)
WDF_interpolated_G1_diffracted_vs_time_tau = WDF_interpolated_G1_diffracted_vs_time_tau_full_spectrum(:,index_min_zoomed_spectrum:index_max_zoomed_spectrum);% 1/(mueV*ps)
WDF_interpolated_diffracted_coherent_vs_time_tau = WDF_interpolated_diffracted_coherent_vs_time_tau_full_spectrum(:,index_min_zoomed_spectrum:index_max_zoomed_spectrum);% 1/(mueV*ps)

%% Spectral energy density and spectral energy flux
% The projection property (see Wiki page linked above) of the WDF(t,omega)
% guarantees that its
% integral over the spectrum omega gives the photon flux as a function of
% time, whereas the integral of WDF over time gives the energy spectral
% density of the flux, which is shown in the spectrometer.
% 
% As we saw above, however, the quantity < b'(t,tau) b(t) > is the sum of two contributions:
%   - A contribution  < b'(t,tau) > < b(t) >, corresponding to the coherent part of the flux.
%   - Complementary, a contribution < b'(t,tau)b(t) > - < b'(t,tau) > < b(t) >, induced by the incoherent
%     part of the flux. This contribution tends towards zero for large delays.

ESD_reflected_photons_vs_omega = real(sum(WDF_interpolated_G1_reflected_vs_time_tau,1))*time_step; % 1/mueV
ESD_coherent_reflected_photons_laser_vs_omega = real(sum(WDF_interpolated_reflected_coherent_vs_time_tau,1))*time_step; % 1/mueV
ESD_transmitted_photons_vs_omega = real(sum(WDF_interpolated_G1_transmitted_vs_time_tau,1))*time_step; % 1/mueV
ESD_coherent_transmitted_photons_laser_vs_omega = real(sum(WDF_interpolated_transmitted_coherent_vs_time_tau,1))*time_step; % 1/mueV
ESD_emitted_photons_vs_omega = real(sum(WDF_interpolated_G1_emitted_vs_time_tau,1))*time_step; % 1/mueV
ESD_coherent_emitted_photons_laser_vs_omega = real(sum(WDF_interpolated_emitted_coherent_vs_time_tau,1))*time_step; % 1/mueV
ESD_diffracted_photons_vs_omega = real(sum(WDF_interpolated_G1_diffracted_vs_time_tau,1))*time_step; % 1/mueV
ESD_coherent_diffracted_photons_laser_vs_omega = real(sum(WDF_interpolated_diffracted_coherent_vs_time_tau,1))*time_step; % 1/mueV

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    Plots     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For plot selection: 
% - g1 vs (t1,t2) : 'g'
% - G1 vs (t1,t2) : 'G'
% - interpolated G1 vs (time, tau): 'I'
% - WDF vs (time, omega) : 'W'
% - photon fluxes : 'F'

% - plot_choice = ['G'];
plot_choice = ['g';'G';'I';'W';'F'];
Plot2levelG1wdfPR;