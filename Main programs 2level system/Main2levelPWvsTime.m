clear
clc
%close all
%% Important note: these paths must be modified if needed
addpath(genpath('..\QotoolboxV015'));
addpath(genpath('..\CQED subprograms'));
addpath(genpath('..\CQED device parameters'))
savepath

% In addition, for the mesolve function to operate the executable files
% (.exe) and batch files (.bat) contained in  '[...] \QotoolboxV015\bin'
% have to be copied to  a folder that is on the Windows system path, in the
% main hard drive where Windows is installed. This can be for example in:
% 'C:\Program Files\Matlab\R2014a\bin'.

% Warning: for the adiabatic version to converge, the tolerance in
% mesolve.m function must be reduced compared to the defaut values. For
% example:
% ode2file('ode_input.dat',L,rho0,t_list,struct('reltol',7e-8,'abstol',8e-8));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Pulsed regime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section indexed "PR" computes the time evolution of the photon flux
% for various fields, and of the exciton occupation, in response to a
% coherent laser pulse with a fixed center frequency, temporal width, and
% average number of photons. One finally verifies that the total photon
% flux has been conserved after the simulation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Choice of full model 'F' or adiabatic model 'A'
model = 'F'; 

%% Experimental conditions
detuning_QD_C_muev = 10; %Detuning between the QD and cavity frequencies, in mueV
detuning_pulse_QD_muev = 0; %Detuning between the pulse central frequency and the QD frequency, in mueV
eta_in = 1; % Injection efficiency for the incoming photons (depends on experimentally-achieved spatial coupling)
Nb_photons = 10;  % Average number of incoming photons in a pulse. This quantity should be multiplied by eta_in to know the number of incoming photons actually coupled to the optical mode
FWHM = 15; %in ps, full width at half-maximum of the incoming Gaussian pulse intensity (unit: ps since angular frequencies are in rad/ps)

%%
%Initialization of parameters, operators, arrays, etc...
Init2levelDeviceParametersOngoingTest;
Init2levelHilbertSpaceAndOperators;
InitLists2levelPRvsTime;

% Definition of the input field in ps^(-1/2), in the form of a fseries (necessary for integrating the master equation)
Standard_deviation_b_in_PR = FWHM/(2*sqrt(log(2))); %Deduced from the properties of a Gaussian function
b_in_fn = fn('gauss',t_delay,Standard_deviation_b_in_PR) * sqrt( eta_in*Nb_photons / ( sqrt(pi) * Standard_deviation_b_in_PR ) ); % square root of the incoming photon number per time unit, in ps-1/2
b_in_vs_time = fsval(b_in_fn,t_list); %scalar array representing b_in vs time

% Initial density matrix before the pulse has started
switch model
    case 'F' %Full model
        psi0 = tensor(Vacuum_state,g_ket); % Initial state: tensorial product of photonic vacuum and QD ground state
        rho0 = psi0*psi0'; % Density matrix corresponding to the initial pure state
    case 'A' % Adiabatic model
        rho0 = g_ket*g_ket'; % Density matrix corresponding to the initial pure state
end

%%% UNUSED HERE: to describe a non-resonant excitation experiment where the
% exciton state |e> is populated at time zero, one can simply change psi0
% and use a very small value for N_in, like 0.0001, to ensure that the
% output fields are almost entirely induced by the initial excitation. In
% such a case we use (e.g. for the full model):
% psi0=tensor(Vacuum_state,e_ket); % Initial state: tensorial product of photonic vacuum and QD excited state
% rho0=psi0*psi0'; % Density matrix corresponding to the initial pure state

%% %%%%%%%% System Hamiltonian and time-dependent operators  %%%%%%%%%%%%
%
% The system Hamiltonian is time-dependent due to the function b_in_fn
% describing the input field b_in(t). 

% In addition, in the case of adiabatic elimination of the cavity mode an
% effective operator $a$ is defined, acting on the QD subspace, based on
% the formula for adiabatic elimination (Eq. 10 of the pdf notes). Since
% this formula depends on b_in(t), we define a time-dependent quantity
% "a_vs_time", which is an array containing, for each time of t_list, the
% corresponding operator "a". In the full model case, to simplify the
% following calculations, we define the same quantity a_vs_time, yet this
% time this array contains the same operator (annihilation operator "a"
% acting on the cavity subspace), replicated for all times of t_list.


switch model 
        case 'A' % Adiabatic model
           
            Delta = 2*(omega_pulse-omega_c)/kappa; %normalized laser detuning appearing in Eq.12 of the pdf notes            
            
            %%%%%%%%% Definition of the adiabatic-model Hamiltonian (which depends on omega_pulse)
            H_PR = (omega_eff-omega_pulse)*sigma_dag*sigma - 1i*sqrt(Gamma_0*eta_top)*...
                    ((1-1i*Delta)^(-1)*b_in_fn*sigma_dag-(1+1i*Delta)^(-1)*b_in_fn'*sigma); % Hamiltonian (eq.13)         
            
            %%% For the redefinition of the operator "a" acting in the QD
            %%% space vs time (Eq. 10 of the pdf notes)
            a_vs_time = -2*(kappa*(1-1i*Delta_QDC))^(-1)*g*sigma*Id_vs_time-2*sqrt(kappa_top)*(kappa*(1-1i*Delta))^(-1)*fsval(b_in_fn,t_list).*Id_vs_time; %annihilation operator a in adiabatic approximation (eq.10), as fseries
            
        case 'F'  % Full model
            
            %%%%%%%%% Definition of the full-model Hamiltonian (which depends on omega_pulse)
            H_PR = (omega_d-omega_pulse)*sigma_dag*sigma...
                 + (omega_c-omega_pulse)*a_dag*a...
                 + 1i*g*(sigma_dag*a-a_dag*sigma)...
                 - 1i*sqrt(kappa_top)*b_in_fn*(a_dag-a); 
            %definition of a_vs_time, even though a is costant in the full
            %model, to reduce the number of "switch" in the following code
             a_vs_time = a*Id_vs_time;
end % end of the "switch model"

%Superoperator associated with the coherent processes (Hamiltonian)
L_coh = -1i * (spre(H_PR) - spost(H_PR)); 

%%%%%%%%% Calculation of the Liouvillian superoperator
Liouvillian  = L_coh + L_incoh; % Total Liouvillian superoperator including both coherent processes (Hamiltonian) and incoherent processes (dissipative jumps)

%%%%%%%%%%% Numerical Integration of the Master Equation %%%%%%%%
% Computation of the density matrix vs time with t_list, i.e. between t_min
% and t_max, requiring a large enough time resolution.
rho_vs_time = mesolve(Liouvillian,rho0,t_list); % second evolution of the system


%%%%%%%%%%%%%%%%%% Definition of the output operators %%%%%%%%%%%%%%%%%%%%%
% These are general formulas for both the adiabatic and full model, depending on "a_vs_time"
b_out_vs_time = b_in_vs_time*Id_vs_time + sqrt(kappa_top)*a_vs_time; % definition of the operator b_out, i.e. the output operator for the reflected light, in ps^(-1/2)
c_out_vs_time = sqrt(kappa_bottom)*a_vs_time; % definition of the operator c_out_vs_time, i.e. the output operator for the transmitted light, in ps^(-1/2)
d_out_vs_time = sqrt(kappa_loss)*a_vs_time; % definition of the operator d_out_vs_time, i.e. the output operator for the diffracted/lost light, in ps^(-1/2)
% e_out = sqrt(gamma_sp)*sigma %  output operator for the light spontaneously emitted outside the cavity mode, in ps^(-1/2), already defined in Init_2level_Hilbert_space_and_operators.m
a_dag_vs_time = a_vs_time';  % definition of describing a_dag(t)

%Calculation of the state population 
expect_sigma_dag_sigma_vs_time = expect(sigma_dag*sigma,rho_vs_time); % List describing <sigma_dag sigma>(t), i.e. the excited state population
expect_sigma_sigma_dag_vs_time = expect(sigma*sigma_dag,rho_vs_time); % List describing <sigma sigma_dag>(t), i.e. the ground state population

% Calculation of the total photon flux as a function of time
flux_injected_photons_vs_time = b_in_vs_time.^2; % total flux of injected photons taking into account  eta_in (so only the photons coupled to the cavity mode), in ps(-1)
flux_reflected_photons_vs_time = real(expect(b_out_vs_time'*b_out_vs_time,rho_vs_time));%flux in ps^(-1)
flux_transmitted_photons_vs_time = real(expect(c_out_vs_time'*c_out_vs_time,rho_vs_time));%flux in ps^(-1)
flux_diffracted_photons_vs_time = real(expect(d_out_vs_time'*d_out_vs_time,rho_vs_time)); %flux in ps^(-1)
flux_emitted_photons_vs_time = real(expect(e_out'*e_out,rho_vs_time));%flux in ps^(-1)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    Plots     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For plot selection: 
% - photon fluxes vs time vs delay  : 'F'
% - occupation probabilities vs time: 'O'

plot_choice = ['F';'O'];
Plot2levelPRvsTime;