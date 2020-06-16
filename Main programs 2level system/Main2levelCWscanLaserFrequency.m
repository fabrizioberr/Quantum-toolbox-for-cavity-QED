clear
clc
%close all
%%%%%%%      Common basis to every programs based "two levels" %%%%%%%%%%%
% Important note: these paths must be modified if needed
addpath(genpath('..\QotoolboxV015'));
addpath(genpath('..\CQED subprograms'));
addpath(genpath('..\CQED device parameters'))

% In addition, for the mesolve function to operate the executable files
% (.exe) and batch files (.bat) contained in  '[...] \QotoolboxV015\bin'
% have to be copied to  a folder that is on the Windows system path, in the
% main hard drive where Windows is installed. This can be for example in:
% 'C:\Program Files\Matlab\R2014a\bin'.

% Warning: for the adiabatic version to converge, the tolerance in
% mesolve.m function must be reduced compared to the defaut values. For
% example:
% ode2file('ode_input.dat',L,rho0,t_list,struct('reltol',7e-8,'abstol',8e-8));

%% %%%%%%%%%%% CW - spectra under stationary resonant excitation %%%%%%%%%%
% This script indexed "CW" considers a fixed incoming power and a variable
% angular frequency for the laser. One calculates the spectral response
% associated to the various fields (reflected, transmitted,
% diffracted/lost, and spontaneously emitted outside the cavity mode). One
% finally verifies the conservation of the total photon flux

% Choice of full model 'F' or adiabatic 'A'
model = 'F';

%%% Experimental conditions
detuning_QD_C_muev = 0; %Detuning between the QD and cavity frequencies, in mueV
eta_in = 1; % Injection efficiency for the incoming photons (depends on experimentally-achieved spatial coupling)
P_in_CW_pW=10;% Incoming continuous-wave power in pW

%Parameters for the calculation of spectra
min_detuning_muev=-100; %minimal detuning (left part of the spectrum), in muev
max_detuning_muev=100; %maximal detuning (right part of the spectrum), in muev
nb_points_spectrum=400;

Init2levelDeviceParametersOngoingTest;
Init2levelHilbertSpaceAndOperators;
InitLists2levelCWscanLaserFrequency;

%Incoming power
P_in_CW = P_in_CW_pW*1e-12;% Incoming power in W %%
b_in_CW = sqrt(eta_in*P_in_CW*1e-24/(hbar*omega_c)); % square root of the photon number per unit time, in ps-1/2
total_flux_injected_photons=abs(b_in_CW)^2; % total flux of incoming photons taking into account  eta_in (so only the photons coupled to the cavity mode), in ps-1

% tic % NB: "tic" is used as a "start" time for the measurement of the
% computing time between "tic" and "toc"

%%%%%%%%%%%%%%  Start the calculation of spectra %%%%%%%%%%%%%
for omega_index=1:nb_points_spectrum   %Loop for frequency scan
    
    omega_laser=omega_laser_list(omega_index);  %Current value of omega_laser, in rad/ps
    
    switch model
        case 'A' % Adiabatic case
            
            Delta = 2*(omega_laser-omega_c)/kappa; %normalized laser detuning appearing in Eq.12 of the pdf notes
            
            %%%%%%%%% Definition of the adiabatic-model Hamiltonian (which depends on omega_laser)
            H_CW = (omega_eff-omega_laser)*sigma_dag*sigma...
                - 1i*sqrt(Gamma_0*eta_top)*(b_in_CW*sigma_dag/(1-1i*Delta)-b_in_CW'*sigma/(1+1i*Delta)); % Adiabatic Hamiltonian
                  
            %%% For the redefinition of the operator "a" acting in the QD space
            % (Eq. 10 of the pdf notes
            a = -2*g*sigma/(kappa*(1-1i*Delta_QDC))-2*sqrt(kappa_top)*b_in_CW*Id/(kappa*(1-1i*Delta)); %annihilation operator a in adiabatic approximation
            
            % UNUSED HERE: Ansatz for the annihilation operator, obtained by taking the time derivative of "a" equal to 0
            % (OK for CW but not for PR (pulsed regime) programs)
            % a = -2*g*sigma/(kappa*(1-1i*Delta))-2*sqrt(kappa_top)*b_in_CW*Id/(kappa*(1-1i*Delta));
            
        case 'F'  % Full model
            
            %%%%%%%%% Definition of the full-model Hamiltonian (which depends on omega_laser)
            H_CW = (omega_d-omega_laser)*sigma_dag*sigma...
                 + (omega_c-omega_laser)*a_dag*a...
                 + 1i*g*(sigma_dag*a-a_dag*sigma)...
                 - 1i*sqrt(kappa_top)*b_in_CW*(a_dag-a);
            
    end % end of the "switch model"
    
    %Superoperator associated to the coherent processes (Hamiltonian)
    L_coh = -1i * (spre(H_CW) - spost(H_CW));
    
    %%%%%%%%% Calculation of the Liouvillian superoperator
    Liouvillian  = L_coh + L_incoh; % Total Liouvillian superoperator including both coherent processes (Hamiltonian) and incoherent processes (dissipative jumps)
    
    %%%%%%%%% Calculation of the density matrix corresponding to the stationary state
    rhoss_CW = steady(Liouvillian);
    
    %%%%%%%%% Definition (or re-definition) of the output operators
    if (model=='A' ||  omega_index==1)  %
        % In the full model, the output flux operators are not
        % frequency-dependent and thus need to be defined only the first
        % time In the adiabatic model the value of "a" is
        % frequency-dependent and the output operators have to be
        % redefined for each frequency
        
        b_out = b_in_CW*Id + sqrt(kappa_top)*a; % definition of the operator b_out, i.e. the output operator for the reflected light, in ps^(-1/2)
        c_out = sqrt(kappa_bottom)*a; % definition of the operator c_out, i.e. the output operator for the transmitted light, in ps^(-1/2)
        d_out = sqrt(kappa_loss)*a; % definition of the operator d_out, i.e. the output operator for the diffracted/lost light, in ps^(-1/2)
        % e_out = sqrt(gamma_sp)*sigma; % definition of the operator e_out, i.e. the output operator for the light spontaneously emitted outside the cavity mode, in ps^(-1/2)
        
        % NB: in the adiabatic model the operators could also have been written directly as:
        % b_out = b_in_CW*Id*(1-2*eta_top/(1-1i*Delta))-sqrt(Gamma_0*eta_top)*sigma/(1-1i*Delta_QDC); %output flux operator,  (eq.12)
        % c_out = -2*g*sqrt(kappa_bottom)*sigma/(kappa*(1-1i*Delta_QDC))-2*sqrt(kappa_top*kappa_bottom)*b_in_CW*Id/(kappa*(1-1i*Delta));
        % d_out = -2*g*sqrt(kappa_loss)*sigma/(kappa*(1-1i*Delta_QDC))-2*sqrt(kappa_top*kappa_loss)*b_in_CW*Id/(kappa*(1-1i*Delta));
        % Such formulas are obtained by directly replacing the value of "a"
        % from the adiabatic model, and are thus equivalent to the above,
        % more general definitions In addition, e_out is independent on the
        % experimental conditions and thus defined in the subprogram
        % "Init_2level_Hilbert_space_and operators.m". It is given here for
        % information and clarity purposes only
        
    end
    
    %%%%%%  Calculation of useful expectation values %%%%%%%%%%%%
    
    % Calculation of the total photon flux as a function of omega_laser
    total_flux_reflected_photons_vs_omega(omega_index) = expect(b_out'*b_out,rhoss_CW); % Total reflected flux  = <b_out_dag b_out>, in ps^(-1)
    total_flux_transmitted_photons_vs_omega(omega_index) = expect(c_out'*c_out,rhoss_CW); % Total transmitted flux = <c_out_dag c_out>, in ps^(-1)
    total_flux_diffracted_photons_vs_omega(omega_index) = expect(d_out'*d_out,rhoss_CW);% Total diffracted/lost flux = <d_out_dag d_out>, in ps^(-1)
    total_flux_emitted_photons_vs_omega(omega_index) = expect(e_out'*e_out,rhoss_CW);% Total flux of spontaneously-emitted photons outside the mode = <e_out_dag e_out>, in ps^(-1)
    
    occupation_excited_state_vs_omega(omega_index) = expect(sigma_dag*sigma,rhoss_CW); % Occupation probability for the excited state |e>
    occupation_ground_state_vs_omega(omega_index) = expect(sigma*sigma_dag,rhoss_CW); % Occupation probability for the ground state |g>
end % end of the frequency scan

% toc
% NB: "toc" is used as a "stop" time for the measurement of the computing time between "tic" and "toc"

% Normalization: coefficients of reflectivity, transmission, diffracted/lost part, and spontaneously-emitted part
R_vs_omega = total_flux_reflected_photons_vs_omega / total_flux_injected_photons;
T_vs_omega = total_flux_transmitted_photons_vs_omega / total_flux_injected_photons;
D_vs_omega = total_flux_diffracted_photons_vs_omega / total_flux_injected_photons;
E_vs_omega = total_flux_emitted_photons_vs_omega / total_flux_injected_photons;

%% %%%%%%%%% Plots %%%%%%%%
% For plot selection:
% - Reflected photons : 'R'
% - Transmitted + diffracted/lost photons : 'T'
% - Photons emitted outside the mode : 'E'
% - Occupation probabilities : 'O'

% plot_choice = ['T'];
plot_choice = ['T';'R';'E';'O'];
Plot2levelCWvsLaserFrequency;
