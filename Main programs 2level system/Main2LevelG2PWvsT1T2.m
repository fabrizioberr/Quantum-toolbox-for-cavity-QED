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
% main hard  drive where Windows is installed. This can be for example in:
% 'C:\Program Files\Matlab\R2014a\bin'.

% Warning: for the adiabatic version to converge, the tolerance in
% mesolve.m function must be reduced compared to the defaut values. For
% example:
% ode2file('ode_input.dat',L,rho0,t_list,struct('reltol',7e-8,'abstol',8e-8));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Intensity correlations in the pulsed regime %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section indexed "g2PR_vs_t1_t2" computes the intensity correlations
% in the pulsed regime, as a function of the detection times t1 and t2 in two
% detectors, for the various optical fields. The normalized function g2(t1,t2)
% is computed, as well as the coincidence maps determining the probability of 
% having double-clicks, one at time t1 and the other at time t_2, during the same
% pulse or for uncorrelated pulses. The conditional occupation probabilities,
% modified by the detection of a first photon at time t1, are also computed, as
% well as the normalized g2 as a function of delay tau = t2 - t1, and the averaged  
% g2(0), i.e. the area of the zero-delay peak of the normalized g2 function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%% Choice of full model 'F' or adiabatic model 'A'
model = 'F'; 

%% Experimental conditions (to be edited)
detuning_QD_C_muev = 0; %Detuning between the QD and cavity frequencies, in mueV
detuning_pulse_QD_muev = 0; %Detuning between the pulse central frequency and the QD frequency, in mueV
eta_in = 1; % Injection efficiency for the incoming photons (depends on experimentally-achieved spatial coupling)
Nb_photons_pulse = 1;  % Average number of incoming photons in a pulse. This quantity should be multiplied by eta_in to know the number of incoming photons actually coupled to the optical mode
FWHM_pulse = 15; %in ps, full width at half-maximum of the incoming Gaussian pulse intensity (unit: ps since angular frequencies are in rad/ps)
t_delay = 2*FWHM_pulse; % Time at which the pulse is maximally intense, so that the computation starts when the pulse has not arrived yet
t_max_ps = (t_delay + 4*FWHM_pulse)*5; % Final time where we stop the computation and plots of time evolutions
nb_points_time = 100; % Time resolution/Number of iterations / <100000 otherwise the integrating the master equation gets difficult (odesolve))
t_min = 0.4*FWHM_pulse; % Initial time considered for the computations and plots of time evolutions


%%
% Initialization of parameters, operators, arrays, etc...
Init2levelDeviceParametersOngoingTest;
Init2levelHilbertSpaceAndOperators;
InitMaps2LevelG2PRvsT1T2;

% Definition of the input field in ps^(-1/2), in the form of a fseries (necessary for integrating the master equation)
Standard_deviation_b_in_PR = FWHM_pulse/(2*sqrt(log(2))); %Deduced from the properties of a gaussian function
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
% The system Hamitonian is time-dependent due to the function b_in_fn
% describing the input field b_in(t). 

% In addition, in the case of adiabatic elimination of the cavity mode an effective
% operator $a$ is defined, acting on the QD subspace, based on the formula for
% adiabatic elimination (Eq. 10 of the pdf notes). Since this formula depends
% on b_in(t), we define a time-dependent quantity "a_vs_time", which is an
% array containing, for each time of t_list, the corresponding operator
% "a". 
% In the full model case, to simplify the following calculations, we define the 
% same quantity a_vs_time, yet this time this array contains the same operator 
% (annihilation operator "a" acting on the cavity subspace), replicated for all
% times of t_list.


switch model 
        case 'A' % Adiabatic model
           
            Delta = 2*(omega_pulse-omega_c)/kappa; %normalized laser detuning appearing in Eq.12 of the pdf notes            
            
            %%%%%%%%% Definition of the adiabatic-model Hamiltonian (in the frame rotating at angular frequency omega_pulse)
            H_PR = (omega_eff-omega_pulse)*sigma_dag*sigma - 1i*sqrt(Gamma_0*eta_top)*...
                    ((1-1i*Delta)^(-1)*b_in_fn*sigma_dag-(1+1i*Delta)^(-1)*b_in_fn'*sigma); % Hamiltonian (eq.13)         
            
            %%% Definition of the time-dependant operator "a_vs_time" acting in the QD subspace (Eq. 10 of the pdf notes)
            a_vs_time = -2*(kappa*(1-1i*Delta_QDC))^(-1)*g*sigma*Id_vs_time-2*sqrt(kappa_top)*(kappa*(1-1i*Delta))^(-1)*fsval(b_in_fn,t_list).*Id_vs_time; %annihilation operator a in adiabatic approximation (eq.10), as fseries
            
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
% d_out_vs_time = sqrt(kappa_loss)*a_vs_time; % UNUSED HERE: definition of the operator d_out_vs_time, i.e. the output operator for the diffracted/lost light, in ps^(-1/2)
% e_out = sqrt(gamma_sp)*sigma %  output operator for the light spontaneously emitted outside the cavity mode, in ps^(-1/2), already defined in Init_2level_Hilbert_space_and_operators.m

% Calculation of the total photon flux as a function of time, for the various optical fields 
flux_injected_photons_vs_time = b_in_vs_time.^2; % total flux of injected photons taking into account  eta_in (so only the photons coupled to the cavity mode), in ps(-1)
flux_reflected_photons_vs_time = real(expect(b_out_vs_time'*b_out_vs_time,density_matrix_vs_time)); %flux in ps^(-1)
flux_transmitted_photons_vs_time = real(expect(c_out_vs_time'*c_out_vs_time,density_matrix_vs_time)); %flux in ps^(-1)
% flux_diffracted_photons_vs_time = real(expect(d_out_vs_time'*d_out_vs_time,rho_vs_time)); %flux in ps^(-1)
flux_emitted_photons_vs_time = real(expect(e_out'*e_out,density_matrix_vs_time)); %flux in ps^(-1)

%Computing the number of reflected/transmitted/emitted photons, by integrating over all times in t_list
Nb_reflected_photons = trapz(t_list,flux_reflected_photons_vs_time);
Nb_transmitted_photons = trapz(t_list,flux_transmitted_photons_vs_time);
% Nb_diffracted_photons = trapz(t_list,flux_diffracted_photons_vs_time);
Nb_emitted_photons = trapz(t_list,flux_emitted_photons_vs_time);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Intensity correlations and conditional probabilities %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (method ensuring the best numerical convergence, using normalized density
% matrices in the mesolve functions)

%%% General notes on calculating two-time correlation functions, valid in
%%% PR (pulsed regime) and CW
%
% Experimentally, the second order correlation function g2(t1,t2) is the
% ratio between two probabilities: 
% -  the probability of detecting a photon at time t2 conditioned by a first
%    photon detected at time t1 
% -  the probability of detecting a photon at time t2 without any information 
%    on previous detections
% 
% On the theory side, the standard definition of a normalized correlation
% function g2(t1,t2) is:
%     g2(t1,t2) = <b'(t1)b'(t2)b(t2)b(t1)> / <b'(t2)b(t2)> <b'(t1)b(t1)>
% 
% We thus need to make the link between both views, and derive a way to
% practically compute such quantities with the physical interpretation in
% mind.


% As a preliminary remark, the theoretical notations above  consider the
% Heisenberg representation, where the operators vary in time but the
% density matrix is constant (equal to its initial value at t=0). For
% example:
%     <b'(t1)b(t1)> = Trace [ b'(t1)b(t1) * rho(0)]
% Fortunately, to go to the Schrodinger evolution we can use a very useful
% rule, which tells us that we can put the time evolution in the density
% matrix (instead of the operators) to calculate any average value:
%     <b'(t1)b(t1)> = Trace [ b'(0)b(0) * rho(t1) ]
% In this Schrodinger representation the operators are constant and thus
% equal to their value at time 0, so that :
%     <b'(t1)b(t1)> = Trace [ b'b * rho(t1) ]
% This is why, to calculate  <b'(t1)b(t1)> we  just have to start at t=0
% and make rho(t) evolve up to time t1 using mesolve, then calculate the
% expectation value of b'b:
%     <b'(t1)b(t1)> = expect ( b'b, rho(t1) )

% Now, to calculate g2(t1,t2) we also also have to switch from the
% Heisenberg definition to a practical quantity we can calculate, i.e. the
% expectation value of some operator on some density matrix. By definition
% in the Heisenberg representation:
%     < b'(t1) b'(t2) b(t2) b(t1) > = Trace[ b'(t1) b'(t2) b(t2) b(t1) * rho(0)]
% But with a circular permutation inside the Trace this can also be seen as:
%     < b'(t1) b'(t2) b(t2) b(t1) > = Trace[ b'(t2) b(t2) * b(t1)rho(0)b'(t1) ]

% This gives us some hint that the correlation function, in the Heisenberg
% representation, can be seen as the expectation value of b'b at time t2,
% starting at time t1 from a fictitious density matrix b(t1)rho(0)b'(t1).
% In the Schrodinger representation, however, the operators b and b' are
% constant and only the density matrix evolves between times 0 and t1, and
% between times t1 and t2. But we can also take into account the effect of
% a photon detection event at time t1, which leads to an abrupt change of
% the system density matrix between two times:
%  -   time t=t1^(-): after evolution between 0 and t1, but just before a
%      click (photon detection event)
%  -   time t=t1^(+): just after a click at time t1
% With this in mind, the idea is to compute < b'(t1) b'(t2) b(t2) b(t1) >
% thanks to three steps:
%  -   evolution from 0 to t1^(-), leading to a density matrix rho(t1)
%      "just before a click"
%  -   modification of the system state due to the detection event at t1,
%      leading to a different density matrix "just after a click" 
%  -   further evolution of the system between times t1^(+) and t2

% However, to physically interpret the results and to get a nice numerical
% convergence, we should not work with the fictitious density matrix
% b(t1)rho(0)b'(t1) in the Heisenberg representation, nor its equivalent in
% the Shrodinger representation b rho(t1) b', since it is not normalized
% (its trace is not unity). Instead we define the real/normalized density
% matrix obtained just after a detection event at time t1:
%      rho(t1, just after a click) =  b  rho(t1) b' / expect [ b'b , rho(t1) ]
% This is a valid density matrix since the denominator is :
%      expect [ b'b , rho(t1) ] = Trace [ b'b * rho(t1) ] = Trace [ b rho(t1) b' ],
% and thus "rho(t1, just after a click)" is well normalized with unity
% trace. This division by expect [ b'b , rho(t1) ] also makes sense since
% this quantity is equal to <b'b(t1)>, i.e. one of the two terms in the
% denominator of g2(t1,t2). From this density matrix at time t1^(+), we can
% then deduce the density matrix at time t2, leading to a density matrix
% "rho(t2, conditioned to a click at t1)". With these notations the
% quantity < b'(t1) b'(t2) b(t2) b(t1) > / < b'b(t1) > is equivalent to :
%      < b'(t1) b'(t2) b(t2) b(t1) > / < b'b(t1) > = Trace[ b'b * rho(t2,
%      conditioned to a click at t1) ]
% And with similar notations the quantity < b'b(t2) > is equivalent to:
%       < b'b(t2) > = Trace [ b'b  *  rho(t2) ]
% So we find that the normalized correlation function g2(t1,t2) is indeed
% the ratio between two quantities:
%  -  the photon flux at time t2, conditioned by a previous photon detection
%     event at time t1 
%  -  the photon flux at time t2, unconditioned
% This is exactly the experimentalist's definition of g2(t1,t2). Note that
% in CW we usually take t1 = 0 and t2 = tau, and we take the stationary
% density matrix state "rhoss" both for rho(t1) and rho (t2), since by
% definition rhoss does not evolve with time: only the density matrix
% "rho(tau, conditioned to a click at time 0)", being different from rhoss,
% does evolve with the delay tau
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%  Calculation of conditional density matrices  %%%%%%%%%%%%%
% In the following the "density_matrix_vs_t1_just_after_click_OPERATOR"s are
% defined for each value of time t1 in t_list. Their N-th element represent
% the normalized density matrix just after a click has occurred at the N-th time.

density_matrix_vs_t1_just_after_click_b_out = b_out_vs_time*density_matrix_vs_time*b_out_vs_time'/flux_reflected_photons_vs_time; % Density matrix just after a reflected photon click at time t1 
density_matrix_vs_t1_just_after_click_c_out = c_out_vs_time*density_matrix_vs_time*c_out_vs_time'/flux_transmitted_photons_vs_time; % Density matrix just after a transmitted photon click at time t1 
density_matrix_vs_t1_just_after_click_e_out = e_out*density_matrix_vs_time*e_out'/flux_emitted_photons_vs_time;% Density matrix just after an emitted (outside the mode) photon click at time t1 

% NB: "tic" is used as a "start" time for the measurement of the computing time between "tic" and "toc"
%tic

%Cycle over all times t1 in t_list, corresponding to the moment where a first click occurred
for t1_index = 1:nb_points_time 
    
    % Both t1 and t2 are values in t_list. However, to compute the
    % normalized density matrices vs t2 after a click at t1, we consider
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
    density_matrix_vs_t2_after_click_b_out_at_t1 = [zero_density_matrix_vs_t2_before_t1  mesolve(Liouvillian,density_matrix_vs_t1_just_after_click_b_out{t1_index},t_list(t1_index:end))]; 
    density_matrix_vs_t2_after_click_c_out_at_t1 = [zero_density_matrix_vs_t2_before_t1  mesolve(Liouvillian,density_matrix_vs_t1_just_after_click_c_out{t1_index},t_list(t1_index:end))];
    density_matrix_vs_t2_after_click_e_out_at_t1 = [zero_density_matrix_vs_t2_before_t1  mesolve(Liouvillian,density_matrix_vs_t1_just_after_click_e_out{t1_index},t_list(t1_index:end))]; 
    
    % Evaluation of g2(t1,t2) as described in the general notes above, for t2 >= t1 (and zero otherwise)
    g2_reflected_vs_t1_t2(t1_index,:) = expect(b_out_vs_time'*b_out_vs_time,density_matrix_vs_t2_after_click_b_out_at_t1)./flux_reflected_photons_vs_time; %Auto-correlation g(2)(tau) for the reflected light
    g2_transmitted_vs_t1_t2(t1_index,:) = expect(c_out_vs_time'*c_out_vs_time,density_matrix_vs_t2_after_click_c_out_at_t1)./flux_transmitted_photons_vs_time; %Auto-correlation g(2)(tau) for the transmitted light
    g2_emitted_vs_t1_t2(t1_index,:) = expect(e_out'*e_out,density_matrix_vs_t2_after_click_e_out_at_t1)./flux_emitted_photons_vs_time; %Auto-correlation g(2)(tau) for the light emitted outside the mode

    % Calculation of the conditional occupation probabilities at time t2 after
     % a photon detection event at t1, for t2 >= t1 (and zero otherwise)
    occupation_ground_vs_t1_vs_t2_after_click_b_out_at_t1(t1_index,:) = expect(sigma*sigma_dag,density_matrix_vs_t2_after_click_b_out_at_t1);
    occupation_excited_vs_t1_vs_t2_after_click_b_out_at_t1(t1_index,:) = expect(sigma_dag*sigma,density_matrix_vs_t2_after_click_b_out_at_t1);
    occupation_ground_vs_t1_vs_t2_after_click_c_out_at_t1(t1_index,:) = expect(sigma*sigma_dag,density_matrix_vs_t2_after_click_c_out_at_t1);     
    occupation_excited_vs_t1_vs_t2_after_click_c_out_at_t1(t1_index,:) = expect(sigma_dag*sigma,density_matrix_vs_t2_after_click_c_out_at_t1);
    occupation_ground_vs_t1_vs_t2_after_click_e_out_at_t1(t1_index,:) = expect(sigma*sigma_dag,density_matrix_vs_t2_after_click_e_out_at_t1);
    occupation_excited_vs_t1_vs_t2_after_click_e_out_at_t1(t1_index,:) = expect(sigma_dag*sigma,density_matrix_vs_t2_after_click_e_out_at_t1);
end

%toc

%% 
%%%%% Completion of previous partially-calculated maps vs t1,t2, to include the case where t2 < t1  

% Due to the symmetry between t1 and t2 (we don't know which detector will
% click first), we have to use properties like g2(t1,t2)=g2(t2,t1) to fill
% the voids in the quantities that we have only partially calculated yet
% (since we systematically considered a zero value when t2 < t1. For each
% map function of t1 and t2, this completion is obtained by adding it to
% its transpose (to replace the zeros at t2 < t1) and dividing by 2 the
% elements along the diagonal (to avoid counting twice the case where
% t2=t1). This is done via an ad-hoc matrix idx below:

idx = ones(nb_points_time)-0.5*diag(ones(nb_points_time,1)); % to divide by 2 the elements along the diagonal

% Completed g2(t1,t2) for the various optical fields
g2_emitted_vs_t1_t2 = real((g2_emitted_vs_t1_t2+g2_emitted_vs_t1_t2.')).*idx; 
g2_reflected_vs_t1_t2 = real((g2_reflected_vs_t1_t2+g2_reflected_vs_t1_t2.')).*idx;
g2_transmitted_vs_t1_t2 = real((g2_transmitted_vs_t1_t2+g2_transmitted_vs_t1_t2.')).*idx;

%Completed conditional occupation probabilities for the excited and ground state
occupation_ground_vs_t1_vs_t2_after_click_b_out_at_t1 = (occupation_ground_vs_t1_vs_t2_after_click_b_out_at_t1+occupation_ground_vs_t1_vs_t2_after_click_b_out_at_t1.').*idx;
occupation_excited_vs_t1_vs_t2_after_click_b_out_at_t1 = (occupation_excited_vs_t1_vs_t2_after_click_b_out_at_t1+occupation_excited_vs_t1_vs_t2_after_click_b_out_at_t1.').*idx;
occupation_ground_vs_t1_vs_t2_after_click_c_out_at_t1 = (occupation_ground_vs_t1_vs_t2_after_click_c_out_at_t1+occupation_ground_vs_t1_vs_t2_after_click_c_out_at_t1.').*idx;
occupation_excited_vs_t1_vs_t2_after_click_c_out_at_t1 = (occupation_excited_vs_t1_vs_t2_after_click_c_out_at_t1+occupation_excited_vs_t1_vs_t2_after_click_c_out_at_t1.').*idx;
occupation_ground_vs_t1_vs_t2_after_click_e_out_at_t1 = (occupation_ground_vs_t1_vs_t2_after_click_e_out_at_t1+occupation_ground_vs_t1_vs_t2_after_click_e_out_at_t1.').*idx;
occupation_excited_vs_t1_vs_t2_after_click_e_out_at_t1 = (occupation_excited_vs_t1_vs_t2_after_click_e_out_at_t1+occupation_excited_vs_t1_vs_t2_after_click_e_out_at_t1.').*idx;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Coincidence maps as a function of t1, t2 %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This part computes quantities such as : 
%   -  < b_out'(t1) b_out'(t2) b_out(t2) b_out(t1) >, which is proportional 
%      to the probability of detecting two reflected photons at time t1 and t2
%      in two detectors, during the same pulse (correlated clicks)
%   -  < b_out'(t1) b_out(t1) > < b_out'(t2) b_out(t2) >, which is proportional 
%      to the probability of detecting two reflected photons at time t1 and t2
%      in two detectors, but for different pulses (uncorrelated clicks)

% Uncorrelated photon_coincidences vs (t1,t2), e.g. <b_out'(t1)b_out(t1)> <b_out'(t2)b_out(t2)>
emitted_photon_coincidences_vs_t1_vs_t2_uncorrelated_pulses = real(expect(e_out'*e_out,density_matrix_vs_time))'*real(expect(e_out'*e_out,density_matrix_vs_time));
reflected_photon_coincidences_vs_t1_vs_t2_uncorrelated_pulses = real(expect(b_out_vs_time'*b_out_vs_time,density_matrix_vs_time))'*real(expect(b_out_vs_time'*b_out_vs_time,density_matrix_vs_time));
transmitted_photon_coincidences_vs_t1_vs_t2_uncorrelated_pulses = real(expect(c_out_vs_time'*c_out_vs_time,density_matrix_vs_time))'*real(expect(c_out_vs_time'*c_out_vs_time,density_matrix_vs_time));

% Correlated photon_coincidences vs (t1,t2), e.g. <b_out'(t1) b_out'(t2) b_out(t2) b_out(t1)>
emitted_photon_coincidences_vs_t1_vs_t2 = g2_emitted_vs_t1_t2.*emitted_photon_coincidences_vs_t1_vs_t2_uncorrelated_pulses;
reflected_photon_coincidences_vs_t1_vs_t2 = g2_reflected_vs_t1_t2.*reflected_photon_coincidences_vs_t1_vs_t2_uncorrelated_pulses;
transmitted_photon_coincidences_vs_t1_vs_t2 = g2_transmitted_vs_t1_t2.*transmitted_photon_coincidences_vs_t1_vs_t2_uncorrelated_pulses;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  Normalized g2(tau)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% As in a standard HBT experiment in the pulsed regime, g2(tau) is obtained
% by defining an histogram integrating all the coincidences corresponding
% to a given delay tau (with tau = t2 - t1 being positive or negative). The
% normalization choice here is that the area of the g2(tau) peak is unity 
% for uncorrelated coincidences, i.e. for all peaks except the zero-delay
% peak.

% Normalized g2 (tau) for the zero-delay peak, corresponding to photons emitted
% within the same pulse. It is obtained from the  correlated coincidence map
% <b'(t1)b'(t2)b(t2)b(t1)>, by integrating over all values corresponding to a 
% given delay tau = t2 - t1, with a time bin t_step, and normalizing by Nb_photons^2.
% Note that tau = 0 for t1=t2, corresponding to  the diagonal elements of the
% g2(t1,t2) map, while non-zero delays are obtained outside the diagonal. 
for j = 1:nb_points_time
    normalized_g2_vs_delay_emitted(j) = real(sum(diag(emitted_photon_coincidences_vs_t1_vs_t2,j)))*t_step/Nb_emitted_photons^2;
    normalized_g2_vs_delay_reflected(j) = real(sum(diag(reflected_photon_coincidences_vs_t1_vs_t2,j)))*t_step/Nb_reflected_photons^2;
    normalized_g2_vs_delay_transmitted(j) = real(sum(diag(transmitted_photon_coincidences_vs_t1_vs_t2,j)))*t_step/Nb_transmitted_photons^2;
end

% Normalized g2 (tau) for the other peaks, corresponding to photons emitted
% within different pulses. It is obtained from the uncorrelated coincidence map
% <b'(t2)b(t2)><b'(t1)b(t1)>, by integrating over all values corresponding to a 
% given delay tau = t2 - t1, with a time bin t_step, and normalizing by Nb_photons^2.
for j = 1:nb_points_time
    normalized_g2_vs_delay_uncorrelated_emitted(j) = real(sum(diag(emitted_photon_coincidences_vs_t1_vs_t2_uncorrelated_pulses,j)))*t_step/Nb_emitted_photons^2;
    normalized_g2_vs_delay_uncorrelated_reflected(j) = real(sum(diag(reflected_photon_coincidences_vs_t1_vs_t2_uncorrelated_pulses,j)))*t_step/Nb_reflected_photons^2;
    normalized_g2_vs_delay_uncorrelated_transmitted(j) = real(sum(diag(transmitted_photon_coincidences_vs_t1_vs_t2_uncorrelated_pulses,j)))*t_step/Nb_transmitted_photons^2;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Mean g2 for the zero delay peak %%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In the pulsed regime, what experimentalists usually call the g2(0) is in
% fact the mean value of g2(tau) for the zero-delay peak, or more precisely
% the area of the normalized g2(tau) curve for the zero-delay peak, as compared
% to the unity area obtained for the other uncorrelated peaks.

% Area of normalized g2(tau) for the zero-delay/correlated peak
mean_g2_zero_delay_peak_reflected_photons = trapz(full_tau_list,[flip(normalized_g2_vs_delay_reflected(2:end)) normalized_g2_vs_delay_reflected]);
mean_g2_zero_delay_peak_transmitted_photons = trapz(full_tau_list,[flip(normalized_g2_vs_delay_transmitted(2:end)) normalized_g2_vs_delay_transmitted]);
mean_g2_zero_delay_peak_emitted_photons = trapz(full_tau_list,[flip(normalized_g2_vs_delay_emitted(2:end)) normalized_g2_vs_delay_emitted]);

% Area of normalized g2(tau) for the zero-delay/correlated peak
mean_g2_uncorrelated_peaks_reflected_photons = trapz(full_tau_list,[flip(normalized_g2_vs_delay_uncorrelated_reflected(2:end)) normalized_g2_vs_delay_uncorrelated_reflected]);
mean_g2_uncorrelated_peaks_transmitted_photons = trapz(full_tau_list,[flip(normalized_g2_vs_delay_uncorrelated_transmitted(2:end)) normalized_g2_vs_delay_uncorrelated_transmitted]);
mean_g2_uncorrelated_peaks_emitted_photons = trapz(full_tau_list,[flip(normalized_g2_vs_delay_uncorrelated_emitted(2:end)) normalized_g2_vs_delay_uncorrelated_emitted]);




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    Plots     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For plot selection: 
% - g2 vs (t1,t2) : 'G'
% - photon coincidences (correlated and not) vs (t1,t2) : 'C'
% - ground state occupations vs (t1,t2) : 'O'
% - photon fluxes vs time & g2 vs delay  : 'F'

% plot_choice = ['G'];
plot_choice = ['G';'C';'O';'F'];
Plot2LevelG2PRvsT1T2;
