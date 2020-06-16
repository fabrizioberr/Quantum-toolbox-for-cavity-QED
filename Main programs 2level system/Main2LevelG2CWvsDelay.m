clear
clc
close all
%% %%%%%      Common basis to every programs based "two levels" %%%%%%%%%%%
% Important note: these paths must be modified if needed
addpath(genpath('..\QotoolboxV015'));
addpath(genpath('..\CQED_subprograms'));
addpath(genpath('..\CQED device parameters'))
savepath

% In addition, for the mesolve function to operate the executable
% files (.exe) and batch files (.bat) contained in  '[...] \QotoolboxV015\bin'
% have to be copied to  a folder that is on the Windows system path, 
% in the main hard drive where Windows is installed. This can be 
% for example in: 'C:\Program Files\Matlab\R2014a\bin'. 

% Warning: for the adiabatic version to converge, the tolerance in
% mesolve.m function must be reduced compared to the defaut values. For
% example:
% ode2file('ode_input.dat',L,rho0,t_list,struct('reltol',7e-8,'abstol',8e-8));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%  g2CW : photon-photon correlations under CW excitation  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This section indexed "g2CW" computes the intensity correlations as a
% function of delay, i.e g2(tau), for the various optical fields. One
% checks that values of g2(0) and g2(infty) correspond to the
% expected ones, and that the photons emitted outside the mode are single
% photons (g2(0)=0). Depending on the detuning large bunchings with
% g2(0)>>1 can also be observed, on the transmitted field for example (cf
% PRL 101, 203602 2008) or in the strong-coupling regime (see for example
% Nature 575, 622-627(2019)).

% Warning: g2CW may not always lead to a proper numerical convergence,
% especially for very low incoming powers: this can be seen when we obtain
% absurd negative values of g(2) and/or noisy values of g2(infty) (instead
% of a smooth convergence to unity), and/or abrupt discontinuities in the
% g2(tau) function. This can usually be solved by adjusting the incoming
% power and/or the Fock space and/or the mesolve function (mainly the
% tolerances, and potentially the calculation algorithm). But it seems that
% the use of normalized density matrices (normalized = unity trace, even
% for conditional density matrices obtained after a detection event) is
% crucial to get the best numerical convergence.

% Choice of full model 'F' or adiabatic 'A'
model = 'F'; 

%%% Experimental conditions
detuning_QD_C_muev = 0; %Detuning between the QD and cavity frequencies, in mueV
eta_in = 1; % Injection efficiency for the incoming photons (depends on experimentally-achieved spatial coupling)
P_in_CW_pW=1000;% Incoming continuous-wave power in pW

Init2levelDeviceParametersOngoingTest;
Init2levelHilbertSpaceAndOperators;

%Incoming power
P_in_CW = P_in_CW_pW*1e-12;% Incoming power in W %%
b_in_CW = sqrt(eta_in*P_in_CW*1e-24/(hbar*omega_c)); % square root of the photon number per unit time, en ps-1/2

% Angular frequency of the incoming CW laser (fixed)
omega_laser_ev = omega_d_ev; % in eV
omega_laser = omega_laser_ev*ev/hbar*1e-12; % in rad/ps    

%Parameters for the calculation of the time dynamics of g(2)(tau)
tau_max=800; % Maximal delay in ps (the minimal delay is fixed to 0 for calculations)
nb_points_time_g2CW=20000; %%% Time resolution/number of iterations for the curves G(2)(tau)
tau_step_g2CW=(tau_max)/(nb_points_time_g2CW-1); % Duration of a time step
tau_list_g2CW = linspace(0,tau_max,nb_points_time_g2CW); % list of all the positive delays considered in the computation and plots
full_tau_list_g2CW = [fliplr(-tau_list_g2CW(2:end)) , tau_list_g2CW ]; % list of all negative and positive delays for the plots

% tic NB: "tic" is used as a "start" time for the measurement of the
% computing time between "tic" and "toc"

switch model 
        case 'A' % Adiabatic case
           
            Delta = 2*(omega_laser-omega_c)/kappa; %normalized laser detuning appearing in Eq.12 of the pdf notes
            
            %%%%%%%%% Definition of the adiabatic-model Hamiltonian (which depends on omega_laser)
            H_CW = (omega_eff-omega_laser)*sigma_dag*sigma...
                 - 1i*sqrt(Gamma_0*eta_top)*(b_in_CW*sigma_dag/(1-1i*Delta)-b_in_CW'*sigma/(1+1i*Delta)); % Hamiltonian (Eq.13 of the pdf notes)
          
            
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
        
%Superoperator associated with the coherent processes (Hamiltonian)
L_coh = -1i * (spre(H_CW) - spost(H_CW)); 

%%%%%%%%% Calculation of the Liouvillian superoperator
Liouvillian  = L_coh + L_incoh; % Total Liouvillian superoperator including both coherent processes (Hamiltonian) and incoherent processes (dissipative jumps)

%%%%%%%%% Calculation of the density matrix corresponding to the stationary state
rhoss_CW = steady(Liouvillian);

%%%%%%%%% Definition of the output operators. These are general formulas
% for both the adiabatic and full model, depending on the choice of "a"
b_out = b_in_CW*Id + sqrt(kappa_top)*a; % definition of the operator b_out, i.e. the output operator for the reflected light, in ps^(-1/2)
c_out = sqrt(kappa_bottom)*a; % definition of the operator c_out, i.e. the output operator for the transmitted light, in ps^(-1/2)
d_out = sqrt(kappa_loss)*a; % definition of the operator d_out, i.e. the output operator for the diffracted/lost light, in ps^(-1/2)
% e_out = sqrt(gamma_sp)*sigma; % definition of the operator e_out, i.e. the output operator for the light spontaneously emitted outside the
% cavity mode, in ps^(-1/2), already defined in "Init_2level_Hilbert_space_and_operators.m"

% NB: in the adiabatic model the operators could also have been written directly as:
% b_out = b_in_CW*Id*(1-2*eta_top/(1-1i*Delta))-sqrt(Gamma_0*eta_top)*sigma/(1-1i*Delta_QDC); %output flux operator,  (eq.12) 
% c_out = -2*g*sqrt(kappa_bottom)*sigma/(kappa*(1-1i*Delta_QDC))-2*sqrt(kappa_top*kappa_bottom)*b_in_CW*Id/(kappa*(1-1i*Delta)); 
% d_out = -2*g*sqrt(kappa_loss)*sigma/(kappa*(1-1i*Delta_QDC))-2*sqrt(kappa_top*kappa_loss)*b_in_CW*Id/(kappa*(1-1i*Delta)); %annihilation operator a in adiabatic approximation
% Such formulas are obtained by directly replacing the value of "a" from
% the adiabatic model, and are thus equivalent to the above, more general
% definitions In addition, e_out is independent on the experimental
% conditions and thus defined in the subprogram
% "Init_2level_Hilbert_space_and operators.m". It is given here for
% information and clarity purposes only

% Calculation of the total photon flux as a function of omega_laser
total_flux_injected_photons=abs(b_in_CW)^2; % total flux of injected photons taking into account  eta_in (so only the photons coupled to the cavity mode), in ps-1
total_flux_reflected_photons=expect(b_out'*b_out,rhoss_CW);%flux in ps^(-1)
total_flux_transmitted_photons=expect(c_out'*c_out,rhoss_CW);%flux in ps^(-1)
total_flux_diffracted_photons=expect(d_out'*d_out,rhoss_CW);%flux in ps^(-1)
total_flux_emitted_photons=expect(e_out'*e_out,rhoss_CW);%flux in ps^(-1)

%%
%%%%%%%%%% Intensity correlations %%%%%%%%%%%
% (method ensuring the best numerical convergence, using normalized density
% matrices in the mesolve functions)

%%% General notes on calculating two-time correlation functions, valid in
%%% PR and CW %%%

% Experimentally, the second order correlation function g2(t1,t2) is the
% ratio between two probabilities:
%   - the probability of detecting a photon at time t2 conditioned by a
%   first photon detected at time t1
%   - the probability of detecting a photon at time t2 without any
%   information on previous detections
% On the theory side, the standard definition of a normalized correlation
% function g2(t1,t2) is:
%     g2(t1,t2) = <b'(t1)b'(t2)b(t2)b(t1)> / <b'(t2)b(t2)> <b'(t1)b(t1)>
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
%      - time t=t1^(-): after evolution between 0 and t1, but just before a
%      click (photon detection event)
%      - time t=t1^(+): just after a click at time t1
% With this in mind, the idea is to compute < b'(t1) b'(t2) b(t2) b(t1) >
% thanks to three steps:
%      - evolution from 0 to t1^(-), leading to a density matrix rho(t1)
%      "just before a click"
%      - modification of the system state due to the detection event at t1,
%      leading to a different density matrix "just after a click"
%      - further evolution of the system between times t1^(+) and t2

% However, to physically interpret the results and to get a nice numerical
% convergence, we should not work with the fictitious density matrix
% b(t1)rho(0)b'(t1) in the Heisenberg representation, nor its equivalent in
% the Shrodinger representation b rho(t1) b', since it is not normalized
% (its trace is not unity). Instead we define the real/normalized density
% matrix obtained just after a detection event at time t1:
%     rho(t1, just after a click) =  b  rho(t1) b' / expect [ b'b , rho(t1) ]
% This is a valid density matrix since the denominator is :
%      expect [ b'b , rho(t1) ] = Trace [ b'b * rho(t1) ] = Trace [ b rho(t1) b' ],
% and thus "rho(t1, just after a click)" is well normalized with unity
% trace. This division by expect [ b'b , rho(t1) ] also makes sense since
% this quantity is equal to <b'b(t1)>, i.e. one of the two terms in the
% denominator of g2(t1,t2). From this density matrix at time t1^(+), we can
% then deduce the density matrix at time t2, leading to a density matrix
% "rho(t2, conditioned to a click at t1)". With these notations the
% quantity < b'(t1) b'(t2) b(t2) b(t1) > / < b'b(t1) > is equivalent to :
%      < b'(t1) b'(t2) b(t2) b(t1) > / < b'b'(t1) > = Trace[ b'b * rho(t2, conditioned to a click at t1) ]
% And with similar notations the quantity < b'b(t2) > is equivalent to:
%       < b'b'(t2) > = Trace [ b'b  *  rho(t2) ]
% So we find that the normalized correlation function g2(t1,t2) is indeed
% the ratio between two quantities:
%      - the photon flux at time t2, conditioned by a previous photon
%      detection event at time t1
%      - the photon flux at time t2, unconditioned
% This is exactly the experimentalist's definition of g2(t1,t2). Note that
% in CW we usually take t1 = 0 and t2 = tau, and we take the stationary
% density matrix state both for rho(t1) and rho (t2), since by definition
% it does not evolve with time: only the density matrix "rho(tau,
% conditioned to a click at time 0)", being different from the stationary
% state's density matrix, does evolve with the delay tau

% Normalized (unity trace) density matrices just after a click at time 0,
% for the various optical fields:
density_matrix_just_after_reflected_photon_detection=b_out*rhoss_CW*b_out'/total_flux_reflected_photons; % Density matrix just after a reflected photon click at time 0 
density_matrix_just_after_transmitted_photon_detection=c_out*rhoss_CW*c_out'/total_flux_transmitted_photons; % Density matrix just after a transmitted photon click at time 0 
density_matrix_just_after_emitted_photon_detection=e_out*rhoss_CW*e_out'/total_flux_emitted_photons;% Density matrix just after an emitted (outside the mode) photon click at time 0 

% Normalized density matrices as a function of the (positive) delay tau,
% conditioned on the detection of a click at time 0
density_matrix_vs_delay_after_reflected_photon_detection=mesolve(Liouvillian,density_matrix_just_after_reflected_photon_detection,tau_list_g2CW); %Density matrix at time tau, conditioned on a reflected photon click at time 0
density_matrix_vs_delay_after_transmitted_photon_detection=mesolve(Liouvillian,density_matrix_just_after_transmitted_photon_detection,tau_list_g2CW); %Density matrix at time tau, conditioned on a transmitted photon click at time 0
density_matrix_vs_delay_after_emitted_photon_detection=mesolve(Liouvillian,density_matrix_just_after_emitted_photon_detection,tau_list_g2CW); %Density matrix at time tau, conditioned on an emitted (outside the mode) photon click at time 0

% Calculation of the normalized auto-correlation functions g2(tau), for
% positive delays and various optical fields
g2_reflected_vs_delay=expect(b_out'*b_out,density_matrix_vs_delay_after_reflected_photon_detection)/total_flux_reflected_photons; %Auto-correlation g(2)(tau) for the reflected light
g2_transmitted_vs_delay=expect(c_out'*c_out,density_matrix_vs_delay_after_transmitted_photon_detection)/total_flux_transmitted_photons; %Auto-correlation g(2)(tau) for the transmitted light
g2_emitted_vs_delay=expect(e_out'*e_out,density_matrix_vs_delay_after_emitted_photon_detection)/total_flux_emitted_photons; %Auto-correlation g(2)(tau) for the light emitted outside the mode

% Calculation of g2(tau) for both negative and positive delays
full_g2_reflected_vs_delay= [fliplr(g2_reflected_vs_delay(2:end) ) g2_reflected_vs_delay ];
full_g2_transmitted_vs_delay= [fliplr(g2_transmitted_vs_delay(2:end) ) g2_transmitted_vs_delay ];
full_g2_emitted_vs_delay= [fliplr(g2_emitted_vs_delay(2:end) ) g2_emitted_vs_delay ];

% Calculation of the conditional occupation probabilities, as a function
% of the delay tau after a photon detection event
occupation_ground_vs_delay_after_reflected_photon_detection=expect(sigma*sigma_dag,density_matrix_vs_delay_after_reflected_photon_detection);
occupation_excited_vs_delay_after_reflected_photon_detection=expect(sigma_dag*sigma,density_matrix_vs_delay_after_reflected_photon_detection);
occupation_ground_vs_delay_after_transmitted_photon_detection=expect(sigma*sigma_dag,density_matrix_vs_delay_after_transmitted_photon_detection);
occupation_excited_vs_delay_after_transmitted_photon_detection=expect(sigma_dag*sigma,density_matrix_vs_delay_after_transmitted_photon_detection);
occupation_ground_vs_delay_after_emitted_photon_detection=expect(sigma*sigma_dag,density_matrix_vs_delay_after_emitted_photon_detection);
occupation_excited_vs_delay_after_emitted_photon_detection=expect(sigma_dag*sigma,density_matrix_vs_delay_after_emitted_photon_detection);

%toc

%%% UNUSED HERE: alternative method with non-normalized conditional density
%%% matrices
% NB: we denote G2(tau) (with a capital "G") the non-normalized intensity
% correlations of the form <b'(0) b'(tau) b(tau) b(0)>, where the output
% fields are in ps-{1/2}. Complementary, we denote g2(tau) the normalized
% intensity correlations <b'(0) b'(tau) b(tau) b(0)>/<b' b><b' b>
 
% rho_vs_delay_if_rho0_equal_a_rhoss_a_dag=mesolve(Liouvillian,a*rhoss_CW*a_dag,tau_list_g2CW);  % Sufficient to calculate correlations for c_out=sqrt(kappa_bottom)a et d_out=sqrt(kappa_loss)a
% rho_vs_delay_if_rho0_equal_sigma_rhoss_sigma_dag=mesolve(Liouvillian,sigma*rhoss_CW*sigma_dag,tau_list_g2CW);  % Sufficient to calculate correlations of e_out=sqrt(gamma_sp) sigma
% rho_vs_delay_if_rho0_equal_rhoss_a_dag=mesolve(Liouvillian,rhoss_CW*a_dag,tau_list_g2CW); % Necessary to calculate correlations of b_out=b_in+sqrt(kappa_top) a
% rho_vs_delay_if_rho0_equal_a_rhoss=mesolve(Liouvillian,a*rhoss_CW,tau_list_g2CW); % Necessary to calculate correlations of b_out=b_in+sqrt(kappa_top) a
% 
% G2_reflected_vs_tau=  kappa_top*expect(b_out'*b_out,rho_vs_delay_if_rho0_equal_a_rhoss_a_dag)...
%                      + abs(b_in_CW)^2*expect(b_out'*b_out,rhoss_CW)...
%                      + sqrt(kappa_top)*b_in_CW*expect(b_out'*b_out,rho_vs_delay_if_rho0_equal_a_rhoss)...
%                      + sqrt(kappa_top)*conj(b_in_CW)*expect(b_out'*b_out,rho_vs_delay_if_rho0_equal_rhoss_a_dag);
% G2_transmitted_vs_tau=kappa_bottom*expect(c_out'*c_out,rho_vs_delay_if_rho0_equal_a_rhoss_a_dag);
% G2_diffracted_vs_tau=kappa_loss*expect(d_out'*d_out,rho_vs_delay_if_rho0_equal_a_rhoss_a_dag);
% G2_emitted_vs_tau=gamma_sp*expect(e_out'*e_out,rho_vs_delay_if_rho0_equal_sigma_rhoss_sigma_dag);
% 
% 
% %%%    Calculation of the G(2) for all delays, through the symmetrization G2(tau) = conj(G2(-tau))
% full_G2_reflected_vs_tau = [fliplr(G2_reflected_vs_tau(2:end) ) G2_reflected_vs_tau ];
% full_G2_transmitted_vs_tau = [fliplr(G_2_transmitted_vs_tau(2:end) ) G_2_transmitted_vs_tau ];
% full_G2_diffracted_vs_tau = [fliplr(G_2_diffracted_vs_tau(2:end) ) G_2_diffracted_vs_tau ];
% full_G2_emitted_vs_tau = [fliplr(G_2_emitted_vs_tau(2:end) ) G_2_emitted_vs_tau ];
% 
% %%%    Calculation of the normalized g(2) for all delays
% full_g2_reflected_vs_delay = full_G2_reflected_vs_tau/(total_flux_reflected_photons)^2;
% full_g2_transmitted_vs_tau = full_G_2_transmitted_vs_tau/(total_flux_transmitted_photons)^2;
% full_g2_diffracted_vs_tau = full_G_2_diffracted_vs_tau/(total_flux_diffracted_photons)^2;
% full_g2_emitted_vs_tau = full_G_2_emitted_vs_tau/(total_flux_emitted_photons^2);


%% %%%%%%%%% Plots %%%%%%%%
% For plot selection: 
% - g2 from reflected photons : 'R'
% - g2 from transmitted + diffracted/lost photons : 'T'
% - g2 from photons emitted outside the mode : 'E'
% - associated conditioned occupation probabilities : 'O'

%plot_choice = ['T'];
plot_choice = ['T';'R';'E';'O'];
Plot2LevelG2CWvsDelay;