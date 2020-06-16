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
% main hard drive where Windows is installed. This can be for example in:
% 'C:\Program Files\Matlab\R2014a\bin'.
 
% Warning: for the adiabatic version to converge, the tolerance in
% mesolve.m function must be reduced compared to the defaut values. For
% example:
% ode2file('ode_input.dat',L,rho0,t_list,struct('reltol',7e-8,'abstol',8e-8));
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% g1SDCW : g1(tau) and spectral densities in CW %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section called "g1SDCW" evaluates the first order temporal coherence
% g(1) as a function of the delay tau, for the various fields. It indicates
% the fractions of the coherent and incoherent contributions to the photon
% flux for each field. It also computes the Fourier transform of the
% incoherent contribution to g(1)(tau) as a function of omega, which gives
% the spectral density of flux for the incoherent part of the optical
% field. It verifies that the integral of the spectral densities
% corresponds to the photon flux (incoherent part only).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% Choice of full model 'F' or adiabatic model 'A'
model = 'F'; 
% Warning: in the full model 'F', the size of the Fock state must be large
% enough to ensure that the last Fock state is negligibly occupied.
% Artifacts can otherwise arise, especially when increasing the incoming
% power.
 
%%% Experimental conditions
detuning_QD_C_muev = 0; %Detuning between the QD and cavity frequencies, in mueV
detuning_laser_QD_muev = 0; %Detuning between the pulse central frequency and the QD frequency, in mueV
eta_in = 1; % Injection efficiency for the incoming photons (depends on experimentally-achieved spatial coupling)
P_in_CW_pW = 10;% Incoming continuous-wave power in pW
 
% Parameters for the evaluation of the temporal evolution
% 
% NB1: the maximum delay "tau_max" will also dictate the frequency
% resolution of the spectra, given by the angular frequency step
% "omega_step". The corresponding angular frequency lists are defined in
% the "Init_lists_..." subprogram (see also below details on the
% calculation and Fast Fourier Transform (FFT) algorithm)
% 
% NB2: one should be careful that tau_max is large enough to include a good
% approximation of "infinite delays" (check that the g1(tau) function has
% had enough time to truly converge), while keeping a number of points
% large enough to ensure a good temporal resolution. This is especially
% important for high input powers where artifacts can appear.

tau_max = 4000; %maximum positive delay in ps
nb_points_delay = 2^13 + 1; % Number of points in the list of positive delays (tau_list). 
% --> This must be of the form 2^N+1 for FFT optimization. For example: 2^13+1=8193
 
% Parameter defining the observed spectral window, in mueV, to avoid
% plotting and calculating spectra over an inadequately large angular
% frequency ranges. NB: should not exceed the size of the full FFT
% spectrum, which depends on the temporal time step and thus on tau_max and
% nb_points_delay.

width_spectral_window_muev = 300; % width of the spectral window to be displayed, centered on omega_laser
 
%%%%%%%%%%% Plots %%%%%%%%
% For plot selection: 
% - g1(tau) : 'G'
% - Spectral densities of flux : 'S'
%
% plot_choice = ['G'];
plot_choice = ['G';'S'];
 
%%
% Initialization of parameters, operators, arrays, etc...
Init2levelDeviceParametersOngoingTest;
Init2levelHilbertSpaceAndOperators;
InitLists2levelG1SDCWvsDelayAndFrequency;
 
% Incoming power
P_in_CW = P_in_CW_pW*1e-12;% Incoming power of the CW laser, in W %%
b_in_CW = sqrt(eta_in*P_in_CW*1e-24/(hbar*omega_c)); % square root of the photon number per unit time, in ps-1/2
flux_injected_photons = abs(b_in_CW)^2; % total flux of injected photons, taking into account eta_in (so only the photons coupled to the cavity mode), in ps-1
 
switch model 
    case 'A' % Adiabatic elimination model
 
        Delta = 2*(omega_laser-omega_c)/kappa; % Normalized laser-cavity detuning
 
        % Hamiltonian including the QD operators, the cavity effect being included by adiabatic elimination
        Hamiltonian_CW = (omega_eff-omega_laser)*sigma_dag*sigma...
                       - 1i*sqrt(Gamma_0*eta_top)*(b_in_CW*sigma_dag/(1-1i*Delta)-b_in_CW'*sigma/(1+1i*Delta)); % Hamiltonian
 
        % Definition of the operator "a" acting in the QD subspace,
        % following adiabatic elimination (Eq. 10 of the pdf notes)
        a = -2*g*sigma/(kappa*(1-1i*Delta_QDC))-2*sqrt(kappa_top)*b_in_CW*Id/(kappa*(1-1i*Delta)); %annihilation operator a in adiabatic approximation
 
        % UNUSED HERE: Ansatz for the annihilation operator, obtained by taking the time derivative of "a" equal to 0 
        % (OK for CW but not for PR (pulsed regime) programs)
        % a = -2*g*sigma/(kappa*(1-1i*Delta))-2*sqrt(kappa_top)*b_in_CW*Id/(kappa*(1-1i*Delta));
 
    case 'F'  % Full model
 
        % Hamiltonian including both QD and cavity operators
        Hamiltonian_CW = (omega_d-omega_laser)*sigma_dag*sigma...
                       + (omega_c-omega_laser)*a_dag*a...
                       + 1i*g*(sigma_dag*a-a_dag*sigma)...
                       - 1i*sqrt(kappa_top)*b_in_CW*(a_dag-a); 
 
end % end of the "switch model"
 
%Superoperator associated with the coherent processes (Hamiltonian)
L_coh = -1i * (spre(Hamiltonian_CW) - spost(Hamiltonian_CW)); 
 
%%%%%%%%% Calculation of the Liouvillian superoperator
Liouvillian  = L_coh + L_incoh; % Total Liouvillian superoperator including both coherent processes (Hamiltonian) and incoherent processes (dissipative jumps)
 
%%%%%%%%% Calculation of the density matrix corresponding to the stationary state
density_matrix_stationary_state = steady(Liouvillian);
 

%%%%%%%%%%%%%%%%%% Definition of the output operators %%%%%%%%%%%%%%%%%%%%%
% These are general formulas for both the adiabatic and full model
b_out = b_in_CW*Id + sqrt(kappa_top)*a; % definition of the operateur b_out, i.e. the output operator for the reflected light, in ps^(-1/2)
c_out = sqrt(kappa_bottom)*a; % definition of the operateur c_out, i.e. the output operator for the transmitted light, in ps^(-1/2)
d_out = sqrt(kappa_loss)*a; % definition of the operateur d_out, i.e. the output operator for the diffracted/lost light, in ps^(-1/2)
% e_out = sqrt(gamma_sp)*sigma; % definition of the operateur e_out, i.e. the output operator for the light spontaneously emitted outside the cavity mode, in ps^(-1/2)

% NB: in the adiabatic model the operators could also have been written
% directly as a contribution from the empty cavity (term proportionnal to
% the identity operator "Id") and a contribution describing QD emission
% (term proportionnal to the decay operator "sigma"):
% b_out = b_in_CW*Id*(1-2*eta_top/(1-1i*Delta))-sqrt(Gamma_0*eta_top)*sigma/(1-1i*Delta_QDC); %output flux operator,  (eq.12) 
% c_out = -2*g*sqrt(kappa_bottom)*sigma/(kappa*(1-1i*Delta_QDC))-2*sqrt(kappa_top*kappa_bottom)*b_in_CW*Id/(kappa*(1-1i*Delta)); 
% d_out = -2*g*sqrt(kappa_loss)*sigma/(kappa*(1-1i*Delta_QDC))-2*sqrt(kappa_top*kappa_loss)*b_in_CW*Id/(kappa*(1-1i*Delta)); %annihilation operator a in adiabatic approximation
%
% Such formulas are obtained by directly replacing the value of "a" from
% the adiabatic model, and are thus equivalent to the above, more general
% definitions. In addition, e_out is independent on the experimental
% conditions and thus defined in the subprogram
% "Init_2level_Hilbert_space_and operators.m". It is given here for
% information and clarity purposes
  

% Total photon flux for the various fields, i.e. expectation values of the form
% < b'b >, including both coherent and incoherent contributions to the optical flux
flux_reflected_photons=real(expect(b_out'*b_out,density_matrix_stationary_state));  %flux in ps^(-1)
flux_transmitted_photons=real(expect(c_out'*c_out,density_matrix_stationary_state)); %flux in ps^(-1)
flux_diffracted_photons=real(expect(d_out'*d_out,density_matrix_stationary_state)); %flux in ps^(-1)
flux_emitted_photons=real(expect(e_out'*e_out,density_matrix_stationary_state)); %flux in ps^(-1)

% Coherent part of the photon flux, i.e. expectation values of the form <
% b' > < b >, corresponding to the fraction of light that has a
% well-defined amplitude and phase (characterized by the complex field
% amplitude < b > ), with respect to the incoming monochromatic laser.
flux_reflected_photons_laser_coherent = abs(expect(b_out,density_matrix_stationary_state))^2; %flux in ps^(-1)
flux_transmitted_photons_laser_coherent = abs(expect(c_out,density_matrix_stationary_state))^2;%flux in ps^(-1)
flux_diffracted_photons_laser_coherent = abs(expect(d_out,density_matrix_stationary_state))^2;%flux in ps^(-1)
flux_emitted_photons_laser_coherent = abs(expect(e_out,density_matrix_stationary_state))^2;%flux in ps^(-1)

% Incoherent part of the photon flux, corresponding to the fraction of the
% optical flux which has no well-defined phase with respect to the incoming
% laser, and thus cannot interfere with it.
flux_reflected_photons_incoh = flux_reflected_photons-flux_reflected_photons_laser_coherent;%flux in ps^(-1)
flux_transmitted_photons_incoh = flux_transmitted_photons-flux_transmitted_photons_laser_coherent;%flux in ps^(-1)
flux_diffracted_photons_incoh = flux_diffracted_photons-flux_diffracted_photons_laser_coherent;%flux in ps^(-1)
flux_emitted_photons_incoh = flux_emitted_photons-flux_emitted_photons_laser_coherent;%flux in ps^(-1)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Evaluation of g1(tau) %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For a given output field operator b, the two-time first-order
% auto-correlation function, g1(t1,t2), also called the "degree of
% first-order coherence", is defined by:
%      g1(t1,t2)= < b'(t2) b(t1) > / sqrt[ <b'b>(t2) <b'b>(t1) ]
% Such a quantity plays a major role in any experiment using an
% interference between two time-delayed components of an optical field (as
% in a Michelson or Mach-Zender setup ensuring the interference of a short
% and long path). It also plays a huge role in the calculation of optical
% spectra (as in photoluminescence or resonance fluorescence experiments),
% due to the Wiener-Khinchine theorem which directly links the optical
% spectra with the Fourier Transform of the first-order autocorrelation
% function.
% 
% To calculate such a quantity, one needs to use the "Quantum Regression
% Theorem". This theorem allows making the link between the Heisenberg
% representation of a two-time correlation function < A(t2) B(t1) >, where
% operators A and B are time-dependant, and the Schrodinger approach that
% we have to use here  (where the density matrix varies). The "recipe" to
% deduce < A(t2)B(t1) > consists in:
%     - Letting the density matrix evolve from time 0 to t1
%     - Replacing rho(t1) by a fictitious density matrix B*rho(t1) 
%     - Computing the evolution of this fictitious density matrix between
%     times t1 and t2 
%- Calculating the expectation value of operator A using this fictitious density matrix
%
% Note that even if B*rho(t1) is not a real/valid density matrix (it's not
% even Hermitian), we can at least make it of the order of unity, to ensure
% an optimal numerical convergence (especially important in the pulsed
% regime where for example the operator b_out is extremely small at the
% beginning of the pulse). Looking at the definition of g1(t1,t2), we see
% that this is readily obtained by taking the operator B as b/sqrt( <b'b>),
% with the consequence that operator A has to be taken equal to b'/sqrt(
% <b'b> ).
%
% NB1: To read more about the Quantum Regression Theorem, and the validity
% of the "recipe": -->
% http://atomoptics-nas.uoregon.edu/~dsteck/teaching/quantum-optics/quantum-optics-notes.pdf
%     (see in particular Sec. 5.7.3, page 199)
% --> "Quantum Noise" by Gardiner & Zoller
%     (see in particular Eq. 5.2.11 and an alternative formulation in Sec.
%     5.2.3)
% --> "Statistical Methods in Quantum Optics 1" by H. J. Carmichael
%     (see in particular Sec. 1.5 up to equations 1.97 and 1. 98)

% NB2: The approach used in the "g2CW_vs_delay" and "g2PR_vs_t1_t2"
% programs, to calculate second-order autocorrelation functions, is more
% focused on the physical/experimental definition of these quantities. It
% also makes use of real/normalized/valid density matrices, contrary to the
% fictitious density matrices used below. But the theory behind is also
% entirely linked to the use of the Quantum Regression Theorem.

% NB3: In the stationary regime, under CW excitation, one can take any
% initial time as time 0, and consider the density matrix of the stationary
% state at this time. Also, the flux <b'b> does not depend on time, hence
% the degree of coherence only depends on the delay tau through:
%      g1(tau)= < b'(tau) b(0) > / <b'b>(0)
% From this formula one can directly see that at zero delay g1(0) has to be
% equal to unity, indeed:
%      g1(0)= < b'b > / < b'b > = 1
% In addition, for very long delays one can guess that the field at delay
% tau=infty is completely uncorrelated from the field at delay zero, and we
% find that g1(infty) has to be equal to the coherent fraction of the
% optical field, corresponding to the fraction of light that has a
% well-defined amplitude and phase with respect to the incoming
% monochromatic laser. Indeed:
%      g1(infty) = < b' > < b > / < b'b > 
% where we used < b'(infty) b(0) > = < b' (infty) > < b(0) >  (uncorrelated
% fields). Conversely, the difference between g1(0) and g(infty) gives the
% incoherent fraction of the optical field, i.e the fraction of intensity
% with no well-defined phase with respect to the incoming laser.


% Evaluation of g1(tau) for positive delays only, starting from the
% stationary-state density matrix at time 0.

tic %Start timer to evaluate the computation time
g_1_reflected_vs_tau = expect(b_out'/sqrt(flux_reflected_photons),mesolve(Liouvillian,b_out/sqrt(flux_reflected_photons)*density_matrix_stationary_state,tau_list)); %Equivalent to <b_out_normalized_dag(t) b_out_normalized(t+tau)> for t--> infinity
g_1_transmitted_vs_tau = expect(c_out'/sqrt(flux_transmitted_photons),mesolve(Liouvillian,c_out/sqrt(flux_transmitted_photons)*density_matrix_stationary_state,tau_list)); %Equivalent to <c_out_normalized_dag(t) c_out_normalized(t+tau)> for t--> infinity
g_1_emitted_vs_tau = expect(e_out'/sqrt(flux_emitted_photons),mesolve(Liouvillian,e_out/sqrt(flux_emitted_photons)*density_matrix_stationary_state,tau_list)); %Equivalent to <e_out_normalized_dag(t) e_out_normalized(t+tau)> for t--> infinity
toc % Stop timer

% Evaluation of g1(tau) for both negative and positive delays (full_tau_list). 
full_g_1_reflected_vs_tau = [fliplr(conj (g_1_reflected_vs_tau(2:end-1))) g_1_reflected_vs_tau ];
full_g_1_transmitted_vs_tau = [fliplr(conj (g_1_transmitted_vs_tau(2:end-1))) g_1_transmitted_vs_tau ];
full_g_1_emitted_vs_tau = [fliplr(conj (g_1_emitted_vs_tau(2:end-1))) g_1_emitted_vs_tau ];
% NB: It is ensured that the total number of points is a power of 2, as required for
% the calculation of spectra using the Fast Fourier Transform algorithm


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Spectral densities of flux %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The optical spectrum of a CW optical field is characterized by its
% "spectral density of flux", S(omega), whose integral  over the whole
% spectrum gives the total photon flux. Computing such a spectrum requires
% using the "Wiener-Khinchine theorem", which states that the spectral
% density of flux, for a field described by operator "b", is the Fourier
% Transform of the non-normalized autocorrelation function vs delay :
%    S(\omega) = (1/2pi)* \int < b'(\tau) b(0) > exp(- i \omega \tau) d\tau
% The normalization of this quantity implies that indeed the integral of
% S(omega) is the total flux:
%       \int S(omega) d\omega = < b'(0) b(0) > = total flux 
% 
% As we saw above, however, the quantity < b'(tau) b > is the sum of two
% contributions:
%   - A contribution  < b' > < b >, corresponding to the coherent part of
%   the flux.
%     For monochromatic (i.e. infinitely coherent) CW laser excitation,
%     this contribution never decays even at infinite delay (constant value
%     < b' > < b >).
%   - Complementary, a contribution < b'(tau)b > - < b' > < b >, induced
%   by the incoherent
%     part of the flux. This contribution tends towards zero for large
%     delays.
%
% Thus, after Fourier transforming < b'(tau) b > to obtain the spectrum:
%   - The first constant contribution leads to a delta function centered on
%   omega_laser,
%     with an area given by < b' > < b >: this  is the coherent part of the
%     optical field, which is monochromatic and oscillating at the same
%     frequency as the laser (with a phase and amplitude given by the
%     complex number < b > ).
%   - The decaying contribution leads to a continuous spectrum, which can
%   not interfere
%     with the incoming laser since it is not monochromatic: this is the
%     spectrum of the incoherent part of the optical field, whose total
%     area is the incoherent flux <b'b> - <b'><b>. Typically, in the
%     weak-coupling regime this incoherent spectrum has the shape of a
%     single Lorentzian peak at low power, but of a Mollow triplet at high
%     power. This is this incoherent part of the spectrum that is
%     calculated here, by Fourier transforming < b'(tau) b > - < b' >< b >
%     as a function of the delay tau.
%
% NB1: To know more about Wiener-Khinchin theorem and the computation of
% optical spectra:
% --> http://atomoptics-nas.uoregon.edu/~dsteck/teaching/quantum-optics/quantum-optics-notes.pdf
%     (see in particular Chap. 2 for Wiener-Khinchin theorem, and Sec. 5.7
%     for resonance fluorescence)
% --> "The Quantum Theory of Light" by Rodney Loudon
%     (see in particular Chap. 3 for classical optics, and Chap. 8 for
%     resonance fluorescence spectra)
% --> "Statistical Methods in Quantum Optics 1" by Howard Carmichael
%     (see in particular Sec. Z.3.4)
%
% NB2: We choose to define a spectral density in terms of the photon energy
% (in mueV), instead of frequency (in ps^-1) or angular frequency (in
% rad/ps). As a flux is measured here in ps^-1 (number of photons per unit
% time), the dimension of our spectral densities has to be in ps^-1 / mueV.
% Indeed each spectral contribution to the total flux is obtained through
% multiplying the spectral density of flux (in ps^-1 / mueV) by the photon
% energy step (in mueV).
%
% NB3: The default technique used here is the Fast Fourier Transform (FFT)
% algorithm, which can provide accurate spectra in a very fast way provided
% we use a number of points that is a power of 2 (same number of points in
% the time and frequency domain). The spectral width of the FFT spectrum
% (here denoted as "FFT_sampling_frequency") and its resolution (here
% related to the step in angular frequency, "omega_step") are fixed by the
% width and resolution used in in the time domain (related to "tau_max" and
% "tau_step" characterizing the "full_tau_list"). In practice, we are
% interested only in a spectral window of width
% "width_spectral_window_muev", and the list of angular frequencies omega
% in "omega_list" (that we use to plot quantities "versus omega" and
% calculate integrals) is just a subset of the full FFT spectrum, centered
% on omega_laser. Also note that the raw FFT spectra have to be adequately
% transformed into real physical quantities:
%    a) Proper normalization of the FFT result has to be ensured, so that
%    the integral of the
%       spectral density indeed corresponds to the optical flux in ps^-1
%       (or more precisely the incoherent component of the optical flux, as
%       discussed above).
%    b) One needs to ensure that the optical spectra are centered on the
%    laser frequency/photon
%       energy, to take into account the  fact that we worked in the
%       rotating frame at this laser frequency. Using the "fftshift"
%       function just ensures that the first half of the FFT spectrum
%       corresponds to negative frequencies, and the second half to
%       positive frequencies. In addition, we ensure that the lists of
%       angular frequencies ("omega_list_full_spectrum" for the complete
%       list in rad/ps, "omega_list" for the selected spectral window in
%       rad/ps, and "omega_list_eV" for the selected energy window in eV)
%       are defined to be centered. on the laser frequency/photon energy.
%    c) The list of delays in "full_tau_list" includes negative and positive
%    delays, with the
%       zero delay in the middle of the list. But the FFT algorithm
%       considers that this is the first point of the list which is the
%       time zero, which leads to a frequency-dependent phase shift along
%       the full FFT spectrum. This phase shift has to be compensated (the
%       spectral density of flux is a real quantity), which is done through
%       a multiplication by "phase_shift_compensation_full_FFT_spectrum".
% --> All the angular frequency lists and quantities related to the
% calculation of the FFT spectra
%     are defined in the "Init_lists_..." subprogram, for clarity.


% Computing the spectral densities of flux over the full frequency spectrum
spectral_density_flux_reflected_photons_incoh_full_spectrum = 1/(2*pi)*fftshift(fft(full_g_1_reflected_vs_tau*flux_reflected_photons-flux_reflected_photons_laser_coherent,nb_points_full_spectrum)./phase_shift_compensation_vs_omega_full_spectrum)/FFT_sampling_frequency * (ev/hbar*1e-18); % spectral density in ps^-1 / muev
spectral_density_flux_transmitted_photons_incoh_full_spectrum = 1/(2*pi)*fftshift(fft(full_g_1_transmitted_vs_tau*flux_transmitted_photons-flux_transmitted_photons_laser_coherent,nb_points_full_spectrum)./phase_shift_compensation_vs_omega_full_spectrum)/FFT_sampling_frequency * (ev/hbar*1e-18); % spectral density in ps^-1 / muev
spectral_density_flux_emitted_photons_incoh_full_spectrum = 1/(2*pi)*fftshift(fft(full_g_1_emitted_vs_tau*flux_emitted_photons-flux_emitted_photons_laser_coherent,nb_points_full_spectrum)./phase_shift_compensation_vs_omega_full_spectrum)/FFT_sampling_frequency * (ev/hbar*1e-18); % spectral density in ps^-1 / mueV

% Spectral densities of flux as a function of omega, limited to the
% selected spectra window (range of angular frequencies defined by "omega
% list")
spectral_density_flux_reflected_photons_incoh_vs_omega = spectral_density_flux_reflected_photons_incoh_full_spectrum(index_min_zoomed_spectrum:index_max_zoomed_spectrum); % spectrum restricted to the selected spectral window, in ps^-1 / mueV
spectral_density_flux_transmitted_photons_incoh_vs_omega = spectral_density_flux_transmitted_photons_incoh_full_spectrum(index_min_zoomed_spectrum:index_max_zoomed_spectrum); % spectrum restricted to the selected spectral window, ps^-1 / mueV
spectral_density_flux_emitted_photons_incoh_vs_omega = spectral_density_flux_emitted_photons_incoh_full_spectrum(index_min_zoomed_spectrum:index_max_zoomed_spectrum);% spectrum restriced to the selected spectral window, in ps^-1 / mueV


%%
% %%% Alternative method: Fourier transform by explicit calculation of the
% Fourier integrals % This method is very inefficient in terms of computing
% time, yet is crucial to allow verifying % the results obtained with the
% FFT algorithm: it must give the same result as the FFT method.
% 
% tic  % Start timer to evaluate the computation 
% 
% for omega_index=1:nb_points_spectrum
%     spectral_density_flux_reflected_photons_incoh_vs_omega(omega_index) = 1/(2*pi)* sum ( (full_g_1_reflected_vs_tau*flux_reflected_photons-flux_reflected_photons_laser_coherent) .* exp(-1i* (omega_list(omega_index) - omega_laser) .*full_tau_list )) * tau_step * (ev/hbar*1e-18) ; % spectral density in ps^-1 / mueV
%     spectral_density_flux_transmitted_photons_incoh_vs_omega(omega_index) = 1/(2*pi)* sum ( (full_g_1_transmitted_vs_tau*flux_transmitted_photons-flux_transmitted_photons_laser_coherent) .* exp(-1i* (omega_list(omega_index) - omega_laser) .*full_tau_list )) * tau_step * (ev/hbar*1e-18); % spectral density in ps^-1 / mueV
%     spectral_density_flux_emitted_photons_incoh_vs_omega(omega_index) = 1/(2*pi)* sum ( (full_g_1_emitted_vs_tau*flux_emitted_photons-flux_emitted_photons_laser_coherent) .* exp(-1i* (omega_list(omega_index) - omega_laser) .*full_tau_list )) * tau_step * (ev/hbar*1e-18); % spectral density in ps^-1 / mueV
% end
% 
% toc   % Stop timer


%% %%%%%%%%% Plots %%%%%%%%

Plot2LevelG1SDCWvsDelayAndFrequency;