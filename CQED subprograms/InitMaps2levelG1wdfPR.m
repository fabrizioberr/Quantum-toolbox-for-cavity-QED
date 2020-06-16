% Parameters of the incoming gaussian pulse
omega_pulse_ev=omega_d_ev+detuning_pulse_QD_muev*1e-6; %center energy of the incoming pulse, in eV
omega_pulse=omega_pulse_ev*ev/hbar*1e-12; %in rad/ps

% Parameters for the computation of time evolutions
t_step=(t_max_ps-t_min)/(nb_points_time-1); % Duration of a time step
t_list=linspace(t_min,t_max_ps,nb_points_time); % list of all the times considered in the computation and plots

% Steps in time and delay for the Wigner Distribution Function (WDF).
% Notice that even though tau = t2-t1, its step here is half the step in
% the density matrix evolution to correctly Fourier transform over
% frequency later, according to the WDF definition. More info in the "main"
% script.

tau_step = t_step/2; %important to avoid aliasing
time_step = t_step/2;

time_list = (t_list(1):time_step:t_list(end));
tau_list = (t_list(1)-t_list(end):tau_step:t_list(end)-t_list(1));
%%
% Parameters for the computation of preliminary time evolution,
% between 0 (long before the pulse) and t_min (time at which we want
% to start plotting and integrating the physical quantities
nb_points_time_before_t_min = 5; % (Low) time resolution for first evolution of the system for initialization
t_list_before_t_min = linspace(0,t_min,nb_points_time_before_t_min); % time array for first evolution of the sistem

% Initialization of qo array of identity operator
Id_vs_time = qo;
for t1_index = 1:nb_points_time
    Id_vs_time{t1_index} = Id;
end

%%% Preallocation of the memory to save computing time

% Initialization of g1 vs (t1,t2)
g1_reflected_vs_t1_t2 = zeros(nb_points_time,nb_points_time);
g1_transmitted_vs_t1_t2 = zeros(nb_points_time,nb_points_time);
g1_diffracted_vs_t1_t2 = zeros(nb_points_time,nb_points_time);
g1_emitted_vs_t1_t2 = zeros(nb_points_time,nb_points_time);

% Initialization of the "zero" density matrix which will be used to fill the 
% conditional density matrices, using 0 values for t2 < t1
zero_density_matrix_vs_t2_before_t1 = qo;
%% Parameters for the calculation of Fourier Transform over tau through the Fast Fourier Transform (FFT) algorithm
% 
% NB1: See comments in the main program for a definition of spectral densities and main concepts involved
%
% NB2: For more information in the matlab FFT function, see: 
%  https://fr.mathworks.com/help/matlab/ref/fft.html
% A discussion on the proper normalization of the FFT signal (which should be performed by
% dividing the fft function's result by the signal's sampling frequency) can be found here:
% https://math.stackexchange.com/questions/636847/understanding-fourier-transform-example-in-matlab
% Such a normalization choice allows respecting the Parseval's theorem, i.e.: 
%       sum_{n=0}^{N-1} |x[n]|^2 = 1/N*sum_{k=0}^{N-1} |X[k]|^2,
% with x[n] the signal and X[k] its discrete Fourier Transform, n and k being positive
% indices between 0 and N-1. Indeed, Parseval's theorem is at the heart of our normalization
% choice that the integral of the spectra density of flux should be the
% total flux. 
% 
% NB3: To get a really fast algorithm, the number of points used in the
% time and frequency domain should be a power of 2. However, in this
% specific program the number of elements is 4*nb_points_time-3, so zero
% padding is performed by the fft-s defined in the main.

FFT_sampling_frequency = 1/t_step; % sampling frequency, also called "sampling rate", in ps^(-1)
nb_points_full_spectrum = 2^nextpow2(length(tau_list)); % ensuring a power of 2 for for optimized FFT performance
omega_step = 2*pi*FFT_sampling_frequency/nb_points_full_spectrum; % angular frequency step in the spectrum, in rad/ps
omega_step_muev = omega_step/ev*hbar/1e-18; % photon energy step in mueV

omega_list_full_spectrum = omega_pulse + (-nb_points_full_spectrum/2:nb_points_full_spectrum/2-1)*omega_step; % array of angular frequencies in rad/ps
omega_list_full_spectrum_muev = omega_list_full_spectrum/ev*hbar/1e-18;
% spectrum centered around omega_pulse, since we work in the frame rotating at this angular frequency)

%% Parameters for the zoomed spectrum, i.e. the list of angular frequencies of interest in the slected spectral window.
%
% NB: This zoomed spectrum is simply a subset of the previous one, between a minimal index
% and a maximal one that are defined below.

index_min_zoomed_spectrum = nb_points_full_spectrum/2+1-round(width_spectral_window_muev/omega_step_muev/2); 
index_max_zoomed_spectrum = nb_points_full_spectrum/2+1+round(width_spectral_window_muev/omega_step_muev/2);

omega_list = linspace(omega_list_full_spectrum(index_min_zoomed_spectrum),omega_list_full_spectrum(index_max_zoomed_spectrum),index_max_zoomed_spectrum-index_min_zoomed_spectrum+1); % list of angular frequencies, in rad/ps over the full FFT spectrum

nb_points_spectrum = length(omega_list);

%% Parameters for the zoomed spectrum expressed in photon energy 

omega_list_ev = omega_list/ev*hbar/1e-12; % selected list of photon energies, in eV

%% Parameters for the compensation of the FFT phase shift compensation
% 
% NB: WVD gives a real quantity. However, as discussed in the main program the fft
% function considers that for each row, the first signal point corresponds to delay 0, while in our case 
% each row contains both negative and positive delays. This is a very general issue arising
% from the fact that MATLAB arrays have only positive indices, so a signal x[n] defined over 
% n =  -N, -N+1, ..., 0, ... N-1, N is treated by the fft function as it were defined over 
% n = 1 , 2 , ... , 2*N+1. Such a translation of the x[n] signal leads, in its Fourier 
% Transform X[k], to a phase shift which linearly increases with the index k. In our case, 
% the phase shift will depend on the angular frequency omega, and has to be compensated by
% by a phase term denoted "phase_shift_compensation_vs_omega_full_spectrum". 
% For more info on the phase shift compensation, see:
% https://www.mathworks.com/matlabcentral/answers/94874-why-is-the-fft-of-an-anti-symmetric-signal-not-correct

k = 0:(nb_points_full_spectrum-1);
%NB 1: the command "repmat" is used since the fft returns a matrix to be
%phaseshifted only along the rows. Therefore, "repmat" is used to generate
%a matrix with identical rows.
%NB 2: "length(tau_list)-1)/2" is the number of negative delay elements for
%each row
phase_shift_compensation_vs_omega_full_spectrum = repmat(exp(-1j*2*pi*((length(tau_list)-1)/2)*k/length(k)),2*nb_points_time-1,1);
