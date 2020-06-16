%M Laser frequency (fixed)

omega_laser_ev = omega_d_ev + detuning_laser_QD_muev*1e-6;
omega_laser = omega_laser_ev*ev/hbar*1e-12; % laser angular frequency in rad/ps, the unit used for calculations


%% Parameters and lists for the evaluation of the temporal evolution as a function of the delay tau 

tau_step = (tau_max)/(nb_points_delay-1); %size of the time step, in ps
tau_list = linspace(0,tau_max,nb_points_delay); % array containg the non-negative time steps
full_tau_list = [fliplr(-tau_list(2:end-1)) , tau_list ];% array containing both positive and negative time steps


%% Parameters for the calculation of spectral densities through the Fast Fourier Transform (FFT) algorithm
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
% time and frequency domain should be a power of 2. 

FFT_sampling_frequency = 1/tau_step; % sampling frequency, also called "sampling rate", in ps^(-1)
nb_points_full_spectrum = 2^nextpow2(length(full_tau_list)); % ensuring a power of 2 for for optimized FFT performance

omega_step = 2*pi*FFT_sampling_frequency/nb_points_full_spectrum; % angular frequency step in the spectrum, in rad/ps
omega_step_muev = omega_step/ev*hbar/1e-18; % photon energy step in mueV

omega_list_full_spectrum = omega_laser + (-nb_points_full_spectrum/2:nb_points_full_spectrum/2-1)*omega_step; % array of angular frequencies in rad/ps
% spectrum centered around omega_laser, since we work in the frame rotating at this angular frequency)


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
% NB: a first-order autocorrelation function is always anti-symmetric in the sense that: 
%          g1(-tau) = g1(tau)*
% This antisymmetry ensures that its Fourier Transform, and thus the corresponding spectral
% density, gives a real physical quantity. However, as discussed in the main program the fft
% function considers that the first signal point corresponds to time 0, while in our case 
% the full_tau_list contains both negative and positive delays, and the antisymmetry point
% (tau=0) is shifted to the middle of the spectrum. This is a very general issue arising
% from the fact that MATLAB arrays have only positive indices, so a signal x[n] defined over 
% n =  -N, -N+1, ..., 0, ... N-1, N is treated by the fft function as it were defined over 
% n = 1 , 2 , ... , 2*N+1. Such a translation of the x[n] signal leads, in its Fourier 
% Transform X[k], to a phase shift which linearly increases with the index k. In our case, 
% the phase shift will depend on the angular frequency omega, and has to be compensated by
% by a phase term denoted "phase_shift_compensation_vs_omega_full_spectrum". 
% For more info on the phase shift compensation, see:
% https://www.mathworks.com/matlabcentral/answers/94874-why-is-the-fft-of-an-anti-symmetric-signal-not-correct

k = 0:(nb_points_full_spectrum-1);
phase_shift_compensation_vs_omega_full_spectrum = exp(-1j*2*pi*(nb_points_full_spectrum/2-1)*k/length(k));


%% Memory preallocation (allows gaining in calculation time)
% NB: The  "_vs_omega" indicates here that this is a list of values related to the different
% values of omega, in rad/ps, in omega_list

spectral_density_flux_reflected_photons_incoh_vs_omega = zeros(1,length(omega_list)); 
spectral_density_flux_transmitted_photons_incoh_vs_omega = zeros(1,length(omega_list));
spectral_density_flux_emitted_photons_incoh_vs_omega = zeros(1,length(omega_list));
