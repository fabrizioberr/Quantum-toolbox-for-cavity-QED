if ~ismember(model,['F'; 'A'])
    f = errordlg('Not valid input for full / adiabatic model selection.','Error');
end

%%%% Physical constants (do not change) %%%%
ev=1.60217646e-19;h=6.626068e-34; hbar=h/(2*pi);

% Parameters: energies of the mode ("c" for "cavity" and of the QD ("d" for
% "dot") in eV. NB: all variable names have to end in "_ev" when they are
% in electron-volts. Afterwards the energies in electron-volts will all be
% converted in angular frequencies in rad/ps for subsequent calculations.
omega_c_ev = 1.329810; % Cavity mode energy
omega_d_ev = omega_c_ev + detuning_QD_C_muev*1e-6; % QD transition energy

% Parameters of the QD (in ueV as indicated by the name ending by "_muev"
g_muev = 17; % Light matter coupling in mueV
gamma_sp_muev = 0.6; %Spontaneous emission of leaky modes (i.e. not in the
                        %cavity mode)
gamma_puredephasing_muev = 0; %Pure dephasing in mueV
gamma_decoherence_muev = gamma_sp_muev/2 + gamma_puredephasing_muev;

%Parameters of the cavity
kappa_muev = 400; % Cavity intensity decay rate, in ueV
eta_top = 0.7; %Extraction efficiency for the top Bragg mirror
eta_bottom = 0.1; %Extraction efficiency for the bottom Bragg mirror
eta_loss = 1-eta_top-eta_bottom; %Losses induced by lateral diffraction or absorption

% Conversions in rad/ps
omega_c = omega_c_ev*ev/hbar*1e-12; %in rad/ps
omega_d = omega_d_ev*ev/hbar*1e-12; %in rad/ps
g = g_muev*10^-6*ev/hbar*1e-12;  %in rad/ps
gamma_decoherence = gamma_decoherence_muev*1e-6*ev/hbar*1e-12;  %in rad/ps
gamma_sp = gamma_sp_muev*1e-6*ev/hbar*1e-12;  %in rad/ps
gamma_puredephasing = gamma_puredephasing_muev*1e-6*ev/hbar*1e-12;  %in rad/ps
kappa = kappa_muev*1e-6*ev/hbar*1e-12;  %in rad/ps
kappa_top = kappa*eta_top;  %in rad/ps
kappa_bottom = kappa*eta_bottom;  %in rad/ps
kappa_loss = kappa*eta_loss;  %in rad/ps

if strcmp(model,'F')
    N=10; %Maximal number of photons in the cavity mode (troncature of the Fock space if complete model is used)
end

% Parameters useful for physical interpretation of data, and required for
% the "adiabatic elimination" model
Delta_QDC = 2*(omega_d-omega_c)/kappa; %normalized QD-cavity detuning
Gamma_0 = 4*g^2/kappa; %Purcell-enhanced emission rate at zero detuning 
Gamma_m = Gamma_0/(1 + Delta_QDC^2); % Purcell-enhanced emission rate 
Gamma_tot = Gamma_m + gamma_sp; %total emission rate
omega_eff = omega_d + 0.5*Gamma_0*Delta_QDC/(1 + Delta_QDC^2); %cavity induced frequency shift
