switch model
    case 'F' 
        text_legend_model='full model';
    case 'A'
        text_legend_model='adiabatic model';
end

% Parameters for the text displayed in figure legends
text_legend_QD = ['g=', num2str(g_muev) ' muev    \gamma_{sp}=' num2str(gamma_sp_muev) 'muev     \gamma^*=' num2str(gamma_puredephasing_muev) 'muev'];
text_legend_cav = ['\kappa=' num2str(kappa_muev) 'muev   \eta_{top}=' num2str(eta_top)  '   C=' num2str(g^2/kappa/gamma_decoherence,3)];
text_legend_P_in = ['P_{in}=' num2str(P_in_CW_pW) 'pW  n_0=' num2str(4*eta_top*abs(b_in_CW)^2/kappa,2) '  n_c=' num2str(gamma_decoherence*gamma_sp/(4*g^2),2)];
text_legend_omega = [' \omega_{laser}-\omega_d=' num2str(detuning_laser_QD_muev) ' muev'];
text_legend = {text_legend_QD,text_legend_cav,text_legend_P_in,text_legend_omega};
% NB: C is the cooperativity and n_c the critical photon number, both depend only on the cavity-QED parameters.
% On the contrary, n_0 depends on the incoming power P_in since it is the  calculated number of intracavity photons
% in the absence of QD (i.e. when g=0). When n_0 is much lower than n_c we are in the weak excitation limit.

%Verification of the maximal photon number in the last Fock state (for full model only)
if model == 'F'
   fprintf(['Full model - Verification - maximal occupation of the last Fock state : ' num2str(abs(expect(occupation_last_Fock_state,density_matrix_stationary_state))) ' \n \n'])
end

% Verification of photon flux conservation
fprintf(['g1SDCW - ' text_legend_model ': Relative error on photon flux conservation: ' num2str(abs(flux_injected_photons-flux_reflected_photons-flux_transmitted_photons-flux_diffracted_photons-flux_emitted_photons)/flux_injected_photons) ' \n \n'])

% Verifications on the extreme values of g1(tau) functions at zero and "infinite" delay 
fprintf(['g1SDCW - ' text_legend_model ': Relative error on g1(0) = 1: ' num2str(abs(1-g_1_reflected_vs_tau(1))) ' \n'])
fprintf(['g1SDCW - ' text_legend_model ': Relative error on g1(0) = 1: ' num2str(abs(1-g_1_transmitted_vs_tau(1))) ' \n'])
fprintf(['g1SDCW - ' text_legend_model ': Relative error on g1(0) = 1: ' num2str(abs(1-g_1_emitted_vs_tau(1))) ' \n  \n'])
fprintf(['g1SDCW - ' text_legend_model ': Relative error on g1(infty) = coherent fraction for reflected photons: ' num2str(abs(flux_reflected_photons_laser_coherent/flux_reflected_photons-g_1_reflected_vs_tau(nb_points_delay))) ' \n'])
fprintf(['g1SDCW - ' text_legend_model ': Relative error on g1(infty) = coherent fraction for transmitted photons: ' num2str(abs(flux_transmitted_photons_laser_coherent/flux_transmitted_photons-g_1_transmitted_vs_tau(nb_points_delay))) ' \n'])
fprintf(['g1SDCW - ' text_legend_model ': Relative error on g1(infty) = coherent fraction for emitted photons: ' num2str(abs(flux_emitted_photons_laser_coherent/flux_emitted_photons-g_1_emitted_vs_tau(nb_points_delay))) ' \n \n'])

% Verification of spectral density normalization for the various fields (incoherent part only). Its integral over
% the whole spectrum must correspond to the incoherent photon flux, considering that each spectral density is measured
% in ps-1/muev, which gives a flux in ps-1 when multiplying by the photon energy step in mueV (omega_step_ev*1e-6) 
fprintf(['g1SDCW - ' text_legend_model ': Relative error on the spectral density normalization - reflected photons (incoherent part) : ' num2str(abs(sum(spectral_density_flux_reflected_photons_incoh_vs_omega*(omega_step_muev))-flux_reflected_photons_incoh)/flux_reflected_photons_incoh) ' \n'])
fprintf(['g1SDCW - ' text_legend_model ': Relative error on the spectral density normalization - transmitted photons (incoherent part) :' num2str(abs(sum(spectral_density_flux_transmitted_photons_incoh_vs_omega*(omega_step_muev))-flux_transmitted_photons_incoh)/flux_reflected_photons_incoh) ' \n'])
fprintf(['g1SDCW - ' text_legend_model ': Relative error on the spectral density normalization - emitted photons (incoherent part) : ' num2str(abs(sum(spectral_density_flux_emitted_photons_incoh_vs_omega*(omega_step_muev))-flux_emitted_photons_incoh)/flux_emitted_photons_incoh) ' \n \n'])


%%%%%%%%%%%%%%%%%%  Plots of g(1)(tau)           %%%%%%%%%%%%%%%%%%%%%%
if ismember('G',plot_choice)
    figure('Name',['g1SDCW - ' text_legend_model ' -  |(g1)| vs tau  - reflected photons - Pin = ' num2str(P_in_CW_pW) ' pW - Pulsation laser = ' num2str(omega_laser_ev) ' eV' ],'NumberTitle','off')
    plot(full_tau_list,abs(full_g_1_reflected_vs_tau),'r');
    xlabel('Delay \tau [ps]'); ylabel('|(g1)| - reflected photons');
    xlim([min(full_tau_list) max(full_tau_list)]);ylim([0 real(max(full_g_1_reflected_vs_tau))]);
    text(full_tau_list(ceil(nb_points_delay/20)),0.2*max(full_g_1_reflected_vs_tau),text_legend,'FontSize',9)
    annotation('textbox',[.15 .65 .2 .2],'String',['Coherent   : ' num2str(abs(100*flux_reflected_photons_laser_coherent/flux_reflected_photons)) ' %'],'FitBoxToText','on');
    annotation('textbox',[.15 .55 .2 .2],'String',['Incoherent : ' num2str(abs(100*flux_reflected_photons_incoh/flux_reflected_photons)) ' %'],'FitBoxToText','on');

    figure('Name',['g1SDCW - ' text_legend_model ' -  |(g1)| vs tau - transmitted photons - Pin = ' num2str(P_in_CW_pW) ' pW - Pulsation laser = ' num2str(omega_laser_ev) ' eV' ],'NumberTitle','off')
    plot(full_tau_list,abs(full_g_1_transmitted_vs_tau),'r');
    xlabel('Delay \tau [ps]'); ylabel('|(g1)| - transmitted photons');
    xlim([min(full_tau_list) max(full_tau_list)]);ylim([0 real(max(full_g_1_transmitted_vs_tau))]);
    text(full_tau_list(ceil(nb_points_delay/20)),0.2*max(full_g_1_transmitted_vs_tau),text_legend,'FontSize',9)
    annotation('textbox',[.15 .65 .2 .2],'String',['Coherent   : ' num2str(abs(100*flux_transmitted_photons_laser_coherent/flux_transmitted_photons)) ' %'],'FitBoxToText','on');
    annotation('textbox',[.15 .55 .2 .2],'String',['Incoherent : ' num2str(abs(100*flux_transmitted_photons_incoh/flux_transmitted_photons)) ' %'],'FitBoxToText','on');

    figure('Name',['g1SDCW - ' text_legend_model ' -real |(g1)| vs tau - emitted photons in leaky modes - Pin = ' num2str(P_in_CW_pW) ' pW - Pulsation laser = ' num2str(omega_laser_ev) ' eV' ],'NumberTitle','off')
    plot(full_tau_list,abs(full_g_1_emitted_vs_tau),'r');
    xlabel('Delay \tau [ps]'); ylabel('|(g1)| - emitted photons');
    xlim([min(full_tau_list) max(full_tau_list)]);ylim([0 real(max(full_g_1_emitted_vs_tau))]);
    text(full_tau_list(ceil(nb_points_delay/20)),0.2*max(full_g_1_emitted_vs_tau),text_legend,'FontSize',9)
    annotation('textbox',[.15 .65 .2 .2],'String',['Coherent   : ' num2str(abs(100*flux_emitted_photons_laser_coherent/flux_emitted_photons)) ' %'],'FitBoxToText','on');
    annotation('textbox',[.15 .55 .2 .2],'String',['Incoherent : ' num2str(abs(100*flux_emitted_photons_incoh/flux_emitted_photons)) ' %'],'FitBoxToText','on');
end
%%%%%%%%%%%%%%%%  Plotting of spectra densities  %%%%%%%%%%%%%%%%%%% 
if ismember('S',plot_choice)

    figure('Name',['g1SDCW - ' text_legend_model ' - Spectral density of flux - Incoherent reflected photons - Pin = ' num2str(P_in_CW_pW) ' pW - Pulsation laser = ' num2str(omega_laser_ev) ' eV' ],'NumberTitle','off')
    plot((omega_list_ev-omega_laser_ev)*1e6,abs(spectral_density_flux_reflected_photons_incoh_vs_omega));
    xlabel('\omega-\omega_{laser}  [muev]');ylabel('Flux spectral density [photons / ps / mueV]');
    text(1e6*(omega_list_ev(ceil(nb_points_spectrum/20))-omega_laser_ev),0.9*real(max(abs(spectral_density_flux_reflected_photons_incoh_vs_omega))),text_legend,'FontSize',9)
    title('Flux spectral density - incoherent part of the reflected field')
    
    figure('Name',['g1SDCW - ' text_legend_model ' - Spectral density of flux - Incoherent transmitted photons - Pin = ' num2str(P_in_CW_pW) ' pW - Pulsation laser = ' num2str(omega_laser_ev) ' eV' ],'NumberTitle','off')
    plot((omega_list_ev-omega_laser_ev)*1e6,abs(spectral_density_flux_transmitted_photons_incoh_vs_omega));
    xlabel('\omega-\omega_{laser}  [muev]');ylabel('Flux spectral density [photons / ps / mueV]');
    text(1e6*(omega_list_ev(ceil(nb_points_spectrum/20))-omega_laser_ev),0.9*real(max(abs(spectral_density_flux_transmitted_photons_incoh_vs_omega))),text_legend,'FontSize',9)
    title('Flux spectral density - incoherent part of the transmitted field')

    figure('Name',['g1SDCW - ' text_legend_model ' - Spectral density of flux - Incoherent emitted photons - Pin = ' num2str(P_in_CW_pW) ' pW - Pulsation laser = ' num2str(omega_laser_ev) ' eV' ],'NumberTitle','off')
    plot((omega_list_ev-omega_laser_ev)*1e6, abs(spectral_density_flux_emitted_photons_incoh_vs_omega));
    xlabel('\omega-\omega_{laser}  [muev]');ylabel('Flux spectral density [photons / ps / mueV]');
    text(1e6*(omega_list_ev(ceil(nb_points_spectrum/20))-omega_laser_ev),0.9*real(max(abs(spectral_density_flux_emitted_photons_incoh_vs_omega))),text_legend,'FontSize',9)
    title('Flux spectral density - incoherent part of the emitted field')
end