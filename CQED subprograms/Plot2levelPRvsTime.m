switch model
    case 'F' 
        text_legend_model='Full model';
    case 'A'
        text_legend_model='Adiabatic model';
end

fprintf([text_legend_model ' - Average number of injected photons per pulse : ' num2str(abs(sum(flux_injected_photons_vs_time)*t_step)) ' \n'])
fprintf([text_legend_model ' - Average number of reflected photons per pulse : ' num2str(abs(sum(flux_reflected_photons_vs_time)*t_step)) ' \n'])
fprintf([text_legend_model ' - Average number of transmitted photons per pulse : ' num2str(abs(sum(flux_transmitted_photons_vs_time)*t_step)) ' \n'])
fprintf([text_legend_model ' - Average number of diffracted photons per pulse : ' num2str(abs(sum(flux_diffracted_photons_vs_time)*t_step)) ' \n'])
fprintf([text_legend_model ' - Average number of spontaneously-emitted photons per pulse : ' num2str(abs(sum(flux_emitted_photons_vs_time)*t_step)) ' \n \n'])

%Verification of the conservation of photon number
fprintf([text_legend_model ' - Verification - relative error on the conservation of photon number : ' num2str(abs(sum(flux_reflected_photons_vs_time+flux_transmitted_photons_vs_time+flux_diffracted_photons_vs_time+flux_emitted_photons_vs_time-flux_injected_photons_vs_time)/sum(flux_injected_photons_vs_time))) ' \n \n'])

%Verification of the maximal photon number in the last Fock state (for full
%model only)
if model == 'F'
    last_Fock_state_ket = basis(N,N);
    occupation_last_Fock_state=tensor(last_Fock_state_ket*last_Fock_state_ket',id_QD);
    expect_occupation_last_Fock_state_vs_time=expect(occupation_last_Fock_state,rho_vs_time); %
    fprintf(['Full model - Verification - maximal occupation of the last Fock state : ' num2str(max(abs(expect_occupation_last_Fock_state_vs_time))) ' \n'])
end

%%%%%%%%%%% Plots %%%%%%%%%%
% Parameters for the text displayed in figure legends
text_legend_QD = ['g=', num2str(g_muev) ' muev    \gamma_{sp}=' num2str(gamma_sp_muev) 'muev     \gamma^*=' num2str(gamma_puredephasing_muev) 'muev'];
text_legend_cav = ['\kappa=' num2str(kappa_muev) 'muev   \eta_{top}=' num2str(eta_top)  '   C=' num2str(g^2/kappa/gamma_decoherence,3)];
text_legend_pulse = ['N_{in}=' num2str(Nb_photons) '       FWHM_{pulse}=' num2str(FWHM) 'ps'];
text_legend = {text_legend_QD,text_legend_cav,text_legend_pulse};
% NB: C is the cooperativity 

if ismember('F',plot_choice)
    figure('Name',['PR - ' text_legend_model ' - Photon flux vs time - Nin = ' num2str(Nb_photons) ' - Pulse FWHM = ' num2str(FWHM) ' ps - pulsation laser = ' num2str(omega_pulse_ev) ' eV'],'NumberTitle','off')
    plot(t_list,real(flux_reflected_photons_vs_time),'r',t_list,real(flux_transmitted_photons_vs_time+flux_diffracted_photons_vs_time),'b',t_list,real(flux_emitted_photons_vs_time),'g')
    legend('Flux of reflected photons','Flux of transmitted + diffracted/lost photons', 'Flux of photons spontaneously-emitted outside the mode')
    xlabel('Time t [ps]');ylabel('Photon flux [ps^{-1}]');
    text(t_list(floor(nb_points_time*0.6)),0.7*max(flux_reflected_photons_vs_time),text_legend,'FontSize',9)
end

if ismember('O',plot_choice)
    figure('Name',['PR - ' text_legend_model ' - Occupation probabilities vs time - Nin = ' num2str(Nb_photons) ' - Pulse FWHM = ' num2str(FWHM) ' ps - pulsation laser = ' num2str(omega_pulse_ev) ' eV'],'NumberTitle','off')
    plot(t_list,real(expect_sigma_dag_sigma_vs_time),'r',t_list,real(expect_sigma_sigma_dag_vs_time),'k')
    legend('Occupation of excited state |e>','Occupation of ground state |g>')
    xlabel('Time t [ps]');ylabel('Occupation probability');
    ylim([0 1])
    text(t_list(floor(nb_points_time*0.6)),0.7,text_legend,'FontSize',9)
end