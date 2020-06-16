switch model
    case 'F' 
        text_legend_model='full model';
    case 'A'
        text_legend_model='adiabatic model';
end

% Parameters for the text displayed in figure legends
text_legend_QD = ['g=', num2str(g_muev) ' muev    \gamma_{sp}=' num2str(gamma_sp_muev) 'muev     \gamma^*=' num2str(gamma_puredephasing_muev) 'muev'];
text_legend_cav = ['\kappa=' num2str(kappa_muev) 'muev   \eta_{top}=' num2str(eta_top)  '   C=' num2str(g^2/kappa/gamma_decoherence,3)];
text_legend_pulse = ['N_{in}=' num2str(Nb_photons_pulse) '       FWHM_{pulse}=' num2str(FWHM_pulse) 'ps'];
text_legend = {text_legend_QD,text_legend_cav,text_legend_pulse};
% NB: C is the cooperativity 
%% g2 vs (t1,t2)
if ismember('G',plot_choice)
    figure('Name',['g2 emitted vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
    surf(flip(t_list),t_list,flip(real(g2_emitted_vs_t1_t2),2))
    xlabel('t_1 [ps]')
    ylabel('t_2 [ps]')
    title('emitted photon g^{(2)}(t_1,t_2)')
    view(2)
    colorbar
    
    figure('Name',['g2 reflected vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
    surf(flip(t_list),t_list,flip(real(g2_reflected_vs_t1_t2),2))
    xlabel('t_1 [ps]')
    ylabel('t_2 [ps]')
    title('reflected photon g^{(2)}(t_1,t_2)')
    view(2)
    colorbar

    figure('Name',['g2 transmitted vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
    surf(flip(t_list),t_list,flip(real(g2_transmitted_vs_t1_t2),2))
    xlabel('t_1 [ps]')
    ylabel('t_2 [ps]')
    title('transmitted photon g^{(2)}(t_1,t_2)')
    view(2)
    colorbar
end
%% Photon coincidences
if ismember('C',plot_choice)
    % photon coincidences uncorrelated
    figure('Name',['Uncorrelated coincidences emitted photons vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
    surf(flip(t_list),t_list,flip(real(emitted_photon_coincidences_vs_t1_vs_t2_uncorrelated_pulses),2))
    xlabel('t_1 [ps]')
    ylabel('t_2 [ps]')
    title('<e_{out}''(t_1) e_{out}(t_1)> <e_{out}''(t_2) e_{out}(t_2)>')
    view(2)
    colorbar

    figure('Name',['Uncorrelated coincidences reflected photons vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
    surf(flip(t_list),t_list,flip(real(reflected_photon_coincidences_vs_t1_vs_t2_uncorrelated_pulses),2))
    xlabel('t_1 [ps]')
    ylabel('t_2 [ps]')
    title('<b_{out}''(t_1) b_{out}(t_1)> <b_{out}''(t_2) b_{out}(t_2)>')
    view(2)
    colorbar

    figure('Name',['Uncorrelated coincidences transmitted photons vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
    surf(flip(t_list),t_list,flip(real(transmitted_photon_coincidences_vs_t1_vs_t2_uncorrelated_pulses),2))
    xlabel('t_1 [ps]')
    ylabel('t_2 [ps]')
    title('<c_{out}''(t_1) c_{out}(t_1)> <c_{out}''(t_2) c_{out}(t_2)>')
    colorbar

    % photon coincidences correlated
    figure('Name',['Correlated coincidences emitted photons vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
    surf(flip(t_list),t_list,flip(real(emitted_photon_coincidences_vs_t1_vs_t2),2))
    xlabel('t_1 [ps]')
    ylabel('t_2 [ps]')
    title('<e_{out}''(t_1) e_{out}''(t_2) e_{out}(t_2) e_{out}(t_1)>')
    view(2)
    colorbar

    figure('Name',['Correlated coincidences reflected photons vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
    surf(flip(t_list),t_list,flip(real(reflected_photon_coincidences_vs_t1_vs_t2),2))
    xlabel('t_1 [ps]')
    ylabel('t_2 [ps]')
    title('<b_{out}''(t_1) b_{out}''(t_2) b_{out}(t_2) b_{out}(t_1)>')
    view(2)
    colorbar

    figure('Name',['Correlated coincidences transmitted photons vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
    surf(flip(t_list),t_list,flip(real(transmitted_photon_coincidences_vs_t1_vs_t2),2))
    xlabel('t_1 [ps]')
    ylabel('t_2 [ps]')
    title('<c_{out}''(t_1) c_{out}''(t_2) c_{out}(t_2) c_{out}(t_1)>')
    view(2)
    colorbar
end
%% conditioned occupations
if ismember('O',plot_choice)
    figure('Name',['Occupation ground after reflected photon vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
    surf(flip(t_list),t_list,flip(real(occupation_ground_vs_t1_vs_t2_after_click_b_out_at_t1),2))
    xlabel('t_1 [ps]')
    ylabel('t_2 [ps]')
    title('Occupation ground after reflected photon')
    view(2)
    colorbar

%     figure('Name',['Occupation excited after reflected photon vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
%     surf(flip(t_list),t_list,flip(real(occupation_excited_vs_t1_vs_t2_after_click_b_out_at_t1),2))
%     xlabel('t_1 [ps]')
%     ylabel('t_2 [ps]')
%     title('Occupation excited after reflected photon')
%     view(2)
%     colorbar

    figure('Name',['Occupation ground after transmitted photon vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
    surf(flip(t_list),t_list,flip(real(occupation_ground_vs_t1_vs_t2_after_click_c_out_at_t1),2))
    xlabel('t_1 [ps]')
    ylabel('t_2 [ps]')
    title('Occupation ground after transmitted photon')
    view(2)
    colorbar

%     figure('Name',['Occupation excited after transmitted photon vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
%     surf(flip(t_list),t_list,flip(real(occupation_excited_vs_t1_vs_t2_after_click_c_out_at_t1),2))
%     xlabel('t_1 [ps]')
%     ylabel('t_2 [ps]')
%     title('Occupation excited after transmitted photon')
%     view(2)
%     colorbar

    figure('Name',['Occupation ground after emitted photon vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
    surf(flip(t_list),t_list,flip(real(occupation_ground_vs_t1_vs_t2_after_click_e_out_at_t1),2))
    xlabel('t_1 [ps]')
    ylabel('t_2 [ps]')
    title('Occupation ground after emitted photon')
    view(2)
    colorbar

%     figure('Name',['Occupation excited after emitted photon vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
%     surf(flip(t_list),t_list,flip(real(occupation_excited_vs_t1_vs_t2_after_click_e_out_at_t1),2))
%     xlabel('t_1 [ps]')
%     ylabel('t_2 [ps]')
%     title('Occupation excited after emitted photon')
%     view(2)
%     colorbar
end
%% photon fluxes
if ismember('F',plot_choice)
    figure('Name',['Photons flux vs t1 - ' text_legend_model],'NumberTitle','off')
    plot(t_list,flux_injected_photons_vs_time,'k','Displayname','incoming')
    hold on
    plot(t_list,flux_reflected_photons_vs_time,'r','Displayname','reflected')
    plot(t_list,flux_transmitted_photons_vs_time,'b','Displayname','transmitted')
    plot(t_list,flux_emitted_photons_vs_time,'g','Displayname','emitted')
    hold off
    title(['Photon flux - Nb =' num2str(Nb_photons_pulse)])
    xlabel('t_1 [ps]')
    ylabel('[1/ps]')
    legend
    text(t_list(ceil(nb_points_time*0.05)),0.95*max(flux_injected_photons_vs_time),text_legend,'FontSize',9)

    figure('Name',['Normalized g2 vs delay, photon emission - ' text_legend_model],'NumberTitle','off')
    plot(full_tau_list,[flip(normalized_g2_vs_delay_emitted(2:end)) normalized_g2_vs_delay_emitted],'g','Displayname', 'correlated')
    hold on
    plot(full_tau_list,[flip(normalized_g2_vs_delay_uncorrelated_emitted(2:end)) normalized_g2_vs_delay_uncorrelated_emitted],'--g','Displayname', 'uncorrelated')
    title('Histogram < g^{(2)}(\tau) > emitted photons')
    xlabel('\tau [ps]')
    ylabel('counts')
    legend('Location','northeast');
    text(full_tau_list(ceil(nb_points_time*0.05)),max(normalized_g2_vs_delay_uncorrelated_emitted)*0.9,text_legend,'FontSize',9)

    figure('Name',['Normalized g2 vs delay, photon reflection - ' text_legend_model],'NumberTitle','off')
    plot(full_tau_list,[flip(normalized_g2_vs_delay_reflected(2:end)) normalized_g2_vs_delay_reflected],'r','Displayname', 'correlated')
    hold on
    plot(full_tau_list,[flip(normalized_g2_vs_delay_uncorrelated_reflected(2:end)) normalized_g2_vs_delay_uncorrelated_reflected],'--r','Displayname', 'uncorrelated')
    title('Histogram < g^{(2)}(\tau) > reflected photons')
    xlabel('\tau [ps]')
    ylabel('counts')
    legend('Location','northeast');
    text(full_tau_list(ceil(nb_points_time*0.05)),max(normalized_g2_vs_delay_uncorrelated_reflected)*0.95,text_legend,'FontSize',9)

    figure('Name',['Normalized g2 vs delay, photon transmission - ' text_legend_model],'NumberTitle','off')
    plot(full_tau_list,[flip(normalized_g2_vs_delay_transmitted(2:end)) normalized_g2_vs_delay_transmitted],'b','Displayname', 'correlated')
    hold on
    plot(full_tau_list,[flip(normalized_g2_vs_delay_uncorrelated_transmitted(2:end)) normalized_g2_vs_delay_uncorrelated_transmitted],'--b','Displayname', 'uncorrelated')
    title('Histogram < g^{(2)}(\tau) > transmitted photons')
    xlabel('\tau [ps]')
    ylabel('counts')
    legend('Location','northeast');
    text(full_tau_list(ceil(nb_points_time*0.05)),max(normalized_g2_vs_delay_uncorrelated_transmitted)*0.95,text_legend,'FontSize',9)
end
%%

% Displaying the mean g2(0), i.e. the area of the correlated HBT peak
fprintf(['PR - ' text_legend_model ' ' ': Area of mean g2(0) from reflected photons: ' num2str(mean_g2_zero_delay_peak_reflected_photons) '\n'])
fprintf(['PR - ' text_legend_model ' ' ': Area of mean g2(0) from transmitted photons: ' num2str(mean_g2_zero_delay_peak_transmitted_photons) '\n'])
fprintf(['PR - ' text_legend_model ' ' ': Area of mean g2(0) from emitted photons: ' num2str(mean_g2_zero_delay_peak_emitted_photons) '\n'])


% Verification of the area of the uncorrelated HBT peaks, that should be
% normalized to unity
fprintf(['PR - ' text_legend_model ' ' ': Relative error on the area of mean g2 for uncorrelated peaks, for reflected photons: ' num2str((-1 + mean_g2_uncorrelated_peaks_reflected_photons)*100 ) '%% \n'])
fprintf(['PR - ' text_legend_model ' ' ': Relative error on the area of mean g2 for uncorrelated peaks, for transmitted photons: ' num2str((-1 + mean_g2_uncorrelated_peaks_transmitted_photons)*100 ) '%% \n'])
fprintf(['PR - ' text_legend_model ' ' ': Relative error on the area of mean g2 for uncorrelated peaks, for emitted photons: ' num2str((-1 + mean_g2_uncorrelated_peaks_emitted_photons)*100 ) '%% \n'])


% Basic verifications on the intensity correlations 
if min(normalized_g2_vs_delay_uncorrelated_transmitted)<0 || min(normalized_g2_vs_delay_reflected)<0 || min(normalized_g2_vs_delay_emitted)<0 
    fprintf(' \n \n !!!!!!!!! Warning: non-physical negative values in g2 vs delay !!!!!!!!!! \n \n');
end