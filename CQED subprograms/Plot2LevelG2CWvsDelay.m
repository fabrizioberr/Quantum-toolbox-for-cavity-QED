switch model
    case 'F' 
        text_legend_model='full model';
    case 'A'
        text_legend_model='adiabatic model';
end

% Parameters for the text displayed in figure legends
text_legend_QD=['g=', num2str(g_muev) ' muev    \gamma_{sp}=' num2str(gamma_sp_muev) 'muev     \gamma^*=' num2str(gamma_puredephasing_muev) 'muev'];
text_legend_cav=['\kappa=' num2str(kappa_muev) 'muev   \eta_{top}=' num2str(eta_top)  '   C=' num2str(g^2/kappa/gamma_decoherence,3)];
text_legend_P_in=['P_{in}=' num2str(P_in_CW_pW) 'pW   n_c=' num2str(gamma_decoherence*gamma_sp/(4*g^2),2)];
text_legend_omega=['\omega_{laser}-\omega_d=' num2str(1e6*(omega_laser_ev-omega_d_ev),3) 'muev    \omega_d-\omega_c=' num2str(1e6*(omega_d_ev-omega_c_ev),3)  'muev'];
text_legend={text_legend_QD,text_legend_cav,text_legend_P_in,text_legend_omega};
% NB: C is the cooperativity and n_c the critical photon number, both depend only on the cavity-QED parameters.
% On the contrary, n_0 depends on the incoming power P_in since it is the  calculated number of intracavity photons
% in the absence of QD (i.e. when g=0). When n_0 is much lower than n_c we are in the weak excitation limit.


%%%%%%%%%%%%%%%%% Verifications %%%%%%%%%%%%%%%%

% Verification on the conservation of photon flux
fprintf(['Verification - relative error on the conservation of photon flux : ' num2str(abs(total_flux_injected_photons-total_flux_reflected_photons-total_flux_transmitted_photons-total_flux_diffracted_photons-total_flux_emitted_photons)/total_flux_injected_photons) ' \n'])

if model=='F'
    %Verification of the maximal photon number in the last Fock state
    expect_occupation_last_Fock_state=expect(occupation_last_Fock_state,rhoss_CW); %
    fprintf(['Verification - occupation of the last Fock state : ' num2str(max(abs(expect_occupation_last_Fock_state))) ' \n \n'])
end


% Basic verifications on the intensity correlations 
% (one also has to verify that the g(2) curves are smooth - irregular curves probably indicate a wrong numerical convergence)
if min(g2_reflected_vs_delay)<0 || min(g2_transmitted_vs_delay)<0 || min(g2_emitted_vs_delay)<0 
    fprintf(' \n \n !!!!!!!!! Warning: non-physical negative values !!!!!!!!!! \n \n');
end

fprintf(['Verification: relative error on g2_reflected(0)=<b_out'' b_out'' b_out b_out>/(<b_out'' b_out>^2) : ' num2str(abs(g2_reflected_vs_delay(1) - (expect(b_out'*b_out'*b_out*b_out,rhoss_CW)/(expect(b_out'*b_out,rhoss_CW)^2) )))  ' \n'])
fprintf(['Verification: relative error on g2_transmitted(0)=<c_out'' c_out'' c_out c_out>/(<c_out'' c_out>^2) : ' num2str(abs(g2_transmitted_vs_delay(1) - (expect(c_out'*c_out'*c_out*c_out,rhoss_CW)/(expect(c_out'*c_out,rhoss_CW)^2) )))  ' \n'])
fprintf(['Verification: relative error on g2_reflected(0)=<e_out'' e_out'' e_out e_out>/(<e_out'' e_out>^2) : ' num2str(abs(g2_emitted_vs_delay(1) - (expect(e_out'*e_out'*e_out*e_out,rhoss_CW)/(expect(e_out'*e_out,rhoss_CW)^2) )))  ' \n'])

fprintf(['Verification: relative error on g2(infty)=1 for reflected photons : ' num2str(abs(g2_reflected_vs_delay(nb_points_time_g2CW)-1)) ' \n'])
fprintf(['Verification: relative error on g2(infty)=1 for transmitted photons : ' num2str(abs(g2_transmitted_vs_delay(nb_points_time_g2CW)-1)) ' \n'])
fprintf(['Verification: relative error on g2(infty)=1 for emitted photons : ' num2str(abs(g2_emitted_vs_delay(nb_points_time_g2CW)-1)) ' \n'])



%%%%%%%%%%%%%%%%%%  Plots of the normalized g(2)(tau) %%%%%%%%%%%%%%%%%%%%%%
if ismember('R',plot_choice)
    
    figure('Name',['g2CW ' text_legend_model ' - g2 vs tau  - reflected photons - Pin = ' num2str(P_in_CW_pW) ' pW - Laser photon energy = ' num2str(omega_laser_ev) ' eV' ],'NumberTitle','off')
    plot(full_tau_list_g2CW,full_g2_reflected_vs_delay,'r');
    xlabel('Delay \tau [ps]'); ylabel('g^{2} - reflected photons');
    xlim([min(full_tau_list_g2CW) max(full_tau_list_g2CW)]);
    text(full_tau_list_g2CW(ceil(nb_points_time_g2CW/20)),0.8*max(full_g2_reflected_vs_delay),text_legend,'FontSize',9)
    
    if ismember('O',plot_choice)
        figure('Name',['g2CW ' text_legend_model ' - occupation probabilities conditioned to a reflected photon detecton - Pin = ' num2str(P_in_CW_pW) ' pW - Laser photon energy = ' num2str(omega_laser_ev) ' eV' ],'NumberTitle','off')
        plot(tau_list_g2CW,occupation_ground_vs_delay_after_reflected_photon_detection,'r',...
             tau_list_g2CW,occupation_excited_vs_delay_after_reflected_photon_detection,'b');
        xlabel('Delay \tau [ps]'); ylabel('Conditionnal occupation probabilities');
        xlim([min(tau_list_g2CW) max(tau_list_g2CW)]); ylim([0 1]);
        text(tau_list_g2CW(ceil(nb_points_time_g2CW/10)),0.7,text_legend,'FontSize',9)
        legend('Occupation of |g> conditioned on a reflected photon detection','Occupation of |e> conditioned on a reflected photon detection','Location','best');
    end
end

if ismember('T',plot_choice)
    figure('Name',['g2CW ' text_legend_model ' - g2 vs tau - transmitted photons - Pin = ' num2str(P_in_CW_pW) ' pW - Pulsation laser = ' num2str(omega_laser_ev) ' eV' ],'NumberTitle','off')
    plot(full_tau_list_g2CW,full_g2_transmitted_vs_delay,'b');
    xlabel('Delay \tau [ps]'); ylabel('g^{2} - transmitted photons');
    xlim([min(full_tau_list_g2CW) max(full_tau_list_g2CW)]);
    text(full_tau_list_g2CW(ceil(nb_points_time_g2CW/20)),0.8*max(full_g2_transmitted_vs_delay),text_legend,'FontSize',9)
    
    if ismember('O',plot_choice)
        figure('Name',['g2CW ' text_legend_model ' - occupation probabilities conditioned to a transmitted photon detecton - Pin = ' num2str(P_in_CW_pW) ' pW - Laser photon energy = ' num2str(omega_laser_ev) ' eV' ],'NumberTitle','off')
        plot(tau_list_g2CW,occupation_ground_vs_delay_after_transmitted_photon_detection,'r',...
             tau_list_g2CW,occupation_excited_vs_delay_after_transmitted_photon_detection,'b');
        xlabel('Delay \tau [ps]'); ylabel('Conditionnal occupation probabilities');
        xlim([min(tau_list_g2CW) max(tau_list_g2CW)]); ylim([0 1]);
        text(tau_list_g2CW(ceil(nb_points_time_g2CW/10)),0.7,text_legend,'FontSize',9)
        legend('Occupation of |g> conditioned on a transmitted photon detection','Occupation of |e> conditioned on a transmitted photon detection','Location','best');
    end
end
if ismember('E',plot_choice)

    figure('Name',['g2CW ' text_legend_model ' - g2 vs tau - photons emitted outside the mode - Pin = ' num2str(P_in_CW_pW) ' pW - Pulsation laser = ' num2str(omega_laser_ev) ' eV' ],'NumberTitle','off')
    plot(full_tau_list_g2CW,full_g2_emitted_vs_delay,'g');
    xlabel('Delay \tau [ps]'); ylabel('g^{2} - photons emitted outside the mode');
    xlim([min(full_tau_list_g2CW) max(full_tau_list_g2CW)]);
    text(full_tau_list_g2CW(ceil(nb_points_time_g2CW/20)),0.3*max(full_g2_emitted_vs_delay),text_legend,'FontSize',9)

    if ismember('O',plot_choice)
        figure('Name',['g2CW ' text_legend_model ' - occupation probabilities conditioned to an emitted photon detecton - Pin = ' num2str(P_in_CW_pW) ' pW - Laser photon energy = ' num2str(omega_laser_ev) ' eV' ],'NumberTitle','off')
        plot(tau_list_g2CW,occupation_ground_vs_delay_after_emitted_photon_detection,'r',...
             tau_list_g2CW,occupation_excited_vs_delay_after_emitted_photon_detection,'b');
        xlabel('Delay \tau [ps]'); ylabel('Conditionnal occupation probabilities');
        xlim([min(tau_list_g2CW) max(tau_list_g2CW)]); ylim([0 1]);
        text(tau_list_g2CW(ceil(nb_points_time_g2CW/10)),0.7,text_legend,'FontSize',9)
        legend('Occupation of |g> conditioned on an emitted photon detection','Occupation of |e> conditioned on an emitted photon detection','Location','best');
    end
end