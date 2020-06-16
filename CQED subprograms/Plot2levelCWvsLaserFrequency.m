switch model
    case 'F' 
        text_legend_model='full model';
    case 'A'
        text_legend_model='adiabatic model';
end

% Parameters for the text displayed in figure legends
text_legend_QD=['g=', num2str(g_muev) ' muev    \gamma_{sp}=' num2str(gamma_sp_muev) 'muev     \gamma^*=' num2str(gamma_puredephasing_muev) 'muev'];
text_legend_cav=['\kappa=' num2str(kappa_muev) '\muev   \eta_{top}=' num2str(eta_top)  '   C=' num2str(g^2/kappa/gamma_decoherence,3)];
text_legend_P_in=['P_{in}=' num2str(P_in_CW_pW) 'pW  n_0=' num2str(4*eta_top*abs(b_in_CW)^2/kappa,2) '  n_c=' num2str(gamma_decoherence*gamma_sp/(4*g^2),2)];
text_legend={text_legend_QD,text_legend_cav,text_legend_P_in};
% NB: C is the cooperativity and n_c the critical photon number, both depend only on the cavity-QED parameters.
% On the contrary, n_0 depends on the incoming power P_in since it is the  calculated number of intracavity photons
% in the absence of QD (i.e. when g=0). When n_0 is much lower than n_c we are in the weak excitation limit.


% Verification of the conservation of total photon flux
fprintf(['CW - ' text_legend_model ' ' ': Maximal relative error on the conservation of photon flux: ' num2str(max(abs(1-R_vs_omega-T_vs_omega-D_vs_omega-E_vs_omega))) ' \n \n'])

%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%


if ismember('R',plot_choice)
    %Color code : Refl in red, Transm + Diffr/lost in blue, Spont. Em. in magenta
    figure('Name',['CW - ' text_legend_model ' ' ' - Reflectivity vs laser photon energy - Pin = ' num2str(P_in_CW_pW) ' pW'],'NumberTitle','off')
    %%% UNUSED HERE: allows comparing total part and coherent part
    % plot(omega_laser_list_ev_CW,real(R_vs_omega),'r',omega_laser_list_ev_CW,real(R_coh_vs_omega),'r--');
    % legend('Fraction of reflected photons','Fraction of coherent reflected photons') 
    plot(detuning_list_muev,real(R_vs_omega),'r');
    xlabel('\omega_{laser}-\omega_{d} [\mueV]'); ylabel('Reflectivity');
    ylim([0 1]);
    text(detuning_list_muev(ceil(nb_points_spectrum*0.55)),0.8,text_legend,'FontSize',9)
    legend('Reflected photons','Location','Northeast') 
end


if ismember('T',plot_choice)
    figure('Name',['CW - ' text_legend_model ' ' ' - Transmission vs laser photon energy - Pin = ' num2str(P_in_CW_pW) ' pW'],'NumberTitle','off')
    %%% UNUSED HERE: allows comparing total part and coherent part
    % plot(omega_laser_list_ev_CW,real(T_vs_omega+D_vs_omega),'b',omega_laser_list_ev_CW,real(T_coh_vs_omega+D_coh_vs_omega),'b--');
    % legend('Fraction of transmitted + diffracted photons','Fraction of coherent transmitted/diffracted photons')
    plot(detuning_list_muev,real(T_vs_omega+D_vs_omega),'b');
    
    xlabel('\omega_{laser}-\omega_{d} [\mueV]'); ylabel('Transmission + diffraction/losses');
    ylim([0 1]);
    text(detuning_list_muev(ceil(nb_points_spectrum*0.55)),0.12,text_legend,'FontSize',9)
    legend('Transmitted + diffracted photons','Location','Northeast') 
end


if ismember('E',plot_choice)
    figure('Name',['CW - ' text_legend_model ' ' ' - Spontaneous emission outside the mode vs laser photon energy - Pin = ' num2str(P_in_CW_pW) ' pW'],'NumberTitle','off')
    % %%% UNUSED HERE: allows comparing total part and coherent part
    % plot(omega_laser_list_ev_CW,real(E_vs_omega),'m',omega_laser_list_ev_CW,real(E_coh_vs_omega),'m--');
    % legend('Fraction of spontaneously emitted photons','Fraction of coherent spontaneously-emitted photons') 
    plot(detuning_list_muev,real(E_vs_omega),'m');
%     xlabel('\omega_{laser}-\omega_{d} [\mueV]'); ylabel('Fraction of photons emitted outside the mode'); 
    ylim([0 1]);
    text(detuning_list_muev(ceil(nb_points_spectrum*0.55)),0.8,text_legend,'FontSize',9)
    legend('Photons emitted outside the mode','Location','Northeast')
end
    


if ismember('O',plot_choice)
    figure('Name',['CW - ' text_legend_model ' ' ' - Occupation probabilities vs laser photon energy - Pin = ' num2str(P_in_CW_pW) ' pW'],'NumberTitle','off')
    plot(detuning_list_muev,real(occupation_excited_state_vs_omega),'r',...
         detuning_list_muev,real(occupation_ground_state_vs_omega),'b');
    xlabel('\omega_{laser}-\omega_{d} [\mueV]'); ylabel('Occupation probability'); 
    ylim([-0.05 1.05]);
    text(detuning_list_muev(ceil(nb_points_spectrum*0.55)),0.7,text_legend,'FontSize',9)
    legend('Occupation state |e>','Occupation state |g>','Location','best')
end