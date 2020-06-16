%Verification of the photon number conservation
fprintf(['RFTR: relative error of the photon number conservation: ' num2str(abs(Nb_reflected_photons+Nb_transmitted_photons+Nb_diffracted_photons+Nb_emitted_photons-Nb_photons_pulse)/Nb_photons_pulse) ' \n'])
%
switch model
    case 'F' 
        text_legend_model='full model';
    case 'A'
        text_legend_model='adiabatic model';
end
%% plotting g1(t1,t2)

if ismember('g',plot_choice)
 
    figure('Name',['g1 reflected vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
    surf(flip(t_list),t_list,flip(abs(g1_reflected_vs_t1_t2),2),'LineStyle','none')
    xlabel('t_1 [ps]')
    ylabel('t_2 [ps]')
    title('reflected |g^{(1)}(t_1,t_2)|')
    view(2)
    colorbar

    figure('Name',['g1 diffracted vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
    surf(flip(t_list),t_list,flip(abs(g1_diffracted_vs_t1_t2),2),'LineStyle','none')
    xlabel('t_1 [ps]')
    ylabel('t_2 [ps]')
    title('diffracted |g^{(1)}(t_1,t_2)|')
    view(2)
    colorbar
    
    figure('Name',['g1 transmitted vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
    surf(flip(t_list),t_list,flip(abs(g1_transmitted_vs_t1_t2),2),'LineStyle','none')
    xlabel('t_1 [ps]')
    ylabel('t_2 [ps]')
    title('transmitted |g^{(1)}(t_1,t_2)|')
    view(2)
    colorbar
    
    figure('Name',['g1 emitted vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
    surf(flip(t_list),t_list,flip(abs(g1_emitted_vs_t1_t2),2),'LineStyle','none')
    xlabel('t_1 [ps]')
    ylabel('t_2 [ps]')
    title('emitted |g^{(1)}(t_1,t_2)|')
    view(2)
    colorbar
end
%% plottig G1 vs(t1,t2)
if ismember('G',plot_choice)

    figure('Name',['G1 reflected vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
    surf(flip(t_list),t_list,flip(abs(G1_reflected_vs_t1_t2),2),'LineStyle','none')
    xlabel('t_1 [ps]')
    ylabel('t_2 [ps]')
    title('|<b_{out}''(t_1) b_{out}(t_2)>|')
    view(2)
    colorbar

    figure('Name',['G1 transmitted vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
    surf(flip(t_list),t_list,flip(abs(G1_transmitted_vs_t1_t2),2),'LineStyle','none')
    xlabel('t_1 [ps]')
    ylabel('t_2 [ps]')
    title('|<c_{out}''(t_1) c_{out}(t_2)>|')
    view(2)
    colorbar
    
    figure('Name',['G1 diffracted vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
    surf(flip(t_list),t_list,flip(abs(G1_diffracted_vs_t1_t2),2),'LineStyle','none')
    xlabel('t_1 [ps]')
    ylabel('t_2 [ps]')
    title('|<d_{out}''(t_1) d_{out}(t_2)>|')
    view(2)
    colorbar

    figure('Name',['G1 emitted vs (t1,t2) - ' text_legend_model],'NumberTitle','off')
    surf(flip(t_list),t_list,flip(abs(G1_emitted_vs_t1_t2),2),'LineStyle','none')
    xlabel('t_1 [ps]')
    ylabel('t_2 [ps]')
    title('|<e_{out}''(t_1) e_{out}(t_2)>|')
    view(2)
    colorbar
end
% Parameters for the text displayed in figure legends
text_legend_QD = ['g=', num2str(g_muev) ' muev    \gamma_{sp}=' num2str(gamma_sp_muev) 'muev     \gamma^*=' num2str(gamma_puredephasing_muev) 'muev'];
text_legend_cav = ['\kappa=' num2str(kappa_muev) 'muev   \eta_{top}=' num2str(eta_top)  '   C=' num2str(g^2/kappa/gamma_decoherence,3)];
text_legend_pulse = ['N_{in}=' num2str(Nb_photons_pulse) '       FWHM_{pulse}=' num2str(FWHM_pulse) 'ps'];
text_legend = {text_legend_QD,text_legend_cav,text_legend_pulse};
% NB: C is the cooperativity 

%% plotting interpolated G1(time,tau)
if ismember('I',plot_choice)
    figure('Name',['interpolated G1 reflected vs time tau - ' text_legend_model ' - Nb_photons = ' num2str(Nb_photons_pulse)],'NumberTitle','off')
    surf(time_list,tau_list, abs(interpolated_G1_reflected_vs_time_tau'),'LineStyle','none')
    axis xy
    colorbar
    xlabel('t, origin at beginning of 1st evolution [ps]','FontSize',14)
    ylabel('\tau [ps]','FontSize',14)
    title('interpolated G1 reflected vs time tau')
    view(2)

    figure('Name',['interpolated G1 transmitted vs time tau - ' text_legend_model ' - Nb_photons = ' num2str(Nb_photons_pulse)],'NumberTitle','off')
    surf(time_list,tau_list, abs(interpolated_G1_transmitted_vs_time_tau'),'LineStyle','none')
    axis xy
    colorbar
    xlabel('t, origin at beginning of 1st evolution [ps]','FontSize',14)
    ylabel('\tau [ps]','FontSize',14)
    title('interpolated G1 transmitted vs time tau')
    view(2)

    figure('Name',['interpolated G1 emitted vs time tau - ' text_legend_model ' - Nb_photons = ' num2str(Nb_photons_pulse)],'NumberTitle','off')
    surf(time_list,tau_list, abs(interpolated_G1_emitted_vs_time_tau'),'LineStyle','none')
    axis xy
    colorbar
    xlabel('t, origin at beginning of 1st evolution [ps]','FontSize',14)
    ylabel('\tau [ps]','FontSize',14)
    title('interpolated G1 emitted vs time tau')
    view(2)

    figure('Name',['interpolated G1 diffracted vs time tau - ' text_legend_model ' - Nb_photons = ' num2str(Nb_photons_pulse)],'NumberTitle','off')
    surf(time_list,tau_list, abs(interpolated_G1_diffracted_vs_time_tau'),'LineStyle','none')
    axis xy
    colorbar
    xlabel('t, origin at beginning of 1st evolution [ps]','FontSize',14)
    ylabel('\tau [ps]','FontSize',14)
    title('interpolated G1 diffracted vs time tau')
    view(2)
end
%% Plotting WDF
if ismember('W',plot_choice)
    figure('Name',['WDF reflected photons - ' text_legend_model ' - Nb_photons = ' num2str(Nb_photons_pulse)],'NumberTitle','off')
    imagesc(time_list,(omega_list_full_spectrum_muev - omega_pulse_ev*1e6),  real(WDF_interpolated_G1_reflected_vs_time_tau_full_spectrum)');
    axis xy
    colorbar
    xlabel('t, origin at beginning of 1st evolution [ps]')
    ylabel('\omega-\omega_p [\mueV]','FontSize',14)
    title('Wigner distribution reflected [1/(\mueV\cdotps)]')

    figure('Name',['WDF reflected photons - zoomed spectrum - ' text_legend_model ' - Nb_photons = ' num2str(Nb_photons_pulse)],'NumberTitle','off')
    imagesc(time_list,(omega_list_ev - omega_pulse_ev)*1e6,  real(WDF_interpolated_G1_reflected_vs_time_tau'));
    axis xy
    colorbar
    xlabel('t, origin at beginning of 1st evolution [ps]')
    ylabel('\omega-\omega_p [\mueV]','FontSize',14)
    title('zoomed Wigner distribution reflected [1/(\mueV\cdotps)]')

    figure('Name',['WDF transmitted photons - ' text_legend_model ' - Nb_photons = ' num2str(Nb_photons_pulse)],'NumberTitle','off')
    imagesc(time_list,(omega_list_full_spectrum_muev - omega_pulse_ev*1e6),  real(WDF_interpolated_G1_transmitted_vs_time_tau_full_spectrum)');
    axis xy
    colorbar
    xlabel('t, origin at beginning of 1st evolution [ps]')
    ylabel('\omega-\omega_p [\mueV]','FontSize',14)
    title('Wigner distribution transmitted [1/(\mueV\cdotps)]')

    figure('Name',['WDF transmitted photons - zoomed spectrum -' text_legend_model ' - Nb_photons = ' num2str(Nb_photons_pulse)],'NumberTitle','off')
    imagesc(time_list,(omega_list_ev - omega_pulse_ev)*1e6,  real(WDF_interpolated_G1_transmitted_vs_time_tau'));
    axis xy
    colorbar
    xlabel('t, origin at beginning of 1st evolution [ps]')
    ylabel('\omega-\omega_p [\mueV]','FontSize',14)
    title('zoomed Wigner distribution transmitted [1/(\mueV\cdotps)]')

    figure('Name',['WDF emitted photons - ' text_legend_model ' - Nb_photons = ' num2str(Nb_photons_pulse)],'NumberTitle','off')
    imagesc(time_list,(omega_list_full_spectrum_muev - omega_pulse_ev*1e6),  real(WDF_interpolated_G1_emitted_vs_time_tau_full_spectrum)');
    axis xy
    colorbar
    xlabel('t, origin at beginning of 1st evolution [ps]')
    ylabel('\omega-\omega_p [\mueV]','FontSize',14)
    title('Wigner distribution emitted [1/(\mueV\cdotps)]')

    figure('Name',['WDF emitted photons - zoomed spectrum - ' text_legend_model ' - Nb_photons = ' num2str(Nb_photons_pulse)],'NumberTitle','off')
    imagesc(time_list,(omega_list_ev - omega_pulse_ev)*1e6,  real(WDF_interpolated_G1_emitted_vs_time_tau'));
    axis xy
    colorbar
    xlabel('t, origin at beginning of 1st evolution [ps]')
    ylabel('\omega-\omega_p [\mueV]','FontSize',14)
    title('zoomed Wigner distribution emitted [1/(\mueV\cdotps)]')

    figure('Name',['WDF diffracted photons - ' text_legend_model ' - Nb_photons = ' num2str(Nb_photons_pulse)],'NumberTitle','off')
    imagesc(time_list,(omega_list_full_spectrum_muev - omega_pulse_ev*1e6),  real(WDF_interpolated_G1_diffracted_vs_time_tau_full_spectrum)');
    axis xy
    colorbar
    xlabel('t, origin at beginning of 1st evolution [ps]')
    ylabel('\omega-\omega_p [\mueV]','FontSize',14)
    title('Wigner distribution diffracted  [1/(\mueV\cdotps)]')

    figure('Name',['WDF diffracted photons - zoomed spectrum - ' text_legend_model ' - Nb_photons = ' num2str(Nb_photons_pulse)],'NumberTitle','off')
    imagesc(time_list,(omega_list_ev - omega_pulse_ev)*1e6,  real(WDF_interpolated_G1_diffracted_vs_time_tau'));
    axis xy
    colorbar
    xlabel('t, origin at beginning of 1st evolution [ps]')
    ylabel('\omega-\omega_p [\mueV]','FontSize',14)
    title('zoomed Wigner distribution diffracted [1/(\mueV\cdotps)]')
end
%% Plotting photon fluxes
if ismember('F',plot_choice)
    figure('Name',['Flux reflected photons vs time - Nb_photons = ' num2str(Nb_photons_pulse)],'NumberTitle','off') 
    plot(t_list, flux_reflected_photons_vs_time,'r', time_list, real(sum(WDF_interpolated_G1_reflected_vs_time_tau_full_spectrum,2))*omega_step_muev,'g--', t_list, flux_reflected_photons_pulse_coherent_vs_time,'k',time_list, real(sum(WDF_interpolated_reflected_coherent_vs_time_tau_full_spectrum,2))*omega_step_muev,'m--');
    xlabel('Time t [ps]');ylabel('Flux reflected photon [ps^{-1}]');
    legend('<b_{out}''b_{out}>','total flux by WDF(t,\omega)','<b_{out}''><b_{out}>','total coherent flux by WDF(t,\omega)');
    text(t_list(floor(nb_points_time*0.6)),0.2*max(flux_reflected_photons_vs_time),text_legend,'FontSize',9)

    figure('Name',['Flux transmitted photons vs time - Nb_photons = ' num2str(Nb_photons_pulse)],'NumberTitle','off') 
    plot(t_list, flux_transmitted_photons_vs_time,'r', time_list, real(sum(WDF_interpolated_G1_transmitted_vs_time_tau_full_spectrum,2))*omega_step_muev,'g--', t_list, flux_transmitted_photons_pulse_coherent_vs_time,'k',time_list, real(sum(WDF_interpolated_transmitted_coherent_vs_time_tau_full_spectrum,2))*omega_step_muev,'m--');
    xlabel('Time t [ps]');ylabel('Flux transmitted photon [ps^{-1}]');
    legend('<c_{out}''c_{out}>','total flux by WDF(t,\omega)','<c_{out}''><c_{out}>','total coherent flux by WDF(t,\omega)');
    text(t_list(floor(nb_points_time*0.6)),0.2*max(flux_transmitted_photons_vs_time),text_legend,'FontSize',9)

    figure('Name',['Flux emitted photons vs time - Nb_photons = ' num2str(Nb_photons_pulse)],'NumberTitle','off') 
    plot(t_list, flux_emitted_photons_vs_time,'r', time_list, real(sum(WDF_interpolated_G1_emitted_vs_time_tau_full_spectrum,2))*omega_step_muev,'g--', t_list, flux_emitted_photons_pulse_coherent_vs_time,'k',time_list, real(sum(WDF_interpolated_emitted_coherent_vs_time_tau_full_spectrum,2))*omega_step_muev,'m--');
    xlabel('Time t [ps]');ylabel('Flux emitted photon [ps^{-1}]');
    legend('<d_{out}''d_{out}>','total flux by WDF(t,\omega)','<d_{out}''><d_{out}>','total coherent flux by WDF(t,\omega)');
    text(t_list(floor(nb_points_time*0.6)),0.2*max(flux_emitted_photons_vs_time),text_legend,'FontSize',9)

    figure('Name',['Flux diffracted photons vs time - Nb_photons = ' num2str(Nb_photons_pulse)],'NumberTitle','off') 
    plot(t_list, flux_diffracted_photons_vs_time,'r', time_list, real(sum(WDF_interpolated_G1_diffracted_vs_time_tau_full_spectrum,2))*omega_step_muev,'g--', t_list, flux_diffracted_photons_pulse_coherent_vs_time,'k',time_list, real(sum(WDF_interpolated_diffracted_coherent_vs_time_tau_full_spectrum,2))*omega_step_muev,'m--');
    xlabel('Time t [ps]');ylabel('Flux diffracted photon [ps^{-1}]');
    legend('<e_{out}''e_{out}>','total flux by WDF(t,\omega)','<e_{out}''><e_{out}>','total coherent flux by WDF(t,\omega)');
    text(t_list(floor(nb_points_time*0.6)),0.2*max(flux_diffracted_photons_vs_time),text_legend,'FontSize',9)


    figure('Name',['Spectral density reflected photons vs omega - Nb_photons = ' num2str(Nb_photons_pulse)],'NumberTitle','off')
    plot((omega_list_ev - omega_pulse_ev)*1e6, real(ESD_reflected_photons_vs_omega) ,'r',(omega_list_ev - omega_pulse_ev)*1e6, real(ESD_coherent_reflected_photons_laser_vs_omega),'g',(omega_list_ev - omega_pulse_ev)*1e6, real(ESD_reflected_photons_vs_omega)-real(ESD_coherent_reflected_photons_laser_vs_omega),'b');
    xlim([min(omega_list_ev - omega_pulse_ev)*1e6 max(omega_list_ev - omega_pulse_ev)*1e6]);
    xlabel('\omega-\omega_p [\mueV]');ylabel('Reflected field spectral density');
    legend('WDF(\omega) - total reflected flux', 'WDF(\omega) coherent flux', 'WDF(\omega) - incoherent flux')
    text(omega_list_ev(1) - omega_pulse_ev,0.7*max(ESD_reflected_photons_vs_omega),text_legend,'FontSize',9)

    figure('Name',['Spectral density transmitted photons vs omega - Nb_photons = ' num2str(Nb_photons_pulse)],'NumberTitle','off')
    plot((omega_list_ev - omega_pulse_ev)*1e6, real(ESD_transmitted_photons_vs_omega) ,'r',(omega_list_ev - omega_pulse_ev)*1e6, real(ESD_coherent_transmitted_photons_laser_vs_omega),'g',(omega_list_ev - omega_pulse_ev)*1e6, real(ESD_transmitted_photons_vs_omega)-real(ESD_coherent_transmitted_photons_laser_vs_omega),'b');
    xlim([min(omega_list_ev - omega_pulse_ev)*1e6 max(omega_list_ev - omega_pulse_ev)*1e6]);
    xlabel('\omega-\omega_p [\mueV]');ylabel('Reflected field spectral density');
    legend('WDF(\omega) - total transmitted flux', 'WDF(\omega) coherent flux', 'WDF(\omega) - incoherent flux')
    text(omega_list_ev(1) - omega_pulse_ev,0.7*max(ESD_transmitted_photons_vs_omega),text_legend,'FontSize',9)

    figure('Name',['Spectral density emitted photons vs omega - Nb_photons = ' num2str(Nb_photons_pulse)],'NumberTitle','off')
    plot((omega_list_ev - omega_pulse_ev)*1e6, real(ESD_emitted_photons_vs_omega) ,'r',(omega_list_ev - omega_pulse_ev)*1e6, real(ESD_coherent_emitted_photons_laser_vs_omega),'g',(omega_list_ev - omega_pulse_ev)*1e6, real(ESD_emitted_photons_vs_omega)-real(ESD_coherent_emitted_photons_laser_vs_omega),'b');
    xlim([min(omega_list_ev - omega_pulse_ev)*1e6 max(omega_list_ev - omega_pulse_ev)*1e6]);
    xlabel('\omega-\omega_p [\mueV]');ylabel('Reflected field spectral density');
    legend('WDF(\omega) - total emitted flux', 'WDF(\omega) coherent flux', 'WDF(\omega) - incoherent flux')
    text(omega_list_ev(1) - omega_pulse_ev,0.7*max(ESD_emitted_photons_vs_omega),text_legend,'FontSize',9)

    figure('Name',['Spectral density diffracted photons vs omega - ' text_legend_model ' - Nb_photons = ' num2str(Nb_photons_pulse)],'NumberTitle','off')
    plot((omega_list_ev - omega_pulse_ev)*1e6, real(ESD_diffracted_photons_vs_omega) ,'r',(omega_list_ev - omega_pulse_ev)*1e6, real(ESD_coherent_diffracted_photons_laser_vs_omega),'g',(omega_list_ev - omega_pulse_ev)*1e6, real(ESD_diffracted_photons_vs_omega)-real(ESD_coherent_diffracted_photons_laser_vs_omega),'b');
    xlim([min(omega_list_ev - omega_pulse_ev)*1e6 max(omega_list_ev - omega_pulse_ev)*1e6]);
    xlabel('\omega-\omega_p [\mueV]');ylabel('Reflected field spectral density');
    legend('WDF(\omega) - total diffracted flux', 'WDF(\omega) coherent flux', 'WDF(\omega) - incoherent flux')
    text(omega_list_ev(1) - omega_pulse_ev,0.7*max(ESD_diffracted_photons_vs_omega),text_legend,'FontSize',9)
end

%% Checking normalization of ESDs
fprintf(['Relative error over normalization of spectral density - reflected field: ' num2str(abs(sum(ESD_reflected_photons_vs_omega)*omega_step_muev-Nb_reflected_photons)/Nb_reflected_photons) ' \n'])
fprintf(['Relative error over normalization of spectral density - transmitted field: ' num2str(abs(sum(ESD_transmitted_photons_vs_omega)*omega_step_muev-Nb_transmitted_photons)/Nb_transmitted_photons) ' \n'])
fprintf(['Relative error over normalization of spectral density - emitted field: ' num2str(abs(sum(ESD_emitted_photons_vs_omega)*omega_step_muev-Nb_emitted_photons)/Nb_emitted_photons) ' \n'])
fprintf(['Relative error over normalization of spectral density - diffracted field: ' num2str(abs(sum(ESD_diffracted_photons_vs_omega)*omega_step_muev-Nb_diffracted_photons)/Nb_diffracted_photons) ' \n'])
