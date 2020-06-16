switch model
    case 'F'
        %%%%%%%%%%%%%%%%%%%% Sub-space "quantum dot"   %%%%%%%%%%%%%%%%%%%%

        %Basis states for the quantum dot (g for "ground", e for "excited")
        g_ket=qo([1;0]); %%% Quantum object "ket" associated to the ground state
        g_bra=g_ket'; %%% Quantum object "bra" associated to the ground state
        e_ket=qo([0;1]); %%% Quantum object "ket" associated to the excited state
        e_bra=e_ket';   %%% Quantum object "bra" associated to the excited state

        %Operators acting in the sub-space of the 2-level QD system (2x2 matrixs)
        id_QD=e_ket*e_bra+g_ket*g_bra; % Quantum object "Identity operator"
        sigma_QD=g_ket*e_bra; % Quantum object "de-excitation of the 2-level system"
        sigma_dag_QD=sigma_QD'; % Quantum object "excitation of the 2-level system"
       
        %%%%%%%%%%%%%%%%%%%% Sub-space "cavity"   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %N denotes the maximal number of photons in the truncated Fock space (see above)
        id_cav = identity(N);  % Quantum object "Identity operator"  (matrix NxN)
        destroy_cav=destroy(N);% Quantum object "annihilation operator"  (matrix NxN)
        Vacuum_state=basis(N,1); % State corresponding to photon vacuum |0> (Nx1)
        
        last_Fock_state_ket = basis(N,N);% State corresponding to the last Fock state |N-1>, to control its occupation probability 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%% Operators acting in the product space cavity-quantum dot (matrices 2Nx2N)
        Id=tensor(id_cav,id_QD); % Tensorial product of identities in the cavity and QD subspace
        sigma = tensor(id_cav,sigma_QD); % De-excitation operator from the excited QD state to the ground QD state
        sigma_dag = sigma'; % De-excitation operator from the ground QD state to the excited QD state
        
        a=tensor(destroy_cav,id_QD); % Annihilation operator for a photon in the cavity mode
        a_dag=a'; % Creation operator for a photon in the cavity mode
        n=a_dag*a; % Intracavity photon number operator

        occupation_last_Fock_state=tensor(last_Fock_state_ket*last_Fock_state_ket',id_QD);% operator for the occupation of the last Fock state, to control that it remains low
        
        
        %%%%%%%%%%%%%%%% Definition of Collapse operators %%%%%%%%%%%%%%%%%%%%%
        C_cav = sqrt(kappa)*a; %Collapse operator for the cavity
        C_sp = sqrt(gamma_sp)*sigma; % Collapse operator for a spontaneous emission outside the cavity mode
        C_pure_dephasing = sqrt(2*gamma_puredephasing)*sigma_dag*sigma; %Collapse operator for pure dephasing

        % Lindblad operators associated to each incoherent process
        L_cav = 1/2 * (2*spre(C_cav)*spost(C_cav') - spre(C_cav'*C_cav) - spost(C_cav'*C_cav));% 1st terme: cavity dumping
        L_sp = 1/2 * (2*spre(C_sp)*spost(C_sp') - spre(C_sp'*C_sp) - spost(C_sp'*C_sp));% 2nd terme: exciton lifetime
        L_pure_dephasing = 1/2 * (2*spre(C_pure_dephasing)*spost(C_pure_dephasing') - spre(C_pure_dephasing'*C_pure_dephasing) - spost(C_pure_dephasing'*C_pure_dephasing));% 3rd term: pure dephasing
        
        % Lindblad operator associated to all the incoherent processes
        L_incoh=L_cav+L_sp+L_pure_dephasing;
        
    case 'A'
        %%%%%%%%%%%%%%%%%%%% space "quantum dot"   %%%%%%%%%%%%%%%%%%%%

        %Basis states for the quantum dot (g for "ground", e for "excited")
        g_ket=qo([1;0]); %%% Quantum object "ket" associated to the ground state
        g_bra=g_ket'; %%% Quantum object "bra" associated to the ground state
        e_ket=qo([0;1]); %%% Quantum object "ket" associated to the excited state
        e_bra=e_ket';   %%% Quantum object "bra" associated to the excited state

        %Operators acting in the sub-space of the 2-level QD system (2x2 matrixs)
        Id=e_ket*e_bra+g_ket*g_bra; % Quantum object "Identity operator"
        sigma=g_ket*e_bra; % Quantum object "de-excitation of the 2-level system"
        sigma_dag=sigma'; % Quantum object "excitation of the 2-level system"
       
        %%%%%%%%%%%%%%%% Definition of Collapse operators %%%%%%%%%%%%%%%%%%%%%
        C_QD = sqrt(Gamma_tot)*sigma; % Collapse operator for a spontaneous emission outside the cavity mode (eq.15)
        C_pure_dephasing = sqrt(2*gamma_puredephasing)*sigma_dag*sigma; %Collapse operator for pure dephasing

        % Lindblad operators associated to each incoherent process
        L_QD = 1/2 * (2*spre(C_QD)*spost(C_QD') - spre(C_QD'*C_QD) ...
            - spost(C_QD'*C_QD));% exciton lifetime
        L_pure_dephasing = 1/2 * (2*spre(C_pure_dephasing)*spost(C_pure_dephasing') - spre(C_pure_dephasing'*C_pure_dephasing) - spost(C_pure_dephasing'*C_pure_dephasing));% 3rd term: pure dephasing
               
        % Lindblad operator associated to all the incoherent processes
        L_incoh=L_QD+L_pure_dephasing;
        
        % NB: in the adiabatic model L_QD includes both the emission outside the mode 
        % and the Purcell-enhanced emission through the cavity mode, hence the 
        % "Gamma_tot" term in the definition of C_QD
        
end

%%%%%%%%%%%  Constant output operator used to describe the field emitted outside the mode 
e_out=sqrt(gamma_sp)*sigma; % Output flux operator in ps^(-1/2)

% NB: the other output flux operators can be power-dependent or frequency-dependent, depending on the model used 
% (full or adiabatic). They are thus defined in the main file. Also note that e_out is equal to C_sp, but we use
% a different notation to insist on its use as an output operator, analogous to the other output operators b_out,
% c_out, and d_out, respectively describing the reflected, transmitted, and diffracted/lost photon field.