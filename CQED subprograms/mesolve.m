% Function developed by Sze Meng Tan. A quantum optics toolbox, 1999
function [rho] = mesolve(L,rho0,t_list)
%%%%%%%%%%%%%%%%%% L : Liouvillien use in the Master Equation
%%%%%%%%%%%%%%%%%% rho 0 : Initial Condition for the solving of the ME
%%%%%%%%%%%%%%%%%% t_list : time list containing each time value where the
%%%%%%%%%%%%%%%%%% Density Matrix will be calculated thanks to the routine

%%%%%%%%%%% Numerical Integration of the Master Equation %%%%%%%%


ode2file('ode_input.dat',L,rho0,t_list,struct('reltol',7e-8,'abstol',8e-7)); % Writes the data into a file
odesolve('ode_input.dat','ode_output.dat'); % Solve/Integration (here Adams method by default)
fid = fopen('ode_output.dat','rb');
rho = qoread(fid,dims(rho0),size(t_list));  
fclose(fid);
%%%%%%%%%%%%%%%%%% rho : list of the Density Matrix calculated for each
%%%%%%%%%%%%%%%%%% time specified in the list





