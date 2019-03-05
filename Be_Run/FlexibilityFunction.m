
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Methode Flexibility Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Eerst twee simulaties uitvoeren met BuiSim
% 1) Run MPC with constant cost signal = minimum heat demand (zoals de
%    toolbox standaard ingesteld staat door Jan) = Q_REF
% 2) Run MPC  with step in price signal = Q_STEP

% Load full reference profile of one year
load('OutputOldEnergy.mat')
load('OutputOldElectricity.mat')
Q_REF_Full_Energy = sum(RefEnergy); % sum(A) telt alle kolomelementen van matrix A op resulterend in 1 rijmatrix
Q_REF_Full_Electricity = sum(RefElectricity);


% Load step response for desired amount of days
BuiInit

% Simulation parameters
SimStart = outdata.SimParam.SimStart;
SimStop = outdata.SimParam.SimStop;
Nsim = outdata.SimParam.run.Nsim; 
Ts = model.plant.Ts;
Time = (1:Nsim)*outdata.model.plant.Ts/3600/24;

%% Energy

% Reference profile and step response profile
Q_REF_Energy = Q_REF_Full_Energy(SimStart:SimStop);
Q_STEP_Energy = sum(outdata.data.U);

% Flexibility Function
FF_Energy = Q_STEP_Energy - Q_REF_Energy;

figure
subplot(3,1,1); 
plot(Time, Q_REF_Energy, 'linewidth', 2);
title('Reference heat demand', 'fontsize', 14)
xlabel('time [days]')
ylabel('Heating power [W]')
grid on
subplot(3,1,2); 
plot(Time, Q_STEP_Energy, 'linewidth', 2);
title('Heating demand after step function', 'fontsize', 14)
xlabel('time [days]')
ylabel('Heating power [W]')
grid on
subplot(3,1,3);
plot(Time, FF_Energy, 'linewidth', 2);
title('Flexibility Function', 'fontsize', 14);
xlabel('time [days]')
ylabel('Heating power [W]')
grid on
 
% ------------------------------------------------------------
% Reduced energy use due to step increase of penalty
% Increased energy use in anticipation of penalty increase
% ------------------------------------------------------------
Reduced_EU = 0; 
Increased_EU = 0; 
    for j = 1:size(FF_Energy,2)
        if FF_Energy(j)<0
            Reduced_EU = Reduced_EU + FF_Energy(j)*Ts/(3.6*10^6);
            %disp(Reduced_EU)
        else
            Increased_EU = Increased_EU + FF_Energy(j)*Ts/(3.6*10^6);
            %disp(Increased_EU)
        end
    end

fprintf('Increased energy use : %.2f kWh\n', Increased_EU)
fprintf('Reduced energy use: %.2f kWh\n', Reduced_EU)

% % --------------------------------------------------
% % Maximum power increase & Maximum power reduction
% % -------------------------------------------------- 

Max_power_red = min(FF_Energy,[],'all');
Max_power_inc = max(FF_Energy,[],'all');    
 
fprintf('Maximum power increase: %.2f W\n', Max_power_inc)
fprintf('Maximum power reduction : %.2f W\n', Max_power_red)



%% Electricity

% Reference profile and step response profile
Q_REF_Electricity = Q_REF_Full_Electricity(SimStart:SimStop);
Q_STEP_Electricity = sum(outdata.data.E);

% Flexibility Function
FF_Electricity = Q_STEP_Electricity - Q_REF_Electricity;

figure
subplot(3,1,1); 
plot(Time, Q_REF_Electricity, 'linewidth', 2);
title('Reference electricity demand', 'fontsize', 14)
xlabel('time [days]')
ylabel('Electric power [W]')
grid on
subplot(3,1,2); 
plot(Time, Q_STEP_Electricity, 'linewidth', 2);
title('Electricity demand after step function', 'fontsize', 14)
xlabel('time [days]')
ylabel('Electric power [W]')
grid on
subplot(3,1,3);
plot(Time, FF_Electricity, 'linewidth', 2);
title('Flexibility Function', 'fontsize', 14);
xlabel('time [days]')
ylabel('Electric power [W]')
grid on
% 
% start_time_penalty = t1; % Tijdslot wanneer penalty signaal verhoogd is
% 
% ------------------------------------------------------------
% Reduced energy use due to step increase of penalty
% Increased energy use in anticipation of penalty increase
% ------------------------------------------------------------
Reduced_ElecU = 0; 
Increased_ElecU = 0; 
    for j = 1:size(FF_Electricity,2)
        if FF_Electricity(j)<0
            Reduced_ElecU = Reduced_ElecU + FF_Electricity(j)*Ts/(3.6*10^6);
            %disp(Reduced_ElecU)
        else
            Increased_ElecU = Increased_ElecU + FF_Electricity(j)*Ts/(3.6*10^6);
            %disp(Increased_ElecU)
        end
    end

fprintf('Increased electricity use : %.2f kWh\n', Increased_ElecU)
fprintf('Reduced electricity use: %.2f kWh\n', Reduced_ElecU)

% % --------------------------------------------------
% % Maximum power increase & Maximum power reduction
% % -------------------------------------------------- 

Max_elec_power_red = min(FF_Electricity,[],'all');
Max_elec_power_inc = max(FF_Electricity,[],'all');    
 
fprintf('Maximum electric power increase: %.2f W\n', Max_elec_power_inc)
fprintf('Maximum electric power reduction : %.2f W\n', Max_elec_power_red)