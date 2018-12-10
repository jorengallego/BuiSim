
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Methode Flexibility Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Eerst twee simulaties uitvoeren met BuiSim
% 1) Run MPC with constant cost signal = minimum heat demand (zoals de
%    toolbox standaard ingesteld staat door Jan) = Q_REF
% 2) Run MPC  with step in price signal (nog bekijken met
%    Jan/Bram/Anke/Alessia hoe dit praktisch kan) = Q_STEP

% Load full reference profile of one year
load('output1.mat')
Q_REF_Full = sum(RefHeat); % sum(A) telt alle kolomelementen van matrix A op resulterend in 1 rijmatrix

% Load step response for desired amount of days
BuiInit

% Simulation parameters
SimStart = outdata.SimParam.SimStart;
SimStop = outdata.SimParam.SimStop;
Nsim = outdata.SimParam.run.Nsim; 
Ts = model.plant.Ts;
Time = (1:Nsim)*outdata.model.plant.Ts/3600/24;

% Reference profile and step response profile
Q_REF = Q_REF_Full(SimStart:SimStop);
Q_STEP = sum(outdata.data.U);

% Flexibility Function
FF = Q_STEP - Q_REF;

figure
subplot(3,1,1); 
plot(Time, Q_REF, 'linewidth', 2);
title('Reference heat demand', 'fontsize', 14)
xlabel('time [days]')
ylabel('Heating power [W]')
subplot(3,1,2); 
plot(Time, Q_STEP, 'linewidth', 2);
title('Heating demand after step function', 'fontsize', 14)
xlabel('time [days]')
ylabel('Heating power [W]')
subplot(3,1,3);
plot(Time, FF, 'linewidth', 2);
title('Flexibility Function', 'fontsize', 14);
xlabel('time [days]')
ylabel('Heating power [W]')
% 
% start_time_penalty = t1; % Tijdslot wanneer penalty signaal verhoogd is
% 
% ------------------------------------------------------------
% Reduced energy use due to step increase of penalty
% Increased energy use in anticipation of penalty increase
% ------------------------------------------------------------
[r c] = size(FF);
Reduced_EU = 0; 
Increased_EU = 0; 
    for j = 1:c
        if FF(j)<0
            Reduced_EU = Reduced_EU + FF(j);
            %disp(Reduced_EU)
        else
            Increased_EU = Increased_EU + FF(j);
            %disp(Increased_EU)
        end
    end
fprintf('Reduced energy use: %.2f kWh\n', Reduced_EU)
fprintf('Increased energy use : %.2f kWh\n', Increased_EU)

% % --------------------------------------------------
% % Maximum power increase & Maximum power reduction
% % -------------------------------------------------- 

Max_power_red = min(FF,[],'all');
Max_power_inc = max(FF,[],'all');    
 
fprintf('Maximum power reduction : %.2f kWh\n', Max_power_red)
fprintf('Maximum power increase: %.2f kWh\n', Max_power_inc)

% % ------------------------------- 
% % Time until maximum reduction   
% % Duration of demand reduction  
% % ------------------------------- 
% 
% [col_max_red] = find(FF == minValue);
% t2 = col_max_red*Ts*60;
% alpha = t2 - t1; % Time until maximum reduction
% 
%         
% % -----------------------------------------------------
% % Time between maximum power and step change penalty
% % Duration of anticipation
% %------------------------------------------------------
% 
% [col_max_inc] = find(FF == maxValue);
% t3 = col_max_inc*Ts*60;
% gamma = t1 - t3; % Time between maximum power and step change penalty
% 
% %%%%%%%%%%%%%%%%%%%%%%%%
% %%% Methode Reynders %%%
% %%%%%%%%%%%%%%%%%%%%%%%%

