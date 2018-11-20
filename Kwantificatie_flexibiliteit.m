
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Methode Flexibility Function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Eerst twee simulaties uitvoeren met BuiSim
% 1) Run MPC with constant cost signal = minimum heat demand (zoals de
%    toolbox standaard ingesteld staat door Jan) = Q_REF
% 2) Run MPC  with step in price signal (nog bekijken met
%    Jan/Bram/Anke/Alessia hoe dit praktisch kan) = Q_STEP

BuiInit

Ts = model.plant.Ts;

Q_REF = sum(outdata.data.U.ref); % sum(A) telt alle kolomelementen van matrix A op resulterend in 1 rijmatrix
Q_STEP = sum(outdata.data.U.step);

FF = Q_STEP - Q_REF;

figure
plot(Ts, FF);
title('Flexibility Function');
ylabel('Heat demand');
xlabel('time');

start_time_penalty = t1; % Tijdslot wanneer penalty signaal verhoogd is

% ------------------------------------------------------------
% Reduced energy use due to step increase of penalty
% Increased energy use in anticipation of penalty increase
% ------------------------------------------------------------

Reduced_EU = 0; 
Increased_EU = 0; 
    for j = 1:size(FF)
        if FF(j)<0
            Reduced_EU = Reduced_EU + FF(j);
        else
            Increased_EU = Increased_EU + FF(j);
        end
    end

% --------------------------------------------------
% Maximum power increase & Maximum power reduction
% -------------------------------------------------- 

Max_power_inc = max(FF,[],'all');    
Max_power_red = min(FF,[],'all'); 

% ------------------------------- 
% Time until maximum reduction   
% Duration of demand reduction  
% ------------------------------- 

[col_max_red] = find(FF == minValue);
t2 = col_max_red*Ts*60;
alpha = t2 - t1; % Time until maximum reduction

        
% -----------------------------------------------------
% Time between maximum power and step change penalty
% Duration of anticipation
%------------------------------------------------------

[col_max_inc] = find(FF == maxValue);
t3 = col_max_inc*Ts*60;
gamma = t1 - t3; % Time between maximum power and step change penalty

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Methode Reynders %%%
%%%%%%%%%%%%%%%%%%%%%%%%

