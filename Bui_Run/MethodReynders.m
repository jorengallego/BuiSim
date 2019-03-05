
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Methode Glenn Reynders %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load full reference profile of one year
load('OutputOldElectricity.mat') 
Q_REF_Full = sum(RefElectricity); % sum(A) telt alle kolomelementen van matrix A op resulterend in 1 rijmatrix


% Load step response for desired amount of days
BuiInit

% Simulation parameters
SimStart = outdata.SimParam.SimStart;
SimStop = outdata.SimParam.SimStop;
Nsim = outdata.SimParam.run.Nsim; 
Ts = model.plant.Ts;
Time = (1:Nsim)*outdata.model.plant.Ts/3600/24;

Q_REF = Q_REF_Full(SimStart:SimStop);
Q_ADR = sum(outdata.data.U); 
Q_verschil_ADR = Q_ADR - Q_REF; % Intervallen nog aanpassen (deze berekening moet maar gaan voor l_ADR
Q_verschil_TOT = Q_ADR - Q_REF;

%% Available storage capacity

C_ADR = 0;  
for j = 1:size(Q_verschil_ADR,2)
    C_ADR = C_ADR + Q_verschil_ADR(j)*Ts/(3.6*10^6);
end

%% Efficiency of the storage process

C_ADR_TOT = 0;  
for j = 1:size(Q_verschil_TOT,2)
    C_ADR_TOT = C_ADR_TOT + Q_verschil_TOT(j)*Ts/(3.6*10^6);
end
Teller = C_ADR_TOT;
Noemer = C_ADR;
Efficiency = 1 - Teller/Noemer;