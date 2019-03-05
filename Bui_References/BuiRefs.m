function references = BuiRefs(model, RefsParam)

% TODO:  FIX THIS 
% TODO: Variable TZoneSetLow not found in Data/preComputedInfrax.mat

if nargin < 2
   RefsParam.Price.variable = 0; %1 =  variable price profile, 0 = fixed to 1
end



if nargin == 0
   model.buildingType = 'Reno';  
   model.pred.ny = 6;
%    addpath('../common_files_all_sims/Data')
end

% TODO: generalize this temporary fix
if strcmp(model.buildingType,'Infrax') 
   model.pred.ny = 19; 
   model.buildingType = 'Reno';
elseif strcmp(model.buildingType,'HollandschHuys') 
   model.pred.ny = 12; 
   model.buildingType = 'Reno';
end
  

bui_path = ['../buildings/', model.buildingType, '/disturbances/'];
fprintf('*** Load references ... \n')

%% comfort boundaries
% function eval for full year comfort boundaries profiles in K
[t_comf, TLow, TUp, TSup, TRefControl] = comfortTemperature(bui_path);

%Zelfgeschreven functie om comfort boundaries makkelijker te kunnen
%aanpassen:
for k = 0:366
    for i = 1:33
        TLow((k*96)+i)= 290.15;
        TUp((k*96)+i)= 299.15;
    end
    for i = 34:93
        TLow((k*96)+i) = 293.15;
        TUp((k*96)+i) = 296.15;
    end
    for i = 94:96
        TLow((k*96)+i) = 290.15;
        TUp((k*96)+i) = 299.15;
    end
end
TLow = TLow(1:35224);
TUp = TUp(1:35224);

% % visualisation of comfort constraints and reference
% plot(t_comf,[TLow, TUp, TLow+(TUp-TLow)/2])
% legend('TLow','TUp','ref')
% figure
% plot(t_comf,TSup)
% legend('supply watter')

references.R = (TLow+(TUp-TLow)/2)*ones(1,model.pred.ny);  % setpoint in K
% references.ref = (TRefControl+2.5)*ones(1,model.pred.ny);  % setpoint in K
references.wb = TLow*ones(1,model.pred.ny); %  below threshold
references.wa = TUp*ones(1,model.pred.ny); %  above threshold
% PMV bonds
references.PMVub(TLow == min(TLow)) = 1;
references.PMVub(TLow > min(TLow)) = 0.5;
references.PMVlb = -references.PMVub;

references.TSup = TSup;   % supply water temperature from HC

%% variable price profiles
Ts = model.plant.Ts;
if RefsParam.Price.variable
%    references.Price = 1+sin(0.01*(1:length(TLow)))'; % variable price profile
   references.EnergyPrice = 0.25 + 0.8*heaviside((1:length(TLow))-(((187*24)+12)*4))'; % + 0.8*heaviside((1:length(TLow))-(((187*24)+14)*4))' ; % [�/kWh] ((#days*24hours) + statpoint step) * #quarters in 1 hour
   %references.EnergyPrice = rectangularPulse(2,3,1:length(TLow));
   references.Price = references.EnergyPrice*Ts/3600/1000; % [�/kW]
%    TODO:  load price profile interface
else
   references.Price = ones(size(TLow));  % standard fixed price 
end

%% COP Heat pump

% Constant COP

references.COP = 3;

% Time varying COP

TSource = 273.5 + 15;
references.COP = zeros(size(TSup,1),model.pred.ny);
for i = 1:model.pred.ny
    references.COP(:,i) = 0.45*(TSup./(TSup - TSource + 10));
end

fprintf('*** Done.\n')
end