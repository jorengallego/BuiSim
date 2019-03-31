function references = BeRefs(model, dist, RefsParam)

% TODO:  FIX THIS 
% TODO: Variable TZoneSetLow not found in Data/preComputedInfrax.mat

if nargin < 3
   RefsParam.Price.variable = 0; %1 =  variable price profile, 0 = fixed to 1
end

if nargin == 0
   model.buildingType = 'Reno';  
   model.pred.ny = 6;
%    addpath('../common_files_all_sims/Data')
end

% TODO: generalize this temporary fix - use Te for infrax
if strcmp(model.buildingType,'Infrax') 
   model.pred.ny = 19; 
   model.buildingType = 'Reno';
end
  
bui_path = ['C:\Users\Joren\Documents\BeSim\buildings\', model.buildingType, '/disturbances/'];
% bui_path = ['../buildings/', model.buildingType, '/disturbances/'];
fprintf('*** Load references ... \n')



if  strcmp(model.buildingType,'HollandschHuys')

    %% cofort zone based on standards ISO7730 and EN15251
    % for more details see Damians PhD page 168

    % Te(i)  - ambient temperature at timestep i
    Te_index = 222;    % TODO: generalize ambient temperature indexing for all models in BuiSim
    Te = dist.d(:,Te_index) - 273.15;
    day_steps = 86400/model.pred.Ts;     % number of steps per day,  24h = 86400 sec
    steps = length(Te);                  % number of datapoints in dataset
    days = floor(steps/day_steps);  %  number of days in dataset

    % Tem(k) - one day average ambient temperature at day k
    for k = 0:days-1
       Tem(k+1) =  mean(Te(1+k*day_steps:(1+k)*day_steps));
    end
    % Tem = movmean(Te,day_steps);

    % TRM =  movmean(Tem,7);
    for k = 1:days
        % Trm(k) - moving mean ambient temperature at day k
        if k == 1
            Trm(k) = Tem(k);
        elseif k == 2
            Trm(k) = (Tem(k) + 0.8*Tem(k-1))/1.8;
        elseif k == 3
            Trm(3) = (Tem(k) + 0.8*Tem(k-1) + 0.6*Tem(k-2))/2.4;    
        elseif k == 4
            Trm(4) = (Tem(k) + 0.8*Tem(k-1) + 0.6*Tem(k-2) + 0.5*Tem(k-3))/2.9;       
        elseif k == 5
            Trm(5) = (Tem(k) + 0.8*Tem(k-1) + 0.6*Tem(k-2) + 0.5*Tem(k-3) + 0.4*Tem(k-4))/3.3;   
        elseif k == 6
            Trm(k) = (Tem(k) + 0.8*Tem(k-1) + 0.6*Tem(k-2) + 0.5*Tem(k-3) + 0.4*Tem(k-4) + 0.3*Tem(k-5))/3.6; 
        else
            Trm(k) = (Tem(k) + 0.8*Tem(k-1) + 0.6*Tem(k-2) + 0.5*Tem(k-3) + 0.4*Tem(k-4) + 0.3*Tem(k-5) + 0.2*Tem(k-6))/3.8;  
        end

        % upper and lower comfort bounds at individual days
        if Trm(k) < 10   % heating season comfort bounds
            TUp(k) = 24;
            TLow(k) = 20;
        elseif Trm(k) > 15  % cooling season comfort bounds
            TUp(k) = 26;
            TLow(k) = 23;
        else  % transient season comfort bounds
            TUp(k) = 24 + 2*(Trm(k)-10)/5;
            TLow(k) = 20+ 3*(Trm(k)-10)/5; 
        end
    end

    % repeat elements of vector to fit the sampling
    WB = repelem(TLow,day_steps)' + 273.15;
    WA = repelem(TUp,day_steps)' + 273.15;


    %%  night setbacks
    % thermal comfort start and end hour
    TCF_start_hour = 7;
    TCF_end_hour = 20;
    % hours to simsteps
    TCF_start_steps = TCF_start_hour*3600/model.plant.Ts;
    TCF_end_steps = TCF_end_hour*3600/model.plant.Ts;

    % therrmal comfort index of active timesteps
    TCF_index = [];
    for k = 0:days-1
       TCF_index = [TCF_index, k*day_steps+TCF_start_steps:k*day_steps+TCF_end_steps];
    end

    % night setback index
    NSB = ones(length(WA),1);
    NSB(TCF_index) = 0;
    % find([1:steps], TCF_index)


    % temperature comfort zone
    references.wb = WB*ones(1,model.pred.ny) - NSB*2;%  below threshold
    references.wa = WA*ones(1,model.pred.ny) + NSB*2; %  above threshold

    % PMV comfort zone
    references.PMVub(WB == min(WB)) = 1;
    references.PMVub(WA > min(WA)) = 0.5;
    references.PMVlb = -references.PMVub;


    %% variable price profiles
    stdprice = RefsParam.Price.stdprice;
    factor = RefsParam.Price.factor;
    day = RefsParam.Price.day-1;
    hour = RefsParam.Price.hour;
    Ts = model.plant.Ts;
    if RefsParam.Price.variable
%        references.Price = 1+sin(0.01*(1:length(WB)))'; % variable price profile
       references.EnergyPrice = stdprice + factor*heaviside((1:length(WB))-(((day*24)+hour)*4))'; % + 0.8*heaviside((1:length(WB))-(((187*24)+14)*4))' ; % [€/kWh] ((#days*24hours) + statpoint step) * #quarters in 1 hour
%        references.EnergyPrice = rectangularPulse(2,3,1:length(WB));
       references.Price = references.EnergyPrice*Ts/3600/1000; % [€/W]
    %    TODO:  load price profile interface
    else
       references.Price = 0.25*Ts/3600/1000*ones(size(WB));  % standard fixed price 
       % references.EnergyPrice = ones(size(WB));
    end
    
    %% COP Heat pump

    % Constant COP

%     references.COP = 3;

    % Time varying COP

    
    Te3d = zeros(size(Te,1),1); % The previous three days average ambient temperature needed for calculation of supply temperature with heating curve
    for i = 1:size(Te,1)
        if i < 3*day_steps+1
            Te3d(i) = (sum(Te(1:i)) + sum(Te(end - 3*day_steps + i:end)))/(3*day_steps);
        else
            Te3d(i) = sum(Te(i-3*day_steps:i))/(3*day_steps);
        end
    end
    
    Tsupply = zeros(size(Te,1),1); % Supply temperature based on heating curve PhD Damien Picard
    for i = 1:size(Tsupply,1)
        if Te3d(i) <= -8
            Tsupply(i) = 29;
        elseif -8 < Te3d(i) <= 15
            Tsupply(i) = 29 - (7/23)*(Te3d(i) + 8);
        elseif 15 < Te3d(i) <= 18
            Tsupply(i) = 22;
        elseif 18 < Te3d(i) <= 30
            Tsupply(i) = 22 - (5/12)*(Te3d(i) - 18);
        else
            Tsupply(i) = 17;
        end
    end 
        
    references.COP = zeros(size(Te,1),model.pred.nu);
    for i = 1:model.pred.nu
        references.COP(:,i) = 0.45*((Tsupply+273.15)./(Tsupply - Te + 10)); % Correlation formula for COP
    end
    references.COP = references.COP';

else
    %% comfort boundaries
    % function eval for full year comfort boundaries profiles in K
    [t_comf, TLow, TUp, TSup, TRefControl] = comfortTemperature(bui_path);

    % % visualisation of comfort constraints and reference
    % plot(t_comf,[TLow, TUp, TLow+(TUp-TLow)/2])
    % legend('TLow','TUp','ref')
    % figure
    % plot(t_comf,TSup)
    % legend('supply watter')
    
%     for k = 0:366
%         for i = 1:33
%             TLow((k*96)+i)= 290.15;
%             TUp((k*96)+i)= 299.15;
%         end
%         for i = 34:93
%             TLow((k*96)+i) = 293.15;
%             TUp((k*96)+i) = 296.15;
%         end
%         for i = 94:96
%             TLow((k*96)+i) = 290.15;
%             TUp((k*96)+i) = 299.15;
%         end
%     end
%     TLow = TLow(1:35224,:);
%     TUp = TUp(1:35224,:);

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
    stdprice = RefsParam.Price.stdprice;
    factor = RefsParam.Price.factor;
    day = RefsParam.Price.day-1;
    hour = RefsParam.Price.hour;
    Ts = model.plant.Ts;
    if RefsParam.Price.variable
    %    references.Price = 1+sin(0.01*(1:length(TLow)))'; % variable price profile
       references.EnergyPrice = stdprice + factor*heaviside((1:length(TLow))-(((day*24)+hour)*4))'; % + 0.8*heaviside((1:length(TLow))-(((187*24)+14)*4))' ; % [€/kWh] ((#days*24hours) + statpoint step) * #quarters in 1 hour
       %references.EnergyPrice = rectangularPulse(2,3,1:length(TLow));
       references.Price = references.EnergyPrice*Ts/3600/1000; % [€/W]
    %    TODO:  load price profile interface
    else
       references.Price = 0.25*Ts/3600/1000*ones(size(TLow));  % standard fixed price 
       % references.EnergyPrice = ones(size(TLow));
    end

    
    %% COP Heat pump

    % Constant COP

%     references.COP = 3;

    % Time varying COP
    Te_index = 41;
    Te = dist.d(:,Te_index);        % Ambient temperature
    Tsupply = 273.15 + 20;                % Supply temperature
%     TSource = 273.5 + 15;
    references.COP = zeros(size(Te,1),model.pred.nu);
    for i = 1:model.pred.nu
        references.COP(:,i) = 0.45*(Tsupply./(Tsupply - Te + 10));
    end
    references.COP = references.COP';
end

fprintf('*** Done.\n')
end