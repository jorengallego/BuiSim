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
       if (k+4)/7 == round((k+4)/7)
           TCF_index = [TCF_index];
       elseif (k+3)/7 == round((k+3)/7)
           TCF_index = [TCF_index];
       else
           TCF_index = [TCF_index, k*day_steps+TCF_start_steps:k*day_steps+TCF_end_steps];
       end
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
        if RefsParam.Price.Belpex
            Belpex2018 = xlsread('BelpexFilter.xlsx');
            references.ElectricityPrice = [];
            for i = 1:8760
                references.ElectricityPrice = [references.ElectricityPrice, Belpex2018(i), Belpex2018(i), Belpex2018(i), Belpex2018(i)];
            end   
            references.ElectricityPrice = references.ElectricityPrice'/1000; % Prices in €/kWh
            references.Price = references.ElectricityPrice*Ts/3600/1000; % [€/W]
        else   
        %        references.Price = 1+sin(0.01*(1:length(WB)))'; % variable price profile
            references.ElectricityPrice = stdprice + factor*heaviside((1:length(WB))-(((day*24)+hour)*4))'; % + 0.8*heaviside((1:length(WB))-(((187*24)+14)*4))' ; % [€/kWh] ((#days*24hours) + statpoint step) * #quarters in 1 hour
        %        references.ElectricityPrice = rectangularPulse(2,3,1:length(WB));
            references.Price = references.ElectricityPrice*Ts/3600/1000; % [€/W]
            %    TODO:  load price profile interface
        end    
    else
       references.ElectricityPrice = 0.25*ones(size(WB)); % [€/kWh]
       references.Price = references.ElectricityPrice*Ts/3600/1000;  % standard fixed price [€/W]
    end
    
 
    %% COP Heat pump
   
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
    
    if RefsParam.HP.ground
        % Ground source heat pump
        c0 = 9.45851; c1 = 0.20159; c2 = -0.17019;
        c3 = 0.00039693; c4 = 0.00092547; c5 = -0.0023720;
        Tground_index = 224;
        Tground = dist.d(:,Tground_index) - 273.15;
        COP = c0 + c1*Tground + c2*Tsupply + c3*(Tground.^2) + c4*(Tsupply.^2) + c5*Tground.*Tsupply;
        references.COP = repmat(COP,1,model.pred.nu)';
    else
        % Air source heat pump
        COP = 0.38*((Tsupply+273.15)./(Tsupply - Te + 10));
        references.COP = repmat(COP,1,model.pred.nu)';
        references.COP(references.COP>6) = 6;
    end
    
    SPF = 13.42*ones(size(Te,1),1);
    references.SPF = repmat(SPF,1,model.pred.nu)';
    
else
    %% comfort boundaries
    % function eval for full year comfort boundaries profiles in K
    %[t_comf, TLow, TUp, TSup, TRefControl] = comfortTemperature(bui_path);

    % % visualisation of comfort constraints and reference
    % plot(t_comf,[TLow, TUp, TLow+(TUp-TLow)/2])
    % legend('TLow','TUp','ref')
    % figure
    % plot(t_comf,TSup)
    % legend('supply watter')
    
%% cofort zone based on standards ISO7730 and EN15251
    % for more details see Damians PhD page 168

    % Te(i)  - ambient temperature at timestep i
    Te_index = 41;    % TODO: generalize ambient temperature indexing for all models in BuiSim
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
    % If you want the building to be kept on a constant temperature
%     WB(:,:) = 22 + 273.15;
%     WA=WB;


    %%  night setbacks
    % thermal comfort start and end hour
    TCF_start1_hour = 7;
    TCF_end1_hour = 9;
    TCF_start2_hour = 17;
    TCF_end2_hour = 22;
    % hours to simsteps
    TCF_start1_steps = TCF_start1_hour*3600/model.plant.Ts;
    TCF_end1_steps = TCF_end1_hour*3600/model.plant.Ts;
    TCF_start2_steps = TCF_start2_hour*3600/model.plant.Ts;
    TCF_end2_steps = TCF_end2_hour*3600/model.plant.Ts;

    % therrmal comfort index of active timesteps
    TCF_index = [];
    for k = 0:days-1
       if (k+4)/7 == round((k+4)/7)
           TCF_index = [TCF_index, k*day_steps+TCF_start1_steps:k*day_steps+TCF_end2_steps];
       elseif (k+3)/7 == round((k+3)/7)
           TCF_index = [TCF_index, k*day_steps+TCF_start1_steps:k*day_steps+TCF_end2_steps];
       else
           TCF_index = [TCF_index, k*day_steps+TCF_start1_steps:k*day_steps+TCF_end1_steps, k*day_steps+TCF_start2_steps:k*day_steps+TCF_end2_steps ];
       end
    end

    % night setback index
    NSB = ones(length(WA),1);
    NSB(TCF_index) = 0;
    % find([1:steps], TCF_index)


    % temperature comfort zone
    references.R = (WB+(WA-WB)/2)*ones(1,model.pred.ny); 
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
        if RefsParam.Price.Belpex
            Belpex2018 = xlsread('BelpexFilter.xlsx');
            references.ElectricityPrice = [];
            for i = 1:8760
                references.ElectricityPrice = [references.ElectricityPrice, Belpex2018(i), Belpex2018(i), Belpex2018(i), Belpex2018(i)];
            end   
            references.ElectricityPrice = references.ElectricityPrice'/1000; % Prices in €/kWh
            references.Price = references.ElectricityPrice*Ts/3600/1000; % [€/W]
        else   
        %        references.Price = 1+sin(0.01*(1:length(WB)))'; % variable price profile
            references.ElectricityPrice = stdprice + factor*heaviside((1:length(WB))-(((day*24)+hour)*4))';% - factor*heaviside((1:length(WB))-(((day*24)+hour+5)*4))' ; % [€/kWh] ((#days*24hours) + statpoint step) * #quarters in 1 hour
        %        references.ElectricityPrice = rectangularPulse(2,3,1:length(WB));
            references.Price = references.ElectricityPrice*Ts/3600/1000; % [€/W]
            %    TODO:  load price profile interface
        end    
    else
       references.ElectricityPrice = 0.25*ones(size(WB)); % [€/kWh]
       references.Price = references.ElectricityPrice*Ts/3600/1000;  % standard fixed price [€/W]
    end

    
    %% COP Heat pump
    
    Tsupply = zeros(size(Te,1),1); % Supply temperature based on heating curve PhD Damien Picard
    for i = 1:size(Tsupply,1)
        if Te(i) <= -7
            Tsupply(i) = 75;
        elseif -7 < Te(i) <= 15
            Tsupply(i) = 68.63636364 - (20/22)*Te(i);
        else
            Tsupply(i) = 55;
        end
    end
    
    if  strcmp(model.buildingType,'Old')
        % High temperature HP 16kW
        c0 = 1.81737; c1 = 0.091224; c2 = 0.043883; 
        c3 = 0.0014346; c4 = -0.00053603; c5 = -0.00073929;
    elseif strcmp(model.buildingType,'Reno')
        % High temperature HP 14kW
        c0 = 1.85031; c1 = 0.10354; c2 = 0.046302; 
        c3 = 0.0016897; c4 = -0.00057315; c5 = -0.00083778;
    elseif strcmp(model.buildingType,'RenoLight')
        if RefsParam.HP.hightemperature
            % High temperature HP 11kW
            c0 = 2.53957; c1 = 0.11909; c2 = 0.028578;
            c3 = 0.0015845; c4 = -0.00045313; c5 = -0.00095888; 
        else    
            % Low temperature HP 6kW
            c0 = 6.57644; c1 = 0.16312; c2 = -0.10941;
            c3 = 0.00046445; c4 = 0.00051668; c5 = -0.0015987;
            for i = 1:size(Tsupply,1)
                if Te(i) <= -7
                    Tsupply(i) = 55;
                elseif -7 < Te(i) <= 15
                    Tsupply(i) = 45.454545 - (30/22)*Te(i);
                else
                    Tsupply(i) = 25;
                end
            end
        end   
        
    end
%     if RefsParam.HP.ground
%         % Ground source heat pump
%         Tground_index = 44;
%         Tground = dist.d(:,Tground_index) - 273.15;
%         COP = c0 + c1*Tground + c2*Tsupply + c3*(Tground^2) + c4*(Tsupply^2) + c5*Tground*Tsupply;
%     else
        % Air source heat pump
    COP = c0 + c1*Te + c2*Tsupply + c3*(Te.^2) + c4*(Tsupply.^2) + c5*Te.*Tsupply;

%     end
    references.COP = repmat(COP,1,model.pred.nu)';
    
    SPF = 15*ones(size(Te,1),1);
    references.SPF = repmat(SPF,1,model.pred.nu)';

    
end

fprintf('*** Done.\n')
end