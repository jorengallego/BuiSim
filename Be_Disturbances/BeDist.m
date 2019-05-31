function dusturb = BeDist(model, DistParam)

if nargin == 0
   model.buildingType = 'Infrax';  
end
if nargin < 2
   DistParam.reload = 0;
end

path = ['C:\Users\Joren\Documents\BeSim\buildings\', model.buildingType];
% path = ['../buildings/', model.buildingType];
disturbanceType = ''; % can be '_lin' if used for linearization validation

if DistParam.reload
        fprintf('*** Load disturbances ... \n')
		% Disturbances
        if strcmp(model.buildingType,'Reno') || strcmp(model.buildingType,'RenoLight') || strcmp(model.buildingType,'Old') 
            [t, v, x0] = disturbances_old(path, 0, 0);
            save([path '/preComputed_matlab/dis' disturbanceType '.mat'], 't', 'v', 'x0');
        else
            [t,  v, inputIndex, dictCtlInputs, dicValVar, dicOutputNameIndex, x0]= disturbances(path,disturbanceType,0, 0);
            save([path '/preComputed_matlab/dis' disturbanceType '.mat'], 't', 'v', 'inputIndex', 'dictCtlInputs', 'dicValVar','dicOutputNameIndex', 'x0');
        end	
		fprintf('*** Done.\n')
else
        fprintf('*** Load disturbances...\n')
		load([path '/preComputed_matlab/dis' disturbanceType '.mat']);
        fprintf('*** Done.\n')
%         TODO:  create mod.mat file also for 6-zone building - connect
%         models in one file
end

        
%% disturbances
if  strcmp(model.buildingType,'HollandschHuys')
%     including fixed ventilation temperature as disturbance for HH
    VenTsup_temp = 20+273.15;
    VenTsup = repmat(VenTsup_temp,[size(v,1),12]);
    dusturb.t = t;
    dusturb.d = [v, VenTsup];
%     load('InternalGainsHH.mat')
    dusturb.d(:,278:289) = 1.3804*dusturb.d(:,278:289);
    
elseif strcmp(model.buildingType,'Reno') || strcmp(model.buildingType,'RenoLight') || strcmp(model.buildingType,'Old')
    HeatGains_temp = [];     %Set internal gains of residential building according to Bram vd Heijde
    for i = 1:96 
        if ismember(i,1:27) || ismember(i,36:67) || ismember(i,88:96)
            HeatGains_temp(i) = 200;
        else
            HeatGains_temp(i) = 835;
        end
    end
    HeatGains = repmat(HeatGains_temp',365,6);
    HeatGains = [HeatGains; HeatGains(1:(size(v,1)-size(HeatGains,1)),:)]; %Weekends nog integreren (zoals bij comfortgrenzen)
    dusturb.t = t;
    dusturb.d = [v, HeatGains];
else
    dusturb.t = t;
    dusturb.d = v;
end


% max and min disturbances
dusturb.dmin = min(dusturb.d, [], 1);
dusturb.dmax = max(dusturb.d, [], 1);

%% TODO: automatic evaluation of statistical properties of datasets
% probability distributions, mean, min, max, etc.


end