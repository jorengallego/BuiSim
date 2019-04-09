Belpex2018 = xlsread('BelpexFilter.xlsx');
BelpexPrice = [];
for i = 1:8760
    BelpexPrice = [BelpexPrice, Belpex2018(i), Belpex2018(i), Belpex2018(i), Belpex2018(i)];
end   
BelpexPrice = BelpexPrice'/1000; % Prices in €/kWh
    