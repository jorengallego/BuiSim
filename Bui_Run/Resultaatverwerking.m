x = [1 2 5 10 50 100];
y1 = [ 0 1.36 4.13 7.68 15.74 18.61 ];
y2 = [ 0 0.72 2.08 3.65 6.98 8.21 ];
y3 = [ 0 661.49 1907.45 3651.31 5684.99 6895.34 ];
y4 = [ 0 1357.45 2392.88 2772.91 3437.28 3575.38 ];

figure
plot(x, y1, 'linewidth', 2)
xlabel('step')
ylabel('Energy flexibility [kWh]')
title('Energy flexibility with Increasing price step (01/01 – 02/01)')
hold on 
plot (x,y2, 'linewidth', 2)
legend({'Increased energy use','reduced energy use'}, 'Location', 'northwest')

figure
plot(x,y3, 'linewidth', 2)
xlabel('step')
ylabel('Power flexibility [W]')
title('Power flexibility with Increasing price step (01/01 – 02/01)')
hold on 
plot (x,y4,'linewidth', 2)
legend({'Maximum power increase','Maximum power reduction'}, 'Location', 'northwest')
