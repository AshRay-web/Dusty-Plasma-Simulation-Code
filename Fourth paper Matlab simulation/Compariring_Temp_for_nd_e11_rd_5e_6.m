load ("OML_e11nd5e-6.mat")
plot(X/DLe, phi, 'LineWidth', 1.5)
hold on

load ("OML_e11nd_rd_5e-6_temp13ev")
plot(X/DLe, phi, 'LineWidth', 1.5)

load ("OML_e11nd_rd_5e-6_temp14ev")
plot(X/DLe, phi, 'LineWidth', 1.5)

%load ("OML_e11nd5e-6temp14ev.mat")
%plot(X/DLe, ni, 'LineWidth', 1.5)


legend('Tempe = 11ev', 'Tempe=13ev',"Tempe=14ev", 'Location', 'SouthEast')
xlabel('X / D_{Le}')
ylabel('Potential')
title("Electron Temperature vs Potential (Dusty Plasma)")
grid on
saveas(gcf, 'Temperature vs Potential (DUSTY PLASMA).png')


