load ("OML_NOB.mat")
plot(X/DLe, phi, 'LineWidth', 1.5)
hold on

load ("OML_e12nd_1e-6.mat")
plot(X/DLe, phi, 'LineWidth', 1.5)



legend('0 B', '0.20 B','Location', 'NorthEast')
xlabel('X / D_{Le}')
ylabel('Potential PHI')
title("Dusty Plasma Magnetic Vs Non Magnetic Field")
grid on
saveas(gcf, 'Magnetic vs Non Magnetic Dusy Plasma.png')

