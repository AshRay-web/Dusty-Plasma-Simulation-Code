load ("OML_e11nd5e-6.mat")
plot(X/DLe, ni, 'LineWidth', 1.5)
hold on

load ("OML_e11nd_rd_3e-6.mat")
plot(X/DLe, ni, 'LineWidth', 1.5)

load ("OML_e11nd_rd_1e-6.mat")
plot(X/DLe, ni, 'LineWidth', 1.5)

legend('rd= 5e-6', 'rd= 3e-6',"rd = 1e-6", 'Location', 'East')
xlabel('X / D_{Le}')
ylabel('Ion Density')
title("Radius of dust particle vs Ion Density")

grid on

saveas(gcf, 'Variation of Radius of Dust particle(Ion_density).png')
