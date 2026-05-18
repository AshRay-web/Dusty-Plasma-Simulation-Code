load ("OML_e11nd5e-6.mat")
plot(X/DLe, phi, 'LineWidth', 1.5)
hold on

load ("OML_e10nd5e-6.mat")
plot(X/DLe, phi, 'LineWidth', 1.5)

load ("OML_e9nd.mat")
plot(X/DLe, phi, 'LineWidth', 1.5)

legend('n_dust = e11', 'n_dust = e10',"ndust = e9", 'Location', "East")
xlabel('X / D_{Le}')
ylabel('Potential')
grid on
title("Variation of Density of Dust in Plasma vs Potential(rd=5e-6)")
saveas(gcf, 'Variation of Dust Density for Plasma_Potential.png')
