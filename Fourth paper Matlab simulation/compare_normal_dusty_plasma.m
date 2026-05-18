load ("new_OML_output.mat")
plot(X/DLe, ni, 'LineWidth', 1.5)
hold on

%load ("new_OML_output.mat")
%plot(X/DLe, phi, 'LineWidth', 1.5)

load ("new_OML_output_temp=1.2.mat")
plot(X/DLe, ni, 'LineWidth', 1.5)


legend('temp=1.33', 'temp=1.2', 'Location', 'SouthEast')
xlabel('X / D_{Le}')
ylabel('Phi Wall Potential')
title("Dust density vs Phi(OML\_NO\_ITER)")
grid on

% ==== make inset plot ====
axes('Position',[0.55 0.55 0.3 0.3])   % [x y width height] normalized coords
box on
hold on

% Replot same data in zoom range
load ("new_OML_output.mat")
plot(X/DLe, phi, 'LineWidth', 1.5)

load ("new_OML_output_temp=1.2.mat")
plot(X/DLe, phi, 'LineWidth', 1.5)

%load ("OML_rd_5e-6_ndust_1e13_B_200e-3_NOITER.mat")
%plot(X/DLe, phi, 'LineWidth', 1.5)

xlim([1.8 3])   % zoomed region in x
ylim([-0.5 -0.25]) % adjust so you see the "interesting" top

%saveas(gcf,"Dustdensity_Vs_Potential(NO_ITER_OML).png")
