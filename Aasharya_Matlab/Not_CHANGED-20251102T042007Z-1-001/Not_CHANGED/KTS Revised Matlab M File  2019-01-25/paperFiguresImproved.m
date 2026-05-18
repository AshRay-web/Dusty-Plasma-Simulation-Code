% GLOBAL STYLE
LW = 1.8;
FS = 12;
FS_title = 14;

% Color palette (consistent across all plots)
cBlue   = [0.10 0.45 0.85];
cRed    = [0.85 0.10 0.10];
cOrange = [0.90 0.50 0.10];
cGreen  = [0.20 0.60 0.30];



% This is a script run it to create all the figures that are in the paper.

data1 = load('DustyPlasma2.mat');    % Plasma output
data2 = load('DustyPlasma1.mat');    % The dusty-plasma output %temp 1.0 eV

data5= load("DustyPlasma5.mat"); % Temp 0.5 ion 
data6= load("DustyPlasma6.mat"); % Temp 2 ion


x = data1.X / data1.DLe;


figure(1);
clf;

% --------- Make figure wide (important) ---------
set(gcf, 'Units', 'pixels', 'Position', [100 100 1200 1200]);  % wide figure
set(gcf, 'PaperPositionMode', 'auto');

% ---------------- Main plot ----------------
plot(x, data1.Ef, 'o', ...
    'LineStyle', 'none', ...
    'LineWidth', LW, ...
    'MarkerSize', 6, ...
    'MarkerFaceColor', 'none');
hold on;

plot(x, data2.Ef, '*', ...
    'LineStyle', 'none', ...
    'LineWidth', LW, ...
    'MarkerSize', 7);
hold off;

xlabel('X/DLe', 'FontSize', FS);
ylabel('Electric Field (E_f)', 'FontSize', FS);
title('Electric Field Comparison', 'FontSize', FS_title);

legend({'Ordinary Plasma', 'Dusty Plasma'}, 'Location', 'best');
grid on;
box on;
set(gca, 'FontSize', FS, 'LineWidth', 1.2);


% ---------------- Save ----------------
saveas(gcf, 'Electric_Field_Comparision_Dusty_Ordinary_Plasma.png');




%% FIGURE 2 – Electron density
figure(2); clf;
plot(x, data1.ne, 'Color', cRed, 'LineWidth', LW); hold on;
plot(x, data2.ne, 'Color', cBlue, 'LineWidth', LW);

xlabel("X/DLe", 'FontSize', FS);
ylabel("Electron Density (n_e)", 'FontSize', FS);
title("Electron Density Comparison", 'FontSize', FS_title);
legend({'Ordinary Plasma', 'Dusty Plasma'}, 'Location','best');
grid on; box on;

saveas(gcf, 'Electron_Density_Comparison.png');



%% FIGURE 3 – Plasma potential
figure(3); clf;
plot(x, data1.phi, 'Color', cRed, 'LineWidth', LW); hold on;
plot(x, data2.phi, 'Color', cBlue, 'LineWidth', LW);

xlabel("X/DLe", 'FontSize', FS);
ylabel("\phi (potential)", 'FontSize', FS);
title("Potential Evolution Comparison", 'FontSize', FS_title);
legend({'Ordinary Plasma','Dusty Plasma'}, 'Location','best');
grid on; box on;

saveas(gcf, 'Potential_Evolution_Comparison.png');



%% FIGURE 7 – Dust charge (temperature comparison)
figure(7); clf;

plot(x, data5.rho , 'Color', cRed, 'LineWidth', LW); hold on;
plot(x, data2.rho , 'Color', cBlue, 'LineWidth', LW);
plot(x, data6.rho , 'Color', cOrange, ...
     'LineWidth', LW, 'MarkerFaceColor', cOrange, 'MarkerSize', 5);

xlabel("X/DLe", 'FontSize', FS);
ylabel("Charge Density", 'FontSize', FS);
title("Charge density near the wall (0 ≤ X/DLe ≤ 4)", 'FontSize', FS_title);
legend({'0.5 eV', '1 eV', '2 eV'}, 'Location','best');

xlim([0 4]);       
grid on; box on;

saveas(gcf, 'ChargeDensity_vs_Ion_Tempreature_wall.png');


%% here try to caputre the space charge density near the wall for different ion tempearture dustty plasma

%% FIGURE 5 – Dusty vs Ordinary rho
figure(5); clf;
plot(x, data1.rho, 'Color', cRed, 'LineWidth', LW); hold on;
plot(x, data2.rho, 'Color', cBlue, 'LineWidth', LW);

xlabel("X/DLe", 'FontSize', FS);
ylabel("Space Charge Density (\rho)", 'FontSize', FS);
title("Charge density near the wall (0 ≤ X/DLe ≤ 4)", 'FontSize', FS_title);
legend({'Ordinary Plasma','Dusty Plasma'}, 'Location','best');
xlim([0 4]); 
grid on; box on;

saveas(gcf, 'Space_Charge_Density_Comparison.png');


%% FIGURE 6 – Bohm sheath velocity (v_x_i)
figure(6); clf;

y1 = 1:size(data1.vxi,1);
y2 = 1:size(data2.vxi,1);

subplot(1,2,1);
contourf(x, y1, data1.vxi, 30, 'LineColor','none');
colorbar; colormap('turbo');
xlabel('X/DLe'); ylabel('Index');
title('v_{xi} – Ordinary Plasma');
axis tight;

subplot(1,2,2);
contourf(x, y2, data2.vxi, 30, 'LineColor','none');
colorbar; colormap('turbo');
xlabel('X/DLe'); ylabel('Index');
title('v_{xi} – Dusty Plasma');
axis tight;

cmin = min([data1.vxi(:); data2.vxi(:)]);
cmax = max([data1.vxi(:); data2.vxi(:)]);
subplot(1,2,1); caxis([cmin cmax]);
subplot(1,2,2); caxis([cmin cmax]);

sgtitle("Comparison of v_{xi} Distributions", 'FontSize', FS_title);

saveas(gcf, 'Comparison_vxi_Distributions.png');
