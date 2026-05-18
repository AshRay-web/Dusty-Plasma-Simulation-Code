% This is a script run it to create all the figures that are in the paper.

data1 = load('DustyPlasma2.mat');    % Plasma output
data2 = load('DustyPlasma1.mat');    % The dusty-plasma output %temp 1.0 eV
data3 = load("DustyPlasma3.mat");    % The argon plasma 
data5= load("DustyPlasma5.mat");     % Temp 0.5 ion 
data6= load("DustyPlasma6.mat");     % Temp 2 ion


x = data1.X / data1.DLe;
%Compairing the change in ion density throughout
figure(1);clf;
plot(x,data1.ni,"r-",'lineWidth',1.5); 

hold on;
plot(x,data2.ni,"b-","lineWidth",1.5) ;
xlabel("X/DLe");ylabel("Ion Density (ni) ");
legend('Ordinary Plasma',"Dusty Plasma",'Location','west');
title("Ion density comparison (ni) ");


%Compairing the change in electron density throughout
figure(2);clf;
plot(x,data1.ne,"r-",'lineWidth',1.5); hold on;
plot(x,data2.ne,"b-","lineWidth",1.5) ;
xlabel("X/DLe");ylabel("Electron Density (ne) ");
legend('Ordinary Plasma',"Dusty Plasma",'Location','west');
title("Electron density comparison (ne) ");

%Compairing the change in the potential
figure(3);clf;
plot(x,data1.phi,"r-","lineWidth",1.5);hold on;
plot(x,data2.phi,"b-","lineWidth",1.5);
xlabel("X/DLe"); ylabel("Phi");
legend('Ordinary Plasma',"Dusty Plasma",'Location','west');
title("Comparision of potential evolution") ;


% Charge accumulate in the dust particle
figure(4);clf;
data3 = load("DustyPlasma3.mat")
plot(x,data3.qdn,"r-","lineWidth",1.5);hold on;
plot(x,data2.qdn,"b-","lineWidth",1.5);
xlabel("X/DLe"); ylabel("Normalized Dust Charge (Qd)");
legend('Argon Plasma',"Hydrogen Plasma",'Location','west');
title("Normalized Dust Charge");

% Charge accumulate in the dust particle
figure(7); clf;

% Use clearer colors
c1 = [0.10 0.45 0.85];   % blue-ish
c2 = [0.85 0.10 0.10];   % red-ish
c3 = [0.90 0.50 0.10];   % orange

plot(x, data5.qdn, 'Color', c2, 'LineWidth', 1.8); hold on;
plot(x, data2.qdn, 'Color', c1, 'LineWidth', 1.8); hold on;
plot(x, data6.qdn, 'Color', c3, 'LineWidth', 1.8, 'MarkerSize', 5, ...
     'MarkerFaceColor', c3);

xlabel("X/DLe", 'FontSize', 12);
ylabel("Normalized Dust Charge (Qd)", 'FontSize', 12);

legend({'0.5 eV', '1 eV', '2 eV'}, 'Location', 'east');

title("Normalized Dust Charge", 'FontSize', 14);

grid on;
box on;

saveas(gcf, 'Charge_Ion_charge.png');


%Space Charge density Comparision
figure(5);clf;
plot(x,data1.rho,"r-",'lineWidth',1.5); hold on;
plot(x,data2.rho,"b-","lineWidth",1.5) ;
xlabel("X/DLe");ylabel("Space Charge Density (rho) ");
legend('Ordinary Plasma',"Dusty Plasma",'Location','west');
title("Compairison of Space Charge density ");

%%%% A PLOT FOR Bohm sheath velocity 
figure(6)

% Create a pseudo y-axis (row index or any physical coordinate)
y1 = 1:size(data1.vxi,1);
y2 = 1:size(data2.vxi,1);

% First subplot (Data 1)
subplot(1,2,1);
contourf(x, y1, data1.vxi, 30, 'LineColor', 'none');  % 30 contour levels
colorbar;
xlabel('x (distance)');
ylabel('Index (pseudo-y)');
title('Data 1: v_x_i distribution');
colormap('turbo');
axis tight;

% Second subplot (Data 2)
subplot(1,2,2);
contourf(x, y2, data2.vxi, 30, 'LineColor', 'none');
colorbar;
xlabel('x (distance)');
ylabel('Index (pseudo-y)');
title('Data 2: v_x_i distribution');
colormap('turbo');
axis tight;

% Optional: consistent color range for fair comparison
cmin = min([data1.vxi(:); data2.vxi(:)]);
cmax = max([data2.vxi(:); data2.vxi(:)]);
subplot(1,2,1); caxis([cmin cmax]);
subplot(1,2,2); caxis([cmin cmax]);

sgtitle('Comparison of v_{xi} for Data1 and Data2');
