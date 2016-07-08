function [ slope_mat] = PhysVMorph2(data_mat, mat_name)
%PhysVMorph plots physiology measurements across many positions in a single
%preparation as function of morphological properties of those same
%positions. The input variable is a compiled data matrix (data_mat) that has the following in each
%column:

% Column 1: Distance from Soma
% Column 2: Diameter in the x-y plane
% Column 3: Neurite Order
% Column 4: Reversal Potential
% Column 5: Mean Response Amplitude at -40 mV

% Each row in the matrix corresponds to a single position. The figure will
% be saved at the filepath. 
% This script also spits out a matrix of the slope values:
% Column 1: Avg Resp V Distance
% Column 2: Avg Resp V Neurite Order
% Column 3: Avg Resp V Diameter
% Column 4: Erev V Distance
% Column 5: Erev V Neurite Order
% Column 6: Erev V Diameter
% Colun 7: Erev V Avg Resp

close all;

load(data_mat); %loads specified data matrix for a given preparation

%slope_mat = zeros(7,2); % Column 1: row is a slope from the corresponding subplot; Column 2: each row is the corresponding MSE 

% Re-name matrix columns for the measured property
dist_vect = mat_name(:,1)
xydiam_vect = mat_name(:,2)
order_vect = mat_name(:,3)
erev_vect = mat_name(:,4)
avgresp_vect = mat_name(:,5)
c_vect = [1 7 10 14 21 31 37 43 48 52 59] 
c = colormap(hsv);
c = c(c_vect(1:length(dist_vect)),:) 

a = 300;

%Plot the colors with position no.
figure(1)
z_vect = zeros(1,length(dist_vect))%length(dist_vect))
scatter(z_vect, 1:1:length(dist_vect), a, c, 'filled', 'LineWidth', 2, 'MarkerEdgeColor','k')
hold on
ylabel('Position No.')
ylim([0 length(dist_vect)])
set(gca,'YDir','reverse')
% rgb2cm
% saveas(gcf, 'PositionColorKey', 'epsc')

%Plot Physiology as a function of Morphology properties
use_vect = find(avgresp_vect > 0); % for the Erev plots only!
%c = linspecer(length(dist_vect)); %linspecer(length(dist_vect(use_vect)));
a = 80;
figure(2)
clf

% MEAN RESPONSE VS DISTANCE
subplot(3,3,1)
scatter(dist_vect', avgresp_vect',a, c, 'filled', 'LineWidth', 2, 'MarkerEdgeColor','k') %(dist_vect(use_vect)', avgresp_vect(use_vect)',a, c, 'filled', 'LineWidth', 2, 'MarkerEdgeColor','k')
xlabel('Distance from Soma (microns)')
xlim([0 800])

%MEAN RESPONSE VS NEURITE ORDER 
subplot(3,3,2)
scatter(order_vect, avgresp_vect, a, c,'filled', 'LineWidth', 2, 'MarkerEdgeColor','k') %(order_vect(use_vect), avgresp_vect(use_vect), a, c,'filled', 'LineWidth', 2, 'MarkerEdgeColor','k')
xlabel('Neurite Order')
xlim([0 30])


% MEAN RESPONSE VS Diameter  
subplot(3,3,3)
scatter(xydiam_vect, avgresp_vect,a, c,'filled', 'LineWidth', 2, 'MarkerEdgeColor','k') %(xydiam_vect(use_vect), avgresp_vect(use_vect),a, c,'filled', 'LineWidth', 2, 'MarkerEdgeColor','k')
xlabel('Diameter in x-y plane (microns)')
xlim([0 15])


%SET Y-AXIS PROPERTIES FOR RESPONSE AMPLITUDE PLOTS
for i = 1:3
    subplot(3,3,i)
    box off
    ylabel('Mean Response (mV)')
    ylim([0 5])
end


%REVERSAL POTENTIAL VS DISTANCE w Linear Fit
use_vect = find(erev_vect < 0)  % remove NaN values
c = c([use_vect], :);
subplot(3,3,4)
scatter(dist_vect(use_vect), erev_vect(use_vect),a, c,'filled', 'LineWidth', 2, 'MarkerEdgeColor','k')
xlabel('Distance from Soma (microns)')
xlim([0 800])
% %Linear Regression
% x = dist_vect(use_vect);
% y = erev_vect(use_vect);
% [R,p] = corrcoef(x,y);
% C = cov(x,y);
% m = R(1,2)*sqrt(C(2,2)/C(1,1));
% b = mean(y) - m*mean(x);
% X = [0 1000];
% Y = m * X + b;
% hold on
% plot(X,Y,'k--', 'LineWidth', 3);
% %Statistics and Such
% Y_fit = m * x + b;
% error = y-Y_fit;
% squarederror = error.*error;% OR, error.^2 is el-by-el raising to 2nd power
% meanerr = mean(squarederror);
% text(400, -75, ['MSE = ' num2str([meanerr])])
% text(400, -80, ['R = ' num2str([R(2,1)])])
% text(400, -85, ['p = ' num2str([p(2,1)])])
% text(400, -90, ['n = ' num2str(length(dist_vect))])
% text(400, -95, ['slope = ' num2str(m)])
% slope_mat(4,1) = m;
% slope_mat(4,2) = meanerr;


%REVERSAL POTENTIAL VS NEURITE ORDER w Linear Fit
subplot(3,3,5)
scatter(order_vect(use_vect), erev_vect(use_vect), a, c,'filled', 'LineWidth', 2, 'MarkerEdgeColor','k')
xlabel('Neurite Order')
xlim([0 30])
% %Linear Regression
% x = order_vect(use_vect);
% y = erev_vect(use_vect);
% [R,p] = corrcoef(x,y);
% C = cov(x,y);
% m = R(1,2)*sqrt(C(2,2)/C(1,1));
% b = mean(y) - m*mean(x);
% X = [0 60];
% Y = m * X + b;
% hold on
% plot(X,Y,'k--', 'LineWidth', 3);
% %Statistics and Such
% Y_fit = m * x + b;
% error = y-Y_fit;
% squarederror = error.*error;% OR, error.^2 is el-by-el raising to 2nd power
% meanerr = mean(squarederror);
% text(25, -75, ['MSE = ' num2str([meanerr])])
% text(25, -80, ['R = ' num2str([R(2,1)])])
% text(25, -85, ['p = ' num2str([p(2,1)])])
% text(25, -90, ['n = ' num2str(length(order_vect(use_vect)))])
% text(25, -95, ['slope = ' num2str(m)])
% slope_mat(5,1) = m;
% slope_mat(5,2) = meanerr;


%REVERSAL POTENTIAL VS DIAMETER w Linear Fit
subplot(3,3,6)
scatter(xydiam_vect(use_vect), erev_vect(use_vect),a, c,'filled', 'LineWidth', 2, 'MarkerEdgeColor','k')
xlabel('Diameter in x-y plane (microns)')
xlim([0 15])
% %Linear Regression
% x = xydiam_vect(use_vect);
% y = erev_vect(use_vect);
% [R,p] = corrcoef(x,y);
% C = cov(x,y);
% m = R(1,2)*sqrt(C(2,2)/C(1,1));
% b = mean(y) - m*mean(x);
% X = [0 60];
% Y = m * X + b;
% hold on
% plot(X,Y,'k--', 'LineWidth', 3);
% %Statistics and Such
% Y_fit = m * x + b;
% error = y-Y_fit;
% squarederror = error.*error;% OR, error.^2 is el-by-el raising to 2nd power
% meanerr = mean(squarederror);
% text(9, -75, ['MSE = ' num2str([meanerr])])
% text(9, -80, ['R = ' num2str([R(2,1)])])
% text(9, -85, ['p = ' num2str([p(2,1)])])
% text(9, -90, ['n = ' num2str(length(xydiam_vect(use_vect)))])
% text(9, -95, ['slope = ' num2str(m)])
% slope_mat(6,1) = m;
% slope_mat(6,2) = meanerr;

%SET Y-AXIS PROPERTIES FOR REVESRAL POTENTIAL PLOTS
for i = 4:6
    subplot(3,3,i)
    box off
    ylabel('Reversal Potential(mV)')
    ylim([-100 -50])
end
%rgb2cm
%saveas(gcf, '878_045PlotMorphVPhys', 'epsc')


% % Plot comparisons across morphology properties
% %c = distinguishable_colors(length(dist_vect));
% c = linspecer(length(dist_vect));
% subplot(4,3,7)
% scatter(dist_vect, xydiam_vect, a, c, 'filled', 'LineWidth', 2, 'MarkerEdgeColor', 'k')
% xlabel('Distance from Soma (microns)')
% ylabel('Diameter in x-y plane (microns)')
% subplot(4,3,8)
% scatter(dist_vect, order_vect, a, c, 'filled', 'LineWidth', 2, 'MarkerEdgeColor', 'k')
% xlabel('Distance from Soma (microns)')
% ylabel('Neurite Order')
% subplot(4,3,9)
% scatter(order_vect, xydiam_vect, a, c, 'filled', 'LineWidth', 2, 'MarkerEdgeColor', 'k')
% xlabel('Neurite Order')
% ylabel('Diameter in x-y plane (microns)')
% rgb2cm
% %saveas(gcf, '878_045PlotMorphVMorph', 'epsc')


% Plot comparisons across physiological properties
%use_vect = find(erev_vect < 0) % remove NaN values
%use_vect = use_vect(find(avgresp_vect > 0))
%c = linspecer(length(erev_vect(use_vect)));
subplot(3,3,7)
scatter(avgresp_vect(use_vect), erev_vect(use_vect), a, c, 'filled', 'LineWidth', 2, 'MarkerEdgeColor', 'k')
% hold on
% for i = 1:length(dist_vect)
%     line([avgresp_vect(i) - std_vect(i) avgresp_vect(i)+ std_vect(i)], [erev_vect(i) erev_vect(i)])
% end
xlabel('Mean Response (mV)')
ylabel('Reversal Potential(mV)')
ylim([-100 -50])
xlim([0 5])
% %Linear Regression
% x = avgresp_vect(use_vect);
% y = erev_vect(use_vect);
% [R,p] = corrcoef(x,y);
% C = cov(x,y);
% m = R(1,2)*sqrt(C(2,2)/C(1,1));
% b = mean(y) - m*mean(x);
% X = [0 5];
% Y = m * X + b;
% hold on
% plot(X,Y,'k--', 'LineWidth', 3);
% %Statistics and Such
% Y_fit = m * x + b;
% error = y-Y_fit;
% squarederror = error.*error;% OR, error.^2 is el-by-el raising to 2nd power
% meanerr = mean(squarederror);
% text(2, -75, ['MSE = ' num2str([meanerr])])
% text(2, -80, ['R = ' num2str([R(2,1)])])
% text(2, -85, ['p = ' num2str([p(2,1)])])
% text(2, -90, ['n = ' num2str(length(avgresp_vect(use_vect)))])
% text(2, -95, ['slope = ' num2str(m)])
% slope_mat(7,1) = m;
% slope_mat(7,2) = meanerr;

rgb2cm
saveas(gcf, 'Plots_878_057hsv_new', 'epsc')


end

