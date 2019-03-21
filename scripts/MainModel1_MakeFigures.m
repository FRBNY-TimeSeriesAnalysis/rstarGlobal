% This file replicates the plots for the "Baseline" model. See figures
% 1, 2, 3, and 4, 5, A1, A2.

%% Plot preliminaries: load, sort, and get quantiles

close all;

addpath('Routines')
load('../results/OutputModel1.mat')

figpath   = '../figures/';
appenpath = [figpath 'appendix/'];

Quant = [0.025 0.160 0.500 0.840  0.975];
M = size(CommonTrends, 3);

set(0,'defaultAxesFontName', 'Times');
set(0,'DefaultTextInterpreter', 'latex')
set(0, 'DefaultAxesFontSize',15)
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth',1.5)
setappdata(0, 'defaultAxesXTickFontSize', 1)
setappdata(0, 'defaultAxesYTickFontSize', 1)


% ------------------------- Trends ---------------------------------------

Rshort_bar    = squeeze(CommonTrends(:, 1,:));
Pi_bar        = squeeze(CommonTrends(:, 2,:));
Ts_bar        = squeeze(CommonTrends(:, 3,:));
Rlong_bar     = Rshort_bar + Ts_bar;

Rshort_bar_us_idio = squeeze(CommonTrends(:, 4,:));
Rshort_bar_de_idio = squeeze(CommonTrends(:, 5,:));
Rshort_bar_uk_idio = squeeze(CommonTrends(:, 6,:));
Rshort_bar_fr_idio = squeeze(CommonTrends(:, 7,:));
Rshort_bar_ca_idio = squeeze(CommonTrends(:, 8,:));
Rshort_bar_it_idio = squeeze(CommonTrends(:, 9,:));
Rshort_bar_jp_idio = squeeze(CommonTrends(:,10,:));

Pi_bar_us_idio = squeeze(CommonTrends(:,11,:));
Pi_bar_de_idio = squeeze(CommonTrends(:,12,:));
Pi_bar_uk_idio = squeeze(CommonTrends(:,13,:));
Pi_bar_fr_idio = squeeze(CommonTrends(:,14,:));
Pi_bar_ca_idio = squeeze(CommonTrends(:,15,:));
Pi_bar_it_idio = squeeze(CommonTrends(:,16,:));
Pi_bar_jp_idio = squeeze(CommonTrends(:,17,:));

Ts_bar_us_idio = squeeze(CommonTrends(:,18,:));
Ts_bar_de_idio = squeeze(CommonTrends(:,19,:));
Ts_bar_uk_idio = squeeze(CommonTrends(:,20,:));
Ts_bar_fr_idio = squeeze(CommonTrends(:,21,:));
Ts_bar_ca_idio = squeeze(CommonTrends(:,22,:));
Ts_bar_it_idio = squeeze(CommonTrends(:,23,:));
Ts_bar_jp_idio = squeeze(CommonTrends(:,24,:));

Rlong_bar_us_idio = Rshort_bar_us_idio + Ts_bar_us_idio;
Rlong_bar_de_idio = Rshort_bar_de_idio + Ts_bar_de_idio;
Rlong_bar_uk_idio = Rshort_bar_uk_idio + Ts_bar_uk_idio;
Rlong_bar_fr_idio = Rshort_bar_fr_idio + Ts_bar_fr_idio;
Rlong_bar_ca_idio = Rshort_bar_ca_idio + Ts_bar_ca_idio;
Rlong_bar_it_idio = Rshort_bar_it_idio + Ts_bar_it_idio;
Rlong_bar_jp_idio = Rshort_bar_jp_idio + Ts_bar_jp_idio;

Rshort_bar_us = repmat(transpose(squeeze(CC(1,1,:))), T, 1) .* Rshort_bar + Rshort_bar_us_idio;
Rshort_bar_de = repmat(transpose(squeeze(CC(2,1,:))), T, 1) .* Rshort_bar + Rshort_bar_de_idio;
Rshort_bar_uk = repmat(transpose(squeeze(CC(3,1,:))), T, 1) .* Rshort_bar + Rshort_bar_uk_idio;
Rshort_bar_fr = repmat(transpose(squeeze(CC(4,1,:))), T, 1) .* Rshort_bar + Rshort_bar_fr_idio;
Rshort_bar_ca = repmat(transpose(squeeze(CC(5,1,:))), T, 1) .* Rshort_bar + Rshort_bar_ca_idio;
Rshort_bar_it = repmat(transpose(squeeze(CC(6,1,:))), T, 1) .* Rshort_bar + Rshort_bar_it_idio;
Rshort_bar_jp = repmat(transpose(squeeze(CC(7,1,:))), T, 1) .* Rshort_bar + Rshort_bar_jp_idio;

Rlong_bar_us = Rlong_bar + Rlong_bar_us_idio;
Rlong_bar_de = Rlong_bar + Rlong_bar_de_idio;
Rlong_bar_uk = Rlong_bar + Rlong_bar_uk_idio;
Rlong_bar_fr = Rlong_bar + Rlong_bar_fr_idio;
Rlong_bar_ca = Rlong_bar + Rlong_bar_ca_idio;
Rlong_bar_it = Rlong_bar + Rlong_bar_it_idio;
Rlong_bar_jp = Rlong_bar + Rlong_bar_jp_idio;

Ts_bar_us    = Ts_bar + Ts_bar_us_idio;
Ts_bar_de    = Ts_bar + Ts_bar_de_idio;
Ts_bar_uk    = Ts_bar + Ts_bar_uk_idio;
Ts_bar_fr    = Ts_bar + Ts_bar_fr_idio;
Ts_bar_ca    = Ts_bar + Ts_bar_ca_idio;
Ts_bar_it    = Ts_bar + Ts_bar_it_idio;
Ts_bar_jp    = Ts_bar + Ts_bar_jp_idio;

Pi_bar_us = repmat(transpose(squeeze(CC(1,2,:))), T, 1) .* Pi_bar + Pi_bar_us_idio;
Pi_bar_de = repmat(transpose(squeeze(CC(2,2,:))), T, 1) .* Pi_bar + Pi_bar_de_idio;
Pi_bar_uk = repmat(transpose(squeeze(CC(3,2,:))), T, 1) .* Pi_bar + Pi_bar_uk_idio;
Pi_bar_fr = repmat(transpose(squeeze(CC(4,2,:))), T, 1) .* Pi_bar + Pi_bar_fr_idio;
Pi_bar_ca = repmat(transpose(squeeze(CC(5,2,:))), T, 1) .* Pi_bar + Pi_bar_ca_idio;
Pi_bar_it = repmat(transpose(squeeze(CC(6,2,:))), T, 1) .* Pi_bar + Pi_bar_it_idio;
Pi_bar_jp = repmat(transpose(squeeze(CC(7,2,:))), T, 1) .* Pi_bar + Pi_bar_jp_idio;

% --------------------------- Sorted trends -------------------------------

sRshort_bar     = sort(Rshort_bar,2);
sPi_bar         = sort(Pi_bar,2);
sTs_bar         = sort(Ts_bar,2);
sRlong_bar      = sort(Rlong_bar,2);

sRshort_bar_us = sort(Rshort_bar_us,2);
sRshort_bar_de = sort(Rshort_bar_de,2);
sRshort_bar_uk = sort(Rshort_bar_uk,2);
sRshort_bar_fr = sort(Rshort_bar_fr,2);
sRshort_bar_ca = sort(Rshort_bar_ca,2);
sRshort_bar_it = sort(Rshort_bar_it,2);
sRshort_bar_jp = sort(Rshort_bar_jp,2);

sRlong_bar_us = sort(Rlong_bar_us,2);
sRlong_bar_de = sort(Rlong_bar_de,2);
sRlong_bar_uk = sort(Rlong_bar_uk,2);
sRlong_bar_fr = sort(Rlong_bar_fr,2);
sRlong_bar_ca = sort(Rlong_bar_ca,2);
sRlong_bar_it = sort(Rlong_bar_it,2);
sRlong_bar_jp = sort(Rlong_bar_jp,2);

sPi_bar_us = sort(Pi_bar_us,2);
sPi_bar_de = sort(Pi_bar_de,2);
sPi_bar_uk = sort(Pi_bar_uk,2);
sPi_bar_fr = sort(Pi_bar_fr,2);
sPi_bar_ca = sort(Pi_bar_ca,2);
sPi_bar_it = sort(Pi_bar_it,2);
sPi_bar_jp = sort(Pi_bar_jp,2);

sRshort_bar_us_idio = sort(Rshort_bar_us_idio,2);
sRshort_bar_de_idio = sort(Rshort_bar_de_idio,2);
sRshort_bar_uk_idio = sort(Rshort_bar_uk_idio,2);
sRshort_bar_fr_idio = sort(Rshort_bar_fr_idio,2);
sRshort_bar_ca_idio = sort(Rshort_bar_ca_idio,2);
sRshort_bar_it_idio = sort(Rshort_bar_it_idio,2);
sRshort_bar_jp_idio = sort(Rshort_bar_jp_idio,2);

sTs_bar_us = sort(Ts_bar_us,2);
sTs_bar_de = sort(Ts_bar_de,2);
sTs_bar_uk = sort(Ts_bar_uk,2);
sTs_bar_fr = sort(Ts_bar_fr,2);
sTs_bar_ca = sort(Ts_bar_ca,2);
sTs_bar_it = sort(Ts_bar_it,2);
sTs_bar_jp = sort(Ts_bar_jp,2);

% ------------------- quantiles of the trends -----------------------------

qRshort_bar     = sRshort_bar(:,ceil(Quant*M));
qPi_bar         = sPi_bar(:,ceil(Quant*M));
qTs_bar         = sTs_bar(:,ceil(Quant*M));
qRlong_bar      = sRlong_bar(:,ceil(Quant*M));

qRshort_bar_us = sRshort_bar_us(:,ceil(Quant*M));
qRshort_bar_de = sRshort_bar_de(:,ceil(Quant*M));
qRshort_bar_uk = sRshort_bar_uk(:,ceil(Quant*M));
qRshort_bar_fr = sRshort_bar_fr(:,ceil(Quant*M));
qRshort_bar_ca = sRshort_bar_ca(:,ceil(Quant*M));
qRshort_bar_it = sRshort_bar_it(:,ceil(Quant*M));
qRshort_bar_jp = sRshort_bar_jp(:,ceil(Quant*M));

qRlong_bar_us = sRlong_bar_us(:,ceil(Quant*M));
qRlong_bar_de = sRlong_bar_de(:,ceil(Quant*M));
qRlong_bar_uk = sRlong_bar_uk(:,ceil(Quant*M));
qRlong_bar_fr = sRlong_bar_fr(:,ceil(Quant*M));
qRlong_bar_ca = sRlong_bar_ca(:,ceil(Quant*M));
qRlong_bar_it = sRlong_bar_it(:,ceil(Quant*M));
qRlong_bar_jp = sRlong_bar_jp(:,ceil(Quant*M));

qPi_bar_us = sPi_bar_us(:,ceil(Quant*M));
qPi_bar_de = sPi_bar_de(:,ceil(Quant*M));
qPi_bar_uk = sPi_bar_uk(:,ceil(Quant*M));
qPi_bar_fr = sPi_bar_fr(:,ceil(Quant*M));
qPi_bar_ca = sPi_bar_ca(:,ceil(Quant*M));
qPi_bar_it = sPi_bar_it(:,ceil(Quant*M));
qPi_bar_jp = sPi_bar_jp(:,ceil(Quant*M));

qTs_bar_us = sTs_bar_us(:,ceil(Quant*M));
qTs_bar_de = sTs_bar_de(:,ceil(Quant*M));
qTs_bar_uk = sTs_bar_uk(:,ceil(Quant*M));
qTs_bar_fr = sTs_bar_fr(:,ceil(Quant*M));
qTs_bar_ca = sTs_bar_ca(:,ceil(Quant*M));
qTs_bar_it = sTs_bar_it(:,ceil(Quant*M));
qTs_bar_jp = sTs_bar_jp(:,ceil(Quant*M));

qRshort_bar_us_idio = sRshort_bar_us_idio(:,ceil(Quant*M));
qRshort_bar_de_idio = sRshort_bar_de_idio(:,ceil(Quant*M));
qRshort_bar_uk_idio = sRshort_bar_uk_idio(:,ceil(Quant*M));
qRshort_bar_fr_idio = sRshort_bar_fr_idio(:,ceil(Quant*M));
qRshort_bar_ca_idio = sRshort_bar_ca_idio(:,ceil(Quant*M));
qRshort_bar_it_idio = sRshort_bar_it_idio(:,ceil(Quant*M));
qRshort_bar_jp_idio = sRshort_bar_jp_idio(:,ceil(Quant*M));



%% Figure 1: Trends in Global and U.S. Real Rates: 1870-2016, Baseline Model

figure

PlotStatesShaded(Year, qRshort_bar)
hold on
plot(Year, qRshort_bar_us(:, 3), ...
    'k:', 'LineWidth', 1.5);  % r-bar US

axis([Year(1) Year(end) -3 6])
% title('$\bar{r}^w_t$ and $\bar{r}_{US, t}$', 'Interpreter', 'latex')
box on

printpdf(gcf, [figpath 'fig1-Model1_Rshortbar-us.pdf'])


%% Figure 2: Trends in Global Real Rates Under Alternative Priors for the 
% Standard Deviation of Innovations to the Trend and Decadal Moving Averages

figure

% Load results using disperse prior for variance of innovations to trends
temp = load('../results/OutputModel1_var01.mat', 'CommonTrends');
CommonTrends_var01 = temp.CommonTrends;
Rshort_bar_var01   = squeeze(CommonTrends_var01(:, 1,:));
sRshort_bar_var01  = sort(Rshort_bar_var01,2);
qRshort_bar_var01  = sRshort_bar_var01(:,ceil(Quant*M));


% Compute real interest rate
rir = [Stir_us-Infl_us, Stir_de-Infl_de, Stir_uk-Infl_uk, ...
    Stir_fr-Infl_fr, Stir_ca-Infl_ca, Stir_it-Infl_it, Stir_jp-Infl_jp];
rir(abs(rir) > 30) = NaN;  % Remove extreme observations

rir_ma = nan(size(rir));
h = 5;  % Centered moving average (plus/minus h years)

for iVar = 1:size(rir,2)
    z_i  = rir(:,iVar);
    rir_ma(:,iVar) = ma_centered(z_i, h);
end

rir_ma_world = nanmean(rir_ma, 2);  % Take cross-sectional average

bands_new(Year, qRshort_bar); hold on;
bands_new(Year, [], qRshort_bar_var01);
bands_new(Year, qRshort_bar);
plot(Year, rir_ma_world, 'LineWidth', 2, 'LineStyle', '-',...
    'Color', 0.75 * [0.9544, 0.0780, 0.1840]);

axis([Year(1) Year(end) -12 12])

printpdf(gcf, [figpath 'fig2-Model1_var01_Rbar-MA.pdf'])


%% Figure 3a: Trends and Observables for Short-Term Real Rates, Baseline Model

figure

PlotStatesShaded(Year, qRshort_bar);
hold on;
p_us = plot(Year, Stir_us - Infl_us, 'k:', 'LineWidth', 1);
p_de = plot(Year, Stir_de - Infl_de, 'b:', 'LineWidth', 1);
p_uk = plot(Year, Stir_uk - Infl_uk, 'c:', 'LineWidth', 1);
p_fr = plot(Year, Stir_fr - Infl_fr, 'y:', 'LineWidth', 1);
p_ca = plot(Year, Stir_ca - Infl_ca, 'r:', 'LineWidth', 1);
p_it = plot(Year, Stir_it - Infl_it, 'g:', 'LineWidth', 1);
p_jp = plot(Year, Stir_jp - Infl_jp, 'm:', 'LineWidth', 1);

axis([Year(1) Year(end) -6 12])
yticks(-5:5:10)

legend([p_us p_de p_uk p_fr p_ca p_it p_jp],...
    {'us', 'de', 'uk', 'fr', 'ca', 'it', 'jp'},...
    'Interpreter', 'latex',...
    'Location',    'SouthOutside',...
    'FontSize',    12,...
    'Orientation', 'horizontal');
legend boxoff;
% title('$\bar{r}^w_t$ and $R_{i,t} - \pi_{i,t}$', 'Interpreter', 'latex')
box on;

printpdf(gcf, [figpath 'fig3a-Model1_Rshortbar-obs.pdf'])


%% Figure 3b: Trends and Observables for Short-Term Real Rates, Baseline Model

figure
hold on

p_us = plot(Year, qRshort_bar_us(:,3), 'k:', 'LineWidth', 2); hold on;
p_de = plot(Year, qRshort_bar_de(:,3), 'b:', 'LineWidth', 2);
p_uk = plot(Year, qRshort_bar_uk(:,3), 'c:', 'LineWidth', 2);
p_fr = plot(Year, qRshort_bar_fr(:,3), 'y:', 'LineWidth', 2);
p_ca = plot(Year, qRshort_bar_ca(:,3), 'r:', 'LineWidth', 2);
p_it = plot(Year, qRshort_bar_it(:,3), 'g:', 'LineWidth', 2);
p_jp = plot(Year, qRshort_bar_jp(:,3), 'm:', 'LineWidth', 2);
plot(Year, qRshort_bar(:, 3), 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');


hline                  = refline(0);
hline.HandleVisibility = 'off';
hline.Color            = 'k';

xlim([Year(1) Year(end)])
ylim([-3 6])
yticks(-2:2:6)

legend([p_us p_de p_uk p_fr p_ca p_it p_jp],...
    {'us', 'de', 'uk', 'fr', 'ca', 'it', 'jp'},...
    'Interpreter', 'latex',...
    'Location',    'SouthOutside',...
    'FontSize',    12,...
    'Orientation', 'horizontal');

%title('$\bar{r}^w_t$ and $\bar{r}_{i,t}$', 'Interpreter','latex')

legend boxoff;
box on

printpdf(gcf, [figpath 'fig3b-Model1_Rshortbar-countries.pdf'])


%% Compute table 1 main model values

Table1 = struct;

Table1(1).Delta_rw = CommonTrends(Year == 2016, 1, :) - CommonTrends(Year == 1980, 1, :);
Table1(1).Years = [1980 2016];

Table1(2).Delta_rw = CommonTrends(Year == 2016, 1, :) - CommonTrends(Year == 1990, 1, :);
Table1(2).Years = [1990 2016];

Table1(3).Delta_rw = CommonTrends(Year == 2016, 1, :) - CommonTrends(Year == 1997, 1, :);
Table1(3).Years = [1997 2016];

clc

for j = 1:size(Table1, 2)
    disp('')
    disp(['-----------------[' num2str(Table1(j).Years) ']------------------'])
    disp(['Median:' num2str(quantile(Table1(j).Delta_rw, .5))])
    disp(['90% posterior coverage: ' ...
        '[' num2str(quantile(Table1(j).Delta_rw, .05)) ', ' num2str(quantile(Table1(j).Delta_rw, .95)) ']'])
end


%% Figure 4a: Trends and Observables for Inflation, Baseline Model

figure;

PlotStatesShaded(Year, qPi_bar); hold on;  % pi-bar world
p_us = plot(Year, Infl_us, 'k:', 'LineWidth', 1); hold on;
p_de = plot(Year, Infl_de, 'b:', 'LineWidth', 1);
p_uk = plot(Year, Infl_uk, 'c:', 'LineWidth', 1);
p_fr = plot(Year, Infl_fr, 'y:', 'LineWidth', 1);
p_ca = plot(Year, Infl_ca, 'r:', 'LineWidth', 1);
p_it = plot(Year, Infl_it, 'g:', 'LineWidth', 1);
p_jp = plot(Year, Infl_jp, 'm:', 'LineWidth', 1);

xlim([1870 2017])
ylim([-3 15])
yticks(0:5:15)

legend([p_us p_de p_uk p_fr p_ca p_it p_jp], ...
    {'us', 'de', 'uk', 'fr', 'ca', 'it', 'jp'}, 'Location', 'southoutside',...
    'Orientation', 'horizontal', 'Box', 'off')
% title('$\bar{\pi}^w_t$ and $\pi_{it}$')

box on

printpdf(gcf, [figpath 'fig4a-Model1_Pibar-obs.pdf'])

%% Figure 4b: Trends and Observables for Inflation, Baseline Model

figure;

median = quantile(CommonTrends(:,2,:), .5, 3);

plotMedian = plot(Year, median, 'k--', 'LineWidth', 2, ...
    'HandleVisibility', 'off');
hold on; box on;

p_us = plot(Year, qPi_bar_us(:, 3), 'k:', 'LineWidth', 1.5); hold on;
p_de = plot(Year, qPi_bar_de(:, 3), 'b:', 'LineWidth', 1.5);
p_uk = plot(Year, qPi_bar_uk(:, 3), 'c:', 'LineWidth', 1.5);
p_fr = plot(Year, qPi_bar_fr(:, 3), 'y:', 'LineWidth', 1.5);
p_ca = plot(Year, qPi_bar_ca(:, 3), 'r:', 'LineWidth', 1.5);
p_it = plot(Year, qPi_bar_it(:, 3), 'g:', 'LineWidth', 1.5);
p_jp = plot(Year, qPi_bar_jp(:, 3), 'm:', 'LineWidth', 1.5);


hline = refline(0);
hline.HandleVisibility = 'off';
hline.Color = 'k';
axis([Year(1) Year(end) -3 15]);

%title('$\bar{\pi}^w_t$ and $\bar{\pi}_{i,t}$', 'Interpreter', 'latex')

legend([p_us p_de p_uk p_fr p_ca p_it p_jp],...
    {'us', 'de', 'uk', 'fr', 'ca', 'it', 'jp'}, 'Location', 'southoutside',...
    'Orientation', 'horizontal', 'Box', 'off')
box on
yticks(0:5:15)

printpdf(gcf, [figpath 'fig4b-Model1_Pibar-countries.pdf'])


%% Figure 5: MY regression fit

fSize = 15;

MY = xlsread('../indata/Data_MY.xlsx');  % Load MY for each country
MY_us = MY(:,2);
MY_de = MY(:,3);
MY_uk = MY(:,4);
MY_fr = MY(:,5);
MY_ca = MY(:,6);
MY_it = MY(:,7);
MY_jp = MY(:,8);
MY_G7 = MY(:,9);

% Compute fitted values
fit_us_MY = [ones(T, 1) MY_us] * regress(qRshort_bar_us(:, 3), [ones(T, 1) MY_us]);
fit_de_MY = [ones(T, 1) MY_de] * regress(qRshort_bar_de(:, 3), [ones(T, 1) MY_de]);
fit_uk_MY = [ones(T, 1) MY_uk] * regress(qRshort_bar_uk(:, 3), [ones(T, 1) MY_uk]);
fit_fr_MY = [ones(T, 1) MY_fr] * regress(qRshort_bar_fr(:, 3), [ones(T, 1) MY_fr]);
fit_ca_MY = [ones(T, 1) MY_ca] * regress(qRshort_bar_ca(:, 3), [ones(T, 1) MY_ca]);
fit_it_MY = [ones(T, 1) MY_it] * regress(qRshort_bar_it(:, 3), [ones(T, 1) MY_it]);
fit_jp_MY = [ones(T, 1) MY_jp] * regress(qRshort_bar_jp(:, 3), [ones(T, 1) MY_jp]);
fit_G7_MY = [ones(T, 1) MY_G7] * regress(qRshort_bar   (:, 3), [ones(T, 1) MY_G7]);

figure
p_wd = PlotStatesShaded(Year, qRshort_bar); hold on;
plot(Year, fit_G7_MY, 'k-')
%title('World')
axis([Year(1) Year(end) -3 6])
set(gca, 'FontSize', fSize)
printpdf(gcf, [figpath 'fig5-Model1_Rbar-countries_MY-fitted-common.pdf'])

figure
p_us = PlotStatesShaded(Year, qRshort_bar_us); hold on;
plot(Year, fit_us_MY, 'k-')
%title('US')
axis([Year(1) Year(end) -3 6])
set(gca, 'FontSize', fSize)
printpdf(gcf, [figpath 'fig5-Model1_Rbar-countries_MY-fitted-us.pdf'])

figure
p_de = PlotStatesShaded(Year, qRshort_bar_de); hold on;
p_de.Color = 'b';
plot(Year, fit_de_MY, 'k-')
%title('Germany')
axis([Year(1) Year(end) -3 6])
set(gca, 'FontSize', fSize)
printpdf(gcf, [figpath 'fig5-Model1_Rbar-countries_MY-fitted-de.pdf'])

figure
p_uk = PlotStatesShaded(Year, qRshort_bar_uk); hold on;
p_uk.Color = 'c';
plot(Year, fit_uk_MY, 'k-')
%title('U.K.')
axis([Year(1) Year(end) -3 6])
set(gca, 'FontSize', fSize)
printpdf(gcf, [figpath 'fig5-Model1_Rbar-countries_MY-fitted-uk.pdf'])

figure
p_fr = PlotStatesShaded(Year, qRshort_bar_fr); hold on;
p_fr.Color = 'y';
plot(Year, fit_fr_MY, 'k-')
%title('France')
axis([Year(1) Year(end) -3 6])
set(gca, 'FontSize', fSize)
printpdf(gcf, [figpath 'fig5-Model1_Rbar-countries_MY-fitted-fr.pdf'])

figure
p_ca = PlotStatesShaded(Year, qRshort_bar_ca); hold on;
p_ca.Color = 'r';
plot(Year, fit_ca_MY, 'k-')
%title('Canada')
axis([Year(1) Year(end) -3 6])
set(gca, 'FontSize', fSize)
printpdf(gcf, [figpath 'fig5-Model1_Rbar-countries_MY-fitted-ca.pdf'])

figure
p_it = PlotStatesShaded(Year, qRshort_bar_it); hold on;
p_it.Color = 'g';
plot(Year, fit_it_MY, 'k-')
%title('Italy')
axis([Year(1) Year(end) -3 6])
set(gca, 'FontSize', fSize)
printpdf(gcf, [figpath 'fig5-Model1_Rbar-countries_MY-fitted-it.pdf'])


figure
p_jp = PlotStatesShaded(Year, qRshort_bar_jp); hold on;
p_jp.Color = 'm';
plot(Year, fit_jp_MY, 'k-')
%title('Japan')
axis([Year(1) Year(end) -3 6])
set(gca, 'FontSize', fSize)
printpdf(gcf, [figpath 'fig5-Model1_Rbar-countries_MY-fitted-jp.pdf'])




%% ---------------------- Appendix figures --------------------------------


%% A1: Trends and Observables for Term Spreads, Baseline Model

set(0, 'DefaultAxesFontSize',15)

f = figure;
filename = 'Tsbar';
h = PlotStatesShaded(Year, qTs_bar);
hold on; box on;
axis([Year(1) Year(end) -2 4]);

p_us = plot(Year, Ltir_us-Stir_us, 'k:', 'LineWidth', 1); hold on;
p_de = plot(Year, Ltir_de-Stir_de, 'b:', 'LineWidth', 1);
p_uk = plot(Year, Ltir_uk-Stir_uk, 'c:', 'LineWidth', 1);
p_fr = plot(Year, Ltir_fr-Stir_fr, 'y:', 'LineWidth', 1);
p_ca = plot(Year, Ltir_ca-Stir_ca, 'r:', 'LineWidth', 1);
p_it = plot(Year, Ltir_it-Stir_it, 'g:', 'LineWidth', 1);
p_jp = plot(Year, Ltir_jp-Stir_jp, 'm:', 'LineWidth', 1);

legend([p_us p_de p_uk p_fr p_ca p_it p_jp],...
    {'us', 'de', 'uk', 'fr', 'ca', 'it', 'jp'},...
    'Interpreter','latex',...
    'Location','SouthOutside',...
    'FontSize',12,'Orientation', 'horizontal'); 

legend boxoff;

%title('$\overline{ts}^w_t$ and $R^L_{i,t} - R_{i,t}$', 'Interpreter', 'latex')

printpdf(gcf, [appenpath 'figa1a-Model1_Tsbar.pdf'])


%% A1: Trends and Observables for Term Spreads, Baseline Model


f = figure;
filename = 'Tsbar-countries';
h = PlotStatesShaded(Year, qTs_bar(:,3));
hold on; box on;
axis([Year(1) Year(end) -2 4]);
p_us = plot(Year, qTs_bar_us(:,3), 'k:', 'LineWidth', 2); hold on;
p_de = plot(Year, qTs_bar_de(:,3), 'b:', 'LineWidth', 2);
p_uk = plot(Year, qTs_bar_uk(:,3), 'c:', 'LineWidth', 2);
p_fr = plot(Year, qTs_bar_fr(:,3), 'y:', 'LineWidth', 2);
p_ca = plot(Year, qTs_bar_ca(:,3), 'r:', 'LineWidth', 2);
p_it = plot(Year, qTs_bar_it(:,3), 'g:', 'LineWidth', 2);
p_jp = plot(Year, qTs_bar_jp(:,3), 'm:', 'LineWidth', 2);

legend([p_us p_de p_uk p_fr p_ca p_it p_jp],...
    {'us', 'de', 'uk', 'fr', 'ca', 'it', 'jp'},...
    'Interpreter','latex',...
    'Location','SouthOutside',...
    'FontSize',12,'Orientation', 'horizontal'); 

legend boxoff;
%title('$\overline{ts}^w_t$ and $\overline{ts}_{i,t}$', 'Interpreter', 'latex')
printpdf(gcf, [appenpath 'figa1b-Model1_Tsbar-countries.pdf'])


%% A2: Country-Specfic Trends r_it and Observables, Baseline Model

fSize = 15;  % Font size

f = figure;

Rshort_country_average = ...
    mean([Stir_us-Infl_us, Stir_de-Infl_de, Stir_uk-Infl_uk,...
          Stir_fr-Infl_fr, Stir_ca-Infl_ca, Stir_it-Infl_it, ...
          Stir_jp-Infl_jp], 2, 'omitnan');
      
h = PlotStatesShaded(Year, qRshort_bar);
hold on; box on; axis([Year(1) Year(end) -3 6]);
%title('World', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa2-Model1_Rshort-countries_trend-idio_obs-average-common.pdf'])

f = figure;
h = PlotStatesShaded(Year, qRshort_bar_us_idio); hold on;
p_us = plot(Year, qRshort_bar_us_idio(:,3), 'k:', 'LineWidth', 2);
plot(Year, Stir_us-Infl_us-Rshort_country_average, 'k:', 'LineWidth', 1);
axis([Year(1) Year(end) -3 6]);
%title('U.S.', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa2-Model1_Rshort-countries_trend-idio_obs-average-us.pdf'])

figure;
h = PlotStatesShaded(Year, qRshort_bar_de_idio); hold on;
h.Color = 'b';
h.LineStyle = '--';
p_de = plot(Year, qRshort_bar_de_idio(:,3), 'b:', 'LineWidth', 2);
plot(Year, Stir_de-Infl_de-Rshort_country_average, 'b:', 'LineWidth', 1);
 axis([Year(1) Year(end) -3 6]);
%title('Germany', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa2-Model1_Rshort-countries_trend-idio_obs-average-de.pdf'])

figure;
h = PlotStatesShaded(Year, qRshort_bar_uk_idio); hold on;
h.Color = 'c';
p_uk = plot(Year, qRshort_bar_uk_idio(:,3), 'c:', 'LineWidth', 2);
plot(Year, Stir_uk-Infl_uk-Rshort_country_average, 'c:', 'LineWidth', 1);
 axis([Year(1) Year(end) -3 6]);
%title('U.K.', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa2-Model1_Rshort-countries_trend-idio_obs-average-uk.pdf'])

figure;
h = PlotStatesShaded(Year, qRshort_bar_fr_idio); hold on;
h.Color = 'y';
p_fr = plot(Year, qRshort_bar_fr_idio(:,3), 'y:', 'LineWidth', 2);
plot(Year, Stir_fr-Infl_fr-Rshort_country_average, 'y:', 'LineWidth', 1);
axis([Year(1) Year(end) -3 6]);
%title('France', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa2-Model1_Rshort-countries_trend-idio_obs-average-fr.pdf'])

figure;
h = PlotStatesShaded(Year, qRshort_bar_ca_idio); hold on;
h.Color = 'r';
p_ca = plot(Year, qRshort_bar_ca_idio(:,3), 'r:', 'LineWidth', 2);
plot(Year, Stir_ca-Infl_ca-Rshort_country_average, 'r:', 'LineWidth', 1);
axis([Year(1) Year(end) -3 6]);
%title('Canada', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa2-Model1_Rshort-countries_trend-idio_obs-average-ca.pdf'])

figure;
h = PlotStatesShaded(Year, qRshort_bar_it_idio); hold on;
h.Color = 'g';
p_it = plot(Year, qRshort_bar_it_idio(:,3), 'g:', 'LineWidth', 2);
plot(Year, Stir_it-Infl_it-Rshort_country_average, 'g:', 'LineWidth', 1);
axis([Year(1) Year(end) -3 6]);
%title('Italy', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa2-Model1_Rshort-countries_trend-idio_obs-average-it.pdf'])

figure;
h = PlotStatesShaded(Year, qRshort_bar_jp_idio); hold on;
h.Color = 'm';
p_jp = plot(Year, qRshort_bar_jp_idio(:,3), 'm:', 'LineWidth', 2);
plot(Year, Stir_jp-Infl_jp-Rshort_country_average, 'm:', 'LineWidth', 1);
axis([Year(1) Year(end) -3 6]);
%title('Japan', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa2-Model1_Rshort-countries_trend-idio_obs-average-jp.pdf'])

