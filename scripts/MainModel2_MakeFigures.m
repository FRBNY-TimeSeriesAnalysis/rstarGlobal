% Model 2
% Generate figure 6, A13-A17

%% Preliminaries

clc;
close all;

load('../results/OutputModel2.mat')
figpath   = '../figures/';
appenpath = [figpath 'appendix/'];

set(0,'defaultAxesFontName', 'Times');
set(0,'DefaultTextInterpreter', 'latex')
set(0, 'DefaultAxesFontSize',15)
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth',1.5)
setappdata(0, 'defaultAxesXTickFontSize', 1)
setappdata(0, 'defaultAxesYTickFontSize', 1)
addpath('Routines');

Quant = [0.025 0.160 0.500 0.840  0.975];
M = size(CommonTrends, 3);

% trends

M_bar         = squeeze(CommonTrends(:, 1,:));
Pi_bar        = squeeze(CommonTrends(:, 2,:));
Ts_bar        = squeeze(CommonTrends(:, 3,:));
Cy_bar        = squeeze(CommonTrends(:, 4,:));
Rshort_bar    = M_bar - Cy_bar;
Rlong_bar     = Rshort_bar + Ts_bar;

Rshort_bar_us_idio = squeeze(CommonTrends(:, 5,:));
Rshort_bar_de_idio = squeeze(CommonTrends(:, 6,:));
Rshort_bar_uk_idio = squeeze(CommonTrends(:, 7,:));
Rshort_bar_fr_idio = squeeze(CommonTrends(:, 8,:));
Rshort_bar_ca_idio = squeeze(CommonTrends(:, 9,:));
Rshort_bar_it_idio = squeeze(CommonTrends(:,10,:));
Rshort_bar_jp_idio = squeeze(CommonTrends(:,11,:));

Pi_bar_us_idio = squeeze(CommonTrends(:,12,:));
Pi_bar_de_idio = squeeze(CommonTrends(:,13,:));
Pi_bar_uk_idio = squeeze(CommonTrends(:,14,:));
Pi_bar_fr_idio = squeeze(CommonTrends(:,15,:));
Pi_bar_ca_idio = squeeze(CommonTrends(:,16,:));
Pi_bar_it_idio = squeeze(CommonTrends(:,17,:));
Pi_bar_jp_idio = squeeze(CommonTrends(:,18,:));

Ts_bar_us_idio = squeeze(CommonTrends(:,19,:));
Ts_bar_de_idio = squeeze(CommonTrends(:,20,:));
Ts_bar_uk_idio = squeeze(CommonTrends(:,21,:));
Ts_bar_fr_idio = squeeze(CommonTrends(:,22,:));
Ts_bar_ca_idio = squeeze(CommonTrends(:,23,:));
Ts_bar_it_idio = squeeze(CommonTrends(:,24,:));
Ts_bar_jp_idio = squeeze(CommonTrends(:,25,:));

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

% sorted trends

sRshort_bar     = sort(Rshort_bar,2);
sPi_bar         = sort(Pi_bar,2);
sTs_bar         = sort(Ts_bar,2);
sCy_bar         = sort(Cy_bar,2);
sM_bar          = sort(M_bar,2);
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

% quantiles of the trends

qRshort_bar     = sRshort_bar(:,ceil(Quant*M));
qPi_bar         = sPi_bar(:,ceil(Quant*M));
qTs_bar         = sTs_bar(:,ceil(Quant*M));
qCy_bar         = sCy_bar(:,ceil(Quant*M));
qM_bar          = sM_bar(:,ceil(Quant*M));
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


%% Figure 6a: r-bar^w

f = figure;

h= bands(Year, qRshort_bar);
box on
axis([Year(1) Year(end) -3 6])
%title('$\overline{r}^w_t$', 'Interpreter', 'latex')

printpdf(gcf, [figpath 'fig6a-Model2_Rbar.pdf'])


%% Figure 6b: r-bar^w and -cy-bar^w_t

figure

tmax = find(Year==1870);

% -qCy_bar is normalized
h = bands(Year, qRshort_bar, -qCy_bar + (qCy_bar(tmax,3) + qRshort_bar(tmax,3)));
hold on
axis([Year(1) Year(end) -3 6])
%title('$\overline{r}^w_t$ and $-\overline{cy}_t^w$', 'Interpreter', 'latex')
box on
printpdf(gcf, [figpath 'fig6b-Model2_Rbar-Cybar.pdf'])


%% Figure 6c: r-bar^w_t and m-bar^w_t

figure

% -qCy_bar is normalized
h = bands(Year, qRshort_bar, qM_bar - (qM_bar(tmax, 3) - qRshort_bar(tmax, 3)));

axis([Year(1) Year(end) -3 6])
%title('$\overline{r}^w_t$ and $\overline{m}_t^w$', 'Interpreter', 'latex')
box on

printpdf(gcf, [figpath 'fig6c-Model2_Rbar-Mbar.pdf'])



%% -------------------------- APPENDIX FIGURES ---------------------------- %%



%% A13: Trends in Global and U.S. Real Rates


figure

PlotStatesShaded(Year, qRshort_bar)
hold on
plot(Year, qRshort_bar_us(:,3), ...
    'k:', 'LineWidth', 1.5);  % r-bar US

axis([Year(1) Year(end) -3 6])
%title('$\overline{r}^w_t$ and $\overline{r}_{US, t}$', 'Interpreter', 'latex')
box on

printpdf(gcf, [appenpath 'figa13-Model2_Rshortbar-us.pdf'])

%% A14a: Trends and Observables for Short-Term Real Rates, Convenience Yield Model


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
%title('$\overline{r}^w_t$ and $R_{i,t} - \pi_{i,t}$', 'Interpreter', 'latex')

printpdf(gcf, [appenpath 'figa14a-Model2_Rshortbar-obs.pdf'])

%% A14b: Trends and Observables for Short-Term Real Rates, Convenience Yield Model

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

%title('$\overline{r}^w_t$ and $\overline{r}_{i,t}$', 'Interpreter','latex')

legend boxoff;
box on;

printpdf(gcf, [appenpath 'figa14b-Model2_Rshortbar-countries.pdf'])

%% A15a: Trends and Observables for Inflation, Convenience Yield Model


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
%title('$\overline{\pi}^w_t$ and $\pi_{it}$')

box on

printpdf(gcf, [appenpath 'figa15a-Model2_pibar-obs.pdf'])

%% A15b: Trends and Observables for Inflation, Convenience Yield Model

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

%title('$\overline{\pi}^w_t$ and $\overline{\pi}_{i,t}$', 'Interpreter', 'latex')

legend([p_us p_de p_uk p_fr p_ca p_it p_jp],...
    {'us', 'de', 'uk', 'fr', 'ca', 'it', 'jp'}, 'Location', 'southoutside',...
    'Orientation', 'horizontal', 'Box', 'off')
box on
yticks(0:5:15)

printpdf(gcf, [appenpath 'figa15b-Model2_pibar-countries.pdf'])

%% A16a: Trends and Observables for Term Spreads, Convenience Yield Model

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

printpdf(gcf, [appenpath 'figa16a-Model2_Rbar-Tsbar.pdf'])

%% A16b: Trends and Observables for Term Spreads, Convenience Yield Model

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

printpdf(gcf, [appenpath 'figa16b-Model2_Rbar-tsbar-countries.pdf'])

%% A17: Country-Specfic Trends r_it and Observables, Convenience Yield Model
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
printpdf(gcf, [appenpath 'figa17-Model2_Rshort-countries_trend-idio_obs-average-common.pdf'])

f = figure;
h = PlotStatesShaded(Year, qRshort_bar_us_idio); hold on;
p_us = plot(Year, qRshort_bar_us_idio(:,3), 'k:', 'LineWidth', 2);
plot(Year, Stir_us-Infl_us-Rshort_country_average, 'k:', 'LineWidth', 1);
axis([Year(1) Year(end) -3 6]);
%title('U.S.', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa17-Model2_Rshort-countries_trend-idio_obs-average-us.pdf'])

figure;
h = PlotStatesShaded(Year, qRshort_bar_de_idio); hold on;
h.Color = 'b';
h.LineStyle = '--';
p_de = plot(Year, qRshort_bar_de_idio(:,3), 'b:', 'LineWidth', 2);
plot(Year, Stir_de-Infl_de-Rshort_country_average, 'b:', 'LineWidth', 1);
 axis([Year(1) Year(end) -3 6]);
%title('Germany', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa17-Model2_Rshort-countries_trend-idio_obs-average-de.pdf'])

figure;
h = PlotStatesShaded(Year, qRshort_bar_uk_idio); hold on;
h.Color = 'c';
p_uk = plot(Year, qRshort_bar_uk_idio(:,3), 'c:', 'LineWidth', 2);
plot(Year, Stir_uk-Infl_uk-Rshort_country_average, 'c:', 'LineWidth', 1);
 axis([Year(1) Year(end) -3 6]);
%title('U.K.', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa17-Model2_Rshort-countries_trend-idio_obs-average-uk.pdf'])

figure;
h = PlotStatesShaded(Year, qRshort_bar_fr_idio); hold on;
h.Color = 'y';
p_fr = plot(Year, qRshort_bar_fr_idio(:,3), 'y:', 'LineWidth', 2);
plot(Year, Stir_fr-Infl_fr-Rshort_country_average, 'y:', 'LineWidth', 1);
axis([Year(1) Year(end) -3 6]);
%title('France', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa17-Model2_Rshort-countries_trend-idio_obs-average-fr.pdf'])

figure;
h = PlotStatesShaded(Year, qRshort_bar_ca_idio); hold on;
h.Color = 'r';
p_ca = plot(Year, qRshort_bar_ca_idio(:,3), 'r:', 'LineWidth', 2);
plot(Year, Stir_ca-Infl_ca-Rshort_country_average, 'r:', 'LineWidth', 1);
axis([Year(1) Year(end) -3 6]);
%title('Canada', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa17-Model2_Rshort-countries_trend-idio_obs-average-ca.pdf'])

figure;
h = PlotStatesShaded(Year, qRshort_bar_it_idio); hold on;
h.Color = 'g';
p_it = plot(Year, qRshort_bar_it_idio(:,3), 'g:', 'LineWidth', 2);
plot(Year, Stir_it-Infl_it-Rshort_country_average, 'g:', 'LineWidth', 1);
axis([Year(1) Year(end) -3 6]);
%title('Italy', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa17-Model2_Rshort-countries_trend-idio_obs-average-it.pdf'])

figure;
h = PlotStatesShaded(Year, qRshort_bar_jp_idio); hold on;
h.Color = 'm';
p_jp = plot(Year, qRshort_bar_jp_idio(:,3), 'm:', 'LineWidth', 2);
plot(Year, Stir_jp-Infl_jp-Rshort_country_average, 'm:', 'LineWidth', 1);
axis([Year(1) Year(end) -3 6]);
%title('Japan', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa17-Model2_Rshort-countries_trend-idio_obs-average-jp.pdf'])

