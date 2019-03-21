%% MainModel1_ReR_MakeFigures.m    Make figures for Baseline model with VAR
%                                  with common trends applied to changes in 
%                                  the real exchange rate.
%
% Corresponds to figures A23 and A24 in the paper.

%% Preliminaries

clc;
close all;

load('../results/OutputReRonly.mat')

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

DQ_bar = squeeze(CommonTrends(:,1,:));

DQ_bar_de_idio = squeeze(CommonTrends(:,2,:));
DQ_bar_uk_idio = squeeze(CommonTrends(:,3,:));
DQ_bar_fr_idio = squeeze(CommonTrends(:,4,:));
DQ_bar_ca_idio = squeeze(CommonTrends(:,5,:));
DQ_bar_it_idio = squeeze(CommonTrends(:,6,:));
DQ_bar_jp_idio = squeeze(CommonTrends(:,7,:));

DQ_bar_de = repmat(transpose(squeeze(CC(1,1,:))), T, 1) .* DQ_bar + DQ_bar_de_idio;
DQ_bar_uk = repmat(transpose(squeeze(CC(2,1,:))), T, 1) .* DQ_bar + DQ_bar_uk_idio;
DQ_bar_fr = repmat(transpose(squeeze(CC(3,1,:))), T, 1) .* DQ_bar + DQ_bar_fr_idio;
DQ_bar_ca = repmat(transpose(squeeze(CC(4,1,:))), T, 1) .* DQ_bar + DQ_bar_ca_idio;
DQ_bar_it = repmat(transpose(squeeze(CC(5,1,:))), T, 1) .* DQ_bar + DQ_bar_it_idio;
DQ_bar_jp = repmat(transpose(squeeze(CC(6,1,:))), T, 1) .* DQ_bar + DQ_bar_jp_idio;

% sorted trends

sDQ_bar = sort(DQ_bar, 2);

sDQ_bar_de = sort(DQ_bar_de,2);
sDQ_bar_uk = sort(DQ_bar_uk,2);
sDQ_bar_fr = sort(DQ_bar_fr,2);
sDQ_bar_ca = sort(DQ_bar_ca,2);
sDQ_bar_it = sort(DQ_bar_it,2);
sDQ_bar_jp = sort(DQ_bar_jp,2);

sDQ_bar_de_idio = sort(DQ_bar_de_idio,2);
sDQ_bar_uk_idio = sort(DQ_bar_uk_idio,2);
sDQ_bar_fr_idio = sort(DQ_bar_fr_idio,2);
sDQ_bar_ca_idio = sort(DQ_bar_ca_idio,2);
sDQ_bar_it_idio = sort(DQ_bar_it_idio,2);
sDQ_bar_jp_idio = sort(DQ_bar_jp_idio,2);

% quantiles of the trends

qDQ_bar = sDQ_bar(:,ceil(Quant*M));

qDQ_bar_de = sDQ_bar_de(:,ceil(Quant*M));
qDQ_bar_uk = sDQ_bar_uk(:,ceil(Quant*M));
qDQ_bar_fr = sDQ_bar_fr(:,ceil(Quant*M));
qDQ_bar_ca = sDQ_bar_ca(:,ceil(Quant*M));
qDQ_bar_it = sDQ_bar_it(:,ceil(Quant*M));
qDQ_bar_jp = sDQ_bar_jp(:,ceil(Quant*M));

qDQ_bar_de_idio = sDQ_bar_de_idio(:,ceil(Quant*M));
qDQ_bar_uk_idio = sDQ_bar_uk_idio(:,ceil(Quant*M));
qDQ_bar_fr_idio = sDQ_bar_fr_idio(:,ceil(Quant*M));
qDQ_bar_ca_idio = sDQ_bar_ca_idio(:,ceil(Quant*M));
qDQ_bar_it_idio = sDQ_bar_it_idio(:,ceil(Quant*M));
qDQ_bar_jp_idio = sDQ_bar_jp_idio(:,ceil(Quant*M));


YRer = [dRer_us, dRer_de, dRer_uk, dRer_fr, dRer_ca, dRer_it, dRer_jp];
YRer(abs(YRer)>30)=NaN;
DQ_country_average = mean(YRer, 2, 'omitnan');


%% Figure A23: Common Trend and Observables for Change in the Real Exchange Rate

figure;

PlotStatesShaded(Year, qDQ_bar);
axis([Year(1) Year(end) -10 10]);

hold on
LW = 1;
p_de = plot(Year, dRer_de, ':b',  'LineWidth', LW);
p_uk = plot(Year, dRer_uk, ':cy', 'LineWidth', LW);
p_fr = plot(Year, dRer_fr, ':y',  'LineWidth', LW);
p_ca = plot(Year, dRer_ca, ':r',  'LineWidth', LW);
p_it = plot(Year, dRer_it, ':g',  'LineWidth', LW);
p_jp = plot(Year, dRer_jp, ':m',  'LineWidth', LW);

% title('$\overline{\Delta q_t^w}$ and $\Delta q_{i,t}$', 'Interpreter', 'latex')

printpdf(gcf, [appenpath 'figa23-ModelReR_DQbar-obs.pdf'])


%% Figure A24: Common Trend and Observables for Changes in the Real Exchange Rate

fSize = 15;

f = figure;
h = PlotStatesShaded(Year, qDQ_bar);
hold on; box on;
axis([Year(1) Year(end) -10 10]);
%title('World')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa24-ModelReR_DQ-countries_trend-idio_obs-average-common.pdf'])

figure
h    = PlotStatesShaded(Year, qDQ_bar_de_idio); hold on;
h.Color = 'b';
plot(Year, dRer_de-DQ_country_average, 'b:', 'LineWidth', 1);
axis([Year(1) Year(end) -10 10]);
%title('Germany')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa24-ModelReR_DQ-countries_trend-idio_obs-average-de.pdf'])

figure
h = PlotStatesShaded(Year, qDQ_bar_uk_idio); hold on;
h.Color = 'c';
plot(Year, dRer_uk-DQ_country_average, 'c:', 'LineWidth', 1);
axis([Year(1) Year(end) -10 10]);
%title('U.K.')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa24-ModelReR_DQ-countries_trend-idio_obs-average-uk.pdf'])

figure
h = PlotStatesShaded(Year, qDQ_bar_fr_idio); hold on;
h.Color = 'y';
plot(Year, dRer_fr-DQ_country_average, 'y:', 'LineWidth', 1);
axis([Year(1) Year(end) -10 10]);
%title('France')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa24-ModelReR_DQ-countries_trend-idio_obs-average-fr.pdf'])

figure
h = PlotStatesShaded(Year, qDQ_bar_ca_idio); hold on;
h.Color = 'r';
plot(Year, dRer_ca - DQ_country_average, 'r:', 'LineWidth', 1);
axis([Year(1) Year(end) -10 10]);
%title('Canada')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa24-ModelReR_DQ-countries_trend-idio_obs-average-ca.pdf'])

figure
h = PlotStatesShaded(Year, qDQ_bar_it_idio); hold on;
h.Color = 'g';
plot(Year, dRer_it - DQ_country_average, 'g:', 'LineWidth', 1);
axis([Year(1) Year(end) -10 10]);
%title('Italy')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa24-ModelReR_DQ-countries_trend-idio_obs-average-it.pdf'])

figure
h = PlotStatesShaded(Year, qDQ_bar_jp_idio); hold on;
h.Color = 'm';
plot(Year, dRer_jp - DQ_country_average, 'm:', 'LineWidth', 1);
axis([Year(1) Year(end) -10 10]);
%title('Japan')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa24-ModelReR_DQ-countries_trend-idio_obs-average-jp.pdf'])

