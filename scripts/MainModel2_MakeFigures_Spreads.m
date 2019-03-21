% MainModel2_MakeFigures_Spreads.m    Creates figures in appendix B9.
%                                     Compares country convenience yields
%                                     to spreads. See figures A25, A26.

%% Preliminaries

%clear;
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

% ---------------------- trends ------------------------------------------

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

% Rshort_bar_us = repmat(transpose(squeeze(CC(1,1,:))), T, 1) .* Rshort_bar + Rshort_bar_us_idio;
% Rshort_bar_de = repmat(transpose(squeeze(CC(2,1,:))), T, 1) .* Rshort_bar + Rshort_bar_de_idio;
% Rshort_bar_uk = repmat(transpose(squeeze(CC(3,1,:))), T, 1) .* Rshort_bar + Rshort_bar_uk_idio;
% Rshort_bar_fr = repmat(transpose(squeeze(CC(4,1,:))), T, 1) .* Rshort_bar + Rshort_bar_fr_idio;
% Rshort_bar_ca = repmat(transpose(squeeze(CC(5,1,:))), T, 1) .* Rshort_bar + Rshort_bar_ca_idio;
% Rshort_bar_it = repmat(transpose(squeeze(CC(6,1,:))), T, 1) .* Rshort_bar + Rshort_bar_it_idio;
% Rshort_bar_jp = repmat(transpose(squeeze(CC(7,1,:))), T, 1) .* Rshort_bar + Rshort_bar_jp_idio;

Rshort_bar_us = Cy_bar - Rshort_bar_us_idio;
Rshort_bar_de = Cy_bar - Rshort_bar_de_idio;
Rshort_bar_uk = Cy_bar - Rshort_bar_uk_idio;
Rshort_bar_fr = Cy_bar - Rshort_bar_fr_idio;
Rshort_bar_ca = Cy_bar - Rshort_bar_ca_idio;
Rshort_bar_it = Cy_bar - Rshort_bar_it_idio;
Rshort_bar_jp = Cy_bar - Rshort_bar_jp_idio;

% ------------------------ sorted trends ----------------------------------

sRshort_bar   = sort(Rshort_bar,2);
sPi_bar       = sort(Pi_bar,2);
sTs_bar       = sort(Ts_bar,2);
sRlong_bar    = sort(Rlong_bar,2);

sRshort_bar_us = sort(Rshort_bar_us,2);
sRshort_bar_de = sort(Rshort_bar_de,2);
sRshort_bar_uk = sort(Rshort_bar_uk,2);
sRshort_bar_fr = sort(Rshort_bar_fr,2);
sRshort_bar_ca = sort(Rshort_bar_ca,2);
sRshort_bar_it = sort(Rshort_bar_it,2);
sRshort_bar_jp = sort(Rshort_bar_jp,2);

% ------------------- quantiles of the trends ----------------------------

qRshort_bar   = sRshort_bar(:,ceil(Quant*M));
qPi_bar       = sPi_bar(:,ceil(Quant*M));
qTs_bar       = sTs_bar(:,ceil(Quant*M));
qRlong_bar    = sRlong_bar(:,ceil(Quant*M));

qRshort_bar_us = sRshort_bar_us(:,ceil(Quant*M));
qRshort_bar_de = sRshort_bar_de(:,ceil(Quant*M));
qRshort_bar_uk = sRshort_bar_uk(:,ceil(Quant*M));
qRshort_bar_fr = sRshort_bar_fr(:,ceil(Quant*M));
qRshort_bar_ca = sRshort_bar_ca(:,ceil(Quant*M));
qRshort_bar_it = sRshort_bar_it(:,ceil(Quant*M));
qRshort_bar_jp = sRshort_bar_jp(:,ceil(Quant*M));

% ---------------------------- Spread ------------------------------------
load('../indata/corpspreadAFE.mat');   % corpspread
corpspread.us = [Year Baa_us - Ltir_us];
corpspread.de = corpspread.Germany;
corpspread.uk = corpspread.UK;
corpspread.fr = corpspread.France;
corpspread.ca = corpspread.Canada;
corpspread.it = corpspread.Italy;
corpspread.jp = corpspread.Japan;


%% Figure A25: Baa-Treasury Spread and Convenience Yield: US


PlotStatesShaded(datenum(Year, 1, 1), qRshort_bar_us);
hold on;
plot(datenum(corpspread.us(:,1), 1, 1), corpspread.us(:, 2), '-k')

xticks(datenum(1880, 1, 1):20*365.25:datenum(2016, 1, 1))
datetick('keepticks')
axis([datenum(Year(1), 1, 1) datenum(Year(end), 1, 1) -1 6])
box on

printpdf(gcf, [appenpath 'figa25-Model2_CY_trend-idio_obs-spread-us.pdf'])


%% Figure A26: Spreads
% Note: This section was written with MATLAB R2017b. Older versions require
% different double-axis functions.

fSize = 15;  % Font size
% ----------------------------- Germany ----------------------------------

figure
yyaxis left

PlotStatesShaded(datenum(Year,1,1), qRshort_bar_de); 
ylim([-1 5])

yyaxis right
plot(corpspread.de(:,1), corpspread.de(:,2), '-b')
box on

xticks(datenum(1990,1,1):5*365.25:datenum(2015,2,1))
datetick('keepticks')
axis([datenum('1990-1-1'), datenum('2016-1-1'), -1, 5])
%title('Germany')
hline = refline(0);
hline.Color = 'k';
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa26-Model2_CY_trend-idio_obs-spread-de_axis.pdf'])


% ----------------------------- France ----------------------------------

figure
yyaxis left
PlotStatesShaded(datenum(Year,1,1), qRshort_bar_fr); 
ylim([-1 5])

yyaxis right
plot(corpspread.fr(:,1), corpspread.fr(:,2), '-y')
box on

xticks(datenum(1990,1,1):5*365.25:datenum(2015,2,1))
datetick('keepticks')
axis([datenum('1990-1-1'), datenum('2016-1-1'), -1, 5])
%title('France')
hline = refline(0);
hline.Color = 'k';
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa26-Model2_CY_trend-idio_obs-spread-fr_axis.pdf'])


% ----------------------------- Italy ----------------------------------
figure
yyaxis left
PlotStatesShaded(datenum(Year,1,1), qRshort_bar_it); 
ylim([-1 5])

yyaxis right
plot(corpspread.it(:,1), corpspread.it(:,2), '-g')
box on

xticks(datenum(1990,1,1):5*365.25:datenum(2015,2,1))
datetick('keepticks')
axis([datenum('1990-1-1'), datenum('2016-1-1'), -1, 5])
%title('Italy')
hline = refline(0);
hline.Color = 'k';
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca, 'FontSize', fSize)

printpdf(gcf, [appenpath 'figa26-Model2_CY_trend-idio_obs-spread-it_axis.pdf'])

% ----------------------------- U.K. ----------------------------------
figure
yyaxis left
PlotStatesShaded(datenum(Year,1,1), qRshort_bar_uk); 
ylim([-1 5])

yyaxis right
plot(corpspread.uk(:,1), corpspread.uk(:,2) / 100, '-cy')  % Change units
box on

xticks(datenum(1990,1,1):5*365.25:datenum(2015,2,1))
datetick('keepticks')
axis([datenum('1990-1-1'), datenum('2016-1-1'), -1, 5])
%title('U.K.')
hline = refline(0);
hline.Color = 'k';
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca, 'FontSize', fSize)

printpdf(gcf, [appenpath 'figa26-Model2_CY_trend-idio_obs-spread-uk_axis.pdf'])

% ----------------------------- Canada ----------------------------------

figure
yyaxis left
PlotStatesShaded(datenum(Year,1,1), qRshort_bar_ca); 
ylim([-1 5])

yyaxis right
plot(corpspread.ca(:,1), corpspread.ca(:,2) / 100, '-r')  % Change units
box on

xticks(datenum(1990,1,1):5*365.25:datenum(2015,2,1))
datetick('keepticks')
axis([datenum('1990-1-1'), datenum('2016-1-1'), -1, 5])
%title('Canada')
hline = refline(0);
hline.Color = 'k';
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
set(gca, 'FontSize', fSize)

% Set figure dimensions
%set(gcf, 'Position', [100 100 1200 900])

printpdf(gcf, [appenpath 'figa26-Model2_CY_trend-idio_obs-spread-ca_axis.pdf'])

