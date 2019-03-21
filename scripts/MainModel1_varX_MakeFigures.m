%% MainModel1_varX_MakeFigures.m    Produce appendix figure A3: This gives
%                                   the Baseline model with alternative 
%                                   specifications for the Variance of innovations to the Trend


%% Plot preliminaries: load, sort, and get quantiles

close all;
addpath('Routines')

resultsPath = '../results/';
figpath   = '../figures/';
appenpath = [figpath 'appendix/'];


set(0,'defaultAxesFontName', 'Times');
set(0,'DefaultTextInterpreter', 'latex')
set(0, 'DefaultAxesFontSize',15)
set(0,'defaultAxesLineStyleOrder','-|--|:', 'defaultLineLineWidth',1.5)
setappdata(0, 'defaultAxesXTickFontSize', 1)
setappdata(0, 'defaultAxesYTickFontSize', 1)

Quant = [0.025 0.160 0.500 0.840  0.975];  % Quantiles

% Extract objects from baseline model (1/100)
temp = load([resultsPath 'OutputModel1.mat'], 'CommonTrends');
Year = load([resultsPath 'OutputModel1.mat'], 'Year'); 
Year = Year.Year;
CommonTrends_var100 = temp.CommonTrends;
M = size(CommonTrends_var100, 3);
Rshort_bar_var100   = squeeze(CommonTrends_var100(:, 1,:));
sRshort_bar_var100  = sort(Rshort_bar_var100,2);
qRshort_bar_var100  = sRshort_bar_var100(:,ceil(Quant*M));

% Extract objects from modified baseline model (1/50)
temp = load([resultsPath 'OutputModel1_var50.mat'], 'CommonTrends');
CommonTrends_var50 = temp.CommonTrends;
M = size(CommonTrends_var50, 3);
Rshort_bar_var50   = squeeze(CommonTrends_var50(:, 1,:));
sRshort_bar_var50  = sort(Rshort_bar_var50, 2);
qRshort_bar_var50  = sRshort_bar_var50(:, ceil(Quant*M));

% Extract objects from modified baseline model (1/25)
temp = load([resultsPath 'OutputModel1_var25.mat'], 'CommonTrends');
CommonTrends_var25 = temp.CommonTrends;
Rshort_bar_var25   = squeeze(CommonTrends_var25(:, 1,:));
sRshort_bar_var25  = sort(Rshort_bar_var25, 2);
qRshort_bar_var25  = sRshort_bar_var25(:, ceil(Quant*M));

% Extract objects from modified baseline model (1/10)
temp = load([resultsPath 'OutputModel1_var10.mat'], 'CommonTrends');
CommonTrends_var10 = temp.CommonTrends;
Rshort_bar_var10   = squeeze(CommonTrends_var10(:, 1,:));
sRshort_bar_var10  = sort(Rshort_bar_var10, 2);
qRshort_bar_var10  = sRshort_bar_var10(:, ceil(Quant*M));

% Extract objects from modified baseline model (1/5)
temp = load([resultsPath 'OutputModel1_var05.mat'], 'CommonTrends');
CommonTrends_var05 = temp.CommonTrends;
Rshort_bar_var05   = squeeze(CommonTrends_var05(:, 1,:));
sRshort_bar_var05  = sort(Rshort_bar_var05, 2);
qRshort_bar_var05  = sRshort_bar_var05(:, ceil(Quant*M));

% Extract objects from modified baseline model (1/2)
temp = load([resultsPath 'OutputModel1_var02.mat'], 'CommonTrends');
CommonTrends_var02 = temp.CommonTrends;
Rshort_bar_var02   = squeeze(CommonTrends_var02(:, 1,:));
sRshort_bar_var02  = sort(Rshort_bar_var02, 2);
qRshort_bar_var02  = sRshort_bar_var02(:, ceil(Quant*M));

% Extract objects from modified baseline model (1/1)
temp = load([resultsPath 'OutputModel1_var01.mat'], 'CommonTrends');
CommonTrends_var01 = temp.CommonTrends;
M                  = size(CommonTrends_var01, 3);
Rshort_bar_var01   = squeeze(CommonTrends_var01(:, 1,:));
sRshort_bar_var01  = sort(Rshort_bar_var01,2);
qRshort_bar_var01  = sRshort_bar_var01(:,ceil(Quant*M));


%% Produce figure
% A3: Alternative Priors ofr the Standard Deviation of Innovaitons to the
% Trend, Baseline Model

fSize = 15;  % Font size

figure
bands_new(Year, qRshort_bar_var100, qRshort_bar_var50); 
hold on;
bands_new(Year, qRshort_bar_var100); hold on;
axis([Year(1) Year(end) -12 12]);
% title('$\sigma_0 = 1/50$', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa3-Model1_var50_Rbar-MA.pdf'])


figure
bands_new(Year, qRshort_bar_var100, qRshort_bar_var25); 
hold on;
bands_new(Year, qRshort_bar_var100); hold on;
axis([Year(1) Year(end) -12 12]);
%title('$\sigma_0 = 1/25$', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa3-Model1_var25_Rbar-MA.pdf'])

figure
bands_new(Year, qRshort_bar_var100, qRshort_bar_var10); 
hold on;
bands_new(Year, qRshort_bar_var100); hold on;
axis([Year(1) Year(end) -12 12]);
%title('$\sigma_0 = 1/10$', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa3-Model1_var10_Rbar-MA.pdf'])

figure
bands_new(Year, qRshort_bar_var100, qRshort_bar_var05); 
hold on;
bands_new(Year, qRshort_bar_var100); hold on;
axis([Year(1) Year(end) -12 12]);
%title('$\sigma_0 = 1/5$', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa3-Model1_var05_Rbar-MA.pdf'])

figure
bands_new(Year, qRshort_bar_var100, qRshort_bar_var02); 
hold on;
bands_new(Year, qRshort_bar_var100); hold on;
axis([Year(1) Year(end) -12 12]);
%title('$\sigma_0 = 1/2$', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa3-Model1_var02_Rbar-MA.pdf'])

figure
bands_new(Year, qRshort_bar_var100, qRshort_bar_var01); 
hold on;
bands_new(Year, qRshort_bar_var100); hold on;
axis([Year(1) Year(end) -12 12]);
%title('$\sigma_0 = 1$', 'Interpreter', 'latex')
set(gca, 'FontSize', fSize)
printpdf(gcf, [appenpath 'figa3-Model1_var01_Rbar-MA.pdf'])

