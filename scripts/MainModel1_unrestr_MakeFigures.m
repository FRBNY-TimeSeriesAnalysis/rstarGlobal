% This file replicates the plots for the "Baseline" model with unrestricted 
% global real rate loadings. See figure A12, table A2.

%% Plot preliminaries: load, sort, and get quantiles


load('../results/OutputModel1_unrestr.mat') 

figpath   = '../figures/';
appenpath = [figpath 'appendix/'];


Quant = [0.025 0.160 0.500 0.840  0.975];
M = size(CommonTrends, 3);
addpath('Routines')

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


%% Figure A12: Trends in Global and U.S. Real Rates: 1870-2016, 
% Unrestricted Version of the Baseline Model

figure

PlotStatesShaded(Year, qRshort_bar)
hold on
plot(Year, qRshort_bar_us(:, 3), ...
    'k:', 'LineWidth', 1.5);  % r-bar US

axis([Year(1) Year(end) -3 6])
%title('$\bar{r}^w_t$ and $\bar{r}_{US, t}$', 'Interpreter', 'latex')
box on

printpdf(gcf, [appenpath 'figa12-Model1_unrestr_Rshortbar-us.pdf'])


