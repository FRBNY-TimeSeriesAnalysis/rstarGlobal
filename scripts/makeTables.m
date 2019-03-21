%% makeTables.m    Created tables 1, A1, A2
% Note: table 1 values come directly from table A1

addpath('Routines')

%% Table A1a: Change in r-bar^w_t and Its Components (baseline), additional details

model      = 'Model1';
model_name = 'Baseline';

load(['../results/Output' model '.mat']);


Quant = [0.025 0.050 0.160 0.500 0.840 0.950 0.975];
M     = size(CommonTrends, 3);
iQ    = ceil(Quant*M);

% trends
Rshort_bar    = squeeze(CommonTrends(:, 1,:));

x = {}; %#ok<*AGROW>

t_start = find(Year == 1980);
t_end   = find(Year == 2016);
dR = Rshort_bar(t_end,:) - Rshort_bar(t_start,:);
p = sum(dR < 0) / M;
temp = sort(dR);
z = temp(iQ);
x = [x {['\makecell{' sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $',...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];

t_start = find(Year == 1990);
t_end   = find(Year == 2016);
dR = Rshort_bar(t_end,:) - Rshort_bar(t_start,:);
p = sum(dR < 0) / M;
temp = sort(dR);
z = temp(iQ);
x = [x {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $',...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];

t_start = find(Year == 1997);
t_end   = find(Year == 2016);
dR = Rshort_bar(t_end,:) - Rshort_bar(t_start,:);
p = sum(dR < 0) / M;
temp = sort(dR);
z = temp(iQ);
x = [x {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $',...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];

row_labels = {'$ \overline{r}^{w}_{t} $'};
style = 'c|c|c|c';
header = {'', '1980-2016', '1990-2016', '1997-2016'};

fid = fopen(['../tables/' strrep(model_name, ' ', '_') '.tex'], 'w');
WriteTeXTable(fid, header, style, [row_labels x], [model_name '\\ \\']);
fclose(fid);


%% Table A1b: Change in r-bar^w_t and Its Components (Convenience Yield), additional details

model = 'Model2';
model_name = 'Convenience yield';

load(['../results/Output' model '.mat']);

Quant = [0.025 0.050 0.160 0.500 0.840 0.950 0.975];
M = size(CommonTrends, 3);
iQ = ceil(Quant*M);

% trends

M_bar         = squeeze(CommonTrends(:, 1,:));
Pi_bar        = squeeze(CommonTrends(:, 2,:));
Ts_bar        = squeeze(CommonTrends(:, 3,:));
Cy_bar        = squeeze(CommonTrends(:, 4,:));
Rshort_bar    = M_bar - Cy_bar;

x = {}; %#ok<*AGROW>

y = {};
t_start = find(Year == 1980);
t_end   = find(Year == 2016);
dR = Rshort_bar(t_end,:) - Rshort_bar(t_start,:);
p = sum(dR < 0) / M;
temp = sort(dR);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $',...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dM = M_bar(t_end,:) - M_bar(t_start,:);
p = sum(dM < 0) / M;
temp = sort(dM);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $', ...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dCy = Cy_bar(t_end,:) - Cy_bar(t_start,:);
p = sum(-dCy < 0) / M;
temp = sort(-dCy);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $', ...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
x = [x y];

y = {}; 
t_start = find(Year == 1990);
t_end   = find(Year == 2016);
dR = Rshort_bar(t_end,:) - Rshort_bar(t_start,:);
p = sum(dR < 0) / M;
temp = sort(dR);
z = temp(iQ);
y = [y; {['\makecell{'...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $',...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dM = M_bar(t_end,:) - M_bar(t_start,:);
p = sum(dM < 0) / M;
temp = sort(dM);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $', ...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dCy = Cy_bar(t_end,:) - Cy_bar(t_start,:);
p = sum(-dCy < 0) / M;
temp = sort(-dCy);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $', ...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
x = [x y];

y = {};
t_start = find(Year == 1997);
t_end   = find(Year == 2016);
dR = Rshort_bar(t_end,:) - Rshort_bar(t_start,:);
p = sum(dR < 0) / M;
temp = sort(dR);
z = temp(iQ);
y = [y; {['\makecell{'...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $', ...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dM = M_bar(t_end,:) - M_bar(t_start,:);
p = sum(dM < 0) / M;
temp = sort(dM);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $',...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dCy = Cy_bar(t_end,:) - Cy_bar(t_start,:);
p = sum(-dCy < 0) / M;
temp = sort(-dCy);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $', ...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
x = [x y];

row_labels = {'$ \overline{r}^{w}_{t} $'; '$ \overline{m}^{w}_{t} $'; '$ -\overline{cy}^{w}_{t} $'};
style = 'c|c|c|c';
header = {'', '1980-2016', '1990-2016', '1997-2016'};

fid = fopen(['../tables/' strrep(model_name, ' ', '_') '.tex'], 'w');
WriteTeXTable(fid, header, style, [row_labels x], [model_name '\\ \\']);
fclose(fid);


%% Table A1c: Change in r-bar^w_t and Its Components (Consumption)


model = 'Model3';
model_name = 'Consumption';

load(['../results/Output' model '.mat']);

Quant = [0.025 0.050 0.160 0.500 0.840 0.950 0.975];
M = size(CommonTrends, 3);
iQ = ceil(Quant*M);

% trends

G_bar         = squeeze(CommonTrends(:, 1,:));
Pi_bar        = squeeze(CommonTrends(:, 2,:));
Ts_bar        = squeeze(CommonTrends(:, 3,:));
Cy_bar        = squeeze(CommonTrends(:, 4,:));
Beta_bar      = squeeze(CommonTrends(:, 5,:));
Gamma_bar     = squeeze(CommonTrends(:, 6,:));
M_bar         = G_bar + Beta_bar;
Rshort_bar    = M_bar - Cy_bar;
Rlong_bar     = Rshort_bar + Ts_bar;
DC_bar        = G_bar + Gamma_bar;

x = {}; %#ok<*AGROW>

y = {};
t_start = find(Year == 1980);
t_end   = find(Year == 2016);
dR = Rshort_bar(t_end,:) - Rshort_bar(t_start,:);
p = sum(dR < 0) / M;
temp = sort(dR);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $', ...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dM = M_bar(t_end,:) - M_bar(t_start,:);
p = sum(dM < 0) / M;
temp = sort(dM);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $',...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dG = G_bar(t_end,:) - G_bar(t_start,:);
p = sum(dG < 0) / M;
temp = sort(dG);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $',...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dBeta = Beta_bar(t_end,:) - Beta_bar(t_start,:);
p = sum(dBeta < 0) / M;
temp = sort(dBeta);
z = temp(iQ);
y = [y; {['\makecell{'...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $', ...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dCy = Cy_bar(t_end,:) - Cy_bar(t_start,:);
p = sum(-dCy < 0) / M;
temp = sort(-dCy);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $',...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dDC = DC_bar(t_end,:) - DC_bar(t_start,:);
p = sum(dDC < 0) / M;
temp = sort(dDC);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $', ...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
x = [x y];

y = {};
t_start = find(Year == 1990);
t_end   = find(Year == 2016);
dR = Rshort_bar(t_end,:) - Rshort_bar(t_start,:);
p = sum(dR < 0) / M;
temp = sort(dR);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $',...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dM = M_bar(t_end,:) - M_bar(t_start,:);
p = sum(dM < 0) / M;
temp = sort(dM);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $', ...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dG = G_bar(t_end,:) - G_bar(t_start,:);
p = sum(dG < 0) / M;
temp = sort(dG);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $',...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dBeta = Beta_bar(t_end,:) - Beta_bar(t_start,:);
p = sum(dBeta < 0) / M;
temp = sort(dBeta);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $',...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dCy = Cy_bar(t_end,:) - Cy_bar(t_start,:);
p = sum(-dCy < 0) / M;
temp = sort(-dCy);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $', ...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dDC = DC_bar(t_end,:) - DC_bar(t_start,:);
p = sum(dDC < 0) / M;
temp = sort(dDC);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $', ...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
x = [x y];

y = {};
t_start = find(Year == 1997);
t_end   = find(Year == 2016);
dR = Rshort_bar(t_end,:) - Rshort_bar(t_start,:);
p = sum(dR < 0) / M;
temp = sort(dR);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $', ...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dM = M_bar(t_end,:) - M_bar(t_start,:);
p = sum(dM < 0) / M;
temp = sort(dM);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $',...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dG = G_bar(t_end,:) - G_bar(t_start,:);
p = sum(dG < 0) / M;
temp = sort(dG);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $', ...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dBeta = Beta_bar(t_end,:) - Beta_bar(t_start,:);
p = sum(dBeta < 0) / M;
temp = sort(dBeta);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $',...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dCy = Cy_bar(t_end,:) - Cy_bar(t_start,:);
p = sum(-dCy < 0) / M;
temp = sort(-dCy);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $', ...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
dDC = DC_bar(t_end,:) - DC_bar(t_start,:);
p = sum(dDC < 0) / M;
temp = sort(dDC);
z = temp(iQ);
y = [y; {['\makecell{' ...
    sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ (%0.2f, %0.2f) $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $ \\\\ $ \\{%0.2f\\} $', ...
    z(4), z(3), z(5), z(2), z(6), z(1), z(7), p) '}']}];
x = [x y];

row_labels = {'$ \overline{r}^{w}_{t} $'; '$ \overline{m}^{w}_{t} $'; ...
    '$ \overline{g}^{w}_{t} $'; '$ \overline{\beta}^{w}_{t} $';...
    '$ -\overline{cy}^{w}_{t} $'; '$ \overline{\Delta} c^{w}_{t} $'};
style = 'c|c|c|c';
header = {'', '1980-2016', '1990-2016', '1997-2016'};

fid = fopen(['../tables/' strrep(model_name, ' ', '_') '.tex'], 'w');
WriteTeXTable(fid, header, style, [row_labels x], [model_name '\\ \\']);
fclose(fid);



%% Table A2: Estimates of lambda^r_i for the Unrestricted Version of the Baseline Model

model = 'OutputModel1_unrestr';

Output.(model) = load(['../results/OutputModel1_unrestr.mat']);
C              = Output.(model).CC; % loading matrix `C` with estimated `lambda` for each draw `j`

% Extract only the common factors

p_quantiles = normcdf([-2 -1 0 1 2]); % quantile probabilities

x = {}; %#ok<*AGROW>
for i = 1:7 % loop across countries
    z = quantile(squeeze(C(i,1,:)), p_quantiles);
    x = [x; {['\makecell{' sprintf('$ %0.2f $ \\\\ $ [%0.2f, %0.2f] $ \\\\ $ \\langle %0.2f, %0.2f \\rangle $', z(3), z(2), z(4), z(1), z(5)) '}']}];
end
row_labels = {'us'; 'de'; 'uk'; 'fr'; 'ca'; 'it'; 'jp'};
style = '|c|c|';
header = {'', '$ \lambda^{\overline{r}_{t}}_{i} $'};

fprintf('\n Table: Estimates of lambda^r_i for the unrestricted versoin of the baseline model\n')
disp('   Country              Median       68% coverage        95% coverage')
disp([row_labels x])

fid = fopen('../tables/MainModel1_unrestr_lambda.tex', 'w');
WriteTeXTable(fid, header, style, [row_labels x]);
fclose(fid);

