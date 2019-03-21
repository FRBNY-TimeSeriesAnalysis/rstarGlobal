%% MainModel1.m    Estimates the "Baseline Model" from Global Trends in 
%                  Interest Rates with 50 degrees of freedom. Results are saved as
%                  "OutputModel1_df50.mat."

%% Initial setup

clear;

filename = '../results/OutputModel1_df50.mat';  % Output filename

addpath Routines

[DATA,TEXT] = xlsread('../indata/DataInflShortLongUpdated.xlsx');
Year = DATA(:,1);

Ndraws = 20000;  % Number of MCMC draws
p = 1;             % Number of lags in the VAR for the cycle;



Mnem = TEXT(2:end);
%1)    'cpi_usa'
%2)    'stir_usa'
%3)    'ltir_usa'
%4)    'cpi_deu'
%5)    'stir_deu'
%6)    'ltir_deu'
%7)    'cpi_uk'
%8)    'stir_uk'
%9)    'ltir_uk'
%10)   'cpi_fr'
%11)   'stir_fr'
%12)   'ltir_fr'
%13)   'cpi_ca'
%14)   'stir_ca'
%15)   'ltir_ca'
%16)   'cpi_it'
%17)   'stir_it'
%18)   'ltir_it'
%19)   'cpi_jp'
%20)   'stir_jp'
%21)   'ltir_jp'


Country = {'US','DE','UK','FR','CA','IT','JP'};
X = DATA(:,2:end-1);

Price_us = X(:,1);
Price_de = X(:,4);
Price_uk = X(:,7);
Price_fr = X(:,10);
Price_ca = X(:,13);
Price_it = X(:,16);
Price_jp = X(:,19);


Stir_us = X(:,2);
Stir_de = X(:,5);
Stir_uk = X(:,8);
Stir_fr = X(:,11);
Stir_ca = X(:,14);
Stir_it = X(:,17);
Stir_jp = X(:,20);


Ltir_us = X(:,3);
Ltir_de = X(:,6);
Ltir_uk = X(:,9);
Ltir_fr = X(:,12);
Ltir_ca = X(:,15);
Ltir_it = X(:,18);
Ltir_jp = X(:,21);



%Inflation rate: national currency
Infl_us   = [NaN;(Price_us(2:end)./Price_us(1:end-1)-1)*100];
Infl_de   = [NaN;(Price_de(2:end)./Price_de(1:end-1)-1)*100];
Infl_uk   = [NaN;(Price_uk(2:end)./Price_uk(1:end-1)-1)*100];
Infl_fr   = [NaN;(Price_fr(2:end)./Price_fr(1:end-1)-1)*100];
Infl_ca   = [NaN;(Price_ca(2:end)./Price_ca(1:end-1)-1)*100];
Infl_it   = [NaN;(Price_it(2:end)./Price_it(1:end-1)-1)*100];
Infl_jp   = [NaN;(Price_jp(2:end)./Price_jp(1:end-1)-1)*100];


Y = [
    Stir_us...
    Stir_de...
    Stir_uk...
    Stir_fr...
    Stir_ca...
    Stir_it...
    Stir_jp...
    Infl_us...
    Infl_de...
    Infl_uk...
    Infl_fr...
    Infl_ca...
    Infl_it...
    Infl_jp...
    Ltir_us...
    Ltir_de...
    Ltir_uk...
    Ltir_fr...
    Ltir_ca...
    Ltir_it...
    Ltir_jp
    ];



Mnem = {
    'Stir_us'...
    'Stir_de'...
    'Stir_uk'...
    'Stir_fr'...
    'Stir_ca'...
    'Stir_it'...
    'Stir_jp'...
    'Infl_us'...
    'Infl_de'...
    'Infl_uk'...
    'Infl_fr'...
    'Infl_ca'...
    'Infl_it'...
    'Infl_jp'...
    'Ltir_us'...
    'Ltir_de'...
    'Ltir_uk'...
    'Ltir_fr'...
    'Ltir_ca'...
    'Ltir_it'...
    'Ltir_jp'
    };




Y(abs(Y)>30)=NaN;   % Eliminate noisy values


[T,n] = size(Y);

T0pre = find(Year==1870);
T1pre = find(Year==1899);
disp(['Avg. and std in the presample: 1954-1959'])
disp([(1:n)' nanmean(Y(T0pre:T1pre,:))' nanstd(Y(T0pre:T1pre,:))'])

disp('mean Stir')
disp(nanmean(nanmean(Y(T0pre:T1pre,1:7))))

disp('mean Infl')
disp(nanmean(nanmean(Y(T0pre:T1pre,8:14))))

disp('mean Ltir')
disp(nanmean(nanmean(Y(T0pre:T1pre,15:end))))


disp('std Stir')
disp(nanmean(nanstd(Y(T0pre:T1pre,1:7))))

disp('std Infl')
disp(nanmean(nanstd(Y(T0pre:T1pre,7:14))))

disp('std Ltir')
disp(nanmean(nanstd(Y(T0pre:T1pre,15:end))))



%% Setup model, initial conditions

T0 = find(Year==1870);
T1 = find(Year==2016);


Y = Y(T0:T1,:);
Year = Year(T0:T1);
y=Y;
[T,n] = size(y);

Ctr =[
    1       1      0%     Stir_us...
    1       1      0%     Stir_de...
    1       1      0%     Stir_uk...
    1       1      0%     Stir_fr...
    1       1      0%     Stir_ca...
    1       1      0%     Stir_it...
    1       1      0%     Stir_jp...
    0       1      0%     Infl_us...
    0       1      0%     Infl_de...
    0       1      0%     Infl_uk...
    0       1      0%     Infl_fr...
    0       1      0%     Infl_ca...
    0       1      0%     Infl_it...
    0       1      0%     Infl_jp...
    1       1      1%     Ltir_us...
    1       1      1%     Ltir_de...
    1       1      1%     Ltir_uk...
    1       1      1%     Ltir_fr...
    1       1      1%     Ltir_ca...
    1       1      1%     Ltir_it...
    1       1      1%     Ltir_jp...
    ];

%Adding country specific trends to real rates
Cadd1                   =    zeros(n,7);
Cadd1(1:7,1:7)          =    eye(7);
Cadd1(15:21,1:7)        =    eye(7);

%Adding the country specific trends in inflation rates
Cadd2              =    zeros(n,7);
Cadd2(1:7,1:7)     =    eye(7);
Cadd2(8:14,1:7)    =    eye(7);
Cadd2(15:21,1:7)   =    eye(7);

%Adding country specific trends to term spread
Cadd3                 =    zeros(n,7);
Cadd3(15:21,1:7)      =    eye(7);



Ctr           = [Ctr Cadd1 Cadd2 Cadd3];
Ccyc          = zeros(n,n*p); 
Ccyc(1:n,1:n) = eye(n);  % Country/series-specific cyclic component
C             = [Ctr Ccyc];


r = size(Ctr,2);    % 3 + nCountries *3, number of non-cyclic components

b0          = zeros(n*p,n);
b0(1:n,1:n) = eye(n)*0;

df0tr = 50;  % Degrees of freedom

%            rs_wrd  pi_wrd        ts_wrd   rs_idio       pi_idio            ts_idio
SC0tr =    ([   1      sqrt(2)       1     1*ones(1,7)    sqrt(2)*ones(1,7)  1*ones(1,7)   ]).^2/100;  % Initial conditions Q
S0tr  =     [   .5     2             1      zeros(1,7)    zeros(1,7)         zeros(1,7)    ]';  % Initial conditions states
P0tr = diag([   1      2             1     1*ones(1,7)/2  2*ones(1,7)/2      1*ones(1,7)/2 ].^2);

%                  stir          infl        ltir
Psi =       (2*[  ones(1,7)   2*ones(1,7)  ones(1,7) ]).^2;

S0cyc = zeros(n*p,1);  % Initialize cyclic component

Atr           = eye(r);
Qtr           = diag(SC0tr);

% Initialize cyclic component
My             = ones(T,1)*nanmean(y); 
yint           = y; 
yint(isnan(y)) = My(isnan(y));
[Trend,Ycyc]   = hpfilter(yint,1000);
[beta, sigma]  = BVAR(Ycyc, p, b0, Psi, .2, 0);  

Acyc                  = zeros(n*p); 
Acyc(n+1:end,1:end-n) = eye(n*(p-1));
Acyc(1:n,:)           = beta';

Qcyc          = zeros(n*p);
Qcyc(1:n,1:n) = (sigma+sigma')/2;  % Symmetric
P0cyc         = dlyap(Acyc,Qcyc);



% Initialize transition matrix
A                   = zeros(r+n*p);
A(1:r,1:r)          = Atr;
A(r+1:end,r+1:end)  = Acyc;

% Initialize variance-covariance matrix of transition equation
Q                  = zeros(r+n*p);
Q(1:r,1:r)         = Qtr;
Q(r+1:end,r+1:end) = Qcyc;


R = eye(n)*1e-12;


%Starting conditions for the Kalman recursion
S0                  = [S0tr;S0cyc];
P0                  = zeros(r+n*p);
P0(1:r,1:r)         = P0tr;
P0(r+1:end,r+1:end) = P0cyc;



tic

% Store MCMC
P_acc  = ones(1,Ndraws)*NaN;
States = ones(T,r+n*p,Ndraws)*NaN;
Trends = ones(T,n,Ndraws)*NaN;
LogLik = ones(1,Ndraws)*NaN;
SS0    = ones(r,Ndraws)*NaN;
Theta  = ones(1,Ndraws);
AA     = ones(r+n*p,r+n*p,Ndraws)*NaN;
QQ     = ones(r+n*p,r+n*p,Ndraws)*NaN;
CC     = ones(n,r+n*p,Ndraws)*NaN;
RR     = ones(n,n,Ndraws)*NaN;

P_acc_var = ones(1,Ndraws)*NaN;
P_acc     = ones(1,Ndraws)*NaN;



% Set prior to paper's lambda (coefficient to inf_wrld)
mean_theta = 1;
std_theta = .5; 


C(1:21,2)        = mean_theta;
C([1:7 15:21],1) = mean_theta;


lambda       = .2;  % Minnesota prior

logML = -inf;  % Initialize log maximum likelihood

%% Begin estimation
% For reference, see the appendix.

for jm = 1:Ndraws  % MCMC draws
    
    % ------------------------- Block 1a ----------------------------------
    
    kf = KF(y, C, R, A, Q, S0, P0);
    loglik = kf.LogLik;
    
    
    theta_old = [];
    theta_old = [theta_old; C(1:7,2)] ;
    
    
    theta_new = theta_old + randn(size(theta_old))*std_theta/5; 
    
    C_new          = C;
    C_new(1:7  ,2) = theta_new(1:7);       % Place new inflation coefficients
    C_new(8:14 ,2) = theta_new(1:7);
    C_new(15:21,2) = theta_new(1:7);
    
    
    kf_new     = KF(y, C_new, R, A, Q, S0, P0);  % Likelihood of new parameters
    loglik_new = kf_new.LogLik;
    
    log_rat    = (loglik_new + sum(log(normpdf(theta_new,mean_theta,std_theta^2)))) ...
               - (loglik     + sum(log(normpdf(theta_old,mean_theta,std_theta^2))));
    p_acc = min(exp(log_rat),1);
    
    
    if rand<=p_acc  % Accept
        C      = C_new;
        loglik = loglik_new;
        kf     = kf_new;
    end;
    
    % ----------------------------- Block 2a ------------------------------
    
    kc   = KC(kf);  % Smoother
    Ytr  = [kc.S0(1:r)'; kc.S(:,1:r)]; % Trend components 
    Ycyc = kc.S(:,r+1:r+n);            % Cyclical component
    
    
    % ----------------------------- Block 2b -------------------------------
    SCtr                = CovarianceDraw(diff(Ytr),df0tr,diag(SC0tr));
    Q(1:r,1:r)          = SCtr;

    
    
    for jp=1:p  % Organize cycle components by including the initial values
        Ycyc = [kc.S0(r+(jp-1)*n+1:r+n*jp)'; Ycyc];
    end
    
    % Estimate cycle component
    [beta,sigma]           = BVAR(Ycyc, p, b0, Psi, lambda, 1);
    Acyc_new               = Acyc;
    Acyc_new(1:n,:)        = beta';
    
    if max(abs(eigs(Acyc_new))) < 1  % Convergence criteria
        
        Qcyc_new               = Qcyc;
        Qcyc_new(1:n,1:n)      = (sigma+sigma')/2;  % Symmetric
        P0cyc_new              = dlyap(Acyc_new, Qcyc_new);
        P0cyc_new              = (P0cyc_new + P0cyc_new')/2;
        Y0cyc                  = kc.S0(r+1:end);
        
        rat = mvnpdf(Y0cyc, Y0cyc*0, P0cyc_new) / mvnpdf(Y0cyc, Y0cyc*0, P0cyc);
        
        p_acc_var = min(rat,1);
        
        
        if rand<=p_acc_var  % Accept
            Acyc  = Acyc_new;
            Qcyc  = Qcyc_new;
            P0cyc = P0cyc_new;
            A(r+1:end,r+1:end)  = Acyc;
            Q(r+1:end,r+1:end)  = Qcyc;
            P0(r+1:end,r+1:end) = P0cyc;
        end;
        
    else
        p_acc_var =0;
    end
    
    % Store MCMC
    States(:,:,jm) = kc.S;
    Trends(:,:,jm) = kc.S(:,1:r)*C(:,1:r)';
    LogLik(jm) = loglik;
    LogML(jm)  = logML;
    SS0(:,jm)  = S0(1:r);
    AA(:,:,jm) = A;
    QQ(:,:,jm) = Q;
    CC(:,:,jm) = C;
    RR(:,:,jm) = R;
    Lambda(jm) = lambda;
    
    P_acc(jm)    = p_acc;
    P_acc_var(jm)= p_acc_var;
    
    if mod(jm,10)==0  % Print to command window
        if jm>1
            if jm <=1000
                disp([num2str(jm),'th draw of ',num2str(Ndraws),'; Elapsed time: ',num2str(toc),' seconds'])
                disp(['Acceptance rate so far: ',num2str(mean(P_acc(1:jm)))])
                disp(['Acceptance rate so far: ',num2str(mean(P_acc_var(1:jm)))])
            elseif jm <10000
                disp([num2str(jm),'th draw of ',num2str(Ndraws),'; Elapsed time: ',num2str(toc),' seconds'])
                disp(['Acceptance rate of the last 1k draws: ',num2str(mean(P_acc(jm-1000+1:jm)))])
                disp(['Acceptance rate of the last 1k draws: ',num2str(mean(P_acc_var(jm-1000+1:jm)))])
            else
                disp([num2str(jm),'th draw of ',num2str(Ndraws),'; Elapsed time: ',num2str(toc),' seconds'])
                disp(['Acceptance rate of the last 10k draws: ',num2str(mean(P_acc(jm-5000+1:jm)))])
                disp(['Acceptance rate of the last 10k draws: ',num2str(mean(P_acc_var(jm-5000+1:jm)))])
            end;
        end;
        
    end;
end

Ndraws = length(Lambda)-1;

skip = 1;
Discard = floor(Ndraws/2);  % Burn-in index

% Remove burn-in draws
States = States(:,:,Discard+1:Ndraws);
Trends = Trends(:,:,Discard+1:Ndraws);
AA     = AA(:,:,Discard+1:Ndraws);
QQ     = QQ(:,:,Discard+1:Ndraws);
CC     = CC(:,:,Discard+1:Ndraws);
RR     = RR(:,:,Discard+1:Ndraws);
LogLik = LogLik(:,Discard+1:Ndraws);
LogML  = LogML(:,Discard+1:Ndraws);
SS0    = SS0(:,Discard+1:Ndraws);


CommonTrends = States(:,1:r,:);
Cycles       = States(:,r+1:r+n,:);


% Compute and store quantiles
Quant = [.025 .16 .5 .84  .975];
sCommonTrends=sort(CommonTrends,3);
sCycles=sort(Cycles,3);
sTrends=sort(Trends,3);

M = size(sCycles,3);
qCommonTrends = sCommonTrends(:,:,floor(Quant*M));
qCycles = sCycles(:,:,floor(Quant*M));
qTrends = sTrends(:,:,floor(Quant*M));


%% Save results

save(filename, '-v7.3')
