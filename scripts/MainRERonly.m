%% MainReRonly.m    Estimates model with only ReR

%% Initial setup


clear;

filename = '../results/OutputReRonly.mat';  % Output filename

addpath Routines

[DATA,TEXT] = xlsread('../indata/DataInflShortRee.xlsx');

Year = DATA(:,1);

Ndraws  =  20000;
p = 1; % Number of lags in the VAR for the cycle;


    
Mnem = TEXT(2:end);
%1)    'cpi_usa'
%2)    'stir_usa'
%3)    'xrusd_usa'
%4)    'cpi_deu'
%5)    'stir_deu'%6)    'xrusd_deu'
%7)    'cpi_uk'
%8)    'stir_uk'
%9)    'xrusd_uk'
%10)   'cpi_fr'
%11)   'stir_fr'
%12)   'xrusd_fr'
%13)   'cpi_ca'
%14)   'stir_ca'
%15)   'xrusd_ca'
%16)   'cpi_it'
%17)   'stir_it'
%18)   'xrusd_it'
%19)   'cpi_jp'
%20)   'stir_jp'
%21)   'xrusd_jp'


Country = {'US','DE','UK','FR','CA','IT','JP'};
X = DATA(:,2:end);

Price_us = X(:,1);
Price_de = X(:,4);
Price_uk = X(:,7);
Price_fr = X(:,10);
Price_ca = X(:,13);
Price_it = X(:,16);
Price_jp = X(:,19);


Stir_us = X(:,2);
Stir_de = X(:,5);
Stir_uk  = X(:,8);
Stir_fr = X(:,11);
Stir_ca  = X(:,14);
Stir_it  = X(:,17);
Stir_jp  = X(:,20);

% Natl. currency for one dollar
Ner_us  = X(:,3);
Ner_de  = X(:,6);
Ner_uk   = X(:,9);
Ner_fr   = X(:,12);
Ner_ca   = X(:,15);
Ner_it   = X(:,18);
Ner_jp   = X(:,21);



%Inflation rate: national currency
Infl_us   = [NaN;(Price_us(2:end)./Price_us(1:end-1)-1)*100];
Infl_de   = [NaN;(Price_de(2:end)./Price_de(1:end-1)-1)*100];
Infl_uk   = [NaN;(Price_uk(2:end)./Price_uk(1:end-1)-1)*100];
Infl_fr   = [NaN;(Price_fr(2:end)./Price_fr(1:end-1)-1)*100];
Infl_ca   = [NaN;(Price_ca(2:end)./Price_ca(1:end-1)-1)*100];
Infl_it   = [NaN;(Price_it(2:end)./Price_it(1:end-1)-1)*100];
Infl_jp   = [NaN;(Price_jp(2:end)./Price_jp(1:end-1)-1)*100];


%%change of nominal exchange rate: positive is depreciation of national currency wrt USD
dNer_us   = [NaN;(Ner_us(2:end)./Ner_us(1:end-1)-1)*100];
dNer_de   = [NaN;(Ner_de(2:end)./Ner_de(1:end-1)-1)*100];
dNer_uk   = [NaN;(Ner_uk(2:end)./Ner_uk(1:end-1)-1)*100];
dNer_fr   = [NaN;(Ner_fr(2:end)./Ner_fr(1:end-1)-1)*100];
dNer_ca   = [NaN;(Ner_ca(2:end)./Ner_ca(1:end-1)-1)*100];
dNer_it   = [NaN;(Ner_it(2:end)./Ner_it(1:end-1)-1)*100];
dNer_jp   = [NaN;(Ner_jp(2:end)./Ner_jp(1:end-1)-1)*100];

%%change of real exchange rate: positive is depreciation of national
%%currency relative to US
dRer_us   = dNer_us - (Infl_us-Infl_us);
dRer_de   = dNer_de - (Infl_de-Infl_us);
dRer_uk   = dNer_uk - (Infl_uk-Infl_us);
dRer_fr   = dNer_fr - (Infl_fr-Infl_us);
dRer_ca   = dNer_ca - (Infl_ca-Infl_us);
dRer_it   = dNer_it - (Infl_it-Infl_us);
dRer_jp   = dNer_jp - (Infl_jp-Infl_us);





Y = [
    dRer_de...
    dRer_uk...
    dRer_fr...
    dRer_ca...
    dRer_it...
    dRer_jp
    ];

Y(abs(Y)>30)=NaN;


%%Some descriptive plots
figure(101)
plotyy(Year,log(Price_de)-log(Price_us),Year,log(Ner_de))
legend('p_{de} - p_{us}','e_{de}')
figure(102)
plotyy(Year,log(Price_uk)-log(Price_us),Year,log(Ner_uk))
legend('p_{uk} - p_{us}','e_{uk}')
figure(103)
plotyy(Year,log(Price_fr)-log(Price_us),Year,log(Ner_fr))
legend('p_{fr} - p_{us}','e_{fr}')
figure(104)
plotyy(Year,log(Price_ca)-log(Price_us),Year,log(Ner_ca))
legend('p_{ca} - p_{us}','e_{ca}')
figure(105)
plotyy(Year,log(Price_it)-log(Price_us),Year,log(Ner_it))
legend('p_{it} - p_{us}','e_{it}')
figure(106)
plotyy(Year,log(Price_jp)-log(Price_us),Year,log(Ner_jp))
legend('p_{jp} - p_{us}','e_{jp}')


%%

T0 = find(Year==1870);
T1 = find(Year==2016);


Y = Y(T0:T1,:);
Year = Year(T0:T1);
y=Y;
[T,n] = size(y);

%       dq_wrd
Ctr =[
    1%dRer_de
    1%dRer_uk
    1%dRer_fr
    1%dRer_ca
    1%dRer_it
    1%dRer_jp
    ];

%Addyng country specific trends to real rates (dq)
Cadd1                   =    eye(n); 

Ctr           = [Ctr Cadd1 ];
Ccyc          = zeros(n,n*p); 
Ccyc(1:n,1:n) = eye(n);
C             = [Ctr Ccyc];

r = size(Ctr,2);

b0          = zeros(n*p,n); 
b0(1:n,1:n) = eye(n)*0;

df0tr = 100;

%             dq_wrd         
SC0tr =    ([  1   1*ones(1,6)     ]).^2/100;
S0tr =      [  0    0*ones(1,6)]';
P0tr = diag([  1   ones(1,6)/2].^2);

%                stir          infl        ltir        
Psi =       (2*[ ones(1,6)   ]).^2;

S0cyc = zeros(n*p,1);

Atr  = eye(r);
Qtr  = diag(SC0tr);

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
A                  = zeros(r+n*p);
A(1:r,1:r)         = Atr;
A(r+1:end,r+1:end) = Acyc;

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

% Set prior to paper's lambda
mean_theta = 1;
std_theta = .5;


lambda = .2;  % Minnesota prior

logML = -inf;

%% Begin estimation
% For reference, see the appendix

for jm = 1:Ndraws  % MCMC draws
    
    % ------------------------- Block 1a ----------------------------------

    kf = KF(y, C, R, A, Q, S0, P0);
    loglik = kf.LogLik;
    
    theta_old = [];
    theta_old = [theta_old; C(1:6,1)] ;  % Adding the loadings of rates to the global m
    
    theta_new = theta_old + randn(size(theta_old))*std_theta/2;
    
    C_new        = C;
    C_new(1:6,1) = theta_new(1:6);  % adjusting loadingns to r_wrd


    kf_new     = KF(y, C_new, R, A, Q, S0, P0);
    loglik_new = kf_new.LogLik;
    
    log_rat = (loglik_new + sum(log(normpdf(theta_new, mean_theta, std_theta^2)))) ...
            - (loglik+      sum(log(normpdf(theta_old, mean_theta, std_theta^2))));
    p_acc = min(exp(log_rat),1);

    if rand<=p_acc
        C      = C_new;
        loglik = loglik_new;
        kf     = kf_new;
    end;
    
    % ----------------------------- Block 2a ------------------------------

    kc   = KC(kf);                     % Smoother
    Ytr  = [kc.S0(1:r)'; kc.S(:,1:r)]; % Trend components 
    Ycyc = kc.S(:,r+1:r+n);            % Cyclical component
    
    % ----------------------------- Block 2b ------------------------------
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
    
    P_acc(jm)= p_acc;
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
            end
        end
        
    end
end;

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

mStates = nanmean(States,3);

save(filename, '-v7.3')


