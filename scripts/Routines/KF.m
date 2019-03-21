%   - y[T x n]:  Series
%   - C[n x ns]: Maps series to states
%   - R[n x n]:  Variance covariance matrix of state equation
%   - A[ns x ns]:  Factor transition equation
%   - Q[ns x ns]:  Variance-covariance matrix of transition equation. 
%   - S0[ns x ns]: Initial states
%   - P0[ns x ns]: MSE matrix
function kf = KF(y,C,R,A,Q,S0,P0);

T = size(y,1);
ns = size(C,2);  % Number of states
Sprev = S0;
Pprev  = P0;

LogLik = 0;

 S = zeros(T,ns)*NaN;
 P = zeros(ns,ns,T)*NaN;
Sf = zeros(T,ns)*NaN;
Pf = zeros(ns,ns,T)*NaN;;
    
    
for t = 1:T  % Loop through time t
    Sft = A*Sprev;       % Forecast for states
    Pft = A*Pprev*A'+Q;  % MSE of forecast
    yt = y(t,:)';
    
    % Remove missing values
    Miss = isnan(yt);
    yt = yt(Miss==0);  
    Ct = C(Miss==0,:);
    Rt = R(Miss==0,Miss==0);
    
    
    yf = Ct*Sft;               % y_t|t-1
    iV  = inv(Ct*Pft*Ct'+Rt);  % inv(E[(y_t - y_t|t-1)(y_t - y_t|t-1)'])
    Gain = Pft*Ct'*iV;         % Kalman gain
    St  = Sft + Gain*(yt-yf);  
    Pt  = Pft - Gain*Ct*Pft;
    LogLik = LogLik + .5*log(det(iV))-.5*(yt-yf)'*iV*(yt-yf)-.5*(2*pi);
    

    S(t,:) = St';
    P(:,:,t) = Pt;
    Sf(t,:) = Sft';
    Pf(:,:,t) = Pft;
    Sprev = St;
    Pprev = Pt;
end;

kf.LogLik = LogLik;
kf.S=S;
kf.P=P;
kf.S0=S0;
kf.P0=P0;
kf.Sf=Sf;
kf.Pf=Pf;
kf.A=A;
kf.Q=Q;
kf.C =C;
kf.R = R;