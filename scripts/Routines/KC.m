function kc = KC(kf);

S  = kf.S;
P  = kf.P;
Sf = kf.Sf;
Pf = kf.Pf;
S0 = kf.S0;
P0 = kf.P0;
A = kf.A;

[T,ns] = size(S);
drS = zeros(T,ns)*NaN;
mu = S(T,:)';
sigma = P(:,:,T)';
sigma = .5*sigma+.5*sigma'+1e-6*eye(ns)*mean(diag(sigma));  % Symmetric, nonzero diagonal
draw = mu + chol(sigma)'*randn(size(mu));
drS(T,:) = draw';

for t=T-1:-1:1
    iPf = inv(Pf(:,:,t+1));
    mu =  S(t,:)'+P(:,:,t)*A'*iPf*(drS(t+1,:)-Sf(t+1,:))';
    sigma  =  P(:,:,t) - P(:,:,t)*A'*iPf*A*P(:,:,t);
    sigma = .5*sigma+.5*sigma'+1e-6*eye(ns)*mean(diag(sigma));
    draw = mu + chol(sigma)'*randn(size(mu));
    drS(t,:) = draw';
end;

iPf = inv(Pf(:,:,1));
mu =  S0+P0*A'*iPf*(drS(1,:)'-Sf(1,:)');
sigma  =  P0 - P0*A'*iPf*A*P0;
sigma = .5*sigma+.5*sigma'+1e-6*eye(ns)*mean(diag(sigma));
draw = mu + chol(sigma)'*randn(size(mu));
drS0 = draw;
kc.S  = drS;
kc.S0 = drS0;
