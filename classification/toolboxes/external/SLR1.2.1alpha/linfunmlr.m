function [f,g,H] = linfunmlr(w, t, ax, X, ixf, ixc, Nclass)
% linear function handle for multinomial logistic regression
% 
% t : integer representing class label. must start from 1
%
% 2009/06/01 OY Revised 
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.

[Nsamp] = size(X,1);

Z = [];
for c = 1 : Nclass
    % calculate classification values
    ixix = find(ixc == c);
    if isempty(ixix)
        Y(:,c) = zeros(Nsamp,1);
    else
        Y(:,c) = X(:,ixf(ixix))*w(ixix);    % Nsamp*Nclass
    end
    % indicator matrix
    tmp=zeros(Nsamp,1);
    tmp(t==c) = 1;
    Z = [Z tmp];  % indicator matrix, Nsamp*Nclass    \
end

% probability matrix
eY = exp(Y); % Nsamp*Nclass
P = eY ./ repmat(sum(eY,2), [1, Nclass]); % Nsamp*Nclass

%%% Value -- negative log-likelihood
fone = sum(Y.*Z,2) - log(sum(eY,2));  % sum over class
f = sum(fone)-1/2*sum(ax.*w.^2);  % sum over Nsamp
f = -f;

%%%% Gradient
if nargout > 1
    del = Z - P; % Nsamp*Nclass
    g = [];
    axax = [];
    for c = 1 : Nclass
        gc = del(:,c); %Nsamp*1
        ixix = find(ixc == c);
        tmp = X(:,ixf(ixix))' * gc - ax(ixix).* w(ixix); 
        g = [g;sum(tmp,2)] ;
        axax = [axax; ax(ixix)];
    end
    g = -g;
end

% %%%% Hessian Information (sparse part only)
if nargout > 2
    H = [];
    for c1 = 1 : Nclass
        Htmp = [];
        for c2 = 1 : Nclass
            if c1 == c2
                D = diag(P(:,c1).*(1-P(:,c1)));
            else
                D = -diag(P(:,c1).*P(:,c2));
            end
            
            %ixix1 = find(ixc == c1);
            %ixix2 = find(ixc == c2);
            %x1 = X(:,ixf(ixix1));  % Nsamp * Nixeff1
            %x2 = X(:,ixf(ixix2));  % Nsamp * Nixeff2
        
            x1 = X(:,ixf(ixc==c1));  % Nsamp * Nixeff1
            x2 = X(:,ixf(ixc==c2));  % Nsamp * Nixeff2
        
%             size(x1)
%             size(x2)
%             size(D)
            
            Htmp = [Htmp, x1'*D*x2];
        end  % c2
        H = [H; Htmp];
    end  % c1
    H = H +diag(axax);
end
