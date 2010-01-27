function [ttr, xtr, tte, xte, g] = gen_simudata(MU, S, Ntr, Nte, varargin)
% Generate simulation data
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.
  
[Nfeat, Nclass] = size(MU);

if ndims(S) == 2,
    SS = repmat(S,[1,1,Nclass]);
else 
    SS = S;
end

Ntot = Ntr+Nte;

% simulation data generation
x = [];
t = [];
for ii = 1 : Nclass
   mutmp = MU(:,ii);
   Stmp  = squeeze(SS(:,:,ii));
   xtmp  = randmn(mutmp, Stmp, ceil(Ntot/Nclass));
   x = [x; xtmp'];
   t = [t; ii*ones(ceil(Ntot/Nclass),1)];
end
   
[ix_tr, ix_te] = separate_train_test(t, Ntr/Ntot);

xtr = x(ix_tr,:);
ttr = t(ix_tr);

xte = x(ix_te,:);
tte = t(ix_te);

if Nclass == 2,
    ttr = ttr - 1;
    tte = tte - 1;
end

g.MU = MU;
g.S  = S;
g.Ntr = Ntr;
g.Nte = Nte;
g.Nclass = Nclass;
g.Nfeat  = Nfeat;
