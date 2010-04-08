function [SC, NCORRECT_TE] = calc_SCval(CVRes, COND, noweight, thres)
% Calculate SC-value (selection count).
%
% -- Syntax
% [N, Ncorrect] = calc_Nval(CVRes)
% [N, NCORRECT_TE] = calc_Nval(CVRes, COND)
% [N, NCORRECT_TE] = calc_Nval(CVRes, COND, noweight)
%
% -- Input
% CVRes : Cross Validation Result structure
% .ix_eff_all  :
% .g           :
% .Ncorrect_te : test correct 
%      or
% .errTable_te : test correct table 
% Cond         : cell array of strings, this is just for diplaying purpose.
% noweight [0] : 1 if percent correct is not used, 0 if percent correct is used. 
%
% -- Ouput
% SC : SC-value
% NCORRECT_TE : percent correct in test (prediction performance)
%
% 2008/04/30 OY support 'noweight'
% 2006/09/12 OY support new output format
% 2006/09/01 OY
%
% Copyright (c) 2009, Okito Yamashita, ATR CNS, oyamashi@atr.jp.


Ncv = length(CVRes);
Nclass = 1;
Nvox = CVRes(1).g.nfeat-1;  % number of voxels

count = zeros(Nclass, Nvox);
NCORRECT_TE = zeros(Ncv,1);

if nargin < 4
    thres = 0;
end

if nargin < 3 
    noweight = 0;
end

if nargin < 2
  for ii = 1 : Nclass
      COND{ii} = '';
  end
end

for ii = 1 :  Ncv
    Res = CVRes(ii);
    if isfield(Res, 'Ncorrect_te'),
        Pcorrect_te = Res.Ncorrect_te/100;
    elseif isfield(Res, 'errTable_te')
        Pcorrect_te = sum(diag(Res.errTable_te))/sum(Res.errTable_te(:));
    else
        Pcorrect_te = sum(diag(Res.errTbale_te))/sum(Res.errTbale_te(:));
    end
        
    for cc = 1 : Nclass
        if ~isempty(COND{cc})
        fprintf('%s     : ', COND{cc});fprintf('%d ', Res.ix_eff_all{cc});
        fprintf('\n');
        end
        
        incr = zeros(1,Nvox);
        if noweight   % no weight
            incr(:, Res.ix_eff_all{cc}) = 1;
        else          % weight
            if Pcorrect_te >= thres  %% thresholding by channce level 
                incr(:, Res.ix_eff_all{cc}) = Pcorrect_te ;
            else
                incr(:, Res.ix_eff_all{cc}) = 0;                
            end
        end
        count(cc,:) = count(cc,:) + incr;
    end
    
    NCORRECT_TE(ii) = Pcorrect_te*100;
end
SC = count;
