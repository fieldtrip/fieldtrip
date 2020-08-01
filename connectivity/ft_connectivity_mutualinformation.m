function output = ft_connectivity_mutualinformation(input, varargin)

% FT_CONNECTIVITY_MUTUALINFORMATION computes mutual information using
% either the information breakdown toolbox (ibtb), as described in Magri
% et al., BMC Neuroscience 2009, 1471-2202, or Robin Ince's Gaussian copula
% based parametric approach (gcmi). The function is a helper function for
% FT_CONNECTIVITYANALYSIS. As a standalone function, it could be used as
% follows:
%
% Use as
%   mi = ft_connectivity_mutualinformation(data, ...)
%
% The input data should be a Nchan x Nobservations matrix.
%
% Additional input arguments come as key-value pairs:
%   method     = string, 'ibtb' (default), or 'gcmi'.
%
% The following arguments pertain to the 'ibtb' method:
%   histmethod = The way that histograms are generated from the data. Possible values
%                are 'eqpop' (default), 'eqspace', 'ceqspace', 'gseqspace'.
%                See the help of the 'binr' function in the ibtb toolbox for more information.
%   numbin     = scalar value. The number of bins used to create the histograms needed for
%                the entropy computations
%   opts       = structure that is passed on to the 'information' function in the ibtb
%                toolbox. See the help of that function for more information.
%   refindx    = scalar value or 'all'. The channel that is used as 'reference channel'.
%
% The output contains the estimated mutual information between all channels and the reference channel(s).
%
% See also FT_CONNECTIVITYANALYSIS

% Copyright (C) 2016 Donders Institute, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

method      = ft_getopt(varargin, 'method',  'ibtb'); % can be gcmi
refindx     = ft_getopt(varargin, 'refindx', 'all', 1);
featureindx = ft_getopt(varargin, 'featureindx'); % this is only function if gcmi is the method
lags        = ft_getopt(varargin, 'lags',    0);  % shift of target w.r.t. source, in samples
sourcelags  = ft_getopt(varargin, 'sourcelags', []);
tra         = ft_getopt(varargin, 'tra',     []); % 1/0-matrix for multivariate combination Nnew x Norg, where Norg = size(input,1)
conditional = istrue(ft_getopt(varargin, 'conditional', false)); % this is only functional if gcmi is the method
combinelags = istrue(ft_getopt(varargin, 'combinelags', false)); % this is only functional if gcmi is the method, with di/dfi
precondition = istrue(ft_getopt(varargin, 'precondition', false)); % this copnorms only once, this is strictly speaking incorrect, but numerical difference are minor

if conditional && ~strcmp(method, 'gcmi')
  error('conditional mutual information can be only computed with ''gcmi'' as method');
end

% check whether the combined options work out
if ~isempty(tra)
  tra = full(tra)>0;
  if strcmp(method, 'ibtb') && ~isequal(tra,eye(size(tra,1))>0)
    error('method ''ibtb'' in combination with a non-identity ''tra'' is not possible');
  end
else
  tra = eye(size(input,1))>0;
end

% ensure that the refindx is numeric, defaults to 1:size(input,1), i.e. do
% all-to-all
if (ischar(refindx) && strcmp(refindx, 'all')) || isempty(refindx)
  refindx = (1:size(tra,1))';
end

% do not allow anything else than a scalar, or 1:nchan as refindx
if numel(refindx)~=1 && numel(refindx)==size(tra,1)
  % ensure column
  refindx = refindx(:);
elseif numel(refindx)~=1 && numel(refindx)~=size(tra,1)
  if ~conditional && isempty(featureindx)
    % this is for plain mi, which currently allows for more than 1 single
    % ref
    refindx = refindx(:)'; % ensure row
  else
    error('this variant of mi can only be computed using a single, or all channels as reference');
  end
end

switch method
  case 'ibtb'
    % check whether the required toolbox is available
    ft_hastoolbox('ibtb', 1);
    
    % set some options
    histmethod = ft_getopt(varargin, 'histmethod', 'eqpop');
    numbin     = ft_getopt(varargin, 'numbin',     10);
    
    % set some additional options that pertain to the algorithmic details of the
    % mutual information computation, see the documentation of ibtb
    opts        = ft_getopt(varargin, 'opts', []);
    opts.nt     = ft_getopt(opts, 'nt', []);
    opts.method = ft_getopt(opts, 'method', 'dr');
    opts.bias   = ft_getopt(opts, 'bias',   'pt');
    
    % deal with NaNs in the input data, e.g. trial boundaries
    finitevals = isfinite(input);
    
    nchans = size(tra,1); 
    n      = size(input, 2);
    output = zeros(nchans, numel(refindx), numel(lags)) + nan;
    
    % for each lag
    for m = 1:numel(lags)
      fprintf('computing mutualinformation for time lag in samples %d\n', lags(m));
      
      % get the samples for the relative shifts
      beg1 = max(0, lags(m))  + 1;
      beg2 = max(0, -lags(m)) + 1;
      n1   = n-abs(lags(m));
        
      end1 = beg1+n1-1;
      end2 = beg2+n1-1;
      
      for p = 1:numel(refindx)
        tmprefdata = nan(sum(tra(refindx(p),:)),n);
        tmprefdata(:, beg1:end1) = input(tra(refindx(p),:), beg2:end2);
        
        finitevals2 = sum(finitevals,1)&sum(isfinite(tmprefdata),1); % this conservatively takes only the non-nan samples across all input data channels
        
        tmpinput    = input(:,finitevals2);
        tmprefdata  = tmprefdata(:,finitevals2);
      
        % discretize signal1
        tmprefdata = binr(tmprefdata, sum(finitevals2), numbin, histmethod);
      
        for k = setdiff(1:size(tmpinput,1),refindx(p))
          signal2 = tmpinput(k,:);
        
          % represent signal2 in bins according to signal1's discretization
          R = zeros(1,3,numbin);
          for j = 1:numbin
            nr         = tmprefdata==j-1;
            opts.nt(j) = sum(nr);
            R(1, 1:opts.nt(j),j) = signal2(nr);
          end
        
          % discretize signal2 and compute mi
          R2 = binr(R, opts.nt', numbin, histmethod);
          output(k,p,m) = information(R2, opts, 'I'); % this computes mutual information
        end
      end
    end
    
  case 'gcmi'
    ft_hastoolbox('gcmi', 1);
    
    % set some options
    cmplx = ft_getopt(varargin, 'complex', 'complex'); % this is only used if data are complex-valued

    % deal with NaNs in the input data, e.g. trial boundaries
    finitevals = isfinite(input);
    
    % verify whether data is complex-valued, check the inputs, and adjust
    % the input data
    if ~all(imag(input(:))==0)
      % a tra deviating from I is currently not supported: ask Robin how to
      % deal with this, if possible at all
      if ~isequal(tra,eye(size(tra,1))>0)
        error('complex-valued input data in combination with multivariate signals is not supported');
      end
      switch cmplx
        case 'complex'
          % tease apart the real/imag parts, treat as 2D-variable, and
          % ensure the nans to behave
          complexrows = sum(imag(input)~=0,2)>0;
          input(~finitevals) = nan+1i.*nan;
          input = cat(1, real(input), imag(input(complexrows,:)));
          tra   = cat(2, tra, tra(:,complexrows));
          finitevals = cat(1, finitevals, finitevals(complexrows,:));
        case 'abs'
          % take the amplitude
          input = abs(input);
        case 'angle'
          % tease apart the real/imag parts, after amplitude normalization,
          % and ensure the nans to behave
          complexrows = sum(imag(input)~=0,2)>0;
          input(~finitevals) = nan+1i.*nan;
          input(complexrows,:) = input(complexrows,:)./abs(input(complexrows,:));
          input = cat(1, real(input), imag(input(complexrows,:)));
          tra   = cat(2, tra, tra(:,complexrows));
          finitevals = cat(1, finitevals, finitevals(complexrows,:));
        otherwise
          error('unsupported value for ''complex''');
      end
    end
  
    if combinelags
      otherlags = lags(2:end);
      lags      = lags(1);
    else
      otherlags = [];
    end
    
    nchans = size(tra,1); 
    n      = size(input, 2);
    if ~isempty(featureindx)
      if isempty(sourcelags)
        sourcelags = lags;
      end
      output = zeros(size(refindx,1), nchans, numel(sourcelags), numel(lags)) + nan;  
    else
      output = zeros(size(refindx,1), nchans, numel(lags)) + nan;
    end
    
    if precondition
      finitevalstmp = sum(isfinite(input))==size(input,1);
      input(:,finitevalstmp) = copnorm(input(:,finitevalstmp)')';
      input(:,~finitevalstmp) = nan;
      finitevals = isfinite(input);
    end
    
    % for each lag if combinelags is false
    for m = 1:numel(lags)
      
      if ~conditional && isempty(featureindx)
        fprintf('computing mutualinformation for time lag in samples %d\n', lags(m));
        % 'normal' mutual information between 2 sets of time series, this
        % allows for negative time lags, as well as a time lag of 0
        
        % get the samples for the relative shifts for the given lag
        beg1 = max(0, lags(m))  + 1;
        beg2 = max(0, -lags(m)) + 1;
        n1   = n-abs(lags(m));
        
        end1 = beg1+n1-1;
        end2 = beg2+n1-1;
        for p = 1:size(refindx,1)
          
          if ~isequal(tra, eye(size(tra,1)))
            tmpsource = nan(sum(tra(refindx(p),:)),n);
            tmpsource(:, beg1:end1) = input(tra(refindx(p),:), beg2:end2);
          else 
            tmpsource = nan(size(refindx,2),n);
            tmpsource(:, beg1:end1) = input(refindx(p,:), beg2:end2);
          end
            
          finitevals2 = sum(finitevals,1)&sum(isfinite(tmpsource),1); % this conservatively takes only the non-nan samples across all input data channels
          
          if ~precondition
            tmptarget  = copnorm(input(:,finitevals2)')';
            tmpsource  = copnorm(tmpsource(:,finitevals2)')';
          else
            tmptarget  = input(:,finitevals2);
            tmpsource  = tmpsource(:,finitevals2);
          end
          
          if ~isequal(tra,eye(size(tra,1)))
            for k = setdiff(1:size(tra,1),refindx(p))
              output(p,k,m) = mi_gg(tmptarget(tra(k,:),:)',tmpsource');%, false, true);
            end
          else
            output(p,:,m) = mi_gg_vec(tmptarget(:,:)',tmpsource',true,true);
          end
        end
      elseif conditional && isempty(featureindx)
        fprintf('computing directed information for time lag in samples %d\n', lags(m));
        
        % condition on the time-lagged version of the target signal, this
        % amounts to what cfg.method = 'di' in ft_connectivityanalysis, it
        % conditions mi between the target and the past of the source 
        % on the past of the target signal
        
        % get the samples for the relative shifts for the given lag
        beg1 = max(0, lags(m))  + 1;
        beg2 = max(0, -lags(m)) + 1;
        n1   = n-abs(lags(m));
        
        end1 = beg1+n1-1;
        end2 = beg2+n1-1;
        
        target_shifted                 = nan(size(input,1),n,numel(otherlags)+1);
        target_shifted(:, beg1:end1,1) = input(:, beg2:end2);
        for k = 1:numel(otherlags)
          % get the samples for the relative shifts for the given lag,
          % accumulate the lags
          otherbeg1 = max(0,  otherlags(k))  + 1;
          otherbeg2 = max(0, -otherlags(k))  + 1;
          n1        = n-abs(otherlags(k));
          
          otherend1 = otherbeg1+n1-1;
          otherend2 = otherbeg2+n1-1;
          target_shifted(:,otherbeg1:otherend1,k+1) = input(:, otherbeg2:otherend2);
        end
        finitevals2    = sum(finitevals,1)>0&sum(sum(isfinite(target_shifted),3),1)>0; % this conservatively takes only the non-nan samples across all input data channels
        
        % the following step is quite expensive computationally, but for
        % the conditioning all shifted versions of the target signal are
        % needed anyway, this bypasses the use of the gcmi toolbox
        if ~precondition
          target = copnorm(input(:,finitevals2)');
          target = bsxfun(@minus,target,mean(target,1));
        else
          target = input(:,finitevals2)';
          target = bsxfun(@minus, target, mean(target,1));
        end
        
        target_shifted = permute(target_shifted(:, finitevals2, :), [2 1 3]);
        if ~precondition
          for k = 1:size(target_shifted,3)
            target_shifted(:,:,k) = copnorm(target_shifted(:,:,k));
          end
          target_shifted = bsxfun(@minus,target_shifted,mean(target_shifted,1));
        end
        
        if isequal(tra, eye(size(tra,1)))
          % compute the covariance between all channels, and their shifted
          % versions only once, and then reorganize into a (Ntarget x
          % Nref) x 3 x 3 matrix
          C = transpose([target, target_shifted])*[target, target_shifted];
          C = C./(size(target,1)-1); 
          
          nt   = size(target,2);
          ns   = numel(refindx);
          
          cT  = diag(C(1:nt,1:nt)); % variance of the target signals
          cTs = diag(C(nt+(1:nt),nt+(1:nt))); % variance of the shifted target signals
          cTTs = diag(C(1:nt, nt+(1:nt))); % covariance between target and shifted target signals
          Cxyz = zeros(ns*nt,3,3);
          for p = 1:numel(refindx)
            ixp = find(tra(refindx(p),:));
            ix  = (p-1)*nt+(1:nt);
            
            Cxyz(ix,1,1) = cT;
            Cxyz(ix,2,1) = C(1:nt, nt+ixp);
            Cxyz(ix,3,1) = cTTs;
            Cxyz(ix,1,2) = Cxyz(ix,2,1);
            Cxyz(ix,2,2) = cTs(ixp);
            Cxyz(ix,3,2) = C(nt+(1:nt), nt+ixp);
            Cxyz(ix,1,3) = Cxyz(ix,3,1);
            Cxyz(ix,2,3) = Cxyz(ix,3,2);
            Cxyz(ix,3,3) = cTs;   
          end
          I = cov2cmi_ggg(Cxyz, size(target,1), true, [1 1 1]);
          output(:,:,m) = reshape(I,[],ns).';
        else
          for p = 1:numel(refindx)
            if ~isequal(tra,eye(size(tra,1)))
              tmpsource  = target_shifted(:,tra(refindx(p),:));
              for k = setdiff(1:size(tra,1),refindx(p))
                output(p,k,m) = cmi_ggg(target(:,tra(k,:)),tmpsource,target_shifted(:,tra(k,:),:), true, true);
              end
            else
              % this part of the code will not be reached, it is kept here
              % to refer to the 'original' implementation, as opposed to
              % the faster version referenced above.
              tmpsource = target_shifted(:,tra(refindx(p),:),:);
              output(p,:,m) = cmi_ggg_vec(target,tmpsource,target_shifted, true, true);
            end
          end
        end
      elseif conditional && ~isempty(featureindx)
        fprintf('computing directed feature information for time lag in samples %d\n', lags(m));
  
        % a featureindx has been specified, this refers to dfi
        % I(A(t1);F|B(t1)) + I(B(t2);F|B(t1)) - I(A(t1) B(t2);F|B(t1))
        if ~isempty(otherlags)
          error('only a single time lag is allowed in dfi');
        end
        
        %if numel(refindx)>1
        %  error('only a single refindx allowed');
        %end
        t_beg1 = max(0,  lags(m)) + 1; 
        t_beg2 = max(0, -lags(m)) + 1;
        t_end1 = n + 1 - t_beg2;
        t_end2 = n + 1 - t_beg1;
        
        % time-lagged version of the target signal,
        % positive lags here mean shifted w.r.t. feature
        target                   = nan(size(input,1),n);
        target(:, t_beg2:t_end2) = input(:, t_beg1:t_end1);
             
        for mm = 1:numel(sourcelags)
          if sourcelags(mm)>=lags(m)
            continue;
          end
          % we have an unshifted feature, a shifted target (with lags(m)),
          % and a shifted source (with sourcelags(mm))
          s_beg1 = max(0,  sourcelags(mm)) + 1; 
          s_beg2 = max(0, -sourcelags(mm)) + 1;
          s_end1 = n + 1 - s_beg2;
          s_end2 = n + 1 - s_beg1;
          
          % feature data
          feature = input(tra(featureindx,:),:);
                    
          % shifted target signals at the time lag of the source
          target_shifted                   = nan(size(input,1),n);
          target_shifted(:, s_beg2:s_end2) = input(:, s_beg1:s_end1);
          
          finitevals2 = sum(finitevals,1)>0&sum(isfinite(target),1)>0&sum(isfinite(target_shifted),1)>0&sum(isfinite(feature),1)>0; % this conservatively takes only the non-nan samples across all input data channels
          
          if ~precondition
            % the following step is quite expensive computationally if it
            % needs to be done each time
            tmptarget  = copnorm(target(:,finitevals2)'); % allow for the original target variable to be kept
            tmptarget  = bsxfun(@minus,tmptarget,mean(tmptarget,1))';
          
            target_shifted = copnorm(target_shifted(:,finitevals2)');
            target_shifted = bsxfun(@minus,target_shifted,mean(target_shifted,1))';
          
            % feature signal
            feature = copnorm(feature(:,finitevals2)');
            feature = bsxfun(@minus,feature,mean(feature,1))';
          else
            tmptarget  = target(:,finitevals2);
            tmptarget  = bsxfun(@minus,tmptarget,mean(tmptarget,2));
          
            target_shifted = target_shifted(:,finitevals2);
            target_shifted = bsxfun(@minus,target_shifted,mean(target_shifted,2));
          
            % feature signal
            feature = feature(:,finitevals2);
            feature = bsxfun(@minus,feature,mean(feature,2));
          end
          
          % time-lagged version of the source signal,
          %source = target_shifted(:,tra(refindx,:));
            
          if ~isequal(tra,eye(size(tra,1)))
            %  not implemented
            error('computation of directed information is not implemented if not all signals are univariate');
          else
            % compute the covariance between all channels, and their shifted
            % versions, and the feature only once, and then reorganize into a (Ntarget x
            % Nref) x 4 x 4 matrix, this bypasses the use of the gcmi
            % toolbox
            dat = cat(1, tmptarget, target_shifted, feature);
            C = dat*transpose(dat);
            %C = transpose([tmptarget, target_shifted, feature])*[tmptarget, target_shifted, feature];
            C = C./(size(tmptarget,1)-1);
            
            nt   = size(tmptarget,1);
            ns   = numel(refindx);
            
            cT  = diag(C(1:nt,1:nt)); % variance of the target signals
            cTs = diag(C(nt+(1:nt),nt+(1:nt))); % variance of the shifted target signals
            cTTs = diag(C(1:nt, nt+(1:nt))); % covariance between target and shifted target signals
            cFT  = C(1:nt, 2*nt+1);
            cFTs = C(nt+(1:nt), 2*nt+1);
            cF   = C(2*nt+1, 2*nt+1); % hard coded only a single feature!
            Cxyz = zeros(ns*nt,4,4);
            for p = 1:numel(refindx)
              ixp = find(tra(refindx(p),:));
              ix  = (p-1)*nt+(1:nt);
              
              Cxyz(ix,1,1) = cT;
              Cxyz(ix,2,1) = C(1:nt, nt+ixp);
              Cxyz(ix,3,1) = cTTs;
              Cxyz(ix,4,1) = cFT;
              Cxyz(ix,1,2) = Cxyz(ix,2,1);
              Cxyz(ix,2,2) = cTs(ixp);
              Cxyz(ix,3,2) = C(nt+(1:nt), nt+ixp);
              Cxyz(ix,4,2) = cFTs(ixp);
              Cxyz(ix,1,3) = Cxyz(ix,3,1);
              Cxyz(ix,2,3) = Cxyz(ix,3,2);
              Cxyz(ix,3,3) = cTs;
              Cxyz(ix,4,3) = cFTs;
              Cxyz(ix,1,4) = Cxyz(ix,4,1);
              Cxyz(ix,2,4) = Cxyz(ix,4,2);
              Cxyz(ix,3,4) = Cxyz(ix,4,3);
              Cxyz(ix,4,4) = cF;
            end
            % compute the three information components, exclude the source
            % and feature 'channels' to avoid potential numerical issues
            I1 = cov2cmi_ggg(Cxyz(:,[1 4 3],[1 4 3]), size(tmptarget,2), true, [1 1 1]); % T,F
            I2 = cov2cmi_ggg(Cxyz(:,[2 4 3],[2 4 3]), size(tmptarget,2), true, [1 1 1]); % S,F
            I3 = cov2cmi_ggg(Cxyz(:,[1 2 4 3],[1 2 4 3]), size(tmptarget,2), true, [2 1 1]);
            output(:,:,mm,m) = reshape(I1+I2-I3,[],ns).'; % equation 4 in Robin Ince's scientific reports paper.
          
%             % compute the three information components, exclude the source
%             % and feature 'channels' to avoid potential numerical issues
%             sel = tra(refindx,:)==0&tra(featureindx,:)==0;
%             I1  = cmi_ggg_vec(source,           feature, target_shifted(:,sel), true, true);
%             I2  = cmi_ggg_vec(tmptarget(:,sel), feature, target_shifted(:,sel), true, true);
%             I3  = cmi_ggg_vec(cat(3, tmptarget(:,sel), repmat(source, 1, sum(sel))), feature, target_shifted(:,sel), true, true);
          end
        end
        
      elseif ~conditional && ~isempty(featureindx)
        fprintf('computing co-information for time lag in samples %d\n', lags(m));
  
        % using a feature without conditioning will lead to the computation
        % of co-information; I(A(t1);F) + I(B(t2);F) - I(A(t1) B(t2);F) ->
        % positive values: redundancy, negative values: synergy
        if any(lags<0)
          error('only lags >=0 allowed');
        end
%         if numel(refindx)>1
%           error('only a single refindx allowed');
%         end
        t_beg1 = max(0,  lags(m)) + 1; 
        t_beg2 = max(0, -lags(m)) + 1;
        
        t_end1 = n+1-t_beg2;
        t_end2 = n+1-t_beg1;
        
        % time-lagged version of the target signal,
        % positive lags here mean shifted w.r.t. feature
        target                   = nan(size(input,1),n);
        target(:, t_beg2:t_end2) = input(:, t_beg1:t_end1);
        
        % if copula is done here, it's not fully correct (because of
        % nanning of shifted values...
        sel = sum(isfinite(target),1)>0;
        target(:,sel) = copnorm(target(:,sel)')';
        target(:,sel) = bsxfun(@minus,target(:,sel),mean(target(:,sel),2));
        
        for mm = 1:numel(sourcelags)
          if sourcelags(mm)>=lags(m)
            continue;
          end
          
          % we have an unshifted feature, a shifted target (with lags(m)),
          % and a shifted source (with lags(mm))
          s_beg1 = max(0,  lags(mm)) + 1; 
          s_beg2 = max(0, -lags(mm)) + 1; 
          
          s_end1 = n+1-s_beg2;
          s_end2 = n+1-s_beg1;
          
          % feature data
          feature = input(tra(featureindx,:),:);
          
          % time-lagged version of the source signal,
          %source                   = nan(sum(tra(refindx,:)),n);
          %source(:, s_beg2:s_end2) = input(tra(refindx,:), s_beg1:s_end1);
          source                   = nan(size(input,1),n);
          source(:, s_beg2:s_end2) = input(:, s_beg1:s_end1);
          
          finitevals2    = sum(finitevals,1)>0&sum(isfinite(target),1)>0&sum(isfinite(source),1)>0&sum(isfinite(feature),1)>0; % this conservatively takes only the non-nan samples across all input data channels
          
          if ~precondition
            % the following step is quite expensive computationally if it
            % needs to be done each time
            tmptarget  = copnorm(target(:,finitevals2)'); % allow for the original target variable to be kept
            tmptarget  = bsxfun(@minus,tmptarget,mean(tmptarget,1));
          
            source = copnorm(source(:,finitevals2)');
            source = bsxfun(@minus,source,mean(source,1));
          
            % feature signal
            feature = copnorm(feature(:,finitevals2)');
            feature = bsxfun(@minus,feature,mean(feature,1));
          else
            tmptarget  = target(:,finitevals2)';
            tmptarget  = bsxfun(@minus,tmptarget,mean(tmptarget,1));
          
            source = source(:,finitevals2)';
            source = bsxfun(@minus,source,mean(source,1));
          
            % feature signal
            feature = feature(:,finitevals2)';
            feature = bsxfun(@minus,feature,mean(feature,1));
          end
          
          if ~isequal(tra,eye(size(tra,1)))
            % tmprefdata  = target_shifted(:,tra(refindx(p),:));
            % for k = setdiff(1:size(tra,1),refindx(p))
            %   output(k,p,m) = cmi_ggg(target(:,tra(k,:)),tmprefdata,target_shifted(:,tra(k,:)), true, false);
            % end
            
            %  not implemented
            error('computation of co-information is not implemented if not all signals are univariate');
          else
            % compute the covariance between all channels, and their shifted
            % versions, and the feature only once, and then reorganize into a (Ntarget x
            % Nref) x 4 x 4 matrix, this bypasses the use of the gcmi
            % toolbox
            C = transpose([tmptarget, source, feature])*[tmptarget, source, feature];
            C = C./(size(tmptarget,1)-1);
            
            nt   = size(tmptarget,2);
            ns   = numel(refindx);
            
            cT  = diag(C(1:nt,1:nt)); % variance of the target signals
            cS  = diag(C(nt+(1:nt),nt+(1:nt))); % variance of the source signals
            %cTS = diag(C(1:nt, nt+(1:nt))); % covariance between target and source signals
            cFT = C(1:nt, 2*nt+1);
            cFS = C(nt+(1:nt), 2*nt+1);
            cF  = C(2*nt+1, 2*nt+1); % hard coded only a single feature!
            Cxy = zeros(ns*nt,3,3);
            for p = 1:numel(refindx)
              ixp = find(tra(refindx(p),:));
              ix  = (p-1)*nt+(1:nt);
              
              Cxy(ix,1,1) = cT;
              Cxy(ix,2,1) = C(1:nt, nt+ixp);
              Cxy(ix,3,1) = cFT;
              Cxy(ix,1,2) = Cxy(ix,2,1);
              Cxy(ix,2,2) = cS(ixp);
              Cxy(ix,3,2) = cFS(ixp);
              Cxy(ix,1,3) = Cxy(ix,3,1);
              Cxy(ix,2,3) = Cxy(ix,3,2);
              Cxy(ix,3,3) = cF;
            end
            % compute the three information components, exclude the source
            % and feature 'channels' to avoid potential numerical issues
            I1 = cov2mi_gg(Cxy(:,[1 3],[1 3]), size(tmptarget,1), true, [1 1]); % T,F
            I2 = cov2mi_gg(Cxy(:,[2 3],[2 3]), size(tmptarget,1), true, [1 1]); % S,F
            I3 = cov2mi_gg(Cxy(:,[1 2 3],[1 2 3]), size(tmptarget,1), true, [2 1]);
            output(:,:,mm,m) = reshape(I1+I2-I3,[],ns).'; % equation 4 in Robin Ince's scientific reports paper.
            
%             % compute the three information components
%             I1 = mi_gg_vec(source,    feature, true, true);
%             I2 = mi_gg_vec(tmptarget, feature, true, true);
%             I3 = mi_gg_vec(cat(3, tmptarget, repmat(source, 1, size(tmptarget,2))), feature, true, true);
%             output(1,:,mm,m) = I1+I2-I3;
          end
        end
      end
    end
  otherwise
end

if size(refindx,1)==1 %&& ~(~conditional && ~isempty(featureindx))
  siz    = [size(output) 1];
  output = reshape(output,[siz(2:end)]);
end

function I = cov2mi_gg(Cxy, N, biascorrect, xy_ind)

% subfunction that computes mutual information, where the covariances have 
% already been computed.
%  Cxy = covariance, MxNsgnxNsgn
%  N = number of samples for covariance computation (needed for bias
%  estimate)
% biascorrect = boolean
% xy_ind = [nx ny nz], vector with dimensionalities of x,y,z
%
% conditional MI is computed between x and y, conditioned on z

xindx = 1:xy_ind(1);
yindx = xindx(end)+(1:xy_ind(2));

% submatrices of joint covariance
Cx  = Cxy(:,xindx, xindx);
Cy  = Cxy(:,yindx, yindx);
Cxy = Cxy(:,[xindx yindx], [xindx yindx]);

chCx   = vecchol(Cx);
chCy   = vecchol(Cy);
chCxy  = real(vecchol(Cxy));

% entropies in nats
% normalisations cancel for mi
HX   = sum(log(vecdiag(chCx)),  2); % + 0.5*Nvarz*log(2*pi*exp(1));
HY   = sum(log(vecdiag(chCy)),  2); % + 0.5*(Nvarx+Nvarz)*log(2*pi*exp(1));
HXY  = sum(log(vecdiag(chCxy)), 2); % + 0.5*(Nvary+Nvarz)*log(2*pi*exp(1));

ln2 = log(2);
if biascorrect
    nX  = numel(xindx);
    nY  = numel(yindx);
    nXY = numel(yindx)+numel(xindx);
    
    psiterms = psi((N - (1:nXY))/2) / 2;
    dterm = (ln2 - log(N-1)) / 2;
    HX   = (HX   - nX*dterm    - sum(psiterms(1:nX)));
    HY   = (HY   - nY*dterm    - sum(psiterms(1:nY)));
    HXY  = (HXY  - nXY*dterm   - sum(psiterms(1:nXY)));
end

% convert to bits
I = (HX + HY - HXY) / ln2;

function I = cov2cmi_ggg(Cxyz, N, biascorrect, xyz_ind)

% subfunction that computes mutual information, conditioned on a third
% variable, where the covariances have already been computed.
%  Cxyz = covariance, MxNsgnxNsgn
%  N = number of samples for covariance computation (needed for bias
%  estimate)
% biascorrect = boolean
% xyz_ind = [nx ny nz], vector with dimensionalities of x,y,z
%
% conditional MI is computed between x and y, conditioned on z

xindx = 1:xyz_ind(1);
yindx = xindx(end)+(1:xyz_ind(2));
zindx = yindx(end)+(1:xyz_ind(3));

% submatrices of joint covariance
Cz  = Cxyz(:,zindx, zindx);
Cyz = Cxyz(:,[yindx zindx], [yindx zindx]);
Cxz = Cxyz(:,[xindx zindx], [xindx zindx]);

chCz   = vecchol(Cz);
chCxz  = vecchol(Cxz);
chCyz  = real(vecchol(Cyz));
chCxyz = real(vecchol(Cxyz));

% entropies in nats
% normalisations cancel for cmi
HZ   = sum(log(vecdiag(chCz)),  2); % + 0.5*Nvarz*log(2*pi*exp(1));
HXZ  = sum(log(vecdiag(chCxz)), 2); % + 0.5*(Nvarx+Nvarz)*log(2*pi*exp(1));
HYZ  = sum(log(vecdiag(chCyz)), 2); % + 0.5*(Nvary+Nvarz)*log(2*pi*exp(1));
HXYZ = sum(log(vecdiag(chCxyz)),2); % + 0.5*(Nvarx+Nvary+Nvarz)*log(2*pi*exp(1));

ln2 = log(2);
if biascorrect
    nZ = numel(zindx);
    nXZ = numel(zindx)+numel(xindx);
    nYZ = numel(zindx)+numel(yindx);
    nXYZ = numel(zindx)+numel(yindx)+numel(xindx);
  
    psiterms = psi((N - (1:nXYZ))/2) / 2;
    dterm = (ln2 - log(N-1)) / 2;
    HZ   = (HZ   - dterm       - sum(psiterms(1:nZ)));
    HXZ  = (HXZ  - nXZ*dterm   - sum(psiterms(1:nXZ)));
    HYZ  = (HYZ  - nYZ*dterm   - sum(psiterms(1:nYZ)));
    HXYZ = (HXYZ - nXYZ*dterm  - sum(psiterms));
end

% convert to bits
I = (HXZ + HYZ - HXYZ - HZ) / ln2;

function out = vecdiag(in)

n = size(in,2);
out = zeros(size(in,1),n);
for k = 1:n
  out(:,k) = in(:,k,k);
end

