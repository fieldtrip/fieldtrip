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

method  = ft_getopt(varargin, 'method',  'ibtb'); % can be gcmi
refindx = ft_getopt(varargin, 'refindx', 'all', 1);
featureindx = ft_getopt(varargin, 'featureindx'); % this is only function if gcmi is the method
lags    = ft_getopt(varargin, 'lags',    0);  % shift of data w.r.t. reference, in samples
featurelags = ft_getopt(varargin, 'featurelags', []);
tra     = ft_getopt(varargin, 'tra',     []); % 1/0-matrix for multivariate combination Nnew x Norg, where Norg = size(input,1)
conditional = istrue(ft_getopt(varargin, 'conditional', false)); % this is only functional if gcmi is the method
combinelags = istrue(ft_getopt(varargin, 'combinelags', false)); % this is only functional if gcmi is the method, with di/dfi
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
if numel(refindx)~=1 && numel(refindx)~=size(tra,1)
  error('mi can only be computed using a single, or all channels as reference');
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
    output = zeros(nchans, numel(refindx), numel(lags)) + nan;
    
    % for each lag if combinelags is false
    for m = 1:numel(lags)
      fprintf('computing mutualinformation for time lag in samples %d\n', lags(m));
      
      % get the samples for the relative shifts for the given lag
      beg1 = max(0, lags(m))  + 1;
      beg2 = max(0, -lags(m)) + 1;
      n1   = n-abs(lags(m));
        
      end1 = beg1+n1-1;
      end2 = beg2+n1-1;
      if ~conditional && isempty(featureindx)
        % 'normal' mutual information between 2 sets of time series
        for p = 1:numel(refindx)
          tmprefdata = nan(sum(tra(refindx(p),:)),n);
          tmprefdata(:, beg1:end1) = input(tra(refindx(p),:), beg2:end2);
          
          finitevals2 = sum(finitevals,1)&sum(isfinite(tmprefdata),1); % this conservatively takes only the non-nan samples across all input data channels
          
          tmpinput    = copnorm(input(:,finitevals2)')';
          tmprefdata  = copnorm(tmprefdata(:,finitevals2)')';
          
          if ~isequal(tra,eye(size(tra,1)))
            for k = setdiff(1:size(tra,1),refindx(p))
              output(k,p,m) = mi_gg(tmpinput(tra(k,:),:)',tmprefdata');%, false, true);
            end
          else
            output(:,p,m) = mi_gg_vec(tmpinput(:,:)',tmprefdata',true,true);
          end
        end
      elseif conditional && isempty(featureindx)
        % condition on the time-lagged version of the target signal
        target_shifted               = nan(size(input,1),n,numel(otherlags)+1);
        target_shifted(:, beg1:end1,1) = input(:, beg2:end2);
        for k = 1:numel(otherlags)
          % get the samples for the relative shifts for the given lag,
          % accumulate the lags
          otherbeg1 = max(0,  otherlags(k))  + 1;
          otherbeg2 = max(0, -otherlags(k))  + 1;
          n1   = n-abs(otherlags(k));
          
          otherend1 = otherbeg1+n1-1;
          otherend2 = otherbeg2+n1-1;
          target_shifted(:,otherbeg1:otherend1,k+1) = input(:, otherbeg2:otherend2);
        end
        
        finitevals2    = sum(finitevals,1)>0&sum(sum(isfinite(target_shifted),3),1)>0; % this conservatively takes only the non-nan samples across all input data channels
        
        % the following step is quite expensive computationally, but for
        % the conditioning all shifted versions of the target signal are
        % needed anyway
        target         = copnorm(input(:,finitevals2)');
        target         = bsxfun(@minus,target,mean(target,1));
        target_shifted = permute(target_shifted(:, finitevals2, :), [2 1 3]);
        for k = 1:size(target_shifted,3)
          target_shifted(:,:,k) = copnorm(target_shifted(:,:,k));
        end
        target_shifted = bsxfun(@minus,target_shifted,mean(target_shifted,1));
        
        for p = 1:numel(refindx)
          if ~isequal(tra,eye(size(tra,1)))
            tmprefdata  = target_shifted(:,tra(refindx(p),:));
            for k = setdiff(1:size(tra,1),refindx(p))
              output(k,p,m) = cmi_ggg(target(:,tra(k,:)),tmprefdata,target_shifted(:,tra(k,:),:), true, true);
            end
          else
            tmprefdata = target_shifted(:,tra(refindx(p),:),:);
            output(:,p,m) = cmi_ggg_vec(target,tmprefdata,target_shifted, true, true);
          end
        end
      elseif conditional && ~isempty(featureindx)
        % a featureindx has been specified, this refers to dfi
        if ~isempty(otherlags)
          error('only a single time lag is allowed in dfi');
        end
        if isempty(featurelags)
          featurelags = lags;
        end
        if ~isequal(featurelags, lags)
          featurebeg1 = max(0,  featurelags)  + 1;
          featurebeg2 = max(0, -featurelags)  + 1;
          n1   = n-abs(featurelags);
          
          featureend1 = featurebeg1+n1-1;
          featureend2 = featurebeg2+n1-1;
        else
          featurebeg1 = beg1;
          featureend1 = end1;
          featurebeg2 = beg2;
          featureend2 = end2;
        end
        
        % feature data
        feature = nan(sum(tra(featureindx,:)),n);
        feature(:, featurebeg1:featureend1) = input(tra(featureindx,:),featurebeg2:featureend2);
        
        % condition on the time-lagged version of the target signal
        target_shifted               = nan(size(input,1),n);
        target_shifted(:, beg1:end1) = input(:, beg2:end2);
        
        finitevals2    = sum(finitevals,1)>0&sum(isfinite(target_shifted),1)>0&sum(isfinite(feature),1)>0; % this conservatively takes only the non-nan samples across all input data channels
        
        % the following step is quite expensive computationally, but for
        % the conditioning all shifted versions of the target signal are
        % needed anyway
        target         = copnorm(input(:,finitevals2)');
        target         = bsxfun(@minus,target,mean(target,1));
        target_shifted = copnorm(target_shifted(:,finitevals2)');
        target_shifted = bsxfun(@minus,target_shifted,mean(target_shifted,1));
        feature        = copnorm(feature(:,finitevals2)');
        feature        = bsxfun(@minus,feature,mean(feature,1));
        %feature        = target_shifted(:,tra(featureindx,:));
        
        for p = 1:numel(refindx)
          if ~isequal(tra,eye(size(tra,1)))
%             tmprefdata  = target_shifted(:,tra(refindx(p),:));
%             for k = setdiff(1:size(tra,1),refindx(p))
%               output(k,p,m) = cmi_ggg(target(:,tra(k,:)),tmprefdata,target_shifted(:,tra(k,:)), true, false);
%             end
             keyboard
          else
            tmprefdata = target_shifted(:,tra(refindx(p),:));
            
            % compute the three information components
            I1 = cmi_ggg_vec(tmprefdata, feature, target_shifted, true, true);
            I2 = cmi_ggg_vec(target,     feature, target_shifted, true, true);
            I3 = cmi_ggg_vec(cat(3, target, repmat(tmprefdata, 1, size(target,2))), feature, target_shifted, true, true);
            output(:,p,m) = I1+I2-I3;
          end
        end
      elseif ~conditional && ~isempty(featureindx)
        % using a feature without conditioning will lead to the computation
        % of co-information; I(A(t1);F) + I(B(t2);F) - I(A(t1) B(t2);F) ->
        % positive values: redundancy, negative values: synergy
        if any(lags<0)
          error('only lags >=0 allowed');
        end
        if numel(refindx)>1
          error('only a single refindx allowed');
        end
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
        
        for mm = 1:numel(lags)
          % we have an unshifted feature, a shifted target (with lags(m)),
          % and a shifted source (with lags(mm))
          s_beg1 = max(0,  lags(mm)) + 1; 
          s_beg2 = max(0, -lags(mm)) + 1; 
          
          s_end1 = n+1-s_beg2;
          s_end2 = n+1-s_beg1;
          
          % feature data
          feature = input(tra(featureindx,:),:);
          
          
          % time-lagged version of the source signal,
          source                   = nan(sum(tra(refindx,:)),n);
          source(:, s_beg2:s_end2) = input(tra(refindx,:), s_beg1:s_end1);
          
          finitevals2    = sum(finitevals,1)>0&sum(isfinite(target),1)>0&sum(isfinite(source),1)>0; % this conservatively takes only the non-nan samples across all input data channels
          
          % the following step is quite expensive computationally
          %target  = copnorm(target(:,finitevals2)');
          %target  = bsxfun(@minus,target,mean(target,1));
          tmptarget  = target(:,finitevals2)';
          tmptarget  = bsxfun(@minus,tmptarget,mean(tmptarget,1));
          source  = copnorm(source(:,finitevals2)');
          source  = bsxfun(@minus,source,mean(source,1));
          feature = copnorm(feature(:,finitevals2)');
          feature = bsxfun(@minus,feature,mean(feature,1));
          
          
          if ~isequal(tra,eye(size(tra,1)))
              %             tmprefdata  = target_shifted(:,tra(refindx(p),:));
              %             for k = setdiff(1:size(tra,1),refindx(p))
              %               output(k,p,m) = cmi_ggg(target(:,tra(k,:)),tmprefdata,target_shifted(:,tra(k,:)), true, false);
              %             end
              keyboard
          else
            % compute the three information components
            I1 = mi_gg_vec(source,    feature, true, true);
            I2 = mi_gg_vec(tmptarget, feature, true, true);
            I3 = mi_gg_vec(cat(3, tmptarget, repmat(source, 1, size(tmptarget,2))), feature, true, true);
            output(:,mm,m) = I1+I2-I3;
          end
        end
      end
    end
  otherwise
end

if numel(refindx)==1 && ~(~conditional && ~isempty(featureindx))
  siz    = [size(output) 1];
  output = reshape(output,[siz([1 3])]);
end
