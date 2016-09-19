function [data, powindx, hasrpt] = univariate2bivariate(data, inparam, outparam, dtype, varargin)

% UNIVARIATE2BIVARIATE is a helper function for FT_CONNECTIVITYANALYSIS
%
% Use as
%   [data, powindx, hasrpt] = univariate2bivariate(data, inparam, outparam, dtype, ...)
% where
%   data        = FieldTrip structure according to dtype (see below)
%   inparam     = string
%   inparam     = string
%   dtype       = string, can be 'freq', 'source', 'raw'
% and additional options come in key-value pairs and can include
%   cmb         = 
%   demeanflag  = 
%   keeprpt     = 
%   sqrtflag    = 

% Copyright (C) 2009-2012, Jan-Mathijs Schoffelen
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

cmb         = ft_getopt(varargin, 'cmb');
demeanflag  = ft_getopt(varargin, 'demeanflag', false);
keeprpt     = ft_getopt(varargin, 'keeprpt',    true);
sqrtflag    = ft_getopt(varargin, 'sqrtflag',   false);

switch dtype
  case 'freq'
    nchan = numel(data.label);
    if isequal(cmb, {'all' 'all'})
      ncmb = nchan^2;
    else
      ncmb = size(cmb,1);
    end
    getpowindx = 0;
    if ncmb==0,
      error('no channel combinations are specified');
    elseif ncmb==nchan^2 || ncmb==(nchan+1)*nchan*0.5,
      dofull = 1;
    else
      dofull = 0;
    end
    
    if strcmp(inparam, 'fourierspctrm') && strcmp(outparam, 'crsspctrm'),
      % fourier coefficients -> cross-spectral density
      if dofull
        data = ft_checkdata(data, 'cmbrepresentation', 'full');
      else
        data = ft_checkdata(data, 'cmbrepresentation', 'sparse', 'channelcmb', cmb);
        getpowindx = 1;
      end
      
    elseif strcmp(inparam, 'powandcsd') && strcmp(outparam, 'crsspctrm'),
      if ~isempty(cmb),
        data = ft_checkdata(data, 'cmbrepresentation', 'sparse', 'channelcmb', cmb);
        % ensure getting powindx later on to prevent crash
        getpowindx = 1;
      else
        % data = ft_checkdata(data, 'cmbrepresentation', 'full');
        % this should not be possible
        error('cannot convert to a full csd representation');
      end
      
    elseif strcmp(inparam, 'fourierspctrm') && strcmp(outparam, 'powcovspctrm'),
      % fourier coefficients -> power covariance
      data = ft_checkdata(data, 'cmbrepresentation', 'sparsewithpow', 'channelcmb', {});
      if sqrtflag, data.powspctrm = sqrt(data.powspctrm); end
      % get covariance by using ft_checkdata
      if demeanflag,
        nrpt = size(data.powspctrm,1);
        mdat = nanmean(data.powspctrm,1);
        data.powspctrm = data.powspctrm - mdat(ones(1,nrpt),:,:,:,:,:);
      end
      data.fourierspctrm = data.powspctrm; % this is necessary for ft_checkdata to work
      data.dimord = ['rpttap',data.dimord(4:end)];
      data = rmfield(data, 'powspctrm');
      data.cumtapcnt(:) = 1;
      data.cumsumcnt(:) = 1;
      if ncmb < (nchan-1)*nchan*0.5,
        data = ft_checkdata(data, 'cmbrepresentation', 'sparse', 'channelcmb', cmb);
      else
        data = ft_checkdata(data, 'cmbrepresentation', 'full');
      end
      data.powcovspctrm = data.crsspctrm;
      data = rmfield(data, 'crsspctrm');
      
    elseif strcmp(inparam, 'powspctrm') && strcmp(outparam, 'powcovspctrm'),
      % power-spectral density -> power covariance
      if sqrtflag, data.powspctrm = sqrt(data.powspctrm); end
      % get covariance by using ft_checkdata
      if demeanflag,
        nrpt = size(data.powspctrm,1);
        mdat = nanmean(data.powspctrm,1);
        data.powspctrm = data.powspctrm - mdat(ones(1,nrpt),:,:,:,:,:);
      end
      data.fourierspctrm = data.powspctrm; % this is necessary for ft_checkdata to work
      data.dimord = ['rpttap',data.dimord(4:end)];
      data = rmfield(data, 'powspctrm');
      data.cumtapcnt(:) = 1;
      data.cumsumcnt(:) = 1;
      if ncmb < (nchan-1)*nchan*0.5,
        data = ft_checkdata(data, 'cmbrepresentation', 'sparse', 'channelcmb', cmb);
      else
        data = ft_checkdata(data, 'cmbrepresentation', 'full');
      end
      data.powcovspctrm = data.crsspctrm;
      data = rmfield(data, 'crsspctrm');
      
    else
      error('unknown conversion from univariate to bivariate representation');
    end % if inparam is fourierspctrm or crsspctrm
    
    if ~isempty(cmb) && (ncmb < (nchan-1)*nchan*0.5 || getpowindx==1),
      powindx = labelcmb2indx(data.labelcmb);
    else
      powindx = [];
    end
    
  case 'source'
    ncmb = numel(cmb);
    
    % the code further down requires this to be a vector with indices
    data = fixinside(data, 'index');
    
    if strcmp(inparam, 'pow') && strcmp(outparam, 'powcov'),
      [nvox,nrpt] = size(data.pow);
      if sqrtflag, data.pow = sqrt(data.pow); end
      if demeanflag,
        mdat = nanmean(data.pow,2);
        data.pow = data.pow - mdat(:,ones(1,nrpt)); % FIXME only works for 1 frequency
      end
      
      if ncmb == size(data.pow,1)
        data.powcov = data.pow * data.pow';
        data.powcovdimord = 'pos_pos';
        powindx = [];
      else
        data.powcov = [reshape(data.pow * data.pow(cmb,:)', [ncmb*nvox 1]); sum(data.pow.^2,2)];
        data = rmfield(data, 'pow');
        data = rmfield(data, 'powdimord');
        % powindx = [nvox+(1:nvox) nvox+(1:nvox); cmb*ones(1,nvox) nvox+(1:nvox)]';
        powindx = [repmat(ncmb*nvox+(1:nvox)',[ncmb 1]) reshape(repmat(ncmb*nvox+cmb(:)', [nvox 1]),[nvox*ncmb 1]); ncmb*nvox+(1:nvox)' ncmb*nvox+(1:nvox)'];
        data.pos = [repmat(data.pos, [ncmb+1 1])];% FIXME come up with something reshape( repmat(data.pos(cmb,:),[nvox 1]);data.pos data.pos];
        data.inside = reshape(repmat(data.inside(:), [1 ncmb+1])+repmat(nvox*(0:ncmb), [nvox 1]), [nvox*(ncmb+1) 1]);
        if ~isempty(data.outside)
          data.outside = reshape(repmat(data.outside(:), [1 ncmb+1])+repmat(nvox*(0:ncmb), [nvox 1]), [nvox*(ncmb+1) 1]);
        end
        data.powcovdimord = 'pos';
        % data.dim(2) = size(data.pos,1);
      end
    
    elseif strcmp(inparam, 'mom') && strcmp(outparam, 'powcov'),
      
      nvox = size(data.pos,1);
      if isfield(data, 'cumtapcnt')
        cumtapcnt = data.cumtapcnt;
      else
        cumtapcnt = ones(size(data.mom{find(data.inside,1,'first')},2),1);
      end
      nrpt = size(cumtapcnt,1);
      
      % make projection matrix to get from mom to pow
      vec = zeros(0,1);
      i1  = zeros(0,1);
      i2  = zeros(0,1);
      for k = 1:nrpt
        i1  = cat(1,i1,ones(cumtapcnt(k),1)*k);
        i2  = cat(1,i2,numel(i2)+(1:cumtapcnt(k))');
        vec = cat(1,vec,ones(cumtapcnt(k),1)./cumtapcnt(k));
      end
      P = sparse(i1,i2,vec);
      
      pow    = nan(nvox,nrpt);
      inside = find(data.inside);
      for k = inside(:)'
        pow(k,:) = P*(abs(data.mom{k}).^2)';
      end
      
      if sqrtflag, pow = sqrt(pow); end
      if demeanflag,
        mdat = nanmean(pow,2);
        pow  = pow - mdat(:,ones(1,nrpt)); % FIXME only works for 1 frequency
      end
      
      if ncmb == size(pow,1)
        data.powcov = pow * pow';
        data.powcovdimord = 'pos_pos';
        powindx = [];
      else
        data.powcov = [reshape(pow * pow(cmb,:)', [ncmb*nvox 1]); sum(pow.^2,2)];
        try,
          data = rmfield(data, 'pow');
          data = rmfield(data, 'powdimord');
        end
        % powindx = [nvox+(1:nvox) nvox+(1:nvox); cmb*ones(1,nvox) nvox+(1:nvox)]';
        powindx = [repmat(ncmb*nvox+(1:nvox)',[ncmb 1]) reshape(repmat(ncmb*nvox+cmb(:)', [nvox 1]),[nvox*ncmb 1]); ncmb*nvox+(1:nvox)' ncmb*nvox+(1:nvox)'];
        data.pos = [repmat(data.pos, [ncmb+1 1])];% FIXME come up with something reshape( repmat(data.pos(cmb,:),[nvox 1]);data.pos data.pos];
        data.inside = reshape(repmat(data.inside(:), [1 ncmb+1])+repmat(nvox*(0:ncmb), [nvox 1]), [nvox*(ncmb+1) 1]);
        if ~isempty(data.outside)
          data.outside = reshape(repmat(data.outside(:), [1 ncmb+1])+repmat(nvox*(0:ncmb), [nvox 1]), [nvox*(ncmb+1) 1]);
        end
        data.powcovdimord = 'pos';
        % data.dim(2) = size(data.pos,1);
      end
    
    elseif strcmp(inparam, 'mom') && strcmp(outparam, 'crsspctrm'),
      % get mom as rpttap_pos_freq matrix
      % FIXME this assumes only 1 freq bin
      sizmom = size(data.mom{data.inside(1)});
      
      if sizmom(1)==1,
        mom = zeros(size(data.pos,1), sizmom(2));
        mom(data.inside, :) = cat(1, data.mom{data.inside});
        
        if keeprpt,
          [nvox, nrpt]   = size(mom);
          data.crsspctrm = transpose([mom.*conj(mom(ones(1,nvox)*cmb,:));abs(mom).^2]);
          data = rmfield(data, 'mom');
          powindx = [nvox+(1:nvox) nvox+(1:nvox); cmb*ones(1,nvox) nvox+(1:nvox)]';
          
          data.pos = [data.pos repmat(data.pos(cmb,:),[nvox 1]);data.pos data.pos];
          data.inside = [data.inside(:); data.inside(:)+nvox];
          data.outside = [data.outside(:); data.outside(:)+nvox];
          data.crsspctrmdimord = 'rpttap_pos';
          
        elseif ncmb<size(mom,2)
          % do it computationally more efficient
          [nvox, nrpt] = size(mom);
          data.crsspctrm = reshape((mom*mom(cmb,:)')./nrpt, [nvox*ncmb 1]);
          tmppow = mean(abs(mom).^2,2);
          data.crsspctrm = cat(1, data.crsspctrm, tmppow);
          tmpindx1 = transpose(ncmb*nvox + ones(ncmb+1,1)*(1:nvox));
          tmpindx2 = repmat(tmpindx1(cmb(:),end), [1 nvox])';
          tmpindx3 = repmat(cmb(:), [1 nvox])'; % expressed in original voxel indices
          powindx  = [tmpindx1(:) [tmpindx2(:);tmpindx1(:,end)]];
          
          data.pos = [repmat(data.pos, [ncmb 1]) data.pos(tmpindx3(:),:); data.pos data.pos];
          data.inside = data.inside(:)*ones(1,ncmb+1) + (ones(length(data.inside),1)*nvox)*(0:ncmb);
          data.inside = data.inside(:);
          data.outside = setdiff((1:nvox*(ncmb+1))', data.inside);
          data = rmfield(data, 'mom');
          data.crsspctrmdimord = 'pos';
        else
          [nvox, nrpt] = size(mom);
          data.crsspctrm = (mom*mom')./nrpt;
          data = rmfield(data, 'mom');
          powindx = [];
          data.crsspctrmdimord = 'pos_pos_freq'; % FIXME hard coded
        end
        
        data.dimord = data.crsspctrmdimord;
        clear mom;
        
      elseif sizmom(1)>1
        % source moments are multivariate
        tmpindx = reshape(1:size(data.pos,1)*sizmom(2), [sizmom(2) size(data.pos,1)]);
        tmpinside = tmpindx(:, data.inside);
        tmpcmb = tmpindx(:, cmb);
        tmpncmb = numel(tmpcmb);
        mom = zeros(sizmom(1), sizmom(2)*size(data.pos,1));
        mom(:, tmpinside(:)) = cat(2, data.mom{data.inside});
        
        if keeprpt,
          error('keeprpt with multivariate dipole moments is not supported');
          % FIXME should this be supported
        elseif tmpncmb<size(mom,2)
          % do it computationally more efficient
          [nrpt,nvox] = size(mom);
          
          % linearly represent the voxel-pair csd matrices
          data.crsspctrm = reshape((transpose(mom)*conj(mom(:,tmpcmb)))./nrpt, [nvox*tmpncmb 1]);
          offsetauto = size(data.crsspctrm,1);
          
          % linearly represent the per voxel csd matrices
          tmppow = zeros(numel(tmpindx), size(tmpindx,1));
          for k = 1:size(tmpindx,2)
            tmppow(tmpindx(:,k), :) = (transpose(mom(:,tmpindx(:,k)))*conj(mom(:,tmpindx(:,k))))./nrpt;
          end
          tmppow = reshape(tmppow, [numel(tmppow) 1]);
          data.crsspctrm = cat(1, data.crsspctrm, tmppow);
          
          tmpindx1 = transpose(ncmb*nvox + ones(ncmb+1,1)*(1:nvox));
          tmpindx2 = repmat(tmpindx1(cmb(:),end), [1 nvox])';
          tmpindx3 = repmat(cmb(:), [1 nvox])'; % expressed in original voxel indices
          powindx = [tmpindx1(:) [tmpindx2(:);tmpindx1(:,end)]];
          
          data.pos = [repmat(data.pos, [ncmb 1]) data.pos(tmpindx3(:),:); data.pos data.pos];
          data.inside = data.inside(:)*ones(1,ncmb+1) + (ones(length(data.inside),1)*nvox)*(0:ncmb);
          data.inside = data.inside(:);
          data.outside = setdiff((1:nvox*(ncmb+1))', data.inside);
          if isfield(data, 'momdimord'),
            data.crsspctrmdimord = ['pos_',data.momdimord(14:end)];% FIXME this assumes dimord to be 'rpttap_...'
          end
          data = rmfield(data, 'mom');
          data = rmfield(data, 'momdimord');
        else
          [nrpt,nvox] = size(mom);
          data.crsspctrm = (transpose(mom)*conj(mom))./nrpt;
          data = rmfield(data, 'mom');
          try, data = rmfield(data, 'momdimord'); end
          powindx = [];
          data.crsspctrmdimord = 'pos_pos_freq'; % FIXME hard coded
        end
        data.dimord = data.crsspctrmdimord;
        clear mom;
        
      end % if sizmom(2)==1 or >1
      
    else
      error('unknown conversion from univariate to bivariate representation');
    end
    
  case 'raw'
    % construct a timelock-like structure that only contains the covariance, see ft_datatype_timelock
    timelock = [];

    if ~strcmp(inparam, 'trial')
      error('incorrect specification of inparam')
    elseif ~strcmp(outparam, 'cov'),
      error('incorrect specification of outparam')
    end
    
    nrpt  = length(data.trial);
    nchan = length(data.label);
    if nrpt==1
      % don't bother to try and keep repetitions
      keeprpt = false;
    end
    
    if keeprpt
        timelock.dimord = 'rpt_chan_chan';
        hasrpt = true;
    else
        timelock.dimord = 'chan_chan';
        hasrpt = false;
    end % if keeprpt

    if isfield(data, 'grad')
      timelock.grad = data.grad;
    end
    if isfield(data, 'elec')
      timelock.elec = data.elec;
    end
    timelock.label = data.label;
    if isfield(data, 'cfg'), timelock.cfg = data.cfg; end
    
    tmpcov   = zeros(nrpt, nchan, nchan);
    nsamples = zeros(1,nrpt);
    for i=1:nrpt
      nsamples(i) = length(data.time{i});
      tmpdat = data.trial{i};
      if demeanflag
        tmpdat = ft_preproc_baselinecorrect(tmpdat);
      end
      tmpcov(i,:,:) = (tmpdat * tmpdat')./(nsamples(i)-1); % use N-1 in the denominator
    end
    
    if ~keeprpt
      % average the covariance over trials
      % the following is not weighted by the number of samples
      %       timelock.cov = reshape(mean(tmpcov, 1), nchan, nchan);
      % the following is weighted by the number of samples
      %       for i=1:nrpt
      %         tmpcov(i,:,:) = tmpcov(i,:,:) * nsamples(i);
      %       end
      %       timelock.cov = reshape(sum(tmpcov,1), nchan, nchan) ./ sum(nsamples);
      % the following is an efficient implementation
      timelock.cov = reshape(nsamples * reshape(tmpcov, nrpt, nchan*nchan) ./ sum(nsamples), nchan, nchan);
    else
      timelock.cov = tmpcov;
    end
    
    if ~isequal(cmb, {'all' 'all'})
      % make a selection of channel combinations
      keyboard
    end
    
    % replace the input raw data with the timelock structure containing the covariance
    data = timelock;
    
  otherwise
    error('unsupported input data type');
end % swith dtype

if ~exist('hasrpt', 'var')
  % in the case of raw data it has already been assigned
  hasrpt = (isfield(data,            'dimord')  && ~isempty(strfind(data.dimord,                'rpt'))) || ...
           (isfield(data, [outparam, 'dimord']) && ~isempty(strfind(data.([outparam,'dimord']), 'rpt')));
end
       
       
% % ----------------------------------------
% function [indx] = labelcmb2indx(labelcmb)
%
% % identify the auto-combinations
% ncmb = size(labelcmb,1);
% indx = zeros(ncmb,2);
%
% label = unique(labelcmb(:));
% nchan = numel(label);
% autoindx = zeros(nchan,1);
% for k = 1:nchan
% sel1 = strcmp(label{k}, labelcmb(:,1));
% sel2 = strcmp(label{k}, labelcmb(:,2));
% autoindx = find(sel1 & sel2);
%
% indx(sel1,1) = autoindx;
% indx(sel2,2) = autoindx;
% end
