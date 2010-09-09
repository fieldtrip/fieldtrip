function [data] = fixcsd(data, desired, channelcmb)

% FIXCSD converts univariate frequency domain data (fourierspctrm) into a bivariate
% representation (crsspctrm), or changes the representation of bivariate frequency
% domain data (sparse/full/sparsewithpow, sparsewithpow only works for crsspctrm or
% fourierspctrm)

% Copyright (C) 2010, Jan-Mathijs Schoffelen, Robert Oostenveld 

if isfield(data, 'crsspctrm') && isfield(data, 'powspctrm')
  current = 'sparsewithpow';
elseif isfield(data, 'powspctrm')
  current = 'sparsewithpow';
elseif isfield(data, 'fourierspctrm') && ~isfield(data, 'labelcmb')
  current = 'fourier';
elseif ~isfield(data, 'labelcmb')
  current = 'full';
elseif isfield(data, 'labelcmb')
  current = 'sparse';
else
  error('Could not determine the current representation of the %s matrix', param);
end

% first go from univariate fourier to the required bivariate representation
if strcmp(current, 'fourier') && strcmp(desired, 'fourier')
  % nothing to do
elseif strcmp(current, 'fourier') && strcmp(desired, 'sparsewithpow')
  dimtok = tokenize(data.dimord, '_');
  if ~isempty(strmatch('rpttap',   dimtok)),
    nrpt = length(data.cumtapcnt);
    flag = 0;
  else
    nrpt = 1;
    flag = 1;
  end
  if ~isempty(strmatch('freq',  dimtok)), nfrq=length(data.freq);      else nfrq = 1; end
  if ~isempty(strmatch('time',  dimtok)), ntim=length(data.time);      else ntim = 1; end
  
  fastflag = all(data.cumtapcnt(:)==data.cumtapcnt(1));

  %create auto-spectra
  nchan     = length(data.label);
  if fastflag
    % all trials have the same amount of tapers
    powspctrm = zeros(nrpt,nchan,nfrq,ntim);
    ntap      = data.cumtapcnt(1);
    for p = 1:ntap
      powspctrm = powspctrm + abs(data.fourierspctrm(p:ntap:end,:,:,:,:)).^2;
    end
    powspctrm = powspctrm./ntap;
  else
    % different amount of tapers
    powspctrm = zeros(nrpt,nchan,nfrq,ntim)+i.*zeros(nrpt,nchan,nfrq,ntim);
    sumtapcnt = [0;cumsum(data.cumtapcnt(:))];
    for p = 1:nrpt
      indx   = (sumtapcnt(p)+1):sumtapcnt(p+1);
      tmpdat = data.fourierspctrm(indx,:,:,:);
      powspctrm(p,:,:,:) = (sum(tmpdat.*conj(tmpdat),1))./data.cumtapcnt(p);
    end
  end

  %create cross-spectra
  if ~isempty(channelcmb),
    ncmb      = size(channelcmb,1);
    cmbindx   = zeros(ncmb,2);
    labelcmb  = cell(ncmb,2);
    for k = 1:ncmb
      ch1 = find(strcmp(data.label, channelcmb(k,1)));
      ch2 = find(strcmp(data.label, channelcmb(k,2)));
      if ~isempty(ch1) && ~isempty(ch2),
        cmbindx(k,:)  = [ch1 ch2];
        labelcmb(k,:) = data.label([ch1 ch2])';
      end
    end

    crsspctrm = zeros(nrpt,ncmb,nfrq,ntim)+i.*zeros(nrpt,ncmb,nfrq,ntim);
    if fastflag
      for p = 1:ntap
        tmpdat1   = data.fourierspctrm(p:ntap:end,cmbindx(:,1),:,:,:);
        tmpdat2   = data.fourierspctrm(p:ntap:end,cmbindx(:,2),:,:,:);
        crsspctrm = crsspctrm + tmpdat1.*conj(tmpdat2);
      end
      crsspctrm = crsspctrm./ntap;
    else
      for p = 1:nrpt
        indx    = (sumtapcnt(p)+1):sumtapcnt(p+1);
        tmpdat1 = data.fourierspctrm(indx,cmbindx(:,1),:,:);
        tmpdat2 = data.fourierspctrm(indx,cmbindx(:,2),:,:);
        crsspctrm(p,:,:,:) = (sum(tmpdat1.*conj(tmpdat2),1))./data.cumtapcnt(p);
      end
    end
    data.crsspctrm = crsspctrm;
    data.labelcmb  = labelcmb;
  end
  data.powspctrm = powspctrm;
  data           = rmfield(data, 'fourierspctrm');
  if ntim>1,
    data.dimord = 'chan_freq_time';
  else
    data.dimord = 'chan_freq';
  end
  
  if nrpt>1,
    data.dimord = ['rpt_',data.dimord];
  end
  
  if flag, siz = size(data.crsspctrm); data.crsspctrm = reshape(data.crsspctrm, siz(2:end)); end
elseif strcmp(current, 'fourier') && strcmp(desired, 'sparse')

  if isempty(channelcmb), error('no channel combinations are specified'); end
  dimtok = tokenize(data.dimord, '_');
  if ~isempty(strmatch('rpttap',   dimtok)),
    nrpt = length(data.cumtapcnt);
    flag = 0;
  else
    nrpt = 1;
    flag = 1;
  end
  if ~isempty(strmatch('freq',  dimtok)), nfrq=length(data.freq); else nfrq = 1; end
  if ~isempty(strmatch('time',  dimtok)), ntim=length(data.time); else ntim = 1; end
  
  ncmb      = size(channelcmb,1);
  cmbindx   = zeros(ncmb,2);
  labelcmb  = cell(ncmb,2);
  for k = 1:ncmb
    ch1 = find(strcmp(data.label, channelcmb(k,1)));
    ch2 = find(strcmp(data.label, channelcmb(k,2)));
    if ~isempty(ch1) && ~isempty(ch2),
      cmbindx(k,:)  = [ch1 ch2];
      labelcmb(k,:) = data.label([ch1 ch2])';
    end
  end

  sumtapcnt = [0;cumsum(data.cumtapcnt(:))];
  fastflag  = all(data.cumtapcnt(:)==data.cumtapcnt(1));
  if fastflag && nrpt>1
    ntap = data.cumtapcnt(1);
    
    % compute running sum across tapers
    for p = 1:ntap
      indx      = p:ntap:nrpt*ntap;
      if p==1
        crsspctrm = data.fourierspctrm(indx,cmbindx(:,1),:,:).*  ...
               conj(data.fourierspctrm(indx,cmbindx(:,2),:,:));
      else
        crsspctrm = data.fourierspctrm(indx,cmbindx(:,1),:,:).*  ...
               conj(data.fourierspctrm(indx,cmbindx(:,2),:,:)) + crsspctrm;
      end
    end
    crsspctrm = crsspctrm./ntap;
  else
    for p = 1:nrpt
      indx    = (sumtapcnt(p)+1):sumtapcnt(p+1);
      tmpdat1 = data.fourierspctrm(indx,cmbindx(:,1),:,:);
      tmpdat2 = data.fourierspctrm(indx,cmbindx(:,2),:,:);
      crsspctrm(p,:,:,:) = (sum(tmpdat1.*conj(tmpdat2),1))./data.cumtapcnt(p);
    end
  end
  data.crsspctrm = crsspctrm;
  data.labelcmb  = labelcmb;
  data           = rmfield(data, 'fourierspctrm');
  data           = rmfield(data, 'label');
  if ntim>1,
    data.dimord = 'chan_freq_time';
  else
    data.dimord = 'chan_freq';
  end
  
  if nrpt>1,
    data.dimord = ['rpt_',data.dimord];
  end

  if flag, siz = size(data.crsspctrm); data.crsspctrm = reshape(data.crsspctrm, siz(2:end)); end
elseif strcmp(current, 'fourier') && strcmp(desired, 'full')

  % this is how it is currently and the desired functionality of prepare_freq_matrices
  dimtok = tokenize(data.dimord, '_');
  if ~isempty(strmatch('rpttap',   dimtok)),
    nrpt = length(data.cumtapcnt);
    flag = 0;
  else
    nrpt = 1;
    flag = 1;
  end
  if ~isempty(strmatch('rpttap',dimtok)), nrpt=length(data.cumtapcnt); else nrpt = 1; end
  if ~isempty(strmatch('freq',  dimtok)), nfrq=length(data.freq);      else nfrq = 1; end
  if ~isempty(strmatch('time',  dimtok)), ntim=length(data.time);      else ntim = 1; end
  nchan     = length(data.label);
  crsspctrm = zeros(nrpt,nchan,nchan,nfrq,ntim);
  sumtapcnt = [0;cumsum(data.cumtapcnt(:))];
  for k = 1:ntim
    for m = 1:nfrq
      for p = 1:nrpt
        %FIXME speed this up in the case that all trials have equal number of tapers
        indx   = (sumtapcnt(p)+1):sumtapcnt(p+1);
        tmpdat = transpose(data.fourierspctrm(indx,:,m,k));
        crsspctrm(p,:,:,m,k) = (tmpdat*tmpdat')./data.cumtapcnt(p);
        clear tmpdat;
      end
    end
  end
  data.crsspctrm = crsspctrm;
  data           = rmfield(data, 'fourierspctrm');
  
  if ntim>1,
    data.dimord = 'chan_chan_freq_time';
  else
    data.dimord = 'chan_chan_freq';
  end
  
  if nrpt>1,
    data.dimord = ['rpt_',data.dimord];
  end   
  
  % remove first singleton dimension
  if flag, siz = size(data.crsspctrm); data.crsspctrm = reshape(data.crsspctrm, siz(2:end)); end

elseif strcmp(current, 'fourier') && strcmp(desired, 'fullfast'),

  dimtok = tokenize(data.dimord, '_');
  nrpt = size(data.fourierspctrm, 1);    
  nchn = numel(data.label);    
  nfrq = numel(data.freq);  
  if ~isempty(strmatch('time',  dimtok)), ntim=numel(data.time); else ntim = 1; end
  
  data.fourierspctrm = reshape(data.fourierspctrm, [nrpt nchn nfrq*ntim]);
  data.fourierspctrm(~isfinite(data.fourierspctrm)) = 0;
  crsspctrm = complex(zeros(nchn,nchn,nfrq*ntim));
  for k = 1:nfrq*ntim
    tmp = transpose(data.fourierspctrm(:,:,k));
    n   = sum(tmp~=0,2);
    crsspctrm(:,:,k) = tmp*tmp'./n(1);
  end
  data           = rmfield(data, 'fourierspctrm');
  data.crsspctrm = reshape(crsspctrm, [nchn nchn nfrq ntim]);
  if isfield(data, 'time'),
    data.dimord = 'chan_chan_freq_time';
  else
    data.dimord = 'chan_chan_freq';
  end

end % convert to the requested bivariate representation

% from one bivariate representation to another
if isequal(current, desired)
  % nothing to do

elseif (strcmp(current, 'full')       && strcmp(desired, 'fourier')) || ...
    (strcmp(current, 'sparse')        && strcmp(desired, 'fourier')) || ...
    (strcmp(current, 'sparsewithpow') && strcmp(desired, 'fourier'))
  % this is not possible
  error('converting the cross-spectrum into a Fourier representation is not possible');

elseif strcmp(current, 'full') && strcmp(desired, 'sparsewithpow')
  error('not yet implemented');
elseif strcmp(current, 'sparse') && strcmp(desired, 'sparsewithpow')
  % convert back to crsspctrm/powspctrm representation: useful for plotting functions etc
  indx     = labelcmb2indx(data.labelcmb);
  autoindx = indx(indx(:,1)==indx(:,2), 1);
  cmbindx  = setdiff([1:size(indx,1)]', autoindx);
  
  if strcmp(data.dimord(1:3), 'rpt')
    data.powspctrm = data.crsspctrm(:, autoindx, :, :);
    data.crsspctrm = data.crsspctrm(:, cmbindx,  :, :);
  else
    data.powspctrm = data.crsspctrm(autoindx, :, :);
    data.crsspctrm = data.crsspctrm(cmbindx,  :, :);
  end 
  data.label    = data.labelcmb(autoindx,1);
  data.labelcmb = data.labelcmb(cmbindx, :);
  
  if isempty(cmbindx)
    data = rmfield(data, 'crsspctrm');
    data = rmfield(data, 'labelcmb');
  end
  
elseif strcmp(current, 'full') && strcmp(desired, 'sparse')
  dimtok = tokenize(data.dimord, '_');
  if ~isempty(strmatch('rpt',   dimtok)), nrpt=numel(data.cumtapcnt); else nrpt = 1; end
  if ~isempty(strmatch('freq',  dimtok)), nfrq=numel(data.freq);      else nfrq = 1; end
  if ~isempty(strmatch('time',  dimtok)), ntim=numel(data.time);      else ntim = 1; end
  nchan    = length(data.label);
  ncmb     = nchan*nchan;
  labelcmb = cell(ncmb, 2);
  cmbindx  = zeros(nchan, nchan);
  k = 1;
  for j=1:nchan
    for m=1:nchan
      labelcmb{k, 1} = data.label{m};
      labelcmb{k, 2} = data.label{j};
      cmbindx(m,j)   = k;
      k = k+1;
    end
  end
  
  % reshape all possible fields
  fn = fieldnames(data);
  for ii=1:numel(fn)
    if numel(data.(fn{ii})) == nrpt*ncmb*nfrq*ntim;
      if nrpt>1,
        data.(fn{ii}) = reshape(data.(fn{ii}), nrpt, ncmb, nfrq, ntim);
      else
        data.(fn{ii}) = reshape(data.(fn{ii}), ncmb, nfrq, ntim);
      end
    end
  end
  % remove obsolete fields
  data           = rmfield(data, 'label');
  try, data      = rmfield(data, 'dof'); end
  % replace updated fields
  data.labelcmb  = labelcmb;
  if ntim>1,
    data.dimord = 'chancmb_freq_time';
  else
    data.dimord = 'chancmb_freq';
  end

  if nrpt>1,
    data.dimord = ['rpt_',data.dimord];
  end

elseif strcmp(current, 'sparsewithpow') && strcmp(desired, 'sparse')

  % this representation for sparse data contains autospectra
  % as e.g. {'A' 'A'} in labelcmb
  if isfield(data, 'crsspctrm'),
    dimtok         = tokenize(data.dimord, '_');
    catdim         = match_str(dimtok, {'chan' 'chancmb'});
    data.crsspctrm = cat(catdim, data.powspctrm, data.crsspctrm);
    data.labelcmb  = [data.label(:) data.label(:); data.labelcmb];
    data           = rmfield(data, 'powspctrm');
  else
    data.crsspctrm = data.powspctrm;
    data.labelcmb  = [data.label(:) data.label(:)];
    data           = rmfield(data, 'powspctrm');
  end
  data = rmfield(data, 'label');

elseif strcmp(current, 'sparse') && strcmp(desired, 'full')
  dimtok = tokenize(data.dimord, '_');
  if ~isempty(strmatch('rpt',   dimtok)), nrpt=numel(data.cumtapcnt); else nrpt = 1; end
  if ~isempty(strmatch('freq',  dimtok)), nfrq=numel(data.freq);      else nfrq = 1; end
  if ~isempty(strmatch('time',  dimtok)), ntim=numel(data.time);      else ntim = 1; end
  
  if ~isfield(data, 'label')
    data.label = unique(data.labelcmb(:));
  end

  nchan     = length(data.label);
  ncmb      = size(data.labelcmb,1);
  cmbindx   = zeros(nchan,nchan);

  for k = 1:size(data.labelcmb,1)
    ch1 = find(strcmp(data.label, data.labelcmb(k,1)));
    ch2 = find(strcmp(data.label, data.labelcmb(k,2)));
    if ~isempty(ch1) && ~isempty(ch2),
      cmbindx(ch1,ch2) = k;
    end
  end

  complete = all(cmbindx(:)~=0);

  fn = fieldnames(data);
  for ii=1:numel(fn)
    if numel(data.(fn{ii})) == nrpt*ncmb*nfrq*ntim;
      if nrpt==1,
        data.(fn{ii}) = reshape(data.(fn{ii}), [nrpt ncmb nfrq ntim]);
      end

      tmpall = nan(nrpt,nchan,nchan,nfrq,ntim);

      for j = 1:nrpt
        for k = 1:ntim
          for m = 1:nfrq
            tmpdat = nan(nchan,nchan);
            indx   = find(cmbindx);
            if ~complete
              % this realizes the missing combinations to be represented as the
              % conjugate of the corresponding combination across the diagonal
              tmpdat(indx) = reshape(data.(fn{ii})(j,cmbindx(indx),m,k),[numel(indx) 1]);
              tmpdat       = ctranspose(tmpdat);
            end
            tmpdat(indx)    = reshape(data.(fn{ii})(j,cmbindx(indx),m,k),[numel(indx) 1]);
            tmpall(j,:,:,m,k) = tmpdat;
          end % for m
        end % for k
      end % for j

      % replace the data in the old representation with the new representation
      if nrpt>1,
        data.(fn{ii}) = tmpall;
      else
        data.(fn{ii}) = reshape(tmpall, [nchan nchan nfrq ntim]);
      end
    end % if numel
  end % for ii

  % remove obsolete fields
  try, data      = rmfield(data, 'powspctrm');  end
  try, data      = rmfield(data, 'labelcmb');   end
  try, data      = rmfield(data, 'dof');        end

  if ntim>1,
    data.dimord = 'chan_chan_freq_time';
  else
    data.dimord = 'chan_chan_freq';
  end

  if nrpt>1,
    data.dimord = ['rpt_',data.dimord];
  end

elseif strcmp(current, 'sparse') && strcmp(desired, 'fullfast')
  dimtok = tokenize(data.dimord, '_');
  if ~isempty(strmatch('rpt',   dimtok)), nrpt=numel(data.cumtapcnt); else nrpt = 1; end
  if ~isempty(strmatch('freq',  dimtok)), nfrq=numel(data.freq);      else nfrq = 1; end
  if ~isempty(strmatch('time',  dimtok)), ntim=numel(data.time);      else ntim = 1; end
  
  if ~isfield(data, 'label')
    data.label = unique(data.labelcmb(:));
  end

  nchan     = length(data.label);
  ncmb      = size(data.labelcmb,1);
  cmbindx   = zeros(nchan,nchan);

  for k = 1:size(data.labelcmb,1)
    ch1 = find(strcmp(data.label, data.labelcmb(k,1)));
    ch2 = find(strcmp(data.label, data.labelcmb(k,2)));
    if ~isempty(ch1) && ~isempty(ch2),
      cmbindx(ch1,ch2) = k;
    end
  end

  complete = all(cmbindx(:)~=0);

  fn = fieldnames(data);
  for ii=1:numel(fn)
    if numel(data.(fn{ii})) == nrpt*ncmb*nfrq*ntim;
      if nrpt==1,
        data.(fn{ii}) = reshape(data.(fn{ii}), [nrpt ncmb nfrq ntim]);
      end

      tmpall = nan(nchan,nchan,nfrq,ntim);

      for k = 1:ntim
        for m = 1:nfrq
          tmpdat = nan(nchan,nchan);
          indx   = find(cmbindx);
          if ~complete
            % this realizes the missing combinations to be represented as the
            % conjugate of the corresponding combination across the diagonal
            tmpdat(indx) = reshape(nanmean(data.(fn{ii})(:,cmbindx(indx),m,k)),[numel(indx) 1]);
            tmpdat       = ctranspose(tmpdat);
          end
          tmpdat(indx)    = reshape(nanmean(data.(fn{ii})(:,cmbindx(indx),m,k)),[numel(indx) 1]);
          tmpall(:,:,m,k) = tmpdat;
        end % for m
      end % for k

      % replace the data in the old representation with the new representation
      if nrpt>1,
        data.(fn{ii}) = tmpall;
      else
        data.(fn{ii}) = reshape(tmpall, [nchan nchan nfrq ntim]);
      end
    end % if numel
  end % for ii

  % remove obsolete fields
  try, data      = rmfield(data, 'powspctrm');  end
  try, data      = rmfield(data, 'labelcmb');   end
  try, data      = rmfield(data, 'dof');        end

  if ntim>1,
    data.dimord = 'chan_chan_freq_time';
  else
    data.dimord = 'chan_chan_freq';
  end

elseif strcmp(current, 'sparsewithpow') && strcmp(desired, 'full')
  % this is how is currently done in prepare_freq_matrices
  data = checkdata(data, 'cmbrepresentation', 'sparse');
  data = checkdata(data, 'cmbrepresentation', 'full');

end
