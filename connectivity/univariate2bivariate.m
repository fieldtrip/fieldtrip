function [data, powindx, hasrpt] = univariate2bivariate(data, inparam, outparam, dtype, varargin)

%if nargin<8, keeprpt    = 1;  end
%if nargin<7, sqrtflag   = 0;  end
%if nargin<6, cmb        = []; end
%if nargin<5, demeanflag = 0;  end

demeanflag = keyval('demeanflag', varargin{:}); if isempty(demeanflag), demeanflag = 0; end
cmb        = keyval('cmb',     varargin{:});
sqrtflag   = keyval('sqrtflag', varargin{:}); if isempty(sqrtflag), sqrtflag = 0;       end
keeprpt    = keyval('keeprt',  varargin{:});  if isempty(keeprpt),  keeprpt  = 1;       end

switch dtype
  case 'freq'
    ncmb  = size(cmb,1);
    nchan = numel(data.label);
    getpowindx = 0;
    if ncmb==0,
      error('no channel combinations are specified');
    elseif ncmb==nchan^2 || ncmb==(nchan+1)*nchan*0.5,
      dofull = 1;
    else
      dofull = 0;
    end
    
    if strcmp(inparam, 'fourierspctrm') && strcmp(outparam, 'crsspctrm'),
      %fourier coefficients -> cross-spectral density
      if dofull
        data    = checkdata(data, 'cmbrepresentation', 'full');
      else
        data    = checkdata(data, 'cmbrepresentation', 'sparse', 'channelcmb', cmb);
      end
    elseif strcmp(inparam, 'powandcsd') && strcmp(outparam, 'crsspctrm'),
      if ~isempty(cmb),
        data    = checkdata(data, 'cmbrepresentation', 'sparse', 'channelcmb', cmb);
        
        %ensure getting powindx later on to prevent crash
        getpowindx = 1;
      else
        %data    = checkdata(data, 'cmbrepresentation', 'full');
        %this should not be possible
        error('cannot convert to a full csd representation');
        
        
      end
    elseif strcmp(inparam, 'fourierspctrm') && strcmp(outparam, 'powcovspctrm'),
      %fourier coefficients -> power covariance
      data = checkdata(data, 'cmbrepresentation', 'sparsewithpow', 'channelcmb', {});
      
      if sqrtflag, data.powspctrm = sqrt(data.powspctrm); end
      
      %get covariance by using checkdata
      if demeanflag,
        nrpt = size(data.powspctrm,1);
        mdat = nanmean(data.powspctrm,1);
        data.powspctrm = data.powspctrm - mdat(ones(1,nrpt),:,:,:,:,:);
      end
      data.fourierspctrm = data.powspctrm; %this is necessary for checkdata to work
      data.dimord        = ['rpttap',data.dimord(4:end)];
      data               = rmfield(data, 'powspctrm');
      data.cumtapcnt(:)  = 1;
      data.cumsumcnt(:)  = 1;
      if ncmb < (nchan-1)*nchan*0.5,
        data    = checkdata(data, 'cmbrepresentation', 'sparse', 'channelcmb', cmb);
      else
        data    = checkdata(data, 'cmbrepresentation', 'full');
      end
      data.powcovspctrm = data.crsspctrm;
      data              = rmfield(data, 'crsspctrm');
    elseif strcmp(inparam, 'powspctrm') && strcmp(outparam, 'powcovspctrm'),
      %power-spectral density -> power covariance
      
      if sqrtflag, data.powspctrm = sqrt(data.powspctrm); end
      
      %get covariance by using checkdata
      if demeanflag,
        nrpt = size(data.powspctrm,1);
        mdat = nanmean(data.powspctrm,1);
        data.powspctrm = data.powspctrm - mdat(ones(1,nrpt),:,:,:,:,:);
      end
      data.fourierspctrm = data.powspctrm; %this is necessary for checkdata to work
      data.dimord        = ['rpttap',data.dimord(4:end)];
      data               = rmfield(data, 'powspctrm');
      data.cumtapcnt(:)  = 1;
      data.cumsumcnt(:)  = 1;
      if ncmb < (nchan-1)*nchan*0.5,
        data    = checkdata(data, 'cmbrepresentation', 'sparse', 'channelcmb', cmb);
      else
        data    = checkdata(data, 'cmbrepresentation', 'full');
      end
      data.powcovspctrm = data.crsspctrm;
      data = rmfield(data, 'crsspctrm');
    else
      error('unknown conversion from univariate to bivariate representation');
    end
    
    if ~isempty(cmb) && (ncmb < (nchan-1)*nchan*0.5 || getpowindx==1),
      powindx = labelcmb2indx(data.labelcmb);
    else
      powindx = [];
    end
  case 'source'
    ncmb = numel(cmb);
    
    if strcmp(inparam, 'pow') && strcmp(outparam, 'powcov'),
      [nrpt,nvox] = size(data.pow);
      if sqrtflag, data.pow = sqrt(data.pow); end
      if demeanflag,
        mdat = nanmean(data.pow,1);
        data.pow = data.pow - mdat(ones(1,nrpt),:); %FIXME only works for 1 frequency
      end
      
      data.powcov = [data.pow .* data.pow(:,ones(1,nvox)*cmb) data.pow.*data.pow];
      data        = rmfield(data, 'pow');
      powindx     = [nvox+(1:nvox) nvox+(1:nvox); cmb*ones(1,nvox) nvox+(1:nvox)]';
      
      data.pos    = [data.pos repmat(data.pos(cmb,:),[nvox 1]);data.pos data.pos];
      data.inside = [data.inside(:); data.inside(:)+nvox];
      data.outside = [data.outside(:); data.outside(:)+nvox];
      data.dim(2) = size(data.pos,1);
    elseif strcmp(inparam, 'mom') && strcmp(outparam, 'crsspctrm'),
      %get mom as rpttap_pos_freq matrix
      %FIXME this assumes only 1 freq bin
      mom = zeros(size(data.mom{data.inside(1)},1), size(data.pos,1));
      mom(:, data.inside) = cat(2, data.mom{data.inside});
      if keeprpt,
        [nrpt,nvox]    = size(mom);
        data.crsspctrm = [mom.*conj(mom(:,ones(1,nvox)*cmb)) abs(mom).^2];
        data           = rmfield(data, 'mom');
        data           = rmfield(data, 'momdimord');
        powindx     = [nvox+(1:nvox) nvox+(1:nvox); cmb*ones(1,nvox) nvox+(1:nvox)]';
        
        data.pos    = [data.pos repmat(data.pos(cmb,:),[nvox 1]);data.pos data.pos];
        data.inside = [data.inside(:); data.inside(:)+nvox];
        data.outside = [data.outside(:); data.outside(:)+nvox];
      elseif ncmb<size(mom,2)
        %do it computationally more efficient
        [nrpt,nvox]    = size(mom);
        
        data.crsspctrm = reshape((transpose(mom)*conj(mom(:,cmb)))./nrpt, [nvox*ncmb 1]);
        tmppow         = mean(abs(mom).^2)';
        data.crsspctrm = cat(1, data.crsspctrm, tmppow);
        tmpindx1       = transpose(ncmb*nvox + ones(ncmb+1,1)*(1:nvox));
        tmpindx2       = repmat(tmpindx1(cmb(:),end), [1 nvox])';
        tmpindx3       = repmat(cmb(:), [1 nvox])'; %expressed in original voxel indices
        powindx        = [tmpindx1(:) [tmpindx2(:);tmpindx1(:,end)]];
        
        data.pos       = [repmat(data.pos, [ncmb 1]) data.pos(tmpindx3(:),:); data.pos data.pos];
        data.inside    = data.inside(:)*ones(1,ncmb+1) + (ones(length(data.inside),1)*nvox)*(0:ncmb);
        data.inside    = data.inside(:);
        data.outside   = setdiff((1:nvox*(ncmb+1))', data.inside);
        if isfield(data, 'momdimord'),
          data.crsspctrmdimord = ['pos_',data.momdimord(14:end)];%FIXME this assumes dimord to be 'rpttap_...'
        end
        data           = rmfield(data, 'mom');
        data           = rmfield(data, 'momdimord');
      else
        [nrpt,nvox]    = size(mom);
        data.crsspctrm = (transpose(mom)*conj(mom))./nrpt;
        data           = rmfield(data, 'mom');
        data           = rmfield(data, 'momdimord');
        powindx        = [];
        data.crsspctrmdimord = 'pos_pos_freq'; %FIXME hard coded
      end
      data.dimord = data.crsspctrmdimord;
      clear mom;
    else
      error('unknown conversion from univariate to bivariate representation');
    end
  otherwise
end

hasrpt  = (isfield(data, 'dimord') && ~isempty(strfind(data.dimord, 'rpt')));

%----------------------------------------
function [indx] = labelcmb2indx(labelcmb)

%identify the auto-combinations
ncmb = size(labelcmb,1);
indx = zeros(ncmb,2);

label = unique(labelcmb(:));
nchan = numel(label);
autoindx = zeros(nchan,1);
for k = 1:nchan
  sel1 = strcmp(label{k}, labelcmb(:,1));
  sel2 = strcmp(label{k}, labelcmb(:,2));
  autoindx = find(sel1 & sel2);
  
  indx(sel1,1) = autoindx;
  indx(sel2,2) = autoindx;
end
