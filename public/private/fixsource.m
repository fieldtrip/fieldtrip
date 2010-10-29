function [output] = fixsource(input, varargin)

% FIXSOURCE converts old style source structures into new style source structures and the
% other way around
%
% Use as:
%   fixsource(input, type)
%    where input is a source structure,
%
% Typically, old style source structures contain
%   avg.XXX or trial.XXX fields
%
% The new style source structure contains:
%   source.pos
%   source.dim (optional, if the list of positions describes a 3D volume
%   source.XXX the old style subfields in avg/trial
%   source.XXXdimord string how to interpret the respective XXX field:
%     e.g. source.leadfield = cell(1,Npos), source.leadfielddimord = '{pos}_chan_ori'
%          source.mom       = cell(1,Npos), source.momdimord       = '{pos}_ori_rpttap'

type   = keyval('type',   varargin{:});
haspow = keyval('haspow', varargin{:});

if isempty(type),   type   = 'old'; end
if isempty(haspow), haspow = 'no';  end

fnames = fieldnames(input);
tmp    = cell2mat(strfind(fnames, 'dimord')); %get dimord like fields
if any(tmp>1),
  current = 'new';
elseif any(tmp==1),
  %don't know what to do yet data is JM's own invention
  current = 'old';
else
  current = 'old';
end

if strcmp(current, type),
  %do nothing
  output = input;
  
  %return
elseif strcmp(current, 'old') && strcmp(type, 'new'),
  %go from old to new

  if isfield(input, 'avg'),
    stuff  = getfield(input, 'avg');
    output = rmfield(input,  'avg');
  elseif isfield(input, 'trial'),
    stuff  = getfield(input, 'trial');
    output = rmfield(input,  'trial');
  else
    %this could occur later in the pipeline, e.g. when doing group statistics using individual subject
    %descriptive statistics
    error('the input does not contain an avg or trial field');
  end

  %-------------------------------------------------  
  %remove and rename the specified fields if present
  removefields = {'xgrid';'ygrid';'zgrid';'method'};
  renamefields = {'frequency' 'freq'; 'csdlabel' 'orilabel'};
  fnames       = fieldnames(output);
  for k = 1:numel(fnames)
    ix = strmatch(fnames{k}, removefields);
    if ~isempty(ix),
      output = rmfield(output, fnames{k});
    end
    ix = strmatch(fnames{k}, renamefields(:,1), 'exact');
    if ~isempty(ix),
      output = setfield(output, renamefields{ix,2}, ...
                        getfield(output, renamefields{ix,1}));
      output = rmfield(output, fnames{k});
    end
  end

  %----------------------------------------------------------------------
  %put the stuff originally in avg or trial one level up in the structure
  fnames       = fieldnames(stuff(1));
  npos         = size(input.pos,1);
  nrpt         = numel(stuff);
  for k = 1:numel(fnames)
    if nrpt>1,
      %multiple trials
      %(or subjects FIXME not yet implemented, nor tested)
      tmp  = getfield(stuff(1), fnames{k});
      siz  = size(tmp);
      if isfield(input, 'cumtapcnt') && strcmp(fnames{k}, 'mom')
        %pcc based mom is orixrpttap
        %tranpose to keep manageable
        for kk = 1:numel(input.inside)
          indx = input.inside(kk);
          tmp{indx} = permute(tmp{indx}, [2 1 3]); 
        end
        nrpttap = sum(input.cumtapcnt);
        sizvox  = [size(tmp{input.inside(1)}) 1];
        sizvox  = [nrpttap sizvox(2:end)];
      elseif strcmp(fnames{k}, 'mom'),
        %this is then probably not a frequency based mom
        nrpttap = numel(stuff);
        sizvox  = [size(tmp{input.inside(1)}) 1];
        sizvox  = [nrpttap sizvox];
      elseif iscell(tmp)
        nrpttap = numel(stuff);
        sizvox  = [size(tmp{input.inside(1)}) 1];
        sizvox  = [nrpttap sizvox];
      end
      
      if siz(1) ~= npos && siz(2) ==npos,
        tmp = transpose(tmp);
      end
      
      if iscell(tmp)
        %allocate memory for cell-array
        tmpall = cell(npos,1);
        for n = 1:numel(input.inside)
          tmpall{input.inside(n)} = zeros(sizvox);
        end
      else
        %allocate memory for matrix
        tmpall = zeros([npos nrpt siz(2:end)]);
      end
      
      cnt = 0;
      for m = 1:nrpt
        tmp = getfield(stuff(m), fnames{k});
        siz = size(tmp);
        if siz(1) ~= npos && siz(2) ==npos,
          tmp = transpose(tmp);
        end
        
        if ~iscell(tmp),
          tmpall(:,m,:,:,:) = tmp;
        else
          for n = 1:numel(input.inside)
            indx   = input.inside(n);
            tmpdat = tmp{indx};
            if isfield(input, 'cumtapcnt') && strcmp(fnames{k}, 'mom'),
              if n==1, siz1 = size(tmpdat,2); end
            else
              if n==1, siz1 = 1; end
            end
            tmpall{indx}(cnt+[1:siz1],:,:,:,:) = tmpdat;
            if n==numel(input.inside), cnt  = cnt + siz1;     end
          end
        end
      end
      output    = setfield(output, fnames{k}, tmpall);
      newdimord = createdimord(output, fnames{k}, 1);
      if ~isempty(newdimord)
        output    = setfield(output, [fnames{k},'dimord'], newdimord);
      end
    
    else
      tmp = getfield(stuff, fnames{k});
      siz = size(tmp);
      if isfield(input, 'cumtapcnt') && strcmp(fnames{k}, 'mom')
        %pcc based mom is orixrpttap
        %tranpose to keep manageable
        for kk = 1:numel(input.inside)
          indx = input.inside(kk);
          tmp{indx} = permute(tmp{indx}, [2 1 3]); 
        end
      end
      if siz(1) ~= npos && siz(2) ==npos,
        tmp = transpose(tmp);
      end
      output    = setfield(output, fnames{k}, tmp);
      newdimord = createdimord(output, fnames{k}); 
      if ~isempty(newdimord)
        output    = setfield(output, [fnames{k},'dimord'], newdimord);
      end
    end
  end
  
  if isfield(output, 'csdlabel')
    output = setfield(output, 'orilabel', getfield(output, 'csdlabel'));
    output = rmfield(output,  'csdlabel');
  end

  if isfield(output, 'leadfield')
    % add dimord to leadfield as well. since the leadfield is not in
    % the original .avg or .trial field it has not yet been taken care of
    output.leadfielddimord = createdimord(output, 'leadfield');  
  end
  
  if isfield(output, 'ori')
    % convert cell-array ori into matrix
    ori = zeros(3,npos) + nan;
    try,
      ori(:,output.inside) = cat(2, output.ori{output.inside});
    catch
      %when oris are in wrong orientation (row rather than column)
      for k = 1:numel(output.inside)
        ori(:,output.inside(k)) = output.ori{output.inside(k)}';
      end
     end
    output.ori = ori;
  end
  current = 'new';
 
elseif strcmp(current, 'new') && strcmp(type, 'old')
  %go from new to old
  error('not implemented yet');
end

if strcmp(current, 'new') && strcmp(haspow, 'yes'), 

  %----------------------------------------------
  %convert mom into pow if requested and possible
  convert = 0;
  if isfield(output, 'mom') && size(output.mom{output.inside(1)},2)==1,
    convert = 1;
  else
    warning('conversion from mom to pow is not possible, either because there is no mom in the data, or because the dimension of mom>1. in that case call ft_sourcedescriptives first with cfg.projectmom');
  end
   
  if isfield(output, 'cumtapcnt')
    convert = 1 & convert;
  else
    warning('conversion from mom to pow will not be done, because cumtapcnt is missing');
  end
    
  if convert,  
    npos = size(output.pos,1);
    nrpt = numel(output.cumtapcnt);
    tmpmom = cat(2,output.mom{output.inside});
    tmppow = zeros(npos, nrpt);
    tapcnt = [0;cumsum(output.cumtapcnt(:))];
    for k = 1:nrpt
      ntap = tapcnt(k+1)-tapcnt(k);  
      tmppow(output.inside,k) = sum(abs(tmpmom((tapcnt(k)+1):tapcnt(k+1),:)).^2,1)./ntap;
    end
    output.pow       = tmppow;
    output.powdimord = ['pos_rpt_freq']; 
  end  

elseif strcmp(current, 'old') && strcmp(haspow, 'yes')
  warning('construction of single trial power estimates is not implemented here using old style source representation');

end
 

%--------------------------------------------------------
function [dimord] = createdimord(output, fname, rptflag);

if nargin==2, rptflag = 0; end

tmp = getfield(output, fname);

dimord = '';
dimnum = 1;
hasori = isfield(output, 'ori'); %if not, this is probably singleton and not relevant at the end

if iscell(tmp) && (size(output.pos,1)==size(tmp,dimnum) || size(output.pos,1)==size(tmp,2))
  dimord = [dimord,'{pos}'];
  dimnum = dimnum + 1;
elseif ~iscell(tmp) && size(output.pos,1)==size(tmp,dimnum)
  dimord = [dimord,'pos'];
  dimnum = dimnum + 1;
end

switch fname
  case 'cov'
    if hasori, dimord = [dimord,'_ori_ori']; end;
  case 'csd'
    if hasori, dimord = [dimord,'_ori_ori']; end;
  case 'csdlabel'
    dimord = dimord;
  case 'filter'
    dimord = [dimord,'_ori_chan']; 
  case 'leadfield'
    %if hasori,
      dimord = [dimord,'_chan_ori'];
    %else
    %  dimord = [dimord,'_chan'];
    %end
  case 'mom'
    if isfield(output, 'cumtapcnt') && sum(output.cumtapcnt)==size(tmp{output.inside(1)},1)
      if hasori,
        dimord = [dimord,'_rpttap_ori'];
      else
        dimord = [dimord,'_rpttap'];
      end
    elseif isfield(output, 'time')
      if rptflag,
        dimord = [dimord,'_rpt'];
        dimnum = dimnum + 1;
      end
      if numel(output.time)==size(tmp{output.inside(1)},dimnum)
        dimord = [dimord,'_ori_time'];
      end
    end
    
    if isfield(output, 'freq') && numel(output.freq)>1,
      dimord = [dimord,'_freq'];
    end    
  case 'nai'
    if isfield(output, 'freq') && numel(output.freq)==size(tmp,dimnum)
      dimord = [dimord,'_freq'];
    end
  case 'noise'
    if isfield(output, 'freq') && numel(output.freq)==size(tmp,dimnum)
      dimord = [dimord,'_freq'];
    end
  case 'noisecsd'
    if hasori, dimord = [dimord,'_ori_ori']; end
  case 'ori'
    dimord = '';
  case 'pow'
    if isfield(output, 'cumtapcnt') && numel(output.cumtapcnt)==size(tmp,dimnum)
      dimord = [dimord,'_rpt'];
      dimnum = dimnum + 1;
    end
    
    if isfield(output, 'freq') && numel(output.freq)>1 && numel(output.freq)==size(tmp,dimnum)
      dimord = [dimord,'_freq'];
    end
    
    if isfield(output, 'time') && numel(output.time)>1 && numel(output.time)==size(tmp,dimnum)
      dimord = [dimord,'_time'];
    end
    
    otherwise
      warning(sprintf('skipping unknown fieldname %s', fname));
      %error(sprintf('unknown fieldname %s', fname));
end
