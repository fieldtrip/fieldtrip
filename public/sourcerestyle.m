function [output] = sourcerestyle(input, type)

% SOURCERESTYLE converts old style source structures into new style source structures and the
% other way around
%
% Use as:
%   sourcerestyle(input, type)
%    where input is a source structure,
%
% Typically, old style source structures contain
%   avg.XXX or trial.XXX fields
%
% The ne wstyle source structure contains:
%   source.pos
%   source.dim (optional, if the list of positions describes a 3D volume
%   source.XXX the old style subfields in avg/trial
%   source.XXXdimord string how to interpret the respective XXX field:
%     e.g. source.leadfield = cell(1,Npos), source.leadfielddimord = '{pos}_chan_ori'
%          source.mom       = cell(1,Npos), source.momdimord       = '{pos}_ori_rpttap'

if nargin==1,
  type = 'old';
end

fnames = fieldnames(input);
tmp    = cell2mat(strfind(fnames, 'dimord')); %get dimord like fields
if any(tmp>1),
  current = 'new';
elseif any(tmp==1),
  %don't know what to do yet data is JM's own invention
else
  current = 'old';
end

if strcmp(current, type),
  %do nothing
  return
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
  
  %remove and rename the specified fields if present
  removefields = {'xgrid';'ygrid';'zgrid';'method'};
  renamefields = {'frequency' 'freq'; 'csdlabel' 'orilabel'};
  fnames       = fieldnames(output);
  for k = 1:numel(fnames)
    ix = strmatch(fnames{k}, removefields);
    if ~isempty(ix),
      output = rmfield(output, fnames{k});
    end
    ix = strmatch(fnames{k}, renamefields(:,1));
    if ~isempty(ix),
      output = setfield(output, renamefields{ix,2}, ...
                        getfield(output, renamefields{ix,1}));
      output = rmfield(output, fnames{k});
    end
  end

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
      if isfield(input, 'cumtapcnt') && strcmp(fnames{k}, 'mom'),
        %pcc based mom is special case: rpttap rather than rpt
        %and single voxel mom = ori x rpttap
        nrpttap = sum(input.cumtapcnt);
        sizvox  = [size(tmp{input.inside(1)}) 1];
        sizvox  = [sizvox(1) nrpttap sizvox(3:end)];
        rptstring = 'rpttap';
      elseif iscell(tmp)
        nrpttap = numel(stuff);
        sizvox  = [size(tmp{input.inside(1)}) 1];
        sozvpx  = [sizvox(1) nrpttap sizvox(3:end)];
        rptstring = 'rpt';
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
            tmpall{indx}(:,cnt+[1:siz1],:,:,:) = tmpdat;
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
 
elseif strcmp(current, 'new') && strcmp(type, 'old')
  %go from new to old
  error('not implemented yet');
end


%-----------------------------------------------
function [dimord] = createdimord(output, fname, rptflag);

if nargin==2, rptflag = 1; end

tmp = getfield(output, fname);

dimord = '';
dimnum = 1;

if iscell(tmp) && size(output.pos,1)==size(tmp,dimnum)
  dimord = [dimord,'{pos}'];
  dimnum = dimnum + 1;
elseif ~iscell(tmp) && size(output.pos,1)==size(tmp,dimnum)
  dimord = [dimord,'pos'];
  dimnum = dimnum + 1;
end

switch fname
  case 'csd'
    dimord = [dimord,'_ori_ori'];
  case 'csdlabel'
    dimord = dimord;
  case 'filter'
    dimord = [dimord,'_ori_chan']; 
  case 'leadfield'
    dimord = [dimord,'_chan_ori'];
  case 'mom'
    if isfield(output, 'cumtapcnt') && sum(output.cumtapcnt)==size(tmp{output.inside(1)},2)
      dimord = [dimord,'_ori_rpttap'];
    elseif isfield(output, 'time') && numel(output.time)==size(tmp{output.inside(1)},2)
      dimord = [dimord,'_ori_time'];
    end
    
    if isfield(output, 'freq')
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
    dimord = [dimord,'_ori_ori'];
  case 'ori'
    %this field is equivalent to a pos-field
    %FIXME should this be matricized (size is not that big)
    dimord = '';
  case 'pow'
    if isfield(output, 'cumtapcnt') && numel(output.cumtapcnt)==size(tmp,dimnum)
      dimord = [dimord,'_rpt'];
      dimnum = dimnum + 1;
    end
    
    if isfield(output, 'freq') && numel(output.freq)==size(tmp,dimnum)
      dimord = [dimord,'_freq'];
    end
  otherwise
    error(sprintf('unknown fieldname %s', fname));
end
