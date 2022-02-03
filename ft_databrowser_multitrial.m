function cfg = ft_databrowser_multitrial(cfg,data)

% FT_DATABROWSER_MULTITRIAL can be used for visual inspection of data. Artifacts that were
% detected by artifact functions (see FT_ARTIFACT_xxx functions where xxx is the type
% of artifact) are marked. Additionally data pieces can be marked and unmarked as
% artifact by manual selection. The output cfg contains the updated specification of
% the artifacts.

% add one option: how many trials shown on the first screen
cfg.ntrialswin = ft_getopt(cfg, 'ntrialswin', 5);
% retrieve events (ft_databrowser code)
cfg.event = ft_getopt(cfg, 'event'); 
if ~isempty(cfg.event)
  % use the events that the user passed in the configuration
  event = cfg.event;
else
  % fetch the events from the data structure in memory
  event = ft_fetch_event(data);
end

%%%%% shifting sampleinfo and event

% store original events
urevent = event;
% samples to shift
eventsample = [event.sample];
sampleinfo = data.sampleinfo;

% time chunks that will be removed
chunksout = [1 sampleinfo(1,1)-1;
  [sampleinfo(1:end-1,2)+1 sampleinfo(2:end,1)-1];
  sampleinfo(end,2)+1 data.hdr.nSamples];

% samples to be subtracted
si2subtract = zeros(size(sampleinfo));
esmp2subtract = zeros(size(eventsample));
for ichunk = 1:size(chunksout,1)
  si2subtract(sampleinfo > chunksout(ichunk,2)) = si2subtract(sampleinfo > chunksout(ichunk,2)) + diff(chunksout(ichunk,:))+1;
  esmp2subtract(eventsample > chunksout(ichunk,2)) = esmp2subtract(eventsample > chunksout(ichunk,2)) + diff(chunksout(ichunk,:))+1;
end

sampleinfo = sampleinfo - si2subtract;
eventsample = eventsample - esmp2subtract;

[event.ursample] = rep2struct([event.sample]);
[event.sample] = rep2struct(eventsample);

%%
boundarysample = sampleinfo(:,1);
boundary = [];
for i = 1:numel(boundarysample)
  boundary(i).type = '';
  boundary(i).value = '';
  boundary(i).sample = boundarysample(i);
  boundary(i).offset = [];
  boundary(i).duration = 0;
  boundary(i).timestamp = boundary(i).sample / data.fsample;
  boundary(i).ursample = NaN;
end
event = [event boundary];
event = sortstruct(event,'sample');

%%%%%%%%%%% now dealing with artfctdef

if isfield(cfg,'artfctdef')
  artnames = fieldnames(cfg.artfctdef);
  for iart = 1:numel(artnames)
    art = cfg.artfctdef.(artnames{iart}).artifact;
    art2subtract = zeros(size(art));
    for ichunk = 1:size(chunksout,1)
      art2subtract(art > chunksout(ichunk,2)) = art2subtract(art > chunksout(ichunk,2)) + diff(chunksout(ichunk,:)) + 1;
    end
    art = art - art2subtract;
    cfg.artfctdef.(artnames{iart}).artifact = art;
  end
end

data.hdr.nSamples = sum(cellfun(@numel,data.time));
data.sampleinfo = sampleinfo;

cfg.event = event;

%%%%%%%%%% forward to databrowser

cfg.continuous = 'yes'; 
cfg.blocksize = (size(cat(2,data.trial{1:cfg.ntrialswin}),2) / data.fsample);

cfgout = ft_databrowser(cfg,data);

%%%%%%%%%% recover previous events
% (these should not have changed in ft_databrowser)
cfgout.event = urevent;

%%%%%%%%%% switch artifacts back to original time frame

if isfield(cfgout,'artfctdef')
  artnames = fieldnames(cfgout.artfctdef);
  for iart = 1:numel(artnames)
    art = cfgout.artfctdef.(artnames{iart}).artifact;
    for ichunk = 1:size(chunksout,1)
      art(chunksout(ichunk,1) < art) = art(chunksout(ichunk,1) < art) + diff(chunksout(ichunk,:)) + 1;
    end
    cfgout.artfctdef.(artnames{iart}).artifact = art;
  end
end

cfg = cfgout;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SUBFUNCTIONS %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [varargout] = rep2struct(varargin)

% [s.target] = rep2struct(dat)
% replicate the value dat into each element of structure s in field target.
% if dat has same nb or elements as s, each element of dat goes into one
% element of s. if dat is more dimensional and doesn't have the same number
% of elements as s, and has same size along dimension 1, then pass each
% slice into s.target.

varargout = cell(1,nargout);
if numel(varargin) == 1
    dat = varargin{1};
    if numel(dat) == nargout
        for i = 1:nargout
            if iscell(dat(i))
                varargout{i} = dat{i};
            else
                varargout{i} = dat(i);
            end
        end
    elseif size(dat,1) == nargout
        for i = 1:nargout
            varargout{i} = dat(i,:);
        end
    else
        for i = 1:nargout
            varargout{i} = dat;
        end
    end
elseif numel(varargin) == nargout
    for i = 1:nargout
        varargout{i} = varargin{i};
    end
else
    error('Wrong number of arguments');
end
return


function S = sortstruct(S, by, direction,sorter)

% S = sortstruct(S, by, direction,sorter)
% sort vector structure by values in field(s) by (may be cell array of
% str).
% direction is 'ascend' or 'descend'
% sorter is a list of unique values that will determine output order.
%
% if by is empty, no sorting occurs.

if isscalar(S) || isempty(S)
    return
elseif ~isvector(S)
    error('input structure must be a vector')
end
s = size(S);
S = S(:);
if nargin < 2
    by = fieldnames(S);
    warning('sortstruct:sortallfields','no sorting fields provided. sorting by all of them...')
end
defifnotexist('direction','ascend');
if ischar(by)
    defifnotexist('sorter',[]);
    by = {by};
    sorter = {sorter};
else
    defifnotexist('sorter',repmat({[]},1,numel(by)));
end
if all(cellfun(@isempty,by))
    return
end
if ischar(S(1).(by{1}))
    L = {S.(by{1})};% list all values
    U = unique(L);% unique values (sorted)
    U = [sorter{1} setxor(sorter{1},U)];
    nS = [];
    for i_u = 1:numel(U)% for each unique value
        tmpS = S(strcmp(L,U{i_u}));% extract them
        if numel(by) > 1
            tmpS = sortstruct(tmpS,by(2:end),direction,sorter(2:end));% sort under
        end
        nS = [nS; tmpS];
    end
elseif isnumeric(S(1).(by{1}))
    L = [S.(by{1})];% list all values
    U = unique(L);% unique values (sorted)
    U = [sorter{1} setxor(sorter{1},U)];
    nS = [];
    for i_u = 1:numel(U)% for each unique value
        tmpS = S(L == U(i_u));% extract them
        if numel(by) > 1
            tmpS = sortstruct(tmpS,by(2:end),direction,sorter(2:end));% sort under
        end
        nS = [nS; tmpS];
    end
else
    error todo
end
S = nS;
S = reshape(S,s);

