function [dat] = ft_fetch_data(data, varargin)

% FT_FETCH_DATA mimics the behaviour of FT_READ_DATA, but for a FieldTrip
% raw data structure instead of a file on disk.
%
% Use as
%   [dat] = ft_fetch_data(data, ...)
%
% See also FT_READ_DATA, FT_FETCH_HEADER, FT_FETCH_EVENT

% Copyright (C) 2009-2013, Jan-Mathijs Schoffelen, Robert Oostenveld
% Copyright (C) 2008, Esther Meeuwissen
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

% check whether input is data
skipcheckdata = ft_getopt(varargin, 'skipcheckdata');
if isempty(skipcheckdata) || skipcheckdata ~= 1
  data = ft_checkdata(data, 'datatype', 'raw', 'hassampleinfo', 'yes');
end

% get the options
hdr          = ft_getopt(varargin, 'header');
begsample    = ft_getopt(varargin, 'begsample');
endsample    = ft_getopt(varargin, 'endsample');
chanindx     = ft_getopt(varargin, 'chanindx');
allowoverlap = ft_getopt(varargin, 'allowoverlap', false);
allowoverlap = istrue(allowoverlap);

if isempty(hdr)
  hdr = ft_fetch_header(data);
end

if isempty(begsample) || isempty(endsample)
  error('begsample and endsample must be specified');
end

if isempty(chanindx)
  chanindx = 1:hdr.nChans;
end

% get trial definition according to original data file
if isfield(data, 'sampleinfo')
  trl = data.sampleinfo;
else
  error('data does not contain a consistent trial definition, fetching data is not possible');
end
trlnum = length(data.trial);

% start with the output data being all NaN
dat = nan(numel(chanindx), endsample-begsample+1);

if trlnum>1,
  % original implementation, used when the input data has multiple trials
  
  trllen = zeros(trlnum,1);
  for trllop=1:trlnum
    trllen(trllop) = size(data.trial{trllop},2);
  end
  
  % check whether data.trial is consistent with trl
  if size(trl,1)~=length(data.trial)
    error('trial definition is not internally consistent')
  elseif any(trllen~=(trl(:,2)-trl(:,1)+1))
    error('trial definition is not internally consistent')
  end
  
  minchan = min(chanindx);
  maxchan = max(chanindx);
  if minchan<1 || maxchan>hdr.nChans
    error('selected channels are not present in the data')
  end
  
  % these are for bookkeeping
  maxsample = max([trl(:,2); endsample]);
  count     = zeros(1, maxsample, 'int32');
  trialnum  = zeros(1, maxsample, 'int32');
  samplenum = zeros(1, maxsample, 'int32');
  
  % determine for each sample in the data where it originates from
  for trllop=1:trlnum
    trlbeg = trl(trllop,1);
    trlend = trl(trllop,2);
    if trlbeg>endsample || trlend<begsample
      % skip this piece, it is not interesting because the requested range falls completely outside
    else
      % make vector with 0= no sample of old data, 1= one sample of old data, 2= two samples of old data, etc
      count(trlbeg:trlend) = count(trlbeg:trlend) + 1;
      % make vector with 1's if samples belong to trial 1, 2's if samples belong to trial 2 etc. overlap/ no data --> Nan
      trialnum(trlbeg:trlend) = trllop;
      % make samplenum vector with samplenrs for each sample in the old trials
      samplenum(trlbeg:trlend) = 1:trllen(trllop);
    end
  end
  
  % overlap --> NaN
  % trialnum(count>1)  = NaN;
  % samplenum(count>1) = NaN;
  
  % make a subselection for the desired samples
  if begsample<1
    count     = count    (1:endsample);
    trialnum  = trialnum (1:endsample);
    samplenum = samplenum(1:endsample);
  else
    count     = count    (begsample:endsample);
    trialnum  = trialnum (begsample:endsample);
    samplenum = samplenum(begsample:endsample);
  end
  
  % check if all samples are present and are not present twice or more
  if any(count>1)
    if ~allowoverlap
      % error('some of the requested samples occur twice in the data');
      % this  can be considered OK if the overlap has exactly identical values
      sel = find(count>1); % must be row vector
      for smplop=sel
        % find in which trials the sample occurs
        seltrl = find(smplop>=trl(:,1) & smplop<=trl(:,2));  % which trials
        selsmp = smplop - trl(seltrl,1) + 1;                 % which sample in each of the trials
        for i=2:length(seltrl)
          % compare all occurences to the first one
          if ~all(data.trial{seltrl(i)}(:,selsmp(i)) == data.trial{seltrl(1)}(:,selsmp(1)))
            error('some of the requested samples occur twice in the data and have conflicting values');
          end
        end
      end
    else
      warning('samples present in multiple trials, using only the last occurence of each sample')
    end
  end
  
  %   if any(count==0)
  %     warning('not all requested samples are present in the data, filling with NaNs');
  %   end
  
  % construct the output data array
  % dat = nan(length(chanindx), length(samplenum));
  % for smplop=1:length(samplenum)
  %   if samplenum(smplop)==0
  %    dat(:, smplop) = nan;
  %   else
  %    dat(:, smplop) = data.trial{trialnum(smplop)}(chanindx,samplenum(smplop));
  %   end
  % end
  
  % the following piece of code achieves the same as the commented code above,
  % but much smaller. rather than looping over samples it loops over the blocks
  % of samples defined by the original trials
  
  utrl = unique(trialnum);
  utrl(~isfinite(utrl)) = 0;
  utrl(utrl==0) = [];
  if length(utrl)==1,
    ok   = trialnum==utrl;
    smps = samplenum(ok);
    dat(:,ok) = data.trial{utrl}(chanindx,smps);
  else
    for xlop=1:length(utrl)
      ok   = trialnum==utrl(xlop);
      smps = samplenum(ok);
      dat(:,ok) = data.trial{utrl(xlop)}(chanindx,smps);
    end
  end
  
else
  % only one trial is present in the input data, so it's quite simple and can be done much faster
  
  % get the indices
  begindx  = begsample - trl(1) + 1;
  endindx  = endsample - trl(1) + 1;
  
  tmptrl = trl([1 2]) - [trl(1) trl(1)]+1; % ignore offset in case it's present
  
  datbegindx = max(1,                     trl(1)-begsample+1);
  datendindx = min(endsample-begsample+1, trl(2)-begsample+1);
  
  %   if begsample<trl(1) || endsample>trl(2)
  %     warning('not all requested samples are present in the data, filling with NaNs');
  %   end
  
  if begsample >= trl(1) && begsample <= trl(2)
    % the begin sample can be found in the available data
    if endsample >= trl(1) && endsample <= trl(2)
      dat(:,datbegindx:datendindx) = data.trial{1}(chanindx,begindx:endindx);
    else
      dat(:, datbegindx:datendindx) = data.trial{1}(chanindx,begindx:tmptrl(2));
    end
  elseif endsample >= trl(1) && endsample <= trl(2)
    % the end sample can be found in the available data
    dat(:, datbegindx:datendindx) = data.trial{1}(chanindx,tmptrl(1):endindx);
  else
    % neither the begin, nor the end sample are in the available data
    dat(:, datbegindx:datendindx) = data.trial{1}(chanindx,tmptrl(1):tmptrl(2));
  end
  
end % if trlnum is multiple or one
