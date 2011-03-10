function [dat] = ft_fetch_data(data, varargin)

% FT_FETCH_DATA mimics the behaviour of FT_READ_DATA, but for a FieldTrip
% raw data structure instead of a file on disk.
%
% Use as
%   [dat] = ft_fetch_data(data, ...)
%
% See also FT_READ_DATA, FT_FETCH_HEADER, FT_FETCH_EVENT

% Copyright (C) 2008, Esther Meeuwissen
% Copyright (C) 2009-2010, Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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
data = ft_checkdata(data, 'datatype', 'raw', 'hastrialdef', 'yes');
    
% get the options
hdr           = keyval('header',        varargin);
begsample     = keyval('begsample',     varargin);
endsample     = keyval('endsample',     varargin);
chanindx      = keyval('chanindx',      varargin);
    
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
  trl    = data.sampleinfo;
else
  error('data does not contain a consistent trial definition, fetching data is not possible');
end
trlnum = length(data.trial);

if trlnum>1,
  % original implementation 
  
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
  maxsample = max(trl(:,2)); 
  count     = zeros(1, maxsample, 'int32');
  trialnum  = zeros(1, maxsample, 'int32');
  samplenum = zeros(1, maxsample, 'int32');
  
  % determine for each sample in the data where it originates from
  for trllop=1:trlnum
    trlbeg = trl(trllop,1);
    trlend = trl(trllop,2);
    if trlbeg>endsample || trlend<begsample
      % this piece of data is not interesting because it falls outside the requested range
      % skip it to speed up the indexing of the trial and sample numbers
      continue
    end
    if trlbeg <= begsample && trlend >= endsample 
        % all data is in this trial!
        % get the indices of the current trial and break the loop
        trlidx = trllop;
        begindx = begsample - trl(trlidx) + 1;
        endindx = endsample - trl(trlidx) + 1;         
        break;
    end
    % make vector with 0= no sample of old data, 1= one sample of old data, 2= two samples of old data, etc
    count(trlbeg:trlend) = count(trlbeg:trlend) + 1;
    % make vector with 1's if samples belong to trial 1, 2's if samples belong to trial 2 etc. overlap/ no data --> Nan
    trialnum(trlbeg:trlend) = trllop;
    % make samplenum vector with samplenrs for each sample in the old trials
    samplenum(trlbeg:trlend) = 1:trllen(trllop);
  end
  
  if exist('trlidx', 'var') && exist('begindx', 'var') && exist('endindx', 'var')
    % fetch the data and return
    dat = data.trial{trlidx}(chanindx,begindx:endindx);
    clear count trialnum samplenum;
  else
      % overlap --> NaN
      %trialnum(count>1)  = NaN;
      %samplenum(count>1) = NaN;

      % make a subselection for the desired samples
      count     = count(begsample:endsample);
      trialnum  = trialnum(begsample:endsample);
      samplenum = samplenum(begsample:endsample);

      % check if all samples are present and are not present twice or more 
      if any(count==0)
        warning('not all requested samples are present in the data, filling with NaNs');
      elseif any(count>1)
        error('some of the requested samples occur twice in the data');
      end

      % construct the output data array
      %dat = nan(length(chanindx), length(samplenum));
      %for smplop=1:length(samplenum)
      %  if samplenum(smplop)==0
      %   dat(:, smplop) = nan; 
      %  else
      %   dat(:, smplop) = data.trial{trialnum(smplop)}(chanindx,samplenum(smplop)); 
      %  end
      %end

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
  end
else
  % only 1 trial present in the input data, so it's quite simple
  % and can be done fast
  
  % check whether the requested samples are present in the input
  if endsample>trl(2) || begsample<trl(1)
    error('some of the requested samples are outside the input data')
  end
  
  % get the indices
  begindx = begsample - trl(1) + 1;
  endindx = endsample - trl(1) + 1;
  
  % fetch the data
  dat = data.trial{1}(chanindx,begindx:endindx);  
end

