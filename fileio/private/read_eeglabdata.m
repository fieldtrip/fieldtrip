% read_eeglabdata() - import EEGLAB dataset files
%
% Usage:
%   >> dat = read_eeglabdata(filename);
%
% Inputs:
%   filename - [string] file name
%
% Optional inputs:
%   'begtrial'  - [integer] first trial to read
%   'endtrial'  - [integer] last trial to read
%   'chanindx'  - [integer] list with channel indices to read
%   'header'    - FILEIO structure header
%
% Outputs:
%   dat    - data over the specified range
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2008-

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2008 Arnaud Delorme, SCCN, INC, UCSD, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: read_eeglabdata.m,v $
% Revision 1.3  2009/07/01 15:35:19  vlalit
% Try several different ways of looking for the data file
%  before giving up (fix suggested by Jakob Scherer)
%
% Revision 1.2  2009/02/27 12:03:09  vlalit
% Arno's fix for a bug reported by Antanas Spokas
%
% Revision 1.1  2009/01/30 04:01:19  arno
% *** empty log message ***
%
% Revision 1.2  2008/04/21 18:45:23  roboos
% fixed bug due to sample/trial selection mixup
% only read the selected trials
%
% Revision 1.1  2008/04/18 14:04:48  roboos
% new implementation by Arno, shoudl be tested
%

function dat = read_eeglabdata(filename, varargin);

if nargin < 1
  help read_eeglabdata;
  return;
end;

header    = keyval('header',     varargin);
begsample = keyval('begsample',  varargin);
endsample = keyval('endsample',  varargin);
begtrial  = keyval('begtrial',   varargin);
endtrial  = keyval('endtrial',   varargin);
chanindx  = keyval('chanindx',   varargin);

if isempty(header)
  header = read_eeglabheader(filename);
end;

if ischar(header.orig.data)
  if strcmpi(header.orig.data(end-2:end), 'set'),
    header.ori = load('-mat', filename);
  else
      
      % assuming that the data file is in the current directory
      fid = fopen(header.orig.data);

      % assuming the .dat and .set files are located in the same directory
      if fid == -1
          pathstr = fileparts(filename);
          fid = fopen(fullfile(pathstr, header.orig.data));
      end

      if fid == -1
          fid = fopen(fullfile(header.orig.filepath, header.orig.data)); %
      end

    if fid == -1, error('Cannot not find data file'); end;

    % only read the desired trials
    if strcmpi(header.orig.data(end-2:end), 'dat')
        dat = fread(fid,[header.nSamples*header.nTrials header.nChans],'float32')';
    else
        dat = fread(fid,[header.nChans header.nSamples*header.nTrials],'float32')';
    end;
    dat = reshape(dat, header.nChans, header.nSamples, header.nTrials);
    fclose(fid);
  end;
else
  dat = header.orig.data;
  dat = reshape(dat, header.nChans, header.nSamples, header.nTrials);
end;

if isempty(begtrial), begtrial = 1; end;
if isempty(endtrial), endtrial = header.nTrials; end;
if isempty(begsample), begsample = 1; end;
if isempty(endsample), endsample = header.nSamples; end;
dat = dat(:,begsample:endsample,begtrial:endtrial);

if ~isempty(chanindx)
  % select the desired channels
  dat = dat(chanindx,:,:);
end
