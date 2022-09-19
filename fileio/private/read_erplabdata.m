% read_erplabdata() - import ERPLAB dataset files
%
% Usage:
%   >> dat = read_erplabdata(filename);
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
% Modified from read_eeglabheader

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

function dat = read_erplabdata(filename, varargin)

if nargin < 1
  help read_erplabdata;
  return;
end

header    = ft_getopt(varargin, 'header');
begsample = ft_getopt(varargin, 'begsample');
endsample = ft_getopt(varargin, 'endsample');
begtrial  = ft_getopt(varargin, 'begtrial');
endtrial  = ft_getopt(varargin, 'endtrial');
chanindx  = ft_getopt(varargin, 'chanindx');

if isempty(header)
  header = read_erplabheader(filename);
end

dat = header.orig.bindata;

if isempty(begtrial), begtrial = 1; end
if isempty(endtrial), endtrial = header.nTrials; end
if isempty(begsample), begsample = 1; end
if isempty(endsample), endsample = header.nSamples; end
dat = dat(:,begsample:endsample,begtrial:endtrial);

if ~isempty(chanindx)
  % select the desired channels
  dat = dat(chanindx,:,:);
end

