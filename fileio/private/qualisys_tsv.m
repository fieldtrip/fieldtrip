function varargout = qualisys_tsv(filename, hdr, begsample, endsample, chanindx)

% QUALISYS_TSV reads motion tracking data from a Qualisys tsv file.
%
% Use as
%   hdr = qualysis_tsv(filename);
%   dat = qualysis_tsv(filename, hdr, begsample, endsample, chanindx);
%   evt = qualysis_tsv(filename, hdr);
%
% See also FT_FILETYPE, FT_READ_HEADER, FT_READ_DATA, FT_READ_EVENT

% Copyright (C) 2019 Robert Oostenveld
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

needhdr = (nargin==1);
needevt = (nargin==2);
needdat = (nargin==5);

if needhdr
  %% read the header
  
  orig = struct();
  fid = fopen_or_error(filename, 'r');
  while true
    line = fgetl(fid);
    if any('0123456789' == line(1))
      break
    else
      tok = strsplit(line);
      key = lower(tok{1});
      if numel(tok)>2
        % keep it as cell-array of strings
        val = tok(2:end);
      elseif numel(tok)==2 && ~isnan(str2double(tok{2}))
        % represent it as a number
        val = str2double(tok{2});
      else
        % represent it as a string
        val = tok{2};
      end
      orig.(key) = val;
    end
  end
  fclose(fid);
  
  % convert it into a FieldTrip header
  hdr = [];
  hdr.label = {};
  for i=1:orig.no_of_markers
    switch orig.data_included
      case '2D'
        hdr.label{end+1} = [orig.marker_names{i} '_x'];
        hdr.label{end+1} = [orig.marker_names{i} '_y'];
      case '3D'
        hdr.label{end+1} = [orig.marker_names{i} '_x'];
        hdr.label{end+1} = [orig.marker_names{i} '_y'];
        hdr.label{end+1} = [orig.marker_names{i} '_z'];
      case '6D'
        % FIXME I don't know how the data is represented, but suspect some quaternion-like representation
        hdr.label{end+1} = [orig.marker_names{i} '_q1'];
        hdr.label{end+1} = [orig.marker_names{i} '_q2'];
        hdr.label{end+1} = [orig.marker_names{i} '_q3'];
        hdr.label{end+1} = [orig.marker_names{i} '_q4'];
        hdr.label{end+1} = [orig.marker_names{i} '_q5'];
        hdr.label{end+1} = [orig.marker_names{i} '_q6'];
    end
  end
  hdr.nChans      = numel(hdr.label);
  hdr.nSamples    = orig.no_of_frames;
  hdr.nSamplesPre = 0; % continuous data
  hdr.nTrials     = 1; % continuous data
  hdr.Fs          = orig.frequency;
  hdr.chantype    = repmat({'motion'}, size(hdr.label));
  hdr.chanunit    = repmat({'unknown'}, size(hdr.label));
  
  % keep the original header details
  hdr.orig = orig;
  
  % return the header details
  varargout = {hdr};
  
elseif needdat
  %% read the data
  
  tab = readtable(filename, 'FileType', 'text', 'Delimiter', '\t');
  dat = table2array(tab(begsample:endsample, chanindx))';
  
  % return the data
  varargout = {dat};
  
elseif needevt
  %% read the events
  
  ft_error('not yet implemented');
  
end

