function smi = read_smi_txt(filename)

% READ_SMI_TXT reads the header information, input triggers, messages
% and all data points from an SensoMotoric Instruments (SMI) *.txt file
%
% Use as
%   smi = read_smi_txt(filename)

% Copyright (C) 2016, Diego Lozano-Soldevilla (CerCo)
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org/
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

smi.header  = {};
smi.label   = {};
smi.dat     = [];
current     = 0;

% read the whole file at once
fid = fopen_or_error(filename, 'rt');
aline = fread(fid, inf, 'char=>char');          % returns a single long string
fclose(fid);

aline(aline==uint8(sprintf('\r'))) = [];        % remove cariage return
aline = tokenize(aline, uint8(sprintf('\n')));  % split on newline

for i=1:numel(aline)
  tline = aline{i};
  
  if numel(tline) && any(tline(1)=='0':'9')
    % if regexp(tline, '^[0-9]')
    tmp = regexp(tline,'[^\t\r\n\f\v]*','match')';
    % exclude SMP sting
    smp = regexp(tmp','SMP');
    smpidx = find(cellfun(@numel,smp)==1);
    
    current = current + 1;
    smi.type{1,current} =  sscanf(tmp{smpidx}, '%s');
    tmp(smpidx,:) = [];
    
    for j=1:size(tmp,1)
      val = sscanf(tmp{j,1}, '%f');
      if isempty(val)
        smi.dat(j,current) =  NaN;
      else
        smi.dat(j,current) =  val;
      end
    end
    
  elseif regexp(tline, '##')
    smi.header = cat(1, smi.header, {tline});
    
    template ='## Sample Rate:';
    if strncmp(tline,template,length(template))
      smi.Fs = cellfun(@str2num,regexp(tline,'[\d]+','match'));
    end
    
    template = '## Number of Samples:';
    if strncmp(tline,template,length(template))
      smi.nsmp = cellfun(@str2num,regexp(tline,'[\d]+','match'));
    end
    
    template = '## Head Distance [mm]:';
    if strncmp(tline,template,length(template))
      smi.headdist = cellfun(@str2num,regexp(tline,'[\d]+','match'));
    end
    	    
  elseif regexp(tline, 'Time*')
    smi.label = cat(1, smi.label, {tline});
    smi.label = regexp(smi.label{:},'[^\t\r\n\f\v]*','match')';
    typeln = regexp(smi.label','Type');
    typeidx = find(cellfun(@numel,typeln)==1);
    
    % delete Type column because only contains strings
    smi.label(typeidx,:) = [];
    % preamble data matrix
    smi.dat = zeros(size(smi.label,1),smi.nsmp);
    
  else
    % all other lines are not parsed
  end
  
end

% remove the samples that were not filled with real data
% smi.dat = smi.dat(:,1:current);

% place the timestamp channel outside of the data
c = find(strcmp(smi.label, 'Time'));
if numel(c)==1
  smi.timestamp = smi.dat(c,:);
  smi.dat(c,:) = [];
  smi.label(c) = [];
end
