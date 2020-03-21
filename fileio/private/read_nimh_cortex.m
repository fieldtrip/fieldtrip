function cortex = read_nimh_cortex(filename, varargin)

% READ_NIMH_CORTEX
%
% Use as
%  cortex = read_nimh_cortex(filename, ...)
%
% Optional input arguments should come in key-value pairs and may
% include
%   begtrial     = number (default = 1)
%   endtrial     = number (default = inf)
%   epp          = read the EPP data, 'yes' or 'no' (default = 'yes')
%   eog          = read the EOG data, 'yes' or 'no' (default = 'yes')
%   feedback     = display the progress on the screen, 'yes' or 'no' (default = 'no')
%
% The output is a structure array with one structure for every trial that was read.

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

% get the optional input arguments
feedback = ft_getopt(varargin, 'feedback', 'no');
begtrial = ft_getopt(varargin, 'begtrial', 1);
endtrial = ft_getopt(varargin, 'endtrial', inf);
% reading the epp and eog data is optional
epp = ft_getopt(varargin, 'epp', 'yes');
eog = ft_getopt(varargin, 'eog', 'yes');

% convert 'yes' or 'no' into boolean
skipepp   = ~istrue(epp);
skipeog   = ~istrue(eog);
clear epp eog

% this will hold the result
cortex  = struct;

% trials are counted with zero-offset
trial = 0;

i_fid = fopen_or_error(filename, 'rb', 'ieee-le');

% read until the end of file or until the specified trial number
while ( ~feof (i_fid) && trial<endtrial )

  % the header of all trials must be read, otherwise the size of the subsequent trial data is not known
  % the header of uninteresting trials is not added to the output
  % the epp and eog data may be skipped completely
  skiptrial = (trial+1)<begtrial;

  length = fread(i_fid, 1, 'ushort');

  % allocate space for the header details
  hd = zeros(1,13);

  if ~isempty(length)
    
    % give some feedback on the screen
    if feedback && (trial+1)>=begtrial
      fprintf('reading trial %d\n', trial+1);
    end

    % Number of bytes for each variable.
    hd(1,1:8)= (fread(i_fid, 8, 'ushort'))';
    
    % Convert bytes to number of float point values (4 bytes apiece in windows and FreeBSD).
    hd(1,5)=hd(1,5)/4;
    % Convert bytes to number of short values (2 bytes apiece in windows and FreeBSD).
    hd(1,6)=hd(1,6)/2;
    hd(1,7)=hd(1,7)/2;
    hd(1,8)=hd(1,8)/2;

    header.cond_no =   hd(1,1);
    header.repeat_no = hd(1,2);
    header.block_no =  hd(1,3);
    header.trial_no =  hd(1,4);
    header.isi_size =  hd(1,5);
    header.code_size = hd(1,6);
    header.eog_size =  hd(1,7);
    header.epp_size =  hd(1,8);

    hd(1,9:10) = (fread(i_fid, 2, 'uchar'))';

    header.kHz_resolution   = hd(1,9);
    header.eye_storage_rate = hd(1,10);

    hd(1,11:13) = (fread(i_fid, 3, 'ushort'))';

    header.expected_response = hd(1,11);
    header.response          = hd(1,12);
    header.response_error    = hd(1,13);

    if skiptrial
      fseek(i_fid,header.isi_size,  'cof');
      fseek (i_fid,header.code_size, 'cof');
    else
      time  = (fread (i_fid,header.isi_size,  'ulong'));
      event = (fread (i_fid,header.code_size, 'ushort'));
    end % if skiptrial
    
    if skipepp || skiptrial
      epp1 = [];
      epp2 = [];
      fseek (i_fid,header.epp_size * 2, 'cof');

    else
      epp = fread (i_fid,header.epp_size, 'short');
      epp1 = zeros(1,header.epp_size/2);
      epp2 = zeros(1,header.epp_size/2);

      % fprintf(o_fid, '\nepp(x)\tepp(y)\n\n');

      for i = 1:2:header.epp_size
        % Must extract the data from the raw epp values.
        % The epp value is made up of 12-bits of data, and 4-bits (the
        % low-order 4 bits) of the channel number.
        % To extract the data, must right shift the raw data by 4 bits (to
        % get rid of the channel number and put the data value in the proper
        % location).  After this conversion, you must still add or subtract
        % 2048, since otherwise the value is not right.  (I think that this
        % is because of the way that matlab handles negative values during
        % the bitwise operations.)
        % These calculations cause the results of ctx2txt.m to be the same as
        % for cortview.exe for the EOG and EPP values.

        s = (i+1)/2;
        epp1(s) = bitshift(epp(i),   -4);
        epp2(s) = bitshift(epp(i+1), -4);

        if (epp1(s) < 0)
          epp1(s) = epp1(s) + 2047;
        else
          epp1(s) = epp1(s) - 2048;
        end;

        if (epp2(s) < 0)
          epp2(s) = epp2(s) + 2047;
        else
          epp2(s) = epp2(s) - 2048;
        end;
      end;
    end; % if skipepp

    if skipeog || skiptrial
      eog = [];
      fseek (i_fid,header.eog_size * 2, 'cof');
    else
      eog = fread (i_fid,header.eog_size, 'short');
      eog = reshape(eog, 2, header.eog_size/2);
    end % if skipeog

    if ~skiptrial
      % collect the headers of all trials in a structure array
      cortex(trial+1).header = header;
      cortex(trial+1).length = length;
      cortex(trial+1).time   = time;
      cortex(trial+1).event  = event;
      cortex(trial+1).epp    = [epp1; epp2];
      cortex(trial+1).eog    = eog;
    end

    trial = trial+1;
  end; % if ~isempty(length)
end;

fclose(i_fid);

% remove the initial empty trials, i.e. the ones that were skipped
cortex = cortex(begtrial:end);

