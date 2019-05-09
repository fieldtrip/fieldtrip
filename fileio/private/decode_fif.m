function [info] = decode_fif(orig)

% DECODE_FIF is a helper function for real-time processing of Neuromag data. This
% function is used to decode the content of the optional neuromag_fif chunk(s).
%
% See also DECODE_RES4, DECODE_NIFTI1, SAP2MATLAB

% Copyright (C) 2013 Arjen Stolk & Robert Oostenveld
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


% check that the required low-level toolbox is available
ft_hastoolbox('mne', 1);

global FIFF
if isempty(FIFF)
  FIFF = fiff_define_constants();
end

if isfield(orig, 'neuromag_header')
  % The binary blob was created on the little-endian Intel Linux acquisition
  % computer, whereas the default for fiff files is that they are stored in
  % big-endian byte order. MATLAB is able to swap the bytes on the fly by specifying
  % 'le" or "be" to fopen. The normal MNE fiff_open function assumes that it is big
  % endian, hence here we have to open it as little endian.
  
  filename = tempname;
  % write the binary blob to disk, byte-by-byte to avoid any swapping between little and big-endian content
  F = fopen_or_error(filename, 'w');
  fwrite(F, orig.neuromag_header, 'uint8');
  fclose(F);
  % read the content of the file using the standard reading functions
  [info, meas] = read_header(filename);
  % clean up the temporary file
  delete(filename);
end

% Typically, at the end of acquisition, the isotrak and hpiresult information
% is stored in the neuromag fiff container which can then (offline) be read by
% fiff_read_meas_info. However, for the purpose of head position monitoring
% (see Stolk et al., Neuroimage 2013) during acquisition, this crucial
% information requires to be accessible online. read_isotrak and read_hpiresult
% can extract information from the additionally chunked (neuromag2ft) files.

if isfield(orig, 'neuromag_isotrak')
  filename = tempname;
  % write the binary blob to disk, byte-by-byte to avoid any swapping between little and big-endian content
  F = fopen_or_error(filename, 'w');
  fwrite(F, orig.neuromag_isotrak, 'uint8');
  fclose(F);
  % read the content of the file using the standard reading functions
  [info.dig] = read_isotrak(filename);
  % clean up the temporary file
  delete(filename);
end

if isfield(orig, 'neuromag_hpiresult')
  filename = tempname;
  % write the binary blob to disk, byte-by-byte to avoid any swapping between little and big-endian content
  F = fopen_or_error(filename, 'w');
  fwrite(F, orig.neuromag_hpiresult, 'uint8');
  fclose(F);
  % read the content of the file using the standard reading functions
  [info.dev_head_t, info.ctf_head_t] = read_hpiresult(filename);
  % clean up the temporary file
  delete(filename);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [info, meas] = read_header(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open and read the file as little endian
[fid, tree] = fiff_open_le(filename); % open as little endian
[info, meas] = fiff_read_meas_info(fid, tree);
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dig] = read_isotrak(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global FIFF

% open the isotrak file (big endian)
% (typically stored in meas_info dir during acquisition, no fif extension required)
[fid, tree] = fiff_open(filename);

% locate the Polhemus data
isotrak = fiff_dir_tree_find(tree,FIFF.FIFFB_ISOTRAK);

dig=struct('kind',{},'ident',{},'r',{},'coord_frame',{});
coord_frame = FIFF.FIFFV_COORD_HEAD;
if length(isotrak) == 1
  p = 0;
  for k = 1:isotrak.nent
    kind = isotrak.dir(k).kind;
    pos  = isotrak.dir(k).pos;
    if kind == FIFF.FIFF_DIG_POINT
      p = p + 1;
      tag = fiff_read_tag(fid,pos);
      dig(p) = tag.data;
    else
      if kind == FIFF.FIFF_MNE_COORD_FRAME
        tag = fiff_read_tag(fid,pos);
        coord_frame = tag.data;
      elseif kind == FIFF.FIFF_COORD_TRANS
        tag = fiff_read_tag(fid,pos);
        dig_trans = tag.data;
      end
    end
  end
end
for k = 1:length(dig)
  dig(k).coord_frame = coord_frame;
end

if exist('dig_trans','var')
  if (dig_trans.from ~= coord_frame && dig_trans.to ~= coord_frame)
    clear('dig_trans');
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dev_head_t, ctf_head_t] = read_hpiresult(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global FIFF

% open the hpiresult file (big endian)
% (typically stored in meas_info dir during acquisition, no fif extension required)
[fid, tree] = fiff_open(filename);

% locate the transformation matrix
dev_head_t=[];
ctf_head_t=[];
hpi_result = fiff_dir_tree_find(tree,FIFF.FIFFB_HPI_RESULT);
if length(hpi_result) == 1
  for k = 1:hpi_result.nent
    kind = hpi_result.dir(k).kind;
    pos  = hpi_result.dir(k).pos;
    if kind == FIFF.FIFF_COORD_TRANS
      tag = fiff_read_tag(fid,pos);
      cand = tag.data;
      if cand.from == FIFF.FIFFV_COORD_DEVICE && ...
          cand.to == FIFF.FIFFV_COORD_HEAD
        dev_head_t = cand;
      elseif cand.from == FIFF.FIFFV_MNE_COORD_CTF_HEAD && ...
          cand.to == FIFF.FIFFV_COORD_HEAD
        ctf_head_t = cand;
      end
    end
  end
end
