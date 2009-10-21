function [dap] = read_mpi_dap(filename)

% READ_MPI_DAP read the analog channels from a DAP file
% and returns the values in microvolt (uV)
%
% Use as
%   [dap] = read_mpi_dap(filename)

% Copyright (C) 2005-2007, Robert Oostenveld
%
% $Log: read_mpi_dap.m,v $
% Revision 1.1  2009/01/14 09:24:45  roboos
% moved even more files from fileio to fileio/privtae, see previous log entry
%
% Revision 1.3  2007/03/19 16:53:09  roboos
% renamed local variable from dat into dap
%
% Revision 1.2  2005/10/06 11:01:24  roboos
% added support for average mode
% only read analogsweep+hdr+data if nanalog>0
% only read spikesweep+hdr+data if nspike>0
%
% Revision 1.1  2005/10/04 16:53:21  roboos
% new implementation
%

fid = fopen(filename, 'rb', 'ieee-le');

% read the file header
filehdr        = readheader(fid);
analog         = {};
analoghdr      = {};
analogsweephdr = {};
spike          = {};
spikehdr       = {};
spikesweephdr  = {};

W = filehdr.nexps;
S = filehdr.nsweeps;

for w=1:filehdr.nexps
  % read the experiment header
  exphdr{w} = readheader(fid);

  if filehdr.nanalog
    if filehdr.analogmode==-1
      % record mode
      for s=1:filehdr.nsweeps
        % read the analogsweepheader
        analogsweephdr{s} = readheader(fid);
        for j=1:filehdr.nanalog
          % read the analog header
          analoghdr{w,s,j} = readheader(fid);
          % read the analog data
          analog{w,s,j} = fread(fid, analoghdr{w,s,j}.datasize, 'int16');
          % calibrate the analog data
          analog{w,s,j} = analog{w,s,j} * 2.5/2048;
        end
      end % for s=1:S
    elseif filehdr.analogmode==0
      % average mode
      s = 1;
      for j=1:filehdr.nanalog
        % read the analog header
        analoghdr{w,s,j} = readheader(fid);
        % read the analog data
        analog{w,s,j} = fread(fid, analoghdr{w,s,j}.datasize, 'int16');
        % calibrate the analog data
        analog{w,s,j} = analog{w,s,j} * 2.5/2048;
      end
    else
      error('unknown analog mode'); 
    end
  end

  if filehdr.nspike
    for s=1:filehdr.nsweeps
      spikesweephdr{s} = readheader(fid);
      for j=1:filehdr.nspike
        % read the spike header
        spikehdr{w,s,j} = readheader(fid);
        % read the spike data
        spike{w,s,j} = fread(fid, spikehdr{w,s,j}.datasize, 'int16');
      end
    end % for s=1:S
  end

end % for w=1:W

dap.filehdr        = filehdr;
dap.exphdr         = exphdr;
dap.analogsweephdr = analogsweephdr;
dap.analoghdr      = analoghdr;
dap.analog         = analog;
dap.spikesweephdr  = spikesweephdr;
dap.spikehdr       = spikehdr;
dap.spike          = spike;

fclose(fid);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hdr = readheader(fid);
% determine the header type, header size and data size
hdr.headertype = fread(fid, 1, 'uchar');
dummy          = fread(fid, 1, 'uchar');
hdr.headersize = fread(fid, 1, 'int16') + 1;    % expressed in words, not bytes
hdr.datasize   = fread(fid, 1, 'int16');        % expressed in words, not bytes
% read the header details
switch hdr.headertype
  case 1 % fileheader
    % fprintf('fileheader at %d\n' ,ftell(fid)-6);
    dummy           = fread(fid, 1, 'uchar')';    % 6
    hdr.nexps       = fread(fid, 1, 'uchar')';    % 7
    hdr.progname    = fread(fid, 7, 'uchar')';    % 8-14
    hdr.version     = fread(fid, 1, 'uchar')';    % 15
    hdr.progversion = fread(fid, 2, 'uchar')';    % 16-17
    hdr.fileversion = fread(fid, 2, 'uchar')';    % 18-19
    hdr.date        = fread(fid, 1, 'int16');     % 20-21
    hdr.time        = fread(fid, 1, 'int16');     % 22-23
    hdr.nanalog     = fread(fid, 1, 'uint8');     % 24
    hdr.nspike      = fread(fid, 1, 'uint8');     % 25
    hdr.nbins       = fread(fid, 1, 'int16');     % 26-27
    hdr.binwidth    = fread(fid, 1, 'int16');     % 28-29
    dummy           = fread(fid, 1, 'int16');     % 30-31
    hdr.nsweeps     = fread(fid, 1, 'int16');     % 32-33
    hdr.analogmode  = fread(fid, 1, 'uchar')';    % 34    "0 for average, -1 for record"
    dummy           = fread(fid, 1, 'uchar')';    % 35
    dummy           = fread(fid, 1, 'int16');     % 36-37
    dummy           = fread(fid, 1, 'int16');     % 38-39
    dummy           = fread(fid, 1, 'int16');     % 40-41
  case 65 % expheader
    % fprintf('expheader at %d\n' ,ftell(fid)-6);
    hdr.time        = fread(fid, 1, 'int16');     % 6-7
    hdr.parallel    = fread(fid, 1, 'int16');     % 8-9
    hdr.id          = fread(fid, 1, 'int16');     % 10-11
    dummy           = fread(fid, 1, 'int16');     % 12-13
    dummy           = fread(fid, 1, 'int16');     % 14-15
    dummy           = fread(fid, 1, 'int16');     % 16-17
  case 129 % analogchannelheader
    % fprintf('analogchannelheader at %d\n' ,ftell(fid)-6);
    dummy           = fread(fid, 1, 'int16');     % 6-7
    hdr.channum     = fread(fid, 1, 'uchar');     % 8
    dummy           = fread(fid, 1, 'uchar');     % 9
    dummy           = fread(fid, 1, 'int16');     % 10-11
    dummy           = fread(fid, 1, 'int16');     % 12-13
    dummy           = fread(fid, 1, 'int16');     % 14-15
    dummy           = fread(fid, 1, 'int16');     % 16-17
  case 161 % spikechannelheader
    % fprintf('spikechannelheader at %d\n' ,ftell(fid)-6);
    dummy           = fread(fid, 1, 'int16');     % 6-7
    hdr.channum     = fread(fid, 1, 'uchar');     % 8
    dummy           = fread(fid, 1, 'uchar');     % 9
    dummy           = fread(fid, 1, 'int16');     % 10-11
    dummy           = fread(fid, 1, 'int16');     % 12-13
    dummy           = fread(fid, 1, 'int16');     % 14-15
    dummy           = fread(fid, 1, 'int16');     % 16-17
  case 137 % analogsweepheader
    % fprintf('analogsweepheader at %d\n' ,ftell(fid)-6);
    hdr.sweepnum    = fread(fid, 1, 'int16');     % 6-7
    hdr.parallel    = fread(fid, 1, 'int16');     % 8-9
    dummy           = fread(fid, 1, 'int16');     % 10-11
    dummy           = fread(fid, 1, 'int16');     % 12-13
    dummy           = fread(fid, 1, 'int16');     % 14-15
    dummy           = fread(fid, 1, 'int16');     % 16-17
  case 169 % spikesweepheader
    % fprintf('spikesweepheader at %d\n' ,ftell(fid)-6);
    hdr.sweepnum    = fread(fid, 1, 'int16');     % 6-7
    hdr.parallel    = fread(fid, 1, 'int16');     % 8-9
    dummy           = fread(fid, 1, 'int16');     % 10-11
    dummy           = fread(fid, 1, 'int16');     % 12-13
    dummy           = fread(fid, 1, 'int16');     % 14-15
    dummy           = fread(fid, 1, 'int16');     % 16-17
  otherwise
    error(sprintf('unsupported format for header (%d)', hdr.headertype));
end
