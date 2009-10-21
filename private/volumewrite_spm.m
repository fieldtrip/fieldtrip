function [Va] = volumewrite_spm(filename, data, transform)

% VOLUMEWRITE_SPM write a anatomical or functional volume to img/hdr file
% using the SPM toolbox
%
% Use as
%   [Va] = volumewrite_spm(filename, data, transform)

% Copyright (C) 2006, Robert Oostenveld
%
% $Log: volumewrite_spm.m,v $
% Revision 1.4  2006/06/07 09:34:19  roboos
% changed checktoolbox into hastoolbox
%
% Revision 1.3  2006/05/08 10:05:50  roboos
% fixed writing of volumes of an unsupported datatype by converting them to double
%
% Revision 1.2  2006/05/03 15:08:12  jansch
% temporary fix in the case of logical input; spm does not know how to handle this
%
% Revision 1.1  2006/01/05 13:43:02  roboos
% New function, is used by all functions that write to hdr/img file
% using SPM (volumenormalise, volumesegment and volumewrite).
%

% check whether the required SPM2 toolbox is available
hastoolbox('SPM2', 1);

if 0
  % not all datatypes are supported by SPM, this can be checked with
  spm_type('double') %           -- Double precision floating point number array
  spm_type('single') %           -- Single precision floating point number array
  spm_type('logical') %          -- Logical array
  spm_type('char') %             -- Character array
  spm_type('cell') %             -- Cell array
  spm_type('struct') %           -- Structure array
  spm_type('int8') %             -- 8-bit signed integer array
  spm_type('uint8') %            -- 8-bit unsigned integer array
  spm_type('int16') %            -- 16-bit signed integer array
  spm_type('uint16') %           -- 16-bit unsigned integer array
  spm_type('int32') %            -- 32-bit signed integer array
  spm_type('uint32') %           -- 32-bit unsigned integer array
  spm_type('int64') %            -- 64-bit signed integer array
  spm_type('uint64') %           -- 64-bit unsigned integer array
end

if isnan(spm_type(class(data)))
  % convert every unsupported data type into double
  data = double(data);
end

typ        = spm_type(class(data));
dim        = size(data);
Va         = [];
Va.mat     = transform;
Va.fname   = filename;
Va.dim     = [dim typ];
Va.n       = 1;
Va.pinfo   = [1 0 0]';
Va.private.hdr.dime.datatype = typ;
Va         = spm_create_vol(Va);
Va         = spm_write_vol(Va,data);

