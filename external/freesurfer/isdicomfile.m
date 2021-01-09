function [res, dcminfo] = isdicomfile(fname)
% [res, dcminfo] = isdicomfile(fname)
%
% Determines whether the given file is dicom by
% trying to read the dicom info. If this fails,
% then res=0 and dcminfo=[]. If successful, then
% res=1, and dcminfo is the result of matlabs
% dicominfo.
%


%
% isdicomfile.m
%
% Original Author: Doug Greve
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%


try
  dcminfo = dicominfo(fname);
  res = 1;
catch
  dcminfo = [];
  res = 0;
end

return;

%%%%%%% This was Anders original code, does not always work %%%%%%%%
fid = fopen(fname,'r');
if fid < 0
  res = 0;
else
  stat = fseek(fid,128,'bof'); % move to DICM string
  tmp = char(fread(fid,4,'uchar')');%'
  res = strcmp(tmp,'DICM');
  fclose(fid);
end
