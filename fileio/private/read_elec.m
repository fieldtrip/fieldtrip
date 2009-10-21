function [el, lab] = read_elec(fn);

% READ_ELEC reads "la/mu" electrode parameters from a MBF electrode file
% which are used to position them on a triangulated surface
%
% [el, lab] = read_elec(filename)
%
% where el = [dhk, la, mu]
% and lab contains the electrode labels (if present)
%
% See also READ_TRI, TRANSFER_ELEC

% Copyright (C) 1998, Robert Oostenveld
%
% $Log: read_elec.m,v $
% Revision 1.1  2009/01/14 09:24:45  roboos
% moved even more files from fileio to fileio/privtae, see previous log entry
%
% Revision 1.2  2003/03/11 15:24:51  roberto
% updated help and copyrights
%

fid = fopen(fn, 'rt');
if fid~=-1

  % read the number of electrodes
  Nel = sscanf(fgetl(fid), '%d'); 
 
  % read the electrode triangle, lambda and mu
  for i=1:Nel
    str = fgetl(fid);
    el(i,:)  = sscanf(str, '%f %f %f')';
    indx = find(str=='!');
    if (indx)
      lab(i,:) = sprintf('%6s', str((indx+1):length(str)));
    end
  end
  fclose(fid);

else
  error('unable to open file');
end

