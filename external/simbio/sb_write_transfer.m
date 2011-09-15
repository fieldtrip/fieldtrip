function sb_write_transfer(transfer,filename)

% sb_write_transfer writes a transfer matrix on disk
%
% Use as
%   sb_write_transfer(matrix,file_name)
% 
% where matrix is the FE-transfermatrix

% Copyright (C) 2011, Johannes Vorwerk

% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailed information

% write .mat file

NumberRows = size(transfer.mat,1);
NumberCols = size(transfer.mat,2);
mat_filename = [filename, '.mat'];
mab_filename = [filename, '.bin'];

try
  fid = fopen(mat_filename,'w');
  fprintf(fid,'Rows:\t%d\n',NumberRows);
  fprintf(fid,'Cols:\t%d\n',NumberCols);
  fprintf(fid,'Bytes:\t%d\n',8);
  fprintf(fid,'FileMatrix:\t%s\n',mab_filename);
  fprintf(fid,'FEMchecksum:\t%s\n',transfer.checksum);
  fclose(fid);

catch
  disp('Error in writing the .mat file')
  rethrow(ME)
end

try
  fid = fopen(mab_filename,'w');
  fwrite(fid, transfer.mat','double');
  fclose(fid);

catch
  disp('Error in writing the .bin file')
  rethrow(ME)
end
