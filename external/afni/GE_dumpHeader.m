function GE_dumpHeader(IFileName, OutFile)
%
% GE_dumpHeader(IFileName,OutFile)
% Dumps much of the header from a GE lx2 or 5.X file to
% OutFile.  If OutFile is not specified, the output goes
% to stdout.
%
% Souheil J. Inati, PhD
% Dartmouth College
% May 2000
% souheil.inati@dartmouth.edu
%

%%%% Call GE_readHeader %%%%
[su_hdr,ex_hdr,se_hdr,im_hdr,pix_hdr] = GE_readHeader(IFileName);

% Open the output file
if nargin == 2
  outid = fopen(OutFile,'w');
else
  outid = 1;  % stdout
end

%%%% Exam Header %%%%
fprintf(outid,'Exam Header Section:\n');
fprintf(outid,GE_dumpExamHeader(ex_hdr));
fprintf(outid,'\n\n');

%%%% Series Header %%%%
fprintf(outid,'Series Header Section:\n');
fprintf(outid,GE_dumpSeriesHeader(se_hdr));
fprintf(outid,'\n\n');

%%%% Image Header %%%%
fprintf(outid,'Image Header Section:\n');
fprintf(outid,GE_dumpImageHeader(im_hdr));

% Close the File if not stdout
if outid ~= 1
  fclose(outid);
end

return
