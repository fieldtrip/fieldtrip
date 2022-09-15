function [err] = writeSPMmat (matin, matout)
%
%   [err] = writeSPMmat (matin, [matout]);
%
%Purpose:
%   Write an ascii version of SPM's .mat files
%   containing a 4x4 transformation matrix.
%   The output of this function can be used
%   with 3dWarp to create AFNI-formatted data from
%   SPM datasets.
%
%Input Parameters:
%   matin : Name of .mat file containing transform
%   matout: (optional) Name of output file.
%           If matout is not specified, the _ASCII
%           extension is added to matin
%
%Output Parameters:
%   err : 0 No Problem
%       : 1 Problems
%
%
%More Info :
%   3dWarp -help
%   https://afni.nimh.nih.gov/ssc/ziad/SUMA/SUMA_doc.htm
%   (search for SPM)
%
%     Author : Ziad Saad
%     Date : Thu Aug 7 11:31:54 EDT 2003
%     SSCC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'writeSPMmat';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

%does the input file exist ?
matin = deblank(matin);
matin = deblank(fliplr(matin));
matin = fliplr(matin);

if (exist(matin) ~= 2),
   fprintf(2,'Error %s:\nFile %s not found.\n', FuncName, matin);
   return;
end

%get its prefix
if (strmatch('tam.', fliplr(matin)) == 1),
   prefixin = matin(1:length(matin)-4);
else
   prefixin = matin;
end

%load it
S = load (matin);
if (~isfield(S,'M')),
   fprintf(2,'Error %s:\nFailed to find M in %s\n', FuncName, matin);
   return;
end

if (size(S.M,2) ~= 4),
   fprintf(2,'Error %s:\nM does not have 4 columns!\nM''s size is %g x %g',...
      FuncName, size(S.M,1), size(S.M,2));
   return;
end

if (size(S.M,1) ~= 3 & size(S.M,1) ~= 4),
   fprintf(2,'Error %s:\nM does not have 3 or 4 rows!\nM''s size is %g x %g',...
      FuncName, size(S.M,1), size(S.M,2));
   return;
end

if (size(S.M,1) == 4),
   if (S.M(4,1) ~= 0 | S.M(4,2) ~= 0 | S.M(4,3) ~= 0 | S.M(4,4) ~= 1),
      beep;
      fprintf(2,'Warning %s:\nAFNI expects the 4th row of M to be [0 0 0 1].\n', FuncName);
      fprintf(2,'This row is currently [%g %g %g %g] and will be ignored by AFNI.\n', ...
         S.M(4,1), S.M(4,2), S.M(4,3), S.M(4,4));
      fprintf(2,'Proceed with caution.\n');
      beep;
   end
end

%write output file
if (nargin == 1),
   matout = [prefixin, '_ASCII.mat'];
end
if (exist(matout) == 2),
   fprintf(2,'Error %s:\nOutput file %s exists already.\nWill not overwrite.\n',...
    FuncName, matout);
   return;
end

fid = fopen(matout,'w');
if (fid < 0),
   fprintf(2,'Error %s:\nFailed in opening %s for writing.\nCheck permissions and disk space.\n',...
      FuncName, matout);
   return;
end
fprintf (fid,'#Conversion of SPM''s %s \n', matin);
fprintf (fid,'# transformation matrix file to ascii format.\n#\n');
fprintf (fid,'#To convert SPM datasets to AFNI''s format, use \n');
fprintf (fid,'# 3dWarp including these options:\n');
fprintf (fid,'# -matvec_in2out %s -matvec_fsl \n', matout);
fprintf (fid,'#See 3dWarp -help for more info.\n#\n');
for (i=1:1:4),
   fprintf (fid,'\t%f\t%f\t%f\t%f\n',...
      S.M(i,1), S.M(i,2), S.M(i,3), S.M(i,4));
end

fclose(fid);
err = 0;
return;

