function [l, sh] = TellAfni_Commands (f, fpart )
%
%   [l, sh] = TellAfni_Commands (f, field)
%
%Purpose:
%   Extract, all the commands from README.driver
%   No guarantee that you'll get all the commands extracted, 
%   Use this as a guide but read the entire README.driver for the full picture
%   
%Input Parameters:
%   f: Full path and name of README.driver file
%      pass no parameters if you want the function to try and
%      locate the file automatically. 
%      Use '' if you want to look for the file interactively.
%   field: Particular command of interest. 
%          If specified, the help section from README.driver 
%          is returned in sh.  
%   
%Output Parameters:
%   l a cell string of the possible parameters
%   sh: help string corresponding to field
%   
%      
%Key Terms:
%   
%More Info :
%    TellAfni
%    TellAfni_Commands
%    Test_TellAfni
%    AFNI's README.driver
%    NewCs
%
%     Author : Ziad Saad
%     Date : Wed Dec 7 15:34:57 EST 2005
%     SSCC/NIMH/ National Institutes of Health, Bethesda Maryland


%Define the function name for easy referencing
FuncName = 'TellAfni_Commands';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

if (nargin == 0), %try to find README.driver
   f = '';
   fpart = '';
elseif (nargin == 1),
   fpart = '';
end

if (isempty(f)),
   [s,w] = unix('locate README.driver');
   if (~isempty(w)), %got a hit
      i = find(isspace(w));
      f = w(1:i(1)-1);
   end
end

if (~isempty(f)),
   if (exist(f,'file') == 2),
      fprintf(1,'Using file %s for info.\n', f);
   else
      fprintf(2,'Failed to find README file %s\n', f);
      return;
   end
else 
   f = '';
end

if (isempty(f)),
   [f, p] = uigetfile('README.driver', 'Pick a README.driver');
   if (f == 0),
      return;
   end
   f = sprintf('%s%c%s',p,filesep,f);
end

fid = fopen(f, 'r');
if (fid < 0),
   fprintf(2,'Failed to read %s\n', f);
   return;
end
   
s = fscanf(fid,'%c');
ns = length(s);
fclose(fid);

[err,inl] = FindChar(s,'NewLine');
ninl = length(inl);

ncm = 50;
cnt = 1;
S = char(zeros(ninl, ncm));
sl = [];
if (ninl),
   for (i=1:1:ninl-1),
      [err,w] = GetWord(s(inl(i):min(inl(i)+ncm, ns)),1);
      w = zdeblank(w);
      nw = length(w);
      if (nw > 1),
         if ( (sum(abs(double(upper(w))-double(w))) == 0) & sum(isletter(w)) > nw/2),
            % fprintf(2,'>%s<\t',w); pause
            S(cnt,1:nw) = w;
            %store location in s
            sl(cnt) = inl(i);
            cnt = cnt + 1;
         end
      end
   end
else
   fprintf(2,'Failed to understand file.\n');
   return;
end

if (cnt < 2),
   fprintf(2,'Failed to understand file, found nothing.\n');
   return;
end

[l, ilu] = unique(cellstr (S(1:cnt-1,:)));
ils = sl(ilu);
 
%do they want particular info ?
sh = '';
if (~isempty(fpart)),
   fpart = upper(zdeblank(fpart));
   cnt = 1;
   fnd = 0;
   while (cnt <= length(l)),
      if (~isempty(findstr(char(l(cnt)),fpart))), %(strncmp(char(l(cnt)),fpart,length(fpart))),
         fnd = cnt;
         %find location of next command
         dils = (ils - ils(cnt)); dils(find(dils<=0)) = max(dils + 1); [jj,imin] = min(dils); 
         %here is the string
         sh = [sh s(ils(cnt)+1:min(ils(imin), length(s)))];
      end
      cnt = cnt + 1;
   end
   if (~fnd),
      fprintf(2,'Failed to find help for %s.\n', fpart);
      return;
   end   
end

return;

