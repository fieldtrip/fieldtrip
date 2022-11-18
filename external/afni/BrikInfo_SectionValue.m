function [err,v, typ] = BrikInfo_SectionValue (BRIKinfo, sSection)
%
%   [err,s, typ] = BrikInfo_SectionValue (sHead, sSection)
%
%Purpose:
%   Get the values of a section in the .HEAD string
%   This function is dedicated for BrikInfo function
%
%Input Parameters:
%   sHead, a string vector containing the .HEAD info
%   sSection, the section name such as: 'DATASET_DIMENSIONS'
%
%
%Output Parameters:
%   err : 0 No Problem
%       : 1 Mucho Problems
%   v a vector (or a string) containing the values in that section
%     v is empty (and err = 1) if the section can't be found
%   typ : is a string to indicate the type of parameter (string or number)
%      that v is
%
%Key Terms:
%
%More Info :
%  see also BrikInfo
%
%
%
%     Author : Ziad Saad
%     Date : Mon Oct 18 13:46:28 CDT 1999


%Define the function name for easy referencing
FuncName = 'BrikInfo_SectionValue';

%Debug Flag
DBG = 1;

%initailize return variables
err = 1;

N_BRIKinfo = length(BRIKinfo);

%itmp = findstr (BRIKinfo, sSection); This fails if sSection appears in places
%like History note and history note appears first.
%Nick OOsterhof suggested the following fix:
htmp = regexp (BRIKinfo, ['name\s*=\s*' sSection]);
if isempty(htmp) | length(htmp) > 1, v = []; return; end
itmp = findstr(BRIKinfo(htmp:N_BRIKinfo),sSection)+htmp-1;
if isempty(itmp), v = []; return; end

%findout where this section ends
inxt = strfind(BRIKinfo(itmp:N_BRIKinfo),'type');
inxt = setdiff(inxt, strfind(BRIKinfo(itmp:N_BRIKinfo),'_type')+1); %JM added to avoid getting stuck on data_type
if isempty(inxt), inxt = N_BRIKinfo; else inxt = inxt(1) + itmp -2; end

%Get the count value
[err,sn, jnk, strt,stp] = GetNextLine(BRIKinfo(itmp:inxt),2);
ic = findstr(sn,'count'); %watch it, spacing sensitive section
if ~isempty(ic)
  %skip the = sign
  sn = sn(ic+5:length(sn));
  ic = findstr(sn, '=');
  sn = sn(ic+1:length(sn));
else
  err = ErrEval(FuncName,'Err_could not find count field');
  return;
end
%get the count value
n = str2num(sn);
if ~n %empty field
  v = [];
  typ = '';
  return;
end

%advance to the next line, after count
itmp = itmp + stp + 1;

%read all values
Svals = BRIKinfo(itmp:inxt); %that's where the values are
inl = find(char(Svals) == 10 | char(Svals) == 13); %find any new lines left in Svals (13 was added to work in DOS (CR and LF))

if (~isempty(inl))
  Svals(inl) = ' '; %replace them by space characters
end

%if (sum(isletter(Svals))), %isletter fails when you have numbers such as 1.4e-3
v = zdeblank(Svals);
if isempty(v)
  typ = 'string';
elseif strcmp(v(1),'''') && strcmp(v(length(v)),'~')
  %This parameter is a string, do not change to numbers
  %remove first and last chars
  v = v(2:length(v)-1);
  typ = 'string';
else
  v = str2num(Svals);
  typ = 'number';
end

%pause

err = 0;
return;

