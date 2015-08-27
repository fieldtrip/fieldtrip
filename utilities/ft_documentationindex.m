function index = ft_documentationindex(filename)

% FT_DOCUMENTATIONINDEX is a function to maintain the online documentation. 
%
% Normal users will not be calling this function, but will rather look at
% http://fieldtrip.fcdonders.nl/reference/index where the output of this
% function can be found.
%
% See FT_DOCUMENTATIONREFERENCE

% Copyright (C) 2008-2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

ft_defaults

[ftver, ftpath] = ft_version;

f1 = dir(fullfile(ftpath, 'ft_*.m'));
f1 = {f1.name}';

funname = f1;

for i=1:length(funname)
  [p, funname{i}, x] = fileparts(funname{i});
end

ncfg  = 0;
index = {};

for j=1:length(funname)
  str = help(funname{j})
  str = tokenize(str, 10);

  % compact the help description
  for i=2:(length(str)-1)
    prevline = str{i-1};
    thisline = str{i  };
    nextline = str{i+1};
    if length(thisline)<5 || length(prevline)<5
      continue
    end
    try
      if ~isempty(regexp(prevline, '^ *cfg')) && isempty(regexp(thisline, '^ *cfg')) && ~all(thisline(1:3)==' ')
        % do not concatenate, this line starts a new paragraph
      elseif ~isempty(regexp(prevline, '^ *cfg')) && all(thisline(1:5)==' ')
        % concatenate a multiline cfg description
        thisline = cat(2, prevline, thisline);
        prevline = '';
      elseif ~all(prevline(1:5)==' ') && ~all(thisline(1:5)==' ') && isempty(regexp(thisline, '^ *cfg'))
        % concatenate the lines of a paragraph
        thisline = cat(2, prevline, thisline);
        prevline = '';
      elseif isempty(regexp(prevline, '^ *cfg')) && ~isempty(regexp(thisline, '^  cfg'))
       % previous line is a paragraph, this line starts with "cfg" but has no extra space in front of it
       % so assume that the cfg is part of the running text in the paragraph and conactenate the lines
       thisline = cat(2, prevline, thisline);
       prevline = '';
      end
    catch
      disp(lasterr);
      disp(thisline);
    end
    str{i-1} = prevline;
    str{i  } = thisline;
    str{i+1} = nextline;
  end
  for i=1:length(str)
    if length(str{i})>1
      % remove double spaces
      dum = findstr(str{i}, '  ');
      str{i}(dum) = [];
    end
    while ~isempty(str{i}) && str{i}(1)==' '
      % remove spaces at the begin of the line
      str{i}(1) = [];
    end
    while ~isempty(str{i}) && str{i}(end)==' '
      % remove spaces at the end of the line
      str{i}(1) = [];
    end
  end
  for i=1:length(str)
    if ~isempty(regexp(str{i}, '^ *cfg.[a-zA-Z0-9_\.]*'))
      ncfg = ncfg+1;
      index{ncfg,1} = funname{j};
      dum = regexp(str{i}, 'cfg.[a-zA-Z0-9_\.]*', 'match');
      index{ncfg,2} = dum{1};
      dum = str{i};
      while length(dum)>0 && dum(1)~=' '
        dum = dum(2:end);
      end
      while length(dum)>0 && (dum(1)=='=' || dum(1)==' ')
        dum = dum(2:end);
      end
      index{ncfg,3} = dum;
      dum1 = index{ncfg,1};
      dum1(end+1:30) = ' ';
      dum2 = index{ncfg,2};
      dum2(end+1:30) = ' ';
    end
  end
end

% add links to reference doc
for i=1:size(index,1)
  index{i,1} = sprintf('[[reference:%s]]', index{i,1});
end

index = sortrows(index(:,[2 3 1]));
index = index(:, [3 1 2]);
count = 0;
for i=2:size(index,1)
  prevfun = index{i-1,1};
  prevcfg = index{i-1,2};
  prevcmt = index{i-1,3};
  thisfun = index{i,1};
  thiscfg = index{i,2};
  thiscmt = index{i,3};

  if strcmp(thiscfg,prevcfg) && strcmp(thiscmt,prevcmt)
    count = count + 1;
    thisfun = [prevfun ', ' thisfun];
    prevfun = '';
    prevcfg = '';
    prevcmt = '';
    index{i  ,1} = thisfun;
    index{i-1,1} = prevfun;
    index{i-1,2} = prevcfg;
    index{i-1,3} = prevcmt;
  end
end
fprintf('merged %d cfg options\n', count);

fid = fopen(filename, 'wb');
currletter = char(96);
fprintf(fid, '====== Index of configuration options ======\n\n');
fprintf(fid, 'A detailed description of each function is available in the [[:reference|reference documentation]].\n\n');
fprintf(fid, 'This index to the reference documentation is automatically generated from the MATLAB code every day. Therefore you should not edit this page manually, since your changes would be overwritten automatically. If you want to suggest corrections to the documentation, please send them by email to the mailing list or to one of the main developers (see [[:contact]]).\n\n');
for i=1:size(index,1)
  fprintf('%s -- %s -- %s\n', index{i,2}, index{i,3}, index{i,1});
  if isempty(index{i,1})
    continue;
  elseif length(index{i,2})<5
    continue;
  end
  thisletter = index{i,2}(5);
  while currletter<thisletter
    currletter = currletter + 1;
    fprintf(fid, '===== %s =====\n\n', upper(char(currletter)));
  end
  fprintf(fid, '** %s ** // %s //\\\\\n', index{i,2}, index{i,1});
  
  % do postprocessing to make sure we don't mess up dokuwiki layout
  % '' is a markup instruction for dokuwiki so escape by replacing it
  % with %%''%%
  index{i,3} = strrep(index{i,3},'''''','%%''''%%');
  
  fprintf(fid, '%s\n\n', index{i,3});
end
fclose(fid);

return

