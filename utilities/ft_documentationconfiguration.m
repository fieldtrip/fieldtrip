function [configuration] = ft_documentationconfiguration(filename)

% FT_DOCUMENTATIONCONFIGURATION is a helper function to maintain the online
% documentation of all configuration options.
%
% Normal users will not be calling this function, but will rather look at
% http://www.fieldtriptoolboxorg/reference/configuration where the output of this
% function can be found.
%
% See also FT_DOCUMENTATIONREFERENCE

% Copyright (C) 2008-2019, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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

subdir = {
  '.'
  'contrib/misc'
  'contrib/nutmegtrip'
  'contrib/spike'
  'external/artinis'
  'realtime/example'
  'realtime/online_eeg'
  'realtime/online_meg'
  'realtime/online_mri'
  };

funname = {};

% find all functions that should be included in the reference documentation
for i=1:length(subdir)
  f = dir(fullfile(ftpath, subdir{i}, 'ft_*.m'));
  f = {f.name}';
  funname = cat(1, funname, f);
end

for j=1:length(funname)
  [p, f, x] = fileparts(funname{j});
  funname{j} = f;
end

ncfg  = 0;
configuration = {};

for j=1:length(funname)
  str = help(funname{j});
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
      str{i} = strrep(str{i}, ' ', '');
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
      configuration{ncfg,1} = funname{j};
      dum = regexp(str{i}, 'cfg.[a-zA-Z0-9_\.]*', 'match');
      configuration{ncfg,2} = dum{1};
      dum = str{i};
      while length(dum)>0 && dum(1)~=' '
        dum = dum(2:end);
      end
      while length(dum)>0 && (dum(1)=='=' || dum(1)==' ')
        dum = dum(2:end);
      end
      configuration{ncfg,3} = dum;
      dum1 = configuration{ncfg,1};
      dum1(end+1:30) = ' ';
      dum2 = configuration{ncfg,2};
      dum2(end+1:30) = ' ';
    end
  end
end

% add links to reference doc
for i=1:size(configuration,1)
  configuration{i,1} = sprintf('[%s](/reference/%s)', configuration{i,1}, configuration{i,1});
end

configuration = sortrows(configuration(:,[2 3 1]));
configuration = configuration(:, [3 1 2]);
count = 0;
for i=2:size(configuration,1)
  prevfun = configuration{i-1,1};
  prevcfg = configuration{i-1,2};
  prevcmt = configuration{i-1,3};
  thisfun = configuration{i,1};
  thiscfg = configuration{i,2};
  thiscmt = configuration{i,3};

  if strcmp(thiscfg,prevcfg) && strcmp(thiscmt,prevcmt)
    count = count + 1;
    thisfun = [prevfun ', ' thisfun];
    prevfun = '';
    prevcfg = '';
    prevcmt = '';
    configuration{i  ,1} = thisfun;
    configuration{i-1,1} = prevfun;
    configuration{i-1,2} = prevcfg;
    configuration{i-1,3} = prevcmt;
  end
end
fprintf('merged %d cfg options\n', count);

fid = fopen(filename, 'wb');
currletter = char(96);

fprintf(fid, '---\n');
fprintf(fid, 'title: Index of all configuration options\n');
fprintf(fid, '---\n');
fprintf(fid, '\n');
fprintf(fid, '# Index of all configuration options \n');
fprintf(fid, '\n');
fprintf(fid, 'A detailed description of each function is available in the [reference documentation](/reference).\n');
fprintf(fid, '\n');

for i=1:size(configuration,1)
  if isempty(configuration{i,1})
    continue;
  elseif length(configuration{i,2})<5
    continue;
  end
  thisletter = configuration{i,2}(5);
  while currletter<thisletter
    currletter = currletter + 1;
    fprintf(fid, '## %s \n\n', upper(char(currletter)));
  end
  fprintf(fid, '**%s** - %s  \n', configuration{i,2}, configuration{i,1});

  % do postprocessing to make sure we don't mess up dokuwiki layout
  % '' is a markup instruction for dokuwiki so escape by replacing it
  % with %%''%%
  % configuration{i,3} = strrep(configuration{i,3},'''''','%%''''%%');

  fprintf(fid, '%s\n\n', configuration{i,3});
end
fclose(fid);

return
