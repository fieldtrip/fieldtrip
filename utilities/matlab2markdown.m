function matlab2markdown(infile, outfile, varargin)

% MATLAB2MARKDOWN converts a MATLAB script or function to Markdown format. All
% comments are converted to text, comment lines starting with %% are converted to
% headings, sections with code are properly formatted, and words that appear in bold,
% italic or monospace are converted.
%
% Use as
%   matlab2markdown(infile, outfile, ...)
%
% If no outfile is specified, it will write it to a .md file with the same name as
% the infile. In case the file exists, it will be written with a numeric suffix.
%
% The best is to provide the full filepath, otherwise it will look for the file within
% the current path.
%
% Optional input arguments can be specified as key-value pairs and can include
%   imagestyle = 'none|inline|jekyll'
%   pageheader = 'none|jekyll'
%   overwrite  = true/false, allow overwriting of the .md file (default = false)
%   highlight  = string, 'matlab', 'plaintext' or '' (default = '')
%   ...
%
% See also MARKDOWN2MATLAB

% Copyright (C) 2018, Sophie Arana and Robert Oostenveld
%
% This file is part of FieldTrip.
%
% megconnectome is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% megconnectome is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with megconnectome.  If not, see <http://www.gnu.org/licenses/>.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following links are useful during development
% https://regexr.com
% https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin>2 && ismember(outfile, {'imagestyle', 'pageheader', 'pagelayout', 'pagetitle', 'pagetags'})
  varargin = [outfile varargin];
  outfile = [];
end

% parse the optional input arguments
% val = ft_getopt(varargin, 'key', default);
imagestyle    = ft_getopt(varargin, 'imagestyle', 'inline');
pageheader    = ft_getopt(varargin, 'pageheader', 'none');
pagelayout    = ft_getopt(varargin, 'pagelayout', '');
pagetitle     = ft_getopt(varargin, 'pagetitle', '');
pagetags      = ft_getopt(varargin, 'pagetags', '');
monospacehelp = ft_getopt(varargin, 'monospacehelp', false); % convert the help in monospace format
overwrite     = ft_getopt(varargin, 'overwrite', false);
highlight     = ft_getopt(varargin, 'highlight', '');

% check input
[inpath, inname, inext] = fileparts(infile);
if isempty(inpath), [inpath, inname, inext] = fileparts(which([inname inext])); end
if ~strcmp(inext,'.m')
  error('please specify a MATLAB file')
end

if nargin < 2 || isempty(outfile)
  outfile = infile;
end

if isempty(highlight)
  highlight = '';
end

% check output
[outpath, outname,outext] = fileparts(outfile);
if isempty(outpath), outpath = inpath; end
if ~strcmp(outext,'.md')
  outext = '.md';
end

% add a suffix to avoid overwriting files
if ~overwrite
  suffix = 1;
  newname = outname;
  while exist(fullfile(outpath,[newname,outext]),'file') == 2
    newname = [outname sprintf('-%0d',suffix)];
    suffix = suffix+1;
  end
  if ~strcmp(outname, newname)
    warning('writing to %s', newname);
  end
  outname = newname;
end

infile = fullfile(inpath,[inname,inext]);
outfile = fullfile(outpath,[outname,outext]);

% read & convert file line by line
infid = fopen(infile,'r');
outfid = fopen(outfile,'w');

if strcmp(pageheader, 'jekyll')
  % print a Jekyll header at the top of the page, like this
  % ---
  % title: Tutorial on EEG analysis
  % layout: default
  % tags: eeg, tutorial
  % ---
  fprintf(outfid, '---\n');
  if ~isempty(pagetitle),   fprintf(outfid, 'title: %s\n', pagetitle); end
  if ~isempty(pagelayout),  fprintf(outfid, 'layout: %s\n', pagelayout); end
  if ~isempty(pagetags),    fprintf(outfid, 'tags: %s\n', pagetags); end
  fprintf(outfid, '---\n');
end

state = 'unknown';    % keeps track of multi-line formatting decisions
index = 0;            % this keeps count of the items in an ordered list
linenumber = 0;       % keep track of the line number

while ~feof(infid)
  
  % reset the index, unless explicitly specified
  reset_index = true;
  % reset the state, unless explicitly specified
  reset_state = true;
  
  line = fgetl(infid);
  linenumber = linenumber + 1;
  if ~ischar(line), break, end
  
  if monospacehelp
    if match(line, '^%') && strcmp(state, 'unknown')
      % the help has started
      state = 'help';
      reset_state = false;
      fprintf(outfid, '```%s\n', highlight);  % start of code block
      fprintf(outfid, '%s\n', formathelp(line));
    elseif match(line, '^%') && strcmp(state, 'help')
      % the help is continuing
      state = 'help';
      reset_state = false;
      fprintf(outfid, '%s\n', formathelp(line));
    elseif ~match(line, '^%') && strcmp(state, 'help')
      % the help has ended
      fprintf(outfid, '```\n'); % end of code block
      state = 'ignore';
      reset_state = false;
    elseif strcmp(state, 'ignore')
      % ignore this line
      state = 'ignore';
      reset_state = false;
    end
    
  elseif isempty(line)
    % keep the complete line as it is
    fprintf(outfid, '%s\n', line);
    
  elseif match(line, '^%%')
    % heading
    [dum,endIndex] = regexp(line, '^%% ');
    remainder = reformat(line((endIndex+1):end));
    fprintf(outfid, '# %s\n', remainder);
    
  elseif match(line, '^% *\* ')
    % unordered list
    [dum,endIndex] = regexp(line, '^% *\* ');
    [dum,level] = regexp(line, '^% *');
    if isempty(level)
      level = 0;
    else
      level = level - 1;
    end
    level = repmat(' ', 1, level);
    remainder = reformat(line((endIndex+1):end));
    fprintf(outfid, '%s- %s\n', level, remainder);
    
  elseif match(line, '^% *# ')
    % ordered list
    [dum,endIndex] = regexp(line, '^% *# ');
    remainder = reformat(line((endIndex+1):end));
    index = index + 1;
    fprintf(outfid, '%d. %s\n', index, remainder);
    % assume that the next line will be part of the ordered list
    reset_index = false;
    
  elseif match(line, '^% *')
    % normal text
    [dum,endIndex] = regexp(line, '^% *');
    remainder = reformat(line((endIndex+1):end));
    fprintf(outfid, '%s\n', remainder);
    
  elseif match(line, '^print.*png')
    % figure, print the code
    fprintf(outfid, '    %s\n', line);
    [startIndex,endIndex] = regexp(line, '[!-~]*.png$');
    image = line(startIndex:endIndex);
    [dum, f] = fileparts(image);
    if strcmp(imagestyle, 'inline')
      % include a default Markdown inline image
      fprintf(outfid, '![%s](%s)\n', f, image);
    elseif strcmp(imagestyle, 'jekyll')
      % include an image as a Jekyll include
      fprintf(outfid, '{%% include image src="%s" %%}\n', image);
    end
    
  elseif ~match(line, '^%')
    % normal code
    fprintf(outfid, '    %s\n', line);
    
  end % if the line matches some pattern
  
  if reset_index
    index = 0;
  end
  
  if reset_state
    state = 'unknown';
  end
  
end % while not feof

fclose(infid);
fclose(outfid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = match(str, pattern)
status = ~isempty(regexp(str, pattern, 'once'));

function str = reformat(str)
str = formatbold(str);
str = formatitalic(str);
str = formatunderline(str);
str = formatmonospace(str);

function str = formathelp(str)
% remove the first %
str = str(2:end);

function str = formatbold(str)
[startIndex,endIndex] = regexp(str, '\*[a-zA-Z0-9]*\*');
for i=1:length(startIndex)
  part1 = str(1:startIndex(i)-1);
  part2 = str(startIndex(i):endIndex(i));
  part3 = str((endIndex(i)+1):end);
  str = [part1 '**' part2(2:end-1) '**' part3];
  % correct the subsequent indices
  startIndex(i+1:end) = startIndex(i+1:end)+2;
  endIndex(i+1:end) = endIndex(i+1:end)+2;
end

function str = formatitalic(str)
% the format is actually the same
[startIndex,endIndex] = regexp(str, '_[a-zA-Z0-9]*_');
for i=1:length(startIndex)
  part1 = str(1:startIndex(i)-1);
  part2 = str(startIndex(i):endIndex(i));
  part3 = str((endIndex(i)+1):end);
  str = [part1 '_' part2(2:end-1) '_' part3];
end

function str = formatunderline(str)
% this is not supported in the conversion from MATLAB live scripts to normal MATLAB code.

function str = formatmonospace(str)
[startIndex,endIndex] = regexp(str, '\|[a-zA-Z0-9]*\|');
for i=1:length(startIndex)
  part1 = str(1:startIndex(i)-1);
  part2 = str(startIndex(i):endIndex(i));
  part3 = str((endIndex(i)+1):end);
  str = [part1 '`' part2(2:end-1) '`' part3];
  % correct the subsequent indices
  startIndex(i+1:end) = startIndex(i+1:end);
  endIndex(i+1:end) = endIndex(i+1:end);
end
