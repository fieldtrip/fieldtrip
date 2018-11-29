function markdown2matlab(infile,outfile,varargin)

% MARKDOWN2MATLAB converts a MarkDown file to a MATLAB script or function. All text
% is converted to comments, headings are converted to comment lines starting with %%
% sections with code are properly formatted, and words that appear in bold, italic or
% monospace are converted.
%
% Use as
%   markdown2matlab(infile, outfile)
%
% If no outfile is specified, it will write it to a .m file with the same name as
% the infile. In case the file exists, it will be written with a numeric suffix.
%
% The best is to provide the full filepath, otherwise it will look for the file within
% the current path.
%
% See also MATLAB2MARKDOWN

% Copyright (C) 2018 Sophie Arana and Robert Oostenveld
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
% You can use https://regexr.com to test/debug regular expressions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check input
[inpath, inname, inext] = fileparts(infile);
if isempty(inpath), [inpath, inname, inext] = fileparts(which([inname inext])); end
if ~strcmp(inext,'.md')
  error('please specify a MarkDown file')
end

if nargin < 2
  outfile = infile;
end

% check output
[outpath, outname,outext] = fileparts(outfile);
if isempty(outpath), outpath = inpath; end
if ~strcmp(outext,'.m')
  outext = '.m';
end

% add a suffix to avoid overwriting files
suffix = 1;
newname = outname;
while exist(fullfile(outpath,[newname,outext]),'file') == 2
  newname = [outname sprintf('-%0d',suffix)];
  suffix = suffix+1;
end
outname = newname;

infile = fullfile(inpath,[inname,inext]);
outfile = fullfile(outpath,[outname,outext]);

% read & convert file line by line
infid = fopen(infile,'r');
outfid = fopen(outfile,'w');

format = 'comment';

while ~feof(infid)
  
  % the format will be reset unless otherwise specified
  reset_format = true;
  
  line = fgetl(infid);
  if ~ischar(line), break, end
  
  if match(line, '^ *[-+*] ')
    % unordered list
    [~,endIndex] = regexp(line, '^ *[-+*] ');
    
    [~,level] = regexp(line, '^ *');
    if isempty(level)
      level = 0;
    end
    level = repmat(' ', 1, level);
    remainder = reformat(line((endIndex+1):end));
    fprintf(outfid, '%%%s* %s\n', level, remainder);
    
  elseif match(line, '^    ')
    % normal code
    [~,endIndex] = regexp(line, '^    ');
    remainder = line((endIndex+1):end);
    fprintf(outfid, '%s\n', remainder);
    format = 'code';
    % assume that the next line will be in the same format
    reset_format = false;
    
  elseif match(line, '^[A-Za-z]')
    % normal text
    remainder = reformat(line);
    fprintf(outfid, '%% %s\n', remainder);
    
  elseif match(line, '^ *[0-9]*\.')
    % ordered list
    [~,endIndex] = regexp(line, '^ *[0-9]*\. ');
    remainder = reformat(line((endIndex+1):end));
    fprintf(outfid, '%% # %s\n', remainder);
    
  elseif match(line, '^#')
    % heading
    [~,endIndex] = regexp(line, '^# *');
    remainder = reformat(line((endIndex+1):end));
    fprintf(outfid, '%%%% %s\n', remainder);
    
  elseif isempty(line)
    % keep the complete line as it is
    if strcmp(format, 'comment')
      fprintf(outfid, '%%\n');
    else
      fprintf(outfid, '\n');
    end
    
  end
  
  if reset_format
    format = 'comment';
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

function str = formatbold(str)
[startIndex,endIndex] = regexp(str, '\*\*[a-zA-Z0-9]*\*\*');
for i=1:length(startIndex)
  part1 = str(1:startIndex(i)-1);
  part2 = str(startIndex(i):endIndex(i));
  part3 = str((endIndex(i)+1):end);
  str = [part1 '*' part2(3:end-2) '*' part3];
  % correct the subsequent indices
  startIndex(i+1:end) = startIndex(i+1:end)-2;
  endIndex(i+1:end) = endIndex(i+1:end)-2;
end

function str = formatitalic(str)
% the format is actually the same
[startIndex,endIndex] = regexp(str, '_[a-z]*_');
for i=1:length(startIndex)
  part1 = str(1:startIndex(i)-1);
  part2 = str(startIndex(i):endIndex(i));
  part3 = str((endIndex(i)+1):end);
  str = [part1 '_' part2(2:end-1) '_' part3];
end

function str = formatunderline(str)
% this is not supported in the conversion from MATLAB live scripts to normal MATLAB code.

function str = formatmonospace(str)
[startIndex,endIndex] = regexp(str, '```[a-z]*```');
for i=1:length(startIndex)
  part1 = str(1:startIndex(i)-1);
  part2 = str(startIndex(i):endIndex(i));
  part3 = str((endIndex(i)+1):end);
  str = [part1 '|' part2(4:end-3) '|' part3];
  % correct the subsequent indices
  startIndex(i+1:end) = startIndex(i+1:end)-4;
  endIndex(i+1:end) = endIndex(i+1:end)-4;
end
