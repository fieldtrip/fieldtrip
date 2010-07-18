function str = guidelines(funname, pat)

% GUIDELINES searches for a contiguous block of commented text and shows
% its contents. It is used to display additional help sections.

if nargin<2
  pat = '[Gg]uidelines';
end

if ispc
  linesep = [13 10];
else
  linesep = 10;
end

if nargout>0
  str = '';
end

funfile = which(funname);
fid = fopen(funfile, 'rt');

blockbeg = false;
blockend = false;

while ~feof(fid)
  line = fgetl(fid);

  if ~blockbeg
    % test whether the block begins
    blockbeg = ~isempty(regexp(line, pat, 'once'));
  else
    % test whether the block ends
    blockend = isempty(line) || ~isequal(line(1), '%');
  end

  if blockbeg && blockend
    break
  end

  if blockbeg && ~blockend
    if nargout>0
      str = cat(2, str, line(2:end), linesep);
    else
      disp(line(2:end));
    end
  end

end

fclose(fid);

