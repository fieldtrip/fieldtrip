function j = struct2json(s, pad)

% STRUCT2JSON formats a MATLAB structure as JSON document
%
% Please use http://jsonlint.com for testing

if nargin<2
  % it should start with 4 spaces for padding
  pad = '    ';
end

eol = '\n';

j = '';
j = [j sprintf('{')];

fn = fieldnames(s);
for i=1:numel(fn)
  key = fn{i};
  val = s.(fn{i});
  
  j = [j sprintf(eol)];
  
  if i<numel(fn)
    % all but the last elements are followed by a comma
    eol = ',\n';
  else
    % the last element is not followed by a comma
    eol = '\n';
  end
  
  j = [j sprintf(pad)];
  switch class(val)
    case 'char'
      assert(size(val,1)==1);
      j = [j sprintf('"%s": "%s"', key, val)];
      
    case {'single', 'double'}
      assert(numel(val)==1);
      j = [j sprintf('"%s": %s', key, num2str(val))];
      
    case 'struct'
      assert(numel(val)==1);
      j = [j sprintf('"%s": ', key)];
      j = [j struct2json(val, [pad '    '])];
      
    case 'logical'
      assert(numel(val)==1);
      if val
        j = [j sprintf('"%s": true', key)];
      else
        j = [j sprintf('"%s": false', key)];
      end
      
    otherwise
      error('unsupported class "%s"', class(val));
  end % switch
  
end % for all fieldnames

j = [j sprintf(eol)];

if isequal(pad, '    ')
  % back at the top level, no padding needed any more
  j = [j '}'];
else
  % not yet at the top level, align the closing bracket with the opening variable
  j = [j sprintf(pad(1:end-4)) '}'];
end




