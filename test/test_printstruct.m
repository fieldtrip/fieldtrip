function test_printstruct

% MEM 2gb
% WALLTIME 00:20:00
% DEPENDENCY printstruct
% DATA no

% the above requirements are quite big, but this script is inherently unpredictable

numtests = 10;

fprintf('generating %d random deep nested structures to test printstruct() serialization\n', numtests);
for k = 1:numtests
  % generate some deep structure
  oldstruct = randomval('struct');

  % get the printed/string version with a new name
  printversion = printstruct('newstruct', oldstruct);

  % do not add it to the existing structure
  clear newstruct

  % eval() it
  eval(printversion);

  % check equality
  % use abstol here because we know all floating point numeric values are
  % generated from standard normal distribution
  [ok,msg] = isalmostequal(oldstruct, newstruct, 'abstol', 1e-6);
  if ok
    fprintf('printstruct() behaves as expected for random structure %d\n', k);
  else
    fprintf('printstruct() DOES NOT behave as expected for random structure %d\n', k);
    fprintf('%s\n', msg{:});
    error('test failed, see above for details');
  end
end

end % function test_printstruct

%%%%%%%%%%%%%%%%%%%% SUBFUNCTION %%%%%%%%%%%%%%%%%%%%
function myval = randomval(type, depth)

if nargin < 2
  depth = 0;
end

if depth > 3 && any(strcmp(type, {'struct' 'cell'}))
  % don't nest too far (unnecessarily slows down the test)
  myval = [];
  return;
end

% the 64-bit int/uint types are not supported by matlab's randi(), so we
% don't test them here
types = {'double' 'double_complex' 'single' 'int8' 'int16' 'int32' 'uint8' 'uint16' 'uint32' 'logical' 'struct' 'cell'};

switch(type)
  case 'struct'
    numfields = 5 + randi(10);
    for k = 1:numfields
      name = sprintf('x%d', randi(intmax));
      myval.(name) = randomval(types{randi(numel(types))}, depth + 1);
    end
  case 'cell'
    numelem = randi(200);
    alldepths = repmat({depth+1}, [1 numelem]);
    alltypes = types(randi(numel(types), 1, numelem));
    % structs and cells cannot be within a cell as far as printstruct is concerned
    alltypes(strcmp(alltypes, 'struct') | strcmp(alltypes, 'cell')) = {'double'};
    myval = cellfun(@randomval, alltypes, alldepths, 'uniformoutput', false);
  case 'logical'
    if rand() < 0.5
      myval = false(randi(100),randi(100));
    else
      myval = true(randi(100),randi(100));
    end
  case 'double_complex'
    siz = [randi(100), randi(100)];
    myval = randn(siz) + 1i .* randn(siz);
  case {'double' 'single'}
    myval = randn(randi(100), randi(100));
  case {'int8' 'int16' 'int32' 'uint8' 'uint16' 'uint32' 'uint64'}
    myval = randi(200, randi(100), randi(100), type);
end

end % function randomval