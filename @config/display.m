function display(x)

% DISPLAY Display function for a config object.

% this flag detemines how much feedback is given
% which is usefull for debugging
fb = false;

if ~isempty(inputname(1))
  fprintf('\n');
  fprintf('%s =\n\n', inputname(1));
end

if fb
  fprintf('----------- value      -----------\n');
end

if numel(x)==0
  disp('1x1 config array with no fields.');
elseif numel(x)==1
  disp(x.value);
else
  siz = sprintf('%dx', size(x)); % construct a string like 1x2x3x, the last 'x' has to be removed
  siz = siz(1:(end-1));
  fprintf('%s config array with fields:\n', siz);
  key = fieldnames(x);
  for i=1:length(key)
    fprintf('    %s\n', key{i});
  end
end

if fb
  fprintf('----------- assignment -----------\n');
  disp(x.assign);
  fprintf('----------- reference  -----------\n');
  disp(x.reference);
end

