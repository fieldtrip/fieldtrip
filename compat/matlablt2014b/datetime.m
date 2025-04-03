function out = datetime(in)

% prior to 2014b datetime does not exist, this function uses the unrecommended date function to return a char array with
% today's date provided the input argument is a char that says 'today', otherwise an error is thrown

if ~ischar(in) || ~isequal(in, 'today')
  error('incorrect input');
else
  out = date;
end

