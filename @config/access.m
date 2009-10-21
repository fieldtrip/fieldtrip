function y = access(x, cmd);

% ACCESS Return the number of accesses (assignments and references) to a CONFIGURATION object.

if nargin==1
  fprintf('----------- values     -----------\n');
  disp(x.value);

  fprintf('----------- original   -----------\n');
  disp(x.original);

  fprintf('----------- assignment -----------\n');
  disp(x.assign);

  fprintf('----------- reference -----------\n');
  disp(x.reference);
else
  switch cmd(1)
    case 'v'
      y = x.value;
    case 'r'
      y = x.reference;
    case 'a'
      y = x.assign;
    case 'o'
      y = x.original;
    otherwise
      error('Incorrect command for accessing a config object.')
  end
end
