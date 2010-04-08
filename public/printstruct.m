%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = printstruct(name, val)
str = [];
if isstruct(val)
  if numel(val)>1
    % this function cannot print struct arrays with multiple elements
    str = sprintf('%s = ''FIXME'';', name);
    return
  else
    % print it as a named structure
    fn = fieldnames(val);
    for i=1:length(fn)
      fv = getfield(val, fn{i});
      switch class(fv)
        case 'char'
          % line = sprintf('%s = ''%s'';\n', fn{i}, fv);
          % line = [name '.' line];
          line = printstr([name '.' fn{i}], fv);
          str  = [str line];
        case {'single' 'double'}
          line = printmat([name '.' fn{i}], fv);
          str  = [str line];
        case {'int8' 'int16' 'int32' 'int64' 'uint8' 'uint16' 'uint32' 'uint64'}
          line = printmat([name '.' fn{i}], fv);
          str  = [str line];
        case 'cell'
          line = printcell([name '.' fn{i}], fv);
          str  = [str line];
        case 'struct'
          line = printstruct([name '.' fn{i}], fv);
          str  = [str line];
        otherwise
          error('unsupported');
      end
    end
  end
elseif ~isstruct(val)
  % print it as a named variable
  switch class(val)
    case 'char'
      str = printstr(name, val);
    case 'double'
      str = printmat(name, val);
    case {'int8' 'int16' 'int32' 'int64' 'uint8' 'uint16' 'uint32' 'uint64'}
      str = printmat(name, val);
    case 'cell'
      str = printcell(name, val);
    otherwise
      error('unsupported');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = printcell(name, val)
str = [];
siz = size(val);
if isempty(val)
  str = sprintf('%s = {};\n', name);
  return;
end
for i=1:prod(siz)
  typ{i} = class(val{i});
end
for i=2:prod(siz)
  if ~strcmp(typ{i}, typ{1})
    warning('different elements in cell array');
    return
  end
end
if all(size(val)==1)
  str = sprintf('%s = { %s };\n', name, printval(val{1}));
else
  str = sprintf('%s = {\n', name);
  for i=1:siz(1)
    dum = '';
    for j=1:siz(2)
      dum = [dum ' ' printval(val{i,j})];
    end
    str = sprintf('%s%s\n', str, dum);
  end
  str = sprintf('%s};\n', str);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = printmat(name, val)
str = [];
siz = size(val);
if any(size(val)==0)
  str = sprintf('%s = [];\n', name);
elseif all(size(val)==1)
  str = sprintf('%s = %s;\n', name, printval(val));
elseif size(val,1)==1
  dum = sprintf('%g ', str, val(:));
  str = sprintf('%s = [%s];\n', name, dum);
else
  str = sprintf('%s = [\n', name);
  for i=1:siz(1)
    dum = sprintf('%g ', val(i,:));
    str = sprintf('%s  %s\n', str, dum);
  end
  str = sprintf('%s];\n', str);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = printstr(name, val)
str = [];
siz = size(val);
if siz(1)~=1
  str = sprintf('%s = \n', name);
  for i=1:siz(1)
    str = [str sprintf('  %s\n', printval(val(i,:)))];
  end
else
  str = sprintf('%s = %s;\n', name, printval(val));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = printval(val)
str = [];
switch class(val)
  case 'char'
    str = sprintf('''%s''', val);
  case 'double'
    str = sprintf('%g', val);
  case {'int8' 'int16' 'int32' 'int64' 'uint8' 'uint16' 'uint32' 'uint64'}
    str = sprintf('%d', val);
  case 'struct'
    str = '''FIXME''';
end

