function j = ft_struct2json(s)

%  FT_STRUCT2JSON

fn = fieldnames(s);
fv = cell(size(fn));
for i=1:numel(fn)
  val = s.(fn{i});
  switch class(val)
    case 'char'
      fv{i} = val;
    case 'double'
      fv{i} = num2str(val);
    otherwise
      error('class %s is not supported\n', type(val));
  end
end
f = cat(1, fn', fv');
j = sprintf('"%s": "%s", ', f{:});
j = j(1:end-2); % remove the last comma and space

