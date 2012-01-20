function db_insert(tablename, s)

% DB_INSERT inserts a structure into a database table. Each field of
% the structure should correspond with one of the fields in the table.
%
% Use as
%   db_insert(tablename, s)

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: db_insert.m,v $
% Revision 1.2  2007/12/03 12:35:24  roboos
% fixed bug when input is double([])
%
% Revision 1.1  2007/11/07 10:50:46  roboos
% created helper functions for easy access to a MySQL database table using a structure for representing the data
%

if isempty(s)
  % the input consists of an empty structure. so nothing to do
  return
end

f  = fieldnames(s);
v = {};
t = {};
l = length(f);
for i=1:l
  v{i} = getfield(s, f{i});
  t{i} = class(v{i});
end

for j=1:numel(s)
  % loop over all elements of a structure array

  cmd = ['INSERT INTO ' tablename ' ( '];

  % these are the fieldnames
  for i=1:l
    cmd = [cmd f{i} ', '];
  end
  cmd = cmd(1:(end-2));  % remove the last ', '
  cmd = [cmd ' ) VALUES ( '];

  % these are the corresponding values
  for i=1:l
    switch t{i}
      case 'char'
        if isempty(v{i})
          cmd = [cmd '""' ];
        else
          cmd = [cmd '"' v{i} '"' ];
        end

      case {'int8', 'int16', 'int32', 'uint8', 'uint16', 'uint32', 'double'}
        if isempty(v{i})
          cmd = [cmd '""' ];
        else
          cmd = [cmd num2str(v{i})];
        end

      otherwise
        error('unsuported data type in structure');
    end
    cmd = [cmd ', '];
  end

  cmd = cmd(1:(end-2));  % remove the last ', '
  cmd = [cmd ' ) ;'];

  % execute the query
  mysql(cmd);

end %for numel(s)

