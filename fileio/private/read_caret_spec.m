function [spec, headerinfo] = read_caret_spec(specfile)

% READ_CARET_SPEC reads in a caret .spec file.
%
% Use as
%   [spec, headerinfo] = read_caret_spec(specfile)
%
% Output arguments:
%   spec       = structure containing per file type the files listed
%   headerinfo = structure containing the specfile header
%
% The file can be an xml-file or an ascii formatted file

% Copyright (C) 2013, Jan-Mathijs Schoffelen

try,
  % read xml-file that contains a description to a bunch of files
  % belonging together
  ft_hastoolbox('gifti', 1);
  g = xmltree(specfile);
  
  % convert into a structure
  s = convert(g);
  
  if isfield(s, 'FileHeader')
    headerinfo = s.FileHeader;
    spec       = rmfield(s, 'FileHeader');
  else
    headerinfo = [];
    spec       = s;
  end
  
  % process the headerinfo
  if ~isempty(headerinfo)
    if isfield(headerinfo, 'Element')
      tmp = headerinfo.Element;
      tmp2 = struct([]);
      for k = 1:numel(headerinfo.Element)
        tmp2(1).(strrep(headerinfo.Element{k}.Name, '-', '_')) = headerinfo.Element{k}.Value;
      end
      headerinfo = tmp2;
    end
  end
  
  % further process the fields in spec
  f = fieldnames(spec);

  for k = 1:numel(f)
    if isempty(strfind(f{k}, 'study_metadata'))
      if iscell(spec.(f{k}))
        tmp = spec.(f{k});
        tmp2 = {};
        for m = 1:numel(tmp)
          tmpx = tmp{m};
          if isstruct(tmpx)
            fn = fieldnames(tmpx)
            for i = 1:numel(fn)
              tmp2{end+1,1} = tmpx.(fn{i});
            end
          end
        end
        spec.(f{k}) = tmp2;
      elseif isstruct(spec.(f{k}))
        tmp = spec.(f{k});
        fn  = fieldnames(tmp);
        tmp2 = {};
        for m = 1:numel(fn)
          tmp2{end+1,1} = tmp.(fn{m});
        end
        spec.(f{k}) = tmp2;
      else
        % don't know what to do with it
        spec = rmfield(spec, f{k});
      end
    else
      % don't know what to do with it
      spec = rmfield(spec, f{k});
    end
  end
      
catch
  
  % process as ASCII-file
  fid  = fopen_or_error(specfile);
  line = 'some text';
  while isempty(strfind(line, 'EndHeader'))
    line = fgetl(fid);
    if isempty(strfind(line, 'BeginHeader')) && isempty(strfind(line, 'EndHeader'))
      tok = tokenize(line, ' ');
      headerinfo.(strrep(tok{1},'-','_')) = tok{2};
    end
  end
  line = fgetl(fid); % empty line
  
  spec = struct([]);
  while 1
    line = fgetl(fid);
    if ~ischar(line), break, end
    tok = tokenize(line, ' ');
    if ~isempty(tok{1})
    if isfield(spec, tok{1})
      spec(1).(tok{1}){end+1,1} = tok{2};
    else
      spec(1).(tok{1}){1} = tok{2};
    end
    end
  end
  fclose(fid);
end

