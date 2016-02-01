function [sS, ok, msgs] = som_set(sS, varargin)

%SOM_SET Create and check SOM Toolbox structs, give values to their fields.
%
% [sS, ok, msgs] = som_set(sS, [field, contents, ...])
%
%   sM              = som_set(sM,'name','SOM#1.1');
%   [dummy,ok,msgs] = som_set(sData);   
%   sT              = som_set('som_topol','msize',[10 10],'lattice','hexa');
%   [sTrain,ok]     = som_set(sTrain,'algorithm','lininit');
%   [sN,ok,msgs]    = som_set('som_norm');
%
% Input and output arguments ([]'s are optional):
%  sS                   the target struct
%              (struct) a SOM Toolbox structure (not visualization struct)
%              (string) structure identifier (see below)
%                       the updated/created structure is returned
%  [field,     (string) field to be given value to (see below)
%   contents]  (varies) the contents for the field
%
%  ok          (vector)  status for each field-contents pair (1=ok)
%  msgs        (cellstr) status string for each field-contents pair (''=ok)
%
%  There can be arbitrarily many field-contents pairs. If there
%  are _no_ field-content pairs, and the first argument is a struct,
%  the fields of the struct are checked for validity.
% 
%  Valid struct and corresponding field identifiers: 
%  'som_map'   : 'codebook', 'labels', 'mask', 'neigh', 'name', 
%                'topol', 'msize, 'lattice', 'shape',
%                'trainhist', 'comp_names', 'comp_norm', 
%  'som_data'  : 'data', 'labels', 'name', 'comp_names', 'comp_norm', 
%                'label_names'
%  'som_topol' : 'msize', 'lattice', 'shape'
%  'som_norm'  : 'method', 'params', 'status'
%  'som_train' : 'algorithm', 'data_name', 'mask', 'neigh', 
%                'radius_ini', 'radius_fin', 'alpha_ini', 'alpha_type', 
%                'trainlen', 'time'
%  'som_grid'  : 'lattice', 'shape', 'msize', 'coord',
%                'line', 'linecolor', 'linewidth', 
%                'marker', 'markersize', 'markercolor', 'surf', 
%                'label', 'labelcolor', 'labelsize'
%                checking given values has not been implemented yet!
%                
% For more help, try 'type som_set' or check out online documentation.
% See also SOM_INFO, SOM_MAP_STRUCT, SOM_DATA_STRUCT, SOM_VS1TO2.

%%%%%%%%%%%%% DETAILED DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% som_set
%
% PURPOSE
%
% Create and set values for fields of SOM Toolbox structs (except
% visualization struct). Can also be used to check the validity of structs.
%
% SYNTAX
%
%  sMap   = som_set('som_map');
%  sData  = som_set(sData);
%  sNorm  = som_set(...,'field',contents,...);
%  [sTopol,ok]      = som_set(sTopol,...);
%  [sTrain,ok,msgs] = som_set('som_train',...);
%
% DESCRIPTION
%
% The function is used to create and set values for fields of SOM
% Toolbox structs, except visualization structs. The given values are
% first checked for validity, and if they are not valid, an error
% message is returned. The function can also be used to check the
% validity of all the fields of the struct by supplying a struct as
% the first and only argument.
% 
% NOTE: Using SOM_SET to create structures does _not_ guarantee that the
% structs are valid (try e.g. sM = som_set('som_map'); som_set(sM)). The
% initial values that the function gives to the fields of the structs are
% typically invalid. It is recommended that when creating map or data 
% structs, the corresponding functions SOM_MAP_STRUCT and SOM_DATA_STRUCT 
% are used instead of SOM_SET. However, when giving values for the fields, 
% SOM_SET tries to guarantee that the values are valid.
%
% If a string is given as the first argument, the corresponding 
% structure is first created and the field-content pairs are then
% applied to it. 
%
% There can be arbitrarily many field-contents pairs. The pairs
% are processed sequentially one pair at a time. For each pair,
% the validity of the contents is checked and the corresponding 
% items in the returned 'ok'-vector and 'msgs'-cellstring are set.
% - if the contents is ok, the status is set to 1 and message to ''
% - if the contents is suspicious, status is set to 1, but a
%   message is produced
% - if the contents is invalid, status is set to 0 and an error
%   message is produced. The contents are _not_ given to the field.
% If there is only one output argument, the status and messages
% for each pair are printed to standard output.
%
% The different field-contents pairs have no effect on each other.
% If a field is given a value multiple times, the last valid one 
% stays in effect.
% 
% In some cases, the order of the given fields is significant.
% For example in the case of 'som_map', the validity of some fields, 
% like '.comp_names', depends on the input space dimension, which is
% checked from the '.data' field (dim = size(sD.data,2) to be specific).
% Therefore, the '.data' field (or '.codebook' field in case of map 
% struct) should always be given a value first. Below is a list of 
% this kind of dependancies:
% 
% som_map:   'comp_names', 'comp_norm', 'msize', 'topol.msize',
%            'labels' and 'mask' depend on 'codebook'
%            new value for 'codebook' should have equal size to the old
%            one (unless the old one was empty)
% som_data:  'comp_names' and 'comp_norm' depend on 'data'
%            new value for 'data' should have equal dimension (size(data,2))
%            as the old one (unless the old one was empty)
% 
% KNOWN BUGS
%
% Checking the values given to som_grid struct has not been
% implemented. Use SOM_GRID function to give the values.
%
% REQUIRED INPUT ARGUMENTS
%
%  sS          The struct.
%     (struct) A SOM Toolbox struct.
%     (string) Identifier of a SOM Toolbox struct: 'som_map', 
%              'som_data', 'som_topol', 'som_norm' or 'som_train'
%   
% OPTIONAL INPUT ARGUMENTS 
%
%  field     (string) Field identifier string (see below).
%  contents  (varies) Value for the field (see below).
%
%  Below is the list of valid field identifiers for the different 
%  SOM Toolbox structs. 
%
%  'som_map' (map struct)
%    'codebook'    : matrix, size [munits, dim] 
%    'labels'      : cell array of strings, 
%                    size [munits, maximum_number_of_labels]
%    'topol'       : topology struct (prod(topol.msize)=munits)
%    'mask'        : vector, size [dim, 1]
%    'neigh'       : string ('gaussian' or 'cutgauss' or 'bubble' or 'ep')
%    'trainhist'   : struct array of train structs
%    'name'        : string
%    'comp_names'  : cellstr, size [dim, 1], e.g. {'c1','c2','c3'}
%    'comp_norm'   : cell array, size [dim, 1], of cell arrays 
%                    of normalization structs
%    Also the following can be used (although they are fields
%    of the topology struct)
%    'msize'       : vector (prod(msize)=munits)
%    'lattice'     : string ('rect' or 'hexa')
%    'shape'       : string ('sheet' or 'cyl' or 'toroid')
%
%  'som_data' (data struct)
%    'data'        : matrix, size [dlen, dim]
%    'name'        : string
%    'labels'      : cell array of strings, 
%                    size [dlen, m]
%    'comp_names'  : cellstr, size [dim, 1], e.g. {'c1','c2','c3'}
%    'comp_norm'   : cell array, size [dim, 1], of cell arrays 
%                    of normalization structs
%    'label_names' : cellstr, size [m, 1]
%
% 'som_topol' (topology struct)
%    'msize'       : vector
%    'lattice'     : string ('rect' or 'hexa')
%    'shape'       : string ('sheet' or 'cyl' or 'toroid')
%
% 'som_norm' (normalization struct)
%    'method'      : string
%    'params'      : varies
%    'status'      : string ('done' or 'undone' or 'uninit')
%
% 'som_train' (train struct)
%    'algorithm'   : string ('seq' or 'batch' or 'lininit' or 'randinit')
%    'data_name'   : string
%    'mask'        : vector, size [dim, 1]
%    'neigh'       : string ('gaussian' or 'cutgauss' or 'bubble' or 'ep')
%    'radius_ini'  : scalar
%    'radius_fin'  : scalar
%    'alpha_ini'   : scalar
%    'alpha_type'  : string ('linear' or 'inv' or 'power')
%    'trainlen'    : scalar
%    'time'        : string
%
% 'som_grid' (grid struct) : checking the values has not been implemented yet!
%    'lattice'     : string ('rect' or 'hexa') or 
%                    (sparce) matrix, size munits x munits
%    'shape'       : string ('sheet' or 'cyl' or 'toroid')
%    'msize'       : vector, size 1x2
%    'coord'       : matrix, size munits x 2 or munits x 3
%    'line'        : string (linespec, e.g. '-', or 'none')
%    'linecolor'   : RGB triple or string (colorspec, e.g. 'k') or 
%                    munits x munits x 3 (sparce) matrix or cell
%                    array of RGB triples 
%    'linewidth'   : scalar or munits x munits (sparce) matrix
%    'marker'      : string (markerspec, e.g. 'o', or 'none') or 
%                    munits x 1 cell or char array of these
%    'markersize'  : scalar or munits x 1 vector
%    'markercolor' : RGB triple or string (colorspec, e.g. 'k')
%    'surf'        : [], munits x 1 or munits x 3 matrix of RGB triples
%    'label'       : [] or munits x 1 char array or 
%                    munits x l cell array of strings 
%    'labelcolor'  : RGB triple or string (colorspec, e.g. 'g' or 'none')
%    'labelsize'   : scalar
%
% OUTPUT ARGUMENTS
% 
%  sS    (struct)  the created / updated struct
%  ok    (vector)  length = number of field-contents pairs, gives
%                  validity status for each pair (0=invalid, 1 otherwise)
%  msgs  (cellstr) length = number of field-contents pairs, gives
%                  error/warning message for each pair ('' if ok)
%
% EXAMPLES
%
% To create a struct:
%  sM  = som_set('som_map');
%  sD  = som_set('som_data');
%  sTo = som_set('som_topol');
%  sTr = som_set('som_train');
%  sN  = som_set('som_norm');
%  sG  = som_set('som_grid');
%
% To check the the contents of a struct: 
%  som_set(sS);
%  [dummy,ok]      = som_set(sS);
%  [dummy,ok,msgs] = som_set(sS);
%
% To give values to fields: 
%  sTo = som_set(sTo,'msize',[10 10],'lattice','hexa','shape','toroid');
%  sM  = som_set('som_map','codebook',rand(100,4),'topol',sTo);
%   
% SEE ALSO
% 
%  som_info         Prints information the given struct.
%  som_map_struct   Create map struct.
%  som_data_struct  Create data struct.
%  som_topol_struct Create topology struct.
%  som_train_struct Create training struct.
%  som_grid         Create and visualize grid struct.
%  som_vs1to2       Conversion from version 1.0 structs to 2.0.
%  som_vs2to1       Conversion from version 2.0 structs to 1.0.

% Copyright (c) 1999-2000 by the SOM toolbox programming team.
% http://www.cis.hut.fi/projects/somtoolbox/

% Version 2.0beta juuso 101199 130300

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create struct if necessary

if ischar(sS), 
  switch sS
   case 'som_map',
    sS=struct('type', 'som_map', ...
              'codebook', [], ...
              'topol', som_set('som_topol'), ...
              'labels', cell(1), ...
              'neigh', 'gaussian', ...
              'mask', [], ...
              'trainhist', cell(1), ...
              'name', '',...
              'comp_names', {''}, ...
              'comp_norm', cell(1));
   case 'som_data', 
    sS=struct('type', 'som_data', ...
              'data', [], ...
              'labels', cell(1), ...
              'name', '', ...
              'comp_names', {''}, ...
              'comp_norm', cell(1), ...
              'label_names', []);
   case 'som_topol',
    sS=struct('type', 'som_topol', ...
              'msize', 0, ...
              'lattice', 'hexa', ...
              'shape', 'sheet');
   case 'som_train',
    sS=struct('type', 'som_train', ...
              'algorithm', '', ...
              'data_name', '', ...
              'neigh', 'gaussian', ...
              'mask', [], ...
              'radius_ini', NaN, ...
              'radius_fin', NaN, ...
              'alpha_ini', NaN, ...
              'alpha_type', 'inv', ...
              'trainlen', NaN, ...
              'time', '');
   case 'som_norm',
    sS=struct('type', 'som_norm', ...
              'method', 'var', ...
              'params', [], ...
              'status', 'uninit');
   case 'som_grid', 
    sS=struct('type','som_grid',...
	      'lattice','hexa',...
	      'shape','sheet',...
	      'msize',[1 1],...
	      'coord',[],...
	      'line','-',...
	      'linecolor',[.9 .9 .9],...
	      'linewidth',0.5,...
	      'marker','o',...
	      'markersize',6,...
	      'markercolor','k',...
	      'surf',[],...
	      'label',[],...
	      'labelcolor','g',...
	      'labelsize',12);    
   otherwise
    ok=0; msgs = {['Unrecognized struct type: ' sS]}; sS = [];
    return;
  end  
  
elseif isstruct(sS) & length(varargin)==0, 
  
  % check all fields
  fields = fieldnames(sS);
  if ~any(strcmp('type',fields)), 
    error('The struct has no ''type'' field.');
  end
  k = 0;
  for i=1:length(fields), 
    contents = getfield(sS,fields{i});
    if ~strcmp(fields{i},'type'), 
      varargin{k+1} = fields{i};
      varargin{k+2} = contents;
      k = k + 2;
    else 
      if ~any(strcmp(contents, ...
        {'som_map','som_data','som_topol','som_train','som_norm'})), 
	error(['Unknown struct type: ' contents]);
      end
    end
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set field values

p = ceil(length(varargin)/2);
ok = ones(p,1);
msgs = cell(p,1);

for i=1:p, 
  field = varargin{2*i-1}; 
  content = varargin{2*i};
  msg = '';
  isok = 0;
  
  si = size(content);
  isscalar = (prod(si)==1);
  isvector = (sum(si>1)==1);
  isrowvector = (isvector & si(1)==1);
  if isnumeric(content), 
    iscomplete = all(~isnan(content(:))); 
    ispositive = all(content(:)>0); 
    isinteger  = all(content(:)==ceil(content(:)));
    isrgb = all(content(:)>=0 & content(:)<=1) & size(content,2)==3;
  end
  
  switch sS.type, 
   case 'som_map',
    [munits dim] = size(sS.codebook);
    switch field, 
     case 'codebook', 
      if ~isnumeric(content), 
	msg = '''codebook'' should be a numeric matrix'; 
      elseif size(content) ~= size(sS.codebook) & ~isempty(sS.codebook), 
	msg = 'New ''codebook'' must be equal in size to the old one.'; 
      elseif ~iscomplete, 
	msg = 'Map codebook must not contain NaN''s.'; 
      else
	sS.codebook = content; isok=1;
      end
     case 'labels', 
      if isempty(content), 
	sS.labels = cell(munits,1); isok = 1;
      elseif size(content,1) ~= munits, 
	msg = 'Length of labels array must be equal to the number of map units.';
      elseif ~iscell(content) & ~ischar(content), 
	msg = '''labels'' must be a string array or a cell array/matrix.';
      else
	isok = 1;
	if ischar(content), content = cellstr(content); 
	elseif ~iscellstr(content), 
	  l = prod(size(content));
	  for j=1:l, 
	    if ischar(content{j}), 
	      if ~isempty(content{j}), 
		msg = 'Invalid ''labels'' array.';
		isok = 0; 
		break;
	      else
		content{j} = ''; 
	      end
	    end
	  end
	end
	if isok, sS.labels = content; end
      end
     case 'topol', 
      if ~isstruct(content), 
	msg = '''topol'' should be a topology struct.'; 
      elseif ~isfield(content,'msize') | ...
	    ~isfield(content,'lattice') | ...
	    ~isfield(content,'shape'), 
	msg = '''topol'' is not a valid topology struct.'; 
      elseif prod(content.msize) ~= munits, 
	msg = '''topol''.msize does not match the number of map units.'; 
      else
	sS.topol = content; isok = 1;
      end
     case 'msize', 
      if ~isnumeric(content) | ~isvector | ~ispositive | ~isinteger, 
	msg = '''msize'' should be a vector with positive integer elements.'; 
      elseif prod(content) ~= munits, 
	msg = '''msize'' does not match the map size.'; 
      else
	sS.topol.msize = content; isok = 1;
      end
     case 'lattice', 
      if ~ischar(content),
	msg = '''lattice'' should be a string'; 
      elseif ~strcmp(content,'rect') & ~strcmp(content,'hexa'),
	msg = ['Unknown lattice type: ' content];
	sS.topol.lattice = content; isok = 1;
      else
	sS.topol.lattice = content; isok = 1;
      end
     case 'shape', 
      if ~ischar(content),
	msg = '''shape'' should be a string';
      elseif ~strcmp(content,'sheet') & ~strcmp(content,'cyl') & ...
	    ~strcmp(content,'toroid'),
	msg = ['Unknown shape type:' content]; 
	sS.topol.shape = content; isok = 1;
      else
	sS.topol.shape = content; isok = 1;
      end
     case 'neigh', 
      if ~ischar(content),
	msg = '''neigh'' should be a string'; 
      elseif ~strcmp(content,'gaussian') & ~strcmp(content,'ep') & ...
	    ~strcmp(content,'cutgauss') & ~strcmp(content,'bubble'),
	msg = ['Unknown neighborhood function: ' content]; 
	sS.neigh = content; isok = 1;
      else
	sS.neigh = content; isok = 1;
      end
     case 'mask', 
      if size(content,1) == 1, content = content'; end
      if ~isnumeric(content) | size(content) ~= [dim 1], 
	msg = '''mask'' should be a column vector (size dim x 1).'; 
      else
	sS.mask = content; isok = 1;
      end
     case 'name', 
      if ~ischar(content), 
	msg = '''name'' should be a string.';
      else 
	sS.name = content; isok = 1;
      end
     case 'comp_names', 
      if ~iscell(content) & ~ischar(content), 
	msg = '''comp_names'' should be a cell string or a string array.'; 
      elseif length(content) ~= dim, 
	msg = 'Length of ''comp_names'' should be equal to dim.'; 
      else
	if ischar(content), content = cellstr(content); end
	if size(content,1)==1, content = content'; end
	sS.comp_names = content;
	isok = 1;
      end        
     case 'comp_norm', 
      if ~iscell(content) & length(content)>0, 
	msg = '''comp_norm'' should be a cell array.'; 
      elseif length(content) ~= dim, 
	msg = 'Length of ''comp_norm'' should be equal to dim.'; 
      else
	isok = 1;
	for j=1:length(content), 
	  if ~isempty(content{j}) & (~isfield(content{j}(1),'type') | ...
				     ~strcmp(content{j}(1).type,'som_norm')), 
	    msg = 'Each cell in ''comp_norm'' should be either empty or type ''som_norm''.';
	    isok = 0; 
	    break; 
	  end
	end
	if isok, sS.comp_norm = content; end
      end        
     case 'trainhist', 
      if ~isstruct(content) & ~isempty(content), 
	msg = '''trainhist'' should be a struct array or empty.';
      else
	isok = 1;
	for j=1:length(content), 
	  if ~isfield(content(j),'type') | ~strcmp(content(j).type,'som_train'), 
	    msg = 'Each cell in ''trainhist'' should be of type ''som_train''.';
	    isok = 0; 
	    break; 
	  end
	end
	if isok, sS.trainhist = content; end      
      end        
     otherwise, 
      msg = ['Invalid field for map struct: ' field]; 
    end
    
   case 'som_data',
    [dlen dim] = size(sS.data);
    switch field,      
     case 'data', 
      [dummy dim2] = size(content);
      if prod(si)==0, 
	msg = '''data'' is empty';
      elseif ~isnumeric(content), 
	msg = '''data'' should be numeric matrix.'; 
      elseif dim ~= dim2 & ~isempty(sS.data), 
	msg = 'New ''data'' must have the same dimension as old one.'; 
      else
	sS.data = content; isok = 1;
      end
     case 'labels', 
      if isempty(content), 
	sS.labels = cell(dlen,1); isok = 1;
      elseif size(content,1) ~= dlen, 
	msg = 'Length of ''labels'' must be equal to the number of data vectors.';
      elseif ~iscell(content) & ~ischar(content), 
	msg = '''labels'' must be a string array or a cell array/matrix.';
      else
	isok = 1;
	if ischar(content), content = cellstr(content); 
	elseif ~iscellstr(content), 
	  l = prod(size(content));
	  for j=1:l, 
	    if ~ischar(content{j}), 
	      if ~isempty(content{j}), 
		msg = 'Invalid ''labels'' array.';
		isok = 0; j
		break;
	      else
		content{j} = ''; 
	      end
	    end
	  end
	end
	if isok, sS.labels = content; end
      end
     case 'name', 
      if ~ischar(content), 
	msg = '''name'' should be a string.';
      else 
	sS.name = content; isok = 1;
      end
     case 'comp_names', 
      if ~iscell(content) & ~ischar(content), 
	msg = '''comp_names'' should be a cell string or a string array.'; 
      elseif length(content) ~= dim, 
	msg = 'Length of ''comp_names'' should be equal to dim.'; 
      else
	if ischar(content), content = cellstr(content); end
	if size(content,1)==1, content = content'; end
	sS.comp_names = content;
	isok = 1;
      end        
     case 'comp_norm', 
      if ~iscell(content) & length(content)>0, 
	msg = '''comp_norm'' should be a cell array.'; 
      elseif length(content) ~= dim, 
	msg = 'Length of ''comp_norm'' should be equal to dim.'; 
      else
	isok = 1;
	for j=1:length(content), 
	  if ~isempty(content{j}) & (~isfield(content{j}(1),'type') | ...
				     ~strcmp(content{j}(1).type,'som_norm')), 
	    msg = 'Each cell in ''comp_norm'' should be either empty or type ''som_norm''.';
	    isok = 0; 
	    break; 
	  end
	end
	if isok, sS.comp_norm = content; end
      end        
     case 'label_names', 
      if ~iscell(content) & ~ischar(content) & ~isempty(content), 
	msg = ['''label_names'' should be a cell string, a string array or' ...
	       ' empty.']; 
      else
	if ~isempty(content), 
	  if ischar(content), content = cellstr(content); end
	  if size(content,1)==1, content = content'; end
	end
	sS.label_names = content;
	isok = 1;
      end        
     otherwise, 
      msg = ['Invalid field for data struct: ' field]; 
    end
    
   case 'som_topol', 
    switch field,      
     case 'msize', 
      if ~isnumeric(content) | ~isvector | ~ispositive | ~isinteger, 
	msg = '''msize'' should be a vector with positive integer elements.'; 
      else
	sS.msize = content; isok=1;
      end
     case 'lattice', 
      if ~ischar(content),
	msg = '''lattice'' should be a string'; 
      elseif ~strcmp(content,'rect') & ~strcmp(content,'hexa'),
	msg = ['Unknown lattice type: ' content]; 
	sS.lattice = content; isok = 1;
      else
	sS.lattice = content; isok = 1;
      end
     case 'shape', 
      if ~ischar(content),
	msg = '''shape'' should be a string';
      elseif ~strcmp(content,'sheet') & ~strcmp(content,'cyl') & ...
	    ~strcmp(content,'toroid'),
	msg = ['Unknown shape type: ' content]; 
	sS.shape = content; isok = 1;
      else
	sS.shape = content; isok = 1;
      end
     otherwise, 
      msg = ['Invalid field for topology struct: ' field]; 
    end
    
   case 'som_train', 
    switch field,      
     case 'algorithm', 
      if ~ischar(content),
	msg = '''algorithm'' should be a string.'; 
      else
	sS.algorithm = content; isok = 1;
      end
     case 'data_name', 
      if ~ischar(content),
	msg = '''data_name'' should be a string'; 
      else
	sS.data_name = content; isok = 1;
      end
     case 'neigh', 
      if ~ischar(content),
	msg = '''neigh'' should be a string'; 
      elseif ~isempty(content) & ~strcmp(content,'gaussian') & ~strcmp(content,'ep') & ...
	    ~strcmp(content,'cutgauss') & ~strcmp(content,'bubble'),
	msg = ['Unknown neighborhood function: ' content]; 
	sS.neigh = content; isok = 1;
      else
	sS.neigh = content; isok = 1;
      end
     case 'mask', 
      if size(content,1) == 1, content = content'; end
      dim = size(content,1); %[munits dim] = size(sS.data); 
      if ~isnumeric(content) | size(content) ~= [dim 1], 
	msg = '''mask'' should be a column vector (size dim x 1).'; 
      else
	sS.mask = content; isok = 1;
      end
     case 'radius_ini', 
      if ~isnumeric(content) | ~isscalar, 
	msg = '''radius_ini'' should be a scalar.'; 
      else
	sS.radius_ini = content; isok = 1;
      end
     case 'radius_fin', 
      if ~isnumeric(content) | ~isscalar, 
	msg = '''radius_fin'' should be a scalar.'; 
      else
	sS.radius_fin = content; isok = 1;
      end
     case 'alpha_ini', 
      if ~isnumeric(content) | ~isscalar,
	msg = '''alpha_ini'' should be a scalar.'; 
      else
	sS.alpha_ini = content; isok = 1;
      end
     case 'alpha_type', 
      if ~ischar(content),
	msg = '''alpha_type'' should be a string'; 
      elseif ~strcmp(content,'linear') & ~strcmp(content,'inv') & ...
	    ~strcmp(content,'power') & ~strcmp(content,'constant') & ~strcmp(content,''),
	msg = ['Unknown alpha type: ' content]; 
	sS.alpha_type = content; isok = 1;
      else
	sS.alpha_type = content; isok = 1;
      end        
     case 'trainlen', 
      if ~isnumeric(content) | ~isscalar, 
	msg = '''trainlen'' should be a scalar.'; 
      else
	sS.trainlen = content; isok = 1;
      end
     case 'time', 
      if ~ischar(content),
	msg = '''time'' should be a string'; 
      else
	sS.time = content; isok = 1;
      end        
     otherwise, 
      msg = ['Invalid field for train struct: ' field]; 
    end
    
   case 'som_norm', 
    switch field,      
     case 'method', 
      if ~ischar(field), 
	msg = '''method'' should be a string.';
      else
	sS.method = content; isok = 1;
      end
     case 'params', 
      sS.params = content; isok = 1;
     case 'status', 
      if ~ischar(content),
	msg = '''status'' should be a string'; 
      elseif ~strcmp(content,'done') & ~strcmp(content,'undone') & ...
	    ~strcmp(content,'uninit'),
	msg = ['Unknown status type: ' content]; 
	sS.status = content; isok = 1;
      else
	sS.status = content; isok = 1;
      end        
     otherwise, 
      msg = ['Invalid field for normalization struct: ' field];
    end

   case 'som_grid', 
    if any(strcmp(field,{'lattice', 'shape', 'msize', 'coord',...
			 'line', 'linecolor', 'linewidth', ...
			 'marker', 'markersize', 'markercolor', 'surf', ... 
			 'label', 'labelcolor', 'labelsize'})), 
      warning('No checking done on field identifier or content.');
      sS = setfield(sS,field,content);       
      isok = 1;
    else
      msg = ['Invalid field for grid struct: ' field];
    end
    
   otherwise, 
    error('Unrecognized structure.'); 
    
  end
  
  msgs{i} = msg;
  ok(i) = isok;
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% return

if nargout < 2, 
  for i=1:p,     
    if ~isempty(msgs{i}), 
      if ~ok(i), fprintf(1,'[Error! '); 
      else fprintf(1,'[Notice ');
      end
      fprintf(1,'in setting %s] ',varargin{2*i-1});
      fprintf(1,'%s\n',msgs{i});
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



