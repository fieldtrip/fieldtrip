function flag=vis_valuetype(value, valid, str);

% VIS_VALUETYPE Used for type checks in SOM Toolbox visualization routines
%
%  flag = vis_valuetype(value, valid, str)
%
%  Input and output arguments:
%   value  (varies) variable to be checked
%   valid  (cell array) size 1xN, cells are strings or vectors (see below)
%   str    (string) 'all' or 'any' (default), determines whether
%                   all or just any of the types listed in argument 'valid' 
%                   should be true for 'value'
%
%   flag   (scalar) 1 or 0 (true or false) 
%
% This is an internal function of SOM Toolbox visualization.  It makes
% various type checks. For example:
%
%  % Return 1 if X is a numeric scalar otherwise 0:
%  f=vis_valuetype(X,{'1x1'});
%
%  % Return 1 if X is a ColorSpec, that is, a 1x3 vector presenting an RGB
%  % value or any of strings 'red','blue','green','yellow','magenta','cyan',
%  % 'white' or 'black' or their shortenings  'r','g','b','y','m','c','w','k': 
%  f=vis_valueype(X,{'1x3rgb','colorstyle'})
%
%  % Return 1 if X is _both_ 10x3 size numeric matrix and has RGB values as rows
%  f=vis_valuetype(X,{'nx3rgb',[10 3]},'all')
%
% Strings that may be used in argument valid: 
%  id             is true if value is 
% 
%  [n1 n2 ... nn] any n1 x n2 x ... x nn sized numeric matrix
%  '1x1'          scalar (numeric)
%  '1x2'          1x2 vector (numeric)
%  'nx1'          any nx1 numeric vector
%  'nx2'              nx2
%  'nx3'              nx3
%  'nxn'          any numeric square matrix
%  'nxn[0,1]'     numeric square matrix with values in interval [0,1]
%  'nxm'          any numeric matrix
%  '1xn'          any 1xn numeric vector
%  '1x3rgb'       1x3 vector v for which all(v>=0 & v<=1), e.g., a RGB code
%  'nx3rgb'       nx3 numeric matrix that contains n RGB values as rows
%  'nx3dimrgb'    nx3xdim numeric matrix that contains RGB values
%  'nxnx3rgb'     nxnx3 numeric matrix of nxn RGB triples
%  'none'         string 'none'
%  'xor'          string 'xor'
%  'indexed'      string 'indexed'
%  'colorstyle'   strings 'red','blue','green','yellow','magenta','cyan','white' 
%                 or 'black', or 'r','g','b','y','m','c','w','k'                 
%  'markerstyle'  any of Matlab's marker chars '.','o','x','+','*','s','d','v',
%                 '^','<','>','p'or 'h'
%  'linestyle'    any or Matlab's line style strings '-',':','--', or '-.'
%  'cellcolumn'   a nx1 cell array
%  'topol_cell'   {lattice, msize, shape} 
%  'topol_cell_no_shape' {lattice, msize}
%  'string'       any string (1xn array of char)  
%  'chararray'    any MxN char array

% Copyright (c) 1999-2000 by the SOM toolbox programming team.
% http://www.cis.hut.fi/projects/somtoolbox/             

% Version 2.0beta Johan 201099 juuso 280800

if nargin == 2
  str='any';
end

flag=0;  
sz=size(value);
dims=ndims(value);

% isnumeric
numeric=isnumeric(value);
character=ischar(value);

% main loop: go through all types in arg. 'valid'
for i=1:length(valid),
  if isnumeric(valid{i}), % numeric size for double matrix
    if numeric & length(valid{i}) == dims,
      flag(i)=all(sz == valid{i});
    else
      flag(i)=0; % not numeric or wrong dimension
    end
  else
    msg=''; % for a error message inside try 
    try 
      switch valid{i}
	
	% scalar
       case '1x1'
	flag(i)=numeric & dims == 2 & sz(1)==1 & sz(2) ==1;
	
	% 1x2 numeric vector
       case '1x2'
	flag(i)=numeric & dims == 2 & sz(1)==1 & sz(2) == 2;
	
	% 1xn numeric vector
       case '1xn'
	flag(i)=numeric & dims == 2 & sz(1) == 1;
	
	% any numeric matrix
       case 'nxm' 
	flag(i)=numeric & dims == 2;
	
	% nx3 numeric matrix 
       case 'nx3'
	flag(i)=numeric & dims == 2 & sz(2) == 3;
	
	% nx2 numeric matrix 
       case 'nx2'
	flag(i)=numeric & dims == 2 & sz(2) == 2;
	
	% nx1 numeric vector
       case 'nx1'
	flag(i)=numeric & dims == 2 & sz(2) == 1;
       
	% nx1xm numric matrix
       case 'nx1xm'
	flag(i)=numeric & dims == 3 & sz(2) == 1;
	
	% nx3 matrix of RGB triples
       case 'nx3rgb'  
	flag(i)=numeric & dims == 2 & sz(2) == 3 & in0_1(value);
	
	% RGB triple (ColorSpec vector)
       case '1x3rgb'
	flag(i) = numeric & dims == 2 & sz(1)==1 & sz(2) == 3 & in0_1(value);
	
	% any square matrix
       case 'nxn'
	flag(i)=numeric & dims == 2 & sz(1) == sz(2);
	
	% nx3xdim array of nxdim RGB triples
       case 'nx3xdimrgb'
	flag(i)=numeric & dims == 3 & sz(2) == 3 & in0_1(value);
	
	% nxnx3 array of nxn RGB triples
       case 'nxnx3rgb'
	flag(i)= numeric & dims == 3 & sz(1) == sz(2) & sz(3) == 3 ...
		 & in0_1(value);
	
	% nxn matrix of values between [0,1]
       case 'nxn[0,1]' 
	
	flag(i)=numeric & dims == 2 & sz(1) == sz(2) & in0_1(value);
	
	% string 'indexed'
       case 'indexed'
	flag(i) = ischar(value) & strcmp(value,'indexed');
	
	% string 'none'
       case 'none'
	flag(i) = character & strcmp(value,'none');
      
	% string 'xor'
       case 'xor'
	flag(i) = character & strcmp(value,'xor');
	
	% any string (1xn char array)
       case 'string'
	flag(i) = character & dims == 2 & sz(1)<=1;
	
	% any char array
       case 'chararray'
	flag(i) = character & dims == 2 & sz(1)>0;
	
	% ColorSpec string
       case 'colorstyle'
	flag(i)=(character &  sz(1) == 1 & sz(2) == 1 & ...
		 any(ismember('ymcrgbwk',value))) | ...
	(ischar(value) & any(strcmp(value,{'none','yellow','magenta',...
		    'cyan','red','green','blue','white','black'})));
	
	% any valid Matlab's Marker
       case 'markerstyle'
	flag(i)=character &  sz(1) == 1 & sz(2) == 1 & ...
		any(ismember('.ox+*sdv^<>ph',value));
	
	% any valid Matlab's LineStyle
       case 'linestyle'
	str=strrep(strrep(strrep(value,'z','1'),'--','z'),'-.','z');
	flag(i)=character & any(ismember(str,'z-:')) & sz(1)==1 & (sz(2)==1 | sz(2)==2);
	
	% any struct
       case 'struct'
	flag(i)=isstruct(value);
	
	% nx1 cell array of strings
       case 'cellcolumn_of_char'
	flag(i)=iscell(value) & dims == 2 & sz(2)==1;  
	try, char(value); catch, flag(i)=0; end
	
	% mxn cell array of strings
       case '2Dcellarray_of_char'  
	flag(i)=iscell(value) & dims == 2; 
	try, char(cat(2,value{:})); catch, flag(i)=0; end
	
	% valid {lattice, msize} 
       case 'topol_cell_no_shape'
	flag(i)=1;
	if ~iscell(value) | length(size(value)) ~= 2 | size(value,2)~=2
	  flag(i)=0;
	else
	  if vis_valuetype(value{1},{'string'}),
	    switch value{1}
	     case { 'hexa','rect'}
	      ;
	     otherwise
	      flag(i)=0;
	    end
	  end
	  if ~vis_valuetype(value{2},{'1xn'}),
	    flag(i)=0;
	  end
	end
	
	% valid {lattice, msize, shape} 
       case 'topol_cell'
	flag(i)=1;
	if ~iscell(value) | length(size(value)) ~= 2 | size(value,2) ~= 3,
	  flag(i)=0;
	else
	  if vis_valuetype(value{1},{'string'}),
	    switch value{1}
	     case { 'hexa','rect'}
	      ;
	     otherwise
	      flag(i)=0;
	    end
	  end
	  if ~vis_valuetype(value{2},{'1xn'})
	    flag(i)=0;
	  end
	  if ~vis_valuetype(value{3},{'string'})
	    flag(i)=0;
	  else
	    switch value{3}
	     case { 'sheet','cyl', 'toroid'}
	      ;
	     otherwise
	      flag(i)=0;
	    end
	  end
	end
       otherwise
	msg='Unknown valuetype!';
      end
    catch 
      % error during type check is due to wrong type of value: 
      % lets set flag(i) to 0
      flag(i)=0; 
    end
    % Unknown indetifier?
    error(msg);
  end
  % set flag according to 3rd parameter (all ~ AND, any ~ OR) 
  if strcmp(str,'all');
    flag=all(flag);
  else
    flag=any(flag);
  end
end


function f=in0_1(value)

f=all(value(:) >= 0 & value(:)<=1);