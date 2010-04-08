% finputcheck() - check Matlab function {'key','value'} input argument pairs
%
% Usage: >> result = finputcheck( varargin, fieldlist );
%        >> [result varargin] = finputcheck( varargin, fieldlist, ... 
%                                                         callingfunc, mode );
% Input:
%   varargin  - 'varargin' argument from a function call using 'key', 'value'
%               argument pairs.
%   fieldlist - A 3- to 5-column cell array, one row per 'key'. The first
%               column contains the key string, the second its type, 
%               the third the accepted value range, and the fourth the 
%               default value.  Allowed types are 'boolean', 'integer', 
%               'real', 'string', 'cell' or 'struct'.  For example:
%                       {'key1' 'string' { 'string1' 'string2' } 'defaultval_key1'}
%                       {'key2' 'int' { minint maxint } 'defaultval_key2'} 
%  callingfunc - Calling function name for error messages. {default: none}.
%  mode        - ['ignore'|'error'] ignore keywords that are either not specified 
%                in the fieldlist cell array or generate an error. 
%                {default: 'error'}.
% Outputs:
%   result     - If no error, structure with 'key' as fields and 'value' as 
%                content. If error this output contain the string error.
%   varargin   - residual varagin containing unrecognized input arguments.
%                Requires mode 'ignore' above.
%
% Note: In case of error, a string is returned containing the error message
%       instead of a structure.
%
% Example:
% 	g = finputcheck(varargin, ...
%               { 'title'         'string'   []       ''; ...
%                 'percent'       'real'     [0 1]    1 ; ...
%                 'elecamp'       'integer'  [1:10]   [] });
% Note: 
%   The 'title' argument should be a string. {no default value}
%   The 'percent' argument should be a real number between 0 and 1. {default: 1}
%   The 'elecamp' argument should be an integer between 1 and 10 (inclusive).
%
%   Now 'g.title' will contain the title arg (if any, else the default ''), etc.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 10 July 2002

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 10 July 2002, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: finputcheck.m,v $
% Revision 1.23  2004/11/06 02:54:06  scott
% a few further small edits to the help msg -sm
%
% Revision 1.22  2004/11/05 15:23:37  arno
% ,sg
%
% Revision 1.21  2004/11/05 04:10:44  scott
% help msg. -sm
%
% Revision 1.20  2004/06/09 16:30:42  arno
% adding or if several types
%
% Revision 1.19  2003/10/29 16:35:57  arno
% msg typo
%
% Revision 1.18  2003/07/30 23:53:58  arno
% debug multiple return values
%
% Revision 1.17  2003/07/26 00:21:17  arno
% allowing cell array for values
%
% Revision 1.16  2003/06/30 02:10:10  arno
% strmatch exact
%
% Revision 1.15  2003/01/31 02:35:38  arno
% debugging lowercase/upercase problem
%
% Revision 1.14  2002/11/20 01:05:44  arno
% take into account duplicate parameters
%
% Revision 1.13  2002/11/18 17:15:18  arno
% adding float arg (=real)
%
% Revision 1.12  2002/11/15 02:16:50  arno
% header for web
%
% Revision 1.11  2002/09/30 15:29:23  arno
% autorizing cell arrays for types
%
% Revision 1.10  2002/09/30 00:42:08  arno
% debug input arguments
%
% Revision 1.9  2002/07/29 18:00:53  arno
% debugging for NaN
%
% Revision 1.8  2002/07/29 17:24:22  arno
% header
%
% Revision 1.7  2002/07/20 19:10:41  arno
% debugging output
%
% Revision 1.6  2002/07/19 17:58:11  arno
% returning non-matched 'key' 'val' arguments
%
% Revision 1.5  2002/07/19 17:46:53  arno
% g empty if no varargin
%
% Revision 1.4  2002/07/19 16:27:14  arno
% adding ignore mode
%
% Revision 1.3  2002/07/10 02:18:32  arno
% header info
%
% Revision 1.2  2002/07/10 02:17:27  arno
% debugging error message passing
%
% Revision 1.1  2002/07/10 01:03:19  arno
% Initial revision
%

function [g, varargnew] = finputcheck( vararg, fieldlist, callfunc, mode )

	if nargin < 2
		help finputcheck;
		return;
	end;
	if nargin < 3
		callfunc = '';
	else 
		callfunc = [callfunc ' ' ];
	end;
    if nargin < 4
        mode = 'do not ignore';
    end;
	NAME = 1;
	TYPE = 2;
	VALS = 3;
	DEF  = 4;
	SIZE = 5;
	
	varargnew = {};
	% create structure
	% ----------------
	if ~isempty(vararg)
        vararg = removedup(vararg);
		for index=1:length(vararg)
			if iscell(vararg{index})
				vararg{index} = {vararg{index}};
			end;
		end;
		try
			g = struct(vararg{:});
		catch
			g = [ callfunc 'error: bad ''key'', ''val'' sequence' ]; return;
		end;
	else 
		g = [];
	end;
	
	for index = 1:size(fieldlist,NAME)
		% check if present
		% ----------------
		if ~isfield(g, fieldlist{index, NAME})
			g = setfield( g, fieldlist{index, NAME}, fieldlist{index, DEF});
		end;
		tmpval = getfield( g, {1}, fieldlist{index, NAME});
		
		% check type
		% ----------
        if ~iscell( fieldlist{index, TYPE} )
            res = fieldtest( fieldlist{index, NAME},  fieldlist{index, TYPE}, ...
                           fieldlist{index, VALS}, tmpval, callfunc );
            if isstr(res), g = res; return; end;
        else 
            testres = 0;
            tmplist = fieldlist;
            for it = 1:length( fieldlist{index, TYPE} )
                if ~iscell(fieldlist{index, VALS})
                     res{it} = fieldtest(  fieldlist{index, NAME},  fieldlist{index, TYPE}{it}, ...
                                           fieldlist{index, VALS}, tmpval, callfunc );
                else res{it} = fieldtest(  fieldlist{index, NAME},  fieldlist{index, TYPE}{it}, ...
                                           fieldlist{index, VALS}{it}, tmpval, callfunc );
                end;
                if ~isstr(res{it}), testres = 1; end;
            end;
            if testres == 0,
                g = res{1};
                for tmpi = 2:length(res)
                    g = [ g 10 'or ' res{tmpi} ];
                end;
                return; 
            end;
        end;
	end;
    
    % check if fields are defined
	% ---------------------------
	allfields = fieldnames(g);
	for index=1:length(allfields)
		if isempty(strmatch(allfields{index}, fieldlist(:, 1)', 'exact'))
			if ~strcmpi(mode, 'ignore')
				g = [ callfunc 'error: undefined argument ''' allfields{index} '''']; return;
			end;
			varargnew{end+1} = allfields{index};
			varargnew{end+1} = getfield(g, {1}, allfields{index});
		end;
	end;


function g = fieldtest( fieldname, fieldtype, fieldval, tmpval, callfunc );
	NAME = 1;
	TYPE = 2;
	VALS = 3;
	DEF  = 4;
	SIZE = 5;
    g = [];
    
    switch fieldtype
     case { 'integer' 'real' 'boolean' 'float' }, 
      if ~isnumeric(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be numeric' ]; return;
      end;
      if strcmpi(fieldtype, 'boolean')
          if tmpval ~=0 & tmpval ~= 1
              g = [ callfunc 'error: argument ''' fieldname ''' must be 0 or 1' ]; return;
          end;  
      else 
          if strcmpi(fieldtype, 'integer')
              if ~isempty(fieldval)
                  if (isnan(tmpval) & ~any(isnan(fieldval))) ...
                          & (~ismember(tmpval, fieldval))
                      g = [ callfunc 'error: wrong value for argument ''' fieldname '''' ]; return;
                  end;
              end;
          else % real or float
              if ~isempty(fieldval)
                  if tmpval < fieldval(1) | tmpval > fieldval(2)
                      g = [ callfunc 'error: value out of range for argument ''' fieldname '''' ]; return;
                  end;
              end;
          end;
      end;  
      
      
     case 'string'
      if ~isstr(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be a string' ]; return;
      end;
      if ~isempty(fieldval)
          if isempty(strmatch(lower(tmpval), lower(fieldval), 'exact'))
              g = [ callfunc 'error: wrong value for argument ''' fieldname '''' ]; return;
          end;
      end;

      
     case 'cell'
      if ~iscell(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be a cell array' ]; return;
      end;
      
      
     case 'struct'
      if ~isstruct(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be a structure' ]; return;
      end;
      
      
     case '';
     otherwise, error([ 'finputcheck error: unrecognized type ''' fieldname '''' ]);
    end;

% remove duplicates in the list of parameters
% -------------------------------------------
function cella = removedup(cella)
    [tmp indices] = unique(cella(1:2:end));
    if length(tmp) ~= length(cella)/2
        fprintf('Warning: duplicate ''key'', ''val'' parameter(s), keeping the last one(s)\n');
    end;
    cella = cella(sort(union(indices*2-1, indices*2)));
