function this = setContinuity(this, cont_table)
% MULTISCALEELECTROPHYSIOLOGYFILE.SETCONTINUITY Set Continuity table
% 
% Syntax:
%   cont_table = setContinuity(this, cont_table)
% Input(s):
%   this        - [obj] MultiscaleElectrophysiologyFile object
%   cont_table  - [table] Continuity table
% 
% Output(s):
%   this        - [obj] MultiscaleElectrophysiologyFile

% Copyright 2019 Richard J. Cui. Created: Mon 05/20/2019  6:35:28.870 PM
% $Revision: 0.1 $  $Date: Mon 05/20/2019  6:35:28.870 PM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(this, cont_table);

this.Continuity = q.cont_table;

end

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(varargin)

% defaults

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addRequired('cont_table', @istable);


% parse and return the results
p.parse(varargin{:});
q = p.Results;

end % function

% [EOF]

