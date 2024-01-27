function this = setSessionPassword(this, password)
% SETSESSIONPASSWORD set Session Password of MEF file
% 
% Syntax:
%   this = setSessionPassword(this, password)
% 
% Input(s):
%   this        - [obj] MultiscaleElectrophysiologyFile object
%   password    - [str] Session password
% 
% Output(s):
%   this        - output MultiscaleElectrophysiologyFile

% Copyright 2019 Richard J. Cui. Created: Mon 05/20/2019 10:02:51.052 PM
% $Revision: 0.1 $  $Date: Mon 05/20/2019 10:02:51.052 PM $
%
% 1026 Rocky Creek Dr NE
% Rochester, MN 55906, USA
%
% Email: richard.cui@utoronto.ca

% =========================================================================
% parse inputs
% =========================================================================
q = parseInputs(this, password);

this.SessionPassword = q.password;

end

% =========================================================================
% subroutines
% =========================================================================
function q = parseInputs(varargin)

% defaults

% parse rules
p = inputParser;
p.addRequired('this', @isobject);
p.addRequired('password', @ischar);

% parse and return the results
p.parse(varargin{:});
q.password = p.Results.password;

end % function

% [EOF]
