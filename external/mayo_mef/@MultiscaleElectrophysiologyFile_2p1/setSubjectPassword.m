function this = setSubjectPassword(this, password)
% SETSUBJECTPASSWORD set Subject Password of MEF file
% 
% Syntax:
%   this = setSubjectPassword(this, password)
% 
% Input(s):
%   this        - [obj] MultiscaleElectrophysiologyFile object
%   password    - [str] Subject password
% 
% Output(s):
%   this        - output MultiscaleElectrophysiologyFile

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
q = parseInputs(this, password);

this.SubjectPassword = q.password;

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
