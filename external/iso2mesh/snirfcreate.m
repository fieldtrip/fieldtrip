function snf = snirfcreate(varargin)
%
%    snf=snirfcreate
%       or
%    snf=snirfcreate(option)
%    snf=snirfcreate('Format',format,'Param1',value1, 'Param2',value2,...)
%
%    Create a empty SNIRF data structure defined in the SNIRF
%    specification: https://github.com/fNIRS/snirf
%
%    author: Qianqian Fang (q.fang <at> neu.edu)
%
%    input:
%        option (optional): option can be ignored. If it is a string with a
%             value 'snirf', this creates a default SNIRF data structure;
%             otherwise, a JSNIRF data structure is created.
%        format: save as option.
%        param/value:   a list of name/value pairs specify
%             additional subfields to be stored under the /nirs object.
%
%    output:
%        snf: a default SNIRF or JSNIRF data structure.
%
%    example:
%        snf=snirfcreate('data',mydata,'aux',myauxdata,'comment','test');
%
%    this file is part of JSNIRF specification: https://github.com/fangq/snirf
%
%    License: GPLv3 or Apache 2.0, see https://github.com/NeuroJSON/jsnirf for details
%

if (nargin == 1)
    snf = jsnirfcreate(varargin{:});
else
    snf = jsnirfcreate('Format', 'snirf', varargin{:});
end
