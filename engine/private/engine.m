% ENGINE is a multi-threaded mex file to manage multiple MATLAB engines
% that can be used for distributed computing. This mex file provides the
% low-level interface.
%
% Use as
%   status = engine('open', poolsize, startcmd)
%   status = engine('close')
%   status = engine('isbusy', num)
%   status = engine('eval',   num, cmd)
%   status = engine('put',    num, name, value)
%   value  = engine('get',    num, name)
%   value  = engine('poolsize')
%            engine('info')
%
% See also ENGPOOL, ENGFEVAL, ENGCELLFUN

error('The mex file is not compiled for your platform (%s)', mexext);

