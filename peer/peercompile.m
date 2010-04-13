function peercompile(action)

% PEERCOMPILE compile the low-level mex file implements that implements
% the TCP/IP, announce and discover services.

% -----------------------------------------------------------------------
% Copyright (C) 2010, Robert Oostenveld
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/
% -----------------------------------------------------------------------

if ~nargin, action = 'all'; end

cwd = pwd;
src = fullfile(fileparts(which(mfilename)),'src');

switch action
    
    case 'all'
        try
            cd(src);
            
            if ispc
                includes     = '-I..\pthreads-win32\include';
                libdirs      = '-L..\pthreads-win32\lib';
                switch upper(computer)
                    case 'PCWIN32'
                        libs = '-lpthreadVC2';
                    case 'PCWIN64'
                        libs = '-lpthreadGC2';
                end
                objext       = '.obj';
            else
                includes     = '-I.';
                libdirs      = '-L.';
                libs         = '-lpthread';
                objext       = '.o';
            end
            
            % create object files
            cfiles = {'announce.c' 'discover.c' 'expire.c' 'extern.c' ...
                'fairshare.c' 'peerinit.c' 'util.c' 'tcpserver.c' 'tcpsocket.c' ...
                'security.c'};
            ofiles = cell(1,numel(cfiles));
            for i=1:numel(cfiles)
                mex(includes,'-c',cfiles{i});
                [p, n]    = fileparts(cfiles{i});
                ofiles{i} = [n objext];
            end
            
            % compile the peer gateway
            mex(includes, libdirs, '-outdir', '..', 'peer.c', ofiles{:}, libs);
            
            % compile the memory profiler
            mex(includes, libdirs, '-outdir', '..', 'memprofile.c', libs);
            
            % compile the delayed exit
            mex(includes, libdirs, '-outdir', '..', 'delayedexit.c', libs);
            
            cd(cwd);
        catch
            cd(cwd);
            rethrow(lasterror);
        end
        
    case 'clean'
        try
            cd(src);
            delete('*.o*');
            cd(cwd);
        catch
            cd(cwd);
            rethrow(lasterror);
        end
        
    case 'distclean'
        peercompile('clean');
        try
            cd(fullfile(src,'..'));
            delete(['*.' mexext]);
            cd(cwd);
        catch
            cd(cwd);
            rethrow(lasterror);
        end
        
    otherwise
        error('Unknown action.');
end
