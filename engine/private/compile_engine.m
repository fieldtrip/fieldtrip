% compile the engine mex file
% note that I am not sure about the pthread libs on 32 and 64 bit windows

switch mexext
  case {'mexmaci64' 'mexmaci' 'mexa64' 'mexglx'}
    mex -DUSE_PTHREADS alternative.c -leng -output engine
    mex ft_getopt.c
    
  case 'mexw64'
    % mex -I.\pthreads-win32\include -L.\pthreads-win32\lib\x64 engine.c -leng -lpthreadVC2
    mex alternative.c -leng -output engine
    
  case 'mexw32'
    % mex -I.\pthreads-win32\include -L.\pthreads-win32\lib\x86 engine.c -leng -lpthreadVC2
    mex alternative.c -leng -output engine
    
  otherwise
    error('unsupported mex platform');
end
