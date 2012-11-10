% compile the engine mex file

% I am not sure about the pthread libs on 32 and 64 bit windows

switch mexext
  case {'mexmaci64' 'mexmaci' 'mexa64' 'mexglx'}
    mex engine.c -leng
    mex ft_getopt.c
    
  case 'mexw64'
    mex -I.\pthreads-win64\include -L.\pthreads-win64\lib engine.c -leng -lpthreadGC
    mex ft_getopt.c
    
  case 'mexw32'
    mex -I.\pthreads-win32\include -L.\pthreads-win32\lib engine.c -leng -lpthreadVC2
    mex ft_getopt.c
    
  otherwise
    error('unsupported mex platform');
end
