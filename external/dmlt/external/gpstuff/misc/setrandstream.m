function prevstream = setrandstream(seed, stream)
%SETRANDSTREAM Set random stream
%
%  Description
%    CURRSTREAM = SETRANDSTREAM()
%    Return current random stream as a random stream object CURRSTREAM.
%
%    PREVSTREAM = SETRANDSTREAM(SEED)
%    Set MATLAB random stream to default stream (Mersenne Twister) with 
%    state/seed SEED and return previous stream PREVSTREAM. This function 
%    takes into account used MATLAB version, and uses correct function 
%    based on that. 
%    
%    PREVSTREAM = SETRANDSTREAM(SEED, STREAM)
%    Set MATLAB random stream to STREAM, with state/seed SEED. Here STREAM
%    is string defining random number generator, e.g. 'mt19937ar' or 
%    'twister' (in new MATLAB versions) for Mersenne Twister.
%
%    PREVSTREAM = SETRANDSTREAM(STREAMOBJ)
%    Set MATLAB random stream to STREAMOBJ. Here STREAMOBJ is random stream
%    object returned from e.g. this function or rng function in new MATLAB
%    versions.
%
%  See also
%    RANDSTREAM, RNG
%   

% Copyright (c) 2012, Ville Tolvanen

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.
  
if nargin>=1 && ~isnumeric(seed) 
  % First argument is random stream object
  stream=seed;
  if str2double(regexprep(version('-release'), '[a-c]', '')) < 2012
    prevstream = RandStream.setDefaultStream(stream);
  else
    prevstream=rng(stream);
  end
else
  if nargin<2
    if nargin<1
      % Get current random stream
      if str2double(regexprep(version('-release'), '[a-c]', '')) < 2012
        prevstream = RandStream.getDefaultStream(stream0);
      else
        prevstream = rng;
      end
      return
    else
      % If stream is not provided, use Mersenne Twister
      stream='mt19937ar';
    end
  end
  if isempty(seed)
    % Default seed
    seed=0;
  end
  if str2double(regexprep(version('-release'), '[a-c]', '')) < 2012
    if ischar(stream)
      stream0 = RandStream(stream,'Seed',seed);
    end
    prevstream = RandStream.setDefaultStream(stream0);
  else
    if ischar(stream)
      prevstream = rng(seed,stream);
    else
      prevstream=rng(stream);
    end
  end
end

end

