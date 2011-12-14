function state=randomseed(setseed)
% randomseed(seedval)
%
% Sets the use of rand/randn/randi to specific state/seed,
% taking into account the different Matlab version specific methods
%
% INPUT
% setseed:   []              does not reset the state, but saves out the state for future use
%            integer         seed value to set to specifc state
%            state vector    state value (vector) output from previous call to setting the state
%
% OUTPUT
% state      vector of current state (or seed only)
%            can be used as input 'setseed' to reset to same state

% Johanna Zumer (2011)
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.

if isempty(setseed) % save out rand state for later use
  if matlabversion(-Inf,'7.3') 
    rand('twister',sum(100*clock)) % can fail otherwise, if first time rand is called per matlab session
    state=rand('twister');
  elseif matlabversion('7.4','7.6') 
    state=rand('twister');
  elseif matlabversion('7.7','7.11')
    stream = RandStream.getDefaultStream;
    state=stream.State;
  elseif matlabversion('7.12',Inf)
    s=rng;
    state=s.State;
  end
else % seedval is actual random seed value set by user, OR saved state
  if matlabversion(-Inf,'7.6') 
    rand('twister',setseed);
    state=rand('twister');
  elseif matlabversion('7.7','7.11')
    if isscalar(setseed)
      stream=RandStream('mt19937ar','Seed',setseed);
    else
      stream=RandStream.getDefaultStream;
      stream.State=setseed;
    end
    RandStream.setDefaultStream(stream);
    state=stream.State;
  elseif matlabversion('7.12',Inf)
    if isscalar(setseed)
      rng(setseed,'twister')
      s=rng;
    elseif isnumeric(setseed)
      s.Type='twister';
      s.State=setseed;
      s.Seed=0;
      rng(s);
    else
      s.Type='Legacy';
      s.Seed='Not applicable';
      s.State=setseed;
      rng(s);
    end
    state=s.State;
  end
end;


