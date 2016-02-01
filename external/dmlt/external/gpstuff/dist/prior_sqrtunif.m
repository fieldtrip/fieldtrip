function p = prior_sqrtunif(varargin)
%PRIOR_SQRTUNIF  Uniform prior structure for the square root of the parameter
%       
%  Description
%    P = PRIOR_SQRTUNIF creates uniform prior structure for the
%    square root of the parameter.
%    
%  See also
%    PRIOR_*
  
% Copyright (c) 2009 Jarno Vanhatalo
% Copyright (c) 2010 Jaakko Riihimäki
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'PRIOR_SQRTUNIFORM';
  ip.addOptional('p', [], @isstruct);
  ip.parse(varargin{:});
  p=ip.Results.p;
  
  if isempty(p)
    init=true;
    p.type = 'Sqrt-Uniform';
  else
    if ~isfield(p,'type') && ~isequal(p.type,'Sqrt-Uniform')
      error('First argument does not seem to be a valid prior structure')
    end
    init=false;
  end
  
  if init
    % set functions
    p.fh.pak = @prior_sqrtunif_pak;
    p.fh.unpak = @prior_sqrtunif_unpak;
    p.fh.lp = @prior_sqrtunif_lp;
    p.fh.lpg = @prior_sqrtunif_lpg;
    p.fh.recappend = @prior_sqrtunif_recappend;
  end

end

function [w, s] = prior_sqrtunif_pak(p)
  w=[];
  s={};
end

function [p, w] = prior_sqrtunif_unpak(p, w)
  w = w;
  p = p;
end

function lp = prior_sqrtunif_lp(x, p)
  lJ  = -log(2*sqrt(x));  % log(1/(2*sqrt(x))) log(|J|) of transformation
  lp  = sum(lJ);
end

function lpg = prior_sqrtunif_lpg(x, p)
  lJg = -1./(2*x);        % gradient of log(|J|) of transformation
  lpg = lJg;
end

function rec = prior_sqrtunif_recappend(rec, ri, p)
% The parameters are not sampled in any case.
  rec = rec;
end
