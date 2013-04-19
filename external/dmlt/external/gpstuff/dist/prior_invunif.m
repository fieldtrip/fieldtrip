function p = prior_invunif(varargin)
%PRIOR_INVUNIF  Uniform prior structure for the inverse of the parameter
%       
%  Description
%    P = PRIOR_INVUNIF creates uniform prior structure for the
%    inverse of the parameter.
%    
%  See also
%    PRIOR_*

% Copyright (c) 2009 Jarno Vanhatalo
% Copyright (c) 2010,2012 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'PRIOR_INVUNIFORM';
  ip.addOptional('p', [], @isstruct);
  ip.parse(varargin{:});
  p=ip.Results.p;
  
  if isempty(p)
    init=true;
    p.type = 'Inv-Uniform';
  else
    if ~isfield(p,'type') && ~isequal(p.type,'Inv-Uniform')
      error('First argument does not seem to be a valid prior structure')
    end
    init=false;
  end
  
  if init
    % set functions
    p.fh.pak = @prior_invunif_pak;
    p.fh.unpak = @prior_invunif_unpak;
    p.fh.lp = @prior_invunif_lp;
    p.fh.lpg = @prior_invunif_lpg;
    p.fh.recappend = @prior_invunif_recappend;
  end
  
end

function [w, s] = prior_invunif_pak(p, w)
  w=[];
  s={};
end

function [p, w] = prior_invunif_unpak(p, w)
  w = w;
  p = p;
end

function lp = prior_invunif_lp(x, p)
  lJ=-log(x)*2;   % log(1/x^2) log(|J|) of transformation
  lp = sum(0 +lJ);
end

function lpg = prior_invunif_lpg(x, p)
  lJg=-2./x;      % gradient of log(|J|) of transformation
  lpg = zeros(size(x)) + lJg;
end

function rec = prior_invunif_recappend(rec, ri, p)
% The parameters are not sampled in any case.
  rec = rec;
end

