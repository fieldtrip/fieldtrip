function p = prior_sqrtinvunif(varargin)
%PRIOR_INVSQRTUNIF  Uniform prior structure for the square root of inverse of the parameter
%       
%  Description
%    P = PRIOR_INVSQRTUNIF creates uniform prior structure for the
%    square root of inverse of the parameter.
%    
%  See also
%    PRIOR_*
  
% Copyright (c) 2009 Jarno Vanhatalo
% Copyright (c) 2010 Jaakko Riihimäki
% Copyright (c) 2010,2012 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'PRIOR_INVSQRTUNIFORM';
  ip.addOptional('p', [], @isstruct);
  ip.parse(varargin{:});
  p=ip.Results.p;
  
  if isempty(p)
    init=true;
    p.type = 'SqrtInv-Uniform';
  else
    if ~isfield(p,'type') && ~isequal(p.type,'SqrtInv-Uniform')
      error('First argument does not seem to be a valid prior structure')
    end
    init=false;
  end
  
  if init
    % set functions
    p.fh.pak = @prior_invsqrtunif_pak;
    p.fh.unpak = @prior_invsqrtunif_unpak;
    p.fh.lp = @prior_invsqrtunif_lp;
    p.fh.lpg = @prior_invsqrtunif_lpg;
    p.fh.recappend = @prior_invsqrtunif_recappend;
  end

end

function [w, s] = prior_invsqrtunif_pak(p)
  w=[];
  s={};
end

function [p, w] = prior_invsqrtunif_unpak(p, w)
  w = w;
  p = p;
end

function lp = prior_invsqrtunif_lp(x, p)
  lJ = -log(2*x.^(3/2));  % log(1/(2*x^(3/2))) log(|J|) of transformation
  lp = -sum(lJ);
end

function lpg = prior_invsqrtunif_lpg(x, p)
  lJg = -3./(2*x);        % gradient of log(|J|) of transformation
  lpg = lJg;
end

function rec = prior_invsqrtunif_recappend(rec, ri, p)
% The parameters are not sampled in any case.
  rec = rec;
end
