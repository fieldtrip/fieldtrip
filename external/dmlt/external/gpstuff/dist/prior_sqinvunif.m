function p = prior_sqinvunif(varargin)
%PRIOR_SQINVUNIF  Uniform prior structure for the square inverse of the parameter
%       
%  Description
%    P = PRIOR_SQINVUNIF creates uniform prior structure for the
%    square inverse of the parameter.
%    
%  See also
%    PRIOR_*

% Copyright (c) 2009 Jarno Vanhatalo
% Copyright (c) 2010,2012 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'PRIOR_SQINVUNIFORM';
  ip.addOptional('p', [], @isstruct);
  ip.parse(varargin{:});
  p=ip.Results.p;
  
  if isempty(p)
    init=true;
    p.type = 'SqInv-Uniform';
  else
    if ~isfield(p,'type') && ~isequal(p.type,'SqInv-Uniform')
      error('First argument does not seem to be a valid prior structure')
    end
    init=false;
  end
  
  if init
    % set functions
    p.fh.pak = @prior_sqinvunif_pak;
    p.fh.unpak = @prior_sqinvunif_unpak;
    p.fh.lp = @prior_sqinvunif_lp;
    p.fh.lpg = @prior_sqinvunif_lpg;
    p.fh.recappend = @prior_sqinvunif_recappend;
  end
  
end

function [w, s] = prior_sqinvunif_pak(p, w)
  w=[];
  s={};
end

function [p, w] = prior_sqinvunif_unpak(p, w)
  w = w;
  p = p;
end

function lp = prior_sqinvunif_lp(x, p)
  lJ = -log(x)*3 + log(2);  % log(-2/x^3) log(|J|) of transformation
  lp = sum(lJ);
end

function lpg = prior_sqinvunif_lpg(x, p)
  lJg = -3./x;              % gradient of log(|J|) of transformation
  lpg = lJg;
end

function rec = prior_sqinvunif_recappend(rec, ri, p)
% The parameters are not sampled in any case.
  rec = rec;
end

