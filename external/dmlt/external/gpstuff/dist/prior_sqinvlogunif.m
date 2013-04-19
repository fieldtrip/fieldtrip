function p = prior_sqinvlogunif(varargin)
%PRIOR_SQINVLOGUNIF  Uniform prior structure for the log of the square inverse of parameter
%       
%  Description
%    P = PRIOR_SQINVLOGUNIF creates uniform prior structure for the
%    log of the square inverse of the parameter.
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
  ip.FunctionName = 'PRIOR_SQINVLOGUNIFORM';
  ip.addOptional('p', [], @(x) isstruct(x) || isempty(x));
  ip.parse(varargin{:});
  p=ip.Results.p;
  
  if isempty(p)
    init=true;
    p.type = 'SqInv-Log-Uniform';
  else
    if ~isfield(p,'type') && ~isequal(p.type,'SqInv-Log-Uniform')
      error('First argument does not seem to be a valid prior structure')
    end
    init=false;
  end
  
  if init
    % set functions
    p.fh.pak = @prior_sqinvlogunif_pak;
    p.fh.unpak = @prior_sqinvlogunif_unpak;
    p.fh.lp = @prior_sqinvlogunif_lp;
    p.fh.lpg = @prior_sqinvlogunif_lpg;
    p.fh.recappend = @prior_sqinvlogunif_recappend;
  end
end

function [w, s] = prior_sqinvlogunif_pak(p)
  w=[];
  s={};
end

function [p, w] = prior_sqinvlogunif_unpak(p, w)
  w = w;
  p = p;
end

function lp = prior_sqinvlogunif_lp(x, p)
  lJ = log(2./x);   % log(|J|) of transformation
  lp = sum(lJ);
end

function lpg = prior_sqinvlogunif_lpg(x, p)
  lJg = -1./x;      % gradient of log(|J|) of transformation
  lpg = lJg;
end

function rec = prior_sqinvlogunif_recappend(rec, ri, p)
% The parameters are not sampled in any case.
  rec = rec;
end
