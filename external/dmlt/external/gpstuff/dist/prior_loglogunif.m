function p = prior_loglogunif(varargin)
%PRIOR_LOGLOGUNIF  Uniform prior structure for the log-log of the parameter
%       
%  Description
%    P = PRIOR_LOGLOGUNIF creates uniform prior structure for the
%    log-log of the parameters.
%    
%  See also
%    PRIOR_*

% Copyright (c) 2009 Jarno Vanhatalo
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'PRIOR_LOGLOGUNIFORM';
  ip.addOptional('p', [], @isstruct);
  ip.parse(varargin{:});
  p=ip.Results.p;
  
  if isempty(p)
    init=true;
    p.type = 'Loglog-Uniform';
  else
    if ~isfield(p,'type') && ~isequal(p.type,'Loglog-Uniform')
      error('First argument does not seem to be a valid prior structure')
    end
    init=false;
  end
  
  if init
    % set functions
    p.fh.pak = @prior_loglogunif_pak;
    p.fh.unpak = @prior_loglogunif_unpak;
    p.fh.lp = @prior_loglogunif_lp;
    p.fh.lpg = @prior_loglogunif_lpg;
    p.fh.recappend = @prior_loglogunif_recappend;
  end

end

function [w, s] = prior_loglogunif_pak(p, w)
  w=[];
  s={};
end

function [p, w] = prior_loglogunif_unpak(p, w)
  w = w;
  p = p;
end

function lp = prior_loglogunif_lp(x, p)
  lp = -sum(log(log(x)) + log(x));     % = log( 1./log(x) * 1./x)
end

function lpg = prior_loglogunif_lpg(x, p)
  lpg = -1./log(x)./x - 1./x;
end

function rec = prior_loglogunif_recappend(rec, ri, p)
% The parameters are not sampled in any case.
  rec = rec;
end
