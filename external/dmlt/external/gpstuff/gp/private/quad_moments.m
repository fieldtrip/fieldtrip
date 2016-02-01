function [m_0, m_1, m_2] = quad_moments(fun, a, b, rtol, atol, minsubs)
% QUAD_MOMENTS Calculate the 0th, 1st and 2nd moment of a given
%              (unnormalized) probability distribution
%
%   [m_0, m_1, m_2] = quad_moments(fun, a, b, varargin) 
%   Inputs:
%      fun  = Function handle to the unnormalized probability distribution
%      a,b  = integration limits [a,b]
%      rtol = relative tolerance for the integration (optional, default 1e-6)
%      atol = absolute tolerance for the integration (optional, default 1e-10)
%               
%   Returns the first three moments:
%      m0  = int_a^b fun(x) dx
%      m1  = int_a^b x*fun(x) dx / m0
%      m2  = int_a^b x^2*fun(x) dx / m0
%
%   The function uses an adaptive Gauss-Kronrod quadrature. The same set of 
%   integration points and intervals are used for each moment. This speeds up 
%   the evaluations by factor 3, since the function evaluations are done only 
%   once.
% 
%   The quadrature method is described by:
%   L.F. Shampine, "Vectorized Adaptive Quadrature in Matlab",
%   Journal of Computational and Applied Mathematics, 211, 2008, 
%   pp. 131-140.

%   Copyright (c) 2010 Jarno Vanhatalo, Jouni Hartikainen
    
% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

    maxsubs = 650;
    
    if nargin < 4
        rtol = 1.e-6;
    end
    if nargin < 5
        atol = 1.e-10;
    end
    if nargin < 6
        minsubs = 10;
    end
    
    rtol = max(rtol,100*eps);
    atol = max(atol,0);
    minsubs = max(minsubs,2); % At least two subintervals are needed
    
    % points and weights
    points15 = [0.2077849550078985; 0.4058451513773972; 0.5860872354676911; ...
        0.7415311855993944; 0.8648644233597691; 0.9491079123427585; ...
        0.9914553711208126];
    points = [-points15(end:-1:1); 0; points15];
    
    w15 = [0.2044329400752989, 0.1903505780647854, 0.1690047266392679, ...
        0.1406532597155259, 0.1047900103222502, 0.06309209262997855, ...
        0.02293532201052922];
    w = [w15(end:-1:1), 0.2094821410847278, w15];
    
    w7 = [0,0.3818300505051189,0,0.2797053914892767,0,0.1294849661688697,0];
    ew = w - [w7(end:-1:1), 0.4179591836734694, w7];
        
    samples = numel(w);
    
    % split the interval.
    if b-a <= 0
        c = a; a = b; b=c;
        warning('The start of the integration interval was less than the end of it.')
    end
    apu = a + (1:(minsubs-1))./minsubs*(b-a);
    apu = [a,apu,b];
    subs = [apu(1:end-1);apu(2:end)];
        
    % Initialize partial sums.
    Ifx_ok = 0;
    Ifx1_ok = 0;
    Ifx2_ok = 0;
    % The main loop
    while true
        % subintervals and their midpoints
        midpoints = sum(subs)/2;   
        halfh = diff(subs)/2;  
        x = bsxfun(@plus,points*halfh,midpoints);
        x = reshape(x,1,[]);
        
        fx = fun(x);
        fx1 = fx.*x;
        fx2 = fx.*x.^2;
        
        fx = reshape(fx,samples,[]);
        fx1 = reshape(fx1,samples,[]);
        fx2 = reshape(fx2,samples,[]);
        
        % Subintegrals.
        Ifxsubs = (w*fx) .* halfh;
        errsubs = (ew*fx) .* halfh;
        Ifxsubs1 = (w*fx1) .* halfh;
        Ifxsubs2 = (w*fx2) .* halfh;

        % Ifx and tol.
        Ifx = sum(Ifxsubs) + Ifx_ok;
        Ifx1 = sum(Ifxsubs1) + Ifx1_ok;
        Ifx2 = sum(Ifxsubs2) + Ifx2_ok;
        tol = max(atol,rtol*abs(Ifx));
        
        % determine the indices ndx of Ifxsubs for which the
        % errors are acceptable and remove those from subs
        ndx = find(abs(errsubs) <= (2/(b-a)*halfh*tol));
        subs(:,ndx) = [];
        if isempty(subs)
            break
        end
        
        % Update the integral.
        Ifx_ok = Ifx_ok + sum(Ifxsubs(ndx));
        Ifx1_ok = Ifx1_ok + sum(Ifxsubs1(ndx));
        Ifx2_ok = Ifx2_ok + sum(Ifxsubs2(ndx));
        
        % Quit if too many subintervals.
        nsubs = 2*size(subs,2);
        if nsubs > maxsubs
            warning('quad_moments: Reached the limit on the maximum number of intervals in use.');
            break
        end
        midpoints(ndx) = []; 
        subs = reshape([subs(1,:); midpoints; midpoints; subs(2,:)],2,[]); % Divide the remaining subintervals in half
    end
    
    % Scale moments
    m_0 = Ifx;
    m_1 = Ifx1./Ifx;
    m_2 = Ifx2./Ifx;
end