function [ssx, y2] = dofit(y, Qd)
%DOFIT Do constrained least squares fit and reduce data for subsequent fits

% Fit y to design matrix in null space
y2 = Qd' * y;            % rotate y into that space: predicted y value
ssx = norm(y2)^2;        % sum of squares explained by fit
