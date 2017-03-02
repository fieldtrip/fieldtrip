function [y] = cellrows(x, rows)

% [Y] = CELLROWS(X, ROWS) outputs cell-array Y with the same dimensionality as X
% Each cell in Y only contains the rows ROWS from the original corresponding cell in X
% 
% X (and Y) should be linear cell-array(s) of matrices for which the number of rows 
% should be the same for all cells

y = cellrowselect(x, rows);
