function shuffledArray = kshuffle(inputArray, DIM)

%
%   shuffledArray = kshuffle(inputArray, DIM)
%
%   Shuffles a matrix by row or column dimension.
%
%   Kian Torab
%   kian.torab@utah.edu
%   Department of Bioengineering
%   University of Utah
%   Version 1.1.0 - March 4, 2010

switch DIM
    case 'row'
        shuffledArray = inputArray(randperm(size(inputArray, 1)), :);
    case 'column'
        shuffledArray = inputArray(:, randperm(size(inputArray, 2)));
    otherwise
        disp('DIM is invalid. See HELP for more information.');
end