% +INTERNAL
%
% Files
%   allInteractions     - Given a cell array of factor names, create all pairwise combinations
%   designMatrix        - Extract a dummy encoded design matix from a linear model .
%   getAllTerms         - Returns a cell array with the names of all terms in the linear model.
%   gMatrix             - Generate a matrix in which each row corresponds to an effect, each column a value
%   interaction         - Create all interaction terms from two dummy coded design matrices.
%   inverseGammaPdf     - The inverse Gamma PDF.
%   mcIntegral          - Monte Carlo integration
%   nWayAnova           - Bayes Factor analysis for an N-Way Anova.
%   rouderS             - The S(g) function of Eq 9 in Rouder et al.
%   scaledInverseChiPdf - The Scaled Inverse Chi-Squared Distribution.
%   simulateLinearModel - Simulate data based on a linear model .
%   zeroSumConstraint   - Impose a zero-sum constraint on a dummy coded predictor matrix X
