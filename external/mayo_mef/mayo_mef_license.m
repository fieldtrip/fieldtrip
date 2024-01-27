function mayo_mef_license()
    % mayo_mef_license show license and setup mayo_mef toolbox

    % Copyright 2020-2023 Richard J. Cui. Created: Sat 03/21/2020 10:35:23.147 PM
    % $Revision: 0.5 $  $Date: Sun 10/08/2023 01:21:01.318 PM $
    %
    % Mayo Foundation for Medical Education and Research
    % Mayo Clinic St. Mary Campus
    % Rochester, MN 55905
    %
    % Email: richard.cui@utoronto.ca (permanent), Cui.Jie@mayo.edu (official)

    % show license
    % ------------
    type mayo_mef_license

    % setup mayo_mef toolbox
    % ----------------------
    % * add current m-file directory and subdirectories to MATLAB path
    addpath(genpath(fileparts(mfilename('fullpath'))));

end % function mayo_mef_license
