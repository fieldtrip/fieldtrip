This is the Snirf package downloaded from https://github.com/fNIRS/snirf_homer3 on 6 August 2020,
corresponding to commit ff9b8e0.

The MATLAB files have been reorganized to make it easier to handle the path with the 
FieldTrip `ft_hastoolbox` function, and the example files have been removed.

    cp ~Downloads/snirf_homer3/Snirf/*        homer3
    cp ~Downloads/snirf_homer3/Utils/*.m      homer3/private
    cp ~Downloads/snirf_homer3/Utils/Hdf5/*   homer3/private

Furthermore, some missing dependencies have been added from the Util directory of https://github.com/BUNPC/Homer3.
