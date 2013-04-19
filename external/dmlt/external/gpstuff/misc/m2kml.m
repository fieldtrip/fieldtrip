% M2KML   Converts GP prediction results to a KML file
%
%  Input: 
%     input_file  - Name of the .mat file containing GP results 
%     cellsize    - Size of the cells in meters
%     output_file - Name of the output file without the file extension! 
%                   If output is not given the name of the input file is used
%                   and the results are written to <input_name>.kmz
% 
%  Assumes that the .mat-file contains (atleast) the following variables:
%    Ef - Logarithm of the relative risk to be displayed
%    X1 - Grid of cells. The size of the data-array is determined from this. 
%    xxii - Indexes of non-zero elements in the cell grid
%  
function m2kml(input_file,cellsize,output_file)
    if nargin < 3
        name = input_file(1:end-4);
        zip_name = [name, '.kmz'];
        output_file = [name, '.kml'];
    else        
        zip_name = [output_file, '.kmz'];
        output_file = [output_file, '.kml'];
    end
    use ge_toolbox
    
    addpath /proj/bayes/software/jmjharti/maps 
    
    load(input_file)
    
    % Form the data for output
    data = zeros(size(X1));
    data = data(:);
    data(:) = NaN;
    data(xxii) = exp(Ef);
    
    cellsize = 5000;
    % is there some mistake in coordinates?
    x=(3050000:cellsize:3750000-1)+cellsize;
    y=(6600000:cellsize:7800000-1)+cellsize;
    
    data = reshape(data,length(y),length(x));
    
    
    % create kml code for polygon overlay
    output = ge_imagesc_ykj(x, y, data, ...
                                 'altitudeMode', 'clampToGround', ...
                                 'transparency', 'ff');

    % Write the KML string to output file
    ge_output(output_file, [output]);
    
    % Zip the KML file
    zip(zip_name,output_file);
    movefile([zip_name,'.zip'],zip_name);
    