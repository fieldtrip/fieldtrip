function val = HDF5_DatasetLoad(gid, name, val0, options)

if ~exist('val0','var')
    val0 = [];
end
if ~exist('options','var')
    options = '';
end

try
    dsetid = H5D.open(gid, name);

    % NOTE: HDF5 stores contiguous muti-dimensional arrays in row-major order.
    % Matlab stores them in row-major order. We want to transpose the loaded data 
    % it back to Matlab's column-major storage order and thus get back the 
    % original array.  
    val = H5D.read(dsetid);    
    val = HDF5_PostProcessing(val, val0, options);

    % val = H5D.read(dsetid);
    H5D.close(dsetid);    
catch
    switch(class(val0))
        case 'char'
            val = '';
        case 'cell'
            val = {};
        otherwise            
            val = [];
    end
end
