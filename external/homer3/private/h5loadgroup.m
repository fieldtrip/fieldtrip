function data = h5loadgroup(loc)

data = struct();

num_objs = H5G.get_num_objs(loc);

%
% Load groups and datasets
%
for j_item = 0:num_objs-1
    objtype = H5G.get_objtype_by_idx(loc, j_item);
    objname = H5G.get_objname_by_idx(loc, j_item);
    
    if objtype == 0
        % Group
        name = regexprep(objname, '.*/', '');
        
        if ~isempty(regexp(name,'^[a-zA-Z].*'))
            group_loc = H5G.open(loc, name);
            try
                sub_data = h5loadgroup(group_loc);
                H5G.close(group_loc);
            catch exc
                H5G.close(group_loc);
                rethrow(exc);
            end
            if isempty(path_parts)
                data = setfield(data, name, sub_data);
            else
                data = sub_data;
                return
            end
        end
        
    elseif objtype == 1
        % Dataset
        name = regexprep(objname, '.*/', '');        
        if ~isempty(regexp(name,'^[a-zA-Z].*'))
            dataset_loc = H5D.open(loc, name);
            try
                % NOTE: HDF5 stores contiguous muti-dimensional arrays in row-major order.
                % Matlab stores them in row-major order. We want to transpose the loaded data
                % it back to Matlab's column-major storage order and thus get back the
                % original array.
                sub_data = H5D.read(dataset_loc, 'H5ML_DEFAULT', 'H5S_ALL','H5S_ALL','H5P_DEFAULT');
                sub_data = HDF5_PostProcessing(sub_data);

                H5D.close(dataset_loc);
            catch exc
                H5D.close(dataset_loc);
                rethrow(exc);
            end
            
            sub_data = fix_data(sub_data);
            
            data = setfield(data, name, sub_data);
        end
    end
    
end


% ---------------------------------------------------------------------
function data = fix_data(data)
% Fix some common types of data to more friendly form.

if isstruct(data)
    fields = fieldnames(data);
    if length(fields) == 2 & strcmp(fields{1}, 'r') & strcmp(fields{2}, 'i')
        if isnumeric(data.r) & isnumeric(data.i)
            data = data.r + 1j*data.i;
        end
    end
end

if isnumeric(data) & ndims(data) > 1
    % permute dimensions
    data = permute(data, fliplr(1:ndims(data)));
end

