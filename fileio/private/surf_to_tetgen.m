function surf_to_tetgen(filename, nodes, elements, element_types, attributes_nodes,attributes_elements)
% This function converts a triangulated mesh in FieldTrip format into a
% surface structure readable by Tetgen software

% Copyright Arno Jannsen, 2011

%Write the node file
node_filename = [filename '.node'];
precision_str = '%u   %.40e  %.40e  %.40e';

N_nodes = size(nodes,1);
if(nargin > 4)
    N_attributes = size(attributes_nodes,2);
    A = [nodes attributes_nodes];
    for i=1:N_attributes
        precision_str = [precision_str '  %.5e'];
    end
    precision_str = [precision_str '\n'];
else 
    N_attributes = 0;   
    A = nodes;
    precision_str = [precision_str '\n'];
end


fid = fopen(node_filename, 'w');
fprintf(fid, [num2str(N_nodes) ' 3 ' num2str(N_attributes) ' 0\n']);
for i=1:(N_nodes)
    fprintf(fid,precision_str,i,A(i,:));
end
fclose(fid);

%Write the poly file
poly_filename = [filename '.poly'];

N_elements = size(elements,1);

fid = fopen(poly_filename, 'w');
fprintf(fid, ['0  3  ' num2str(N_attributes) '  0\n']);
fprintf(fid, [num2str(N_elements) ' 1\n']);
for i=1:(N_elements)
    switch element_types(i)
        case 302
                fprintf(fid,['1  0  302  # ' num2str(i) '\n']);
                fprintf(fid,'3  %u  %u  %u\n',elements(i,:));
        otherwise
            ft_error('currently, only type 302 elements (triangles) are supported')
    end             
end
fprintf(fid,'0\n');
if(nargin > 5)
    N_regions = size(attributes_elements,1);
    if(size(attributes_elements,2) < 5)
        attributes_elements = [attributes_elements, -1*ones(N_regions,1)];
    end
    fprintf(fid,[num2str(N_regions) '\n']);
    for i=1:N_regions
        fprintf(fid,'%u  %.3f  %.3f  %.3f  %u  %.2f\n',i,attributes_elements(i,:));
    end
else
   fprintf(fid,'0\n');
end
fclose(fid);

end
