%% Test the reciprocity data routine
infile = '../data/SPM/recipdata-gray.h5';
[data,compress,gridlocs,node_sizes,voxel_sizes] = hh_read_recipdata(infile);
decompress = gridlocs(:,1) + gridlocs(:,2) * node_sizes(1) + ...
    gridlocs(:,3) * node_sizes(2) * node_sizes(1) + 1;

%% Test the lead field matrices
sol = zeros(length(compress),1);
tic;
for nnode = 1:length(decompress)
    L = hh_leadfield_new(data,compress,node_sizes, ...
                         decompress(nnode));
% $$$     L = hh_leadfield(data,compress,node_sizes, gridlocs(nnode,1),gridlocs(nnode,2),gridlocs(nnode,3));
    sol(decompress(nnode)) = ~isempty(L);
end
toc;
fprintf('The number of valid leadfield matrices: %d\n',sum(sol > 0));

%% Test the results
X = zeros(node_sizes(1),node_sizes(2),node_sizes(3));
X(:) = sol;

%% Plot the solution image
imagesc(X(:,:,190))
