function [order,coords,points] = perform_embedding(points,options)

% RUN_EMBEDDING   Perform Spectral Embedding
%
%   Perform Spectral Embedding
%
%   SYNTAX
%       [ORDER,COORDS,POINTS] = PERFORM_EMBEDDING(POINTS,OPTIONS)
%

% $Id: perform_embedding.m 4 2009-08-15 21:10:35Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-08-15 17:10:35 -0400 (Sam, 15 aoû 2009) $
% $Revision: 4 $

if options.disp_log
    disp('---- Running spectral embedding');
end

points = double(points);

global PointCloud

if nargin<3
    options.null = 0;
end

if ~isfield(options, 'mu_lapacian')
    options.mu_lapacian = 1/2;
end
mu_lapacian = options.mu_lapacian;

[npoints,space_dimension] = size(points);
% disp(['Embedding ', num2str(npoints), ' points in dimension ', num2str(space_dimension)]);
if options.disp_log
    disp(['Nb trials : ',num2str(npoints)]);
    disp(['Nb time samples : ',num2str(space_dimension)]);
end

%% Centering data
if options.CENTER_DATA % Centering data
    points = points - repmat(mean(points,2),1,size(points,2));
end

%% Normalizing data
if options.NORMALIZE_DATA % Sphering data
    points = points ./ repmat(sqrt(diag(points*points')),1,size(points,2));
end

%% Computing distances
if options.disp_log
    disp('Computing distances');
end
distances = compute_distances(points,options.distance);

%% Computing weights
if options.disp_log
    disp('Computing weights');
end
if options.USE_ADAPTIVE_SIGMA
    weights = compute_weights_smart(distances,options.ratio);
else
    weights = compute_weights(distances,options.sigma);
end

%% Computing eigenvalues/eigenvectors
if options.disp_log
    disp('Computing eigenvalues/eigenvectors');
end
sums = sum(weights);
if options.mu_lapacian
    weights = diag(1./sums.^(mu_lapacian))*weights*diag(1./sums.^(mu_lapacian));
    sums = sum(weights);
end

diago_method = 2;
switch diago_method
    case 1 % std method using the 'eig' function
        [vecs,vals] = eig(weights,diag(sums));
    case 2 % iterative method using the 'eigs' function
        options.tol = 1e-7; options.disp = 0; options.isreal = 1; options.issym = 1; options.maxit = 50;
        % options.tol = 1e-6; options.disp = 0; options.isreal = 1; options.issym = 1; options.maxit = 20;

        if options.verbose
            nb_vectors = 4;
        else
            nb_vectors = 2;
        end
        try
            options.disp = 0;
            [vecs,vals,flag] = eigs(double(weights),diag(double(sums)),nb_vectors,'LM',options);
        catch
            if options.disp_log
                warning('Could not estimate eigenvectors')
            end
            vals = zeros([1,nb_vectors]);
            vecs = zeros([npoints,nb_vectors]);
        end
otherwise
    error('Unknown diagonalisation method');
end

vals = diag(vals);
[vals,perm] = sort(-vals); vals = -vals; % sorting eigenvalues
vecs = vecs(:,perm); % applying the permutation to the eigenvectors

% get the coordinates in the small vector space
coords = vecs(:,2:end); % embedding in R^3 and skip first eigenvalue

if options.verbose
    %% scale coords with eigenvalues
    coords(:,1) = vals(2)*coords(:,1);
    coords(:,2) = vals(3)*coords(:,2);
    coords(:,3) = vals(4)*coords(:,3);
end

%% reordering with second eigenvalues
[trash,order] = sort(-coords(:,1));
weights_reordered = weights(order,order);
clear trash

if options.disp_log
    disp('-- Embedding done !');
end

if 0
    %% plotting eigenvalues
    nb_vp_display = 10;
    figure
    plot(vals(1:min(nb_vp_display,length(vals))));
    ylim([0,1.1]);
    title('Eigenvalues')
    figure
    imagesc(weights_reordered);
    title('Weights ordered with vp2')

    figure
    clf
    PointCloud = coords';
    ClickA3DPoint(PointCloud);
    hold on; % so we can highlight the clicked points without clearing the point cloud
    if options.CONNECT_EMBEDDED_POINTS
        plot3(PointCloud(1,:), PointCloud(2,:), PointCloud(3,:), 'k'); % visualize the point cloud
        % hold off
    end
    xlabel('vp 2');
    ylabel('vp 3');
    zlabel('vp 4');
    grid on
    title('Points embedded in 3D (Click on a point!)');
    camorbit(-45,-45)
    grid on

    if 0
        % plot coords separately
        figure
        subplot(3,1,1)
        plot(coords(order,1)');
        title('coord 1');
        subplot(3,1,2)
        plot(coords(order,2)');
        title('coord 2');
        subplot(3,1,3)
        plot(coords(order,3)');
        title('coord 3');
    end

    if 1
        figure
        imagesc(points(order,:))
        % colorbar
        xlabel('Time (ms)')
        ylabel('Trial')
        title('Data Ordered')
    end
end

