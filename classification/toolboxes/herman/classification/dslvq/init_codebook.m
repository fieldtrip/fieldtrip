
function [codebook label] = init_codebook(data,labels,methodname,methodparams)

% INIT_CODEBOOK implements the process of initializing codebook for DSLVQ or LVQ 
%
%  Use as
%    1) [codebook label] = init_codebook(data,labels,methodname,methodparams)
%    2) As part of LVQ object (see lvq.m)
%
%  INPUT:
%         data   - training data (features)
%         labels - class labels
%         methodname, methodparams  - initialization method and its parameters
%                                    'randinit1','randinit2' - two approaches to random initialization
%                                    'eigenlininit' - utilising eigen structure of the data 
%                                           (projection of the data poiints sampled from a linear grid)    
%                                    'subclustinit' - subtractive clustering
%                                    'kmeansinit'   - k-means clustering
%                                    'fcminit'      - fuzzy c-means clustering
%         
%  OUTPUT:
%         codebook, label - trained codebook vectors with their labels
%
% SEE ALSO
%       dslvq_train.m
%       lvq1_train.m
%       lvq3_train.m
%       lvq.m

%  Pawel Herman, 2009


%  methodname    -  methodparams
%    'randinit1'     -  [numcods]
%    'randinit2'     -  [numcods]
%    'eigenlininit'  -  gridsize, e.g. for 3D case it could be [10 10 10]
%    'fcminit'       -  [numcods randout] 2nd param equal 1 implies random initialisation of labels 
%    'kmeansinit'    -  [numcods randout] 2nd param equal 1 implies random initialisation of labels
%    'subclustinit'  -  [radius  randout] 2nd param equal 1 implies random initialisation of labels

classes = unique(labels);

[ndata ndim] = size(data);

if nargin<3
    methodname   = 'randinit2';
    methodparams = 5*length(classes);
end
if nargin==3 || (nargin==4 && isempty(methodparams))
    switch methodname
        case 'randinit1',
            methodparams = 5*length(classes);
        case 'eigenlininit',
            methodparams = repmat([3],1,size(data,2));
        case 'fcminit',
            methodparams = [5 1];            
        case 'kmeansinit',
            methodparams = [5 1];
        case 'subclustinit',
            methodparams = [0.5 1];
        otherwise,
            methodname   = 'randinit2';
            methodparams = 5*length(classes);
    end
end


switch  methodname

    case 'randinit1',

        numcods =  methodparams(1);

        codebook = randn(numcods,ndim);
        codebook = (sqrtm(cov(data,1))*codebook' + mean(data)'*ones(1,numcods))';

        label = classes(floor(1+length(classes)*rand(size(codebook,1),1)));

    case 'randinit2',

        numcods =  methodparams(1);

        codebook = rand(numcods,ndim);

        max_val = max(data);     min_val = min(data);
        for i=1:ndim
            codebook(:,i) = (max_val(i) - min_val(i)) * codebook(:,i) + min_val(i);
        end

        label = classes(floor(1+length(classes)*rand(size(codebook,1),1)));

    case 'eigenlininit',

        gridsize =  methodparams;
        numcods = cumprod(gridsize);
        numcods = numcods(end);

        A = zeros(ndim);
        D = data;
        cdb = mean(D);
        D = D - ones(ndata,1)*cdb;

        for i=1:ndim,
            for j=i:ndim,
                c = D(:,i).*D(:,j);
                A(i,j) = sum(c)/length(c); A(j,i) = A(i,j);
            end
        end

        nev = ndim;

        [EV,S]   = eig(A);
        eigval   = diag(S);
        [aux,index]  = sort(eigval);
        eigval  = eigval(flipud(index));
        EV       = EV(:,flipud(index));
        EV       = EV(:,1:nev);
        eigval  = eigval(1:nev);

        EV = (sqrt(eigval)*ones(1,ndim)) .* EV ./ sqrt(ones(nev,1)*diag(EV*EV')');  %normalization

        grid = {};
        for i=1:ndim
            grid{1,i} = {linspace(-1,1,gridsize(i))};
        end
        init_codebook = myvariation(grid{:});
        %
        %         max_val = max(init_codebook);
        %         min_val = min(init_codebook);
        %         init_codebook =  (init_codebook - ones(numcods,1)*min_val) ./ (ones(numcods,1)*(max_val-min_val));
        %         init_codebook(abs(init_codebook)==Inf) = 0.5;
        %         init_codebook = (init_codebook - 0.5)*2;

        codebook = ones(numcods,1) * cdb;
        for i = 1:numcods
            for j = 1:nev
                codebook(i,:) = codebook(i,:) + init_codebook(i,j) * EV(:,j)';
            end
        end

        label = classes(floor(1+length(classes)*rand(size(codebook,1),1)));


    case 'kmeansinit',

        if  methodparams(2)
            [idx,C] =  kmeans(data, methodparams(1));
            codebook = C;
            label = classes(floor(1+length(classes)*rand(size(codebook,1),1)));
        else
            [idx,C] =  kmeans([data labels], methodparams(1));
            codebook = C(:,1:end-1);
            aux = abs(C(:,end)*ones(1,length(classes)) - (classes * ones(1,length(C(:,end))))');
            [val ind] = min(aux,[],2);
            label = classes(ind);
        end

    case 'fcminit',

        if  methodparams(2)
            [C] =  fcm(data, methodparams(1));
            codebook = C;
            label = classes(floor(1+length(classes)*rand(size(codebook,1),1)));
        else
            [C] =  fcm([data labels], methodparams(1));
            codebook = C(:,1:end-1);
            aux = abs(C(:,end)*ones(1,length(classes)) - (classes * ones(1,length(C(:,end))))');
            [val ind] = min(aux,[],2);
            label = classes(ind);
        end
        
        
    case 'subclustinit',

        if  methodparams(2)
            [C] =  subclustering(data, methodparams(1));
            codebook = C;
            label = classes(floor(1+length(classes)*rand(size(codebook,1),1)));
        else
            [C] =  subclustering([data labels], methodparams(1));
            codebook = C(:,1:end-1);
            aux = abs(C(:,end)*ones(1,length(classes)) - (classes * ones(1,length(C(:,end))))');
            [val ind] = min(aux,[],2);
            label = classes(ind);
        end
end