function  [blocks, Xu, tstblocks] = set_PIC(x, dims, blocksize, indtype, visualize, xtest)
% SET_PIC     Set the inducing inputs and blocks for two dimensional input data
%
%  [blocks, Xu] = set_PIC(x, dims, blocksize, indtype, visualize)
%     x = input matrix, where each row is length two co-ordinate vector, 
%     dims = dimensions of the area that needs to be covered with inducing inputs 
%           and block. The format is, e.g., dims = [1    60     1    35];
%     blocksize = length of the block edge (blocks are square)
%     indtype = 'corners' (inducing inputs on the corners of the blocks)
%               'corners+1xside' (inducing inputs on the corners of the blocks and 
%                                 one on each block boundary)
%               'corners+2xside' (inducing inputs on the corners of the blocks and 
%                                 two on each block boundary)
%     visualize = 1 for plottin the data, blocks and inducing inputs
%                 0 for no plot
%
%     Returns:  
%     blocks = cell array with vectors appointing the block indeces
%     Xu = matrix of inducing inputs
%
%   [blocks, Xu, tstblocks] = set_PIC(x, dims, blocksize, indtype, visualize, xtest)
%     takes also test inputs 'xtest' and returns:
%     tstblocks = cell array with vectors appointing the block indeces for xtest
%
%     NOTE! The function works for test data only if the data is inside a square
%      area limeted by min(x) and max(x).
%
%   BUGS:
%    - The function works only if the data is distributed in the upper right quadrant
%      The reason for this is that the function was written for spatial data where this 
%      is the case. 

% Copyright (c) 2007-2008 Jarno Vanhatalo

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

% Find the indeces also for test data
if nargin >= 6
    also_test = 1;
    indextest={};
else
    also_test = 0;
end

ly = dims(2)/blocksize;
lx = dims(4)/blocksize;

b1 = linspace(0,dims(4),lx);
b2 = linspace(0,dims(2),ly);
index={}; indexsize = []; 
for i1=1:length(b1)-1
    for i2=1:length(b2)-1
        ind = 1:size(x,1);
        ind = ind(: , b1(i1)<=x(ind',1) & x(ind',1) < b1(i1+1));
        ind = ind(: , b2(i2)<=x(ind',2) & x(ind',2) < b2(i2+1));
        if ~isempty(ind)
            index{end+1} = ind';
            indexsize(end+1) = length(ind);
            if also_test == 1 
                indtest = 1:size(xtest,1);
                indtest = indtest(: , b1(i1)<=xtest(indtest',1) & xtest(indtest',1) < b1(i1+1));
                indtest = indtest(: , b2(i2)<=xtest(indtest',2) & xtest(indtest',2) < b2(i2+1));
                indextest{length(index)} = indtest';
            end 
        end
    end       
end
minblock = min(indexsize); maxblock = max(indexsize); avgblock = mean(indexsize);

switch indtype
  case 'corners'
    [xii,yii]=meshgrid(b1,b2);
    xii=[xii(:) yii(:)];
    qm=min(sqrt(bsxfun(@minus,x(:,1),xii(:,1)').^2+bsxfun(@minus,x(:,2),xii(:,2)').^2));
    qii=qm<=blocksize/2;
    Xu = unique(xii(qii,:), 'rows');
  case 'corners+1xside'
    b12 = linspace((b1(2)-b1(1))/2,dims(4)+(b1(2)-b1(1))/2,lx);
    b22 = linspace((b2(2)-b2(1))/2,dims(2)+(b2(2)-b2(1))/2,ly);

    [xii,yii]=meshgrid(b1,b2);
    xii=[xii(:) yii(:)];
    [xii2,yii2]=meshgrid(b12,b2);
    xii2=[xii2(:) yii2(:)];
    xii = [xii;xii2];
    [xii2,yii2]=meshgrid(b1,b22);
    xii2=[xii2(:) yii2(:)];
    xii = [xii;xii2];
    qm=min(sqrt(bsxfun(@minus,x(:,1),xii(:,1)').^2+bsxfun(@minus,x(:,2),xii(:,2)').^2));
    %qii=qm<=blocksize/4;sum(qii);
    qii=qm<=blocksize/4;
    Xu = unique(xii(qii,:), 'rows');
  case 'corners+2xside'
    b12 = linspace((b1(2)-b1(1))/3,dims(4)+(b1(2)-b1(1))/3,lx);
    b22 = linspace((b2(2)-b2(1))/3,dims(2)+(b2(2)-b2(1))/3,ly);
    b13 = linspace(2*(b1(2)-b1(1))/3,dims(4)+2*(b1(2)-b1(1))/3,lx);
    b23 = linspace(2*(b2(2)-b2(1))/3,dims(2)+2*(b2(2)-b2(1))/3,ly);
        
    [xii,yii]=meshgrid(b1,b2);
    [xii2a,yii2a]=meshgrid(b12,b2);
    [xii2b,yii2b]=meshgrid(b1,b22);
    [xii3a,yii3a]=meshgrid(b13,b2);
    [xii3b,yii3b]=meshgrid(b1,b23);
    xii = [xii(:) yii(:) ; xii2a(:) yii2a(:) ; xii2b(:) yii2b(:) ; xii3a(:) yii3a(:) ; xii3b(:) yii3b(:)];
    qm=min(sqrt(bsxfun(@minus,x(:,1),xii(:,1)').^2+bsxfun(@minus,x(:,2),xii(:,2)').^2));
    %qii=qm<=blocksize/4;sum(qii);
    qii=qm<=blocksize/6;
    Xu = unique(xii(qii,:), 'rows');
end

numu = length(Xu);
% $$$ D = 0;
% $$$ for i1 = 1:size(Xu,2);
% $$$     D = D + bsxfun(@minus,Xu(:,i1), Xu(:,i1)').^2;
% $$$ end
% $$$ meandist = mean(mean(sqrt(D)));

% Include the blocks with too few data points into larger ones
go = 1; i1=1;
while go == 1
    % If the block is too small include it in an other block
    if length(index{i1}) < (blocksize^2/4)
        meandist=[];
        if i1>1
            others = 1:i1-1;
        else
            others = [];
        end
        if i1 < length(index)
            others = [others i1+1:length(index)];
        end;
        for i2 = 1:length(others)
            D = 0;
            for i3 = 1:size(x,2);
                D = D + bsxfun(@minus,x(index{i1},i3), x(index{others(i2)},i3)').^2;
            end
            meandist(i2) = mean(mean(sqrt(D)));
        end
        sortedmeandist = sort(meandist);
        minmean = others(find(meandist == min(meandist)));
        if max(size(minmean)) > 1 
            minmean = minmean(1);
        end
        index{i1} =  sort([index{i1};index{minmean}]);
% $$$         % Set the block into the nearest (in averige) if the total block size is 
% $$$         % less than 1.5 x max blocksize
% $$$         if length(index{i1})+length(index{minmean}) < 1.25*maxblock
% $$$             index{i1} =  sort([index{i1};index{minmean}]);
% $$$             % Else take the second nearest and so on until found reasonable block
% $$$         elseif length(index{i1})+length(index{others(find(meandist == sortedmeandist(2)))}) < 1.25*maxblock 
% $$$             minmean = others(find(meandist == sortedmeandist(2)));
% $$$             index{i1} =  sort([index{i1};index{minmean}]);
% $$$             % Else take the third nearest and so on until found reasonable block
% $$$         elseif length(index{i1})+length(index{others(find(meandist == sortedmeandist(3)))}) < 1.25*maxblock 
% $$$             minmean = others(find(meandist == sortedmeandist(3)));
% $$$             index{i1} =  sort([index{i1};index{minmean}]);
% $$$             % Else take the fourth nearest and so on until found reasonable block
% $$$         else length(index{i1})+length(index{others(find(meandist == sortedmeandist(4)))}) < 1.25*maxblock 
% $$$             minmean = others(find(meandist == sortedmeandist(4)));
% $$$             index{i1} =  sort([index{i1};index{minmean}]);
% $$$         end
        if minmean==1
            index = {index{2:end}};
        elseif minmean==length(index)
            index = {index{1:end-1}};
        else
            index = {index{1:minmean-1} index{minmean+1:end}};
        end
        
        % Set also the test indeces
        if also_test == 1
            indextest{i1} =  sort([indextest{i1};indextest{minmean}]);
            if minmean==1
                indextest = {indextest{2:end}};
            elseif minmean==length(indextest)
                indextest = {indextest{1:end-1}};
            else
                indextest = {indextest{1:minmean-1} indextest{minmean+1:end}};
            end          
        end        
    else
        i1 = i1+1;
    end
    if i1 > length(index)
        go = 0;
    end
end

for i1 = 1:length(index)
    indexsize2(i1) = length(index{i1});
end
minblock = min(indexsize2); maxblock = max(indexsize2); avgblock = mean(indexsize2);


if visualize == 1
        
     col = {'b.','g.','c.','m.','y.','k.','b.','g.','c.','m.';
          'g.','c.','m.','y.','k.','b.','g.','c.','m.','b.';
          'c.','m.','y.','k.','b.','g.','c.','m.','b.','g.';
          'm.','y.','k.','b.','g.','c.','m.','b.','g.','c.';
          'y.','k.','b.','g.','c.','m.','b.','g.','c.','m.';
          'k.','b.','g.','c.','m.','b.','g.','c.','m.','y.';
          'b.','g.','c.','m.','b.','g.','c.','m.','y.','k.';};

    
    figure 
    set(gcf,'units','centimeters');
    hold on
    for i1=1:length(index)
        plot(x(index{i1},1),x(index{i1},2),'.', 'Color', [i1./length(index) abs(sin(i1)) 1/i1])
    end
    plot(Xu(:,1),Xu(:,2), 'r*', 'MarkerSize', 8, 'LineWidth', 2)

    tmpy = 1:dims(2);
    tmpx = 1:dims(4);
    for i1 = 1:length(b1)
        plot(repmat(b1(i1), 1, length(tmpy)), tmpy,'k--')
    end
    for i1 = 1:length(b2)
        plot(tmpx, repmat(b2(i1), 1, length(tmpx)),'k--')
    end
    axis equal
    xlim([0 dims(4)])
    ylim([0 dims(2)])
    set(gcf, 'pos', [11.3896   15.9581   15   25])
    
    S = sprintf('%d blocks, with %.f data points in averige and %d/%d at most/least \n %d inducing inputs.', ...
                length(index), avgblock, maxblock, minblock, numu);
    H = text(-5, -5, S);
    set(H, 'FontSize', 14);
end


if also_test == 1
    for i=1:length(indextest)
        if min(size(indextest{i})) == 0
            indextest{i} = [];
        end        
    end
    tstblocks = indextest;
end

blocks = index;
