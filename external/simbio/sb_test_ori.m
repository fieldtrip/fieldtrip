function err = sb_test_ori(pnt,elem)
err = 1;
if(size(elem,2) == 4)
    det = sum(cross(pnt(elem(:,2),:)-pnt(elem(:,1),:),pnt(elem(:,4),:)-pnt(elem(:,1),:),2).*(pnt(elem(:,3),:)-pnt(elem(:,1),:)),2);
    if sum(det <= 0) > 0
        err = 0;
    end
elseif(size(elem,2) == 8)
    det1 = sum(cross(pnt(elem(:,6),:)-pnt(elem(:,1),:),pnt(elem(:,8),:)-pnt(elem(:,1),:),2).*(pnt(elem(:,5),:)-pnt(elem(:,1),:)),2);
    det2 = sum(cross(pnt(elem(:,3),:)-pnt(elem(:,1),:),pnt(elem(:,6),:)-pnt(elem(:,1),:),2).*(pnt(elem(:,2),:)-pnt(elem(:,1),:)),2);
    det3 = sum(cross(pnt(elem(:,8),:)-pnt(elem(:,1),:),pnt(elem(:,3),:)-pnt(elem(:,1),:),2).*(pnt(elem(:,4),:)-pnt(elem(:,1),:)),2);
    det4 = sum(cross(pnt(elem(:,8),:)-pnt(elem(:,3),:),pnt(elem(:,6),:)-pnt(elem(:,3),:),2).*(pnt(elem(:,7),:)-pnt(elem(:,3),:)),2);
    if sum(det1 <= 0) || sum(det2 <= 0) || sum(det3 <= 0) || sum(det4 <= 0)
        err = 0;
    end
else
    error('Invalid number of nodes per element!');
end
end