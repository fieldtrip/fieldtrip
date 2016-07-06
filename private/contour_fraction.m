function [pnt] = contour_fraction(cnt1, cnt2, fraction)

% CONTOUR_FRACTION

ncnt = size(cnt1,1);

% determine the total length of the contour
tot_l = 0;
for i=1:ncnt
  tot_l = tot_l + pntdist(cnt1(i,:), cnt2(i,:));
end

frac_l = fraction * tot_l;

% propagate along the contour untill we get at the desired fraction
sum_l = 0;
for i=1:ncnt
  seg_l = pntdist(cnt1(i,:), cnt2(i,:));
  if (sum_l+seg_l)>=frac_l
    % the desired point lies on this segment
    la = frac_l - sum_l;
    vec = cnt2(i,:)-cnt1(i,:);
    pnt = cnt1(i,:) + la * vec/norm(vec);
    return
  else
    sum_l = sum_l + seg_l;
    sum_f = sum_l/tot_l;
  end
end

pnt = [nan nan nan];

