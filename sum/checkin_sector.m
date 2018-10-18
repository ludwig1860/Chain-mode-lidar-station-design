function [check_mask] = checkin_sector(ctr,sec_points, theta, rho)
% checkin_sector function check points exists in circle sector or not.
% Input:
%       ctr: 2-elements vector, center of circle
%       sec_points: 1-by-4 vector, two sector end points (x1,y1,x2,y2)
%       theta: m-by-n matrix, angle 
%       rho: m-by-n matrix, length
%
% Outpur:
%       check_mask: m-by-n matrix, logical matrix, pixels in sector are
%       assigned as 1.
%
% Date: 11/10/2018,  Linyuan Li

[ctr_x, ctr_y] = deal(ctr(:,1), ctr(:,2));
[pts_01_xy, pts_02_xy] = deal(sec_points(1:2),sec_points(3:4));
[pts_01_x, pts_01_y] = deal(pts_01_xy(:,1), pts_01_xy(:,2));
[pts_02_x, pts_02_y] = deal(pts_02_xy(:,1), pts_02_xy(:,2));

[vec_01_x, vec_02_x] = deal(pts_01_x-ctr_x, pts_02_x-ctr_x);
[vec_01_y, vec_02_y] = deal(pts_01_y-ctr_y, pts_02_y-ctr_y);

% polar coordinate, basis is ctr
[vec_01_theta,radi] = cart2pol(vec_01_x,vec_01_y);
[vec_02_theta,radi] = cart2pol(vec_02_x,vec_02_y);


vec_01_theta = -rad2deg(vec_01_theta);  % anti-clockwise direction
vec_01_theta(vec_01_theta<0) = vec_01_theta(vec_01_theta<0)+360;

vec_02_theta = -rad2deg(vec_02_theta);  % anti-clockwise direction
vec_02_theta(vec_02_theta<0) = vec_02_theta(vec_02_theta<0)+360;

sort_theta = sort([vec_01_theta,vec_02_theta]);
max_limit = sort_theta(2);
min_limit = sort_theta(1);

% checkin pixels inside sector or not
if max_limit - min_limit <= 180
    checkin_idx = find(theta<=max_limit & theta>=min_limit & rho<=radi);
elseif max_limit - min_limit > 180
    checkin_idx_01 = find(theta>=0 & theta<=min_limit & rho<=radi);
    checkin_idx_02 = find(theta<=360 & theta>=max_limit & rho<=radi);
    checkin_idx = [checkin_idx_01;checkin_idx_02];
end

check_mask = false(size(theta));
check_mask(checkin_idx) = 1;

end