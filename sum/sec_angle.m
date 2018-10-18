function [sec_ang] = sec_angle(ctr_xy, pts_01_xy, pts_02_xy)
% sec_angle function calculate the sector angle of a circle,
%       given centre and two end points
% Input:
%       ctr_xy: m-by-2 matrix, centres of circles
%       pts_01_xy: m-by-2 matrix, x,y of point A
%       pts_02_xy: m-by-2 matrix, x,y of point B
% Output:
%       sec_ang: m-by-1 vector, the sector angle of circle
% Date: 05/10/2018,  Linyuan Li

[ctr_x, ctr_y] = deal(ctr_xy(:,1), ctr_xy(:,2));
[pts_01_x, pts_01_y] = deal(pts_01_xy(:,1), pts_01_xy(:,2));
[pts_02_x, pts_02_y] = deal(pts_02_xy(:,1), pts_02_xy(:,2));

[vec_01_x, vec_02_x] = deal(pts_01_x-ctr_x, pts_02_x-ctr_x);
[vec_01_y, vec_02_y] = deal(pts_01_y-ctr_y, pts_02_y-ctr_y);

[vec_01_theta,~] = cart2pol(vec_01_x,vec_01_y);
[vec_02_theta,~] = cart2pol(vec_02_x,vec_02_y);

sec_ang = rad2deg(abs(vec_01_theta-vec_02_theta));

if sec_ang>180
    sec_ang = 360-sec_ang;
else
%     sec_ang = sec_ang;
end

end