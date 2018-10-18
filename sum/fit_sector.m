function [two_end_points, sec_ang] = fit_sector(ctr, radi, edge_points)
% checkin_sector function check points exists in circle sector or not.
% Input:
%       ctr: 2-elements vector, center of circle
%       radi: 2-elements vector, radius of circle
%       edge_points: n-by-2 matrix, edge points
%
% Output:
%       two_end_points: 2-by-2 matrix, two end points coordinate
%       sec_ang: scalar, angle of sector
%
% Date: 11/10/2018,  Linyuan Li

[edge_pt_x,edge_pt_y] = deal(edge_points(:,1),edge_points(:,2));

% image coordinate, (ctr(1), ctr(2)) is basis
x_ctr_coor = edge_pt_x - ctr(1);
y_ctr_coor = edge_pt_y - ctr(2);

% convert cartesian coordinate to polar coordinate
[theta, rho] = cart2pol(x_ctr_coor,y_ctr_coor);
theta = -rad2deg(theta);  % anti-clockwise direction
theta(theta<0) = theta(theta<0)+360;
rho = round(rho);

% find end points
radi_idx = rho>=radi-2 & rho<=radi+1;
theta_candi = theta(radi_idx);

if length(find(radi_idx==1))<2
    sec_ang = 0;
    two_end_points = [];
else
    
    min_idx = find(theta_candi==min(theta_candi));
    max_idx = find(theta_candi==max(theta_candi));
    
    if max(theta_candi) - min(theta_candi)<=180
        sec_ang = max(theta_candi) - min(theta_candi);
    else
        sec_ang = (360-max(theta_candi))+(min(theta_candi)-0);
    end
    
    edge_points_candi = edge_points(radi_idx,:);
    end_pt_01 = edge_points_candi(min_idx(1),:);
    end_pt_02 = edge_points_candi(max_idx(1),:);
    
    two_end_points = [end_pt_01, end_pt_02];
end

end