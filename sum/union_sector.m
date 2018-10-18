function [union_sec_ang] = union_sector(ctr, radi, sectors_points, imageSize)
% union_sector function compute union of multiple sectors
% Input:
%       ctr: 2-elements vector, center of circle
%       radi: 2-elements vector, radius of circle
%       sectors_points: n-by-4 vector, n pair of sector end points (x1,y1,x2,y2)
%
% Outpur:
%       union_sec_ang: scalar, the union of sectors
%
% Date: 12/10/2018,  Linyuan Li

% if isempty(sectors_points)
%     union_sec_ang = 0;
% else

numPair = size(sectors_points,1);

if numPair == 1
    pts_01_xy = sectors_points(1:2);
    pts_02_xy = sectors_points(3:4);
    union_sec_ang = sec_angle(ctr, pts_01_xy, pts_02_xy);
else
    
    check_mask_all = false(imageSize(1),imageSize(2),numPair);
    for pp = 1:numPair
        sec_points = sectors_points(pp,:);
        
        [theta, rho] = coor_transform(ctr,imageSize);
        
        [check_mask] = checkin_sector(ctr,sec_points, theta, rho);
        check_mask_all(:,:,pp) = check_mask;
    end
    check_mask_final = sum(check_mask_all,3);
    sec_area = length(find(check_mask_final>0));
    union_sec_ang = 360*sec_area/(pi*radi.^2);

end