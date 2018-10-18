function [vis_mask] = visible_mask(imageSize,obsvr_point,xtng_circles_pair,ytng_circles_pair )
% visible_mask function find the visible area given an observation point
% Input:
%       imageSize: 1-by-2 vector, image size
%       obsvr_point: 2-elements vector, the orw-column coordinate of
%                    observation point
%       xtng_circles_pair: m-by-2 matrix, row of tangent points for each circle,
%       ytng_circles_pair: m-by-2 matrix, column tangent points for each circle,
%
% Outpur:
%       vis_mask: size = imageSize, matrix, visible area mask
% Date: 02/10/2018,  Linyuan Li


% imageSize = [size(m_plot,1),size(m_plot,2)];
numCircles = length(xtng_circles_pair);
% obsvr_point = [640,640];
%% Reshape
xtng_circles_all = reshape(xtng_circles_pair',...
    [2*length(xtng_circles_pair), 1]);
ytng_circles_all = reshape(ytng_circles_pair',...
    [2*length(ytng_circles_pair), 1]);

%% Line equation
numTangePoints = length(xtng_circles_all);
lines_eq = zeros(numTangePoints,3);

for i=1:numTangePoints
    coefficients = polyfit([obsvr_point(1), xtng_circles_all(i)],...
        [obsvr_point(2), ytng_circles_all(i)], 1);
    a = coefficients (1);
    b = coefficients (2);
    lines_eq(i,:) = [a,-1,b] ;
end
% tx = 1:1:500; ty = a*tx+b;
% plot(tx,ty,'r','linewidth',5); hold on
% scatter([obsvr_point(1), xtng_circles(i)], [obsvr_point(2), ytng_circles(i)],50,'c','filled'); hold on

%% Intersection points of lines in image and image border
border_points = lineToBorderPoints(lines_eq,imageSize);
[pts_side1,pts_side2] = deal(border_points(:,1:2),border_points(:,3:4));

   % Choose the intersection point located at the same direction with tangent poing
[border_side1_theta,~] = cart2pol(pts_side1(:,1)-obsvr_point(1),...
    pts_side1(:,2)-obsvr_point(2)); % side 1 vector
border_side1_deg = rad2deg(border_side1_theta);

[border_side2_theta,~] = cart2pol(pts_side2(:,1)-obsvr_point(1),...
    pts_side2(:,2)-obsvr_point(2)); % side 2 vector
border_side2_deg = rad2deg(border_side2_theta);

[tng_theta,~] = cart2pol(xtng_circles_all-obsvr_point(1),...
    ytng_circles_all-obsvr_point(2));  % tangent vector
tng_deg = rad2deg(tng_theta);

    % determine the border point along the direction of tangent vector
pts_use = zeros(numTangePoints,2);
for j = 1:numTangePoints

    if abs(border_side1_deg(j) - tng_deg(j))<0.1
        pts_use(j,:) = pts_side1(j,:);
    elseif abs(border_side2_deg(j) - tng_deg(j))<0.1
        pts_use(j,:) = pts_side2(j,:);
    end

%     line([pts_use(j,1),obsvr_point(1)],[pts_use(j,2),obsvr_point(2)]);hold on
end


   % reshape
pts_use_pair_x = transpose(reshape(pts_use(:,1),[2,length(pts_use)/2]));
pts_use_pair_y = transpose(reshape(pts_use(:,2),[2,length(pts_use)/2]));

%% Visible mask
vis_mask = true(imageSize);
for k = 1:numCircles
    c = [xtng_circles_pair(k,1),pts_use_pair_x(k,1),pts_use_pair_x(k,2),xtng_circles_pair(k,2)];  % X - column
    r = [ytng_circles_pair(k,1),pts_use_pair_y(k,1),pts_use_pair_y(k,2),ytng_circles_pair(k,2)];  % Y - row
    BW = roipoly(vis_mask,c,r);
    vis_mask(BW==1) = 0;
end
% img = imread('station_placement_optimize_test_example02.png');
% figure;imshow(imoverlay(img, vis_mask, 'cyan'),'InitialMagnification',67); hold on
end
