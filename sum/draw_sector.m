function draw_sector(ctr_xy, pts_01_xy, pts_02_xy, color)
% draw_sector function draw a sector 
%
% Input:
% Input:
%       ctr_xy: 1-by-2 vector, centres of circles
%       pts_01_xy: 1-by-2 vector, x,y of point A
%       pts_02_xy: 1-by-2 vector, x,y of point B
%       color: 1-by-3 vector, color vector
%

% Date:08/10/2018, Linyuan Li

% end angle
[x0,y0] = deal(ctr_xy(1), ctr_xy(2));
r = pdist([ctr_xy;pts_01_xy],'euclidean');

[pts_01_x, pts_01_y] = deal(pts_01_xy(1), pts_01_xy(2));
[pts_02_x, pts_02_y] = deal(pts_02_xy(1), pts_02_xy(2));

[vec_01_x, vec_02_x] = deal(pts_01_x-x0, pts_02_x-x0);
[vec_01_y, vec_02_y] = deal(pts_01_y-y0, pts_02_y-y0);

[angle1,~] = cart2pol(vec_01_x,vec_01_y);
[angle2,~] = cart2pol(vec_02_x,vec_02_y);

% select angle less than 180 degrees
[angle1, angle2] = deal(rad2deg(angle1),rad2deg(angle2));
angle1(angle1<0) = angle1(angle1<0)+360;
angle2(angle2<0) = angle2(angle2<0)+360;

if abs(angle2-angle1)<180
    % do nothing
elseif angle2-angle1>180
    angle2 = angle2-360;
elseif angle1-angle2>180
    angle1  = angle1-360;
end

t = linspace(angle1/180*pi, angle2/180*pi, 100);

x = x0+r*cos(t);
y = y0+r*sin(t);

h = fill([x0 x],[y0 y],color);
set(h,'facealpha',.3);

end

