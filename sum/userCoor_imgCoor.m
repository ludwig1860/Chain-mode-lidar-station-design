function [scene_img, centers_imgCoor, radius_imgCoor] = ...
    userCoor_imgCoor(centers_userCoor, radius_userCoor,tScale, times)
% userCoor_imgCoor function converts user coordinate to image coordinate
% input:
%       centers_userCoor: m-by-2 matrix, center of circle in user coordinat
%       radius_userCoor: m-by-1 vector, radius of circle in user coordinate
%       times: scalar, scale convertion between user coordinate and meter
%       tScale: scalar, scale convertion between user coordinate and image
%               coordinate.
% output:
%       centers_imgCoor: m-by-2 matrix, center of circle in user coordinat
%       radius_imgCoor:  m-by-1 vector, radius of circle in user coordinate
%
% date: 28/09/2018, Linyuan Li

x = centers_userCoor(:,1);
y = centers_userCoor(:,2);
r = radius_userCoor/times;

% translate coordinate

xn = x-min(x);
yn = y-min(y);
yn = max(yn)-yn;

max_plot_range = ceil(max([max(xn), max(yn)]));
min_plot_range = ceil(min([max(xn), max(yn)]));

min_frame_range_width = max_plot_range*tScale +200;
min_frame_range_height = min_plot_range*tScale +200;
scene_img = 255*ones(min_frame_range_height, min_frame_range_width);

x_img = tScale*xn + 100;
y_img = tScale*yn + 100;
r_img = tScale*r;

centers_imgCoor = [x_img, y_img];
radius_imgCoor = r_img;

end


% %display
% figure;
% imshow(img); hold on
% viscircles([x_img,y_img],r_img,...
%     'Color',[1 0 0]); hold on
% 
% scatter(x_scan_img,y_scan_img,100,'pk','filled')
% centers_all = [x_img,y_img];
% radius_all = r_img;
% 
% print('wanpeng_plot.png','-dpng');