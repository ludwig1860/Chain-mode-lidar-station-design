function [stt_optl] = station_pick(centers_all,radius_all,...
    sum_vis, obscured_tree_pos, varargin)
% station_pick function find optimal scan station
% input:
%       centers_all: m-by-2 matrix, center of circle
%       radius_all: m-by-1 vector, radius of circle
%       sum_vis: k-by-p vector, visibility map
%       obscured_tree_pos: n-by-2 matrix, obscured tree centers
% output:
%       stt_optl: 1-by-2 matrix, position of station
%
% date: 02/10/2018, Linyuan Li

numCircles = length(centers_all);
imgSize = size(sum_vis);

% envelope
x = centers_all(:,1);
y = centers_all(:,2);
k = convhull(x,y);

BW = roipoly(sum_vis, x(k), y(k));
sum_vis(BW==0)=0;

% buffer area
x = 1:imgSize(1);
y = 1:imgSize(2);
[xx,yy] = meshgrid(x,y);

for i = 1:numCircles
    %circlePixels is a 2D "logical" array.
    mask = hypot(xx - centers_all(i,1),...
        yy - centers_all(i,2)) <= 3.7*radius_all(i);
    sum_vis(mask) = 0;
end

% position

if nargin==3
    stt_optl_ind = find(sum_vis==max(sum_vis(:)));
    [stt_optl_r,stt_optl_c] = ind2sub(size(sum_vis), stt_optl_ind);
    stt_optl = [stt_optl_c(1),stt_optl_r(1)];
    
elseif nargin==4
    % convex hull of obscured tree
    [tree_x,tree_y] = deal(obscured_tree_pos(:,1),...
        obscured_tree_pos(:,2));
    K = convhull(tree_x,tree_y);
    BW = roipoly(sum_vis,tree_x(K),tree_y(K));
    sum_vis_sub = sum_vis.*BW;
    %
    
    stt_optl_ind = find(sum_vis_sub==max(sum_vis_sub(:)));
    numCandid = numel(stt_optl_ind);
    [stt_optl_r,stt_optl_c] = ind2sub(size(sum_vis), stt_optl_ind);
    stt_optl = [stt_optl_c(round(numCandid/2)),...
        stt_optl_r(round(numCandid/2))];
    
end

end