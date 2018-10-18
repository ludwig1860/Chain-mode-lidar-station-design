function [xtng_circles, ytng_circles] = pt_circ_tangent(centers, radius, point_XY)
% pt_circ_tangent function calculate the tangent point of a circle, given a point
% input:
%       centers: m-by-2 matrix, center of circle
%       radius: m-by-1 vector, radius of circle
%       points_XY: 1-by-2 matrix, position of point
% output:
%       xtng_circles: m-by-2 matrix, x coordinate of two tangent point
%       ytng_circles: m-by-2 matrix, y coordinate of two tangent point
%
% date: 28/09/2018, Linyuan Li

numCircles = length(centers);
[xtng_circles,ytng_circles] = deal(zeros(numCircles, 2));

for i=1:numCircles
    
    XY = transpose(point_XY);
    ctr = [centers(i,1); centers(i,2)];
    r = radius(i);
    
    % two tangent points of circle
    syms xh yh
    Q = [xh;yh];
    
        % Q lies on the circle.
    E1 = r^2 - [1 1]*(ctr - Q).^2 == 0;
    
        % Pythagoras applies
    E2 = [1 1]*(XY - ctr).^2 == r^2 + [1 1]*(Q - XY).^2;
    res = solve(E1,E2);
    
    % results
    xtng = double(vpa(res.xh));
    ytng = double(vpa(res.yh));
    
    xtng_circles(i,:) = real(xtng');
    ytng_circles(i,:) = real(ytng');
end

end