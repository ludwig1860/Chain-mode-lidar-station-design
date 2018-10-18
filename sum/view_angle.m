function [fov_angle] = view_angle(TRI,x,y)
% view_angle function calculate the occulusion angle of a circle, given a point
% Input:
%       TRI: m-by-3 matrix, point index belongs to a triangle
%            A row of TRI contains indices into the vectors x and y that define a single triangle
%       x: m-by-1 vector, x coordinat of point
%       y: m-by-1 vector, y coordinat of point
%
% Outpur:
%       view_angle: m-by-1 vector, the fov angle of each circle
% Date: 29/09/2018,  Linyuan Li

numCircles = length(TRI);

% tangent vector
tng_vec_x = x(1:end-1) - x(end);
tng_vec_y = y(1:end-1) - y(end);

tri_tng_vec_x = transpose(reshape(tng_vec_x,2,numCircles));
tri_tng_vec_y = transpose(reshape(tng_vec_y,2,numCircles));

% angle calculation 
fov_angle = atan2d(tri_tng_vec_x(:,1).*tri_tng_vec_y(:,2)-tri_tng_vec_y(:,1).*tri_tng_vec_x(:,2),...
    tri_tng_vec_x(:,1).*tri_tng_vec_x(:,2)+tri_tng_vec_y(:,1).*tri_tng_vec_y(:,2));

fov_angle = abs(fov_angle);

end