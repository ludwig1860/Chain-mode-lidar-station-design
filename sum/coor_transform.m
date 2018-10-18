function [theta, rho] = coor_transform(ctr,imgSize)
% coor_transform function convert image coordinate to polar coordinate, the
% basis of later coordinate is the center of circle
% Input:
%       ctr: 2-elements vector, center of circle
%       imgSize: 2-elements vector, size of image(m-by-n)
%
% Outpur:
%       theta: m-by-n matrix, angle 
%       rho: m-by-n matrix, length
%
% Date: 11/10/2018,  Linyuan Li

[numRow, numCol] = deal(imgSize(1),imgSize(2));

% image coordinate, (0,0) is basis
x_img = repmat(1:1:numCol,numRow,1);
y_img = repmat(transpose(1:1:numRow),1,numCol);

% image coordinate, (ctr(1), ctr(2)) is basis
x_ctr_coor = x_img - ctr(1);
y_ctr_coor = y_img - ctr(2);

% convert cartesian coordinate to polar coordinate
[theta, rho] = cart2pol(x_ctr_coor,y_ctr_coor);

theta = -rad2deg(theta);  % anti-clockwise direction

theta(theta<0) = theta(theta<0)+360;

end