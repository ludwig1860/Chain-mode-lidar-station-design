function plot_fill_circle(c, r, color)
% plot_fill_circle function plot colored circles
% Input:
%       c: m-by-2 matrix, centers of circles
%       r: m-by-1 vector, radius of circles
%       color: 1-by-3 vector, color vector
%
% Date: 30/09/2018,  Linyuan Li

pos = [c-r 2*r 2*r];

for i = 1:size(pos,1)
    rr = rectangle('Position',pos(i,:),'Curvature',[1 1], 'FaceColor', color, 'Edgecolor','k');
    hold on
end

axis equal

end