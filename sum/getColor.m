function color = getColor(value, min, max)
% getColor function set defaut colormap and extract a color
% 
% Input:
%       value: scalar, range: 0~max-1
%       min: scalar, 1
%       max: scalar, user defined
% Example:
%       color = getColor(4, 1, 10);

 set(groot,'DefaultFigureColormap', parula);

 cmap = get(groot,'defaultfigurecolormap');
 [m, ~] = size(cmap);
 
 row = round((value/(max-min))*(m-1)) + 1;
 color = cmap(row, :);  
 
end

