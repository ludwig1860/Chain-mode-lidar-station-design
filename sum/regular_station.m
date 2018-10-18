
% given scan stations, compute visibility

warning('off');

%% Read tree position image
img = imread('manual_simu_plot _row_plant.png');
% img = 255*ones(3800,3800);
imageSize = [size(img,1),size(img,2)];

%% Detect tree position and tree DBH in this image

disp('Echo_1: detecting trees.....');

[centers_all,radius_all] = imfindcircles(img,[1 50]);

% [centers_all,radius_all] = ...
%     deal(importdata('centers_all_wp.mat'),...
%     importdata('radius_all_wp.mat'));

numCircles = length(centers_all);

[centers_all,radius_all] = deal(round(centers_all),...
    round(radius_all)); % pixel row and column
[x_centers_all, y_centers_all] = deal(centers_all(:,1),centers_all(:,2));


%% Visibility map
disp('Echo_2: computing the visibility.....');
scan_wan_img = [370,486;1737,486;2785,486;...
    370,1422;1737,1422;2785,1422;...
    370,2605;1737,2605;2785,2605];
[x_scan_w_img,y_scan_w_img] = deal(scan_wan_img(:,1),scan_wan_img(:,2));
numScans = length(x_scan_w_img);

% The visibility map - for choosing the first station
vis_mask_overlay = zeros(size(img,1), size(img,2), numScans);
parfor k = 1:numScans
    
    obsvr_point = [x_scan_w_img(k), y_scan_w_img(k)];
    
    % calculate tangent points of circle given point
    [xtng_circles_obr, ytng_circles_obr] = pt_circ_tangent(centers_all, ...
        radius_all, obsvr_point);
    
    % visible area
    vis_mask = visible_mask(imageSize,obsvr_point,...
        xtng_circles_obr,ytng_circles_obr);
    vis_mask_overlay(:,:,k) = vis_mask;
end

sum_vis = sum(vis_mask_overlay,3);

figure; imagesc(sum_vis);
viscircles(centers_all,radius_all,...
    'Color',[0.5 0.25 0]); hold on
axis equal
axis off
hcb = colorbar;
set(hcb,'YTick',0:1:5)
ylabel(hcb, 'Observation times');

id_circle_text = cellstr(num2str( transpose(1:1:numCircles)));
text(x_centers_all-15, y_centers_all, ...
    id_circle_text, 'color','r'); hold on

for t = 1:numScans
    scatter(scan_wan_img(t,1), scan_wan_img(t,2),...
        300,'p','MarkerFaceColor','k'); 
    hold on
    text(scan_wan_img(t,1)-20, scan_wan_img(t,2)-40, ...
        ['Scan ',num2str(t)],...
        'color','k');
    hold on
end
hold off
% print('wanpeng_plot_5scans.tiff','-dtiff','-r300');