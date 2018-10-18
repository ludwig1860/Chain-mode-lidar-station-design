warning('off');

%% stations
stations = [370,486;1737,486;2785,486;...
    370,1422;1737,1422;2785,1422;...
    370,2605;1737,2605;2785,2605];

%% Read tree position data

disp('Echo_0: reading trees dataset.....');

img = imread('manual_simu_plot _row_plant.png');  % mode 1: input image
% 
% trees_data = ...
%     xlsread('D:\linyuan_work\06_TLS station placement_UAV guide\Examples\Wanpeng_field_plot\plot_example_01_CD_10.xlsx');

% [centers_userCoor, radius_userCoor] = deal(trees_data(:,2:3), trees_data(:,4));
% times = 1;
% tScales = 100;


%% Detect tree position and tree DBH in this image

disp('Echo_1: detecting trees.....');
% 
% [img, centers_all,radius_all] = ...
%     userCoor_imgCoor(centers_userCoor, radius_userCoor, tScales, times);  % mode 2: input positions
imageSize = [size(img,1),size(img,2)];

[centers_all,radius_all] = imfindcircles(img,[6 50]);  % mode 1: input image

numCircles = length(centers_all);

[centers_all,radius_all] = deal(round(centers_all),...
    round(radius_all)); % pixel row and column
[x_centers_all, y_centers_all] = deal(centers_all(:,1),centers_all(:,2));

h1 = figure; imshow(img);
viscircles(centers_all,radius_all,...
    'Color',[0.5 0.25 0]); hold on

id_circle_text = cellstr(num2str( transpose(1:1:numCircles)));
text(x_centers_all-15, y_centers_all, ...
    id_circle_text, 'color','r','fontsize',15); hold on

%% 

disp('Echo_3: searching the rest scan stations.....');
centers = centers_all;
radius = radius_all;

intsec_point_pair_station = {};
tree_ang_station = [];

for gg=1:size(stations,1)
   
    point_XY = stations(gg,:);

    % calculate tangent points of circle given point
    
    [xtng_circles, ytng_circles] =...
        pt_circ_tangent(centers, radius, point_XY);
    
    % create triangle index
    TRI = [transpose(reshape(1:2*numCircles,2,numCircles)),...
        (2*numCircles+1)*ones(numCircles,1)];
    x = [reshape(transpose(xtng_circles),2*numCircles,1); point_XY(1)];
    y = [reshape(transpose(ytng_circles),2*numCircles,1); point_XY(2)];
%     triplot(TRI,x,y); hold on
    
    % fov angle of each circle calculation
    fov_angle = view_angle(TRI,x,y);
    
    % calculate visible tree angle and end points
    [tree_ang, intsec_point_pair_xy] =...
        tree_angle(TRI,x,y,fov_angle,centers_all,radius_all,imageSize);
    
    intsec_point_pair_station = ...
        [intsec_point_pair_station,intsec_point_pair_xy];
    tree_ang_station = [tree_ang_station, tree_ang];
    
    
        % draw
    color = [rand(1),rand(1),rand(1)];
    for n=1:numCircles
        numPairs = size(intsec_point_pair_xy{n},1);
        if numPairs>0
            for m=1:numPairs
                ctr_xy = centers_all(n,:);
                pts_01_xy = intsec_point_pair_xy{n}(m,1:2);
                pts_02_xy = intsec_point_pair_xy{n}(m,3:4);
                draw_sector(ctr_xy, pts_01_xy, pts_02_xy, color); hold on
            end
        end
    end
    
    
end
    %

%% plot

disp('Echo_4: mapping visibility for scan stations.....');


numStations = size(stations,1);
for t = 1:numStations
    scatter(stations(t,1), stations(t,2),...
        300,'p','MarkerFaceColor',getColor(t-1,1,numStations));
    hold on
    text(stations(t,1)-20, stations(t,2)-40, ...
        ['Scan ',num2str(t)],...
        'color','k');
    hold on
end

print(h1,'manual_regular_10CD_vis_tree.tiff','-dtiff','-r300');                                                                                                          
