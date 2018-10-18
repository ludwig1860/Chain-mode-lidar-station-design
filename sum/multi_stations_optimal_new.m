warning('off');
%% user setting
user_vis_ang = 150;
trees_vis_rate = 0.90;


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


%% Multiple optimal station calculation
disp('Echo_2: searching the first optimal scan station.....');

[stations,vis_tree_id_station,...
    obscure_tree_id_station] = deal([]);

% The visibility map - for choosing the first station
vis_mask_overlay = zeros(size(img,1), size(img,2), numCircles);
parfor k = 1:numCircles
    
    obsvr_point = [x_centers_all(k), y_centers_all(k)];
    
    centers_others = centers_all;
    centers_others(k,:) = [];
    radius_others = radius_all;
    radius_others(k,:) = [];
    
    % calculate tangent points of circle given point
    [xtng_circles_obr, ytng_circles_obr] = pt_circ_tangent(centers_others, ...
        radius_others, obsvr_point);
    
    % visible area
    vis_mask = visible_mask(imageSize,obsvr_point,...
        xtng_circles_obr,ytng_circles_obr);
    vis_mask_overlay(:,:,k) = vis_mask;
    
end

sum_vis = round(100*sum(vis_mask_overlay,3)./numCircles);

% first station placement
[stt_optl] = station_pick(centers_all,radius_all, sum_vis);

first_optl_station = stt_optl;
stations = [stations;first_optl_station];

%% second, third and .... station

disp('Echo_3: searching the rest scan stations.....');

point_XY = first_optl_station;
centers = centers_all;
radius = radius_all;
intsec_point_pair_station = {};
tree_ang_station = [];
ii = 1;

while 1
    
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
%     
            % draw
%     color = [rand(1),rand(1),rand(1)];
%     for n=1:numCircles
%         numPairs = size(intsec_point_pair_xy{n},1);
%         if numPairs>0
%             for m=1:numPairs
%                 ctr_xy = centers_all(n,:);
%                 pts_01_xy = intsec_point_pair_xy{n}(m,1:2);
%                 pts_02_xy = intsec_point_pair_xy{n}(m,3:4);
%                 draw_sector(ctr_xy, pts_01_xy, pts_02_xy, color); hold on
%             end
%         end
%     end
    
    
    intsec_point_pair_station = ...
        [intsec_point_pair_station,intsec_point_pair_xy];
    tree_ang_station = [tree_ang_station, tree_ang];
    
    if ii==1
        % determine loop
        union_sec_ang = tree_ang;

        vis_circle_id = find(union_sec_ang >= user_vis_ang);
        occul_circle_id = find(union_sec_ang < user_vis_ang);
        obscure_tree_id_station{ii} = occul_circle_id;
        vis_tree_id_station{ii} = vis_circle_id;
        if numCircles - numel(vis_circle_id) < (1-trees_vis_rate)*numCircles
            break
        end
        
    else
        
        
        % ??????? ??????????????
        union_sec_ang = -1*ones(numCircles,1);
        for kk = 1:numCircles
            sectors_points = -1*ones(ii,4);
            for jj = 1:ii
                if ~isempty(intsec_point_pair_station{kk,jj})
                    sectors_points(jj,:)= ...
                        intsec_point_pair_station{kk,jj};
                end
            end
            
            sectors_points(sectors_points==-1) = [];
            sectors_points = reshape(sectors_points,...
                [numel(sectors_points)/4, 4]);
            
            if isempty(sectors_points)
                union_sec_ang(kk)=0;
            else
                ctr = centers_all(kk,:);
                radi = radius_all(kk);
                union_sec_ang(kk) = ...
                    union_sector(ctr, radi, sectors_points, imageSize);
            end
        end
        
        % determine loop
        
        vis_circle_id = find(union_sec_ang >= user_vis_ang);
        occul_circle_id = find(union_sec_ang < user_vis_ang);
        obscure_tree_id_station{ii} = occul_circle_id;
        vis_tree_id_station{ii} = vis_circle_id;
        
        obs_intsc = obscure_tree_id_station{1};

        for hh = 2:ii
            obs_intsc = intersect(obscure_tree_id_station{hh}, obs_intsc);
        end
        
        if length(obs_intsc)< (1-trees_vis_rate)*numCircles
            break
        end
        
    end
    
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
    
    
    % find the next optimal station
    vis_mask_overlay_sub = vis_mask_overlay(:,:,occul_circle_id);
    sum_vis_sub = ...
        round(100*sum(vis_mask_overlay_sub,3)./numel(occul_circle_id));
    
    
    [stt_optl] = station_pick(centers_all,radius_all,...
        sum_vis_sub, centers_all(occul_circle_id,:));
    
    point_XY = stt_optl;
    stations = [stations;point_XY];
    
    ii = ii+1;
end
stations = stations(1:end-1,:);
vis_tree_id_station(end) = [];
obscure_tree_id_station(end) = [];

%% visibility map for scan stations

disp('Echo_4: mapping visibility for scan stations.....');

numStations = size(stations,1);
vis_scan_overlay = zeros(imageSize(1),imageSize(2),numStations);
for t = 1:numStations
    vis_mask_scan = visible_mask(imageSize,stations(t,:),...
        xtng_circles,ytng_circles);
    vis_scan_overlay(:,:,t) = vis_mask_scan;
    
end

sum_vis_scan = sum(vis_scan_overlay,3);

h2 = figure; imagesc(sum_vis_scan);
viscircles(centers_all,radius_all,...
    'Color',[0.5 0.25 0]); hold on
axis equal
axis off
hcb = colorbar;
set(hcb,'YTick',0:1:numStations)
ylabel(hcb, 'Observation times');

id_circle_text = cellstr(num2str( transpose(1:1:numCircles)));
text(x_centers_all-15, y_centers_all, ...
    id_circle_text, 'color','r'); hold on

%% plot

for t = 1:numStations
    scatter(stations(t,1), stations(t,2),...
        800,'p','MarkerFaceColor',[1,0,0]);
    hold on
    text(stations(t,1)-20, stations(t,2)-40, ...
        ['Scan ',num2str(t)],...
        'color','k');
    hold on
end
% figure
for t = numStations:-1:1
    plot_fill_circle(centers_all(vis_tree_id_station{t},:), ...
        radius_all(vis_tree_id_station{t}), getColor(t-1,1,numStations));
    hold on
end
% axis square
hold off
print(h1,'manual_regular_10CD_vis_tree_02.tiff','-dtiff','-r300');
    print(h2,'manual_regular_10CD_stations_02.tiff','-dtiff','-r300');                                                                                                                
% clear
% LSY_fig('scan stations locations')