function [ tree_ang_final, intsec_point_pair_xy_final] =...
    tree_angle(TRI,x,y,fov_angle,centers_all,radius_all,imageSize)
% tree_angle function calculate the occulusion angle of a circle,
%       given a point
% Input:
%       TRI: m-by-3 matrix, point index belongs to a triangle
%            A row of TRI contains indices into the vectors x and y
%            that define a single triangle
%       x: (m+1)-by-1 vector, x coordinat of point
%       y: (m+1)-by-1 vector, y coordinat of point
%       fov_angle: m-by-1 vector, fov angle of a circle, triangle angle
%       centers_all:  m-by-2 matrix, center of circles
%
% Output:
%       tree_ang_final: n-by-1 cell, the sector of tree observed
%       intsec_point_pair_xy_final: n-by-1 cell, the intersection point between line
%                       and circle, the sector observed has two end point.
% Call:
%       sec_angle.m
% Date: 04/10/2018,  Linyuan Li

% Initialize data
numCircles = length(TRI);
tree_ang = cell(numCircles,1);
intsec_point_pair_xy = cell(numCircles,1);

[x_ctr, y_ctr] = deal(centers_all(:,1),centers_all(:,2));

x_tan = transpose(reshape(x(1:end-1),2,numCircles));
y_tan = transpose(reshape(y(1:end-1),2,numCircles));

% Tangent vector: view point -- tangent point
tng_vec_x = x(1:end-1) - x(end);
tng_vec_y = y(1:end-1) - y(end);

% cartesian to polar coordinate
[tng_vec_theta,tng_vec_rho] = cart2pol(tng_vec_x,tng_vec_y);
tng_vec_theta = -rad2deg(tng_vec_theta);  % anti-clockwise direction

tri_tng_vec_theta = transpose(reshape(tng_vec_theta,2,numCircles));
tri_tng_vec_rho = transpose(reshape(tng_vec_rho,2,numCircles));

% Rotate polar coordinate to fit for human
tri_tng_vec_theta(tri_tng_vec_theta<0) = ...
    tri_tng_vec_theta(tri_tng_vec_theta<0)+360;

% sort tangent vector angle
[tri_tng_vec_theta, sort_idx_i] = sort(tri_tng_vec_theta, 2);

% two matrices that are related and I need to sort one matrix,
%while sorting the corresponding matrix in the same way.
[tan_pts_x_order, tan_pts_y_order] = deal(zeros(numCircles,2));
for kk = 1:numCircles
    tan_pts_x_order(kk,:)  = x_tan(kk,sort_idx_i(kk,:));
    tan_pts_y_order(kk,:)  = y_tan(kk,sort_idx_i(kk,:));
end

% start and end vector
tng_theta_start = tri_tng_vec_theta(:,1); % min
tng_theta_end = tri_tng_vec_theta(:,2);  % max

[tan_pts_x_start,tan_pts_x_end] =...
    deal(tan_pts_x_order(:,1), tan_pts_x_order(:,2));
[tan_pts_y_start,tan_pts_y_end] =...
    deal(tan_pts_y_order(:,1), tan_pts_y_order(:,2));

% max_A - min_B
tng_theta_end_repl = repmat(tng_theta_end,1,numCircles);
tng_theta_start_repl = transpose(repmat(tng_theta_start,1,numCircles));

phi = abs(tng_theta_end_repl - tng_theta_start_repl); % max_A - min_B
verify = phi - phi';  % to verify which element is 'max_A - min_B'
phi(verify<0) = 100; % flag


%% three cases: no intersection, intersection: full overlap, intersection: partial overlap

omega = repmat(fov_angle,1,numCircles) +...
    transpose(repmat(fov_angle,1,numCircles));  % (max_A-min_A)+(max_B-min_B)
omega(verify<0) = 0; % flags

decider = phi - omega;
decider(1:numCircles+1:end) = 100;  % diagram element

% case 1(no intersection): max_A - min_B > (max_A-min_A)+(max_B-min_B)
no_intsc_ind = find(decider>=0);


% case 2 and 3  - intersection
intsc_ind = find(decider<0);

occul_circle_id = zeros(length(intsc_ind), 1);
[intsc_pt_pair_x_,intsc_pt_pair_y_] = deal(zeros(length(intsc_ind),2));
for i = 1:numel(intsc_ind)
    
    [intsc_r, intsc_c] = ind2sub([numCircles,numCircles],intsc_ind(i));
    
    
    sign_compr_01 = tng_theta_end(intsc_r)-tng_theta_end(intsc_c);
    sign_compr_02 = tng_theta_start(intsc_r)-tng_theta_start(intsc_c);
    
    % which circle is occulused
    if tri_tng_vec_rho(intsc_r) > tri_tng_vec_rho(intsc_c)
        occul_circle_id(i) = intsc_r;
        vis_temp_id = intsc_c;
        
    elseif tri_tng_vec_rho(intsc_r) < tri_tng_vec_rho(intsc_c)
        occul_circle_id(i) = intsc_c;
        vis_temp_id = intsc_r;
        
    end
    
    % case 2 (full overlap): sign(max_A - max_B) = -sign(min_A - min_B)
    if sign(sign_compr_01) == -sign(sign_compr_02)
        % do nothing
        
        % case 3 (partial overlap): sign(max_A - max_B) = sign(min_A - min_B)
    elseif sign(sign_compr_01) == sign(sign_compr_02)
        condi = abs(tng_theta_end(occul_circle_id(i)) - tng_theta_start(vis_temp_id))> ...
            abs(tng_theta_end(vis_temp_id) - tng_theta_start(occul_circle_id(i)));
        if condi
            if x(end)==tan_pts_x_end(vis_temp_id)
                slope = inf;
                intercpt = x(end);
                
            elseif y(end)==tan_pts_y_start(vis_temp_id)
                slope = 0;
                intercpt = y(end);
            else
                coefficients = polyfit([x(end),tan_pts_x_end(vis_temp_id)],...
                    [y(end),tan_pts_y_end(vis_temp_id)], 1);
                slope = coefficients (1);
                intercpt = coefficients (2);
            end
            [xout,yout] = linecirc(slope,intercpt,...
                centers_all(occul_circle_id(i),1),...
                centers_all(occul_circle_id(i),2),...
                radius_all(occul_circle_id(i)));

            
            if  all(isnan(xout))||all(isnan(yout))
                sec_ang = 0;
                [pts_01_xy,pts_02_xy] = deal([]);
            else
                
                d1 = pdist([x(end),y(end);xout(1),yout(1)],'euclidean');
                d2 = pdist([x(end),y(end);xout(2),yout(2)],'euclidean');
                
                if d1<d2
                    intsc_pt_pair_x_(i,1) = xout(1);
                    intsc_pt_pair_y_(i,1) = yout(1);
                elseif d1>d2
                    intsc_pt_pair_x_(i,1) = xout(2);
                    intsc_pt_pair_y_(i,1) = yout(2);
                end
                
                intsc_pt_pair_x_(i,2) = tan_pts_x_end(occul_circle_id(i));
                intsc_pt_pair_y_(i,2) = tan_pts_y_end(occul_circle_id(i));
                
                
                
                pts_01_xy = [intsc_pt_pair_x_(i,1),intsc_pt_pair_y_(i,1)];
                pts_02_xy = [intsc_pt_pair_x_(i,2),intsc_pt_pair_y_(i,2)];
                [sec_ang] = sec_angle(centers_all(intsc_c,:), pts_01_xy, pts_02_xy);
                
            end
            
            tree_ang{occul_circle_id(i)} = [tree_ang{occul_circle_id(i)};sec_ang];
            intsec_point_pair_xy{occul_circle_id(i)} = [intsec_point_pair_xy{occul_circle_id(i)};[pts_01_xy,pts_02_xy]];
            
            
        else
            if x(end)==tan_pts_x_start(vis_temp_id)
                slope = inf;
                intercpt = x(end);
                
            elseif y(end)==tan_pts_y_start(vis_temp_id)
                slope = 0;
                intercpt = y(end);
            else
                
                coefficients = polyfit([x(end),tan_pts_x_start(vis_temp_id)],...
                    [y(end),tan_pts_y_start(vis_temp_id)], 1);
                slope = coefficients (1);
                intercpt = coefficients (2);
            end
            [xout,yout] = linecirc(slope,intercpt,...
                centers_all(occul_circle_id(i),1),...
                centers_all(occul_circle_id(i),2),...
                radius_all(occul_circle_id(i)));
            
            if all(isnan(xout))||all(isnan(yout))
                sec_ang = 0;
                [pts_01_xy,pts_02_xy] = deal([]);
            else
                
                d1 = pdist([x(end),y(end);xout(1),yout(1)],'euclidean');
                d2 = pdist([x(end),y(end);xout(2),yout(2)],'euclidean');
                
                if d1<d2
                    intsc_pt_pair_x_(i,1) = xout(1);
                    intsc_pt_pair_y_(i,1) = yout(1);
                elseif d1>d2
                    intsc_pt_pair_x_(i,1) = xout(2);
                    intsc_pt_pair_y_(i,1) = yout(2);
                end
                
                intsc_pt_pair_x_(i,2) = tan_pts_x_start(occul_circle_id(i));
                intsc_pt_pair_y_(i,2) = tan_pts_y_start(occul_circle_id(i));
                
                pts_01_xy = [intsc_pt_pair_x_(i,1),intsc_pt_pair_y_(i,1)];
                pts_02_xy = [intsc_pt_pair_x_(i,2),intsc_pt_pair_y_(i,2)];
                [sec_ang] = sec_angle(centers_all(intsc_c,:), pts_01_xy, pts_02_xy);
            end
            tree_ang{occul_circle_id(i)} = [tree_ang{occul_circle_id(i)};sec_ang];
            intsec_point_pair_xy{occul_circle_id(i)} = [intsec_point_pair_xy{occul_circle_id(i)};[pts_01_xy,pts_02_xy]];
            
        end
        
    end
    
end

unique_occul_circle_id = unique(occul_circle_id);

% sector of visible circles
comp_vis_circle_id = 1:1:numCircles;
comp_vis_circle_id(unique_occul_circle_id) = [];


comp_vis_circle_x_tan = x_tan(comp_vis_circle_id,:);
comp_vis_circle_y_tan = y_tan(comp_vis_circle_id,:);

x_c_vis_cir = x_ctr(comp_vis_circle_id);
y_c_vis_cir = y_ctr(comp_vis_circle_id);

sec_ang = sec_angle([x_c_vis_cir,y_c_vis_cir],...
    [comp_vis_circle_x_tan(:,1),comp_vis_circle_y_tan(:,1)],...
    [comp_vis_circle_x_tan(:,2),comp_vis_circle_y_tan(:,2)]);


for ttt = 1:length(sec_ang)
    tree_ang{comp_vis_circle_id(ttt)}  = sec_ang(ttt);
    ptsv_01_xy = [x_tan(comp_vis_circle_id(ttt),1),y_tan(comp_vis_circle_id(ttt),1)];
    ptsv_02_xy = [x_tan(comp_vis_circle_id(ttt),2),y_tan(comp_vis_circle_id(ttt),2)];
    intsec_point_pair_xy{comp_vis_circle_id(ttt)} = [ptsv_01_xy,ptsv_02_xy];
end


%% the overlap of tree visible parts
intsec_point_pair_xy_final = cell(numCircles,1);
tree_ang_final = zeros(numCircles,1);
for mm = 1:numCircles
    if isempty(intsec_point_pair_xy{mm})
        continue
    end
    
    ctr = centers_all(mm,:);
    radi = radius_all(mm);
    numPair = size(intsec_point_pair_xy{mm},1);
    check_mask_all = false(imageSize(1),imageSize(2),numPair);
    
    % compute mask for checking points exists in circle sector or not
    
    if numPair == 1
        intsec_point_pair_xy_final{mm} = intsec_point_pair_xy{mm};
        tree_ang_final(mm) = tree_ang{mm};
    else
        
        
        for pp = 1:numPair
            sec_points = intsec_point_pair_xy{mm}(pp,:);
            
            [theta, rho] = coor_transform(ctr,imageSize);
            
            [check_mask] = checkin_sector(ctr,sec_points, theta, rho);
            check_mask_all(:,:,pp) = check_mask;
        end
        check_mask_final = logical(prod(check_mask_all,3));
        
        % fit sector
        in_idx = find(check_mask_final==1);
        if isempty(in_idx)
            intsec_point_pair_xy_final{mm} = [];
            tree_ang_final(mm) = 0;
        else
            
            [in_r,in_c] = ind2sub(imageSize, in_idx);
            
            tf = collinear([in_r,in_c], 1e-14); % determine colinear
            
            if tf==false
                K= convhull(in_c, in_r);
                edge_points = [in_c(K),in_r(K)];
                
                [two_end_points, sec_ang] = fit_sector(ctr, radi, edge_points);
            else
                sec_ang = 0;
                two_end_points = [];
            end
            intsec_point_pair_xy_final{mm} = two_end_points;
            tree_ang_final(mm) = sec_ang;
        end
    end
    
end

% figure
% imshow(check_mask_final);
% hold on
% scatter(in_c(K),in_r(K),'.g');
% hold on
% scatter(ctr(1),ctr(2),'.r');
% hold on
%
% draw_sector(ctr, end_pt_01, end_pt_02, [1 0 0]);


end


