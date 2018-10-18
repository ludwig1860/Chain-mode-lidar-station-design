function [occul_circle_id, cumla_intsec_ang] = occulu_angle(TRI,x,y,fov_angle)
% occulu_angle function calculate the occulusion angle of a circle,
%       given a point
% Input:
%       TRI: m-by-3 matrix, point index belongs to a triangle
%            A row of TRI contains indices into the vectors x and y 
%            that define a single triangle
%       x: m-by-1 vector, x coordinat of point
%       y: m-by-1 vector, y coordinat of point
%       fov_angle: m-by-1 vector, fov angle of a circle, triangle angle
%
% Outpur:
%       occul_circle_id: m-by-1 vector, flag of tree - occulusion or not
%       cumla_intsec_ang: scalar, sum of all occulusion angle
% Date: 30/09/2018,  Linyuan Li

numCircles = length(TRI);

% tangent vector
tng_vec_x = x(1:end-1) - x(end);
tng_vec_y = y(1:end-1) - y(end);

% cartesian to polar coordinate
[tng_vec_theta,tng_vec_rho] = cart2pol(tng_vec_x,tng_vec_y);
tng_vec_theta = -rad2deg(tng_vec_theta);  % anti-clockwise direction

tri_tng_vec_theta = transpose(reshape(tng_vec_theta,2,numCircles));
tri_tng_vec_rho = transpose(reshape(tng_vec_rho,2,numCircles));

% universal polar coordinate
tri_tng_vec_theta(tri_tng_vec_theta<0) = ...
    tri_tng_vec_theta(tri_tng_vec_theta<0)+360;

% sort tangent vector angle
tri_tng_vec_theta = sort(tri_tng_vec_theta, 2);

% start and end vector
tng_theta_start = tri_tng_vec_theta(:,1); % min
tng_theta_end = tri_tng_vec_theta(:,2);  % max

% max_A - min_B
tng_theta_end_repl = repmat(tng_theta_end,1,numCircles);
tng_theta_start_repl = transpose(repmat(tng_theta_start,1,numCircles));

phi = abs(tng_theta_end_repl - tng_theta_start_repl); % max_A - min_B
verify = phi - phi';
phi(verify<0) = 100; % flag

phi(phi>180) = 360 - phi(phi>180);

%% three cases: no intersection, intersection: full overlap, intersection: partial overlap

omega = repmat(fov_angle,1,numCircles) +...
    transpose(repmat(fov_angle,1,numCircles));
omega(verify<0) = 0; % flag


decider = phi - omega;
decider(1:numCircles+1:end) = 100;  % diagram element

 % case 1(no intersection): max_A - min_B > (max_A-min_A)+(max_B-min_B)
 no_intsc_ind = find(decider>0);   
 
 % case 2 and 3  - intersection
 intsc_ind = find(decider<0);
 
 cumla_intsec_ang = 0;
 occul_circle_id = zeros(length(intsc_ind), 1);
 for i = 1:numel(intsc_ind)
     
      [intsc_r, intsc_c] = ind2sub([numCircles,numCircles],intsc_ind(i));
     
      sign_compr_01 = tng_theta_end(intsc_r)-tng_theta_end(intsc_c);
      sign_compr_02 = tng_theta_start(intsc_r)-tng_theta_start(intsc_c);
      
     % case 2 (full overlap): sign(max_A - max_B) = -sign(min_A - min_B)
     if sign(sign_compr_01) == -sign(sign_compr_02)
         
        intsc_ang_i = max([fov_angle(intsc_r), fov_angle(intsc_c)]);
        cumla_intsec_ang = intsc_ang_i + cumla_intsec_ang;
        
     % case 3 (partial overlap): sign(max_A - max_B) = sign(min_A - min_B)
     elseif sign(sign_compr_01) == sign(sign_compr_02)
        
        intsc_ang_i = min([abs(tng_theta_end(intsc_r) -...
            tng_theta_start(intsc_c)),...
            abs(tng_theta_end(intsc_c) - tng_theta_start(intsc_r))]);
        cumla_intsec_ang = intsc_ang_i + cumla_intsec_ang;
     end
     
     % which circle is occulused
     if tri_tng_vec_rho(intsc_r) > tri_tng_vec_rho(intsc_c)
         occul_circle_id(i) = intsc_r;
         
     elseif tri_tng_vec_rho(intsc_r) < tri_tng_vec_rho(intsc_c)
         occul_circle_id(i) = intsc_c;
         
     end
     
 end
 
 occul_circle_id = unique(occul_circle_id);

end

