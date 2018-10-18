% test_for wan peng

plot_corner = [-20.107,-10.139;...
    -10.061,19.916;...
    20.364,10.030;...
    10.070,-20.157];

scan_wan = [0,0,;-9.077,-3.864;...
    -4.736,11.048;...
    8.710,7.334;...
    2.961,-11.540];
[x_scan, y_scan] = deal(scan_wan(:,1),scan_wan(:,2));

tree_pos = xlsread('D:\linyuan_work\06_TLS station placement_UAV guide\Examples\plot_example_01.xlsx');
times = 1;

x = tree_pos(:,2);
y = tree_pos(:,3);
r = tree_pos(:,4)/times;

%% translate coordinate

xn = x-min(x);
yn = y-min(y);
yn = max(yn)-yn;

max_plot_range = ceil(max([max(xn), max(yn)]));
min_plot_range = ceil(min([max(xn), max(yn)]));

min_frame_range_width = max_plot_range*100 +200;
min_frame_range_height = min_plot_range*100 +200;
img = 255*ones(min_frame_range_height, min_frame_range_width);

x_img = 100*xn + 100;
y_img = 100*yn + 100;
r_img = 100*r;


xn_scan = x_scan - min(x);
yn_scan = y_scan - min(y);
yn_scan = max(yn)-yn_scan;

x_scan_img = 100*xn_scan+100;
y_scan_img = 100*yn_scan+100;
scan_wan_img = [x_scan_img,y_scan_img];

%% display
figure;
imshow(img); hold on
viscircles([x_img,y_img],r_img,...
    'Color',[1 0 0]); hold on

scatter(x_scan_img,y_scan_img,100,'pk','filled')
centers_all = [x_img,y_img];
radius_all = r_img;

print('wanpeng_plot.png','-dpng');