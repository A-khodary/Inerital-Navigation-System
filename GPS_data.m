%% Load GPS data
gps_mod = load('malaga-urban-dataset-extract-04_all-sensors_GPS.txt');
% gps_select = [stamp latitude longitude altitude heading]
gps_select = [gps_mod(:,1) gps_mod(:,2) gps_mod(:,3) gps_mod(:,4) gps_mod(:,8)];
n_gps = size(gps_select,1);

stamp_gps = gps_select(:,1);

lat_v = gps_select(:,2);

lon_v = gps_select(:,3);

alt_v = gps_select(:,4);

heading_v =gps_select(:,5);
