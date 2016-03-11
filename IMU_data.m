%% Load IMU data
path(path,'functions')

imu_mod = load('malaga-urban-dataset-extract-04_all-sensors_IMU.txt');
% imu_select = [stamp orientation_euler angular_velocity linear_acceleration ]
imu_select = [imu_mod(:,1) imu_mod(:,11:13) imu_mod(:,5:7) imu_mod(:,2:4)];
n_imu = size(imu_select,1);

stamp_imu = imu_select(:,1);

e_imu_v = imu_select(:,2:4)';

ang_vel_v = imu_select(:,5:7)';

lin_accel_v  = imu_select(:,8:10)';

