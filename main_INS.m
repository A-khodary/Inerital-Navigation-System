clear all
close all
clc

path(path,'functions')

n = 3100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Syncronization data parameters
%  GPS
deltaGPS = 0;
stepGPS = 100;
% IMU
deltaIMU  = 78;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GPS and IMU data
GPS_data
IMU_data
%% Angle conversion
r2d = 180/pi;
d2r = pi/180;

%% Filter selection
opt = 0;
if opt == 0
    filter =  'Extended Kalman filter'; 
end

%% Plot range
range = 100;

%% Simulation or actual experiment selection
% state = 1 simulation, state = 0  actual
state = 0;
%% System initial conditions
fIMU = 100;                       % Hz
TIMU = 1/fIMU;                    % Sample time [s].

fGPS = 1;                        % Hz
TGPS = 1/fGPS;                    % Sample time [s].

g = -9.81;                        % gravity[m/s^2].
%% Correlation time [s].
%  Rate gyros
tauR(1,1) = 626.8115;
tauR(2,1) = 6468.0515;
tauR(3,1) = 602.5784; 
% Accelerometers
tauA(1,1) = 1438.6558;
tauA(2,1) = 3807.8042;
tauA(3,1) = 1883.8307;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Weighting of the Kalman filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variances of the Q matrice
% AHRS - Attitude and heading reference system
varQ1(1,1) = 0.0001;
varQ1(2,1) = 0.0001;
varQ1(3,1) = 0.01;
varQ1(4,1) = 0.0001;
varQ1(5,1) = 0.0001;
varQ1(6,1) = 0.0001;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INS aided GPS - Inertial Navigation system aided GPS
varQ2(1,1) = 0.00059539;
varQ2(2,1) = 0.00062077;
varQ2(3,1) = 0.001261;
varQ2(4,1) = 1.5151e-6;
varQ2(5,1) = 2.181e-6;
varQ2(6,1) = 6.853e-6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variances of the R matrice
% AHRS - Attitude and heading reference system
varR1(1,1) = 6;
varR1(2,1) = 6;
varR1(3,1) = 6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INS aided GPS - Inertial navigation system aided GPS
varR2(1,1) = 2.591128492879236e-008^2;
varR2(2,1) = 2.591128492879236e-008^2;
varR2(3,1) = 0.00033^2;
varR2(4,1) = 0.5144^2;
varR2(5,1) = 0.5144^2;
varR2(6,1) = 0.5144^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P1 = eye(6);
P2 = eye(9)*1e-10;
% load P_matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sumDeltaz1= 0;
sumDeltaz2= 0;
magDec =  0.3604127;      % Magnetic declination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GPScont = 0;
for k=1:n
    %% System time
    k
    t = (k-1)*TIMU;
    t_v(k) = t;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% IMU reading
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Euler angles - IMU
    eangles.phi     = e_imu_v(1,k+deltaIMU);
    eangles.theta     = e_imu_v(2,k+deltaIMU);
    eangles.psi     = e_imu_v(3,k+deltaIMU);
    eangles.phi     = eangles.phi;
    eangles.theta   = -eangles.theta;
    eangles.psi     = -normalize_angle_f(eangles.psi +(magDec + pi/2),-pi);
    % Changing quaternion 
    qIMU = euler2quat_f(eangles);
    qIMU_v(:,k) =  qIMU;
    % Linear acceleration
    a = lin_accel_v(:,k+deltaIMU);
    a(2:3,1) = -a(2:3,1);  % Changing axis
    a_a = a/g;    
    a_a_v(:,k) = a_a;
    % Angular velocity
    omega_g = ang_vel_v(:,k+deltaIMU); 
    omega_g(2:3,1) = -omega_g(2:3,1); % Changing axis
    omega_g_v(:,k) = omega_g;
    % Roll and pitch angles - accelerometers
    thetaAccel =  normalize_angle_f(atan2(-a_a(1),sqrt(a_a(2)^2 + a_a(3)^2)),-pi); % Pitch accelerometer 
    thetaAccel_v(k) = thetaAccel;                                                  % Pitch accelerometer vector
    phiAccel   =  normalize_angle_f(atan2(-a_a(2),-a_a(3))+pi,-pi);                % Roll accelerometer
    phiAccel_v(k) = phiAccel;                                                      % Roll accelerometer vector

    % Euler angles normalization
    psiIMU          = normalize_angle_f(eangles.psi,-pi);   % Yaw IMU
    psiIMU_v(k)     = psiIMU;                               % Yaw IMU vector
    thetaIMU        = normalize_angle_f(eangles.theta,-pi); % Pitch IMU
    thetaIMU_v(k)   = thetaIMU;                             % Pitch IMU vector
    phiIMU          = normalize_angle_f(eangles.phi,-pi);   % Roll IMU
    phiIMU_v(k)     = phiIMU;                               % Roll IMU vector
    if k == 1
        qhat = euler2quat_f(eangles);
        eangleshat = eangles;
        
        psiHat   = eangleshat.psi;
        thetaHat = eangleshat.theta;
        phiHat   = eangleshat.phi;
        
        phiGyro    = phiHat;
        psiGyro    = psiHat;
        thetaGyro  = thetaHat;    
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% GPS reading
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Yaw angle obtained by the GPS
    if mod(k-1,stepGPS) == 0
        GPSstate = 1;           % There is GPS data
        GPScont = GPScont +1;   
        tGPS = t;
        tGPS_v(GPScont) =  t;
        % Reading GPS data
        psiGPS              = normalize_angle_f(heading_v(GPScont + deltaGPS),-pi);    % Yaw GPS
        deltapsi = psiGPS - eangles.psi;
        
        psiGPS_v(GPScont)   = psiGPS;                                                  % Yaw GPS vector
        latGPS              = lat_v(GPScont + deltaGPS);
        latGPS_v(GPScont)   = latGPS; 
        lonGPS              = lon_v(GPScont + deltaGPS);
        lonGPS_v(GPScont)   = lonGPS;
        altGPS              = alt_v(GPScont + deltaGPS);
        altGPS_v(GPScont)   = altGPS;
                   
        p = [latGPS;lonGPS;altGPS];
        p_v(:,GPScont) = p;

        if GPScont == 1
            deltaT = 0;
            dlat = 0;
            dlon = 0;
            dalt = 0;
            
            latGPSprevious = latGPS;
            lonGPSprevious = lonGPS;
            altGPSprevious = altGPS;
        else
            %deltaT = TGPS;
            deltaT = (stamp_gps(GPScont + deltaGPS) - stamp_gps(GPScont + deltaGPS -1))/10^9;
            dlat = (latGPS - latGPSprevious)/deltaT;
            dlon = (lonGPS - lonGPSprevious)/deltaT;
            dalt = (altGPS - altGPSprevious)/deltaT;
            
            latGPSprevious = latGPS;
            lonGPSprevious = lonGPS;
            altGPSprevious = altGPS;
        end

        lambda = latGPS;
        h =  altGPS ;

        Rlambda = Rlambda_m(lambda);
        Rphi = Rphi_m(lambda);

        vCalc = [(Rlambda + h)  0 0;
            0                                   (Rphi + h)*cos(lambda) 0;
            0                                   0 -1]* [dlat; dlon; dalt];
        vCalc_v (:,k) = vCalc;

        v = [vCalc(1);vCalc(2);vCalc(3)];

        v_v(:,GPScont) = v;

        pIG = p;
        vIG = v;  
    else 
        GPSstate = 0; % There is not GPS data
    end
    
    %% Filters   
    if k == 1
        %Bias initial conditions
        bg = zeros(3,1);  % Bias of the rate gyros
        ba = zeros(3,1);  % Bias of the accelerometers
        
        % State initial conditions
%         x1 = [phit; theta; psi; bg];
        
        x1 = [eangles.phi; eangles.theta; eangles.psi; bg]
        x2 = [p; v; ba];
    end
    % AHRS - Attitude and heading reference system
    
        [sys_a] = sysEuler_f(omega_g, x1, varQ1, varR1, g, tauR, TIMU);
        if opt == 0

            [x1p, P1p] = ekf_p (x1,P1,sys_a);
        end
        
        z1 = [eangles.phi; eangles.theta; eangles.psi + deltapsi];
        
        if opt == 0
            delta_z1 = normalize_angle_f(z1-sys_a.H*x1p,-pi);
            [x1,P1] = ekf_u(x1p, P1p, delta_z1, sys_a,1); 
        end
        

        eangleshat.phi     = normalize_angle_f(x1(1),-pi);
        eangleshat.theta   = normalize_angle_f(x1(2),-pi);
        eangleshat.psi     = normalize_angle_f(x1(3),-pi);
        bg = x1(4:6);

        % Quaternion 
        qhat = euler2quat_f(eangleshat); 
        qhat = qhat/norm(qhat);       
        % INS aided GPS - Inertial Navigation system aided GPS
        if GPSstate == 1
            [pINS,vINS] = GPS_IMU(x2(1:3), x2(4:6), qhat, a-ba, g, TIMU);
%           [p_INS,v_INS] = INS(x2(1:3,1),v, x1(1:4), a_a-ba, g,T);
%           [p_IG,v_IG] = GPS_IMU(p_IG,v_IG, x1(1:4), a-ba, g, T);
            x2(1:3) = pINS;
            x2(4:6) = vINS;
        end
        if GPSstate == 0;
            %%%%%%%%%%%%%%%%%%%%%% Position Reference System %%%%%%%%%%%%%%%%%%%
            sysP = syspos_m(p, v, qhat, varQ2, varR2, tauA, TIMU);
            if opt == 0
                [xp2, Pp2] = ekf_p (x2,P2,sysP);          % Prediction
            end
%            hhat2 = xp2(1:6);
            [pHat,vHat] = INS(x2(1:3),x2(4:6), qhat, a_a-ba, g,TIMU);
            hhat2 = [pHat;vHat];
            
            z2(1:6,1)  = [p; v];               % Position observation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Position reference system observation error
            deltaz2    = z2 - hhat2;
            deltaz2(1) = subtr_ang(z2(1),hhat2(1));
            deltaz2(2) = subtr_ang(z2(2),hhat2(2));
            
            sumDeltaz2 = deltaz2'*deltaz2*TIMU + sumDeltaz2;
            deltaz2_v(:,k) = deltaz2;   % Observation position error
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if opt == 0
                [x2,P2] = ekf_u(xp2, Pp2, deltaz2, sysP,1);   % Filtering 
            end
        end
        p_IG_v(:,k) = pIG;
        v_IG_v(:,k) = vIG;
        % Accelerometer bias estimated
        ba = x2(7:9);
        ba_v(:,k) = ba;
        latHat_v(k) = x2(1);
        lonHat_v(k) = x2(2);
        hHat_v(k)   = x2(3);
        vnHat_v(k) = x2(4);
        veHat_v(k) = x2(5);
        vdHat_v(k) = x2(6);
        % Rate gyros bias estimated
        bg_v(:,k) = bg;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        psiHat     = normalize_angle_f(eangleshat.psi,   -pi);
        thetaHat   = normalize_angle_f(eangleshat.theta, -pi);
        phiHat     = normalize_angle_f(eangleshat.phi,   -pi);
        psiHat_v(k)   = psiHat;
        thetaHat_v(k) = thetaHat;
        phiHat_v(k)   = phiHat;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if k >1
%             % Bias correction
%             phiGyro    = normalize_angle_f((omega_g(1) - bg(1))*TIMU + phiGyro,-pi);
%             thetaGyro  = normalize_angle_f((omega_g(2) - bg(2))*TIMU + thetaGyro,-pi);
%             psiGyro    = normalize_angle_f((omega_g(3) - bg(3))*TIMU + psiGyro,-pi);
            % Whithout bias correction
            phiGyro    = normalize_angle_f((omega_g(1))*TIMU + phiGyro,-pi);
            thetaGyro  = normalize_angle_f((omega_g(2))*TIMU + thetaGyro,-pi);
            psiGyro    = normalize_angle_f((omega_g(3))*TIMU + psiGyro,-pi);
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        phiGyro_v(k)   = phiGyro;
        thetaGyro_v(k) = thetaGyro;
        psiGyro_v(k)   = psiGyro;
        
       
end
p_GPS_v = p_v;
v_GPS_v = v_v;
plot_fig
drawnow
