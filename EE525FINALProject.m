
%% 2D Kalman filter
% This code needs the Px, Py, and Time data already in the workspace
% do this by double-clicking DATA2D.mat
clc
close all

last = size(Px); % simulation bound
W = 1; % process noise
deltaT = 1; % delta t initialization
Hk = [1 0 0 0;0 0 1 0]; % measurement-state matrix

% can try different Vx, Vy values, and Pk values
XkApriori = [Px(1) 0  Py(1) 0]'; % estimation of the state vector prior to k observation
% PkApriori = [[1 1 1 1] % error covariance matrix prior to k observation
%              [1 1 1 1]
%              [1 1 1 1]
%              [1 1 1 1]];
PkApriori = [[10000 0 0 0] % error covariance matrix prior to k observation
             [0 10000 0 0]
             [0 0 10000 0]
             [0 0 0 10000]];

Rk = [1 0 ; 0 1]; % measurement noise

Xk = zeros(4, last(1) - 1)'; % estimation of the state vector after k observation

for k = 1:(last(1) - 1) % iterate through time, measurement and time update
    %get measurements and deltaT
    zk = [Px(k);Py(k)]; % observation / measurement
    deltaT = (Time(k+1) - Time(k)) * 3600*24; % convert to seconds

    %Measurement Update
    K_k = PkApriori*Hk' * (Hk*PkApriori*Hk' + Rk)^-1; % Kalman Gain
    Xk(k,:) = XkApriori + K_k*(zk - (Hk*XkApriori));
    Pk = (eye(4,4) - K_k*Hk)*PkApriori; % error covariance after k observation

    %Time update; recalculate phi and Q based on delta T
    Phi = [[1 deltaT 0 0] % state transition matrix
           [0 1 0 0]
           [0 0 1 deltaT]
           [0 0 0 1]];

    Q = [[(W/3)*deltaT^3 (W/2)*deltaT^2 0 0] % model process covariance matrix
        [(W/2)*deltaT^2 W*deltaT 0 0]
        [0 0 (W/3)*deltaT^3 (W/2)*deltaT^2]
        [0 0 (W/2)*deltaT^2 W*deltaT]];

    XkApriori = Phi*Xk(k,:)';
    PkApriori = Phi*Pk*Phi' + Q;


end

%Do not plot the last sample data; it is ignored due to the lack of a
%deltaT
time = Time(1:end-1);
px = Px(1:end-1);
py = Py(1:end-1);

subplot(2,2,1);
plot(time,px, time, Xk(:,1));
ylim([87650 87750]);
title("Position X");
legend("Measured", "Filtered");

subplot(2,2,2);
plot(time,py, time,Xk(:,3));
title("Position Y");
legend("Measured", "Filtered");

subplot(2,2,3);
plot(time,Xk(:,2));
title("Velocity x");
legend("Filtered");

subplot(2,2,4);
plot(time,Xk(:,4));
title("Velocity y");
legend("Filtered");

%% 3D Kalman Filter

%This code needs the Px3D, Py3D, Pz3D and Time data already in the workspace
clc
close all

last = size(Px3D);
W = 1;
deltaT = 1;
Hk = [[1 0 0 0 0 0]
      [0 0 1 0 0 0]
      [0 0 0 0 1 0]];

XkApriori = [Px3D(1) 0  Py3D(1) 0 Pz3D(1) 0]';
PkApriori = [[1 0 0 0 0 0]
             [0 1 0 0 0 0]
             [0 0 1 0 0 0]
             [0 0 0 1 0 0]
             [0 0 0 0 1 0]
             [0 0 0 0 0 1]];

Rk = [[1 0 0]
      [0 1 0]
      [0 0 1]];

Xk = zeros(6, last(1) - 1)';

for k = 1:(last(1) - 1)
    %get measurements and deltaT
    zk = [Px3D(k);Py3D(k);Pz3D(k)];
    deltaT = (Time3D(k+1) - Time3D(k)) * 3600*24; % convert to seconds

    %Measurement Update
    K_k = PkApriori*Hk' * (Hk*PkApriori*Hk' + Rk)^-1;
    Xk(k,:) = XkApriori + K_k*(zk - (Hk*XkApriori));
    Pk = (eye(6,6) - K_k*Hk)*PkApriori;

    %Time update; recalculate phi and Q based on delta T
    Phi = [[1 deltaT 0 0 0 0] 
           [0 1 0 0 0 0]
           [0 0 1 deltaT 0 0]
           [0 0 0 1 0 0]
           [0 0 0 0 1 deltaT]
           [0 0 0 0 0 1]];

    Q = [[(W/3)*deltaT^3 (W/2)*deltaT^2 0 0 0 0]
        [(W/2)*deltaT^2 W*deltaT 0 0 0 0]
        [0 0 (W/3)*deltaT^3 (W/2)*deltaT^2 0 0]
        [0 0 (W/2)*deltaT^2 W*deltaT 0 0]
        [0 0 0 0 (W/3)*deltaT^3 (W/2)*deltaT^2]
        [0 0 0 0 (W/2)*deltaT^2 W*deltaT]];

    XkApriori = Phi*Xk(k,:)';
    PkApriori = Phi*Pk*Phi' + Q;


end

%Do not plot the last sample data; it is ignored due to the lack of a
%deltaT
time3D = Time3D(1:end-1);
px3D = Px3D(1:end-1);
py3D = Py3D(1:end-1);
pz3D = Pz3D(1:end-1);

subplot(2,3,1);
plot(time3D,px3D, time3D, Xk(:,1));
% ylim([87650 87750]);
title("Position X");
legend("Measured", "Filtered");

subplot(2,3,2);
plot(time3D,py3D, time3D,Xk(:,3));
title("Position Y");
legend("Measured", "Filtered");

subplot(2,3,3);
plot(time3D,pz3D, time3D, Xk(:,5));
% ylim([87650 87750]);
title("Position Z");
legend("Measured", "Filtered");


subplot(2,3,4);
plot(time3D,Xk(:,2));
title("Velocity x");
legend("Filtered");

subplot(2,3,5);
plot(time3D,Xk(:,4));
title("Velocity y");
legend("Filtered");

subplot(2,3,6);
plot(time3D,Xk(:,6));
title("Velocity z");
legend("Filtered");


