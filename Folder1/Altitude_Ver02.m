%% Flat-Earth 6-Dof Simulation + 4th RK
% Just inner loop, to hovering

clc; clear all; close all;
%% Initialization
Rad2Deg  = 180/pi;
Deg2Rad = pi/180;

dt = 0.01; % Time step
Tf = 60; % Total simulation time

g = [0 0 9.80665]';
l = 0.225; % [m] dist. between CoM and rotor
MASS  = 0.468; % [kg]
Ix = 4.856*10^(-3) ;
Iy = 4.856*10^(-3);
Iz = 8.801*10^-3;  
J = diag([Ix, Iy, Iz]);

k_M = 1.140*10^(-7); % [N/(rad^2/s^2)] Thrust coefficient
k_F = 2.980*10^(-6); % [Nm/(rad^2/s^2)] Thrust coefficient
base_throttle = - MASS * g(3);

Prop.MASS = MASS;
Prop.J = J;
Prop.g = g;
Prop.k_F = k_F;
Prop.k_M = k_M;
Prop.l = l;

I_err_alt = 0;
ki_alt = 0.01;
kp_alt = 0.1;
kd_alt = 0.3;

P = zeros(3, Tf/dt + 1); % Pos. vec. in NED coord. sys.
V = zeros(3, Tf/dt + 1); % Vel. vec. in body fixed coord. sys.
EA = zeros(3, Tf/dt + 1); % Euler angles: phi, theta, and psi.
AngVel = zeros(3, Tf/dt + 1); % p, q, r
Q = zeros(4, Tf/dt + 1); % Quaternion
V_n = zeros(3, Tf/dt + 1); % Vel. vec. in NED coord. sys.

% Initial values
P0 = [0 0 0]';
P(:,1) = P0;

V0 = [0 0 0]';
V(:,1) = V0;

EA0 = [0 0 0]'*Deg2Rad;
EA(:,1) = EA0;

q0 = cos(EA0(1)/2)*cos(EA0(2)/2)*cos(EA0(3)/2) + sin(EA0(1)/2)*sin(EA0(2)/2)*sin(EA0(3)/2);
q1 = sin(EA0(1)/2)*cos(EA0(2)/2)*cos(EA0(3)/2) - cos(EA0(1)/2)*sin(EA0(2)/2)*sin(EA0(3)/2);
q2 = cos(EA0(1)/2)*sin(EA0(2)/2)*cos(EA0(3)/2) + sin(EA0(1)/2)*cos(EA0(2)/2)*sin(EA0(3)/2);
q3 = cos(EA0(1)/2)*cos(EA0(2)/2)*sin(EA0(3)/2) - sin(EA0(1)/2)*sin(EA0(2)/2)*cos(EA0(3)/2);

Q0 = [q0 q1 q2 q3]';
Q(:,1) = Q0;

AngVel0 = [0 0 0]*Deg2Rad';
AngVel(:,1) = AngVel0;

% To Store Value
Alt_err = zeros(1, Tf/dt + 1);
Rol_err = zeros(1, Tf/dt + 1);
Pit_err = zeros(1, Tf/dt + 1);
store_cmd_alt = zeros(1, Tf/dt + 1);

%% Main Part(Iteration)
for m = 1 : Tf/dt
    Cn2b = [Q(1,m)^2 + Q(2,m)^2 - Q(3,m)^2 - Q(4,m)^2      2*(Q(2,m)*Q(3,m) + Q(1,m)*Q(4,m))                 2*(Q(2,m)*Q(4,m) - Q(1,m)*Q(3,m));
            2*(Q(2,m)*Q(3,m) - Q(1,m)*Q(4,m))              Q(1,m)^2 - Q(2,m)^2 + Q(3,m)^2 - Q(4,m)^2     2*(Q(3,m)*Q(4,m) + Q(1,m)*Q(2,m));
            2*(Q(2,m)*Q(4,m) + Q(1,m)*Q(3,m))              2*(Q(3,m)*Q(4,m) - Q(1,m)*Q(2,m))             Q(1,m)^2 - Q(2,m)^2 - Q(3,m)^2 + Q(4,m)^2];
    
    EA(1,m) = atan2(Cn2b(2,3), Cn2b(3,3));
    EA(2,m) = -asin(Cn2b(1,3));
    EA(3,m) = atan2(Cn2b(1,2), Cn2b(1,1));
    Prop.Cn2b = Cn2b;

    V_n(:,m) = (Cn2b)'*V(:,m);

    %% Sensor modeling (Navigaion)
    Alt_sensor = P(3, m);
    
    %% Guidance
    Alt_ref = -1;

    %% Control for Altitude
    Alt_err(m) = Alt_ref - Alt_sensor;

    I_err_alt = I_err_alt + Alt_err(m) * dt; % 적분 오차 업데이트
    if m > 1
        D_err_alt = (Alt_err(m) - Alt_err(m-1)) / dt; % 오차 미분
    else
        D_err_alt = 0; 
    end

    cmd_alt = kp_alt * Alt_err(m) + kd_alt * D_err_alt + ki_alt * I_err_alt;
    

    store_cmd_alt(m) = cmd_alt;
    %% Main Control
    u1 = base_throttle + cmd_alt; % F_z
    u2 = 0; % L
    u3 = 0; % M
    u4 = 0; % N
    
    U = [u1, u2, u3, u4]';

    % F_z를 그대로 dynamics에 제어 입력을 사용
    % (기존 코드에서 모터의 각속도를 이용해 원하는 추력을 올바르게 생성하는 것을 확인하였음)
    CONTROL = U; 
    
    %% Numerical integration (dynamics, kinematics)
    STATE = [P(:,m); V(:,m); AngVel(:,m); Q(:,m)];
    
    DEL_X_1 = Dynamics_GyL_Ver01(STATE,                 CONTROL, Prop)*dt;
    DEL_X_2 = Dynamics_GyL_Ver01(STATE + (1/2)*DEL_X_1, CONTROL, Prop)*dt;
    DEL_X_3 = Dynamics_GyL_Ver01(STATE + (1/2)*DEL_X_2, CONTROL, Prop)*dt;
    DEL_X_4 = Dynamics_GyL_Ver01(STATE + DEL_X_3,       CONTROL, Prop)*dt;

    STATE_AFTER = STATE + (1/6)*(DEL_X_1 + 2*DEL_X_2 + 2*DEL_X_3 + DEL_X_4); % 4-th RK
    
    % Update the State
    P(:,m + 1) = STATE_AFTER(1:3, 1);
    V(:,m + 1) = STATE_AFTER(4:6, 1);
    AngVel(:,m + 1) = STATE_AFTER(7:9, 1);
    Q(:,m + 1) = STATE_AFTER(10:13, 1);
    
    if m == Tf/dt 
        Cn2b = [Q(1,m+1)^2 + Q(2,m+1)^2 - Q(3,m+1)^2 - Q(4,m+1)^2      2*(Q(2,m+1)*Q(3,m+1) + Q(1,m+1)*Q(4,m+1))                 2*(Q(2,m+1)*Q(4,m+1) - Q(1,m+1)*Q(3,m+1));
            2*(Q(2,m+1)*Q(3,m+1) - Q(1,m+1)*Q(4,m+1))              Q(1,m+1)^2 - Q(2,m+1)^2 + Q(3,m+1)^2 - Q(4,m+1)^2     2*(Q(3,m+1)*Q(4,m+1) + Q(1,m+1)*Q(2,m+1));
            2*(Q(2,m+1)*Q(4,m+1) + Q(1,m+1)*Q(3,m+1))              2*(Q(3,m+1)*Q(4,m+1) - Q(1,m+1)*Q(2,m+1))             Q(1,m+1)^2 - Q(2,m+1)^2 - Q(3,m+1)^2 + Q(4,m+1)^2];
    
        EA(1,m+1) = atan2(Cn2b(2,3), Cn2b(3,3));
        EA(2,m+1) = -asin(Cn2b(1,3));
        EA(3,m+1) = atan2(Cn2b(1,2), Cn2b(1,1));

        V_n(:,m+1) = (Cn2b)'*V(:,m+1);
    end
end

Alt_err(end) = Alt_err(end-1); % 임시로 설정
%% Result
Realtime = 0:dt:Tf;
LineWidth1 = 1.5;

figure;
subplot(3,1,1);
plot(Realtime, P(1,:), 'b', 'LineWidth',LineWidth1);
ylabel('X[m]')
grid on;
title('Position in X direction');

subplot(3,1,2);
plot(Realtime, P(2,:), 'r', 'LineWidth',LineWidth1);
ylabel('Y[m]')
grid on;
title('Position in Y direction');

subplot(3,1,3);
plot(Realtime, P(3,:), 'r', 'LineWidth',LineWidth1);
ylabel('Z[m]')
xlabel('Time[sec]')
grid on;
title('Position in Z direction');

figure;
plot(Realtime, Alt_err, 'b', 'LineWidth', LineWidth1);
ylabel('Alt err[m]')
grid on;
title('Altitude Error');

figure;
plot(Realtime, store_cmd_alt, 'b', 'LineWidth', LineWidth1);
grid on;
ylabel('CMD ALT')
title('Commanded Altitude');

figure;
plot(Realtime, P(3,:), 'r', 'LineWidth', LineWidth1);
set(gca, 'YDir','reverse'); % y축 반전
ylabel('Z[m]');
xlabel('Time[sec]');
grid on;
title('Altitude over Time');
