%% Flat-Earth 6-Dof Simulation + 4th RK
% Outer loop design.

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

Rmat_inv = [-1/(4*k_F), -1/(2*sqrt(2)*l*k_F),  1/(2*sqrt(2)*l*k_F),  1/(4*k_M);
                -1/(4*k_F), -1/(2*sqrt(2)*l*k_F), -1/(2*sqrt(2)*l*k_F), -1/(4*k_M);
                -1/(4*k_F),  1/(2*sqrt(2)*l*k_F), -1/(2*sqrt(2)*l*k_F),  1/(4*k_M);
                -1/(4*k_F),  1/(2*sqrt(2)*l*k_F),  1/(2*sqrt(2)*l*k_F), -1/(4*k_M)];

w1 = 0;
w2 = 0;
w3 = 0;
w4 = 0;
w = zeros(4, Tf/dt + 1); % speed of rotor
w(:,1) = [w1; w2; w3; w4];

Prop.MASS = MASS;
Prop.J = J;
Prop.g = g;
Prop.k_F = k_F;
Prop.k_M = k_M;
Prop.l = l;

% position control gain
kp_pos = 1; 
ki_pos = 0;
kd_pos = 2;

% altitude control gain
I_err_alt = 0; % 적분 오차 초기화
kp_alt = 0.1;
ki_alt = 0.001;
kd_alt = 0.4;

% rolling control gain
I_err_rol = 0; % 적분 오차 초기화
kp_rol = 0.1;
ki_rol = 0.001;
kd_rol = 0.4;

% pitching control gain
I_err_pit = 0; % 적분 오차 초기화
kp_pit = 0.1;
ki_pit = 0.001;
kd_pit = 0.4;

% yawing control gain
I_err_yaw = 0; % 적분 오차 초기화
kp_yaw = 0.1;
ki_yaw = 0.001;
kd_yaw = 0.3;

P = zeros(3, Tf/dt + 1); % Pos. vec. in NED coord. sys.
V = zeros(3, Tf/dt + 1); % Vel. vec. in body fixed coord. sys.
EA = zeros(3, Tf/dt + 1); % Euler angles: phi, theta, and psi.
AngVel = zeros(3, Tf/dt + 1); % p, q, r
Q = zeros(4, Tf/dt + 1);

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


Pos_err = zeros(2, Tf/dt + 1);
Alt_err = zeros(1, Tf/dt + 1);
Rol_err = zeros(1, Tf/dt + 1);
Pit_err = zeros(1, Tf/dt + 1);
Yaw_err = zeros(1, Tf/dt + 1);
store_cmd_pos_roll = zeros(1,Tf/dt + 1);
store_cmd_pos_pitch = zeros(1,Tf/dt + 1);
store_cmd_y = zeros(1,Tf/dt + 1);
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
    % Sampling ratio를 고려하지 않고 진행
    
    Pos_sensor = P(1:2, m);
    Alt_sensor = P(3, m);
    
    
    EA_sensor = EA(1:2,m);
    Yaw_sensor = EA(3,m);

    %% Guidance
    % Reference Values
    Pos_ref = [30; -30]; % X, Y position reference

    Alt_ref = -40;
    
    Yaw_ref_deg = 30;
    Yaw_ref = Yaw_ref_deg * Deg2Rad;
    
    %% Control for Position(Outer loop)
    Pos_err(:,m) = Pos_ref - Pos_sensor;

    if m > 1
        D_err_pos = (Pos_err(:, m) - Pos_err(:, m-1)) / dt; % 오차 미분
    else
        D_err_pos = [0; 0]; 
    end

    cmd_x = kp_pos * Pos_err(1,m) + kd_pos * D_err_pos(1);
    cmd_y = kp_pos * Pos_err(2,m) + kd_pos * D_err_pos(2);
    
    store_cmd_y(m) = cmd_y;
    cmd_pos = Cn2b * [cmd_x; cmd_y; 0];
    cmd_pos_pitch = - cmd_pos(1); % x
    cmd_pos_roll = cmd_pos(2); % y 

    store_cmd_pos_roll(m) = cmd_pos_roll; % y
    store_cmd_pos_pitch(m) = cmd_pos_pitch; % x
    %% Attitude_Control_Roll
    
    Rol_err(m) = cmd_pos_roll * Deg2Rad - EA_sensor(1);
    
    I_err_rol = I_err_rol + Rol_err(m) * dt; % 적분 오차 업데이트

    if m > 1
        D_err_rol = (Rol_err(m) - Rol_err(m-1)) / dt; % 오차 미분
    else
        D_err_rol = 0; 
    end

    cmd_rol = kp_rol * Rol_err(m) + ki_rol * I_err_rol + kd_rol * D_err_rol;

    %% Attitude_Control_Pitch

    Pit_err(m) = cmd_pos_pitch * Deg2Rad - EA_sensor(2);

    I_err_pit = I_err_pit + Pit_err(m) * dt; % 적분 오차 업데이트

    if m > 1
        D_err_pit = (Pit_err(m) - Pit_err(m-1)) / dt; % 오차 미분
    else
        D_err_pit = 0; 
    end

    cmd_pit = kp_pit * Pit_err(m) + ki_pit * I_err_pit + kd_pit * D_err_pit;

    %% Attitude_Control_Yaw
    Yaw_err(1,m) = Yaw_ref - Yaw_sensor;
    Yaw_err(m) = wrapTo180(Yaw_ref - Yaw_sensor);

    I_err_yaw = I_err_yaw + Yaw_err(m) * dt; % 적분 오차 업데이트

    if m > 1
        D_err_yaw = wrapTo180(Yaw_err(m) - Yaw_err(m-1)) / dt; % 오차 미분
    else
        D_err_yaw = 0; 
    end

    cmd_yaw = kp_yaw * Yaw_err(m) + ki_yaw * I_err_yaw + kd_yaw * D_err_yaw;

    %% Control for Altitude
    Alt_err(1,m) = Alt_ref - Alt_sensor;
    I_err_alt = I_err_alt + Alt_err(m) * dt; % 적분 오차 업데이트

    if m > 1
        D_err_alt = (Alt_err(m) - Alt_err(m-1)) / dt; % 오차 미분
    else
        D_err_alt = 0; 
    end

    cmd_alt = kp_alt * Alt_err(m) + ki_alt * I_err_alt + kd_alt * D_err_alt;

    %% Main Control
    
    u1 = cmd_alt + base_throttle; % F_z
    u2 = cmd_rol; % L
    u3 = cmd_pit; % M
    u4 = cmd_yaw; % N

    U = [u1, u2, u3, u4]';

    w(:, m+1) = sqrt(Rmat_inv * U);
    CONTROL(1:4, 1) = w(:,m+1); % AngVel.

    %% Numerical integration (dynamics, kinematics)
    STATE = [P(:,m); V(:,m); AngVel(:,m); Q(:,m)];
    
    DEL_X_1 = Dynamics_GyL_Ver03(STATE,                 CONTROL, Prop)*dt;
    DEL_X_2 = Dynamics_GyL_Ver03(STATE + (1/2)*DEL_X_1, CONTROL, Prop)*dt;
    DEL_X_3 = Dynamics_GyL_Ver03(STATE + (1/2)*DEL_X_2, CONTROL, Prop)*dt;
    DEL_X_4 = Dynamics_GyL_Ver03(STATE + DEL_X_3,       CONTROL, Prop)*dt;

    STATE_AFTER = STATE + (1/6)*(DEL_X_1 + 2*DEL_X_2 + 2*DEL_X_3 + DEL_X_4); % 4-th RK
    
    P(:,m + 1) = STATE_AFTER(1:3, 1);
    V(:,m + 1) = STATE_AFTER(4:6, 1);
    AngVel(:,m + 1) = STATE_AFTER(7:9, 1);
    Q(:,m + 1) = STATE_AFTER(10:13, 1);
    Q(:,m + 1) = Q(:,m + 1) / norm(Q(:,m + 1));
    
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
%% Result
Realtime = 0:dt:Tf;

LineWidth1 = 1.5;

figure;
subplot(3,1,1); grid on;
plot(Realtime, P(1,:), 'b', 'LineWidth',LineWidth1);
ylabel('X[m]')
subplot(3,1,2); grid on;
plot(Realtime, P(2,:), 'r', 'LineWidth',LineWidth1);
ylabel('Y[m]')
subplot(3,1,3); grid on;
plot(Realtime, P(3,:), 'r', 'LineWidth',LineWidth1);
ylabel('Z[m]')
xlabel('Time[sec]')

figure;
subplot(3,1,1); grid on;
plot(Realtime, EA(1,:) * Rad2Deg, 'b', 'LineWidth',LineWidth1);
ylabel('Roll(deg)')
subplot(3,1,2); grid on;
plot(Realtime, EA(2,:) * Rad2Deg, 'r', 'LineWidth',LineWidth1);
ylabel('Pitch[deg]')
subplot(3,1,3); grid on;
plot(Realtime, EA(3,:) * Rad2Deg, 'r', 'LineWidth',LineWidth1);
ylabel('Yaw[rad]')
xlabel('Time[sec]')

figure;
subplot(4,1,1); grid on;
plot(Realtime, w(1,:), 'b', 'LineWidth',LineWidth1);
ylabel('w1(rad/s)')
subplot(4,1,2); grid on;
plot(Realtime, w(2,:), 'r', 'LineWidth',LineWidth1);
ylabel('w2(rad/s)')
subplot(4,1,3); grid on;
plot(Realtime, w(3,:), 'r', 'LineWidth',LineWidth1);
ylabel('w3(rad/s)')
subplot(4,1,4); grid on;
plot(Realtime, w(4,:), 'r', 'LineWidth',LineWidth1);
ylabel('w4(rad/s)')
xlabel('Time[sec]')

figure;
plot(Realtime, Alt_err, 'b', 'LineWidth', LineWidth1);
ylabel('Alt err[m]')

figure;
plot(Realtime, Yaw_err * Rad2Deg, 'b', 'LineWidth', LineWidth1);
ylabel('Yaw err[m]')


