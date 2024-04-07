%% Flat-Earth 6-Dof Simulation + 4th RK

clc; clear all; close all;
%% Initialization
Rad2Deg  = 180/pi;
Deg2Rad = pi/180;

dt = 0.1; % Time step
Tf = 100; % Total simulation time

g = [0 0 9.80665]';

k_F = 0.1; % [N/(rad^2/s^2)] Thrust coefficient
k_M = 7.8263*10^(-4); % [Nm/(rad^2/s^2)] Thrust coefficient
max_throttle = 10;
base_throttle = - max_throttle * 0.1291174078117;

l = 0.062; % [m] dist. between CoM and rotor
MASS  = 0.068; % [kg]
Ix = 6.85*10^(-5);
Iy = 9.2*10^(-5);
Iz = 1.366*10^(-4);
J = diag([Ix, Iy, Iz]);

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

I_err_alt = 0; % 적분 오차 초기화
ki_alt = 0.01;
kp_alt = 0.6;
kd_alt = 0.4;

D_err_alt = 0; 

dt_gc = 0.1; % Time step for guidance and control unit.
sampling_ratio = dt_gc/dt;

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
% EA(:,1) = EA0;

q0 = cos(EA0(1)/2)*cos(EA0(2)/2)*cos(EA0(3)/2) + sin(EA0(1)/2)*sin(EA0(2)/2)*sin(EA0(3)/2);
q1 = sin(EA0(1)/2)*cos(EA0(2)/2)*cos(EA0(3)/2) - cos(EA0(1)/2)*sin(EA0(2)/2)*sin(EA0(3)/2);
q2 = cos(EA0(1)/2)*sin(EA0(2)/2)*cos(EA0(3)/2) + sin(EA0(1)/2)*cos(EA0(2)/2)*sin(EA0(3)/2);
q3 = cos(EA0(1)/2)*cos(EA0(2)/2)*sin(EA0(3)/2) - sin(EA0(1)/2)*sin(EA0(2)/2)*cos(EA0(3)/2);

Q0 = [q0 q1 q2 q3]';
Q(:,1) = Q0;

AngVel0 = [0 0 0]*Deg2Rad';
AngVel(:,1) = AngVel0;

H_err = zeros(1, Tf/dt + 1);

%% Main Part(Iteration)
for m = 1 : Tf/dt
    % DCM from NED to body-fixed coordinate system : 3-2-1
    R1 = R1_GyL(EA(1, m));
    R2 = R2_GyL(EA(2, m));
    R3 = R3_GyL(EA(3, m));
    Cn2b = R1*R2*R3; % n2b = NED to body
    
    Cn2b = [Q(1,m)^2 + Q(2,m)^2 - Q(3,m)^2 - Q(4,m)^2  2*(Q(2,m)*Q(3,m) + Q(1,m)*Q(4,m))          2*(Q(2,m)*Q(4,m) - Q(1,m)*Q(3,m));
        2*(Q(2,m)*Q(3,m) - Q(1,m)*Q(4,m))              Q(1,m)^2 - Q(2,m)^2 + Q(3,m)^2 - Q(4,m)^2  2*(Q(3,m)*Q(4,m) + Q(1,m)*Q(2,m));
        2*(Q(2,m)*Q(4,m) + Q(1,m)*Q(3,m))              2*(Q(3,m)*Q(4,m) - Q(1,m)*Q(2,m))          Q(1,m)^2 - Q(2,m)^2 - Q(3,m)^2 + Q(4,m)^2];
    
    EA(1,m) = atan2(Cn2b(2,3), Cn2b(3,3));
    EA(2,m) = -asin(Cn2b(1,3));
    EA(3,m) = atan2(Cn2b(1,2), Cn2b(1,1));

    Prop.Cn2b = Cn2b;

    V_n(:,m) = (Cn2b)*V(:,m);

    % Currnet time
    t_curr = (m - 1)*dt;

    %% Sensor modeling (Navigaion)
    % Sampling ratio를 고려하지 않고 진행
    
    H_sensor = P(3, m);
    
    %% Guidance
    H_ref = -40;
    H_err(1,m) = H_ref - H_sensor;
    
    %% Control for Altitude
    I_err_alt = I_err_alt + H_err(1,m) * dt; % 적분 오차 업데이트

    if m > 1
        D_err_alt = (H_err(1,m) - H_err(1,m-1)) / dt; % 오차 미분
    else
        D_err_alt = 0; 
    end

    cmd_alt = kp_alt * H_err(1,m) + ki_alt * I_err_alt + kd_alt * D_err_alt;
    
    
    %% Control for roll & pitch
    


    %% Main Control
    
    u1 = base_throttle + cmd_alt; % 스로틀 조절
    u2 = 0; % L
    u3 = 0; % M
    u4 = 0; % N

    U = [u1, u2, u3, u4]';

    Rmat_inv = [-1/(4*k_F), -1/(2*sqrt(2)*l*k_F),  1/(2*sqrt(2)*l*k_F),  1/(4*k_M);
                -1/(4*k_F), -1/(2*sqrt(2)*l*k_F), -1/(2*sqrt(2)*l*k_F), -1/(4*k_M);
                -1/(4*k_F),  1/(2*sqrt(2)*l*k_F), -1/(2*sqrt(2)*l*k_F),  1/(4*k_M);
                -1/(4*k_F),  1/(2*sqrt(2)*l*k_F),  1/(2*sqrt(2)*l*k_F), -1/(4*k_M)];

    w(:, m) = sqrt(Rmat_inv * U);

    CONTROL(1:4, 1) = w(:,m); % AngVel.

    %% Numerical integration (dynamics, kinematics)
    STATE = [P(:,m); V(:,m); AngVel(:,m); Q(:,m)];
    
    DEL_X_1 = Dynamics_GyL2(STATE,                 CONTROL, Prop)*dt;
    DEL_X_2 = Dynamics_GyL2(STATE + (1/2)*DEL_X_1, CONTROL, Prop)*dt;
    DEL_X_3 = Dynamics_GyL2(STATE + (1/2)*DEL_X_2, CONTROL, Prop)*dt;
    DEL_X_4 = Dynamics_GyL2(STATE + DEL_X_3,       CONTROL, Prop)*dt;

    STATE_AFTER = STATE + (1/6)*(DEL_X_1 + 2*DEL_X_2 + 2*DEL_X_3 + DEL_X_4); % 4-th RK
    
    P(:,m + 1) = STATE_AFTER(1:3, 1);
    V(:,m + 1) = STATE_AFTER(4:6, 1);
    AngVel(:,m + 1) = STATE_AFTER(7:9, 1);
    Q(:,m + 1) = STATE_AFTER(10:13, 1);

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
plot(Realtime, H_err, 'b', 'LineWidth', LineWidth1);
ylabel('H_err[m]')
