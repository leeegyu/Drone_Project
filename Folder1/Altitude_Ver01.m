clc; clear all; close all;

m = 0.468; 
g = 9.80665;

K_p = 0.1; 
K_i = 0.01; 
K_d = 0.3;

% Transfer function of drone z-direction dynamics
s = tf('s');
G = 1 / (m * s^2);

% Controller T.F.
C = K_p + K_i / s + K_d * s;

% Feedback loop T.F.
T = feedback(C * G, 1);

% Time Response
t = 0:0.01:60; 
[y, t] = step(T, t); 

% Plot Altitude Response
figure;
plot(t, y, 'b', 'LineWidth', 2);
grid on;
xlabel('Time [s]');
ylabel('Altitude');
title('Altitude Response for PID Controller');
