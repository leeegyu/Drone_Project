clc; clear all; close all;

m = 0.468; 

K_p = 0.1; 
K_i = 0.01; 
K_d = 0.3;

% transfer function of drone z-direction dynamics
s = tf('s');
G = 1 / (m * s^2);

% Controller T.F.
C = K_p + K_i / s + K_d * s;

% Feedback loop T.F.
T = feedback(C * G, 1);
T
poles = pole(T);
zeros = zero(T);

figure;
plot(real(poles), imag(poles), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
plot(real(zeros), imag(zeros), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
xlabel('Real Part');
ylabel('Imaginary Part');
title('Pole, Zero of the Closed-Loop Sys');
legend('Pole', 'Zero');
xlim([-2 2]);
ylim([-2 2]);

% Time Rsponse
t = 0:0.01:60; 
[y, t] = step(T, t); 

figure;
plot(t, y, 'b', 'LineWidth', 2);
grid on;
xlabel('Time');
ylabel('Altitude');
title('Altitude Response for PID Controller');
