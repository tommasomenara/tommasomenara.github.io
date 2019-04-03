% Spring-Mass-Damper modeling

clear all
close all
clc

m = 1; % mass [kg]
k = 1; % spring constant [N/m]
b = 0.2; % damping constant [Ns/m]  
F = 1; % input force [N]

A = [0 1; -k/m -b/m]; % system matrix
B = [0 1/m]'; % input matrix
C = [1 0]; % output matrix
D = [0]; % feedforward matrix

%% Equivalent representations of the system in time and frequency domains

s = tf('s');
sys1 = 1/(m*s^2+b*s+k)

sys1 = ss(A,B,C,D) % state-space representation of the system

[b a] = ss2tf(A,B,C,D) % state-space to transfer function representation

%% Impulse, unit and step response

figure(1);
step(sys1);

figure(2);
impulse(sys1);

%% change some parameters, e.g. damping

m = 1; % mass [kg]
k = 1; % spring constant [N/m]
b = 0.4; % damping constant [Ns/m]  (try also 0 and 2 to see undamped and overdamped cases)
F = 1; % input force [N]

A = [0 1; -k/m -b/m]; % system matrix
B = [0 1/m]'; % input matrix
C = [1 0]; % output matrix
D = [0]; % feedforward matrix

sys2 = ss(A,B,C,D) % state-space representation of the system

figure(3)
step(sys1,'b',sys2,'r')
legend('b = 0.2','b = 0.4')

%% integrate system

t = 0:0.1:60;
u = F*ones(length(t),1);

figure(4)
lsim(sys1,u,t)