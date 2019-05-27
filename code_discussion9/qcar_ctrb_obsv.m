% Simulate quarter-car suspension model.

%%%%%%%%%%%%%%%%%%%%%%%%
%   Tommaso Menara     %
%      UCR-ME121       %
%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

% Initialize parameters:
k1 = 16e3; % spring stiffness [N/m]
k2 = 16e4; % tire stiffness [N/m]
c1 = 1e3; % 
m = 45; % unsprung mass [kg]
M = 250; % sprung mass [kg]
g = 9.81; % gravity acceleration
W = M*g; % weight
F = 0.01*W; % active suspension force
                        
%% Linear state space model 

% system matrix
A = [0 1 0 0; 
    -k1/M -c1/M k1/M c1/M;
    0 0 0 1;
    k1/m c1/m -(k1+k2)/m -c1/m];

% input matrix
B = [0 F/m 0 -F/m;
    0 0 0 1e-2*k2/m]';

% output matrix
C = [1 0 0 0];

% feedforward matrix
D = [0 0];

qcar_input1 = ss(A,B(:,1),C,D(:,1)); % state space with input F
qcar_input2 = ss(A,B(:,2),C,D(:,2)); % state space with input road

Kdc = dcgain(qcar_input2); % compute ss value from input 2
Kr = 1/Kdc; % compute gain to scale input to follow step correctly
qcar_input2_scaled = ss(A,B(:,2)*Kr,C,D(:,2)); % state space with B matrix scaled

figure(1)
step(qcar_input2_scaled) % step response to second input

eig(A) % eigenvalues of A

%% controller design

CC = ctrb(A, B(:,1)); % controllability matrix
rank(CC)

p = [-10.4354+50i ...
    -10.4354-50i ...
    -1.6757+7.5142i ...
    -1.6757-7.5142i]; % desired poles
K = place(A, B(:,1), p) % gain matrix

qcar_CL = ss(A-B(:,1)*K, B(:,2), C, []) % close loop system

Kdc2 = dcgain(qcar_CL);
Kr2 = 1/Kdc2;
qcar_CL_scaled = ss(A-B(:,1)*K, B(:,2)*Kr2, C, []) % scaled-input close loop system

figure(1)
hold on
step(qcar_CL_scaled) % step response to first input

%% Luenberger state observer design

O = obsv(A,C); % observability matrix
rank(O)

q = [-11+50i ...
    -11-50i ...
    -7+1i ...
    -7-1i]; % desired observer poles

L = place(A', C', q).' % observer gain matrix
eig(A-L*C)

%% full-state feedback with an observer

At = [ A-B(:,1)*K             B(:,1)*K
       zeros(size(A))    A-L*C ];

Bt = [    B(:,2)*Kr2
       zeros(size(B,1),1) ];

Ct = [ C    zeros(size(C)) ];

x0 = rand(4,1);

sys_full = ss(At, Bt, Ct, 0); % state-space of full state feedback with observer

t = 0:0.01:5; % time vector for simulation
figure(2)
[y, t, x] = lsim(sys_full, zeros(size(t)), t, [x0; x0]); % simulate system

n = 4; % state dimension
x_tilde = x(:,n+1:end); % error
x = x(:,1:n); % state
x_hat = x - x_tilde; % state estimate

% Save state variables explicitly to aid in plotting
x1 = x(:,1); x2 = x(:,2); x3 = x(:,3); x4=x(:,4);
x1_hat = x_hat(:,1); x2_hat = x_hat(:,2); x3_hat = x_hat(:,3); x4_hat = x_hat(:,4);
plot(t,x1,'-r',t, x1_hat,':r',t, x2,'-b',t, x2_hat,':b',...,
    t,x3,'-g',t, x3_hat,':g',t, x4, '-m',t,x4_hat, ':m', 'Linewidth', 2)

%% LQR design

Q = eye(n); % do not penalize the state
% Q = diag([100 100 100 100]); % penalize state changes
R = 1; % do not penalize control energy
% R = .1; % more control energy allowed

[K_LQR, P, e] = lqr(qcar_input1, Q, R); % compute LQR gain

qcar_LQR = ss(A-B(:,1)*K_LQR, B(:,2), C, []) % feeback with LQR

Kdc3 = dcgain(qcar_LQR);
Kr3 = 1/Kdc3;
qcar_LQR_scaled = ss(A-B(:,1)*K_LQR, B(:,2)*Kr3, C, []) % scaled-input close loop system

figure(1)
step(qcar_LQR_scaled)
legend('sys', 'closed-loop sys', 'LQR')