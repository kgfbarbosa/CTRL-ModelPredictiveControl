% =========================================================================
% Title       : Code 5 – State-Space Simulation with Input Delay (d = 4)
% Author      : K.G.F. Barbosa
% Course      : Model Predictive Control (MPC) – 2025.1
% Date        : July 10, 2025
% License     : GNU General Public License v3.0 (GPL-3.0)
% Repository  : https://github.com/kgfbarbosa/CTRL-ModelPredictiveControl
% =========================================================================
%
% Description:
% This script simulates Example 4.4 from the book:
%
%   "Model Predictive Control" by Eduardo F. Camacho and Carlos Bordons
%
% The system is implemented using **state-space model** in MATLAB/Simulink.
% The key feature in this simulation is the presence of a fixed 
% **input delay**:
%   - Input delay d = 4 samples
%
% Purpose:
%   - Analyze the effect of fixed delays in predictive control design
%   - Test the performance of GPC under delayed input actuation
%
% Requirements:
%   - MATLAB with Simulink
% -------------------------------------------------------------------------


% Initialize IDE
clc;        % Clear command window
clear;      % Clear all variables in the memory
close all;  % Close all open figures


% Model
B = [0 0.2713];
A = [1 -0.8351];
C = [1 -0.8 0.12];
d = 4;                  % I/O delay
Ts = 0.1;               % Sampling time
z = tf('z', Ts);        % z transform operator
Gz = tf(B, A, Ts, 'variable', 'q^-1');  % Transfer function
Gz = z^-d*Gz;


% Controller parameters
Ny = 10;                % Prediction horizon
Nu = 10;                % Control horizon
lambda = 100;           % Control weight parameter
delta = 1;              % Output error tracking weight parameter
R = lambda*eye(Nu);     % Control weight matrix
Q = delta*eye(Ny);      % Output error tracking weight matrix


% Carima model
denDelta = [1 -1];
Delta = tf(1, denDelta, Ts, 'variable', 'z^-1');
Atilde = conv(A, denDelta);             % A integrated
    

% Adjust vectors length to be the same
maxLength = max([length(Atilde)-1, length(B)-1, length(C)-1]);
Atilde = [Atilde, zeros(1, maxLength-(length(Atilde)-1))];
A = [A, zeros(1, maxLength-(length(A)-1))];
B = [B, zeros(1, maxLength-(length(B)-1))];
C = [C, zeros(1, maxLength-(length(C)-1))];


% State-space formulation (Observable canonical form)
ssA = [-Atilde(2:end)' [eye(length(Atilde)-2); zeros(1,length(Atilde)-2)]];
ssB = [B(1:end-1)'];
ssD = [C(2:end)'-Atilde(2:end)'];
ssH = [1 zeros(1,maxLength-1)];


% Prediction Matrices

    % Creates G matrix
G = zeros(Ny,Nu);
for i=1:Ny
    G(i,1) = ssH*ssA^(i-1)*ssB;
end
for i=2:Ny
    for j=2:Nu
        G(i,j) = G(i-1,j-1);
    end
end

    % Creates F matrix
F = zeros(Ny, length(Atilde)-1);
for i=1:Ny
    F(i,:) = ssH*ssA^i;
end

    % Creates E matrix
E = zeros(Ny, 1);
for i=1:Ny
    E(i,1) = ssH*ssA^(i-1)*ssD;
end

    % Creates K gain vector
K = (G'*Q*G + R)\(G'*Q);
K = K(1,:);

    % Creates Kr reference gain
Kr = sum(K);

    % Creates KE and KF
if d == 0
    KE = K*E;
else
    KE = 0;
end
KF = K*F;


% Simulation
Tsim = 120*Ts;
Tref = 0;
refInitialValue = 0;
refFinalValue = 1;
Tdist = 0.4*Tsim;
distInitialValue = 0;
distFinalValue = -0.3;
Tnoise = 0.8*Tsim;
noiseInitialValue = 0;
noiseFinalValue = -0.3;
sim = sim('../MPCSimulink');

% Plot
figure
subplot(3,1,1)
stairs(sim.y,'b', linewidth=1.5)
hold on
stairs(sim.ref,'--r', linewidth=1.5)
xlabel('discrete time')
ylabel('y(k)')
axis("tight")
hold off
grid on
subplot(3,1,2)
stairs(sim.u,'r', linewidth=1)
xlabel('discrete time')
ylabel('u(k)')
axis("tight")
grid on
subplot(3,1,3)
stairs(sim.du,'r', linewidth=1)
xlabel('discrete time')
ylabel('Δu(k)')
axis("tight")
grid on