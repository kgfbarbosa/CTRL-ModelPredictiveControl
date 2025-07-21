% =========================================================================
% Title       : Code 2.2 – Simulation of Example 4.11 from Julio's MPC Book
% Author      : K.G.F. Barbosa
% Course      : Model Predictive Control (MPC) – 2025.1
% Date        : July 18, 2025
% License     : GNU General Public License v3.0 (GPL-3.0)
% Repository  : https://github.com/kgfbarbosa/CTRL-ModelPredictiveControl
% =========================================================================
%
% Description:
% This script simulates Example 4.11 from the book:
%
%   "Controle Preditivo" by Julio E. Normey-Rico and Edson M. Hemerly
%
% The simulation is implemented in MATLAB and Simulink and aims to replicate
% the behavior of the predictive controller discussed in the example.
%
% Purpose:
%   - Apply MPC techniques to constrained systems
%   - Develop and simulate predictive controllers for academic understanding
%
% Requirements:
%   - MATLAB with Simulink
% -------------------------------------------------------------------------


% Initialize IDE
clc;        % Clear command window
clear;      % Clear all variables in the memory
close all;  % Close all open figures

% System
B = [0.3];                      % Polynomial A of the system (numerator)
A = [1, -0.9];                  % Polynomial B of the system (denominator)
C = [1, -0.4];                  % Polynomial C
d = 0;                          % System delay
Ts = 1;                         % Sample time
z = tf('z', Ts);                % Discrete z transform variable
Gz = tf(B, A, Ts);              % System transfer function
Gz = z^-d*Gz;                   % Add delay

% GPC parameters
N = 3;                          % Prediction horizon
lambda = 0.8;                   % Control weight variable
R = lambda * eye(N);            % Control weight matrix
numDelta = 1;                   % Integrator numerator
denDelta = [1 -1];              % Integrator denominator
Atilde = conv(A, denDelta);     % ATilde = A(z^-1) * delta operator
e = zeros(1, N);
F = [];
Gpartial = [];
impulse = [C, zeros(1, length(Atilde)-length(C))];
G_temp = zeros(N, N);   % Initialize temporary G matrix without delay
G = zeros(N, N);  % delayed G matrix

% Creates E vector, F matrix and Gpartial
for i = 1:2*N
    [q, r] = deconv(impulse, Atilde);  % Impulse/Atilde
    impulse = [r(2:end), 0];           % Impulse update for next step
    e(i) = q;                          % Store 1st quocient coef
    F = [F; r(2:end)]                 % F matrix
end

% G_temp = zeros(N, N);   % Initialize temporary G matrix without delay
for i = 1:N
    EB = conv(e, B);
    % Save the (n+1)-th element of the convolution (corresponds to y(k+n))
    if length(EB) >= i+1
        Gpartial = [Gpartial; EB(i+1)] % Auxiliar to create G
    else
        Gpartial = [Gpartial; 0]
    end

    EBbyC = deconv(EB, C);
    n_elements = min(length(EBbyC), N - i + 1);
    G_temp(i:i + n_elements - 1, i) = EBbyC(1:n_elements);
end
e = e(:,1:N);
F = F(1:N,:);

% G = zeros(N, N);  % delayed G matrix
if d < N
    G(d+1:end, :) = G_temp(1:N-d, :);
end


K0 = (G' * G + R) \ G';                % Cost function
K = K0(1, :);
Kr = sum(K);
KGp = K * Gpartial;
KF = K * F;

% Show results
disp('E = '); disp(e);
disp('F = '); disp(F);
disp('G = '); disp(G);
disp("(G'G+R)G' = "); disp(K0);

% Simulation
Tsim = 30 * Ts;
Tref = 0;
refInitialValue = 0;
refFinalValue = 1;
sim = sim('diophantineMPCSimulink');

% Plot
figure(1); hold on;
plot(sim.time, sim.y, 'b', LineWidth = 1.5);
plot(sim.time, sim.ref, '--r', LineWidth = 1);
plot(sim.time, sim.du, 'k', LineWidth = 1.5);
grid on;
xlabel('Discrete time');
title("Julio's MPC book - example 4.11");
legend('y(t)', 'u(t)', 'du(t)');
hold off;