% =========================================================================
% Title       : Code 1 – Diophantine Equation and Disturbance Modeling
% Author      : K.G.F. Barbosa
% Course      : Model Predictive Control (MPC) – 2025.1
% Date        : July 1, 2025
% License     : GNU General Public License v3.0 (GPL-3.0)
% Repository  : https://github.com/kgfbarbosa/CTRL-ModelPredictiveControl
% =========================================================================
%
% Description:
% This script solves a control design problem using a Diophantine equation
% approach. Specifically, it computes the controller polynomials Phi and f,
% as well as the disturbance model for a given linear system.
%
% Problem Setup:
%   - System parameters: a = 0.9, b = 0.2
%   - Prediction horizon: N = 5
%
% The code follows the structure of standard MPC formulation with
% disturbance modeling and is intended for academic purposes.
% -------------------------------------------------------------------------


% Initialize IDE
clc; clear; close all;

% --------------------------------------
% Model Parameters
% --------------------------------------
% System: y(k) = a * y(k-1) + b * u(k-1)
a = 0.9;        % A coefficient (auto-regressive)
b = 0.2;        % B coefficient (control)
N = 5;          % Prediction horizon
A = [1, -a];    % Polynomial A(z^-1) = 1 - a*z^-1
B = b;          % Polynomial B(z^-1) is scalar (gain only)


% --------------------------------------
% Diophantine equation: get F and G
% --------------------------------------
    % Diophantine equation solves:
    % 1 = A(z^-1) * F(z^-1) + z^-N * G(z^-1)
    % where:
    % A = A(z^-1) vector coefficients, ex: [1, -a];
    % N = Prediction horizon, ex: 5;
    
    % Initializes n_a auxiliar variable
na = length(A) - 1;
    % Calculates F
F = zeros(1, N+1);  % N+1 terms for F
F(1) = 1;
for k = 2:N+1
    F(k) = -A(2:end) * F(k-1:-1:max(k-na,1)).';
end
    % Calculate G(1) with the remaining error term
G = zeros(1, N);
G(1) = -A(2:end) * F(N:-1:max(N-na+1,1)).';


% --------------------------------------
% Build dynamic matrix Phi and free response f
% --------------------------------------
    % F: F(z^-1) coefficients vector
    % B: B(z^-1) gain, scalar
    % N: Prediction horizon
    
    % Builds Phi
Phi = zeros (N, N);
for i = 1:N
    for j = 1:i
        Phi(i, j) = F(i - j + 1) * B;
    end
end
    % f vector (free response)
f_vec = F(2:N+1).';     % Transpose for column vector


% --------------------------------------
% Disturbance model (assumed)
% --------------------------------------
C = A;  % assuming C(z^-1) = A(z^-1)


% --------------------------------------
% Results display
% --------------------------------------
disp("Phi matrix:");
disp(Phi);
disp("Free response vector f:");
disp(f_vec);
disp("Disturbance model C(z^-1):");
disp(C);
