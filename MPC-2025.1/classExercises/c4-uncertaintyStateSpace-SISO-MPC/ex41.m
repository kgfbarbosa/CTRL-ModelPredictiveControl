% =========================================================================
% Title       : Code 4 – Robust MPC Simulation of Camacho's Example 4.4
% Author      : K.G.F. Barbosa
% Course      : Model Predictive Control (MPC) – 2025.1
% Date        : July 8, 2025
% License     : GNU General Public License v3.0 (GPL-3.0)
% Repository  : https://github.com/kgfbarbosa/CTRL-ModelPredictiveControl
% =========================================================================
%
% Description:
% This script simulates Example 4.4 from the book:
%
%   "Model Predictive Control" by Eduardo F. Camacho and Carlos Bordons
%
% The system is modeled using a **state-space representation** with an
% **input delay**.
% Key characteristics of this simulation include:
%   - Input delay: d = 5
%   - Uncertain delay handling
%   - Robustness index analysis
%   - Parameter λ (lambda) adjusted to 110 to ensure GPC convergence
%
% Purpose:
%   - Explore the effect of input delay uncertainty on predictive control
%   - Apply robust GPC (Generalized Predictive Control) tuning strategies
%
% Requirements:
%   - MATLAB with Simulink
% -------------------------------------------------------------------------


% Initialize IDE
clc;        % Clear command window
clear;      % Clear all variables in the memory
close all;  % Close all open figures


% Delay uncertainty
% define uncertainty in percentage (50 = 50%)
delay_unc = 20;         % delay uncertainty value
% define uncertainty parameters
dU = (1+delay_unc/100); % delay uncertainty (1 + uncertainty percentage)
% disclaimer: the delay always must be an integer, because of the use of z 


% Model
B = [0.4 0.6];
A = [1 -0.8];
C = [1 -0.8 0.12];
d = 5*dU;               % I/O delay
Ts = 1;                 % Sampling time
z = tf('z', Ts);        % z transform operator
Gz = tf(B, A, Ts, 'variable', 'q^-1');  % Transfer function
Gz = z^-d*Gz;


% Controller parameters
Ny = 3;                 % Prediction horizon
Nu = 3;                 % Control horizon
lambda = 110;           % Control weight parameter
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


% Robustness index
% Multiplicative uncertainty
wz = logspace(-1,log10(pi/Ts),1000);
dP = abs(freqresp((Gz - Gz*z^-1),wz))./abs(freqresp(Gz*z^-1,wz));
dP = squeeze(dP(1,1,:));

aux = log10(pi/Ts);
n = 1000;
w = logspace(-2,aux,n);
Ir = zeros(1,n);
for k = 1:n
    zw = freqresp(z,w(k));
Gzw = zw^(-1)*(0.4+0.6*zw^-1)/(1 - 0.8*zw^-1);
V = K*E + (K*F-K*E*ssH)*((zw*eye(size(ssA))-ssA+ssD*ssH )\ssD);
ir = (Gzw*V\((K*F)*((zw*eye(size(ssA))-ssA)\ssB)+1))*(1-zw^(-1));
Ir(:,k) = abs(ir);
end


% Simulation
Tsim = 100*Ts;         % 100 samples
Tref = 0;
refInitialValue = 0;
refFinalValue = 1;
Tdist = 0.5*Tsim;
distInitialValue = 0;
distFinalValue = -0.5;
Tnoise = 0.8*Tsim;
noiseInitialValue = 0;
noiseFinalValue = -0.5;
sim = sim('../MPCSimulink');

% Plot
figure(1)
subplot(3,1,1)
stairs(sim.y,'b', linewidth=1.5)
xlabel('discrete time')
ylabel('y(k)')
axis("tight")
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

figure(2)
loglog(w, dP, '--k', linewidth=1)
hold on
loglog(w, Ir, 'b', linewidth=1.5)
% grid on
% axis ([1e-2 1e1 1e-2 1e3])
xlabel('frequency ω')
legend('Uncertainty δP(ω)','Robustness index Ir(ω)')