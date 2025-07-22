% =========================================================================
% Title       : Code 7 – Unconstrained MPC Simulation of Camachoos MIMO 
%               Case Study (Ex. 7.6)
% Author      : K.G.F. Barbosa
% Course      : Model Predictive Control (MPC) – 2025.1
% Date        : July 13, 2025
% License     : GNU General Public License v3.0 (GPL-3.0)
% Repository  : https://github.com/kgfbarbosa/CTRL-ModelPredictiveControl
% =========================================================================
%
% Description:
% This script simulates the **MIMO case study (Example 7.6)** from page 203
% of the book:
%
%   "Model Predictive Control" by Eduardo F. Camacho and Carlos Bordons
%
% The system is modeled using a **state-space formulation** in MATLAB. This
% version uses **unconstrained predictive control** with QuadProg, meaning
% no bounds are imposed on control inputs or system outputs.
%
% Key Features:
%   - Multi-Input Multi-Output (MIMO) dynamic system
%   - Example 7.6 from Camacho & Bordons (pg. 203)
%   - Unrestricted Generalized Predictive Control (GPC)
%
% Purpose:
%   - Explore GPC implementation for realistic multivariable systems
%
% Requirements:
%   - MATLAB
% -------------------------------------------------------------------------


% Initialize IDE
clc;        % Clear command window
clear;      % Clear all variables in the memory
close all;  % Close all open figures

% Process
Ts = 0.1;
z = tf('z',Ts);
G1 = tf(0.1133,[1.783 4.48 1]);
G2 = tf(0.9222,[2.071 1]);
G3 = tf(0.3378,[0.361 1.09 1]);
G4 = tf(-0.321,[0.104 2.463 1]);
GG = [G1 G2; G3 G4];
GGd = c2d(GG,Ts,'zoh');  % Discretization
GGd = [GGd(1,1)*z^(-14),GGd(1,2)*z^(-1);GGd(2,1)*z^(-6),GGd(2,2)*z^(-19)]

[b1,a1] = tfdata(GGd(1,1),'v');
[b2,a2] = tfdata(GGd(1,2),'v');
[b3,a3] = tfdata(GGd(2,1),'v');
[b4,a4] = tfdata(GGd(2,2),'v');
c1 = [1 0 0];
c2 = [1 0 0];

delta = [1 -1];
Delta = tf([1],delta,Ts,'Variable','z^-1','IODelay',0); % Delta
Delta = [Delta 0; 0 Delta];
a1t = conv(delta,conv(a1,a2));
a2t = conv(delta,conv(a3,a4));
b1t = conv(b1,a2);
b2t = conv(b2,a1);
b3t = conv(b3,a4);
b4t = conv(b4,a3);

% adjusting the vectors to have the same length
Max = max([length(a1t)-1,length(a2t)-1,length(b1t)-1,length(b2t)-1,length(b3t)-1,length(b4t)-1,length(c1)-1,length(c2)-1]);
a1t = [a1t zeros(1,Max-(length(a1t)-1))];
a2t = [a2t zeros(1,Max-(length(a2t)-1))];
b1t = [b1t zeros(1,Max-(length(b1t)-1))];
b2t = [b2t zeros(1,Max-(length(b2t)-1))];
b3t = [b3t zeros(1,Max-(length(b3t)-1))];
b4t = [b4t zeros(1,Max-(length(b4t)-1))];
c1 = [c1 zeros(1,Max-(length(c1)-1))];
c2 = [c2 zeros(1,Max-(length(c2)-1))];

N = 150;  % Prediction Horizon
Nu = 30;  % Control Horizon

% Ponderation matrices
Q = 1*eye(2*N);
R = 10*eye(2*Nu);

% State-Space formulation observable cannonical form
% Augmented Matrices
A1 = [-a1t(2:end)' [eye(length(a1t)-2); zeros(1,length(a1t)-2)]];
A2 = [-a2t(2:end)' [eye(length(a2t)-2); zeros(1,length(a2t)-2)]];
A = [A1 zeros(size(A1)); zeros(size(A2)) A2]
B1=[b1t(2:end)]';
B2=[b2t(2:end)]';
B3=[b3t(2:end)]';
B4=[b4t(2:end)]';
B=[B1 B2; B3 B4];
D1 = [c1(2:end)'-a1t(2:end)'];
D2 = [c2(2:end)'-a2t(2:end)'];
D = [D1 zeros(size(D1)); zeros(size(D2)) D2]
H = [1 zeros(1,Max-1)];
H = [H zeros(size(H)); zeros(size(H)) H]

% Prediction Matrices
G = zeros(2*N,2*Nu);
for i=1:N
    for j=1:Nu
        if i>=j
            G(2*i-1:2*i,2*j-1:2*j) = H*A^(i-j)*B;
        end
    end
end

F = zeros(N * size(H,1), size(A,2));
for i = 1:N
    row_idx = (i-1)*size(H,1)+1 : i*size(H,1);
    F(row_idx, :) = H * A^i;
end

disp('a1t');size(a1t)
disp('b1t');size(b1t)
disp('A1');size(A1)
disp('A');size(A)
disp('B');size(B)
disp('D');size(D)
disp('H');size(H)
disp('G');size(G)
disp('F');size(F)
disp('Q');size(Q)
disp('R');size(R)

    % Creates E matrix (noise/disturbance contribution)
E = zeros(N*size(H,1), size(D,2));
for i=1:N
    row_idx = (i-1)*size(H,1)+1 : i*size(H,1);
    E(row_idx, :) = H*A^(i-1)*D;
end

    % Creates GPC quadractic cost function matrix
HH = 2*(G'*Q*G + R);

% Simulation
    % Initial conditions
Tsim = 100*Ts;              % 100 samples
k_qtd = ceil(Tsim/Ts);      % Total simulation time
ref1Value = 0.5;
ref2Value = 0.3;
f = zeros(2*N,1);          % System free response vector
r = zeros(2*N,1);        % Reference vector
y0 = zeros(2,1);            % Initial output vector
du0 = zeros(2,1);           % Initial control variation vector
u0 = zeros(2,1);            % Initial control vector
x0 = zeros(length(A),1);    % Initial state vector
e0 = zeros(2,1);            % Initial error vector

    % Filters internal states (system startes at rest)
[~,z01] = filter(b1(2:end),a1,0);
[~,z02] = filter(b2(2:end),a2,0);
[~,z03] = filter(b3(2:end),a3,0);
[~,z04] = filter(b4(2:end),a4,0);

    % Simulation loop
for k=1:k_qtd
    
        % System update (controller internal model)
    y(:,k) = y0;
    x(:,k) = A*x0 + B*du0 + D*e0;
    e(:,k) = y(:,k) - H*x(:,k); % error between real and model outputs
    
        % Free response
    for i=1:2*N
        f(i,1)=F(i,:)*x(:,k) + E(i,:)*e(:,k); % prediction of the system
        % behavior without applying control variation (du = 0).
    end
    
        % Reference update
    if k>(k_qtd/2)
        ref1Value = 0.4;
    end
    r = repmat([ref1Value; ref2Value], ceil(2*N/2), 1);
    r = r(1:2*N);
    ref(:,k) = r;

        % Linear cost vector calculation (Quadratic programming
        % equation)
    qpe = 2*(f-r)'*Q*G;
    
        % QuadProg solution
    opt = optimoptions("quadprog","Display","off");
    Acon = []; Bcon = []; % unrestricted
    solution = quadprog(HH,qpe,Acon,Bcon,[],[],[],[],[],opt);

        % Control and real output signals update
    du(:,k) = solution(1:2);
    u(:,k) = u0 + du(:,k); % control variation
    
        % Simulates real plant
    [y1,z01] = filter(b1(2:end),a1,u(1,k),z01);
    [y2,z02] = filter(b2(2:end),a2,u(2,k),z02);
    [y3,z03] = filter(b3(2:end),a3,u(1,k),z03);
    [y4,z04] = filter(b4(2:end),a4,u(2,k),z04);
    
        
        % Outputs
    Y1 = y1+y2;
    Y2 = y3+y4;

        % State update for the next iteration
    du0 = du(:,k);
    u0 = u(:,k);
    x0 = x(:,k);
    e0 = e(:,k);
    y0 = [Y1; Y2];
end

% Plot
t = 1:1:k_qtd;

figure(1);

subplot(2,1,1);
hold on; grid on;
plot(t,ref(1,:),'--k','LineWidth',1);
plot(t,ref(2,:),'--k','LineWidth',1);
plot(t,y(1,:),'r','LineWidth',1);
plot(t,y(2,:),'b','LineWidth',1);
xlabel('Samples'); ylabel('y'); title('Output');

subplot(2,1,2);
hold on; grid on;
plot(t,u(1,:),'r','LineWidth',1);
plot(t,u(2,:),'b','LineWidth',1);
xlabel('Samples'); ylabel('u'); title('Control signal');