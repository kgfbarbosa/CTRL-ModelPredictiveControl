% =========================================================================
% Title       : Code 6.2 – Restricted State-Space Simulation of Camacho's
%               Book MIMO Example from page 139
% Author      : K.G.F. Barbosa
% Course      : Model Predictive Control (MPC) – 2025.1
% Date        : July 14, 2025
% License     : GNU General Public License v3.0 (GPL-3.0)
% Repository  : https://github.com/kgfbarbosa/CTRL-ModelPredictiveControl
% =========================================================================
%
% Description:
% This script simulates the **MIMO example** from page 139 of the book:
%
%   "Model Predictive Control" by Eduardo F. Camacho and Carlos Bordons
%
% The system is modeled in MATLAB using a **state-space representation**,
% and the predictive controller includes **constraints** on control actions
% and outputs (i.e., constrained MPC), using QuadProg.
%
% Key Features:
%   - Multi-Input Multi-Output (MIMO) system
%   - State-space modeling approach
%   - Constrained predictive control
%
% Purpose:
%   - Apply constrained MPC to a MIMO system
%   - Evaluate how constraints affect performance and feasibility
%
% Requirements:
%   - MATLAB
% -------------------------------------------------------------------------


% Initialize IDE
clc;        % Clear command window
clear;      % Clear all variables in the memory
close all;  % Close all open figures


% Model
Ts = 0.03;  % Sampling time
G1 = tf(1,[0.7 1]); G1d = c2d(G1, Ts, 'zoh');
G2 = tf(5,[0.3 1]); G2d = c2d(G2, Ts, 'zoh');
G3 = tf(1,[0.5 1]); G3d = c2d(G3, Ts, 'zoh');
G4 = tf(2,[0.4 1]); G4d = c2d(G4, Ts, 'zoh');


% Controller parameters
Ny = 3;                 % Prediction horizon
Nu = 2;                 % Control horizon
lambda = 0.05;          % Control weight parameter
delta = 1;              % Output error tracking weight parameter
R = lambda*eye(2*Nu);   % Control weight matrix
Q = delta*eye(2*Ny);    % Output error tracking weight matrix


% Carima model
denDelta = [1 -1];      % discrete differentiation polynomial
Delta = tf(1, denDelta, Ts, 'variable', 'z^-1');
                        % integrator


% Polynomials b(z^-1), a(z^-1) and c(z^-1) {disturbance rejection}
[b1, a1] = tfdata(G1d, 'v');
[b2, a2] = tfdata(G2d, 'v');
[b3, a3] = tfdata(G3d, 'v');
[b4, a4] = tfdata(G4d, 'v');
c1 = conv([1 0], [1 0]);
c2 = conv([1 0], [1 0]);

% Common denominator model
a1t = conv(denDelta, conv(a1,a2));  % delta*a1*a2
a2t = conv(denDelta, conv(a3,a4));  % delta*a3*a4
b1t = conv(b1,a2);                  % b1*a2
b2t = conv(b2,a1);                  % b2*a1
b3t = conv(b3,a4);                  % b3*a4
b4t = conv(b4,a3);                  % b4*a3

% Adjust vectors length to be the same
maxLength = max([length(a1t)-1, length(a2t)-1, length(b1t)-1, ...
    length(b3t)-1,length(b4t)-1, length(c1)-1, length(c2)-1]);


% Padding the vectors with zeros on the right
a1t = [a1t zeros(1, maxLength-(length(a1t)-1))];
a2t = [a2t zeros(1, maxLength-(length(a2t)-1))];
b1t = [b1t zeros(1, maxLength-(length(b1t)-1))];
b2t = [b2t zeros(1, maxLength-(length(b2t)-1))];
b3t = [b3t zeros(1, maxLength-(length(b3t)-1))];
b4t = [b4t zeros(1, maxLength-(length(b4t)-1))];
c1 = [c1 zeros(1, maxLength-(length(c1)-1))];
c2 = [c2 zeros(1, maxLength-(length(c2)-1))];


% State-space formulation (Observable canonical form)

    % Creates A coeficients (state polynomial)
A1 = [-a1t(2:end)' [eye(length(a1t)-2); zeros(1,length(a1t)-2)]];
A2 = [-a2t(2:end)' [eye(length(a2t)-2); zeros(1,length(a2t)-2)]];
    
    % Creates B coeficients (input polynomial)
B1 = [b1t(2:end)]';
B2 = [b2t(2:end)]';
B3 = [b3t(2:end)]';
B4 = [b4t(2:end)]';
    
    % Creates D coeficients (integrator polynomial)
D1 = [c1(2:end)'-a1t(2:end)'];
D2 = [c2(2:end)'-a2t(2:end)'];

    % Creates H coeficients (output polynomial))
H = [1 zeros(1,maxLength-1)];

    % Creates augmented matrices
A = [A1 zeros(size(A1)); zeros(size(A2)) A2];
B = [B1 B2; B3 B4];
D = [D1 zeros(size(D1)); zeros(size(D2)) D2];
H = [H zeros(size(H)); zeros(size(H)) H];


% Prediction matrices

    % Creates G matrix (prediction matrix)
G = zeros(2*Ny, 2*Nu);
for i=1:Ny
    for j=1:Nu
        if i>=j
            G(2*i-1:2*i, 2*j-1:2*j) = H*A^(i-j)*B;
        end
    end
end

    % Creates F matrix (initial state contribution)
F = zeros(Ny*size(H,1), size(A,2));
for i=1:Ny
    row_idx = (i-1)*size(H,1)+1 : i*size(H,1);
    F(row_idx, :) = H*A^i;
end

    % Creates E matrix (noise/disturbance contribution)
E = zeros(Ny*size(H,1), size(D,2));
for i=1:Ny
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
f = zeros(2*Ny,1);          % System free response vector
r = zeros(2*Ny,1);        % Reference vector
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
    for i=1:2*Ny
        f(i,1)=F(i,:)*x(:,k) + E(i,:)*e(:,k); % prediction of the system
        % behavior without applying control variation (du = 0).
    end
    
        % Reference update
    if k>(k_qtd/2)
        ref1Value = 0.4;
    end
    r = repmat([ref1Value; ref2Value], ceil(2*Ny/2), 1);
    r = r(1:2*Ny);
    ref(:,k) = r;

        % Linear cost vector calculation (Quadratic programming
        % equation)
    qpe = 2*(f-r)'*Q*G;


        % Constraints
    du_min = -0.2;
    du_max = 0.5;
    u_min = -0.1;
    u_max = 0.5;
    y_min = 0;
    y_max = 0.6;

        % Constraint matrices
    Acon = [ G;               % Output upper bound
            -G;               % Output lower bound
             eye(2, 2*Nu);    % u_max
            -eye(2, 2*Nu);    % u_min
             eye(2, 2*Nu);    % du_max
            -eye(2, 2*Nu)];   % du_min

    Bcon = [y_max - f;
            -y_min + f;
            (u_max - u0);
            -(u_min - u0);
            du_max * ones(2,1);
            -du_min * ones(2,1)];


        % QuadProg solution
    opt = optimoptions("quadprog","Display","off");
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

%
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