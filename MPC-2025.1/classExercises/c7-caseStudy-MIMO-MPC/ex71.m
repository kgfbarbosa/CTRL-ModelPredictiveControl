% =========================================================================
% Title       : Code 7 – Unconstrained MPC Simulation of Camacho's MIMO 
%               Case Study (Ex. 7.6)
% Author      : K.G.F. Barbosa
% Course      : Model Predictive Control (MPC) – 2025.1
% Date        : July 13, 2025 (edited in August 6, 2025 to "final version")
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
% For me is the final code using Quadprog to SISO and MIMO delayed systems
% even when added uncertainties to it.
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
%
% Parameters used by me to reproduce the example:
%      -- controller parameters --
%   Ny = 30;                    % Prediction horizon
%   Nu = 3;                     % Control horizon
%   lambda = 1;                 % Control weight parameter
%   delta = 0.8;                % Output error tracking weight parameter
%      -- polynomial c for disturbance rejection --
%   c{i} = conv([1 0], [1 0]);  % Example: double integrator
%      -- simulation time --
%   Tsim = 60;
%      -- reference setup --
%   ref1 is 1 from 5 to 35 seconds
%   ref2 is 0.7 from 15 to 45 seconds
%       only to try to match with the book's graphs
%      -- restrictions setup --
%   u_min = -2.75;
%   u_max = 2.75;
%   y_min = 0;
%   y_max = 1;
% -------------------------------------------------------------------------


% -- Initialize IDE --
clc;                     % Clear command window
clear;                   % Clear all variables in the memory
close all;               % Close all open figures

% -- Define uncertainties --
gain_unc = 0;            % gain uncertain value
delay_unc = 0;           % delay uncertain value
sample_unc = 0;          % sample uncertain value

% -- Define uncertainty parameters --
gU = (1+gain_unc/100);   % gain uncertainty (1 + uncertainty percentage)
dU = (1+delay_unc/100);  % delay uncertainty (1 + uncertainty percentage)
sU = (1+sample_unc/100); % sample uncertainty (1 + uncertainty percentage)

% -- Input and output numbers; sample time value and 's' and 'z'
% parameters --
n_in = 2; n_out = 2;
s = tf('s'); z = tf('z');
Ts = 0.05*sU;


% System
    % -- System delay --
dC = [0.715, 0.0; 0.299, 0.94]*dU;
dD = round(dC*(1/Ts));
    % -- System TF matrix --
G = [ ...
    (0.1133*gU)/(1.783*s^2+4.48*s+1), (0.9222*gU)/(2.071*s+1);
    (0.3378*gU)/(0.361*s^2+1.09*s+1), (-0.321*gU)/(0.104*s^2+2.463*s+1);];

    % -- System discretization (ZOH) --
Gd = tf(zeros(n_out,n_in));
for i = 1:n_out
    for j = 1:n_in
        Gd(i,j) = c2d(G(i,j), Ts, 'zoh');
    end
end
    % -- System delay implementation --
GGd = tf(zeros(n_out,n_in));
for i = 1:n_out
    for j = 1:n_in
        GGd(i,j) = Gd(i,j)*z^(-dD(i,j));
    end
end

% Controller
    % -- Controller parameters --
Ny = 30;                     % Prediction horizon
Nu = 3;                     % Control horizon
lambda = 1;                  % Control weight parameter
delta = 0.8;                   % Output error tracking weight parameter
% In GPC, to “see” the effect of control in the future, the prediction
% horizon Ny must be greater than the system's largest delay, because if 
% the horizon ends before the delay, the controller cannot anticipate the
% impact of the current control action on the output.
% 
% Therefore:
% For GPC to work well, it's necessary that Ny ≥ the largest delay
% in samples.
% 
% The control horizon Nu defines for how many future steps the controller
% can adjust the control actions. It [doesn't need to be as large as ' ...
% Ny — in practice:]
% 
% Small Nu (e.g., 10–20) → smoother control.
% Large Nu (e.g., ≥ Ny) → more aggressive control (because the controller
% tries to correct everything at once).

R = lambda * eye(n_in * Nu); % Control weight matrix
Q = delta * eye(n_out * Ny); % Output error tracking weight matrix

% Carima model
denDelta = [1 -1];           % Discrete differentiation polynomial
Delta = tf(1, denDelta, Ts, 'variable', 'z^-1');
                             % Integrator


% Polynomials
    % -- b(z^-1), a(z^-1) and c(z^-1) {disturbance rejection} --
    % -- Extracting coeficients --
for i = 1:n_out
    for j = 1:n_in
        [b{i,j}, a{i,j}] = tfdata(GGd(i,j), 'v');
    end
end
    % -- Initialize vectors --
a_t = cell(n_out,1);         % Output common denominators
b_t = cell(n_out, n_in);     % Adjusted numerators
c = cell(n_out,1);           % Disturbance rejection vector
    % -- Create common denominator output --
for i = 1:n_out
    a_temp = denDelta;
    for j = 1:n_in
        a_temp = conv(a_temp, a{i,j});
    end
    a_t{i} = a_temp;
end
    % -- Create adjusted b for every output-input --
for i = 1:n_out
    for j = 1:n_in
        % Another channels denominator from the same output (except j)
        denom_others = 1;
        for k = 1:n_in
            if k ~= j
                denom_others = conv(denom_others, a{i,k});
            end
        end
        b_t{i,j} = conv(b{i,j}, denom_others); % b_ij*another denominators
    end
end

    % -- Create polynomial c for disturbance rejection --
for i = 1:n_out
    c{i} = conv([1 0], [1 0]);  % Example: double integrator
    % c{i} = [1 -0.98]; % Integrator
end

    % -- Find max Length --
maxLength = 0;
for i = 1:n_out
    maxLength = max(maxLength, length(a_t{i}) - 1);
    maxLength = max(maxLength, length(c{i}) - 1);
    for j = 1:n_in
        maxLength = max(maxLength, length(b_t{i,j}) - 1);
    end
end
    % -- Padding --
for i = 1:n_out
        % Polynomial a_t's padding (one for output)
    a_t{i} = [a_t{i}, zeros(1, maxLength - (length(a_t{i}) - 1))];
        % Polynomial c's padding (one for output)
    c{i} = [c{i}, zeros(1, maxLength - (length(c{i}) - 1))];
        % Polynomial b_t's padding (one for input-output)
    for j = 1:n_in
        b_t{i,j} = [b_t{i,j}, zeros(1, maxLength-(length(b_t{i,j}) - 1))];
    end
end


% State-space formulation (Observable canonical form)
    % -- Initialize matrices --
A = []; B = []; D = []; H = [];

for i = 1:n_out
    % -- Create matrix A --
    N = maxLength;   % System's total order
    a_vec = a_t{i};  % Polynomial a_i(z) = 1 + a1*z^-1 + a2*z^-2 + ...
    c_vec = c{i};
    Ai = zeros(N);                    % Initialize Ai with zeros
    for row = 1:N
        Ai(row,1) = -a_vec(row+1);    % First column: -a2, -a3, ...
        if row < N
            Ai(row, row+1) = 1;       % Displaced identity
        end
    end
    
    % -- Create matrix D --
    Di = c_vec(2:end) - a_vec(2:end);
    if length(Di) < maxLength
        Di = [Di, zeros(1, maxLength - length(Di))];
    end
    Di = Di';
    
    % -- Create matrix H --
    Hi = [1, zeros(1, maxLength-1)];
    
    % -- Create matrix B --
    Bi_row = [];
    for j = 1:n_in
        b_vec = b_t{i,j};
        b_rest = b_vec(2:end);
        if length(b_rest) < maxLength
            b_rest = [b_rest, zeros(1, maxLength - length(b_rest))];
        end
        Bi_row = [Bi_row, b_rest'];
    end

    % -- Accumulates in the final matrices --
    A = blkdiag(A, Ai);
    D = blkdiag(D, Di);
    H = blkdiag(H, Hi);
    B = [B; Bi_row];
end


% Prediction matrices
% -- Create matrix G (impact of future controls on future outputs) --
G = zeros(n_out*Ny, n_in*Nu);
for i = 1:Ny
    for j = 1:Nu
        if i >= j
            G((i-1)*n_out+1:i*n_out, (j-1)*n_in+1:j*n_in) = H * A^(i-j)* B;
        end
    end
end

% -- Create matrix F (impact of initial states on future outputs) --
F = zeros(n_out*Ny, size(A,2));
for i = 1:Ny
    F((i-1)*n_out+1:i*n_out, :) = H * A^i;
end

% -- Create matrix E (impact of disturbance on future outputs) --
E = zeros(n_out*Ny, size(D,2));
for i = 1:Ny
    E((i-1)*n_out+1:i*n_out, :) = H * A^(i-1) * D;
end

% -- Quadratic cost function --
HH = 2*(G'*Q*G + R);

% Simulation
Tsim = 60;                    % Number of samples

k_qtd = ceil(Tsim/Ts);        % Total simulation time
f = zeros(n_out*Ny,1);        % System free response vector
r = zeros(n_out*Ny,1);        % Reference vector
y0 = zeros(n_out,1);          % Initial output vector
du0 = zeros(n_in,1);          % Initial control variation vector
u0 = zeros(n_in,1);           % Initial control vector
x0 = zeros(length(A),1);      % Initial state vector
e0 = zeros(n_out,1);          % Initial error vectors

% -- Filters internal states (system startes at rest) --
z = cell(n_out, n_in);
for i = 1:n_out
    for j = 1:n_in
        [~,z{i,j}] = filter(b{i,j}(2:end), a{i,j},0);
    end
end

for k=1:k_qtd
    y(:,k) = y0;
    x(:,k) = A*x0 + B*du0 + D*e0;
    e(:,k) = y(:,k) - H*x(:,k);
    
    % -- Free response --
    for i=1:n_out*Ny
       f(i,1) = F(i,:)*x(:,k) + E(i,:)*e(:,k); 
    end
    % -- Future reference --
    if k > 35/Ts
        r1 = 0;
    elseif k < 5/Ts
        r1 = 0;                     % Reference 1
    else
        r1 = 1;
    end
    if k > 45/Ts
        r2 = 0;
    elseif k < 15/Ts
        r2 = 0;
    else
        r2 = 0.7;                   % Reference 2
    end
    r = repmat([r1; r2], Ny, 1);    % Repeat to Ny steps
    ref(:,k) = r;
    e_real(:,k) = ref(1:n_out,k) - y(:,k); % Error used in ISE
    
    % -- Quadratic equation --
    qpe = 2*(f-r)'*Q*G;
    % -- QP solver --
    opt = optimoptions('quadprog','Display','off');
    
    
    % Constraints
    % -- Constraints values --
    % du_min = -1; % Non disturbance restrictions,
    % du_max = 1;  % comment disturbance add below
    % u_min = -0.2;
    % u_max = 0.3;
    % y_min = -0.1;
    % y_max = 1.1;
    du_min = -inf; % Disturbance restrictions,
    du_max = inf;  % comment non disturbance add below
    u_min = -2.75;
    u_max = 2.75;
    y_min = 0;
    y_max = 1;
    % -- Constraint matrices --
    Acon = [ G;                      % Output upper bound
            -G;                      % Output lower bound
             eye(n_in, n_in*Nu);     % u_max
            -eye(n_in, n_in*Nu);     % u_min
             eye(n_in, n_in*Nu);     % du_max
            -eye(n_in, n_in*Nu) ];   % du_min
    Bcon = [y_max - f;
            -y_min + f;
            (u_max - u0);
            -(u_min - u0);
            du_max * ones(n_in,1);
            -du_min * ones(n_in,1)];

    [sol,~,exitflag] = quadprog(HH,qpe,Acon,Bcon,[],[],[],[],[],opt);
    exitflag;
    if isempty(sol)
        warning("Quadprog didn't return a solution at iteration %d", k);
        break;
    end
    if exitflag <= 0
        warning("Quadprog failed at iteration %d with exitflag = %d",...
            k, exitflag);
        break;
    end

    du(:,k) = sol(1:n_in);                  % Control signal
    u(:,k) = u0 + du(:,k);                  % Control signal variation
    
    % Create disturbance vectors
    d1 = 0.3;                 % Constant
    d2 = 0;
    if k*Ts >= 15             % When time >= 15 s
        d2 = 0.3;
    end
    % u_with_dist = u(:,k) + [d1; d2]; % Disturbance add

    u_with_dist = u(:,k);            % Non disturbance add
                                       % both in control signal

    Y = zeros(n_out,1);
    for i_out = 1:n_out
        Y(i_out) = 0;
        for j_in = 1:n_in
            [y_tmp, z{i_out,j_in}] = filter(b{i_out,j_in}(2:end), ...
                a{i_out,j_in}, u_with_dist(j_in), z{i_out,j_in});
            Y(i_out) = Y(i_out) + y_tmp;
        end
    end
    du0 = du(:,k);
    u0 = u(:,k);
    x0 = x(:,k);
    e0 = e(:,k);
    y0 = Y;
end


% Plot
% t = 1:1:k_qtd;
k_final = k;          % Last valid iteration number
% t = 1:k_final;
t = (0:k_final-1)*Ts; % Show time instead of samples
figure(1);

% Performance Index calculation
% ISE - Integral of Squared Error index
ISE_1 = sum(e_real(1,1:k_final).^2);
ISE_2 = sum(e_real(2,1:k_final).^2);
% TV - Total Variation index
TV_1 = sum(abs(diff(u(1,1:k_final))));
TV_2 = sum(abs(diff(u(2,1:k_final))));
% Show results
fprintf('ISE1 = %.4f | ISE2 = %.4f\n', ISE_1, ISE_2);
fprintf('TV1 = %.4f | TV2 = %.4f\n', TV_1, TV_2);

subplot(2,2,1);
hold on; grid on;
plot(t, ref(1,1:k_final), '--r', 'LineWidth', 1);
% plot(t, ref(2,1:k_final), '--k', 'LineWidth', 1);
% plot(t, ref(3,1:k_final), '--k', 'LineWidth', 1);
plot(t,y(1,:),'b','LineWidth',1);
% plot(t,y(2,:),'b','LineWidth',1);
xlabel('Samples'); ylabel('y1'); title('Output');
legend('Ref', 'GPC');

subplot(2,2,3);
hold on; grid on;
plot(t,u(1,:),'b','LineWidth',1);
% plot(t,u(2,:),'b','LineWidth',1);
% plot(t,u(3,:),'g','LineWidth',1);
xlabel('Samples'); ylabel('u1'); title('Control signal');

subplot(2,2,2);
hold on; grid on;
% plot(t, ref(1,1:k_final), '--k', 'LineWidth', 1);
plot(t, ref(2,1:k_final), '--r', 'LineWidth', 1);
% plot(t, ref(3,1:k_final), '--k', 'LineWidth', 1);
% plot(t,y(1,:),'r','LineWidth',1);
plot(t,y(2,:),'b','LineWidth',1);
xlabel('Samples'); ylabel('y2'); title('Output');
% legend('Ref', 'GPC');

subplot(2,2,4);
hold on; grid on;
% plot(t,u(1,:),'r','LineWidth',1);
plot(t,u(2,:),'b','LineWidth',1);
% plot(t,u(3,:),'g','LineWidth',1);
xlabel('Samples'); ylabel('u2'); title('Control signal');