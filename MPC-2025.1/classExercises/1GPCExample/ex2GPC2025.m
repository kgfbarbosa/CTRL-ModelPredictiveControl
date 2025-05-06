clc; clear; close all;

A = [1, -0.9];
B = [0.3];      % B deve ser vetor
C = [1 -0.4];
d = 1;          % delay
Ts = 1;         % sampling time

na = length(A);
z = tf('z', Ts);
Delta = 1/(1 - z^-1);
[numDelta, denDelta] = tfdata(Delta, 'v');
ATilde = conv(A, denDelta);

F = []; Gp = []; G = [];
aux = [1, zeros(1,na)];
N = 3;
lambda = 0.8;
R = lambda*eye(N);

numGz = B;
denGz = A;
Gz = tf(numGz, denGz, Ts, 'Variable', 'q^-1', 'IODelay', d);

E = zeros(1, N); % pre-alocacao para E
for n = 1:N
    [q,r] = deconv(aux, ATilde);
    F = [F; r(2:end)];
    E(n) = q;
    aux = [r(2:end), 0];
    aux2 = conv(E,B); % E vetor com coeficientes escalar
    if length(aux2) >= n+1
        Gp = [Gp; aux2(n+1)];
    else
        Gp = [Gp; 0];
    end
end
Ga = aux2(1:n)';
G = toeplitz(Ga, [Ga(1), zeros(1,N-1)]);

aux3 = (G'*G+R)\G';
K = aux3(1,:);
Kr = sum(K);
KGp = [0, K*Gp];
KF = K*F;
