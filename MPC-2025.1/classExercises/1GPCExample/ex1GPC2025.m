clc; clear; close all;

A = [1, -0.8];
B = [0.4, 0.6];
C = 1;
d = 0;

na = length(A);
Ts = 1;
z = tf('z', Ts);
Delta = 1/(1 - z^-1);
[numDelta, denDelta] = tfdata(Delta, 'v');
ATilde = conv(A, denDelta);
F = []; Gp = []; G = [];
aux = [1, zeros(1,na)];
N = 3;
lambda = 0.8; R = 0.8*eye(N);
numGz = B;
denGz = A;
Gz = tf(numGz,denGz,Ts,'Variable','q^-1','IODelay',1);

for n = 1:N
    [q,r] = deconv(aux,ATilde);
    F = [F;r(2:end)];
    E(n) = q;
    aux = [r(2:end), 0];
    aux2 = conv(E,B); % cresce a cada iteracao mas inicialmente pode ser
    % bem curto. if garante que Gp sempre tenha N elementos, mesmo que
    % aux2(n+1) ainda nao exista.
    if length(aux2) >= n+1
        Gp = [Gp; aux2(n+1)];
    else
        Gp = [Gp; 0]; % Preenche com zero se aux2 nao tem n+1 elementos
    end
end
Ga = aux2(1:n)';
G = toeplitz(Ga, [Ga(1), zeros(1,N-1)]);

aux3 = (G'*G+R)\G';
K = aux3(1,:);
Kr = sum(K);
KGp = [0, K*Gp];
KF = K*F;
