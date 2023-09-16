function T = fkine_compass(q)
l = 1;
M = eye(4);
v1 = zeros(3,1);
w1 = [1;0;0];

q2 = [0;0;l];
w2 = [1;0;0];
v2 = cross(-w2,q2);

S1 = [[skew(w1); zeros(1,3)] [v1;0]];
S2 = [[skew(w2); zeros(1,3)] [v2;0]];

T = expm(S1*q(1))*(expm(S2*q(2))*M);