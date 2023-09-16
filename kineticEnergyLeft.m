function K = kineticEnergyLeft(q,qd)
global m mH 
A = [cos(sum(q))-cos(q(1)) cos(sum(q));
     sin(sum(q))-sin(q(1)) sin(sum(q))];
M = zeros(2,2);
M(1,1) = mH + (3/2)*m-m*cos(q(2));
M(1,2) = (1/4)*m*(1-2*cos(q(2)));
M(2,1) = (1/4)*m*(1-2*cos(q(2)));
M(2,2) = (1/4)*m;
M = inv(M);
D = A*inv(A'*M*A)*A';
K = -0.5*(qd'*D*qd);