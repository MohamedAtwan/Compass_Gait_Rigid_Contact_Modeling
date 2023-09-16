function tau = inverseDribbelDynamics(q,qd,qdd)
global N k m mH g T gamma Q
M = zeros(2,2);
M(1,1) = mH + (3/2)*m-m*cos(q(2));
M(1,2) = (1/4)*m*(1-2*cos(q(2)));
M(2,1) = (1/4)*m*(1-2*cos(q(2)));
M(2,2) = (1/4)*m;

c11 = 0.0;
% c12 = 0.5*qd(1)*(-0.5*m*sin(q(2))-m*sin(q(2))-0) + 0.5*qd(2)*(-0.5*m*sin(q(2))-0.5*m*sin(q(2))-0);
% c21 = 0.5*qd(1)*(0 + 0 + m*sin(q(2))) + 0.5*qd(2)*(0 + 0 +0.5*m*sin(q(2)));
c12 = ((3/4)*qd(1)+qd(2)/2)*m*sin(q(2));
c21 = -(qd(1)/2 + qd(2)/4)*m*sin(q(2));
c22 = 0.0;

C = [c11 c12; c21 c22];

dvdq = zeros(2,1);
% dvdq(1) = -m*g*0.5*sin(q(1)-gamma) - mH*g*sin(q(1)-gamma) - m*g*sin(q(1)-gamma) + 0.5*m*g*sin(q(1)+q(2)-gamma);
dvdq(1) = -((3/2)*m*g+mH*g)*sin(q(1)-gamma)+(1/2)*m*g*sin(q(1)+q(2)-gamma);
dvdq(2) = m*g*0.5*sin(q(1)+q(2)-gamma);

tau = M*qdd+C*qd+dvdq;

