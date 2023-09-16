function [c,ceq] = dribbelConstraints(b)
global N k m mH g T gamma Q q4
bm = b;
b = zeros(2,k+1);
b(1,:) = modifyPolynomial(bm(1,1:k+1),k);
b(2,:) = modifyPolynomial(bm(1,k+2:end),k);
% dist = 0.1;
% y = -dist*cos(gamma);
% z = dist*sin(gamma);
q0 = b(:,end);
q0d = b(:,end-1);
jacc1 = jacobianLeft(q0);
ts = ones(k+1,1);
qT = b*ts;
qTd = [modifyPolynomial(polyder(b(1,:)),k); modifyPolynomial(polyder(b(2,:)),k)]*ts;
q = qT;
MPM = [(-m-2*m*cos(2*q(1))+4*(m+mH)*cos(4*q(1)))/(3*m+4*mH-2*m*cos(4*q(1))),                 (-2*m*cos(2*q(1)))/(3*m+4*mH-2*m*cos(4*q(1)));
       (8*(m+mH)*(1+2*cos(2*q(1)))*sin(q(1))^2)/(3*m+4*mH-2*m*cos(4*q(1))),                  (-m+2*m*cos(2*q(1)))/(3*m+4*mH-2*m*cos(4*q(1)))];
[qT,jacc] = relabelled_coords(qT);
jacc2 = jacobianLeft(qT);
c2 = forwardKinematicsLeft([0;0;q]);
ceq = [c2(2)-0.2;q0(2)+2*q0(1); q(2)+2*q(1);jacc1*q0d-jacc2*MPM*qTd];%q0d-jacc*(MPM*qTd)];
cz = zeros(N,1);
t = linspace(0,T,N);
for tm = 1:N
    ts = ones(k+1,1);
    for i = k:-1:0
        ts(k-i+1) = t(tm)^i;
    end
    q = b*ts;
    if tm == 1
        bd = [modifyPolynomial(polyder(b(1,:)),k); modifyPolynomial(polyder(b(2,:)),k)];
        qd = bd*ts;
        bdd = [modifyPolynomial(polyder(bd(1,:)),k); modifyPolynomial(polyder(bd(2,:)),k)];
        qdd = bdd*ts;
        c2dd = unconstrained_acceleration(q,qd,qdd);
    end
    c2 = forwardKinematicsLeft([0;0;q]);
    cz(tm) = c2(3);
end
c = [-cz;-c2dd];%; -sin(q(1))*qTd(1)+sin(qT(1)+qT(2))*(qTd(1)+qTd(2))+1.2]; %#[-sin(q0(1)+sin(q0(1)+q0(2))); 0];