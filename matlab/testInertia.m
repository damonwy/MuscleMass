l = 2;          % length of rigid body 1 and 2
m = 1;          % mass of rigid body 2
r = l/10;       % muscle offset by r from the proximal joint
theta1 = 0.22;
theta2 = 1;

%%%%%%%%%%%%%%%%%%%%%%%%% Inertia of Rigid Body 1
syms q1(t1, t2, s)
q0 = [0 0];
q1 = l * [cos(t1) sin(t1)];
qm = (1 - s) * q0 + s * q1;
J1 = jacobian(qm, [t1, t2]);
JTJ1 = J1' * J1;
I1(t1, t2) = m * int(JTJ1, s, 0, 1);
i1_ = I1(theta1, theta2);
i1 = m / 3 * [l^2 0; 0 0];
norm(double(i1_) - i1)

%%%%%%%%%%%%%%%%%%%%%%%%% Inertia of Rigid Body 2
syms p0(t1, t2, s)                                      % t1, t2 are angles of joint 1 and 2
p0 = l * [cos(t1) sin(t1)];                             % start point of rigid body 2
p1(t1, t2, s) = p0 + l * [cos(t1 + t2), sin(t1 + t2)];  % end point of rigid body 2
pm = (1 - s) * p0 + s * p1;                             % material point on rigid body 2
J2 = jacobian(pm, [t1, t2]);                            % differentiating pm wrt t1, t2
%simplify(J);

JTJ2 = J2' * J2;
I2 = m * int(JTJ2, s, 0, 1);                              % the inertia function of rigid body 2                             

% the inertia of rigid body 2 computed by symbolic func
i2_ = I2(theta1, theta2);

% my derivation 
i2 = m * l^2 / 3.0 * [4 + 3*cos(theta2) 1 + 1.5*cos(theta2); 1 + 1.5 * cos(theta2) 1];
norm(double(i2_)- i2)

% formula from paper
i22 = m/3.0 * [1 + 3 * l^2 + 3 * l * cos(theta2) 1 + 1.5 * l * cos(theta2); 1 + 1.5 * l *cos(theta2) 1];
norm(double(i2_) - i22)

%%%%%%%%%%%%%%%%%%%%%%%%% Inertia of Mucsle

wr(t1, t2, s) = p0 + r * [cos(t1 + t2), sin(t1 + t2)];  % muscle end point
wm = s * wr;                                            % muscle material point
Jm_ = jacobian(wm, [t1, t2]);                           % Jacobian 

syms Jm(t1, t2, s)
Jm = s * [-l * sin(t1) - r * sin(t1 + t2) -r * sin(t1 + t2); l * cos(t1) + r * cos(t1 + t2) r * cos(t1 + t2)]; % Jacobian written in the paper
% Jm_ = Jm checked

JTJm = Jm_' * Jm_;
Im(t1, t2) = m * int(JTJm, s, 0, 1);
im_ = Im(theta1, theta2);
im = m/3 * [l^2 + r^2 + 2*l*r*cos(theta2) r^2 + l*r*cos(theta2); r^2 + l*r*cos(theta2) r^2];
norm(double(im_) - im)
