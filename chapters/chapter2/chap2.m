%% Foundations
%% 2.1.1 Relative pose
%% Coordinate frames
%% Pose Graphs
%% Summary
%% Working in Two Dimensions (2D)
%% Orientation in Two Dimensions
%% 2D Rotation Matrix
R = rotm2d(0.3)
plottform2d(R);
det(R)
det(R * R)
syms theta real
R = rotm2d(theta)
simplify(R * R)
det(R)
simplify(ans)
%% Matrix Exponential for Rotation
R = rotm2d(0.3);
L = logm(R)
S = skew2vec(L)
X = vec2skew(2)
skew2vec(X)
expm(L)
expm(vec2skew(S))
%% Pose in Two Dimensions
%% 2D Homogeneous Transformation Matrix
rotm2d(0.3)
tformr2d(0.3)
TA = trvec2tform(1, 2) * tformr2d(30, 'deg')
axis([0 5 0 5]); % new plot with both axes from 0 to 5
plottform2d(TA, frame="A", color="b");
T0 = trvec2tform(0, 0);
plottform2d(T0, frame="0", color="k");  % reference frame
TB = trvec2tform(2, 1)
plottform2d(TB, frame="B", color="r");
TAB = TA * TB
plottform2d(TAB, frame="AB", color="g");
TBA = TB * TA;
plottform2d(TBA, frame="BA", color="c");
P = [3; 2];  % column vector
plotpoint(P, "ko", label="P");
inv(TA) * [P; 1]
h2e(inv(TA) * e2h(P))
%% Rotating a Coordinate Frame
axis([-5 4 -1 5]);
T0 = trvec2tform(0, 0);
plottform2d(T0, frame="0", color="k");
TX = trvec2tform(2, 3);
plottform2d(TX, frame="X", color="b");
TR = tformr2d(2);
plottform2d(TR * TX, framelabel="RX", color="g");
plottform2d(TX * TR, framelabel="XR", color="g");
C = [3; 2];
plotpoint(C, "ko", label="C");
TC = trvec2tform(C) * TR * trvec2tform(-C)
plottform2d(TC * TX, framelabel="XC", color="r");
%% Matrix exponential for Pose
L = logm(TC)
S = skewa2vec(L)'
expm(vec2skewa(S))
X = vec2skewa([1 2 3])
skewa2vec(X)
%% 2D Twists
S = Twist2d.UnitRevolute(C)
expm(vec2skewa(2 * S.S))
S.exp(2)
S.pole'
S = Twist2d.UnitPrismatic([0, 1])
S.exp(2)
T = trvec2tform(3, 4) * tformr2d(0.5)
S = Twist2d(T)
S.w
S.pole'
S.exp(1)
%% Working in Three Dimensions (3D)
%% Orientation in Three Dimensions
%% 3D Rotation Matrix
R = rotmx(pi / 2)
plottform(R);
animtform(R)
R = rotmx(pi / 2) * rotmy(pi / 2)
plottform(R);
rotmy(pi / 2) * rotmx(pi / 2)
%% Three-Angle Representations
R = rotmz(0.1) * rotmy(-0.2) * rotmz(0.3);
R = eul2rotm([0.1 -0.2 0.3], "ZYZ")
gamma = rotm2eul(R, "ZYZ")
R = eul2rotm([0.1 0.2 0.3], "ZYZ")
gamma = rotm2eul(R, "ZYZ")
eul2rotm(gamma, "ZYZ")
R = eul2rotm([0.1 0 0.3], "ZYZ")
rotm2eul(R, "ZYZ")
R = eul2rotm([0.1 0.2 0.3], "ZYX")
gamma = rotm2eul(R, "ZYX")
R = eul2rotm([0.1 0.2 0.3], "XYZ")
gamma = rotm2eul(R, "ZYX")
%tripleangle
%% Singularities and Gimbal Lock
%% Two-Vector Representation
a = [-1 0 0];
o = [1 1 0];
R = oa2rotm(o, a)
%% Rotation about an Arbitrary Vector
R = eul2rotm([0.1 0.2 0.3]);
rotm2axang(R)
[x,e] = eig(R)
theta = angle(e(2,2))
R = axang2rotm([1 0 0 0.3])
%% Matrix Exponential for Rotation
R = rotmx(0.3)
L = logm(R)
S = skew2vec(L)
expm(L)
expm(vec2skew(S))
R = rotmx(0.3);
R = expm(0.3 * vec2skew([1, 0, 0]));
X = vec2skew([1, 2, 3])
skew2vec(X)
%% Unit Quaternions
q = quaternion(rotmx(0.3), "rotmat", "point")
q = q * q;
q.conj
q * q.conj
q.rotmat("point")
q.compact
q.rotatepoint([0 1 0])
%% Pose in Three Dimensions
%% Homogeneous Transformation Matrix
T = trvec2tform(2, 0, 0) * tformrx(pi / 2) * trvec2tform(0, 1, 0)
plottform(T);
tform2rotm(T)
tform2trvec(T)'
%% Matrix exponential for Pose
T = trvec2tform(2, 3, 4) * tformrx(0.3)
L = logm(T)
S = skewa2vec(L)
expm(vec2skewa(S))
X = vec2skewa([1, 2, 3, 4, 5, 6])
skewa2vec(X)
%% 3D Twists
S = Twist.UnitRevolute([1 0 0], [0 0 0])
expm(0.3 * vec2skewa(S.S));
S.exp(0.3)
S = Twist.UnitRevolute([0, 0, 1], [2, 3, 2], 0.5);
X = trvec2tform(3, 4, -4);
hold on
for theta=[0:0.3:15]
  plottform(S.exp(theta) * X, style="rgb", width=2)
end
L = S.line();
L.plot('k:', linewidth=2);
S = Twist.UnitPrismatic([0 1 0])
S.exp(2)
T = trvec2tform(1, 2, 3) * eul2tform([0.3, 0.4, 0.5]);
S = Twist(T)
S.w'
S.pole'
S.pitch
S.theta
%% Vector-Quaternion Pair
%% Unit Dual Quaternion
%% Advanced Topics
%% Pre- and Post-Multiplication of Transforms
%% Rotating a Frame Versus Rotating a Point
%% Efficiency of Representation
%% Distance Between Orientations
q1 = quaternion(rotmx(pi / 2), "rotmat", "point");
q2 = quaternion(rotmz(-pi / 2), "rotmat", "point");
q1.dist(q2)
%% Normalization
R = eye(3, 3);
det(R) - 1
for i=1:100
  R = R * eul2rotm([0.2, 0.3, 0.4]);
end
det(R) - 1
R = tformnorm(R);
det(R) - 1
q = q.normalize();
%% Understanding the Exponential Mapping
%% More About Twists
S = Twist.UnitRevolute([1, 0, 0], [0, 0, 0])
S.S
S.v
S.w
S.se3()
expm(0.3 * S.se3())
S.exp(0.3)
S2 = S * S
S2.printline(orient="axang", unit="rad")
line = S.line()
line.plot("k:", linewidth=2);
T = trvec2tform(1, 2, 3) * eul2tform([0.3, 0.4, 0.5]);
S = Twist(T)
S / S.theta
S.unit();
S.exp(0)
S.exp(1)
S.exp(0.5)
%% Configuration Space
%% Using the Toolbox
R = rotmx(0.3)  % create SO(3) matrix as a MATLAB matrix
whos R
R = SO3.Rx(0.3)  % create SO3 object
whos R
R.double()
R = SO3(rotmx(0.3));  % convert SO(3) matrix
R = SO3.Rz(0.3);  % rotation about z-axis
R = SO3.RPY(10, 20, 30, unit="deg");  % from ZYX roll-pitch-yaw angles
R.rpy();  % convert to roll-pitch-yaw angles
R.eul();  % convert to ZYZ-Euler angles
R.printline();  % compact single-line print
TA = SE2(1, 2) * SE2(30, unit="deg")
whos TA
TA = SE2(1, 2, 30, unit="deg");
TA.R
TA.t
TA.plot(frame="A", color="b");
TA.printline()
P = [3; 2];
TA.inv() * P
R = SO3.Rx([0:0.1:0.5]);
whos R
R * [0; 1; 0]
%% Wrapping Up
%% Further Reading
%% Exercises
