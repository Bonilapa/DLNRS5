clear all
syms fi s dfi ds ddfi dds L1 L2 L d2 m1 m2 I1 I2 g real
d2 = L2/2;
h = 10;
%% Kinetic energy
% Jacobians found in classic approach on paper
Jv1 = [
    0 0
    0 0
    0 0];
Jv2 = [
    -(s+d2)*sin(fi) cos(fi)
    (s+d2)*cos(fi) sin(fi) 
    0 0];
Jw1 = [
    0 0
    0 0
    1 0];
Jw2 = [
    0 0
    0 0
    1 0];
R1 = [
cos(fi) -sin(fi) 0
sin(fi) cos(fi) 0
0 0 1];
R2 = R1;

q = [fi; s];
D1 = m1 * Jv1' * Jv1 + Jw1' * R1 * I1 * R1' * Jw1;
D2 = m2 * Jv2' * Jv2 + Jw2' * R2 * I2 * R2' * Jw2;

%% Potential energy
D = D1 + D2;
D = simplify(D);

% assume the first link is under the surface on a height h
P1 = m1 * g * h;
P2 = m2 * g * (h + (s+d2)*sin(fi));
P = P1 + P2;

G1 = diff(P, fi);
G2 = diff(P, s);
G = [G1; G2];

%% Coriolis force
dq = [dfi; ds];
ddq = [ddfi; dds];
C = Coriolis(D, q, dq, 2);
C = simplify(C);

%% Final equation
tor = D * ddq + C * dq + G;
D(fi, s) = subs(D, {m1, m2, I1, I2, L2}, {2 2 1 2 1});
C(fi, s, dfi, ds) = subs(C*dq, {m1, m2, I1, I2, L2}, {2 2 1 2 1});
G(fi, s) = subs(G, {m1 m2 I1 I2 L2 g}, {2 2 1 2 1 9.81});
D(0, 0);
C(0, 0, 0, 0);
G(0, 0);

fi_0 = pi/2;
s_0 = 0;
dfi_0 = 0;
ds_0 = 0;
ddfi_0 = 0;
dds_0 = 0;
dt = 0.01;
U = [0; 5.81*2]; % stable vertical position with fi = pi/2
n = 300;

for i = 1:n
    fip(i) = fi_0;
    sp(i) = s_0;
    dfip(i) = dfi_0;
    dsp(i) = ds_0;
    ddfip(i) = ddfi_0;
    ddsp(i) = dds_0;
    ddq = inv(D(fi_0, s_0)) * (U-C(fi_0, s_0, dfi_0, ds_0) - G(fi_0, s_0));
    ddfi_0 = ddq(1);
    dds_0 = ddq(2);
    dfi_0 = dfip(i) + double(ddq(1) * dt);
    ds_0 = dsp(i) + double(ddq(2) * dt);
    fi_0 = fip(i) + dfi_0 * dt;
    s_0 = sp(i) + ds_0*dt;
end
t = 0:0.1:(0.1*(n-1));
figure()
plot(t, fip, 'k', 'linewidth', 2)
hold on
plot(t, sp, 'r', 'linewidth', 2)
title('Position vs time')
legend('fi' , 's')
figure()
plot(t, dfip, 'k', 'linewidth', 2)
hold on
plot(t, dsp, 'r', 'linewidth', 2)
title('Velocity vs time')
legend('fi' , 's')
figure()
plot(t, ddfip, 'k', 'linewidth', 2)
hold on
plot(t, ddsp, 'r', 'linewidth', 2)
title('Acceleration vs time')
legend('fi' , 's')