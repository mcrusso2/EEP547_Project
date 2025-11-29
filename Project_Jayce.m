%% Minseg Robot Control Project EEP 547
% Team Names: Jayce Gaddis, Matthew Russo

%% Minseg Parameters
%Initializing
clear;
close all;
clc;

%Defining parameters
g = 9.81;
syms alpha; % Angle between minseg pendulum and vertical axis
syms x; % traveling distance of the wheel
syms L; % Distance between wheel center and pendulum CoM
syms mp; % Mass of Pendulum
syms Ip; % Moment of inertia at reference
syms mw; % Mass of the wheel
syms rw; % Radius of wheel
syms Icmw; % Moment of inertia at center of mass of wheel
syms kt; % Torque constant in Nm/A
syms kb; % Back EMF constant in Vs/rad
syms R;  % Resistance in ohms


%% Step 1: Setup System Matrices
% Simplifying A and B matrices by defining numerator and denumerator

den1_a = Icmw*(Ip+L^2*mp) + (L^2*mp*mw + Ip*(mp+mw))*rw^2;
num1_a = g*L*mp*(Icmw + (mp+mw)*rw^2);
num2_a = kt*(Icmw + rw*(mw*rw + mp*(L + rw)));
num3_a = g*L^2*mp^2*rw^2;
num4_a = kt*(Ip + L*mp*(L+rw));

A21 = num1_a/den1_a;
A22 = -kb*num2_a/(R*den1_a);
A24 = -kb*num2_a/(R*rw*den1_a);
A41 = num3_a/den1_a;
A42 = -kb*rw*num4_a/(R*den1_a);
A44 = -kb*num4_a/(R*den1_a);

%Defining the matrices A, B, C, and D
A_b = [0 1 0 0;
      A21 A22 0 A24;
      0 0 0 1;
      A41 A42 0 A44];
B_b = [0;
    -num2_a/(R*den1_a);
    0
    -num4_a*rw/(R*den1_a)];
C_b = eye(size(A_b));
D_b = [0;0;0;0];

%% Step 2: Measurements of Minseg
% Measure the physical parameters of our MinSeg system. The _kt_, _Kb_, and
% _R_ values were provided for us. The rest of the parameters were
% calculated using various methods. While lengths and masses were measured
% with rulers and scales, the hardest to calculate was the moment of
% inertias. 
% 
% The wheel's inertia was calculated using a ramp and finding the
% acceleration of the wheel. That acceleration was used to backcalculated
% the the inertia. The pendulum's inertia was calculated by using an
% equation with relation to its length and mass.

% Without Battery			
%  Length to CoG (with wheels)	L	95.83	mm
%  Length to CoG	L	100	mm
mp_m = 201/1000;
L_m = 100/1000;


Ip_m = mp_m*L_m^2;
% With Battery			
%  Length to CoG (with wheels)	L	111.78	mm	
mp_m_b = 339/1000;
L_m_b = 111.78/1000;
measured_period_with_batt = 0.91;
Ip_m_b = mp_m_b * g * L_m_b * (measured_period_with_batt / (2 * pi))^2;
%fprintf("Measured Moment of Inertia: %.3e kg·m²\n", Ip_m_b);
fprintf("Pendulum measured Moment of Inertia No Batt: %.3e kg·m²\n", Ip_m);

% Wheel parameters
dw_m = 42;                % wheel diameter in mm
rw_m = dw_m / 2 / 1000;   % radius in meters

mw_kg = 18 / 1000;        % mass in kg
time_down_ramp_s = 1.89;  % time to roll distance (s)
ramp_angle_deg = 3;       % incline angle (degrees)
ramp_length_m = 0.508;    % ramp distance (m)

% Compute acceleration down the ramp
a = 2 * ramp_length_m / time_down_ramp_s^2;   % m/s^2

% Compute wheel moment of inertia about center (two wheels total)
Icmw_m = 2 * mw_kg * rw_m^2 * ((g * sind(ramp_angle_deg)) / a - 1);
fprintf('Wheel inertia Icmw = %.3e kg·m²\n', Icmw_m);

%  Motor			
%   Torque Constant	kt	0.3233	Nm/A
%   Back EMF	Kb	0.4953	Vs/rad
%   Resistance	R	5.2628	ohms
kt_m = 0.3233;
kb_m = 0.4953;
R_m = 5.2628;

% Build substitution pairs
sub_pairs = [ mp  Ip   L   mw   rw   Icmw   kt    kb    R ;
              mp_m  Ip_m  L_m  mw_kg rw_m  Icmw_m   kt_m  kb_m  R_m ];

sub_pairs_batt = [ mp  Ip   L   mw   rw   Icmw   kt    kb    R ;
              mp_m_b  Ip_m_b  L_m_b  mw_kg rw_m  Icmw_m   kt_m  kb_m  R_m ];

% Substitute measured values into A and B
A_b = subs(A_b, sub_pairs(1,:), sub_pairs(2,:));
B_b = subs(B_b, sub_pairs(1,:), sub_pairs(2,:));

A_b = double(vpa(A_b,14));
B_b = double(vpa(B_b,14));
disp("A= ");disp(A_b);
disp("B= ");disp(B_b);
			
%% Step 3: Transfer Function
% Transfer function matrix of the linearized system. The system possesses
% four transfer functions, all sharing the same denominator.

[tfnum,tfden] = ss2tf(A_b,B_b,C_b,D_b);

tf1 = tf(tfnum(1,:),tfden);
tf2 = tf(tfnum(2,:),tfden);
tf3 = tf(tfnum(3,:),tfden);
tf4 = tf(tfnum(4,:),tfden);

%% Step 4: Characteristic Polynomial and Eigenvalues
% The characteristic polynomial and eigenvalues of matrix A.

charPoly = charpoly(A_b);
rootPoly = roots(charPoly);
eigA_b = eig(A_b);

disp("Characteristic Polynomial: "); disp(charPoly);
disp("Roots: ");disp(rootPoly);
disp("Eigenvalues: ");disp(eigA_b);

%% Step 5: Check for Stability
% Is the system asymptotically stable? Is it marginally stable? It is
% neither, the system is unstable. The system has four eignvalues, three are
% negative, while one is positive. This positive eignvalue is the reason the
% system in unstable.

if eigA_b < 0
    fprintf('The system is asymptotically stable\n')
else
    fprintf('The system is not asymptotically stable\n')
end

%% Step 6:  Check Transfer Function for Stability
% The system is also not BIBO stable. For this system, the eignvalues and
% the poles are the same value. Since a single pole is positive, the system
% is not BIBO stable.

poles_tf1 = pole(tf1);

if poles_tf1 < 0
    fprintf('The system is BIBO stable\n')
else
    fprintf('The system is not BIBO stable\n')
end


%% Step 7: Check for Controllabilty
% To find the controllability of the system, a control matrix must be
% calculated. Comparing the rank of this matrix with the length of the A
% matrix will determine controllability. 
%
% $$ Uncontrollable - ctrb_rank < length(A) $$

ctrb_b = ctrb(A_b,B_b);
ctrb_b_rank = rank(ctrb_b);

if ctrb_b_rank < length(A_b)
    fprintf('The system is not controllable\n')
else
    fprintf('The system is controllable\n')
end

%% Step 8: Check for Observability
% To find the observability of the system, a observation matrix must be
% calculated. Comparing the rank of this matrix with the length of the A
% matrix will determine observability. 
%
% $$ Unobservable - obsv_rank < length(A) $$

obsv_b = obsv(A_b,C_b);
obsv_b_rank = rank(obsv_b);

if obsv_b_rank < length(A_b)
    fprintf('The system is not observable\n')
else
    fprintf('The system is observable\n')
end

%% Step 9:  Transform into Canonical Form
% Transforming the state-space matrixes into CCF and OCF.
sys = ss(A_b, B_b, C_b, D_b);
CCF = compreal(sys,"c");
CCF.A = CCF.A';
CCF.B = flip(CCF.B);
CCF.C(1,:) = flip(tfnum(1,2:5));
CCF.C(2,:) = flip(tfnum(2,2:5));
CCF.C(3,:) = flip(tfnum(3,2:5));
CCF.C(4,:) = flip(tfnum(4,2:5));

OCF = ss(CCF.A',CCF.C',CCF.B',CCF.D');

disp('A Obs:');disp(OCF.A);
disp('A Cont:');disp(CCF.A);

%% Step 10: Designing State Estimator with Pole Placement
% Develop a closed-loop state estimator (full-dimensional observer) for the open-loop system (no feedback
% yet) such that the poles of the observer are stable, and the dynamics of the observer is at least 6-8 times
% faster than the dynamics of the linearized model.
%  Include in your report the value of the estimator gain, L.

% Try ~8 times faster than system
desiredPoles = [-24+0.1j,-24-0.1j, -20, -10]*1.5;

% Calculate the estimator gain L using pole placement
L = place(A_b', C_b', desiredPoles)';
disp('Estimator Gain L:'); disp(L);

% Verify the Pole placement
A_obs = A_b - L*C_b;
disp("New Poles: ");disp(eig(A_obs));

%% Step 11: State Estimator Simulink
% Develop a Simulink model of the linearized system (open-loop system with full-dimensional observer
% designed above). Add the state estimator derived in previous step to your Simulink model and set the initial
% conditions of the state estimator to 𝑥M0 = [0 0 0 0]. Simulate the behavior of the system when a unit-step
% u(t) = 1, t ≥ 0 is applied at the input. Plot the estimated state-variables and output variables on the same
% graph.
sim("StateEstimator_project");
figure(1);
plot(t, x, 'LineWidth', 2.5);
hold on;
plot(t,xhat, '*-');
grid on;
xlabel('time (sec)');
title('Observer Step Input');

%% Step 12: Feedback via Pole Placement
% Consider the case when the linearized system is stabilized by using feedback control. Using the pole
% placement method, via MATLAB, develop a proportional controller such that the poles of the closed-loop
% system are stable, and the dynamics of the closed-loop model is at least 4-6 times faster than the dynamics
% of the open-loop model. Include in your report the value of the proportional gain, K.
desiredPropGainPoles = [-21, -17, -14, -10];
K_pg = place(A_b, B_b, desiredPropGainPoles);

disp("Step 12- Proportional Gain Matrix: "); disp(K_pg);

%% Step 13: Closing the Feedback Loop
% Derive the state-space representation of the closed-loop system. Find the characteristic polynomial and the
% eigenvalues of the closed-loop system. Is this closed-loop system asymptotically stable?
Acl = A_b - B_b*K_pg;
charPoly_cl = charpoly(Acl);
disp("Closed Loop Charactaristic Polynomial: ");disp(charPoly_cl);

eigenValues_cl = eig(Acl);
disp("Closed Loop Eigenvalues: ");disp(eigenValues_cl);

if eigenValues_cl < 0
    fprintf('The system is asymptotically stable\n')
else
    fprintf('The system is not asymptotically stable\n')
end

%% Step 14: Feedback Simulink Model and Step Response
% Develop a Simulink model of the linearized closed-loop system (no estimator) when the output of the
% system equals state variables 𝑦 = 𝑥=[ 𝛼 𝛼̇ 𝑥 𝑥̇]/. Add the proportional controller developed above to the
% Simulink model and simulate the response of the closed-loop system when a unit step u(t) = 1, t ≥ 0 is
% applied. Plot all the outputs on the same graph.
% Create a simulink model of the system with full state feedback and the
% gains calculated previously.
sim("ProportionalGainController_NO_ESTIMATOR_project");
figure(2);
plot(t, y, 'LineWidth',2.5);
grid;
title("Proportional Gain Controller (No State Estimator)");
legend('x_p_o_s', 'x_v_e_l', 'a_p_o_s', 'a_v_e_l')
xlabel('Time (sec)')



%% Step 15: Estimator With Feedback
% Combine the proportional feedback controller with the state estimator from 4.3 and 4.4 in Simulink.
% Simulate the system in to analyze its performance. Try using the estimator for states not measured (this will
% change your y(t)). Plot the error function and discuss your results.
sim("ProportionalGainController_w_Estimator_project");
figure(3);
plot(t, y, 'LineWidth',2.5);
grid;
title("Proportional Gain Controller (State Estimator)");
legend('x_p_o_s', 'x_v_e_l', 'a_p_o_s', 'a_v_e_l');
xlabel('Time (sec)');


%% Step 16: PID Tuning
% A) Demonstrate the feedback control system using an LQR controller.
% B) Demonstrate the feedback control system using a PID controller (You could try PI, PD, and PID.).
% Using the model below, we can use MATLAB's built in PID functions to give
% us our proportional, integral, and derivative gains. Additionally we can
% calculate the poles and zeros to see if the system with PID feedback is
% stable.



%% Step 17: LQR Tuning
% Demonstrate the feedback control system using an LQR controller. 
% Show that the LQR controller balances your MinSeg robot. 
% Show you separate Simulink model. 
% This can be demonstrated in a video and live during the project presentation.
R = [75];

Q = diag([10000, 1, 500000, 5000]);
KLQR = lqr(A_b, B_b, Q, R);
disp("KLQR: ");disp(KLQR);

lqr_poles = eig(A_b - B_b*KLQR)

K_pg = KLQR;

sim("ProportionalGainController_NO_ESTIMATOR_project");
figure(4);
plot(t, y, 'LineWidth',2.5);
grid;
title("LQR (State Estimator)");
legend('x_p_o_s', 'x_v_e_l', 'a_p_o_s', 'a_v_e_l');
xlabel('Time (sec)');

%% Step 18: Sonar Sensor
% Use a sonar sensor to implement some type of feedback control. 
% One option would be to implement an alert system on the balancing MinSeg robot, 
% which generates sound as the distance between MinSeg robot and an object is less 
% than certain threshold by using an ultrasonic sensor and a speaker. 
% As the balancing robot approaches an object, the alert system should 
% increase the volume, frequency, or tone of an alerting sound. 
% This part is open for you to be creative to use the sonar sensor. 
% Demonstrate this in a video or live during the project presentation.




