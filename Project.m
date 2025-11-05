%% Minseg Robot Control Project EEP 547 
% Team Names: Matthew Russo & Jayce

%% Minseg Parameters
%Initializing
clear;
close all;
clc;

%Defining parameters
%
%
%
%
syms Ip Ic mp mw L r g;
syms kt kb R V;
syms Z;
kt = 0.3233; %Nm/A
kb = 0.4953; %Vs/rad
R = 5.2628; %ohms

r = 0.021; %m (~42mm diameter)
L = 0.197; %m (from center orange to end of black board)
g = 9.8; % meters/sec^2

mp = 0.326; %Kgrams
mw = 0.01417; %Kgrams
Ic = 0.5*mw*(r^2); %1/2 * mw *radius^2
%%%%%%^do we multiply by 2 because 2 wheels????
Ip = mp*(L^2)/3; %1/3 * mp * L^2


Z = Ic*Ip + Ic*L^2*mp + Ip*mp*r^2 + Ip*mw*r^2 + L^2*mp*mw*r^2;
%% Step 1: Setup System Matrices
% Simplifying A and B matrices by defining numerator and denumerator

den1_b = [R*Z];
num1_b = [];
num2_b = [-(kt*(Ic + mp*r^2 + mw*r^2 + L*mp*r))];
num3_b = [];
num4_b = [-(kt*r*(mp*L^2 + mp*r*L + Ip))];

num21_a = (L*g*mp*(Ic + mp*r^2 + mw*r^2));
num22_a = -(kb*kt*(Ic + mp*r^2 + mw*r^2 + L*mp*r));
num23_a = [];
num24_a = -(kb*kt*(Ic + mp*r^2 + mw*r^2 + L*mp*r));
num41_a = (L^2*g*mp^2*r^2);
num42_a = -(kb*kt*r*(mp*L^2 + mp*r*L + Ip));
num43_a = [];
num44_a = -(kb*kt*(mp*L^2 + mp*r*L + Ip));

%Defining the matrices A, B, C, and D
A_b = [0 1 0 0; num21_a/Z num22_a/(R*Z) 0 num24_a/(R*r*Z) ; 0 0 0 1; num41_a/Z num42_a/(R*Z) 0 num44_a/(R*Z)];
B_b = [0; num2_b/den1_b ; 0; num4_b/den1_b];
C_b = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
D_b = [0; 0; 0; 0];

%% Step 2: Measurements of Minseg
% Measure the physical parameters of your MinSeg as shown in Figure 2. Find creative ways to measure the
% weight in grams (home kitchen scale, post office, grocery store, etc.). Create a table to list the values of
% your measurement in SI unit. Set reference values of (kt, Kb, R) as (kt, Kb, R) = (0.3233 Nm/A, 0.4953
% Vs/rad, 5.2628 ohms).
% Example:
% We used a ramp to measure the wheels acceleration, then backcalculated
% the moment of inertia. Likewise, we found the pendulum's fundamental
% frequency, which gave us the period. We then used that to backcalculate
% the moment of inertia of the pendulum. These steps were excluded from the
% report to reduce length.

% Without Battery			
%  Length to CoG (with wheels)	L	95.83	mm
%  Length to CoG	L	100	mm
%  Diameter of Wheel	dw	40	mm
%  radius of Wheel	rw	20	mm		

% With Battery			
%  Length to CoG (with wheels)	L	111.78	mm
%  diameter of wheel	dw	40	mm
%  radius of wheel	rw	20	mm			
%  Motor			
%   Torque Constant	kt	0.3233	Nm/A
%   Back EMF	Kb	0.4953	Vs/rad
%   Resistance	R	5.2628	ohms
			
%% Step 3: Transfer Function
% Find the transfer function matrix of the linearized system.

[tfnum,tfden] = ss2tf(A_b,B_b,C_b,D_b);

tf1 = tf(tfnum(1,:),tfden);
tf2 = tf(tfnum(2,:),tfden);
tf3 = tf(tfnum(3,:),tfden);
tf4 = tf(tfnum(4,:),tfden);

%% Step 4: Characteristic Polynomial and Eigenvalues
% Find the characteristic polynomial and eigenvalues of matrix A.

charPoly = poly(A_b);
rootPoly = pole(tf1); %same for all 4 tf
eigA_b = eig(A_b);

%% Step 5: Check for Stability
% Is the system asymptotically stable? Is it marginally stable? Explain why.

if eigA_b < 0
    fprintf('The system is asymptotically stable\n')
else
    fprintf('The system is not asymptotically stable\n')
end

%% Step 6:  Check Transfer Function for Stability
% Find the poles of the transfer function. Is the system BIBO stable? Explain why.

poles_tf1 = pole(tf1);

if poles_tf1 < 0
    fprintf('The system is BIBO stable\n')
else
    fprintf('The system is not BIBO stable\n')
end


%% Step 7: Check for Controllabilty
% Find the controllability matrix of the linearized system. What is the rank of the controllability matrix? Is the
% linearized system controllable?

ctrb_b = ctrb(A_b,B_b);
ctrb_b_rank = rank(ctrb_b);

if ctrb_b_rank < length(A_b)
    fprintf('The system is not controllable\n')
else
    fprintf('The system is controllable\n')
end

%% Step 8: Check for Observability
% Analyze the observability of the linearized system with the output vector as 𝑦 = 𝑥 =[ 𝛼 𝛼̇ 𝑥 𝑥̇]/ . Is the
% linearized system observable?

obsv_b = obsv(A_b,C_b);
obsv_b_rank = rank(obsv_b);

if obsv_b_rank < length(A_b)
    fprintf('The system is not observable\n')
else
    fprintf('The system is observable\n')
end

%% Step 9:  Transform into Canonical Form
% Transform the matrices into CCF and OCF. Note matlab returns the
% observable canonical form when given companion.






%% Step 10: Designing State Estimator with Pole Placement
% Develop a closed-loop state estimator (full-dimensional observer) for the open-loop system (no feedback
% yet) such that the poles of the observer are stable, and the dynamics of the observer is at least 6-8 times
% faster than the dynamics of the linearized model. Include in your report the value of the estimator gain, L.




%% Step 11: State Estimator Simulink
% Develop a Simulink model of the linearized system (open-loop system with full-dimensional observer
% designed above). Add the state estimator derived in previous step to your Simulink model and set the initial
% conditions of the state estimator to 𝑥M0 = [0 0 0 0]. Simulate the behavior of the system when a unit-step
% u(t) = 1, t ≥ 0 is applied at the input. Plot the estimated state-variables and output variables on the same
% graph.




%% Step 12: Feedback via Pole Placement
% Consider the case when the linearized system is stabilized by using feedback control. Using the pole
% placement method, via MATLAB, develop a proportional controller such that the poles of the closed-loop
% system are stable, and the dynamics of the closed-loop model is at least 4-6 times faster than the dynamics
% of the open-loop model. Include in your report the value of the proportional gain, K.



%% Step 13: Closing the Feedback Loop
% Derive the state-space representation of the closed-loop system. Find the characteristic polynomial and the
% eigenvalues of the closed-loop system. Is this closed-loop system asymptotically stable?




%% Step 14: Feedback Simulink Model and Step Response
% Develop a Simulink model of the linearized closed-loop system (no estimator) when the output of the
% system equals state variables 𝑦 = 𝑥=[ 𝛼 𝛼̇ 𝑥 𝑥̇]/. Add the proportional controller developed above to the
% Simulink model and simulate the response of the closed-loop system when a unit step u(t) = 1, t ≥ 0 is
% applied. Plot all the outputs on the same graph.
% Create a simulink model of the system with full state feedback and the
% gains calculated previously. 
 




%% Step 15: Estimator With Feedback
% Combine the proportional feedback controller with the state estimator from 4.3 and 4.4 in Simulink.
% Simulate the system in to analyze its performance. Try using the estimator for states not measured (this will
% change your y(t)). Plot the error function and discuss your results.



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




%% Step 18: Sonar Sensor
% Use a sonar sensor to implement some type of feedback control. 
% One option would be to implement an alert system on the balancing MinSeg robot, 
% which generates sound as the distance between MinSeg robot and an object is less 
% than certain threshold by using an ultrasonic sensor and a speaker. 
% As the balancing robot approaches an object, the alert system should 
% increase the volume, frequency, or tone of an alerting sound. 
% This part is open for you to be creative to use the sonar sensor. 
% Demonstrate this in a video or live during the project presentation.




