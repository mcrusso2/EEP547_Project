% housekeeping
clear all
close all
clc

% Digital filter used to smooth the wheel velocity curve
% load('dfFilter1.mat')

% NXT parameters use by the Simulink model
encoder_counts=720;   % number of counts (if using quad encoding)
RPM_MAX = 170;        % spec sheet max RPM

% constants or conversion factors:
RADSEC2RPM = 60/(2*pi);       % radians/sec to RPM
RAD2C=encoder_counts/(2*pi);  % conversion from radians to counts
C2RAD=1/RAD2C; 
C2DEG=360/encoder_counts;
R2D=180/pi;
D2R=1/R2D;
g=9.8;

% Wheel Info
Rw = 0.021336;    % radius of wheel in m (small wheel), under load

% hardware
Vsupply=9;  % 9 volts when on batter, 5 volts on USB
DCB2V=Vsupply/255;   % Duty cycle bits to PWM (volts)
V2DCB=1/DCB2V;       % Volts to duty cycle in bits.


% Gain derived from measurements, using EE547 model and Mathematica's LQR,
% but with the default MinSeg motor parameters.
%KLQR=[10.0, 53.4938, -93.1855, -18.0756];     % 56mm wheel .345kg Q=100


% Gain derived from measurements, using EE547 model and Mathematica's LQR
% and actual measurements.
% KLQR=[12.2474, 55.528, -104.241, -20.7600];   % 56mm wheel .381kg Q=150
%KLQR=[10.9545, 54.285, -95.4137, -18.6104];   % 56mm wheel .381kg Q=120
%KLQR=[10.0000, 53.412, -94.9204, -18.2348];   % 56mm wheel .381kg Q=100
%KLQR=[10.0000, 53.412, -96.1486, -19.4129];   % 56mm wheel .381kg Q=100 (With viscous damping)

% KLQR=[ 9.4868, 52.960, -92.8867, -17.6796];   % 56mm wheel .381kg Q=90
% KLQR=[ 9.4868 ];
    
%KLQR=[ 9.2195, 52.730, -91.8432, -17.3941];   % 56mm wheel .381kg Q=85
%KLQR=[ 7.07107, 51.008, -83.9058, -15.2105];   % 56mm wheel .381kg Q=50
% KLQR=[ 4.47214, 49.2633, -75.7192, -12.9462];   % 56mm wheel .416kg Q=20

% Gain from pole placement, 1X
%KLQR=[28.885, 61.468, -85.655,  -13.872];   % 56mm wheel .381kg
% Gain from pole placement, 1X, closer to the real axis
% KLQR=[ 8.1799, 48.261, -92.6513, -17.0587 ];   % 56mm wheel .381kg.  Less jitter, but more horizontal movement.

% Poles placed  {-1065, -3.6+0.5i, -0.4, -3.6-0.5i}.
% KLQR=[  7.7404, 47.163, -92.5425, -16.8075 ];   % 56mm wheel .381kg. 

% Poles placed [-21, -17, -14, -10];
% KLQR=-1*[ -20.9615   -2.5068   51.3725   37.8606];

%LQR Selected 
% KLQR=[81.6497   74.2481   -93.3999  -13.3790];
KLQR=[ 70.7107   73.7136 -108.1757  -16.8574];

% Quite stable, not so much oscillation

%Poles placed {-1065, -3.8+0.1, -0.6, -3.8-0.1i}.  
%KLQR=[ 13.2162, 52.289,-102.8783, -18.1050 ];   % 56mm wheel .381kg.
% Eventually goes unstable through slow oscillations

% combine some constants:
ES=-C2DEG*D2R;

% automatic calibration of gyro offset
% discrete 1st order filter recurrence relation - 
% discrete-time implementation of a low-pass filter is 
% the exponentially-weighted moving average
% alpha=.5 (time contant is equal to sampling period)
% alpha<.5 (time constant is larger than sampling period) tau ~ TS/alpha

% .006 requires 800 samples to get to steardy state (see gyro filter design)
% 1st order system take 3Tau to get to .95 ss value, 3.9Tau, 98, 4.56tau 99
% tau ~ .0025/.006 = .4167 = 3*.4167 = 1.25 to get to .95  - we need .99
% tau ~ .0025/.006 = .4167 = .4167*.456 = 1.9sec to get to 99% ss
a_go = .006;  % alpha for initial gyro offset calibration (2 sec for ss)

% Start time for balancing - Gyro calibartion time in second:
tstart=2;

%% MPU5060
%TS=.005 % Minimum simulation time
TS=0.005; % Simulation time with motor controller
%TS=.005 % fastest with default mpu5060 library settings... and filter set to DLPF set to 4
GyroS=250/32768;   % MPU5060 set to a maximum rate of 250 deg/second and 16-bit sampling
GS=GyroS*D2R;      % Convert from degrees to radians
tstart=.6;
return

% Low pass filter for output:
a_ov = 0.7;  % on USB
a_ov = 0.5; % on battery
% a_ov = 1;  % no filter.


