function [ rr, vv ] = kep2car( a, e, i, RAAN, omega, theta, mu )

%Transform from Keplerian elements to cartesian coordinates 
% -------------------------------------------------------------------------
%% Input Arguments:  
% [a] semi-major axis [km]
% [e] eccentricity  [-]
% [i] inclination   [rad]
% [RAAN] right ascension of the ascending node    [rad]
% [omega] argument of periapsis   [rad]
% [theta] true anomaly   [rad]
% [mu] gravitational parameter  [km^3/s^2]

% -------------------------------------------------------------------------
% Output arguments: 
% [rr] : [3x1] position vector [km]
% [vv] : [3x1] velocity vector [km/s]

% 1. Compute position and velocity vectors in the perifocal frame

p=a*(1-e^2);          					% Semilatus rectum
r=p/(1+e*cos(theta)); % Magnitude of the position
rrpf=r*[cos(theta), sin(theta), 0]';             	% Position vector expressed in perifocal coordinate system
vvpf=(mu/p)^0.5*[-sin(theta), e+cos(theta), 0]'; 	% Velocity vector expressed in perifocal coordinate system

% 2. Define rotation matrices from perifocal to internal frame

R_RAAN = [cos(RAAN), sin(RAAN), 0; -sin(RAAN), cos(RAAN), 0;0 0 1]'; 		% Rotation matrix around axis 3
R_i=[1 0 0;0 cos(i) sin(i);0 -sin(i) cos(i)]';                       		% Rotation matrix around axis 1
R_omega = [cos(omega), sin(omega), 0; -sin(omega), cos(omega), 0;0 0 1]'; 	% Rotation matrix around axis 3
Rpf_ge= R_RAAN*R_i*R_omega; 				% Global rotation matrix from perifocal to inertial

% 3. Rotate the vectors

rr=Rpf_ge*rrpf; 					% Position vector in cartesian
vv=Rpf_ge*vvpf; 					% Velocity vector in cartesian

end
