% function [ a, e, i, RAAN, omega, theta ] = car2kep( rv, vv, mu )
function X = c2k( rv, vv, mu )

% Transformation from cartesian coordinates to Keplerian elements 
% -------------------------------------------------------------------------
% Input arguments:
% rr [3x1] position vector[km]
% vv [3x1]  velocity vector [km/s]
% mu [1x1] gravitational parameter [km^3/s^2]

% -------------------------------------------------------------------------
% Output arguments: 
% a [1x1] semi-major axis [km]
% e [1x1] eccentricity [-]
% i [1x1] inclination [rad]
% RAAN [1x1] Right Ascension of ascending node [rad]
% omega [1x1] argument of the perigee [rad]
% theta [1x1] true anomaly [rad]

%% Procedure to obtain orbital elements from state vector
% 1 and 2. Determine the magnitude of r and v
r = norm(rv);
v = norm(vv);

% 3 and 4. Angular momemtum (vector and magnitude)
hv = cross(rv,vv);
h = norm(hv);

% 5. Inclination
i = acos(dot(hv,[0 0 1])/h);

% 6 and 7. Line of nodes direction

Nv = cross([0 0 1],hv);
N = norm(Nv); % Magnitude
Nv = Nv/N;

% 8. Right ascencion of ascending node
if Nv(2)>=0
    RAAN = acos(Nv(1));
else
    RAAN = 2*pi-acos(Nv(1));
end

% 9 and 10. Eccentricity vector
vr = dot(vv,rv)/r; % Radial velocity
ev = 1/mu*((v^2-mu/r)*rv-r*vr*vv);
e = norm(ev);

% 11. Argument of the perigee

if ev(3)>=0
    omega = acos(dot(Nv,ev)/e);
else
    omega = 2*pi-acos(dot(Nv,ev)/e);
end


% 12. True anomaly
if dot(rv,vv)>=0
    d = dot(ev,rv)/e/r;
    if abs(d - 1) < 1e-10
        d=1;
    end
    theta = acos(d);
else
    d = dot(ev,rv)/e/r;
    if abs(d - 1) < 1e-10
        d=1;
    end
    theta = 2*pi-acos(d);
end

% 13. Specific Energy and Semimajor axis
E = norm(vv)^2/2-mu/norm(rv);   % Specific Energy
a = -mu/2/E; 			% Semi major axis

X =  [a, e, i, RAAN, omega, theta ] ;
end





        