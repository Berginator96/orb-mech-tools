function DT = kep_tof (a,e,th_i,th_f, mu)

% Compute the time of flight along a keplerian orbit between two angular positions
% -------------------------------------------------------------------------
% Input arguments:
% a [1x1] semi-major axis [km]
% e [1x1] eccentricity [-]
% th_i [1x1] initial angular position (true anomaly) along the orbit [rad]
% th_f [1x1] final angular position (true anomaly) along the orbit [rad]
% mu [1x1] gravitational parameter [km^3/s^2]

% -------------------------------------------------------------------------
% Output arguments: 
% DT [1x1] time of flight [s] 


% 1. Eccentric Anomaly and time law functions 

E=@(th) 2*atan( sqrt((1-e)/(1+e))*tan(th/2));
t=@(th) sqrt(a^3/mu)*(E(th)-e*sin(E(th)));


% 2. Set angles correctly

th_i=normAngle(th_i);
th_f=normAngle(th_f);
if (th_i==2*pi)
th_i=0;
end
if (th_f==0)
th_f=2*pi;
end


% 3. Compute all times of flight

t1=t(th_i);
t2=t(th_f);

DT=t2-t1;

if (DT<0)
T=2*pi*sqrt(a^3/mu);
DT=T+DT;



function th=normAngle (t)

% Normalize angle between 0 and 2*pi rad
% -------------------------------------------------------------------------
% Input arguments:
% t [1x1] angle [rad]
% -------------------------------------------------------------------------
% Output arguments: 
% th [1x1] normalized angle [rad]


if (t<0)
    t=2*pi+t;
end
while t>2*pi
    t=t-2*pi;
end
th=t;
end


end
