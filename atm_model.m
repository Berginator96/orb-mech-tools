function rho=atm_model(h)

% Compute exponential atmospheric density model. Source: Wertz, 1979,820
% -------------------------------------------------------------------------
%% Input Arguments:  
% h [1x1] altitude [km]
%
% -------------------------------------------------------------------------
% Output arguments: 
% rho: density [kg/km^3]
%

% data is structured in the form : [minimum interval altitude [km] ,maximum interval altitude [km], base altitude [km], nominal density [kg/m^3], scale height [km]]

data=[0, 25, 0, 1.225, 7.249;
    25, 30, 25, 3.899e-2, 6.349;
    30, 40, 30, 1.774e-2, 6.682;
    40, 50, 40, 3.972e-3, 7.554;
    50, 60, 50, 1.057e-3, 8.382;
    60, 70, 60, 3.206e-4, 7.714;
    70, 80, 70, 8.770e-5, 6.549;
    80, 90, 80, 1.905e-5, 5.799;
    90, 100, 90, 3.396e-6, 5.382;
    100, 110, 100, 5.297e-7, 5.877;
    110, 120, 110, 9.661e-8, 7.263;
    120, 130, 120, 2.438e-8, 9.473;
    130, 140, 130, 8.484e-9, 12.636;
    140, 150, 140, 3.845e-9, 16.149;
    150, 180, 150, 2.070e-9, 22.523;
    180, 200, 180, 5.464e-10, 29.740;
    200, 250, 200, 2.789e-10, 37.105;
    250, 300, 250, 7.248e-11, 45.546;
    300, 350, 300, 2.418e-11, 53.628;
    350, 400, 350, 9.158e-12, 53.298;
    400, 450, 400, 3.725e-12, 58.515;
    450, 500, 450, 1.585e-12, 60.828;
    500, 600, 500, 6.967e-13, 63.822;
    600, 700, 600, 1.454e-13, 71.835;
    700, 800, 700, 3.614e-14, 88.667;
    800, 900, 800, 1.170e-14, 124.64;
    900, 1000, 900, 5.245e-15, 181.05;
    1000, 1e4, 1000, 3.019e-15, 268.00];


% Identify the row index of interest

index=1;

while h>=data(index,2)
    
 index=index+1;
 
end      

% Compute density with exponential formula

h_base=data(index,3);
rho_0=data(index,4);
H=data(index,5);
rho=rho_0*exp(-(h-h_base)/H)*10^9;

end
    