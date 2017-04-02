%% ME 643 : Project 1
%
%% Solve kinematics with Vector Loops
%

clear all;%#ok
close all;

%% Machine Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r2      = 0.007;                        % AoA    
r5      = 0.02;                         % BoB
r4      = 0.1;                          % BC
rab     = 0.045;                        % AoBo
H       = 0.102;                        % H
h       =  .010;                        % Rod Height
w       =  .010;                        % Rod Width
rho     =  1070;                        % rho       - kg / m^3
I4      =  (rho*r4*w*h*(r4^2))/3;       % Mass Moment of Interia piece 4
omega   =  2 * pi / 60 ;                % omega     - rad per sec
m       =  0.45;                        % mass      - Kg
k       =  175;                         % k         - N / m
ry      =  H - rab; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Transient and Solution Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.1;
t4_high = 2;
t4_low = 0.8;
diff_t = Inf;
min_diff = 0.00000001 * dt;

END = 60;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Loop Solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
for i = 1:END/dt
    
    t(i) = dt * (i - 1);
    theta2(i) = omega * t(i); 
    
    t4_high = 2;
    t4_low = 0.8;
    diff_t = Inf;
    
    while diff_t > min_diff
        t4 = (t4_high + t4_low)/2;
        
        loop = 100 * (r2 * cos( theta2(i) ) - cot( t4 ) * ( ry - r2 * sin( theta2(i) ) )...
            + r4 * cos( t4 ) - r5 * cos(...
            asin( 1/r5 * ( H - r4 * sin( t4 ) ) ) ) );
        
        diff_t = loop;
        
        if diff_t > 0
            t4_low = t4;
        else
            t4_high = t4;
        end
        diff_t = abs(diff_t);
    end
    
    theta4(i) = t4;
    theta5(i) = asin( 1/r5 * (H - r4 * sin( theta4(i) ) ) );
    r41(i) = csc(theta4(i)) * (ry - r2 * sin(theta2(i)));
    rx(i) = r4 * cos( theta4(i) ) - r5 * cos( theta5(i) );
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Angular Velocities and accelerations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Omega4     =    smooth( diff(theta4)/dt, 25)';
Omega5     =    smooth( diff(theta5)/dt, 25)';
alpha4     =    smooth( diff(Omega4)/dt, 25)';
alpha5     =    smooth( diff(Omega5)/dt, 25)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Kinematics at Joint C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
velocityC      =  smooth(diff(rx)/dt,25)';
accelerationC  =  smooth(diff(velocityC)/dt,25)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Truncate Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = t(3:end);
theta2 = theta2(3:end);
theta4 = theta4(3:end);
theta5 = theta5(3:end);
r41 = r41(3:end);
rx = rx(3:end);

velocityC = velocityC(2:end);
Omega4 = Omega4(2:end);
Omega5 = Omega5(2:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Forces at Joint C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Force64x  =  k * rx - m * accelerationC;
Force64y  =  Force64x .* tan( theta4 );

Force64mag    =  ((Force64x .^2 + (Force64y .^2))).^(1/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Distance between Middle Slider and Joint B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rBAy = r4 - r41 .* sin( theta4 );
rBAx = r4 - r41 .* cos( theta4 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Forces at Joint A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Force32 = -1 ./ (r4 - r41) .* (...
    I4 * accelerationC + Force64x .* r4 .* sin( theta4 )...
    + Force64y .* r4 .* cos( theta4 ) );

Force32x   =  Force32 .* sin( theta4 );
Force32y   =  Force32 .* cos( theta4 );
Force32mag =  abs( Force32 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Forces at Joint B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Force54 = (I4 .* alpha4 + Force32x .* rBAy + Force32y .* rBAx ) ./ (r4 .* (cos((theta5) .* sin(theta4)) + (sin(theta5) .* cos(theta4))));
Force54x = Force54 .* cos(theta5);
Force54y = Force54 .* sin(theta5);
Force54mag = abs(Force54);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Forces at Joint Bo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Reaction2x     =      -Force54x;
Reaction2y     =      -Force54y;
Reaction2Mag   =   abs(Force54);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Forces at Joint Ao %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Reaction5x     =      -Force32x;
Reaction5y     =      -Force32y;
Reaction5Mag   =   abs(Force32);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






