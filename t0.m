%% ME 643 : Project 1
%
%% Solve kinematics with Vector Loops
%

clear all;%#ok
close all;

r2 = 0.007;             % AoA       - m
r5 = 0.02;              % BoB
r4 = 0.1;               % BC
rab = 0.045;            % AoBo
H = 0.102;              % H

omega = 2 * pi / 60;    % omega     - rad per sec
m = 0.45;               % mass      - Kg
k = 175;                % k         - N / m
rho = 1070;             % rho       - kg / m^3
ry = H - rab;

t4 = 0.8:0.00001:2;
dt = 0.01;

END = 100;

% syms t4

% Set variable
tic
for i = 1:END/dt
    
    t(i) = dt * (i - 1);
    
    theta2(i) = omega * t(i);
    
%     assume(t4 > 0.9 & t4 < 2)
% 
%     loop = abs( r2 * cos( theta2(i) ) - cot( t4 ) * ( ry - r2 * sin( theta2(i) ) )...
%         + r4 * cos( t4 ) - r5 * cos(...
%         asin( 1/r5 * ( H - r4 * sin( t4 ) ) ) ) ) == 0;
% 
%     [sol, params, cond] = solve(loop, t4, 'returnconditions', true);
%     
%     assume(cond)
%     
%     interval = [sol >= 0.9, sol < 2];
%     
%     Sol = solve(interval, params);

    loop = r2 * cos( theta2(i) ) - cot( t4 ) * ( ry - r2 * sin( theta2(i) ) )...
            + r4 * cos( t4 ) - r5 * cos(...
            asin( 1/r5 * ( H - r4 * sin( t4 ) ) ) );
    
    [~, index] = min( abs(loop) );

    theta4(i) = t4(index);

    theta5(i) = asin( 1/r5 * (H - r4 * sin( theta4(i) ) ) );
    
    r41(i) = csc(theta4(i)) * (ry - r2 * sin(theta2(i)));
    
    rx(i) = r4 * cos( theta4(i) ) - r5 * cos( theta5(i) );
end
toc

%% Save data to .mat file
%

save data_lowRes.mat

%% Plotting
%

close all;

fig = figure();
fig2 = figure();

ax1 = axes(fig, 'next', 'add', 'position', [0.12, 0.6, 0.8, 0.35],...
    'fontsize', 14, 'xgrid', 'on', 'ygrid', 'on');
ax2 = axes(fig, 'next', 'add', 'position', [0.12, 0.12, 0.8, 0.35],...
    'fontsize', 14, 'xgrid', 'on', 'ygrid', 'on');

ax3 = axes(fig2, 'next', 'add',...
    'fontsize', 14, 'xgrid', 'on', 'ygrid', 'on');

plot(ax1, t, rx)
plot(ax1, t, r41)
plot(ax2, t, theta4)
plot(ax2, t, theta5)
plot(ax3, t, theta2)

legend(ax1, 'r_x', 'r_{4,1}', 'location', 'best');
legend(ax2, '\theta_4', '\theta_5', 'location', 'best');
legend(ax3, '\theta_2', 'location', 'best');

xlabel(ax1, 'time (s)');
xlabel(ax2, 'time (s)');
xlabel(ax3, 'time (s)');

ylabel(ax1, 'distance (m)');
ylabel(ax2, 'angle (rad)');
ylabel(ax3, 'angle (rad)');



