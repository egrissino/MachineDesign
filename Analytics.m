%% ME 643 Project 1
%
%% 
%

clear all;%#ok
close all;

load data_UHR.mat

amp = 0.5 * (max(rx) - min(rx));
offset = max(rx) - amp;

phase = 0;
error = Inf;

while error > 0.0001
    phase = phase + 0.0001;
    C = amp * sin( omega * t + phase) + offset;
    error = abs( sum(C - rx) );
    
    if phase >= 2 * pi
        break
    end
end
    

C = amp * sin( omega * t + phase) + offset;


fig = figure();
ax = axes(fig, 'next', 'add', 'fontsize', 14);

plot(ax, t, rx);
% plot(ax, t, C);

xlabel(ax, 'Time (s)');
ylabel(ax, 'Distance (m)');

title('Analytical r_x vs. Numerical r_x');

legend('Numerical', 'Theoretical', 'location', 'best');
grid on

%% Velocity of C
%

vx = diff(rx) ./ dt;

fig2 = figure();
ax2 = axes(fig2, 'next', 'add');

plot(ax2, t(1:length(vx)), vx);

ax = diff(vx) ./ dt;

fig3 = figure();
ax3 = axes(fig3, 'next', 'add');

plot(ax3, t(1:length(ax)), ax);

