% drawline_thick(ax, ch1, ch2, color, lw)
%   ax    - axes handle
%   ch1, ch2 - channel indices (1-based, into the 60-channel layout)
%   color - [R G B] or [R G B A]
%   lw    - line width
function drawline_thick(ax, ch1, ch2, color, lw)

load systems.mat;
t = systems(18).layout;
t([5 10 21 27],:) = [];

plot(ax, [t(ch1,1) t(ch2,1)], [t(ch1,2) t(ch2,2)], ...
    'Color', color, 'LineWidth', lw);
end
