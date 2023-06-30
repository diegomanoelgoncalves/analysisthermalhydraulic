function y = interpolate(x1, y1, x2, y2, x)
% Calculate corresponding y-coordinate
y = y1 + (y2-y1)/(x2-x1) * (x-x1);