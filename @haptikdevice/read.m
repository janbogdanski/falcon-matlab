function [orientation,velocity,buttons] = read(h)
% READ Get data from device
%
% [orientation,velocity,buttons] = read(h)
%
[orientation,velocity,buttons] = haptik_matlab(4,h.id);
