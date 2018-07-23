%% f: problem data
%% x is the matrix containing the coordinates of the volume
%% quadrature points in an element.
%%
%% size(x) == (3,N_vol_pts_E)
%%
%% y is the vector s.t. y(i) == f(x(:,i))
function [y] = f(x)
	N_vol_pts_E = size(x,2):
	y = sin(x(1,:).*x(2,:)).+x(3,:);

