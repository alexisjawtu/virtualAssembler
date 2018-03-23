%% are_coplanar: function description
function [bit] = are_coplanar(P)
	N 		= cross(P(:,2)-P(:,1),P(:,3)-P(:,1));
	const 	= N.'*P(:,1);
	bit 	= (abs(const - N.'*P(:,4)) < 10e-4);
endfunction