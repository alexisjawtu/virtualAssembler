%% are_coplanar: function description
function [bit] = counter_clockwise(P)
	bit 	= (dot(P(:,2)-P(:,1),P(:,4)-P(:,1)) == 0);
endfunction