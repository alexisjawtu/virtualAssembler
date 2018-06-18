%% assmbl_tetr: function description
function [outputs] = assmbl_tetr(arg, normalFacesE, measFacesE)
	pts_of_faces        = zeros(3,3,4);
  	pts_of_faces(:,:,1) = [.5 .5 0; 0 .5 .5; 0 0 0];
  	pts_of_faces(:,:,2) = [0 0 0; .5 .5 0; 0 .5 .5];
  	pts_of_faces(:,:,3) = [.5 .5 0; 0 0 0; 0 .5 .5];
  	pts_of_faces(:,:,4) = [.5 .5 0; 0 .5 .5; .5 0 .5];

    face_quad_coef = [measFacesE;measFacesE;measFacesE]/3;



	outputs = 1;
endfunction
