%% assmbl_prism: function description
function [outputs] = assmbl_prism(arg, normalFacesE, measFacesE)
	prismatic_face_points          = zeros(3,9,5);
	prismatic_face_points(:,:,1)   = [0 .5 1 0 .5  1 0 .5 1; 0 0  0  0  0  0 0 0 0; 0 0  0 .5 .5 .5 1 1 1];
	prismatic_face_points(:,:,2)   = [1 .5 0 1 .5 0 1 .5 0; 0 .5 1 0 .5 1 0 .5 1; 0 0 0 .5 .5 .5 1 1 1];
	prismatic_face_points(:,:,3)   = [0 0 0 0 0 0 0 0 0; 0 .5 1 0 .5 1 0 .5 1; 0 0 0 .5 .5 .5 1 1 1];
	prismatic_face_points(:,1:3,4) = [.5 .5 0; 0 .5 .5; 0 0 0];
	prismatic_face_points(:,1:3,5) = [.5 .5 0; 0 .5 .5; 1 1 1];
	prismatic_face_points


    face_quad_coef   = [(prism_face_coefs*measFacesE(1:3))/36, [repmat(measFacesE(4:5)/3,3,1);zeros(6,2)]];



	outputs = 1;
endfunction
