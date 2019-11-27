%% check_rectangles: function description
%% checks whether the rectangle is:  a--b || c--d
%% with 4,a,b,c,d being a line in faces.txt 
function [bool] = check_rectangles(faces,vertices)

	corre = 0;
	corre2 = 0;
	corre3 = 0;
	corre4 = 0;
	corre5 = 0;
	corre6 = 0;
	tot = 0;

	F = size(faces,1);
	for f = 1:F
		if faces(f,1) == 4
			tot += 1;
			if (norm( vertices(faces(f,2),:) - vertices(faces(f,3),:) ) == norm(vertices(faces(f,4),:) - vertices(faces(f,5),:)))
				corre6 += 1;
			end
			if (norm( vertices(faces(f,5),:) - vertices(faces(f,3),:) ) == norm(vertices(faces(f,4),:) - vertices(faces(f,2),:)))
				corre5 += 1;
			end
			if (norm( vertices(faces(f,4),:) - vertices(faces(f,3),:) ) == norm(vertices(faces(f,5),:) - vertices(faces(f,2),:)))
				corre4 += 1;
			end
			if (norm( vertices(faces(f,4),:) - vertices(faces(f,3),:) ) > norm(vertices(faces(f,3),:) - vertices(faces(f,2),:)))
				corre3 += 1;
			end
			if (norm( vertices(faces(f,4),:) - vertices(faces(f,3),:) ) > norm(vertices(faces(f,4),:) - vertices(faces(f,2),:)))
				corre2 += 1;
			end
			if cross (vertices(faces(f,4),:)-vertices(faces(f,2),:), vertices(faces(f,5),:)-vertices(faces(f,3),:)) == 0
				corre +=1;
			end
	
		end
	end

	Tot = tot*ones(1,6);
	aux = [corre corre2 corre3 corre4 corre5 corre6];
	
	if all(Tot == aux) ~= 1
		'check_rectangle diff'
		Tot
		aux
	end

	bool = all(Tot == aux);

endfunction