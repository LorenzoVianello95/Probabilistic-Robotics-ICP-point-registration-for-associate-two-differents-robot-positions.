
#alternative: hungarian method? min column, best friend, lonely best friend
function [x_visualize,y_visualize, new_landmarks, new_measurements] = data_association(landmarks, measurements,gating_tau, lonely_best_friend_gamma_threshold, w)
	N= size(landmarks,1); # number of landsmark and number of measurements
	M=size(measurements,1);
	x_visualize=[];		y_visualize=[]; 	new_landmarks=[];	 new_measurements=[];
	A=A_matrix(landmarks, measurements, w);

	#now associate the measurement to the most promising landmark
	# proposed associations will be stored in a [Mx3] matrix composed
	# in this way
	#
	#	[measurement id, proposed landmark id , association matrix value a_mn] 	
	#
	# we will populate such a matrix with the associations surviving
	# the gating heuristic for each step
	# 
	# In the best friends and lonely best friend heuristics we will keep
	# the association as is in case of success, otherwise we will put
	# an id=0 to that measurment, meaning that the association is doubtful

	#configuration heuristics
	#gating_tau                         = 1;#1e-3
  	#lonely_best_friend_gamma_threshold = 1e-3;#1e-3;

	#1. Gating
  	for m=1:M

		#return the min and index on the 'm-th' row
		[a_mn, min_index] = min(A(m,:));

    		#if the association passes the gate
		if(a_mn < gating_tau)	
		
			#add the possible association - creating the associations vector
      			#[measurement id, proposed landmark id , association matrix value a_mn] 
			associations(end+1,:) = [m, min_index, a_mn];
		endif
  	endfor

	#associations that survived the gating
	number_of_gated_associations = size(associations, 1);

	#2. Best friends
	for i=1:number_of_gated_associations
		a_mn                 = associations(i, 3);
		proposed_landmark_id = associations(i, 2);

   	 #compute column minimum
		min_on_column = min(A(:, proposed_landmark_id));

   	 #if the association is not the minimum in the column
		if(a_mn != min_on_column)
			associations(i, 2) = 0; #discard association, it is doubtful
		endif
	endfor
	
	#3. Lonely best friend
  	number_of_valid_associations = 0;
	if(M > 1)
		for i=1:number_of_gated_associations
			a_mn                 = associations(i, 3);
			measurement_id       = associations(i, 1);
			proposed_landmark_id = associations(i, 2);

      			#this association is doubtful, skip evaluation
			if(proposed_landmark_id == 0)
				continue;
			endif

			#obtain second best(aka min) value of the row
			ordered_row          = unique(A(measurement_id,:));
			second_row_min_value = ordered_row(2);

			#obtain second best(aka min) value of the column
			ordered_col          = unique(A(:,proposed_landmark_id));
			second_col_min_value = ordered_col(2);

			#check if the association is ambiguous
			if( (second_row_min_value - a_mn) < lonely_best_friend_gamma_threshold ||
		 	    (second_col_min_value - a_mn) < lonely_best_friend_gamma_threshold )

        		#discard association, it is doubtful
				associations(i,2) = 0;
			else
        
        			#we have found a valid association!
        			++number_of_valid_associations;
      			endif
		endfor
	endif

	#assign the associations to the observations
	for i=1:number_of_gated_associations
		if associations(i,2)!=0
			x_visualize=[x_visualize;landmarks(associations(i,2),1);measurements(associations(i,1),1);nan];	
			y_visualize=[y_visualize;landmarks(associations(i,2),2);measurements(associations(i,1),2);nan];
			
			#vengono restituiti come 2*n quindi ogni colonna landmark o meas
			new_landmarks= [new_landmarks, landmarks(associations(i,2),:)'];
			new_measurements=[new_measurements,measurements(associations(i,1),:)'];
			
		endif
	endfor

  printf("valid associations: %u / measurements: %u / landmarks: %u\n", number_of_valid_associations, M, N);

endfunction


#build A matrix that contain likelihood for each landmark/measurement pair 
#per il momento considero solo indici vicini perchè translazioni e rotazionipiccole 
#un idea sarebbe quella di valutare di più misurazioni centrali perchè laterali potrebbero perdersi 
function A = A_matrix(landmarks, measurements, w)
	N= size(landmarks,1); # number of landsmark and number of measurements
	M= size(measurements,1);
	A=ones(max(N,M),max(N,M))*1000;

	sigma_zx = eye(2,2)*0.01;
	sigma_nn = sigma_zx;
	omega_nn = inv(sigma_nn);
	
	if w>0
		intorno=w;
	else
		intorno= max(N,M);
	endif
	for n= 1:N
		#for m= max(n-intorno,1):min(n+intorno,M)	#IF WANT FASTER DATA ASSOCIATION UNCOMMENT, VALID ONLY IN SUCCESSIVE FRAMES
		for m= 1:M
			landmark_position = landmarks(n,:);
			measurement_position = measurements(m,:);
			#a_mn= (landmark_position(1,1)-measurement_position(1,1))^2+(landmark_position(1,2)-measurement_position(1,2))^2
			a_mn= (landmark_position-measurement_position)*omega_nn*(landmark_position-measurement_position)';
			#c=1+(abs(m-M/2)+abs(n-N/2))/M;
			#a_mn*= c;
			A(m,n)= a_mn;
		endfor
	endfor
endfunction

