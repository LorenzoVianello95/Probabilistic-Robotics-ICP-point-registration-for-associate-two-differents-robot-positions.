#CREA G2O FILE

source "./data_loader.m"
source "./laser2map.m"
source "./data_association.m"
source "common_stuff.m"
source "icp_2d_manifold.m"

fast= false; # se true fa solo le prime 50 iterazioni

global observations= loadData('laser_scans.txt', fast);
odometry= loadOdometry('out.txt', fast);
printf("Observations extracted")

#FUNCTION THAT GIVEN A TRANSFORM AND LANDMARKS MOVE THE LANDMARK IN BASE OF THE TRANSFORM
function newpositions= movePoints(transform, points)
	newpositions= [];
	for index=1:size(points,1)
		old_x=points(index, 1);
		old_y= points(index, 2);
		old_position=  [old_x; old_y; 1];
		new_p= transform*old_position;
		newpositions=[newpositions; [new_p(1) new_p(2)]];
	endfor 
endfunction

#FUNCTION THAT BUILD LOOP CLOSURE
function edge= calculateEdge(index_from, index_to, pose_from, pose_to, landsmark_positions_from, landsmark_positions_to_b, w)
	transform= inv(v2t(pose_from))*v2t(pose_to);
	landsmark_positions_to= movePoints(transform, landsmark_positions_to_b);
	subplot(3,1,1);
	plot(landsmark_positions_to(:,1), landsmark_positions_to(:,2), 'b.', landsmark_positions_from(:,1), landsmark_positions_from(:,2), 'g.');
	pause(0.1)
	global gating_tau; 	global gamma_threshold;
	[x_visualize,y_visualize, P, Z] = data_association(landsmark_positions_from, landsmark_positions_to,gating_tau, gamma_threshold, w);
	if size(P,2)>20
		subplot(3,1,2);
		plot (x_visualize,y_visualize, 'c', P(1,:),P(2,:),'g.',Z(1,:),Z(2,:),'b.');
		pause(0.1)

		global num_iterations;	 global DAMPING;	   global kernel_threshold;	global X_guess;
		[X_current, chi_evolution, inliers_evolution]=doICPManifold(X_guess, P, Z, num_iterations, DAMPING, kernel_threshold);
		edge= t2v(inv(X_current)*transform); #ICP RETURN POSE OF WORLD IN ROBOT FRAME SO WE NEED TO INVERT IT
		
	else 
		edge= [1000; 1000;1000]; 	#RISULTATO FINTO SE NON HO ABBASTANZA ASSOCIAZIONI
	endif
endfunction

MIN_ANGLE=-1.5708;
MAX_ANGLE=1.5708;
ANGLE_INCREMENT=0.00436332;
MAX_RANGE=30;
global num_iterations=100;
global DAMPING=0;#1e-3;
global kernel_threshold=1e-3;

global gating_tau      = 1;#1.5;
global gamma_threshold = 1e-3;#1e-2;#1e-3;

position_x=0;
position_y=0;
orientation=0;

pose=[0; 0; 0];

global X_guess=eye(3);

trajectory_x=[];
trajectory_y=[];
trajectory_theta=[];

gT_trajectory_x=[];
gT_trajectory_y=[];
gT_trajectory_theta=[];

step=5;

er=0;
number_iteration=0


dataOut = fopen('graph.g2o', 'w+');

vertex_counter=0;

nodes=[];
for i=1:step:size(observations, 1)
	i
	if (false && vertex_counter-1>=0)
		landsmark_positions=laser2map(observations(i,:),MIN_ANGLE, MAX_ANGLE, ANGLE_INCREMENT, MAX_RANGE);
		landsmark_positions_0=laser2map(observations(max(1,i-step),:),MIN_ANGLE, MAX_ANGLE, ANGLE_INCREMENT, MAX_RANGE);
		[x_visualize,y_visualize, P, Z] = data_association(landsmark_positions_0, landsmark_positions,gating_tau, gamma_threshold,10);
		[X_current, chi_evolution, inliers_evolution]=doICPManifold(X_guess, P, Z, num_iterations, DAMPING, kernel_threshold);
		[X_current_b, chi_evolution_b, inliers_evolution_b]=doICPManifold(X_guess, Z, P, num_iterations, DAMPING, kernel_threshold);
	#DAL MOMENTO CHE ICP 'SOTTOVALUTA' IL RISULTATO LO HO UN PO ACRESCIUTO AGGIUNGENDO IL RISULTATO OTTENUTO FACENDO DA Z A L		
		X_current=v2t(t2v(X_current)+ t2v(inv(X_current_b)));
		pose= t2v(v2t(pose)*inv(X_current));
	endif

	#NODES
	fprintf(dataOut, 'VERTEX_SE2 ');
	gT_pose= odometry(i,:);
	#if i>1700
	#	gT_pose(2)-=0.01;
	#endif
	node_info=[vertex_counter, gT_pose, i];#pose', i];		#nodes info è composto da idnode, x,y, theta, index_odometry
	nodes= [nodes; node_info];
	fprintf(dataOut, mat2str(node_info(1:4))(2:end-1));
	fprintf(dataOut, '\n');
	
	#LASERS
	fprintf(dataOut, 'ROBOTLASER1 0 -1.5708 3.14159 0.00436332 30 0.1 0 721 ');
	fprintf(dataOut, mat2str(observations(i,:))(2:end-1));
	fprintf(dataOut, ' 0 -45.0844 -10.9803 1.54375 -45.0844 -10.9803 1.54375 0.000000 0.000000 0.000000 0.000000 0.000000 -1.000000 hostname -1.000000\n');
	
	
	#EDGES SUCCESSIVI
	if( vertex_counter-1>=0)
		gT_pose_precedent_step=odometry(max(1,i-step),:);
		X_guess=(inv(v2t(gT_pose))*v2t(gT_pose_precedent_step));
		x= t2v(inv(X_guess));
		x_inv= t2v(X_guess);
		fprintf(dataOut, 'EDGE_SE2 ');
		fprintf(dataOut, cstrcat(num2str(vertex_counter-1)," ",num2str(vertex_counter)," "));
		fprintf(dataOut, mat2str(x')(2:end-1));
		fprintf(dataOut, ' 10000 0 0 10000 0 1000\n'); #questi valori si rifanno all'altro dataset...
		fprintf(dataOut, 'EDGE_SE2 ');
		fprintf(dataOut, cstrcat(num2str(vertex_counter)," ",num2str(vertex_counter-1)," "));
		fprintf(dataOut, mat2str(x_inv')(2:end-1));
		fprintf(dataOut, ' 10000 0 0 10000 0 1000\n'); #questi valori si rifanno all'altro dataset... 
	endif
	
	#EDGES LOOPCLOUSER
	if(  mod(vertex_counter,5)==0 && vertex_counter-1>=0)# && i>1700)
		# lookfor nodes close to actual node but not close in index
		for (index= 1:size(nodes,1))
			n= nodes(index, :);
			if abs(n(1)-vertex_counter)>100
				if sqrt((n(2)-node_info(2))^2+(n(3)-node_info(3))^2)<0.5
					gT_pose_precedent_step=odometry(max(1,n(5)),:);
					subplot(3,1,3);
					plot(odometry(:,1), odometry(:,2), '.b',"markersize", 0.7,gT_pose_precedent_step(1),gT_pose_precedent_step(2), 'og', gT_pose(1), gT_pose(2), 'or' );
					landsmark_positions_to=laser2map(observations(i,:),MIN_ANGLE, MAX_ANGLE, ANGLE_INCREMENT, MAX_RANGE);
					landsmark_positions_from=laser2map(observations(max(1,n(5)),:),MIN_ANGLE, MAX_ANGLE, ANGLE_INCREMENT, MAX_RANGE);
					#X_guess=(inv(v2t(gT_pose))*v2t(gT_pose_precedent_step));
					#x= t2v(inv(X_guess));
					#x_inv= t2v(X_guess);
					edge= calculateEdge(n(1), vertex_counter, gT_pose_precedent_step, gT_pose, landsmark_positions_from, landsmark_positions_to, -1);
					if edge(1)<900
						fprintf(dataOut, 'EDGE_SE2 ');
						fprintf(dataOut, cstrcat(num2str(n(1))," ",num2str(vertex_counter)," "));
						fprintf(dataOut, mat2str(edge')(2:end-1));
						fprintf(dataOut, ' 10000 0 0 10000 0 1000\n'); #questi valori si rifanno all'altro dataset...
					#edge= calculateEdge(vertex_counter, n(1),  gT_pose, gT_pose_precedent_step, landsmark_positions_to, landsmark_positions_from, -1);
						edge= t2v(inv(v2t(edge)));
						fprintf(dataOut, 'EDGE_SE2 ');
						fprintf(dataOut, cstrcat(num2str(vertex_counter)," ",num2str(n(1))," "));
						fprintf(dataOut, mat2str(edge')(2:end-1));
						fprintf(dataOut, ' 10000 0 0 10000 0 1000\n'); #questi valori si rifanno all'altro dataset... 
					else
						printf("DISCARTED\n");	
						break;
					endif
				endif
			endif
		endfor
	endif

	vertex_counter+=1;
	
	end

