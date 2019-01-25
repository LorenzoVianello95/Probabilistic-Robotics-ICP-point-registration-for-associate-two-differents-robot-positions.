source "./data_loader.m"
source "./laser2map.m"
source "./data_association.m"
#source "./ICP.m"
#source "icp_2d.m"
source "common_stuff.m"
source "icp_2d_manifold.m"

fast= true; # se true fa solo le prime 50 iterazioni

global observations= loadData('laser_scans.txt', fast);
odometry= loadOdometry('out.txt', fast);
printf("Observations extracted")

MIN_ANGLE=-1.5708;
MAX_ANGLE=1.5708;
ANGLE_INCREMENT=0.00436332;
MAX_RANGE=30;
num_iterations=100;
DAMPING=0;
kernel_threshold=1e-2;

gating_tau      = 1.5;
gamma_threshold = 1e-2;#1e-3;

position_x=0;
position_y=0;
orientation=0;

pose=[0; 0; 0];

X_guess=eye(3);

trajectory_x=[];
trajectory_y=[];
trajectory_theta=[];

gT_trajectory_x=[];
gT_trajectory_y=[];
gT_trajectory_theta=[];

step=5;

er=0;
number_iteration=0

localization=true;
if localization
	for i=1:step:size(observations, 1)
		i
		landsmark_positions=laser2map(observations(i,:),MIN_ANGLE, MAX_ANGLE, ANGLE_INCREMENT, MAX_RANGE);
		landsmark_positions_0=laser2map(observations(max(1,i-step),:),MIN_ANGLE, MAX_ANGLE, ANGLE_INCREMENT, MAX_RANGE);
		[x_visualize,y_visualize, P, Z] = data_association(landsmark_positions_0, landsmark_positions,gating_tau, gamma_threshold, 10);
		#plot(landsmark_positions(:,1), landsmark_positions(:,2), 'r.');
		#plot(x_visualize,y_visualize, 'b', P(1,:),P(2,:),'g.',Z(1,:),Z(2,:),'r.');

		gT_pose= odometry(i,:);
		gT_pose_precedent_step=odometry(max(1,i-step),:);
		
		X_guess=inv(v2t(gT_pose))*v2t(gT_pose_precedent_step)

		[X_current, chi_evolution, inliers_evolution]=doICPManifold(X_guess, P, Z, num_iterations, DAMPING, kernel_threshold);
		[X_current_b, chi_evolution_b, inliers_evolution_b]=doICPManifold(X_guess, Z, P, num_iterations, DAMPING, kernel_threshold);
	#DAL MOMENTO CHE ICP 'SOTTOVALUTA' IL RISULTATO LO HO UN PO ACRESCIUTO AGGIUNGENDO IL RISULTATO OTTENUTO FACENDO DA Z A L		
		X_current=v2t(t2v(X_current)+ t2v(inv(X_current_b)));
	
		er+= norm(t2v(X_guess)-t2v(X_current));
		number_iteration+=1;
		printf("mean error since now %d \n", er/number_iteration)
		
		pose= t2v(v2t(pose)*inv(X_current));	#ICP RETURN POSE OF WORLD IN ROBOT FRAME SO WE NEED TO INVERT IT
				
		trajectory_x=[trajectory_x;pose(1)];
		trajectory_y=[trajectory_y;pose(2)];
		trajectory_theta=[trajectory_theta; pose(3)];

		gT_trajectory_x=[gT_trajectory_x;gT_pose(1)];
		gT_trajectory_y=[gT_trajectory_y;gT_pose(2)];
		gT_trajectory_theta=[gT_trajectory_theta; gT_pose(3)];

		if true
			subplot (3, 3, 1)
			plot (trajectory_x,trajectory_y,'o');
			title ("xy trajectory");
			subplot (3, 3, 2)
			plot (gT_trajectory_x,gT_trajectory_y,'o');
			title ("xy trajectory gT");
			subplot(3,3,3)
			plot (x_visualize,y_visualize, 'b', P(1,:),P(2,:),'g.',Z(1,:),Z(2,:),'r.');
			title("data association");

			subplot (3, 3, 4)
			plot (trajectory_theta,'o');
			title("trajectory theta");
			subplot (3, 3, 5)
			plot (trajectory_x,'o');
			title("trajextory x");
			subplot (3, 3, 6)
			plot (trajectory_y,'o');
			title("trajectory y");

			subplot (3, 3, 7)
			plot (gT_trajectory_theta,'or');
			title("trajectory theta gt");
			subplot (3, 3, 8)
			plot (gT_trajectory_x,'or');
			title("trajectory x gt");
			subplot (3, 3, 9)
			plot (gT_trajectory_y,'or');
			title("trajectory y gt");
		endif

		pause(0.001);
		end
	pause	
endif
