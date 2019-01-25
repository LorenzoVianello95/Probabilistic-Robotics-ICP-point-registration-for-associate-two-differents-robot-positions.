
%laser lecture are the 721 laser scan
function landsmark_positions= laser2map(laserlecture,min_angle, max_angle, delta_angle, max_range)
	lands_position=[];

	for i= 1:721
		l_depth= laserlecture(1,i);
		l_ang=(i-1)*delta_angle+min_angle;
		land_pos= [l_depth*cos(l_ang) l_depth*sin(l_ang)];
		#DISCART LECTURE TOO FAR
		if l_depth<max_range
			lands_position=[lands_position; land_pos];
		endif
	end
	landsmark_positions=lands_position;

endfunction
