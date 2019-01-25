#PRENDE I VALORI DAL FILE OSSIA LE LETTURE LASER

function observations= loadData(filepath, fast)
	ob=[];
	ang=[];
	%open the file
	fid = fopen(filepath, 'r');
	counter=1;
	
	while true
		%get current line
		c_line = fgetl(fid);
		
		%stop if EOF
		if c_line == -1
			break;
		end
	
		if counter> 50 && fast
			break;
		end
		counter+=1

		elements = strsplit(c_line,' ');
		
		ob=[ob;extractMeasurement(elements)];
	end
	observations=ob;
	anglesInfo=ang;

endfunction

#LEGGE LE ODOMETRIE DAL FILE 
function trajectory= loadOdometry(filepath, fast)
	pos=[];
	fid = fopen(filepath, 'r');

	counter=1;
	while true
		%get current line
		c_line = fgetl(fid);
		
		%stop if EOF
		if c_line == -1
			break;
		end
	
		if counter> 50 && fast
			break;
		end
		counter+=1

		%Split the line using space as separator
		elements = strsplit(c_line,' ');
		
		pos=[pos;extractOdometry(elements)];
	end
	trajectory= pos;
endfunction


function out = extractMeasurement(elements)
  meas=[];
  for i=22:742
	meas(end+1)=str2double(elements{i});
	end
  out =meas;
endfunction

function out = extractOdometry(elements)
	x= str2double(elements{13});
	y= str2double(elements{14});
	#qz is expressed in quaternion unity
	q_z= str2double(elements{18});
		
	#TRANSFORM QUATERNION IN R,P,Y
	q_w= sqrt(1-q_z^2);
	sin_yaw= 2*q_w*q_z;
	cos_yaw= 1-2*q_z^2;
	yaw= atan2(sin_yaw, cos_yaw);
	
  	od=[x, y, yaw];
  	out =od;
endfunction

function out = extractAngleInformations(elements)
  min_angle = str2double(elements{16});
  max_angle = str2double(elements{17});
  angle_increment = str2double(elements{18});
  out =[min_angle,max_angle,angle_increment];
endfunction


