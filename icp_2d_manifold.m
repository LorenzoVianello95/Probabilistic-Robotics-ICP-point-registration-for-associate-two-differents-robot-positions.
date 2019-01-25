source "common_stuff.m"

function [e,J]=errorAndJacobianManifold(X,p,z)
 t=X(1:2,3);
 R=X(1:2,1:2);
 z_hat=R*p+t;
 e=z_hat-z;

 J=zeros(2,3);
 J(1:2,1:2)=eye(2);
 J(1:2,3)=[-z_hat(2),
	   z_hat(1)]';
endfunction;

function [chi,inl, X]=icp2dManifold(X,P,Z)
  chi=0; %cumulative chi2
  inl=0;
  threshold=1;#1e-5;#todo da abbassare
  H=zeros(3,3);  b=zeros(3,1); %accumulators for H and b
  for(i=1:size(P,2))
     p=P(:,i); z=Z(:,i); % fetch point and measurement
     [e,J]=errorAndJacobianManifold(X,p,z); %compute e and J for the point
     if(e'*e>threshold)
	e*=sqrt(threshold/(e'*e));
     else
	inl+=1;
     endif
     H+=J'*J;            %assemble H and B
     b+=J'*e;
     chi+=e'*e;          %update cumulative error
  endfor
  dx=-H\b;               %solve the linear system
  X=v2t(dx)*X;
  inl/=size(P,2);
endfunction

function [X_current,chi_evolution, inliers_evolution]=testICP2DManifold(X_current, P, Z, num_iterations)
  chi_evolution=zeros(num_iterations,1);
  inliers_evolution=zeros(num_iterations,1);
  #X_current=eye(3);
  for (i=1:num_iterations)
      [chi, inl, X_current]=icp2dManifold(X_current,P,Z);
      chi_evolution(i)=chi;
      inliers_evolution(i)=inl;
  endfor;
endfunction

function [X, chi_stats, num_inliers]=doICPManifold(X_guess, P, Z, num_iterations, damping, kernel_threshold)
  X_guess=eye(3);
  X=X_guess;
  chi_stats=zeros(1,num_iterations);
  num_inliers=zeros(1,num_iterations);
  for (iteration=1:num_iterations)
    H=zeros(3,3);
    b=zeros(3,1);
    chi_stats(iteration)=0;
    for (i=1:size(P,2))
      [e,J] = errorAndJacobianManifold(X, P(:,i), Z(:,i));
      chi=e'*e;
      if (chi>kernel_threshold)
	 e*=sqrt(kernel_threshold/chi);
	 chi=kernel_threshold;
      else
	num_inliers(iteration)++;
      endif;
      chi_stats(iteration)+=chi;
      H+=J'*J;
      b+=J'*e;
    endfor
    H+=eye(3)*damping;
    dx=-H\b;
    X=v2t(dx)*X;
  endfor
endfunction

