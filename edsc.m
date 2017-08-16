function [Missrate, C, grp] = edsc(X,s,lambda,affine,outlier,Dim,alpha)
% Efficient Dense Subspace Clustering (EDSC)
% Input: X -- data matrix
%        s -- groundtruth labels by [1 1 2 2 3 3...]
%        lambda -- tradeoff parameter
%        affine -- affine constraint or not, true or false, false by
%                  default
%        outlier -- contain outliers or not, true or false, false by
%                  default
%        Dim -- maximal intrinsic dimension of each subspace, 4 by default
%        alpha -- 4 by default, see section 5 of our paper
% Output: Missrate -- missclassification rate
%         C -- Coefficient matrix, see our paper
%         grp -- cluster result after ncuts
% copyright@ Pan Ji, Australian National University, pan.ji@anu.edu.au
% If you use this code in your research, please cite our paper:
% Ji, P., Salzmann, M., Li, H.: Efficient dense subspace clustering. In:
% IEEE WACV, PP 461-468 (2014)
if(nargin<5)
	outlier = false;
end
if(nargin<4)
	affine = false;
end
if(nargin<3)
	lambda = 100;
	lambda1 = 100;
	lambda2 = 50;
elseif(length(lambda)==2)
	lambda1 = lambda(1);
	lambda2 = lambda(2);
elseif(length(lambda)==1)
	lambda1 = lambda;
	lambda2 = lambda/2;
end

[d,N] = size(X);
K = max(s);

if(~outlier)
	lxtx = lambda*(X'*X);
	if(~affine) %linear
		C = (eye(N)+lxtx)\lxtx;
	else %affine
		% There may be some numerical differences between these two
		% methods. Theoretically, they are equavilent.
		%% a solution got by Lagrange multiplier method:
% 		tic
% 		B = (eye(N)+lxtx)\lxtx;
% 		A = (eye(N)+lxtx)\ones(N,1);
% 		lhand = (A'*A) + lambda*((X*A)'*(X*A)) - 2*(ones(1,N)*A);
% 		rhand = B'*A + lambda*(X'-B'*X')*X*A - (B'*ones(N,1)-ones(N,1));
% 		Y = rhand/lhand;
% 		C = B - A*Y';
% 		t1 = toc;
        %%		
		% a direct implementation	
		tic
		One = ones(1,N);
		V0 = eye(N);
		V = null(One);
		lhand = V'*V+V'*lxtx*V;
		rhand = -V'*V0;
		phi = lhand\rhand;
		C = V0+V*phi;	
		time = toc;
	end
else
	lxtx = lambda1*(X'*X);
	if(~affine) %linear
		% ALM
		%parameters
		epsilon = 1e-8; rho = 1e-6; maxIter = 1e4; eta = 3; max_rho = 1e10;
		%initialize
		C = zeros(N,N);
		D = zeros(d,N);
		E = zeros(d,N);
		Y = zeros(d,N);
		iter = 0;
		while(iter<maxIter)
			iter = iter+1;
			%update D
			lhs = lambda1*eye(N)+rho*(eye(N)-C)*(eye(N)-C');
			rhs = lambda1*(X-E)-Y*(eye(N)-C');
			D = rhs/lhs;
			%update C
			dtd = D'*D;
			C = (eye(N)+rho*dtd)\(D'*Y+rho*dtd);
			%update E
			xmd = X-D;
			E = max(0,xmd - lambda2/lambda1)+min(0,xmd + lambda2/lambda1);
			
			leq1 = D-D*C;			
			stpC = max(max(abs(leq1)));
			if(iter == 1 || mod(iter,50)==0 || stpC<epsilon)
				disp(['iter ' num2str(iter) ',rho=' num2str(rho,'%2.1e') ',stopALM=' num2str(stpC,'%2.3e')]);
			end
			if(stpC<epsilon)
				break;
			else
				Y = Y + rho*leq1;				
				rho = min(max_rho,rho*eta);
			end
		end
	else %affine
		disp(['Fill in the right code!!']);
	end
end
grp = postProC(C,K,Dim,alpha);
grp = bestMap(s,grp);
Missrate = sum(s(:) ~= grp(:)) / length(s);

end










