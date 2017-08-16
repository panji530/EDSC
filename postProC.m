function s = post_proC(C,K,d,alpha)
if(nargin<4)
	alpha = 4;
end
if(nargin<3)
	d = 4;
end
N = length(C);

C = (C + C')/2;
r = d*K+1;
opts.u0 = ones(N,1);
[U,S,~] = svds(C,r,'largest',opts);
S = diag(S);
U = U*diag(sqrt(S));
U = normr(U);
Z = U*U';
L = Z.^alpha;

% ncut
L = (L + L')/2;
grp = ncutW(L,K); % get your ncutW funtion via http://www.cis.upenn.edu/~jshi/software/

s = zeros(N,1);
for i = 1:K
    s = s+grp(:,i)*i;
end

end