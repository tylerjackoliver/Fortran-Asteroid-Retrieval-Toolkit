function [sn,un,cn,Ws,Wu,Wc] = eigGet(A,discrete) ;

%        [sn,un,cn,Ws,Wu,Wc] = eigGet(A,discrete) ;
%
% Compute the eigenvalues and eigenvectors of the matrix A spanning the 
% three local subspaces <Es,Eu,Ec> where A is MxM and s+u+c=M, which
% locally approximate the invariant manifolds <Ws,Wu,Wc>
%
% if discrete   time system, discrete = 1
% if continuous time system, discrete = 0
%
% Shane Ross (revised 2.19.04)

sn=[];
un=[];
cn=[];
Ws=[];
Wu=[];
Wc=[];

% arbitrary small displacement for use in discrete case 
delta = 1.e-4;  % <== maybe needs tweeking?
[M,dum]=size(A);
[V,D]  =eig(A); % obtain eigenvectors (V) and eigenvalues (D) of matrix A
V = cleanUpMatrix(V) ;
D = cleanUpMatrix(D) ;
s=0;u=0;c=0;
for k = 1:M            
	if     discrete == 1             % discrete time system
       		if     abs(D(k,k))-1 < -delta 
	     		s=s+1;
			sn(s)   = D(k,k);
			jj=1; while abs(V(jj,k)) == 0, jj=jj+1; end
               	 	Ws(:,s) = V(:,k)./V(jj,k);% stable   (s dimensional)
			Ws(1:M,s) = RemoveInfinitesimals(Ws(:,s));
 		elseif abs(D(k,k))-1 >  delta 
			u=u+1;
			un(u)   = D(k,k);
			jj=1; while abs(V(jj,k)) == 0, jj=jj+1; end
                	Wu(:,u) = V(:,k)./V(jj,k);% unstable (u dimensional)
			Wu(1:M,u) = RemoveInfinitesimals(Wu(:,u));
        	else
			c=c+1;
			cn(c)   = D(k,k);
			jj=1; while abs(V(jj,k)) == 0, jj=jj+1; end
                	Wc(:,c) = V(:,k)./V(jj,k);% center   (c dimensional)
			Wc(1:M,c) = RemoveInfinitesimals(Wc(:,c));
        	end
	elseif discrete == 0	         % continuous time system
        	if     real(D(k,k)) < 0
	     		s=s+1;
			sn(s)   = D(k,k);
			jj=1; while abs(V(jj,k)) == 0, jj=jj+1; end
               	 	Ws(:,s) = V(:,k)./V(jj,k);% stable   (s dimensional)
			Ws(1:M,s) = RemoveInfinitesimals(Ws(:,s));
        	elseif real(D(k,k)) > 0
			u=u+1;
			un(u)   = D(k,k);
			jj=1; while abs(V(jj,k)) == 0, jj=jj+1; end
                	Wu(:,u) = V(:,k)./V(jj,k);% unstable (u dimensional)
			Wu(1:M,u) = RemoveInfinitesimals(Wu(:,u));
        	else
			c=c+1;
			cn(c)   = D(k,k);
			jj=1; while abs(V(jj,k)) == 0, jj=jj+1; end
                	Wc(:,c) = V(:,k)./V(jj,k);% center   (c dimensional)
			Wc(1:M,c) = RemoveInfinitesimals(Wc(:,c));
                end
	end
end

function A=RemoveInfinitesimals(A);

% Remove all entries with absolute value smaller than TOL
TOL=1.e-14;

for k=1:length(A),
	if abs(real(A(k))) < TOL,
		A(k)=A(k) - real(A(k));
	end
	if abs(imag(A(k))) < TOL,
		A(k)=A(k) - i*imag(A(k));
	end
end
