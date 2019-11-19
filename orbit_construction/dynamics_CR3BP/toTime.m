function I = toTime(t,tf)
%        I = toTime(t,tf);
% find the index I of the element in matrix t which 
% has an absolute value	closest to tf

t=abs(t);

for k=1:length(tf)
	T=abs( tf(k)*ones(size(t)) - t );
	[dum,II]=min(T);
	I(k)=II;
end
