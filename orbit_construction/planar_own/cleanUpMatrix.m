function A = cleanUpMatrix(A) ;

%        A = cleanUpMatrix(A) ;
%
% Remove all entries in matrix A with absolute value smaller than TOL
% where TOL is set inside cleanUpMatrix.m

TOL=1.e-14;

for k=1:length(A),
  for l = 1:length(A)
    if abs(real(A(k,l))) < TOL,
      A(k,l)=i*imag(A(k,l)) ;
    end
    if abs(imag(A(k,l))) < TOL,
      A(k,l)=real(A(k,l)) ;
    end
  end
end

