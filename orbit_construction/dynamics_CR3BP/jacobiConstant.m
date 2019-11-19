function jc = jacobiConstant(x, t, massParameter)
    m1 = 1-massParameter;
    m2 = massParameter;
    r13= sqrt((x(1)+m2)^2+x(2)^2+x(3)^2);      % r: distance to m1, LARGER MASS
    r23= sqrt((x(1)-1+massParameter)^2+x(2)^2+x(3)^2);      % R: distance to m2, smaller mass
%     omega = 0.5*(x(1)^2+x(2)^2) + (1-massParameter)/r13 ...
%         + massParameter/r23;
    u = 2*(1-massParameter)/r13 + 2*massParameter/r23+x(1)^2+x(2)^2;
    v = x(4)^2+x(5)^2+x(6)^2;
    jc = u - v;
%     fileId = fopen('jacobi.txt', 'a+');r
%     fprintf(fileId,'%2.8f \t %2.8f\n', t,r jc);
%     fclose(fileId);
end