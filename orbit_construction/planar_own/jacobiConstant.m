function jc = jacobiConstant(x, ~, mu)

    m2 = mu;
    r13= sqrt((x(1)+m2)^2+x(2)^2+x(3)^2);                                   % r: distance to m1, LARGER MASS
    r23= sqrt((x(1)-1+mu)^2+x(2)^2+x(3)^2);                                 % R: distance to m2, smaller mass

    u = 2*(1-mu)/r13 + 2*mu/r23+x(1)^2+x(2)^2;
    v = x(4)^2+x(5)^2+x(6)^2;
    jc = u - v;

end