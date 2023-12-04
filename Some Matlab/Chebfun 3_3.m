L = chebop(-1, 1);
n=10;
L.op = @(x,y) (1-x^2)*diff(y,2)-2*x*diff(y)+n*(n+1)*y;
L.lbc = (-1)^n; L.rbc =1; u=L\0;

L = chebop(-1, 1);
n=40;
L.op = @(x,y) (1-x^2)*diff(y,2)-2*x*diff(y)+n*(n+1)*y;
L.lbc = (-1)^n; L.rbc =1; v=L\0;

subplot(1,2,1)
plot(u);
title('BVP Solution P_{40}')

subplot(1,2,2)
plot(legpoly(40));
title('Legendre Polynomial P_{40}')

a=u'*u;
b=v'*v;
c=u'*v;