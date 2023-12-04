L = chebop(0, 20);
L.op = @(t,y) diff(y) - log(1-abs(y))*y/2-3*cos(t);
L.lbc = 0; u = L\0;
du = chebfun(diff(u));
ddu = chebfun(diff(du));
plot (ddu), grid on
y=chebfun(u,[1,20]);
dy=chebfun(du,[1,20]);
ddy=chebfun(ddu,[1,20]);
c=max(u);
xlabel('t') 
ylabel('diff^2(y)') 
