L=chebop(-4*pi,4*pi); L.op=@(x,y) 0.003/16*diff(y,2)-sin(x)*y;
L.lbc=0; L.rbc=0;
y=L\1; plot(y)
axis([-4*pi,4*pi -40 40])