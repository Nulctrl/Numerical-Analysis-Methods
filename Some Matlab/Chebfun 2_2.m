L = chebop(0, 40);
p=-1.5;
L.op = @(t,y) diff(y,2) +y*abs(y)^(p-1);
L.lbc = [1;0.5*1i]; u = L\0;

plot (u), grid on
axis([-1 1 -1 1]),axis square.
xlabel('t') 
ylabel(y) 
