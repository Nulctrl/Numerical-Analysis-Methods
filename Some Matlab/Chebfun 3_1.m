ep = [0.5];
k=8;
x = chebfun('x');
F = (abs(x)<ep)/(2*ep); 
L = chebop(@(x,y) diff(y,2),[-1 1],0,0);
M = chebop(@(x,y) -F(x)*y);
[V,D] = eigs(L,M,k);
for i=1:8
    if V{i}(0.75)<0
        VV{i} = V{i}*(-1);   
    else
        VV{i} = V{i}; 
    end
end

subplot(3,3,1)
plot(VV{1})
title('Eigenfunction 1')

subplot(3,3,2)
plot(VV{2})
title('Eigenfunction 2')

subplot(3,3,3)
plot(VV{3})
title('Eigenfunction 3')

subplot(3,3,4)
plot(VV{4})
title('Eigenfunction 4')

subplot(3,3,5)
plot(VV{5})
title('Eigenfunction 5')

subplot(3,3,6)
plot(VV{6})
title('Eigenfunction 6')

subplot(3,3,7)
plot(VV{7})
title('Eigenfunction 7')

subplot(3,3,8)
plot(VV{8})
title('Eigenfunction 8')