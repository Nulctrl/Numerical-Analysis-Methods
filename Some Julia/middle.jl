using LinearAlgebra
A=[1 0 1/2; 0 3 1/4; 1/2 1/4 6]
x=[1,1,1]
y=inv((A-3.01I))*x
w=y/norm(y)

for i=1:100
w=inv((A-3.01I))*w/norm(inv((A-3.01I))*w)
end

y=1/w+s

print(y)