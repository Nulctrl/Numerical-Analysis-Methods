import LinearAlgebra
f1(x)=log(x)
f2(x)=x*log(x)
f3(x)=x^2*log(x)
f4(x)=x^3*log(x)
n=320
result1=0
result2=0
result3=0
result4=0


for i=1:n
    global result1+=f1(i/n)
end
global result1/=n
global result1-=f1(1)/2n
global result1=-1-result1
global result1*=n


for i=1:n
   global result2+=f2(i/n)
end
global result2/=n
global result2-=f2(1)/2n
global result2=-0.25-result2
global result2*=n


for i=1:n
   global result3+=f3(i/n)
end
global result3/=n
result3-=f3(1)/2n
result3=-1/9-result3
result3*=n


for i=1:n
   global result4+=f4(i/n)
end
global result4/=n
global result4-=f4(1)/2n
global result4=-0.0625-result4
global result4*=n


A=[f1(1/n) f1(2/n) f1(3/n) f1(4/n);f2(1/n) f2(2/n) f2(3/n) f2(4/n);f3(1/n) f3(2/n) f3(3/n) f3(4/n);f4(1/n) f4(2/n) f4(3/n) f4(4/n)]

b=[result1; result2; result3; result4]

x=A^(-1)*b
#[f1(1/n), f1(2/n), f1(3/n), f1(4/n)]   [c1]   [result1]
#[f2(1/n), f2(2/n), f2(3/n), f2(4/n)] * [c2] = [result2]
#[f3(1/n), f3(2/n), f3(3/n), f3(4/n)] * [c3] = [result3]
#[f4(1/n), f4(2/n), f4(3/n), f4(4/n)]   [c4]   [result4]
println(x)
trueresult1=0
for i=1:n
   global trueresult1+=f1(i/n)
end
global trueresult1/=n
global trueresult1-=f1(1)/2n
for j=1:4
   global trueresult1+=x[j]*f1(j/n)/n
end
print(trueresult1)