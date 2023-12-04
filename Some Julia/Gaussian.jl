function P(n, x)
if n==0
    return 1
end
if n==1
    return x
end
    a=1
    b=x
    for i=1:n-1
        c=a
        a=b
        b=(2i+1)*b*x/(i+1)-c*i/(i+1)
    end
    return b
end

function Pd(n,x)
if n==0
    return 0
end
if n==1
    return 1
end
a=1
    for i=1:n-1
        a=x*a+(i+1)*P(i,x)
    end
return a
end
        
  
function myroots2(n)
    
roots=[]

for k=1:n
    node=cos((2k-1)*pi/2n)
    push!(roots,node)
end

f(x)=P(n,x)
fp(x)=Pd(n, x)
roots2=Float64[]
    for i in 1:n
        a=roots[i]
        for m in 1:1000
            b=a
            a=a-f(a)/fp(a)
            if abs((a-b))<1e-17
                break
            end
        end
        push!(roots2, a)
    end
    
    return roots2
end

import LinearAlgebra

function weights(n)
A=zeros(n,n)
    for i=1:n
         for j=1:n
            A[i,j]=cos((i-1)*acos(myroots2(n)[j]))
        end
    end

results=[]
push!(results,2.0)
push!(results,0.0)
    for i=2:n-1
        result=(cos((i+1)acos(1))/(i+1)-cos((i-1)acos(1))/(i-1)-cos((i+1)acos(-1))/(i+1)+cos((i-1)acos(-1))/(i-1))/2
        push!(results,result)
    end
weights=A^(-1)*results  
    
return(weights)
    
end

trueval1=10*(1-exp(-2pi))/101
phi1(t)=pi*(exp(-(pi+pi*t)))*sin(10*(pi+pi*t))

answer=0.0
for n=3:100
    weights3=weights(n)
    nodes3=myroots2(n)
    answer=0.0
    for i=1:n
        answer+=weights3[i]*phi1(nodes3[i])
    end
    
    if (abs(answer-trueval1)<1e-8)
        println(answer)
        println(n)
        break
    else
        println(-n)
    end
end