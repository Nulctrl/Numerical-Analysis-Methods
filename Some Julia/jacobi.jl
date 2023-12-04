using LinearAlgebra

N=20
A=Matrix{Float64}(I,N,N)
for a=1:N
    for b=1:N
          A[a,b]=sqrt(a^2+b^2)
    end
end
A1=A

R1=Matrix{Float64}(I,N,N)

for i=1:100000

B=copy(A1)
for f=1:N
B[f,f]=0.0
end

p=argmax(broadcast(abs,B))[1]
q=argmax(broadcast(abs,B))[2]
fi=(1/2)*atan(2A1[p,q],A1[q,q]-A1[p,p])
R=Matrix{Float64}(I,N,N)
R[p,p]=cos(fi)
R[p,q]=sin(fi)
R[q,p]=-sin(fi)
R[q,q]=cos(fi)

sum=0
for n=1:N
for m=1:N
if(m!=n)
sum+=(A1[n,m])^2
end
end
end
sum=sqrt(sum)

if(sum<1e-8)
print(i)
println("times")
println(sum)
break
end

global A1=transpose(R)*A1*R
global R1*=R
end

s=[]
for z=1:N
push!(s,A1[z,z])
end
    
u=R1

d=sort(s)
for g=1:5
println(d[N-g+1])

end
println(u)