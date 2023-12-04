using LinearAlgebra

function Matrixn(n)

    A=Matrix{Float64}(I,n,n)
    for i = 1:n-1
        A[i,i] = sqrt(i)
        A[i,i+1] = (1+i)^(0.03)
    end
    A[n,n] = sqrt(n)
    return A
end

function QRmethod(A::Matrix)

n=size(A)[1]
H1=A'*A
eigval=eigvals(H1)
s=sqrt.(eigval)

s1=[]    
    for i=1:1e6
        Q,R=qr(H1-I)
        H2=R*Q+I
        H1=H2
        s1=sqrt.(diag(H1))
        err1=(sort(s)[n]-sort(s1)[n])/s[n]
        err2=(sort(s)[n]-sort(s1)[n])/s[n]
        if err1<1e-6 && err2<1e-6
            break
        end
    end
    return s1
end
                             
x=QRmethod(Matrixn(25))
print(x)                    

a1=Matrixn(10)
a2=Matrixn(20)
a3=Matrixn(40)
a4=Matrixn(80)
a5=Matrixn(160)
a6=Matrixn(320)
a7=Matrixn(640)
println(@time QRmethod(a1))
println(@time QRmethod(a2))
println(@time QRmethod(a3))
println(@time QRmethod(a4))
println(@time QRmethod(a5))
println(@time QRmethod(a6))
println(@time QRmethod(a7))

using Plots
x=rand(7)
y=rand(7)

for i=1:7
    x[i]=10*2^(i-1)    
end

using Plots
x=rand(7)
y=rand(7)

for i=1:7
    x[i]=10*2^(i-1)    
end

y[1]=(0.001604)^(1/3)
y[2]=(0.012711)^(1/3)
y[3]=(0.047803)^(1/3)
y[4]=(0.474650)^(1/3)
y[5]=(3.755898)^(1/3)
y[6]=(42.712512)^(1/3)
y[7]=(460.706094)^(1/3)

plot(x,y, label="x versus y^(1/3)", legend=:topleft)                                                                                                                                                                               