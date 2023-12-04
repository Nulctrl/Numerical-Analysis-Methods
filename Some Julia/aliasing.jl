function q(x)
answer=0
for k=-1e7:1e7
answer+=1/(1+k^2)exp(im*k*x)
end
return answer
end

a=[q(0)]

for i=1:15
push!(a,q(i*pi/8))
end

using FFTW
b=fft(a)/16

answer=0

for j=1:1e6
global answer+=coef(16*j)*2
end

c=[answer+1]

for i=1:15
global answer=0
for j=-1e6:1e6
global answer+=coef(i+16*j)
end
push!(c,answer)
end

print(b-c)
