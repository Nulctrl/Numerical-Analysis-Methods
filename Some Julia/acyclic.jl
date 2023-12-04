using FourierTools
f=[1/32]
for i=1:31
push!(f,1/32)
end

for i=1:351
push!(f,0)
end

g=[0.0]
for i=1:191
push!(g,0)
end

for i=1:64
push!(g,1/64)
end

for i=1:127
push!(g,0)
end

a=conv(f,g)

print(a)

x=rand(383)
for i=1:383
    x[i]=i
end

plot(x,a,label="Acyclic Convolution f*g",legend=:topleft )