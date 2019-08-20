push!(LOAD_PATH, string(pwd(),"/src"))
using ModifiedStokesSolver
using PyPlot
using LaTeXStrings

n = 500
zmax = 5
z = linspace(1e-8, zmax, n)
T1S = zeros(n)
T2S = zeros(n)
T3S = zeros(n)
T1L = zeros(n)
T2L = zeros(n)
T3L = zeros(n)

for i=1:n
    T1S[i], T2S[i], T3S[i], T1L[i], T2L[i], T3L[i] = ModifiedStokesSolver.Ti_split(z[i])
end

figure(1)
clf()
plot(z, T1S, label="T1S(z)")
plot(z, T2S, label="T2S(z)")
plot(z, T3S, label="T3S(z)")
plot(z, T1L, label="T1L(z)")
plot(z, T2L, label="T2L(z)")
plot(z, T3L, label="T3L(z)")

xlim(0, zmax)
legend()
xlabel("z")
grid("on")
show()


figure(2)
clf()
f = [1.0, 0.0]
alpha = 1.0
x = linspace(-3.0, 3.0, 10000)
u = map( x -> stokeslet([x, 0.0], f, alpha)[1], x)
plot(x, u, label=L"S(\alpha|x-y|)")
axis([-3,3,0,0.5])
xlabel(L"\alpha |x-y|")
legend()
