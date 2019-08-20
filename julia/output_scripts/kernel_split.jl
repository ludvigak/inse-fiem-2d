push!(LOAD_PATH, string(pwd(),"/src"))
using ModifiedStokesSolver
using PyPlot
using LaTeXStrings

n = 500
zmax = 10
z = linspace(1e-3, zmax, n)
T1S = zeros(n)
T2S = zeros(n)
T3S = zeros(n)
T1L = zeros(n)
T2L = zeros(n)
T3L = zeros(n)


for i=1:n
    T1S[i], T1L[i], T2S[i], T2L[i], T3S[i], T3L[i] = ModifiedStokesSolver.Ti_split(z[i])
end

close("all")
fs = (5,3)
figure(1, figsize=fs)
clf()
semilogy(z, T1S,  label=L"|\mathcal{T}_1{ }^S(z)|")
semilogy(z, -T1L, label=L"|\mathcal{T}_1{ }^L(z)|")
semilogy(z, -T2S, label=L"|\mathcal{T}_2{ }^S(z)|")
semilogy(z, T2L,  label=L"|\mathcal{T}_2{ }^L(z)|")
semilogy(z, -T3S, label=L"|\mathcal{T}_3{ }^S(z)|")
semilogy(z, T3L,  label=L"|\mathcal{T}_3{ }^L(z)|")
xlim(0, zmax)
xlabel(L"z")
grid("on")
legend()
tight_layout(0)


figure(2, figsize=fs)
clf()
z = logspace(-2, 2, n)
T1 = zeros(n)
T2 = zeros(n)
T3 = zeros(n)
for i=1:n
    T1[i], T2[i], T3[i] = ModifiedStokesSolver.Ti(z[i])
end
loglog(z,  abs.(T1), label=L"|\mathcal{T}_1(z)|")
loglog(z,  abs.(T2), label=L"|\mathcal{T}_2(z)|")
loglog(z,  abs.(T3), label=L"|\mathcal{T}_3(z)|")
loglog(z, 1./z.^2, "--", label=L"z^{-2}")
loglog(z, 1./z.^4, "--", label=L"z^{-4}")
loglog(z, 1./z.^6, "--", label=L"z^{-6}")
xlabel(L"z")
legend()
grid()
tight_layout(0)


function write_fig(num, name)
    figure(num)
    savefig(name)
    println("Saved $name")
end

filename = joinpath(dirname(@__DIR__),"figs", "Ti_asymp.pdf")
write_fig(2, filename)
run(`cp $filename ../docs/paper/figs`)

filename = joinpath(dirname(@__DIR__),"figs", "Ti_split_asymp.pdf")
write_fig(1, filename)
run(`cp $filename ../docs/paper/figs`)

