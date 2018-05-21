using PyPlot
include("data.jl")
h = linspace(400,2000);
a = exp(R) +h*1000;
T = 2*π*sqrt.(a.^3/exp(μ));
c_W = log(1)
T_g = 3*exp(μ)*exp(c_W)./(2*a.^3);
T_m = B;
L = T_g.*T;
ρ_M =11e4
ρ_P2 = 9/0.1

h[26]
T_g[26]
hm = T.*T_g/4*sqrt(2)/2;
hm[11]
pulses = 3*365*5
F = ρ_P2*hm;
mass = F*pulses/(exp(I_sp)*exp(G));
(ρ_P2*T.*T_g/4*sqrt(2)/2*3*365*5/(exp(I_sp)*exp(G)))[26]
ll = log(ρ_P2) + log.(T) + log.(T_g) + log(1/4*sqrt(2)/2*3*365) + log(exp_Lt_min) - I_sp - G;
exp.(ll)[26]
mass[26]

Isc = 1/12*5000*(0.1*0.1*0.3*(0.3^2-0.1^2))
D_W = 0.021
D_W =(L/(0.5*exp(rho_W)*exp(w_W))).^(1/4);
m_W = exp(rho_W)*D_W.^2;
#0.5*m_W*D_W^2*exp(w_W) #torque
B = 2*exp(M_B)./a.^3;
T_g[1]
B[1]
rho_M = log(0.005)
m_M = ρ_M*T_g;
