γ_1 = 0.464835

α_1 = log(0.4009298)
μ = log(3.986005e14) # m^3/s^2 standard gravitational parameter
R = log(6378e3) # m
η_A = log(0.29)
Q = log(1367) # W/m^2
EN = log(40)
D_r = log(5.3) # m
L = log(9.772372209558107)
k = log(1.38064852e-23) # J/K Boltzman constant R/Na
T_s = log(135) # K
B = log(8) # bit
N = log(2e3) # pixel width
η = log(0.55)
λ_v = log(500e-9) # m
f = log(2.2e9) # Hz
c = log(2.998e8) # m/s
ρ_A = log(10) # kg/m^2
ρ_p = log(100) # kg/m^1.5
ρ_T = log(2) # kg/m^1.5
ρ_P = log(500e3) # kg*m
ρ_b = log(0.1e-3) # kg/J
P_l = log(5) # W
X_r = log(20) # m
m_c = log(0.2) # kg required weight
η_S = log(0.2) # 20% of total mass is structural mass

λ_c = c - f
G_r = η + 2*(log(π) + D_r - λ_c)

h_min = log(1e5) # orbit altitude
h_max = log(2e6)

d_min = log(0.1) # daylight fraction of orbit
d_max = log(0.9)

e_min = log(0.1) # eclipse fraction of orbit
e_max = log(0.9)

exp_RhR_min = 1/(exp(h_min - R) + 1)
exp_RhR_max = 1/(exp(h_max - R) + 1)

g_min = log(acos(exp_RhR_min)/pi) # ground station viewing fraction of orbit
g_max = log(acos(exp_RhR_max)/pi)
;
