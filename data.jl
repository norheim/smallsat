# γ_1 = 0.464835
# α_1 = log(0.4009298)


μ = log(3.986005e14) # m^3/s^2 standard gravitational parameter
R = log(6378e3) # m
#η_A = log(0.29)
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
#ρ_A = log(10) # kg/m^2
#ρ_p = log(100) # kg/m^1.5
#ρ_T = log(2) # kg/m^1.5
ρ_P = log(500e3) # kg*m
#ρ_b = log(0.1e-3) # kg/J
P_l = log(5) # W
X_r = log(20) # m
m_c = log(0.2) # kg required weight
m_max = log(10)
η_S = log(0.2) # 20% of total mass is structural mass
C_D = log(2.2) # coefficient of drag
G = log(9.81) # m/s^2 gravitational constant
I_sp = log(536) # s impulse specific constant
w_W = log(1000) # radians/s for maximum angular velocity of reaction wheel
c_W = log(1000) # constant to maque gravity gradient significant enough to count
rho_M = log(1) # density for mass of magnetorquers
rho_W = log(100) # density for mass for reaction wheels

λ_c = c - f
G_r = η + 2*(log(π) + D_r - λ_c)

exp_Lt_min = 3 # years lifetime

exp_Ln_min = 0.5
exp_Ln_max = 25

exp_Lp_min = 0.01
exp_Lp_max = 10.0

h_min = log(1e5) # orbit altitude
h_max = log(2e6)

d_min = log(0.1) # daylight fraction of orbit
d_max = log(0.9)

#e_min = log(0.1) # eclipse fraction of orbit
#e_max = log(0.9)

exp_RhR_min = 1/(exp(h_min - R) + 1)
exp_RhR_max = 1/(exp(h_max - R) + 1)

g_min = log(acos(exp_RhR_min)/pi) # ground station viewing fraction of orbit
g_max = log(acos(exp_RhR_max)/pi)


n_T = 2 # transmitter catalog
D_Ti = log.([0.07, 0.14])
m_Ti = log.([0.053, 0.300])
P_Ti = log.([10, 10])
#G_Ti = log.([10, 44.7])

n_b = 5 # battery catalog
E_bi = log.([138600, 144000, 144000, 165600, 1607040])
m_bi = log.([0.270, 0.310, 0.355, 0.710, 3.95])

n_p = 2 # payload catalog
D_pi = log.([0.04, 0.1])
m_pi = log.([0.080, 2])

n_A = 2 # solar panel catalog
ρ_Ai = log.([5.35, 1.05])
η_Ai = log.([0.305, 0.268])


# piecewise linear table
h_ρi = log(1e3) + log.([100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1250, 1500])
ρi = log.([4.79e-07, 1.81e-09, 2.53e-10, 6.24e-11, 1.95e-11, 6.98e-12, 2.72e-12, 1.13e-12, 4.89e-13, 2.21e-13, 1.04e-13, 5.15e-14, 2.72e-14, 1.55e-14, 9.63e-15, 6.47e-15, 4.66e-15, 3.54e-15, 2.79e-15, 1.11e-15, 5.21e-16])
Hi = log(1e3) + log.([5.9, 25.5, 37.5, 44.8, 50.3, 54.8, 58.2, 61.3, 64.5, 68.7, 74.8, 84.4, 99.3, 121, 151, 188, 226, 263, 296, 408, 516])
;
