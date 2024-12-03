# Imports
import numpy as np
import scipy.integrate as integrate
import GEOMODEL as GEOMODEL

# Constants
frequencies = np.logspace(0,5,100000)
omega = 2*np.pi*frequencies

def BScoatingBM():

    # Constants
    kB = 1.38e-23  # Boltzmann constant in J/K
    T = 300  # Temperature in K
    phi_coating = 2.1388e-4  # Loss angle for coating
    elastic_energy = 4.9098e-13
    
    # Make NSR for noise Curve
    BS_coating_BM_noise = np.sqrt((8 * kB * T * elastic_energy*phi_coating) / (np.pi*2*frequencies))
    
    return BS_coating_BM_noise


def BSsubstrateBM():
    
    # Constants
    kB = 1.38*(10**-23)
    T=300
    
    # Make ASD of strain for noise
    BS_substrate_BM_noise = np.sqrt((1.97*(10**-9))*((8*kB*T*(10**-8))/(np.pi*2*frequencies)))
    
    return BS_substrate_BM_noise


def BSthermorefractive():
    
    # Constants
    kB = 1.38*(10**-23)
    kap = 1.38
    T=300
    beta=8.5*(10**-6)
    a = (8*(10**-2))/(np.cos(29*(np.pi/180)))
    eta=1.23
    k=8.56*(10**6)
    r=0.71*(10**-2)
    C=746
    p=2200
    lth = np.sqrt(kap/(C*p*2*np.pi*frequencies))
    
    # Make ASD of strain for noise
    BS_substrate_TR_noise = np.sqrt(((4*kB*kap*(T**2)*(beta**2)*a*(eta + (eta**-1)))/(np.pi*((C*p*(r**2)*(2*np.pi*frequencies))**2)*2*(eta**2)))*(1+((2*(k**2)*(r**2)*eta)/((eta + (eta**-1))*(1+((2*k*lth)**4))))))
    
    return BS_substrate_TR_noise


def TMcoatingBM():
    
    # Constants
    kB = 1.38e-23  # Boltzmann constant in J/K
    T = 300  # Temperature in K
    phi_coating = 5*(10**-5)  # Loss angle for coating
    phi_coatingt = 2.4*(10**-4)
    h = 5.32e-6
    Y = 72*(10**9)
    Yt = 140*(10**9)
    r_o = (25e-3)/np.sqrt(2)
    sigma = 0.17
    sigmat=0.23
    elastic_energy = (h/(4*np.pi*(r_o**2)))*((1/Y)*(((1+sigma)*(1-2*sigma))/(1-sigma)) + (1/Y)*((((1+sigma)**2)*((1-2*sigma)**2))/(1-(sigma**2))))
    elastic_energyt = (h/(4*np.pi*(r_o**2)))*((1/Yt)*(((1+sigmat)*(1-2*sigmat))/(1-sigmat)) + (Yt/(Y**2))*((((1+sigma)**2)*((1-2*sigma)**2))/(1-(sigmat**2))))
    
    # Make NSR for noise Curve
    TM_coating_BM_noise1 = (8 * kB * T * elastic_energy*phi_coating) / (np.pi*2*frequencies)
    TM_coating_BM_noise2 = (8 * kB * T * elastic_energyt*phi_coatingt) / (np.pi*2*frequencies)
    
    TM_coating_BM_noise_sum = np.sqrt(TM_coating_BM_noise1**2 + TM_coating_BM_noise2**2)
    TM_coating_BM_noise = np.sqrt(TM_coating_BM_noise_sum)
    
    return TM_coating_BM_noise


def TMcoatingTE():
    
    # Constants
    kB = 1.38e-23  # Boltzmann constant in J/K
    kap = 0.5
    kapt = 0.6
    T=300
    alpha = 5.1*(10**-7) #0.00000309 
    alphat = 3.6*(10**-6)
    h = 2.66e-7
    Y = 72*(10**9)
    Yt = 140*(10**9)
    C=746
    Ct=306
    p=2200
    pt=6850
    r_o = (25e-3)/np.sqrt(2)
    sigma = 0.17
    sigmat = 0.23
    M = (1/(1-sigma))*((-(1+sigma)/Y) + (((1+sigma)*(2*sigma-1))/(Y)))
    Mt = (1/(1-sigmat))*((-(1+sigmat)/Yt) + (((1+sigma)*(2*sigma-1))/(Y)))
    
    # Make NSR for noise Curve
    TM_coating_TE_noise1 = ((np.sqrt(2)*(alpha**2)*(Y**2)*(M**2)*kB*(T**2)*(h**2))/(np.pi*np.sqrt(kap*p*C)*(r_o**2)*np.sqrt(2*np.pi*frequencies)))*90 #had to scale by extra factor of 3 (20 layers * 1.5-cause folded arm mirror- * 3)
    TM_coating_TE_noise2 = ((np.sqrt(2)*(alphat**2)*(Yt**2)*(Mt**2)*kB*(T**2)*(h**2))/(np.pi*np.sqrt(kapt*pt*Ct)*(r_o**2)*np.sqrt(2*np.pi*frequencies)))*90
    
    TM_coating_TE_noise_sum=np.sqrt(TM_coating_TE_noise1**2 + TM_coating_TE_noise2**2)
    TM_coating_TE_noise = np.sqrt(TM_coating_TE_noise_sum)
    
    return TM_coating_TE_noise


def TMsubstrateBM():

    # Constants
    kB = 1.38e-23  # Boltzmann constant in J/K
    T = 300  # Temperature in K
    phi_coating = 10**-8  # Loss angle for coating
    Y = 72*(10**9)
    r_o = (25e-3)/np.sqrt(2)
    sigma = 0.17
    elastic_energy = (1-(sigma**2))/(2*np.sqrt(2*np.pi)*Y*r_o)
    
    # Make NSR for noise Curve
    TM_substrate_BM_noise = np.sqrt((8 * kB * T * elastic_energy*phi_coating) / (np.pi*2*frequencies))
    
    return TM_substrate_BM_noise


def TMsubstrateTE():
    
    # Constants
    kB = 1.38e-23  # Boltzmann constant in J/K
    kap = 1.38
    T=300
    alpha = 5.1*(10**-7)
    ks = 2464
    C=746
    p=2200
    r_o = (25e-3)/np.sqrt(2)
    sigma = 0.17
    
    # Make NSR for noise Curve
    TM_substrate_TE_noise = np.sqrt(((8*kap*(alpha**2)*((1+sigma)**2)*kB*(T**2))/(np.sqrt(2*np.pi)*(p**2)*(C**2)*(r_o**3)*((np.pi*2*frequencies)**2)))*(1 + ((ks*r_o)/(np.sqrt(2*np.pi)))))
    
    return TM_substrate_TE_noise

def ETMsubstrateTE():

    # Constants
    kB = 1.38e-23  # Boltzmann constant in J/K
    kap = 1.38
    T=300
    alpha = 5.1*(10**-7)
    ks = 2464
    C=746
    p=2200
    r_o = (25e-3)/np.sqrt(2)
    sigma = 0.17
    
    
    # Make NSR for noise Curve
    ETM_substrate_TE_noise = np.sqrt((8*kap*(alpha**2)*((1+sigma)**2)*kB*(T**2))/(np.sqrt(2*np.pi)*(p**2)*(C**2)*(r_o**3)*((np.pi*2*frequencies)**2)))
    
    return ETM_substrate_TE_noise


def pendulumModeTE():
    
    # Constants
    kB = 1.38e-23  # Boltzmann constant in J/K
    T=300
    omega_o = 2*np.pi*GEOMODEL.adjust_model().py.fz
    phi= 10**-8
    m = 5.6
    
    # Make NSR for noise Curve
    pendulum_Mode_TE_noise = np.sqrt((4*(omega_o**2)*kB*T*phi)/(m*omega*((omega_o**2 - omega**2)**2 + (omega_o**4) * (phi**2))))
    
    return pendulum_Mode_TE_noise

def violinModeTE():
    
    # Constants
    kB = 1.38e-23  # Boltzmann constant in J/K
    T=300
    g = 9.81
    rho = 2200 #density of the fiber material 
    m = 5.6 
    L = 0.28 #length of the pendulum 
    P = (m*g)/4 #tension 
    R = 112.5e-6#fiber radius 
    rho_l = rho*np.pi*(R**2) #material linear mass density
    E = 7.2*(10**10) #young's modulus
    I = (np.pi*(R**4))/4 #bending moment of inertia
    d_s = 100e-6 #dissipation depth
    D = 0.00000084084816 #coefficient of thermal diffusion
    A = 0.000000039 # cross section area of fiber
    omega_o = 2*np.pi*GEOMODEL.adjust_model().py.fz
    alpha = 5.1*(10**-7)
    beta = 2*(10**-4)
    u_o = P/(A*E) #static strain
    tau_R = (4*(R**2))/(13.55*D)
    C_v = 1.64 * (10**6) #heat capacity per volume
    phi_bulk = 1.8*(10**-8)
    
    def Dn_1(n):
        return (2/(L*np.sqrt((np.sqrt((P**2)+4*E*I*rho_l*(omega**2))+P)/(2*E*I))))*(1 + (4 + (((n*np.pi)**2)/2))*(1/(L*np.sqrt((np.sqrt((P**2)+4*E*I*rho_l*(omega**2))+P)/(2*E*I)))))
    
    def omega_n(n):
        omega_n = ((n*np.pi)/L)*np.sqrt(P/(rho_l))*(1 + (2/L)*np.sqrt((E*I)/P) + (4 + (((n*np.pi)**2)/2))*(2/(L**2))*((E*I)/P))
        return omega_n
    
    def phi_fiber(n):
        phi_nonlin = (E*((alpha - u_o*beta)**2)*T*omega*tau_R)/(C_v*(1+((omega*tau_R)**2)))
        return (Dn_1(n))*((1+((8*d_s)/R))*phi_bulk+phi_nonlin)
    def phi_n(n):
        return phi_fiber(n)*Dn_1(n)
    def displacement_noise(n):
         return np.sqrt((8*(omega_o**2)*kB*T*phi_fiber(n))/((m*omega)*(((omega_n(n)**2 - (omega**2))**2) + (omega_n(n)**4) * (phi_n(n)**2))))
    
    # Make NSR for noise Curve
    total_noise_squared = np.zeros_like(omega)
    for n in range(1, 33):  # Evaluate for the first 32 modes
        total_noise_squared += displacement_noise(n)**2
    
    total_noise = np.sqrt(total_noise_squared)
    return total_noise