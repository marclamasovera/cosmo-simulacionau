import numpy as np
import matplotlib.pyplot as plt
## PARAMETERS 

c = 299792.458 # km/s
mpc = 3.085677581e19 # km
spy = 3.15576e7 # segons en un any
doesnt_converge = 1e6


class ship:
    def __init__(self, Omega_m, Omega_r, Omega_Lambda,is_light,initial_velocity,target_distance=1.0):
        self.H_0 = 70
        self.Omega_m = Omega_m
        self.Omega_r = Omega_r
        self.Omega_Lambda = Omega_Lambda
        self.Omega_k = 1 - Omega_m - Omega_r - Omega_Lambda
        self.is_light = is_light
        self.initial_velocity = initial_velocity
        self.target_distance = target_distance

    def Hubble_parameter(self,a):
        term = self.Omega_r / a**4 + self.Omega_m / a**3 + self.Omega_Lambda + self.Omega_k / a**2
        try:
            return 1.0 * np.sqrt(term)
        except ValueError as e:
            print(f"Error in Hubble_parameter: {e}")
            return np.nan


    def peculiar_velocity(self,a : float) -> float:
        if self.is_light:
            return self.initial_velocity
        else:
            v0 = self.initial_velocity
            gamma0 = 1/np.sqrt(1 - v0**2)
            u = gamma0 * v0 / a
            return u / np.sqrt(1 + u**2)

    def comoving_distances(self,b) -> float:
        dx = 0.01
        a_start = 1
        n = (int((b-a_start)/dx) + 1)
        a = np.linspace(a_start,b,n)
        f = self.peculiar_velocity(a) / (a**2 * self.Hubble_parameter(a))     #int_a_start^b v(a)/(a^2 H(a)) da
        comoving_distance = np.trapezoid(f,a)
        return comoving_distance 

    def maximum_distance_reacheable(self):
        return self.comoving_distances(doesnt_converge) 

    def find_scale_factor(self,
                          tolerance = 1e-5,
                          doesnt_converge=doesnt_converge,
                          max_iter=100):

        a_low = 1.0
        a_high = 2.0

        while self.comoving_distances(a_high) - self.target_distance < 0:
            a_high *= 2.0
            if a_high > doesnt_converge:
                print(f"The maximum distance is {self.comoving_distances(a_high)} (c/H_0) ")
                raise ValueError("La nau no acaba mai d'arribar a destinació. Una mica com Rodalies")

        for trial in range(max_iter):
            a_mid = (a_high + a_low)/2
            mid_value = self.comoving_distances(a_mid) - self.target_distance

            if abs(mid_value) < tolerance or (a_high - a_low)/2.0 < tolerance:
                print(f"Converged after {trial} iterations")
                return a_mid

            # Bolzano
            if mid_value < 0:
                a_low = a_mid
            else:
                a_high = a_mid

        return (a_high+a_low)/2.0



    def Universe_age_at_point(self):
        a_reach = self.find_scale_factor()
        dx = 0.01
        n = int((a_reach)/dx + 1) # a_start = 0
        a = np.linspace(1e-13,a_reach,n)# posem 1e-13 per evitar singularitat a a=0
        f = 1/(a*self.Hubble_parameter((a))) #int_0^{a_reach} 1/(a H(a)) da


        age = np.trapezoid(f,a)
        
        return age * (mpc/(self.H_0 * spy)) # tornem les unitats a anys

omegas_m = np.linspace(0.1,2.9,150)
omegas_l = np.linspace(2.9,0.1,150)
factors = []
distances = []

for o_m,o_l in zip(omegas_m,omegas_l):
    ship_1 = ship(Omega_Lambda=o_l,
                Omega_m=o_m,
                Omega_r=0.0,
                is_light=False,
                initial_velocity=1/2 #(c=1))
    )
    print(f" Thinking... Don't stall, I'm still here")
    d_max = ship_1.maximum_distance_reacheable()
    distances.append(d_max)
    print("Finished for now, let's go to the next one!")

    try:
        a = ship_1.find_scale_factor()
        print(f"Scale factor when the ship reaches its destination: {a:.2e}")
    except ValueError as e:
        print(e)
        a = np.inf

    factors.append(a)

fig,ax = plt.subplots()

ax.scatter(omegas_m,distances)
ax.set_xlabel(r'$\Omega_m$')
ax.set_ylabel(r'Maximum distance')
ax.set_title('Maximum distance vs. Matter density')
plt.show()
