import numpy as np
from scipy.integrate import quad
## PARAMETERS 

c = 299792.458 # km/s
mpc = 3.085677581e19 # km
spy = 3.15576e7 # segons en un any
doesnt_converge = 1e6
A_MAX = 1e12


class ship:
    def __init__(self, Omega_m, Omega_r, Omega_Lambda, is_light, initial_velocity, target_distance=1.0):
        self.H_0 = 70
        self.Omega_m = Omega_m
        self.Omega_r = Omega_r
        self.Omega_Lambda = Omega_Lambda
        self.Omega_k = 1 - Omega_m - Omega_r - Omega_Lambda
        self.is_light = is_light
        self.initial_velocity = initial_velocity
        self.target_distance = target_distance

    def Hubble_parameter(self, a):         
        term = self.Omega_r / a**4 + self.Omega_m / a**3 + self.Omega_Lambda + self.Omega_k / a**2
        return np.sqrt(term)

    def peculiar_velocity(self, a):
        if self.is_light:
            return 1.0
        v0 = self.initial_velocity
        u = (v0 / np.sqrt(1 - v0**2)) / a    # gamma0 v0 / a
        return u / np.sqrt(1 + u**2)

    def comoving_distances(self, b, a_start=1.0):
         #int_a_start^b v(a)/(a^2 H(a)) da
        f = lambda lna: self.peculiar_velocity(np.exp(lna)) / (np.exp(lna) * self.Hubble_parameter(np.exp(lna)))
        return quad(f, np.log(a_start), np.log(b), limit=200)[0]

    def bisect_distance(self, a_start, target, tolerance=1e-8, max_iter=200):
        #bolzano
        a_low, a_high = a_start, A_MAX
        dmax = self.comoving_distances(a_high, a_start)
        if dmax < target:
            print(f"The maximum distance is {dmax} (c/H_0)")
            raise ValueError("La nau no acaba mai d'arribar a destinació. Una mica com Rodalies")
        for _ in range(max_iter):
            a_mid = np.sqrt(a_low * a_high)
            if self.comoving_distances(a_mid, a_start) < target:
                a_low = a_mid
            else:
                a_high = a_mid
            if (a_high - a_low) / a_mid < tolerance:
                break
        return np.sqrt(a_low * a_high)

    def find_scale_factor(self):
        return self.bisect_distance(1.0, self.target_distance)

    def age_at(self, a_end):

        f = lambda lna: 1.0 / self.Hubble_parameter(np.exp(lna))
        age = quad(f, -40, np.log(a_end), limit=200)[0]
        return age * mpc / (self.H_0 * spy)

    def Universe_age_at_point(self):
        #int_0^{a_reach} 1/(a H(a)) da
        return self.age_at(self.find_scale_factor())

    def coming_back(self):
        a_arr = self.find_scale_factor()
        light = ship(self.Omega_m, self.Omega_r, self.Omega_Lambda, True, 1.0, self.target_distance)
        a_rec = light.bisect_distance(a_arr, self.target_distance)
        return a_arr, a_rec, self.age_at(a_rec), a_rec / a_arr - 1

        
################################################################################

# sections i, j
# for omega_m in [0.3,1.0]:
#     ship_humans = ship(Omega_m=omega_m, Omega_r=0.0, Omega_Lambda=0.0,is_light=False,initial_velocity=0.5,target_distance=1.0)
#     ship_light = ship(Omega_m=omega_m, Omega_r=0.0, Omega_Lambda=0.0,is_light=True,initial_velocity=1.0,target_distance=1.0)

#     ship_light_age = ship_light.Universe_age_at_point()
#     print(f"Age of ship with light: {ship_light_age:2e}, omega_m: {omega_m}")
#     try:
#         ship_humans_age = ship_humans.Universe_age_at_point()
#     except ValueError as e:
#         print(f"Error for Omega_m={omega_m}: {e} ")
#         continue
#     print(f"Age of ship with humans: {ship_humans_age:2e},omega_m: {omega_m}")

#section l
# for omega_m in [0.3, 1.0]:
#     humans = ship(omega_m, 0.0, 0.0, False, 0.5, 1.0)
#     try:
#         a_arr, a_rec, age_rec, z = humans.coming_back()
#         print(f"Omega_m={omega_m}: a_arr={a_arr:.2f}, a_rec={a_rec:.2f}, "
#               f"age at reception={age_rec:.3e} yr, z={z:.3f}")
#     except ValueError:
#         print(f"Omega_m={omega_m}: the ship never arrives, so there is no return signal")

#section n

Omega_m = 0.3
Omega_r = 0.0
Omega_Lambda = 0.7

trial_ship = ship(Omega_m, Omega_r, Omega_Lambda, False, 0.5, 1.0)

try:
    print(f"Age of the universe at the point of arrival: {trial_ship.Universe_age_at_point():.3e} yr")
except ValueError as e:
    print(f"Error: {e}")