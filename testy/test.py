import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

z=0.001
r2=2
a=0

quad = integrate.quad

lowBound = abs(z)
highBound = np.sqrt(z**2+r2**2)

def integrandFinite(u):
    exp1 = np.e**(-u*u)
    exp2 = np.e**(2*r2*np.sqrt(u*u-z*z))
    frac1 = (u*u-3*z*z)/(u**4)
    frac2 = (-np.sqrt(u*u-z*z)+r2)**a/np.sqrt(u*u-z*z)
    return exp1*exp2*frac1*frac2

def integrandInfinite(u):
    exp1 = np.e**(-u*u)
    exp2 = np.e**(-2*r2*np.sqrt(u*u-z*z))
    frac1 = (u*u-3*z*z)/(u**4)
    frac2 = (np.sqrt(u*u-z*z)+r2)**a/np.sqrt(u*u-z*z)
    return exp1*exp2*frac1*frac2

def integrandZFinite(z):
    def localIntegrandFinite(u):
        exp1 = np.e**(-u*u)
        exp2 = np.e**(2*r2*np.sqrt(u*u-z*z))
        frac1 = (u*u-3*z*z)/(u**4)
        frac2 = (-np.sqrt(u*u-z*z)+r2)**a/np.sqrt(u*u-z*z)
        return exp1*exp2*frac1*frac2
    return quad(localIntegrandFinite,lowBound,highBound)[0]

def integrandZInfinite(z):
    def localIntegrandInfinite(u):
        exp1 = np.e**(-u*u)
        exp2 = np.e**(-2*r2*np.sqrt(u*u-z*z))
        frac1 = (u*u-3*z*z)/(u**4)
        frac2 = (np.sqrt(u*u-z*z)+r2)**a/np.sqrt(u*u-z*z)
        return exp1*exp2*frac1*frac2
    return quad(localIntegrandInfinite,lowBound,np.inf)[0]

resultFinite = quad(integrandFinite,lowBound,highBound)
resultInfinite = quad(integrandInfinite,lowBound,np.inf)

resultZFinite = quad(integrandZFinite,0.4,10)
resultZInfinite = quad(integrandZInfinite,0.4,10)

#print(resultFinite)
#print(resultInfinite)
#print("--------------\nZIntegration:\n--------------")
#print(resultZFinite)
#print(resultZInfinite)

middle_array = []
val_array = []
def fullIntegral(L):
    def outer(z2):
        print("done outer")
        def middle(r2):
            print("d mid")
            def inner(u, z1):
                z = abs(z1-z2)
                u+=z
                exp1 = np.e**(-u*u)
                exp2 = np.e**(-2*r2*np.sqrt(u*u-z*z))
                frac1 = (u*u-3*z*z)/(u**4)
                frac2 = 1/np.sqrt(u*u-z*z)
                return exp1*exp2*frac1*frac2*np.exp(-(z1-z2)*(z1-z2))

            lower_u_bound = lambda z1: abs(z1 - z2)
            return integrate.dblquad(inner, 0, L, lower_u_bound, np.infty)[0] * np.exp(-r2 * r2)
        midres = quad(middle,0,np.inf)[0]
        middle_array.append(midres)
        val_array.append(z2)
        return midres
    return quad(outer,0,L)[0]



middle_array2=[]
val_array2 = []
def fullIntegrallow(L):
    def outer(z2):
        print("done")
        def middle(r2):
            print("d mid")
            def inner(u, z1):
                z = abs(z1 - z2)
                exp1 = np.e**(-u*u)
                exp2 = np.e**(+2 * r2 * np.sqrt(u*u - z*z))
                frac1 = (u*u - 3 * z*z) / (u**4)
                frac2 = 1 / np.sqrt(u*u - z*z)
                return exp1 * exp2 * frac1 * frac2 * np.exp(-(z1 - z2) * (z1 - z2))

            lower_u_bound = lambda z1: abs(z1 - z2)
            upper_u_bound = lambda z1: np.sqrt(r2**2 + (z1 - z2)**2)
            return integrate.dblquad(inner, 0, L, lower_u_bound, upper_u_bound)[0] * np.exp(-r2 * r2)
        
        midres = quad(middle,0,np.inf)[0]
        middle_array.append(midres)
        val_array2.append(z2)
        return midres
    
    return integrate.quad(outer, 0, L)[0]

result = fullIntegral(1)
print(result)

result2 = fullIntegrallow(1)
print(result2)

with open("integration_results.txt", "w") as file:
    file.write("Integration Result 1: " + str(result) + "\n")
    file.write("Integration Result 2: " + str(result1) + "\n")