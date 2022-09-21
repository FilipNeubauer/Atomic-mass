import cmath
from turtle import color
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import cm
from scipy.optimize import curve_fit, fsolve

import numpy as np

mp = 938.2720882
mn = 939.5654205
me = 0.51099895

pole_magicky = np.zeros((161, 111))

z_pole = []
n_pole = []
b_pole = []

m_pole = []


pole = np.zeros((161, 111))

pole_hmotnosti = np.zeros((161, 111))

f = open("mass_1.mas20.txt", "r")
lines = f.readlines()

all_nuclids = []

binding_energy_list = []

# for i in lines:
#     print(i)
for number, line in enumerate(lines):
    if number < 36:
        continue
    nuclid = []
    n = int(line[6:10])
    nuclid.append(n)

    

    z = int(line[11:15])

    

    nuclid.append(z)
    a = int(line[15:20])
    nuclid.append(a)

    try:
        binding_energy = float(line[57:68]) / 1000
    except:
        continue

    nuclid.append(binding_energy)

    all_nuclids.append(nuclid)

    binding_energy_list.append([binding_energy, n, z])

    if a >= 3:

        z_pole.append(z)

        b_pole.append(binding_energy)

        n_pole.append(n)

    


    pole[n, z] = binding_energy

    if n != 0 or z != 0:
        celkova_hmotnost = mp * z + mn * n - binding_energy * a
        pole_hmotnosti[n, z] = celkova_hmotnost

    if celkova_hmotnost != 0 and pole_hmotnosti[n-2, z]:
        pole_magicky[n, z] = - celkova_hmotnost + mn*2 + pole_hmotnosti[n-2, z]





#print(pole)
sorted_list = sorted(binding_energy_list, key=lambda x: x[0], reverse=True)
#print(sorted_list[0])

# plt.imshow(pole.transpose(), origin="lower", interpolation="none")
# plt.colorbar(label="MeV")
# plt.xlabel("N")
# plt.ylabel("Z")
# plt.show()


#print(pole_hmotnosti)

hmotnost_helia = pole_hmotnosti[2, 2]
#print(hmotnost_helia)

pole_rozpadu = np.zeros((161, 111))

for n in range(1, 161):
    for z in range(1, 111):
        hmotnost_rozpadajiciho = pole_hmotnosti[n, z]
        if hmotnost_rozpadajiciho == 0:
            continue
        
        
        hmotnost_konecnyho = pole_hmotnosti[n-2, z-2]

        if hmotnost_konecnyho == 0:
            continue

        if hmotnost_rozpadajiciho > hmotnost_konecnyho + hmotnost_helia:
            pole_rozpadu[n, z] = 1
        else:
            pole_rozpadu[n, z] = -1



# -----první obrázek, rovnice-------
cmap_alfa = ["yellow", "purple", "white"]
plt.imshow(pole_rozpadu.transpose(), origin="lower", interpolation="none", cmap=cmap_alfa)
#plt.colorbar(label="MeV")
#plt.title("Alfa rozpad")
plt.xlabel("N")
plt.ylabel("Z")
plt.show()



# hmotnost protonu 938.2720882 MeV
# neutron 939.5654205 MeV

# -----beta rozpad, porovnání s nudat, rovnice-------
cmap_beta = ListedColormap(["white", "black", "#e58dc4", "#71cade", "orange"])


pole_rozpadu_beta = np.zeros((161, 111))

for n in range(161):
    for z in range(111):
        vysledek = 0
 
        if z+1 >= 110 or pole_hmotnosti[n,z] == 0 or pole_hmotnosti[n-1, z+1] == 0 or pole_hmotnosti[n+1, z-1] == 0:
            vysledek = -1
        else:
            hmotnostRozpadajiciho = pole_hmotnosti[n, z]
            hmotnostBmVysledku = pole_hmotnosti[n-1, z+1] + me
            hmotnostBpVysledku = pole_hmotnosti[n+1, z-1]
            betap = 0
            betam = 0
            if hmotnostRozpadajiciho > hmotnostBmVysledku:
                betam = 1
            if hmotnostRozpadajiciho + me > hmotnostBpVysledku:
                betap = 1
            if betam and betap:
                vysledek = 3
            else:
                if betam:
                    vysledek = 1
                if betap:
                    vysledek = 2
 
        pole_rozpadu_beta[n,z] = vysledek



            


plt.imshow(pole_rozpadu_beta.transpose(), origin="lower", interpolation="none", cmap=cmap_beta)
#plt.colorbar()
plt.xlabel("N")
plt.ylabel("Z")
plt.show()

# beta = 1
# beta+ = -1







# ------formule-----
def formule(x, a_v, a_s, a_c, a_A, a_p):
    z, n = x
    a = n + z
    b = a_v - a_s * a**(-1/3) - a_c * z**2 * a**(-4/3) - a_A * ((n-z)/a)**2
    delta = a_p * a**(-3/2)


    if type(z) is int:
        if z % 2 == 1 and n % 2 == 1:
            b -= delta
        elif z % 2 == 0 and n % 2 == 0:
            b += delta
    else:
        for i in range(len(z)):
            if z[i] % 2 == 1 and n[i] % 2 == 1:
                b[i] -= delta[i]
            elif z[i] % 2 == 0 and n[i] % 2 == 0:
                b[i] += delta[i]

    # if z % 2 == 1 and n % 2 == 1:
    #     b -= delta
    # elif z % 2 == 0 and z % 2 == 0:
    #     b += delta

    return b





vysledek = curve_fit(formule, (z_pole, n_pole), b_pole)

print(vysledek[0])


pole_rozdilu = np.zeros((161, 111))

# print(len(z_pole))
# print(len(n_pole))
# print(len(b_pole))

for i in range(len(z_pole)):

    n = n_pole[i]
    z = z_pole[i]
    b = b_pole[i]

    rozdil = b - formule((z, n), vysledek[0][0], vysledek[0][1], vysledek[0][2], vysledek[0][3], vysledek[0][4])
    pole_rozdilu[n, z] = rozdil

# print(pole_rozdilu)

# ------ dalsi obrazek-----
plt.imshow(pole_rozdilu.transpose(), origin="lower", interpolation="none", vmin=-0.2, vmax=+0.2, cmap=cm.seismic)
#plt.imshow(pole_rozdilu.transpose(), origin="lower")
#plt.colorbar()
plt.xlabel("N")
plt.ylabel("Z")
plt.show()

# magicky_cisla = [8, 20, 28, 50, 82, 126]

# for i in magicky_cisla:
#     plt.plot([i, i], [0, 110], color="darkgrey")
#     if i != 126:
#         plt.plot([0, 160], [i, i], color="grey")


# plt.show()


# print(pole_magicky)

plt.imshow(pole_magicky.transpose(), origin="lower", interpolation="none", vmin=0, vmax=30)
plt.colorbar()
plt.xlabel("N")
plt.ylabel("Z")
plt.show()
