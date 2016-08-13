import numpy as np
from elastopy import dynamics
from elastopy import gmsh
from elastopy import plotter
from elastopy import data
from diffuspy import steadystate
from diffuspy import plotter as plottert
from diffuspy import transient


mesh_name = 'dam-2mat'

model = gmsh.Parse(mesh_name)
# plotter.model(model, nodes_label=True)
# plotter.show()
material = data.Collect()
s = list(model.surf.keys())

# Concrete Dam
material.cndtvt[s[0]] = 0.89     # J/(m s C)
# http://www.pietrangeli.it/assets/uploads/pubblicazioni/file/536a0225156d5.pdf
material.spcfht[s[0]] = 820      # J/kg C
material.dnsty[s[0]] = 2500     # kg/m3
material.E[s[0]] = 21*1e+9      # N/m2
material.nu[s[0]] = 0.2
material.thrmexp[s[0]] = 9.8*1e-6  # 1/C

# Foundation Basalt
# http://www.simetric.co.uk/si_materials.htm
material.dnsty[s[1]] = 3011
material.cndtvt[s[1]] = 5
material.thrmexp[s[1]] = 5.4*1e-6
# http://ceae.colorado.edu/~amadei/CVEN5768/PDF/NOTES5.pdf
material.E[s[1]] = 62.6*1e+9
material.nu[s[1]] = 0.33
# http://www.endmemo.com/chem/specificheatsearch.php?q=Basalt%20Rock
material.spcfht[s[1]] = 840


def b_force(x1, x2, t=1):
    return np.array([0.0,
                     (x2 >= 70.0) * -material.dnsty[s[0]] * 10 +
                     (x2 < 70.0) * -material.dnsty[s[1]] * 10])

gma = 9820.0


def trac_bc(x1, x2, t=1):
    return {('line', 7): [(x2 <= 186.13)*(gma*(186.13-x2)), 0.0],
            ('line', 8): [(17.07 * gma - 116.13 * gma) * (x2 - 70.0) /
                          (169.06 - 70.0) + 116.13 * gma,
                          -((17.07 * gma - 116.13 * gma) * (x2 - 70.0) /
                            (169.06 - 70.0) + 116.13 * gma)],
            ('line', 9): [0.0, - 116.13 * gma]}


def displ_bc(x1, x2):
    return {('line', 0): [0.0, 0.0],
            ('line', 10): [0.0, 'free'],
            ('line', 1): [0.0, 'free']}


def b_heat(x1, x2, t=1):
    if t == 0:
        return 0
    else:
        return (x2 >= 70) * 0.05 # 330e3 * material.dnsty[s[0]] * np.exp(-0.3*t) * 0.3

print(b_heat(1, 70, 80))
# x, y = [], []
# for t in range(0, 30, 1):
#     x.append(t)
#     y.append(b_heat(1, 70, t))
# import matplotlib.pyplot as plt
# plt.plot(x, y, '-b')
# plotter.show()

def flux_bc(x1, x2, t=1):
    return {}


def temp_bc(x1, x2, t=1):
    return {0: 10,
            1: 10/70*x2 + 10,
            2: 20 + 15*np.sin(2*np.pi/(24*60*60)*(t-9*60*60)),
            3: 20 + 15*np.sin(2*np.pi/(24*60*60)*(t-9*60*60)),
            4: 20 + 15*np.sin(2*np.pi/(24*60*60)*(t-9*60*60)),
            5: 20 + 15*np.sin(2*np.pi/(24*60*60)*(t-9*60*60)),
            6: 20 + 15*np.sin(2*np.pi/(24*60*60)*(t-9*60*60)),
            7: 20 + 15*np.sin(2*np.pi/(24*60*60)*(t-9*60*60)),
            8: (10/(169-70)*(x2 - 70)) + 10 + 10*np.sin(2*np.pi /
                                                   (24*60*60)*(t-9*60*60)),
            9: 10,
            10: 10}

# import matplotlib.pyplot as plt
# import numpy as np
# x = np.linspace(0, 24)
# y = np.sin(2 * np.pi/(24)*(x-9))
# plt.plot(x, y, '-b')
# plt.xlim(0,24)
# plt.show()

# T = steadystate.solver(model, material, b_heat, flux_bc, temp_bc, t=12*60*60)
# plottert.contour(model, T)
# plotter.show()

T0 = 15

t_int = 24*60*60
dt = t_int/240

# frames = transient.solver(model, material, b_heat, flux_bc,
#                           temp_bc, T0, t_int, dt, vmin=0., vmax=100.)
# Tani = plottert.anime(frames, dt)
# plotter.show()

U, SIG = dynamics.solver(model, material, b_force,
                         trac_bc, displ_bc, t_int, dt,
                         b_heat=b_heat, flux_bc=flux_bc,
                         temp_bc=temp_bc, T0=T0)

from elastopy import stress

spmax2, spmin2 = [], []
time = []
for time_index in range(int(t_int/dt)+1):
    spmaxl =stress.principal_max(SIG[:, 0, time_index],
                                 SIG[:, 1, time_index], SIG[:, 2, time_index])
    spminl = stress.principal_min(SIG[:, 0, time_index],
                                  SIG[:, 1, time_index], SIG[:, 2, time_index])
    spmax2.append(spmaxl[105])
    spmin2.append(spminl[105])
    time.append(time_index*dt/(60*60))
import matplotlib.pyplot as plt
# plt.plot(time, spmax2, '-b')
plt.plot(time, spmin2, '-g')
plt.show()

# time = [2, 8, 12, 16, 20]

# for ti in time:
#     time_index = int((t_int/dt)/24*ti)

#     # spmax = stress.principal_max(SIG[:, 0, time_index],
#     #                              SIG[:, 1, time_index], SIG[:, 2, time_index])
#     # print(np.amax(spmax), np.amin(spmax))
#     # plotter.stresses(model, SIG[:, :, time_index], ftr=1e6, spmax=True,
#     #                  vmin=-3.5, vmax=2.5, lev=20)
#     # import matplotlib.pyplot as plt
#     # plt.gca().set_title('Time: '+str(ti)+'h')

#     spmin = stress.principal_min(SIG[:, 0, time_index],
#                                  SIG[:, 1, time_index], SIG[:, 2, time_index])
#     print(np.amax(spmin), np.amin(spmin))
#     plotter.stresses(model, SIG[:, :, time_index], ftr=1e6, spmin=True,
#                      vmin=-8.3, vmax=-1.2, lev=20)
#     plt.gca().set_title('Time: '+str(ti)+'h')

#     # plotter.model_deformed(model, U[:, time_index], magf=1000, color='k', ele=True)
#     # plt.gca().set_title('Time: '+str(ti)+'h')

#     # T = steadystate.solver(model, material, b_heat, flux_bc,
#     #                        temp_bc, t=ti*60*60)
#     # plottert.contour(model, T, title='Time: '+str(ti)+'h', vmin=0, vmax=50)
# plotter.show()


# plotter.stress_animation(SIG, model, t_int, dt, spmax=True,
#                          name="stre-cilamce1.gif", vmin=-50, vmax=150,
#                          interval=100, ftr=1e6)

# plotter.displ_animation(U, model, t_int, dt, name='displ-cilamce1.gif', magf=100,
#                         interval=100)
