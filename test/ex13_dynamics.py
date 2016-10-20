import numpy as np
from elastopy import dynamics
from elastopy import gmsh
from elastopy import plotter
from elastopy import data

mesh_name = 'patch'

model = gmsh.Parse(mesh_name)

material = data.Collect()
s = list(model.surf.keys())

material.cndtvt[s[0]] = 1.
material.spcfht[s[0]] = 1.
material.dnsty[s[0]] = 0.1
material.E[s[0]] = 1000
material.nu[s[0]] = 0.2
material.dnsty[s[0]] = 1.
material.thrmexp[s[0]] = 1e-6


def b_force(x1, x2, t=1):
    return np.array([0.0,
                     0.0])


def trac_bc(x1, x2, t=1):
    return {('line', 3): [0., 0.0],
            ('line', 1): [0., 0.0]}


def displ_bc(x1, x2):
    return {('node', 0): [0.0, 0.0],
            ('node', 1): ['free', 0.0]}


def b_heat(x1, x2, t=1):
    return 0


def flux_bc(x1, x2, t=1):
    return {1: 0.0,
            2: 0.0,
            0: 0.0}


def temp_bc(x1, x2, t=1):
    return {3: 1.}

t_int = 60*60*10
dt = t_int/5

U, SIG = dynamics.solver(model, material, b_force,
                         trac_bc, displ_bc, t_int, dt)
                         # b_heat=b_heat, flux_bc=flux_bc, temp_bc=temp_bc)

plotter.stress_animation(SIG, model, t_int, dt, spmax=True,
                         name="stre-nobc.gif", vmin=-.5, vmax=.5,
                         interval=100)

plotter.displ_animation(U, model, t_int, dt, name='displ-nobc.gif', magf=10000,
                        interval=100)

"""
PLOT SPMAX, SPMIN, T AND U FIELD AT SPECIFIC TIMES
"""
time = np.linspace(0, 9, 10)  # in [h]


T_point = []
time_point = []
for ti in time:
    time_index = int((t_int/dt)/t_int*(60*60*ti))

    """
    PRINCIPAL STRESS MAXIMUM
    """
    # spmax = stress.principal_max(SIG[:, 0, time_index],
    #                              SIG[:, 1, time_index],
    #                              SIG[:, 2, time_index])
    # print(np.amax(spmax)/1e6, np.amin(spmax)/1e6)
    # plotter.stresses(model, SIG[:, :, time_index], ftr=1e6, spmax=True,
    #                  vmin=-5.9, vmax=5.5, lev=20, title='Time: '+str(ti)+'h')

    """
    PRINCIPAL STRESS MINIMUM
    """
    # spmin = stress.principal_min(SIG[:, 0, time_index],
    #                              SIG[:, 1, time_index],
    #                              SIG[:, 2, time_index])
    # print(np.amax(spmin)/1e6, np.amin(spmin)/1e6)
    # plotter.stresses(model, SIG[:, :, time_index], ftr=1e6, spmin=True,
    #                  vmin=-12.2, vmax=0, lev=20, title='Time: '+str(ti)+'h')

    # DISPLACEMENT
    plotter.model_deformed(model, U[:, time_index], magf=1000, color='r')
    plt.gca().set_title('Time: '+str(ti)+'h')

    # # TEMPERATURE
    # T = steadystate.solver(model, material, b_heat, flux_bc,
    #                        temp_bc, t=ti*60*60)
  #   T_point.append(T[105])
#     time_point.append(ti)
#     plottert.contour(model, T, title='Time: '+str(ti)+'h', vmin=0, vmax=90)

    plt.savefig('displ'+str(ti)+'h.pdf', bbox_inches='tight')

# fig, ax = plt.subplots()
# ax.plot(time_point, T_point, '-b')
# ax.set_xlabel('Time (h)')
# ax.set_ylabel('Temperature (C)')
# plt.savefig('T_evo.pdf', bbox_inches='tight')

plt.show()
