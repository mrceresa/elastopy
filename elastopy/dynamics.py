from elastopy import boundary
from elastopy import stiffness
from elastopy import load
from elastopy import traction
from elastopy import stress
from elastopy import mass
from elastopy import plotter
from elastopy import newmark
import numpy as np
import matplotlib.pyplot as plt


def solver(model, material, b_force, trac_bc, displ_bc, t_int, dt,
           EPS0=0, t=1, gamma=0.5, beta=0.25, U0=0, V0=0, plot_u=False,
           plot_s=False, vmin=None, vmax=None, magf=1):
    """Solves the elasto-dynamics problem

    """
    # Number of steps
    N = int(t_int/dt)

    # Generate initial values
    U_t, V_t = initial_values(U0, V0, model)
    K, M, P = generate(model, material, b_force, trac_bc, EPS0, t=0)
    A_t = np.linalg.solve((M + dt**2 * beta * K), (P - K @ U_t))
    SIG = stress.recovery(model, material, U_t, EPS0)

    # PLOT U_updt and SIG here
    if plot_s is True:
        fig_s, ax_s = plotter.initiate()
        # plot stress at t=0
        im_s, s_range = plotter.stresses_dyn(model, SIG, ax_s, spmax=True)
        te_s = ax_s.text(.5, 0, "Time: 0", ha='center', va='baseline',
                         transform=ax_s.transAxes)
        frm_s = []
        frm_s.append(im_s + [te_s])
        vrange = []
        vrange.append(s_range)  # s_range = [max, min]

    if plot_u is True:
        fig_u, ax_u = plotter.initiate()
        frm_u = []
        # plot deformed shape at t=0
        im_u = plotter.model_deformed_dyn(model, U_t, ax_u, ele=True,
                                          magf=magf)
        te_u = ax_u.text(.5, 1, "Time: 0", ha='center', va='top',
                         transform=ax_u.transAxes)
        frm_u.append([im_u, te_u])

    for n in range(1, N):
        t = n*dt
        K, M, P = generate(model, material, b_force, trac_bc, EPS0, t)

        c0, c2, c3, c6, c7 = newmark.coefficients(dt, beta, gamma)

        # Compute matrices of the system K_ U_updt = P_
        K_ = K + c0 * M
        P_ = P + M @ (c0 * U_t + c2 * V_t + c3 * A_t)
        K_m, P_m = boundary.dirichlet(K_, P_, model, displ_bc)

        # Calculate at time t + dt
        U_updt = np.linalg.solve(K_m, P_m)
        A_updt = c0*(U_updt - U_t) - c2*V_t - c3*A_t
        V_updt = V_t + c6*A_t + c7*A_updt
        SIG = stress.recovery(model, material, U_updt, EPS0)

        # Update the variables
        U_t = U_updt
        V_t = V_updt
        A_t = A_updt

        if plot_s is True:
            im_s, s_range = plotter.stresses_dyn(model, SIG, ax_s, spmax=True)
            te_s = ax_s.text(.5, 0, "Time: "+str(round(t, 3)), ha='center',
                             va='baseline',
                             transform=ax_s.transAxes)
            frm_s.append(im_s + [te_s])
            frm, fig = frm_s, fig_s
            vrange.append(s_range)  # s_range = [max, min]

        if plot_u is True:
            im_u = plotter.model_deformed_dyn(model, U_t, ax_u, ele=True,
                                              magf=magf)
            te_u = ax_u.text(0.5, 1, "Time: "+str(round(t, 3)), ha='center',
                             va='top',
                             transform=ax_u.transAxes)
            frm_u.append([im_u, te_u])
            frm, fig = frm_u, fig_u

    if plot_s is True:
        vrange = np.array(vrange)
        print('Min and Max: ', np.amin(vrange), np.amax(vrange))

        # Change the colorbar range
        sm = plt.cm.ScalarMappable(cmap='plasma',
                                   norm=plt.Normalize(vmin=vmin, vmax=vmax))
        # fake up the array of the scalar mappable. Urgh...
        sm._A = []
        cbar = plt.colorbar(sm)
        cbar.set_label(r'Stress')

    return frm, fig


def initial_values(U0, V0, model):
    """Return the initial values if given

    """
    if np.size(U0) == 1:
        U = np.zeros(model.ndof)
        V = np.zeros(model.ndof)
    else:
        U = U0
        V = V0
    return U, V


def generate(model, material, b_force, trac_bc, EPS0, t):
    """Generate the matrices and vectors

    """
    K = stiffness.K_matrix(model, material, t)

    M = mass.M_matrix(model, material, t)

    Pb = load.Pb_vector(model, b_force, t)

    Pt = traction.Pt_vector(model, trac_bc, t)

    Pe = load.Pe_vector(model, material, EPS0, t)

    P = Pb + Pt + Pe

    return K, M, P
