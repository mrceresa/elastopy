def coefficients(dt, beta, gamma):
    """Generates the coefficients used in the procedure

    """
    c0 = 1/(beta * dt**2)
    c2 = 1/(beta * dt)
    c3 = 1/(2 * beta) - 1
    c6 = dt*(1 - gamma)
    c7 = gamma*dt
    return c0, c2, c3, c6, c7
