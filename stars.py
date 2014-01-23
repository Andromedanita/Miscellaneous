import numpy as np
import astropy.units as u
from astropy.constants import G, h, c, m_e, m_p
from astropy.table import Table


Knr = 1./20. * (3/np.pi)**(2./3.) * h**2 / m_e
Ker = 1./8. * (3/np.pi)**(1./3.) * h * c


def mche_r(r, y):
    """Mass conservation and hydrostatic eq. and with radius as indep. var.

    d M_r
    ----- = 4 pi r^2 rho
     d r

    d P     G M_r rho
    --- = - ---------
    d r        r^2

    Parameters
    ----------
    r : float
      independent variable: current radius
    y : tuple of two dependent variables
      current enclosed mass, density

    Returns
    -------
    dmdr, dddr : tuple of two derivatives
      mass conservation and hydrostatic equilibrium derivatives
    """
    m, d = y[0], abs(y[1])  # need absolute value since d<0 at last step
    # mass conservation
    dMdr = 4.*np.pi*r**2*d
    # hydrostatic equilibrium
    dPdr = -G*m*d/r**2
    P, dPdd = eos(d)
    dddr = dPdr/dPdd
    return dMdr, dddr


def polytrope(d):
    """
    Polytropic equation of state: P = K * d ** gamma

    Density is converted to a number density by dividing by mu_e * m_p

    """
    P = K * (d/mu_e/m_p) ** gamma
    dPdd = gamma*P/d
    return P, dPdd


def rk4(x, dx, y, f):
    """Runge-Kutta 4th order integration"""
    add = lambda y1, y2: [_y1 + _y2 for _y1, _y2 in zip(y1, y2)]
    mul = lambda dx, y: [dx * _y for _y in y]
    k1 = f(x, y)
    k2 = f(x+0.5*dx, add(y, mul(0.5*dx, k1)))
    k3 = f(x+0.5*dx, add(y, mul(0.5*dx, k2)))
    k4 = f(x+dx, add(y, mul(dx, k3)))
    ynew = add(y, mul(dx/6., add(k1, add(mul(2.,k2),
                                         add(mul(2.,k3), k4)))))
    ynew = [_ynew.to(_yold.unit) for _ynew, _yold in zip(ynew, y)]
    return x + dx, ynew


def integrate(dc, dr=1.e-4*u.Rsun):
    """
    Integrate a given equation of state for a given central density
    """
    # starting r, M_r, d
    r = [0.*dr.unit]
    mr = [0.*u.Msun]
    d = [dc]
    # use polytropic expansion for first step: theta = 1-1/6*xi**2
    # do decompose to avoid rounding in exponents in the next line
    Pc = eos(dc)[0].decompose()
    alpha = np.sqrt(gamma * Pc / (4.*np.pi*G*(gamma-1.)*dc**2)).to(dr.unit)
    r.append(dr)
    xi = dr / alpha
    theta = 1. - 1./6. * xi**2
    n = 1./(gamma-1.)
    d.append(dc*theta**n)
    mr.append((dc*alpha**3*4./3.*np.pi*xi**3*(1.-n*xi**2/10.)).to(u.Msun))
    while d[-1] > 0:
        y = [mr[-1], d[-1]]
        rnew, ynew = rk4(r[-1], dr, y, mche_r)
        r.append(rnew)
        mr.append(ynew[0])
        d.append(ynew[1])

    # last step will have overshot to d<0; interpolate to d=0
    fraction = -d[-1]/(d[-2]-d[-1])
    r[-1] = r[-2] + fraction * (r[-1] - r[-2])
    mr[-1] = mr[-2] + fraction * (mr[-1] - mr[-2])
    d[-1] = 0. * d[-1].unit
    # u.Quantity makes this a single array rather than a list
    return u.Quantity(r), u.Quantity(mr), u.Quantity(d)


if __name__ == '__main__':
    eos = polytrope
    gamma = 5./3.
    mu_e = 2.
    K = Knr
    r, mr, d = integrate(1.e6*u.g/u.cm**3)
    star = Table([r, mr, d], names=('r', 'm', 'd'))
    import matplotlib.pylab as plt
    plt.plot(star['r'], star['d'])
    plt.show()
