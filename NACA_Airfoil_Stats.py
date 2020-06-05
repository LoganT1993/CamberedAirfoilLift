import math
import numpy as np
from numpy import pi
import matplotlib as mat
import scipy
import sympy as sym


def LiftNACA(c, m, p, alpha, rho, V):

    # This function returns the lift experienced by airfoils with and without camber

    # c = chord length
    # M = first digit, so 5 is 0.05. Max camber relative to the chord length
    # P = 2nd digit, % of distance to point of maximum camber relative to the chord length (ex: 4 is 40%, or 0.4).
    # alpha = angle of attack, degrees

    c = float(c)    # The program messes up if it gets an int instead of a float for the value of c
    x = sym.Symbol('x')
    z = sym.Symbol('z')
    th = sym.Symbol('theta')
    P = 0.1*p   # Must be between 0 and 1
    M = 0.01*m   # Must be between 0 and .095

    z1 = 2*M
    z2 = P**2
    z3 = (1-P)**2

    try:   # If cambered, gotta use the hard formulas with those integrals

        dz = (z1/z2)*(P - x) + (z1/z3)*(P - x)   # Split into two parts
        xth = (1/2)*(1 - sym.cos(th))   # Coordinate transformation from x to theta

        dth1 = (z1/z2)*(P-((1/2)*(1-sym.cos(th))))   # Now in terms of theta, first half
        dth2 = (z1/z3)*(P-((1/2)*(1-sym.cos(th))))   # Now in terms of theta, second half

        # So we have the expression. Now we need to integrate it with respect to theta from 0 to the position of max camber,
        # then from the position of max camber to pi.

        # Now we need to convert the limits of integration. Max camber occurs at 2P = 1-cos(th) -> th = cos^-1(1-2P)

        lim = sym.acos(1-2*P)

        # The angle of attack of zero lift is calculated using the following formula:
        # [alpha_0L = (-1/pi) integral[(dz/dx)(cos(theta)-1) dtheta]

        cosi = (sym.cos(th) - 1)   # Shorthand to express this
        th_1a = sym.Mul(dth1, cosi)   # Concatenate dz/dx and the (cos(theta) - 1) terms, first half
        th_1b = sym.Mul(dth2, cosi)    # Concatenate dz/dx and the (cos(theta) - 1) terms, second half
        alpha_0a = -(1 / pi) * sym.integrate(th_1a, (th, 0, lim))   # Integrate from 0 to the first limit
        alpha_0b = -(1 / pi) * sym.integrate(th_1b, (th, lim, pi))    # Integrate from the first limit to pi
        alpha_0Lrad = alpha_0a+alpha_0b    # Angle of attack for zero lift, radians
        alpha_0Ldeg = (180/pi)*alpha_0Lrad   # Angle of attack for zero lift, degrees

        th_2a = sym.Mul(dth1, sym.cos(th))   # This is the equation that goes into A1, first half
        th_2b = sym.Mul(dth2, sym.cos(th))   # This is the equation that goes into A1, second half
        A1 = (2 / pi) * sym.integrate(th_2a, (th, 0, lim)) + (2 / pi) * sym.integrate(th_2b, (th, lim, pi))   # rads
        th_3a = sym.Mul(dth1, sym.cos(2*th))   # This is the equation that goes into A2, first half
        th_3b = sym.Mul(dth2, sym.cos(2*th))   # This is the equation that goes into A2, first half
        A2 = (2 / pi) * sym.integrate(th_3a, (th, 0, lim)) + (2 / pi) * sym.integrate(th_3b, (th, lim, pi))   # rads

        cl = 2*pi*((pi/180)*alpha - alpha_0Lrad)   # Coefficient of lift: [angle of attack - angle of attack for 0 lift], rads

        LiftForce = 0.5*cl*c*rho*(V**2)   # Lift, in Newtons
        cm4 = (pi / 4) * (A2 - A1)   # Coefficient of moment about the quarter-chord
        cmle = -(cl/4)   # Coefficient of moment about the leading edge
        x_cp = (c/4)*(1 + (pi/cl)*(A1-A2))   # Center of pressure, in meters

        print("The lift per unit span is", f"{LiftForce:.1f}", "Newtons.")
        print("The coefficient of moment about the quarter chord is", f"{cm4:.3f}")
        print("The coefficient of lift is", f"{cl:.3f}")
        print("The coefficient of moment about the leading edge is", f"{cmle:.3f}")
        print("The center of pressure is", f"{x_cp:.3f}", "meters.")
        print('The angle of attack at zero degrees is', f"{alpha_0Lrad:.3f}", "rad")
        print('The angle of attack at zero degrees is', f"{alpha_0Ldeg:.3f}", "degrees")


        return f"{LiftForce:.1f}", f"{cm4:.3f}", f"{cl:.3f}", f"{cmle:.3f}", f"{x_cp:.3f}", f"{alpha_0Lrad:.3f}", f"{alpha_0Ldeg:.2f}"

    except:   # If symmetric, use appropriate equations per thin airfoil theory

        cl = 2*pi*((pi/180)*alpha)
        cm4 = -(cl/4)
        LiftForce = 0.5*cl*c*rho*(V**2)   # Lift, in Newtons
        print("The lift per unit span is", f"{LiftForce:.1f}", "Newtons.")
        print("The coefficient of lift is", f"{cl:.3f}")
        print("The coefficient of moment about the leading edge is", f"{cm4:.3f}")

        return f"{LiftForce:.1f}", f"{cm4:.3f}", f"{cl:.3f}"

#--------------------

cLength = 1.00   # Chord length, m
Vinf = 10   # Speed of the aircraft, m/s
rho = 1.23   # Density of the air, kg/m^3
alpha = 3   # Angle of the airfoil in the air, in degrees

Lift = LiftNACA(cLength, 4, 4, alpha, rho, Vinf)