"""trafic model, assignment of numerical mooc module 3.
Sympy is used to solve a set of equations
"""
import sympy
import numpy as np
import matplotlib.pyplot as plt
#~ sympy.init_printing()

# just following the ipython notebook...
u_max, u_star, rho_max, rho_star, A, B = sympy.symbols('u_max u_star rho_max rho_star A B')

eq1 = sympy.Eq( 0, u_max*rho_max*(1 - A*rho_max - B*rho_max**2) ) # eq (5) in ipython notebook
eq2 = sympy.Eq( 0, u_max*(1 - 2*A*rho_star - 3*B*rho_star**2) )   # eq (6)
eq3 = sympy.Eq( u_star, u_max*(1 - A*rho_star - B*rho_star**2) )  # eq (7)

eq4 = sympy.Eq(eq2.lhs - 3*eq3.lhs, eq2.rhs - 3*eq3.rhs)
#~ print(eq4.simplify())

rho_sol = sympy.solve(eq4, rho_star)[0]
#~ print("rho_star: ", rho_sol)
B_sol = sympy.solve(eq1, B)[0]
#~ print(B_sol)
quadA = eq2.subs([(rho_star, rho_sol), (B,B_sol)])
#~ print(quadA.simplify())
A_sol = sympy.solve(quadA, A)
#~ print(A_sol[0])
#~ print(A_sol[1])

model_params = {u_star: 1.5, u_max:2.0, rho_max:15.0} # model parameters changed for assignment 
aval = A_sol[0].evalf(subs=model_params) 
print("A:", aval)
model_params[A] = aval
bval = B_sol.evalf(subs=model_params)
print("B:", bval)
model_params[B] = bval
rho_sol = sympy.solve(eq2, rho_star)[0] # rho_star for F' = 0
#~ print(rho_sol)
rho_star_val = rho_sol.evalf(subs=model_params)
print("rho_star:", rho_star_val)

# graphical verification
rho = np.linspace(0., model_params[rho_max], 100)
plt.plot(rho, model_params[u_star]*rho*(1. - aval*rho - bval*rho**2))
plt.xlabel(r'$\rho$')
plt.ylabel('F')
plt.axvline(x=rho_star_val, c='k')
plt.show()