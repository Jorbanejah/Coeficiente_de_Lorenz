import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter

#                           **Definiciones(RK4 y sistema de Lorenz)**
def metodo_Rk(f,t0,X0,t_final,h):

    n = round((t_final-t0)/h)
    t=np.linspace(t0,t_final,n+1)
    m = len(X0)
    x = np.zeros((n+1, m))
    x[0] = X0
    for i in range(n):
        k1 = h * f(x[i])
        k2 = h * f(x[i] + k1/2)
        k3 = h * f(x[i] + k2/2)
        k4 = h * f(x[i] + k3)
        x[i+1] = x[i] + (k1 + 2*k2 + 2*k3 + k4) / float(6)
    return t, x

def sistema_lorentz(x,sigma,beta,rho):
 
    dx_dt=sigma*(x[1]-x[0])
    dy_dt=x[0]*(rho-x[2])-x[1]
    dz_dt=x[0]*x[1]-beta*x[2]
    return np.array([dx_dt,dy_dt,dz_dt])

#                           **Condiciones iniciales, parámetros y pertubaciones**
sigma = 10 #    10  16  14 
beta = 8/3 #    8/3  4  3   
rho = 28   #    99.86 45.92 35  
X0 = np.array([0.,1.,0.])
X01 = np.array([0.,1.,0.01])
t0 = 0
t_final = 200
h = 0.01

# Cálculo del sistema de Lorenz
t, x = metodo_Rk(lambda x: sistema_lorentz(x, sigma, beta, rho), t0, X0, t_final, h)
t1, x1 = metodo_Rk(lambda x: sistema_lorentz(x, sigma, beta, rho), t0, X01, t_final, h)

#Exportamos las trayectorias a archivos CSV para utilizarlas en Coef_liapunov.py
np.savetxt("trayectoria_x_con_t.csv", np.column_stack((t, x)), delimiter=",", header="t,X,Y,Z", comments='')
np.savetxt("trayectoria_x1_con_t.csv", np.column_stack((t1, x1)), delimiter=",", header="t,X,Y,Z", comments='')


#                           **Visualización de los resultados**

fig, ax = plt.subplots(subplot_kw={'projection': '3d'})

lorenz_plt = ax.plot([], [], [], lw=2, color = 'red')[0] 
lorenz_plt1 = ax.plot([], [], [], lw=2, color='orange')[0] 

ax.set_xlim([-20, 20])
ax.set_ylim([-30, 30])
ax.set_zlim([0, 50])

# Cada frame es un paso de tiempo en el sistema de Lorenz
def animacion(frame):

    # Trayectoria principal
    x_current = x[:frame+1, 0]
    y_current = x[:frame+1, 1]
    z_current = x[:frame+1, 2]
    lorenz_plt.set_data(x_current, y_current)
    lorenz_plt.set_3d_properties(z_current)

    # Trayectoria perturbada
    x_current1 = x1[:frame+1, 0]
    y_current1 = x1[:frame+1, 1]
    z_current1 = x1[:frame+1, 2]
    lorenz_plt1.set_data(x_current1, y_current1)
    lorenz_plt1.set_3d_properties(z_current1)
   
    return lorenz_plt, lorenz_plt1

animation = FuncAnimation(fig, animacion, frames=range(2,len(x),5), interval=5, blit=False)

plt.title('Atractor de Lorenz')
plt.show()