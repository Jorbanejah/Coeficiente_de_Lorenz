import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
# 

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

#                           **Definimos los diferentes constantes y resolvemos el sistema de Lorenz**.
#Algunos valores de sigma, beta y rho son los que se suelen usar para el atractor de Lorenz.
sigma = 10 #    10  16  14 
beta = 8/3 #    8/3  4  3   
rho = 28   #    99.86 45.92 35  
X0 = np.array([0.,1.,0.])
t0 = 0
t_final = 200
h = 0.01

t, x = metodo_Rk(lambda x: sistema_lorentz(x, sigma, beta, rho), t0, X0, t_final, h)#Solucionamos el sistema por el metodo RK4 y metemos los resultados en x

#                           **Empezamos la animación**.

fig, ax = plt.subplots(subplot_kw={'projection': '3d'})

lorenz_plt = ax.plot([], [], [], lw=2)[0] #Inicializamos la línea del gráfico 3D
#Ponemos los límites de los ejes
ax.set_xlim(np.min(x[:,0]), np.max(x[:,0]))
ax.set_ylim(np.min(x[:,1]), np.max(x[:,1]))
ax.set_zlim(np.min(x[:,2]), np.max(x[:,2]))

# Cada frame es un paso de tiempo en el sistema de Lorenz
def animacion(frame):

    x_current = x[:frame+1, 0]
    y_current = x[:frame+1, 1]
    z_current = x[:frame+1, 2]

    lorenz_plt.set_data(x_current, y_current)
    lorenz_plt.set_3d_properties(z_current)
    return lorenz_plt,

animation = FuncAnimation(fig, animacion, frames=range(0,len(x),5), interval=5, blit=False)

plt.title('Atractor de Lorenz')
plt.show()