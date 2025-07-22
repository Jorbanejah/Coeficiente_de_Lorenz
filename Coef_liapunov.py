#Primer intento de crear un programa para el cálculo del coeficiente de Liapunov para sistemas caoticos:
#Atractor de Lorentz
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.linalg import solve
from mpl_toolkits import mplot3d
from functools import reduce
from scipy import stats
from sklearn.linear_model import LinearRegression

#Condiciones iniciales, numero de pasos, perturbación.
X01=np.array([0.,1.,0.])#Condicion inicial.
perturbacion=1e-15
X02=np.array([0.,1+perturbacion,0.])#Condicion inicial con una leve perturbación de 1*10^-15
t0=100
t_final=200
h=0.001

def metodo_Rk(f,t0,X0,t_final,h):
        n = round((t_final-t0)/h)
        t=np.linspace(t0,t_final,n+1)
        m = len(X0)
        x = np.zeros((n+1, m), dtype='float64')
        x[0] = X0
        for i in range(n):
            k1 = h * f(x[i])
            k2 = h * f(x[i] + k1/2)
            k3 = h * f(x[i] + k2/2)
            k4 = h * f(x[i] + k3)
            x[i+1] = x[i] + (k1 + 2*k2 + 2*k3 + k4) / float(6)
        return t, x
def sistema_lorentz(x):
        sigma=10
        rho=28
        beta=8/3
        dx_dt=sigma*(x[1]-x[0])
        dy_dt=x[0]*(rho-x[2])-x[1]
        dz_dt=x[0]*x[1]-beta*x[2]
        return np.array([dx_dt,dy_dt,dz_dt])
    
def comparador_matriz(tol,matriz,parada):
        filas, columnas = matriz.shape
        j=0
        while j<columnas:
            for i in range(filas):
                if abs(matriz[i,j]) >=tol:
                    parada[j]=i
                if j==2:
                    return(parada)
                break
        j=j+1
        
def modulo(matriz,t):
        filas=len(t)
        vector=np.zeros(len(t))
        for i in range(filas):
             vector[i]=math.sqrt(matriz[i,0]**2 + matriz[i,1]**2 + matriz[i,2]**2)
        return vector

def liapunov(t0,X01,X02,t_final,h,perturbacion):
    #Cálculo de los puntos recorridos por el atractor y su perturvación
    t,x=metodo_Rk(sistema_lorentz,t0,X01,t_final,h)
    t,y=metodo_Rk(sistema_lorentz,t0,X02,t_final,h)
    #Vemos donde la distancia entre puntos iguales (mismos tiempos) es suficientemente grande 
    # como para notarse
    z=np.zeros((len(x),len(X01)))
    z=(y-x)
    tol=1e-4 #Es un ejemplo, seguramente esta no sea la tolerancia optima para poder sacar el coef.
    parada=np.array([0,0,0])
    parada=comparador_matriz(tol,z,parada) #No entiendo porque no va.
    print(parada)
    # Linealizamos la ecuación que teniamos.
    deltaln=np.zeros(len(t))
    delta=modulo(z,t)
    for i in range(len(t)):
        deltaln[i]=np.log(delta[i])
    #Dado que no todos los puntos sirven, ya que tenemos una zona más o menos estable al final, y otra, al inicio, 
    # que creo!! que se tiene que quitar porque la diferencia no es sufienciente. Por eso,
    # calculo con un algoritmo de medias el punto estable, para el máximo, y cogo el valor 
    # de la parada más alto para el mínimo.

    LV=int(len(t)/1000) #Longitud de mi vector de medias
    media= np.zeros(LV)
    for l in range(LV):
        k=0
        suma=0
    while k<=1000:
        for i in (len(t)-l*1000-1,len(t)-l*1000-1001):
            k=k+1
            suma=suma+deltaln[i]
            if k==1000:
                media[l]=suma/k

    mayor=0 #Para el siguiente algoritmo y ver cual es la diferencia máxima
    valoritermax=0#donde se guarda la iteracion maxima que estamos buscando

    for k in range(LV-1):
        if (media[k]-abs(media[k+1]))>mayor:
            mayor = media[k]-media[k+1]
            valoritermax=len(t)-k*1000
    print(valoritermax)

    #Por último, hacemos una regresión lineal.
    slope1, intercept1, r1, p1, std_err1 = stats.linregress(t[parada[1] : valoritermax], deltaln[parada[1] : valoritermax])
    print('La pendiente es:',slope1)
    print('La ordenada en el origen es: ', intercept1, 'que al hacer el logaritmo neperiano sale: ', np.exp(intercept1))
    print('La correlación de la regresión: ', r1)
    thorizonte=1/slope1*np.log(10**(-14)/perturbacion)
    print(thorizonte)
    x=np.linspace(110,150,1000)
    plt.plot(t[10000 : 50000], deltaln[10000 : 50000])
    plt.plot(x,slope1*x+intercept1)
    plt.axhline(thorizonte)
    plt.xlabel('t')
    plt.ylabel('||ln(delta(t))||')
    plt.grid(True)
    plt.show()
    return slope1

f = liapunov(t0,X01,X02,t_final,h,perturbacion)
print (f)