#Programa_1: ver como cambia una función de la forma dx/dt=-ax____ Aproximando la derivada por un euler explicito
import sys
import numpy as np
import matplotlib.pyplot as plt

def f(t, x, a):
    return a*x-x**3


#Defino el método de euler explicito donde f es la función a la que esta igualada la derivada, x0 el valor inicial, contiene los valores del tiempo y
# n número de pasos:
def euler_explicito(x0,t,f,a):
    n=len(t)#número de pasos
    dt = t[1] - t[0] # Tamaño del paso de tiempo
    x = np.zeros(n) #Almacenar la solución
    x[0] = x0 # Establecezco el valor inicial

    for i in range(1,n):
        x[i] = x[i-1] + f(t[i-1],x[i-1], a)*dt # Método de Euler explícito
    return x

t = np.linspace(0, 5, 1000) # Tiempos de 0 a 5, divididos en 100 pasos
# Llamo a la función:
a=float(input("Que constante 'a' quieres? Positiva, negativa o cero (escribela): "))

# a=int(a)
# if (a==0):
#     x0=input("Cual es la condición inicial del programa: ")
#     def f(y,t):
#         return 0
# elif(a<0):
#     x0=input("Cual es la condición inicial del programa: ")
#     def f(y,t):
#         return -y
# elif (a>0):
#     x0=input("Cual es la condición inicial del programa: ")
#     def f(y,t):
#         return y
x0=float(input("Cual es la condición inicial del programa: "))



solucion1 = euler_explicito(x0, t, f, 1)
solucion2 = euler_explicito(-x0, t, f, 1)


# Graficar la solución
plt.plot(t, solucion1, label='Solución aproximada')
plt.plot(t, solucion2, label='Solución aproximada')
plt.axhline(y=0)
plt.xlabel('Tiempo')
plt.ylabel('Valor de y')
plt.title('Método de Euler Explícito')
plt.legend()
plt.grid(True)
plt.show()
