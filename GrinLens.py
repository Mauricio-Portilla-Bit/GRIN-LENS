import matplotlib.pyplot as plt
import numpy as np

# [A:  1.30 ] [B:  1 ]
# xf = 47 / C = 0.24
# xf = 57 / C = 0.78
# xf = 80 / C = 0.64
# xf = 45 / C = 0.39 / x0 = 0
# [A:  1.66 ] [B:  1.4 ] [C:  0.27 ] [XF:  70 ]
# Min [A:  2.0 ] [B:  1.0 ] [C:  0.28 ] [XF:  40 ] [SE CRUZA ANTES:  False ]: [FOCO EN POSICIÓN:  False ]: [MIN:  1.5721012970809944 ]
# Min [A:  1.0 ] [B:  1.0 ] [C:  0.29 ] [XF:  40 ] [SE CRUZA ANTES:  False ]: [FOCO EN POSICIÓN:  False ]: [MIN:  2.2787364547591054 ]


x0 = 0
xf = 40
y0 = -10
yf = 10
dx = 0.01

x = np.arange(x0,xf,dx)
y = np.arange(y0,yf,dx)
X,Y = np.meshgrid(x,y)

A = 1
B = 1
C = 0.29
# Mapa de índices de refracción
N = A/np.cosh((C*Y)) + B
indice = lambda y: A/np.cosh(C*y) + B

# Definir a los Rayos
r_shape = np.arange(-50,80+xf,dx)

r1 = np.full(r_shape.shape, fill_value= -5.0, dtype="float")
r2 = np.full(r_shape.shape, fill_value= -4.0, dtype="float")
r3 = np.full(r_shape.shape, fill_value= -3.0, dtype="float")
r4 = np.full(r_shape.shape, fill_value= -2.0, dtype="float")
r5 = np.full(r_shape.shape, fill_value= -1.0, dtype="float")
r6 = np.full(r_shape.shape, fill_value= 0.0, dtype="float")

pend_inicial = 0.15

r1 = r1 + r_shape*pend_inicial
r2 = r2 + r_shape*pend_inicial
r3 = r3 + r_shape*pend_inicial
r4 = r4 + r_shape*pend_inicial
r5 = r5 + r_shape*pend_inicial
r6 = r6 + r_shape*pend_inicial

# Refractar los rayos
def refractar (r):

    r_nueva = np.zeros(r.shape)

    for k in range(len(r)):

        # Antes de la intersección
        if x0 > r_shape[k]:
            r_nueva[k] = r[k]

        # En la primera intersección
        if x0 == round(r_shape[k], 2):
            pendiente = (r_nueva[k] - r_nueva[k-1])/dx
            k1 = np.tan(np.arcsin((1 * np.sin(np.arctan(pendiente))) / indice(r_nueva[k-1])))
            k2 = np.tan(np.arcsin((1 * np.sin(np.arctan(pendiente))) / indice(r_nueva[k-1] + k1*dx/2)))
            k3 = np.tan(np.arcsin((1 * np.sin(np.arctan(pendiente))) / indice(r_nueva[k-1] + k2*dx/2)))
            k4 = np.tan(np.arcsin((1 * np.sin(np.arctan(pendiente))) / indice(r_nueva[k-1] + k3*dx)))
            r_nueva[k] = r_nueva[k-1] + (dx/6)*(k1 + 2*k2 + 2*k3 + k4)

        # En el lente de GRIN
        if x0 < r_shape[k] and xf > r_shape[k]:
            pendiente = (r_nueva[k-1] - r_nueva[k-2])/dx
            k1 = np.tan(np.arcsin((indice(r_nueva[k-1]) * np.sin(np.arctan(pendiente))) / indice(r_nueva[k-1])))
            k2 = np.tan(np.arcsin((indice(r_nueva[k-1]) * np.sin(np.arctan(pendiente))) / indice(r_nueva[k-1] + k1*dx/2)))
            k3 = np.tan(np.arcsin((indice(r_nueva[k-1]) * np.sin(np.arctan(pendiente))) / indice(r_nueva[k-1] + k2*dx/2)))
            k4 = np.tan(np.arcsin((indice(r_nueva[k-1]) * np.sin(np.arctan(pendiente))) / indice(r_nueva[k-1] + k3*dx)))
            r_nueva[k] = r_nueva[k-1] + (dx/6)*(k1 + 2*k2 + 2*k3 + k4)

        # En la salida del lente
        if xf == round(r_shape[k], 2):
            r_nueva[k] = r_nueva[k-1] + dx*np.tan(np.arcsin((indice(r_nueva[k-1]) * np.sin(np.arctan(pendiente))) / 1))
            r_nueva[k+1] = r_nueva[k] + dx*np.tan(np.arcsin((indice(r_nueva[k]) * np.sin(np.arctan(pendiente))) / 1))
            pendiente = (r_nueva[k+1] - r_nueva[k])/dx

        # Después de la intersección
        if xf < r_shape[k]:
            r_nueva[k] = r_nueva[k-1] + pendiente*dx

    return r_nueva

r1 = refractar(r1)
r2 = refractar(r2)
r3 = refractar(r3)
r4 = refractar(r4)
r5 = refractar(r5)
r6 = refractar(r6)

plt.plot(r6, r_shape, label='label6')
plt.plot(r5, r_shape, label='label5')
plt.plot(r4, r_shape, label='label4')
plt.plot(r3, r_shape, label='label3')
plt.plot(r2, r_shape, label='label2')
plt.plot(r1, r_shape, label='label1')

plt.title("GRIN LENS (N = SECH(CY))")
plt.xlabel("X")
plt.ylabel("Y")

plt.pcolormesh(Y,X,N, cmap="plasma")
#plt.pcolormesh(X,Y,N, cmap="plasma")
cax = plt.axes([0.905, 0.1, 0.02, 0.8])
plt.colorbar(cax=cax)
plt.show()


