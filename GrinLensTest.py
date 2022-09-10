import numpy as np

def obtenerGradiente():

    error = 1
    a_var = 0
    b_var = 0
    c_var = 0
    x_var = 0
    dummy_min = -1

    rayox_inicial = -50
    rayox_final = 70

    while (error > 0.05):

        x0 = 0
        xf = 40 + x_var
        y0 = -9
        yf = 9
        dx = 0.01

        x = np.arange(x0, xf, dx)
        y = np.arange(y0, yf, dx)

        X, Y = np.meshgrid(x, y)

        A = 1.00 + a_var
        B = 1.00 + b_var
        C = 0.0 + c_var
        # Mapa de índices de refracción
        N = A / np.cosh((C * Y)) + B
        indice = lambda y: A / np.cosh(C * y) + B

        # Definir a los Rayos
        r_shape = np.arange(rayox_inicial, 120+xf, dx)

        r1 = np.full(r_shape.shape, fill_value=-5.0, dtype="float")
        r2 = np.full(r_shape.shape, fill_value=-4.0, dtype="float")
        r3 = np.full(r_shape.shape, fill_value=-3.0, dtype="float")
        r4 = np.full(r_shape.shape, fill_value=-2.0, dtype="float")
        r5 = np.full(r_shape.shape, fill_value=-1.0, dtype="float")
        r6 = np.full(r_shape.shape, fill_value=0.0, dtype="float")
        r7 = np.full(r_shape.shape, fill_value=1.0, dtype="float")
        r8 = np.full(r_shape.shape, fill_value=2.0, dtype="float")

        pend_inicial = 0.15

        r1 = r1 + r_shape * pend_inicial
        r2 = r2 + r_shape * pend_inicial
        r3 = r3 + r_shape * pend_inicial
        r4 = r4 + r_shape * pend_inicial
        r5 = r5 + r_shape * pend_inicial
        r6 = r6 + r_shape * pend_inicial
        r7 = r7 + r_shape * pend_inicial
        r8 = r8 + r_shape * pend_inicial

        # Refractar los rayos
        def refractar(r):

            r_nueva = np.zeros(r.shape)

            for k in range(len(r)):
                if x0 > r_shape[k]:
                    r_nueva[k] = r[k]
                if x0 == round(r_shape[k], 2):
                    pendiente = (r_nueva[k] - r_nueva[k - 1]) / dx
                    k1 = np.tan(np.arcsin((1 * np.sin(np.arctan(pendiente))) / indice(r_nueva[k - 1])))
                    k2 = np.tan(np.arcsin((1 * np.sin(np.arctan(pendiente))) / indice(r_nueva[k - 1] + k1 * dx / 2)))
                    k3 = np.tan(np.arcsin((1 * np.sin(np.arctan(pendiente))) / indice(r_nueva[k - 1] + k2 * dx / 2)))
                    k4 = np.tan(np.arcsin((1 * np.sin(np.arctan(pendiente))) / indice(r_nueva[k - 1] + k3 * dx)))
                    r_nueva[k] = r_nueva[k - 1] + (dx / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

                if x0 < r_shape[k] and xf > r_shape[k]:
                    pendiente = (r_nueva[k - 1] - r_nueva[k - 2]) / dx
                    k1 = np.tan(np.arcsin((indice(r_nueva[k - 1]) * np.sin(np.arctan(pendiente))) / indice(r_nueva[k - 1])))
                    k2 = np.tan(np.arcsin(
                        (indice(r_nueva[k - 1]) * np.sin(np.arctan(pendiente))) / indice(r_nueva[k - 1] + k1 * dx / 2)))
                    k3 = np.tan(np.arcsin(
                        (indice(r_nueva[k - 1]) * np.sin(np.arctan(pendiente))) / indice(r_nueva[k - 1] + k2 * dx / 2)))
                    k4 = np.tan(np.arcsin(
                        (indice(r_nueva[k - 1]) * np.sin(np.arctan(pendiente))) / indice(r_nueva[k - 1] + k3 * dx)))
                    r_nueva[k] = r_nueva[k - 1] + (dx / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

                    # r_nueva[k] = r_nueva[k-1] + dx*np.tan(np.arcsin((indice(r_nueva[k-1]) * np.sin(np.arctan(pendiente))) / indice(r_nueva[k-1] + pendiente*dx)))

                if xf == round(r_shape[k], 2):
                    r_nueva[k] = r_nueva[k - 1] + dx * np.tan(
                        np.arcsin((indice(r_nueva[k - 1]) * np.sin(np.arctan(pendiente))) / 1))
                    r_nueva[k + 1] = r_nueva[k] + dx * np.tan(
                        np.arcsin((indice(r_nueva[k]) * np.sin(np.arctan(pendiente))) / 1))
                    pendiente = (r_nueva[k + 1] - r_nueva[k]) / dx

                if xf < r_shape[k]:
                    r_nueva[k] = r_nueva[k - 1] + pendiente * dx

            return r_nueva

        r1 = refractar(r1)
        r2 = refractar(r2)
        r3 = refractar(r3)
        r4 = refractar(r4)
        r5 = refractar(r5)
        r6 = refractar(r6)
        r7 = refractar(r7)
        r8 = refractar(r8)

        #1) Diferencia entre los radios
        rs = abs(r1-r2) + abs(r1-r3) + abs(r1-r4) + abs(r1-r5) + abs(r1-r6) + abs(r1-r7) + abs(r1-r8)

        #2) Verificar que sucede un cruce en una distancia de 70 micrómetros
        posicion_del_foco = False
        if abs(rs[int((abs(rayox_inicial) + xf + 70) / dx )]) < 5:
            posicion_del_foco = True

        #3) Verificar que no existan cruces previo a la intersección
        cruce_previo = False
        for i in range(1, int((abs(rayox_inicial) + xf + 70)/dx) - 1):
            if (r8[i] < r7[i]):
                cruce_previo = True
                break
            if (r7[i] < r6[i]):
                cruce_previo = True
                break
            if (r6[i] < r5[i]):
                cruce_previo = True
                break
            if (r5[i] < r4[i]):
                cruce_previo = True
                break
            if (r4[i] < r3[i]):
                cruce_previo = True
                break
            if (r3[i] < r2[i]):
                cruce_previo = True
                break
            if (r2[i] < r1[i]):
                cruce_previo = True
                break

        print("Min", "[A: ", round(A,2) ,"]","[B: ", round(B,2) ,"]","[C: ", round(C,2) ,"]", "[XF: ", int(xf) ,"]", "[SE CRUZA ANTES: ", cruce_previo ,"]:", "[FOCO EN POSICIÓN: ", posicion_del_foco ,"]:", "[MIN: ", min(rs) ,"]")

        # Modificación de Variables a partir del análisis
        if (cruce_previo == False and posicion_del_foco == True):
            error = 0
            return [A,B,round(C,3),xf]

        # Disminución del Gradiente
        up = False

        # Si aumenta, avanzar más rápido
        if (dummy_min < min(rs)):
            up = True

        if (min(rs) > 50 or up == True):
            c_var = c_var + 0.15
        elif (min(rs) > 30):
            c_var = c_var + 0.10
        elif (min(rs) > 20):
            c_var = c_var + 0.08
        elif (min(rs) > 10):
            c_var = c_var + 0.05
        elif (min(rs) > 5):
            c_var = c_var + 0.02
        elif (min(rs) > 3):
            c_var = c_var + 0.01
        else:
            c_var = c_var + 0.01

        # Variación de Parámetros
        if c_var >= 1.5:
            c_var = 0
            x_var = x_var + 1

        if x_var == 100:
            break

        dummy_min = min(rs)

    return -1

print(obtenerGradiente())
print("End")
