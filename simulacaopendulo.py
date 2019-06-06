from numpy import *
import matplotlib.pyplot as plt
from vpython import *
seterr(all='raise')

def RungeKutta():
    global flag
    flag         = False
    variaveis[0] = array([r0, rponto0, -theta0, -thetaponto0])  #Theta negativo para ficar de acordo com o desenho do problema
    canvas(visible=True,center=vector(1,-2,0))
    m1 = sphere(radius = 0.1, color = vector(1,0.752,0.329)) #Massa que oscila
    corda1=cylinder(pos=vector(0,0,0),axis=m1.pos,radius=0.01,color=vector(1,1,1))
    polia1 = ring(pos=vector(0.1,0,0),axis=vector(0,0,1),radius=0.1,thickness=0.01, color=0.5*vector(1,1,1))
    cylinder(pos=polia1.pos+vector(0,polia1.radius,0),axis=vector(2,0,0),radius=0.01,color=corda1.color)
    polia2 = ring(pos=polia1.pos+vector(2,0,0),axis=vector(0,0,1),radius=0.1,thickness=0.01, color=0.5*vector(1,1,1))
    n=6 #Número de raias na polia
    for i in range(0,n):
        cylinder(pos=polia1.pos,axis=vector(polia1.radius*cos(i*2*pi/n),polia1.radius*sin(i*2*pi/n),0),radius=0.01,color=0.5*vector(1,1,1))
        cylinder(pos=polia2.pos,axis=vector(polia2.radius*cos(i*2*pi/n),polia2.radius*sin(i*2*pi/n),0),radius=0.01,color=0.5*vector(1,1,1))
    m2 = sphere(pos = polia2.pos+vector(polia2.radius,-10*r0,0),radius = 0.1, color = m1.color) #Massa que só se move na vertical
    corda2=cylinder(pos=polia2.pos+vector(polia2.radius,0,0),axis=m2.pos-polia2.pos-vector(polia2.radius,0,0),radius=0.01,color=corda1.color)
    trajetoria = curve(color=color.cyan)
    for i in range(indices - 1):
        rate(500)
        try:
            if(((theta0 == 0 and thetaponto0 == 0) or theta0 == pi) and i == 0):
                print('Aviso: Equilíbrio estável ou instável escolhido como valor inicial.')
            if(variaveis[i][0] <= (limInf * r0)):
                break
            m1.pos = variaveis[i][0]*vector(sin(variaveis[i][2]), -cos(variaveis[i][2]),0 )
            corda1.axis=m1.pos
            m2.pos.y = m1.pos.mag - 11*r0
            corda2.axis=m2.pos-polia2.pos-vector(polia2.radius,0,0)
            trajetoria.append(m1.pos)
            k1 = precisao * EqMovimento(variaveis[i], tempos[i])
            k2 = precisao * EqMovimento(variaveis[i] + k1 / 2, tempos[i] + precisao/2)
            k3 = precisao * EqMovimento(variaveis[i] + k2 / 2, tempos[i] + precisao/2)
            k4 = precisao * EqMovimento(variaveis[i] + k3, tempos[i] + precisao)

            variaveis[i + 1] = variaveis[i] + k1/6 + k2/3 + k3/3 + k4/6
        except FloatingPointError:
            flag = True
            break

def EqMovimento(variaveis, tempos):
    r    = variaveis[0]
    rponto = variaveis[1]
    theta     = variaveis[2]
    thetaponto  = variaveis[3]

    rpontoponto = ((r / (1 + mu)) * ((thetaponto) ** 2)) + (((g * cos(theta)) - (g * mu)) / (1 + mu))
    thetapontoponto  = - ((g * sin(theta)) / r) - (2 * ((rponto) * (thetaponto)) / r)

    return array([rponto, rpontoponto, thetaponto, thetapontoponto])

#Variáveis do algorítmo.
precisao      = 0.001
tmax   = 5
limInf = 0.01

indices = int(tmax / precisao)
tempos   = linspace(0, (indices - 1) * precisao, indices)

#Variáveis físicas.
r0         = 0.5
rponto0      = 0
theta0     = -pi/4
thetaponto0  = 0
g          = 9.8
mu         = 2.5 #Razão das massas m2/m1

#Variáveis do algorítmo de Runge-Kutta
variaveis = zeros([indices, 4], dtype=float)
RungeKutta()


#Construir os gráficos
plt.figure()
ax = plt.subplot(111, projection='polar')
ax.set_theta_zero_location("S")
ax.plot(variaveis[:,2], variaveis[:,0], color='r', linewidth=1)
plt.title('Oscilações para $\mu = %.3f$, $r_0 = %.3f$ $m$, $\\theta_0 = %.3f^\degree$, $t_{max} = %.2f$ $s$' % (mu, r0, theta0 * 180 / pi, tmax), y = 1.06)
ax.grid(True)
if(flag == False):
    print('Simulacao concluida com sucesso.')
    plt.show()
