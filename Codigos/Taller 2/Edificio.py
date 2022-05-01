#%% TALLER 2 ANALISIS ESTRUCTURAL AVANZADO
##
## Franklin Andres Lizarazo Muñoz
##
##-------------------------------------------------------------------------
#%% Importo las librerias

import pandas as pd
import numpy as n
import matplotlib.pyplot as pl

#%% Funciones
def forces(D,x):
    global p, P, j,XY,IJ,ba
    
    gl      = file[f'Portico{x}_g']
    XY      = file[f'Portico{x}_n']
    element = file[f'Portico{x}_e']

    IJ      = element[['I','J']]
    gdl     = pd.DataFrame([3*IJ.I-3, 3*IJ.I-2, 3*IJ.I-1,3*IJ.J-3, 3*IJ.J-2, 3*IJ.J-1],\
                        index=['1','2','3','4','5','6']).transpose()
        
    j  = int(XY.count().unique())
    ba = int(element.count().unique())
    P, p = ba*[None], ba*[None]

    for e in range(ba):
        De   = D[gdl.loc[e+1]]
        p[e] = n.reshape( kep[x-1][e]@Tep[x-1][e]@De, [2,3]).transpose()
        P[e] = n.reshape( Kep[x-1][e]            @De, [2,3]).transpose()
    
def graf(D):
    D = n.reshape(D,[j,3])

    pl.figure()
    XYdef = XY.copy()
    XYdef.X += fac*D[:,0]
    XYdef.Y += fac*D[:,1]

    for e in range(1,ba+1):
        if p[e-1][0][1]>0:
            col='b'
        else:
            col='r'
        xx=XY.loc[IJ.I.loc[e]],XY.loc[IJ.J.loc[e]]
        pl.plot([xx[0][0], xx[1][0]],[xx[0][1], xx[1][1]],color=(0.5,0.5,0.5,.3),linewidth=1,linestyle='dashed')

        xx=XYdef.loc[IJ.I.loc[e]],XYdef.loc[IJ.J.loc[e]]
        pl.plot([xx[0][0], xx[1][0]],[xx[0][1], xx[1][1]],color=col) 

def Excel(p,D,x):
    D=n.reshape(D,[j,3])
    p=n.reshape(p,[ba*3,2])
    D=pd.DataFrame(D,columns=['Desplazamientos x','Desplazamientos en y', 'Giros'])
    p=pd.DataFrame(p,columns=['Nodo1_Local','Nodo2_Local'])
    with pd.ExcelWriter(f'Desplazamientos{x}.xlsx', engine='xlsxwriter') as writer:
        p.to_excel(writer, sheet_name=f'Esfuerzos{x}')
        D.to_excel(writer, sheet_name=f'Deplazamientos{x}')                           

#%% Leo los datos y asigno las variables

n_porticos = 3;
file=pd.read_excel('Portico1.xlsx',sheet_name=None)
datos = file['Datos']
datos.set_index('Dato', inplace=True)

pisos    = int(datos.Valor.loc['Pisos'])
K_ss_Acu = pisos*[None]
K_ps_Acu = pisos*[None]
C_Acu    = pisos*[None]
Kep      = pisos*[None]
kep      = pisos*[None]
Tep      = pisos*[None]

for u in range(1,n_porticos+1):
    
    element = file[f'Portico{u}_e']
    element.set_index('Elemento', inplace=True)
    
    gl = file[f'Portico{u}_g']
    gl.set_index('gdl', inplace=True)
    
    XY = file[f'Portico{u}_n']
    XY.set_index('Nodo', inplace=True)
    
    #%% Geometría
     
    E   = element.E
    A   = element.A
    ba  = int(element.count().unique())
    j   = int(XY.count().unique())
    ngdl = 3*j
    I  = element.In
    
    #%% Defino la longitud de cada barra
    
    length = element.longitud
    
    #%% Defino los angulos de cada elemento
    
    beta   = element.angulo
    
    #%% Defino las cargas en el grado de libertad dado
    
    P      = gl.carga
    
    #%% Se define la topología de los elementos
    
    IJ     = element[['I','J']]
    
    #%% Se define la matriz con los grados de libertad de cada barra
    
    gdl = pd.DataFrame([3*IJ.I-3, 3*IJ.I-2, 3*IJ.I-1,3*IJ.J-3, 3*IJ.J-2, 3*IJ.J-1],\
                      index=['1','2','3','4','5','6']).transpose()
    
    #%% Se calcula la matriz k local
    
    K = n.zeros([ngdl,ngdl])
    k_acu, T_acu, K_acu = ba*[None], ba*[None], ba*[None]
    
    for e in range(1,ba+1):
        Ea           = E.loc[e]*A.loc[e]/length[e]
        Ei, Ei2, Ei3 = E.loc[e]*I.loc[e]/length[e], E.loc[e]*I.loc[e]/length[e]**2, E.loc[e]*I.loc[e]/length[e]**3
        
        eta, mu      = n.cos(n.deg2rad(beta[e])), n.sin(n.deg2rad(beta[e]))
        
        ke = n.array([[ Ea,       0,      0, -Ea,       0,      0],\
                      [  0,  12*Ei3,  6*Ei2,   0, -12*Ei3,  6*Ei2],\
                      [  0,   6*Ei2,   4*Ei,   0,  -6*Ei2,   2*Ei],\
                      [-Ea,       0,      0,  Ea,       0,      0],\
                      [  0, -12*Ei3, -6*Ei2,   0,  12*Ei3, -6*Ei2],\
                      [  0,   6*Ei2,   2*Ei,   0,  -6*Ei2,   4*Ei]])
        k_acu[e-1] = ke
    
        Te = n.array([[eta,  mu,  0,   0,   0, 0],\
                      [-mu, eta,  0,   0,   0, 0],\
                      [  0,   0,  1,   0,   0, 0],\
                      [  0,   0,  0, eta,  mu, 0],\
                      [  0,   0,  0, -mu, eta, 0],\
                      [  0,   0,  0,   0,   0, 1]])
        T_acu[e-1] = Te
    
        Ke               = Te.transpose()@ke@Te
        K_acu[e-1]       = Ke
        ge               = list(gdl.loc[e])
        K[n.ix_(ge,ge)] += Ke
    
    #%% Grados de libertad restringidos y libres
    
    libertad = gl.Libertad
    
    a  = list(libertad[libertad == 1].index-1)
    b  = list(libertad[libertad != 1].index-1)
    
    gp = list(gl.g_principal[gl.g_principal == 1].index)
    s  = n.setdiff1d(b,gp)
    
    Kaa = K[n.ix_(a,a)];  Kab = K[n.ix_(a,b)]
    Kba = K[n.ix_(b,a)];  Kbb = K[n.ix_(b,b)]
    
    Pb  = P[list(libertad[libertad==0].index)]
    
    #%% Se resuelve el sistema de ecuaciones
    
    Db   = n.linalg.solve(Kbb,Pb)
    D    = n.zeros(ngdl)
    D[b] = Db
    
    Pa   = Kab@Db
    
    #%% Fuerzas al interior de cada elemento
    
    P, p = ba*[None], ba*[None]
    
    Kep[u-1] = K_acu
    kep[u-1] = k_acu
    Tep[u-1] = T_acu
    
    for e in range(ba):
        De   = D[gdl.loc[e+1]]
        p[e] = n.reshape( k_acu[e]@T_acu[e]@De, [2,3]).transpose()
        P[e] = n.reshape( K_acu[e]         @De, [2,3]).transpose()
        
    #%% Matriz condensada
    K_pp = K[n.ix_(gp,gp)];
    K_ps = K[n.ix_(gp,s)];
    K_sp = K[n.ix_(s,gp)];
    K_ss = K[n.ix_(s,s)];
    C    = K_pp-K_ps@n.linalg.inv(K_ss)@K_sp
    
    K_ss_Acu[u-1] = K_ss
    K_ps_Acu[u-1] = K_ps
    C_Acu[u-1]    = C
    
    #%% Grafico
    pl.figure(u)
    pl.spy(K)
    fac=5000

#%% Se ensambla el Edificio

A_p = [[0,1,6],[0,1,2],[0,1,-2],[0,1,-6],[1,0,0],[1,0,-4],[1,0,4]]

AA = n.zeros([30,9]);
for k in range(7+1):
    for i in range(3+1):
        ik=(k-1)*3+i-1;
        g_i=[3*i-3, 3*i-2, 3*i-1];
        AA[n.ix_([ik,ik,ik],g_i)] = A_p[k-1][:];

A_1=AA[0:3][:];
A_2=AA[3:6][:];
A_3=AA[6:9][:];
A_4=AA[9:12][:];
A_5=AA[12:15][:];
A_6=AA[15:18][:];
A_7=AA[18:21][:];

#% Se genera la matriz de rigidez de todo el edificio: 

S = A_1.transpose()@C_Acu[0]@A_1 + A_2.transpose()@C_Acu[1]@A_2+\
    A_3.transpose()@C_Acu[0]@A_3 + A_4.transpose()@C_Acu[1]@A_4+\
    A_5.transpose()@C_Acu[2]@A_5 + A_6.transpose()@C_Acu[2]@A_6+\
    A_7.transpose()@C_Acu[2]@A_7

#%%
    
# =============================================================================
#
#  Metodo de la fuerza horizontal equivalente
# 
# =============================================================================

z = datos.Valor.loc['Altura']
zi = datos.Valor.loc['Altura_piso']
peso_losa = datos.Valor.loc['Peso_losa']
Ba = datos.Valor.loc['Base']
L = datos.Valor.loc['Largo']

Aa ,Av, Fa, Fv, I, R = [.25, .25, 1.5, 1.5, 1, 7]
C = 0.048
a = 0.85

# Periodo aproximado de vibración
T  = C*z**a
Tc = 0.48*Av*Fv/Aa/Fa
Tl = 2.4*Fv

# Aceleración espectral
if T < Tc:
    Sa = 2.5*Aa*Fa*T/Aa/Fa
elif T >= Tl:
    Sa = 1.2*Av*Fv*I/R/T
else:
    Sa = 1.2*Av*Fv*Tl*I/R/T**2


if T < 0.5:
    k = 1
elif T >= 2.5:
    k = 2 
else:
    k= 0.75+0.5*T

M = peso_losa*Ba*L/9.81*n.ones(3)
m = M
mz = zi**k*m

h = n.zeros(int(pisos))
for i in range(int(pisos)):
    h[i]=m[i]*zi**k*M[2]*Sa*9.81/mz.sum()

ex, ey = [0.05*Ba, 0.05*L]

B = 0.3*ex*h+ey*h

H = n.zeros(9)
g_x=n.arange(0,7,3)
g_y=n.arange(1,8,3)
g_b=n.arange(2,9,3)

H[g_x] = 0.3*h
H[g_y] = h
H[g_b] = B

U = n.linalg.inv(S)@H

#%% Desplazamientos principales (moviemientos horizontales) de los pórticos:

D_p1 = A_1@U
D_p2 = A_2@U
D_p3 = A_3@U
D_p4 = A_4@U
D_p5 = A_5@U
D_p6 = A_6@U
D_p7 = A_7@U

#%% Desplazamientos secundarios

D_s1 = -n.linalg.inv(K_ss_Acu[0])@(K_ps_Acu[0].transpose()@D_p1)
D_s2 = -n.linalg.inv(K_ss_Acu[1])@(K_ps_Acu[1].transpose()@D_p2)
D_s3 = -n.linalg.inv(K_ss_Acu[0])@(K_ps_Acu[0].transpose()@D_p3)
D_s4 = -n.linalg.inv(K_ss_Acu[1])@(K_ps_Acu[1].transpose()@D_p4)
D_s5 = -n.linalg.inv(K_ss_Acu[2])@(K_ps_Acu[2].transpose()@D_p5)
D_s6 = -n.linalg.inv(K_ss_Acu[2])@(K_ps_Acu[2].transpose()@D_p6)
D_s7 = -n.linalg.inv(K_ss_Acu[2])@(K_ps_Acu[2].transpose()@D_p7)
             
#%% Desplazamientos totales de los porticos

gl = file['Portico1_g']
b  = list(gl.Libertad[gl.Libertad != 1].index-1)
gp = list(gl.g_principal[gl.g_principal == 1].index)
s  = n.setdiff1d(b,gp)

#Portico 1
D1 = n.zeros(36)
D1[gp] = D_p1
D1[s]  = D_s1

forces(D1,1)
graf(D1)
Excel(p,D1,1)

#Portico 3
D3 = n.zeros(36)
D3[gp] = D_p3
D3[s]  = D_s3

forces(D3,1)
graf(D3)
Excel(p,D3,3)

gl = file['Portico2_g']
b  = list(gl.Libertad[gl.Libertad != 1].index-1)
gp = list(gl.g_principal[gl.g_principal == 1].index)
s  = n.setdiff1d(b,gp) 

#Portico 2
D2 = n.zeros(36)
D2[gp] = D_p2
D2[s]  = D_s2

forces(D2,2)
graf(D2)
Excel(p,D2,2)

#Portico 4
D4 = n.zeros(36)
D4[gp] = D_p4
D4[s]  = D_s4

forces(D4,2)
graf(D4)
Excel(p,D4,4)

gl = file['Portico3_g']
b  = list(gl.Libertad[gl.Libertad != 1].index-1)
gp = list(gl.g_principal[gl.g_principal == 1].index)
s  = n.setdiff1d(b,gp)

#Portico 5
D5 = n.zeros(48)
D5[gp] = D_p5
D5[s]  = D_s5

forces(D5,3)
graf(D5)
Excel(p,D5,5)

#Portico 6
D6 = n.zeros(48)
D6[gp] = D_p6
D6[s]  = D_s6

forces(D6,3)
graf(D6)
Excel(p,D6,6)

#Portico 7
D7 = n.zeros(48)
D7[gp] = D_p7
D7[s]  = D_s7

forces(D7,3)
graf(D7)
Excel(p,D7,7)

#%%

# =============================================================================
#
# Metodo de Dinamico
#
# =============================================================================

mr = m*(8^2+12^2)/12
M  = sum(m)

T  = 0.047*z**0.9
Masa=n.diag([m[0], m[0], mr[0], m[1], m[1], mr[1], m[2], m[2], mr[2]])

Phi, Omega = [file['Phi'].to_numpy(), file['Omega'].to_numpy()]

T  = 2*n.pi/Omega
T0 = 0.1*Av*Fv/(Aa*Fa)
Tc = 0.48*Av*Fv/(Aa*Fa)
Tl = 2.4*Fv
Sa = n.zeros(n.shape(T))

for i in range(9):
    if T[i]<=T0:
        Sa[i] = 2.5*Aa*Fa*I*(0.4+0.6*T[i]/T0)/R
    elif T[i] <= Tc:
        Sa[i] = 2.5*Aa*Fa*I/R
    elif T[i] > Tl:
        Sa[i] = 1.2*Av*Fv*Tl*I/(R*T[i]**2)
    else:
        Sa[i] = 1.2*Av*Fv*I/(R*T[i])

Sa *= 9.81

rx      = n.zeros(9);
ry      = n.zeros(9);
rx[g_x] = 1; 
ry[g_y] = 0.3;
r       = rx + ry; 

F = n.zeros(n.shape(Phi)) 
q = n.zeros(n.shape(T)) 

for i in range(9):
    q[i]    = Phi[:][i].transpose()@Masa@r/(Phi[:][i].transpose()@Masa@Phi[:][i])
    F[:][i] = n.reshape(n.reshape(Masa@Phi[:][i],[9,1])@q[i],[9,1])@Sa[i]

#%% Fuerzas finales - Método RCSC (NSR-10 A.5.4.4)

Fsrss=n.zeros([9,1])
for k in range(9):
    for i in range(9):
        Fsrss[k] = Fsrss[k] + F[k][i]**2

Fsrss = n.sqrt(Fsrss); #Fuerzas arrojadas por el método dinámico

Fx          = Fsrss[g_x]                    # Fuerzas dinámicas en X
Fy          = Fsrss[g_y]                    # Fuerzas dinámicas en Y
Bsrss       = Fsrss[g_b]                    # Pares torsores
Ba_srss     = Fx*ey + Fy*ex               # Torsión accidental adicional
Btotal_srss = Bsrss + Ba_srss
Facc        = n.zeros([9,1])
Facc[g_b]   = Ba_srss                       
Fsrss      += Facc

Usrss = n.linalg.inv(S)@Fsrss;              # Desplazamientos, método RCSC
 

#%% Desplazamientos principales (moviemientos horizontales) de los pórticos:

D_psrss1 = A_1@Usrss;
D_psrss2 = A_2@Usrss;
D_psrss3 = A_3@Usrss;
D_psrss4 = A_4@Usrss;
D_psrss5 = A_5@Usrss;
D_psrss6 = A_6@Usrss;
D_psrss7 = A_7@Usrss;

#%% Desplazamientos secundarios:

D_ssrss1 = -n.linalg.inv(K_ss_Acu[0])@(K_ps_Acu[0].transpose()@D_psrss1);
D_ssrss2 = -n.linalg.inv(K_ss_Acu[1])@(K_ps_Acu[1].transpose()@D_psrss2);
D_ssrss3 = -n.linalg.inv(K_ss_Acu[0])@(K_ps_Acu[0].transpose()@D_psrss3);
D_ssrss4 = -n.linalg.inv(K_ss_Acu[1])@(K_ps_Acu[1].transpose()@D_psrss4);
D_ssrss5 = -n.linalg.inv(K_ss_Acu[2])@(K_ps_Acu[2].transpose()@D_psrss5);
D_ssrss6 = -n.linalg.inv(K_ss_Acu[2])@(K_ps_Acu[2].transpose()@D_psrss6);
D_ssrss7 = -n.linalg.inv(K_ss_Acu[2])@(K_ps_Acu[2].transpose()@D_psrss7);

#%%Vector de desplazamientos de cada portico 
 
gl = file['Portico1_g']
bsrss  = list(gl.Libertad[gl.Libertad != 1].index-1)
psrss  = list(gl.g_principal[gl.g_principal == 1].index)
ssrss  = n.setdiff1d(bsrss,psrss)

#Portico 1 
D_1srss           = n.zeros([36,1])
D_1srss[psrss]    = D_psrss1
D_1srss[ssrss]    = D_ssrss1

forces(D_1srss,1)
graf(D_1srss)
Excel(p,D_1srss,11)
                        
#Portico 3
D_3srss           = n.zeros([36,1])
D_3srss[psrss]    = D_psrss3
D_3srss[ssrss]    = D_ssrss3

forces(D_3srss,1)
graf(D_3srss)
Excel(p,D_3srss,33)
 

gl = file['Portico2_g']
bsrss  = list(gl.Libertad[gl.Libertad != 1].index-1)
psrss  = list(gl.g_principal[gl.g_principal == 1].index)
ssrss  = n.setdiff1d(bsrss,psrss)

#Portico 2
D_2srss           = n.zeros([36,1])
D_2srss[psrss]    = D_psrss2
D_2srss[ssrss]    = D_ssrss2

forces(D_2srss,2)
graf(D_2srss)
Excel(p,D_2srss,22)
                       
#Portico 4
D_4srss           = n.zeros([36,1])
D_4srss[psrss]    = D_psrss4
D_4srss[ssrss]    = D_ssrss4

forces(D_4srss,2)
graf(D_4srss)
Excel(p,D_4srss,44)

gl = file['Portico3_g']
bsrss  = list(gl.Libertad[gl.Libertad != 1].index-1)
psrss  = list(gl.g_principal[gl.g_principal == 1].index)
ssrss  = n.setdiff1d(bsrss,psrss)

#Portico 5
D_5srss           = n.zeros([48,1])
D_5srss[psrss]    = D_psrss5
D_5srss[ssrss]    = D_ssrss5

forces(D_5srss,3)
graf(D_5srss)
Excel(p,D_5srss,55)
                       
#Portico 6
D_6srss           = n.zeros([48,1])
D_6srss[psrss]    = D_psrss6
D_6srss[ssrss]    = D_ssrss6

forces(D_6srss,3)
graf(D_6srss)
Excel(p,D_6srss,66)

#Portico 7
D_7srss           = n.zeros([48,1])
D_7srss[psrss]    = D_psrss7
D_7srss[ssrss]    = D_ssrss7

forces(D_7srss,3)
graf(D_7srss)
Excel(p,D_7srss,77)

Resultado=pd.DataFrame(n.concatenate((n.reshape(H,[9,1]),Fsrss),axis=1),columns=['Metodo Fuerza Equivalente','Metodo Dinamico'])
Resultado.to_excel('Resultado.xlsx')
