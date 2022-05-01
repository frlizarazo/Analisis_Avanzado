%% TALLER 1 ANALISIS ESTRUCTURAL AVANZADO
%
% Franklin Andres Lizarazo Muñoz
%
%-------------------------------------------------------------------------

%% Geometría General
E  = 2e5;              % kPa
A  = 0.005;            % m^2
ba = 21;               % Numero de barras
j  = 12;               % Numero de nodos
I  = (0.05*0.1^3)/12;  % m^4 Inercia

% Se define la longitud de cada barra manualmente [en m]
le = [5 5 5 5 5 5 6.4 5.5 5.04 5.04 5.5 6.4 4 6.4 6.3 8.04 7 8.04 6.3 6.4 4];

% Se define la inclinación de cada barra manualmente [en °]
beta  = [0 0 0 0 0 0 39 25 8 -8 -25 -39 90 -39 90 -52 90 52 90 39 -90];

% Se define la topología de los elementos
IJ = [ 1  2;  2  3;  3  4; 4  5; 5  6; 6  7; 1 12; 12 11; 11 10; 10  9; 
       9  8; 8  7; 2 12; 12 3; 3 11; 11  4; 4 10; 4  9; 5  9; 5  8; 8  6];

% Se define la posición de cada uno de los nodos
XY = [ 0   0;  5   0; 10   0; 15   0; 20   0; 25   0; 30   0; 25   4;
      20 6.3; 15   7; 10 6.3;  5   4];

%% ------------------------------------------------------------------------
% Uniones Articuladas
%--------------------------------------------------------------------------
%% Cargas y Grados de Libertad

% Se ubican las cargas en los nodos especificados en la dirección dada[en kN]
%     x1 y1 ...
P  = [0  0  0 -2 0 -2 0 -2 0 -2 0 -2 0 0 0 -5 0 -10 0 -15 0 -10 0 -5]';

% Se define la matriz con los grados de libertad de cada barra
gdl=zeros(ba,4);
for e=1:ba
   gdl(e,:)=[2*IJ(e,1)-1 2*IJ(e,1) 2*IJ(e,2)-1 2*IJ(e,2)];
end

%% Calcular la matriz k local

% Se definen unos acumuladores para usar datos mas adelante
[eta_acu, mu_acu, K_acu, k_acu, T_acu] = deal(cell(1,ba));

% Se separa el espacio de la matriz K global
K = zeros(j*2);

% Se calculan los ke locales, se guardan en los acumuladores y se ensabla
% la matriz global K
for e=1:ba
    ke   = E*A/le(e)*[ 1 0 -1 0
                       0 0  0 0
                      -1 0  1 0
                       0 0  0 0]; k_acu{e} = ke;
                 
    %Se extrae el angulo del elemento e del vector de angulos
    eta  = cosd(beta(e)); eta_acu{e} = eta;      
    mu   = sind(beta(e)); mu_acu{e}  = mu;
    
    %Se define la matriz de transformación
    Te   = [ eta  mu   0   0
             -mu eta   0   0
               0   0 eta  mu
               0   0 -mu eta]; T_acu{e} = Te;
           
    Ke   = Te'*ke*Te; K_acu{e} = Ke;
    ge   = gdl(e,:);
    
    K(ge,ge) = K(ge,ge)+Ke;   
end

%% Grados de Libertad Restringidos y libres
a = [1 2 13 14];   b = setdiff(1:2*j,a);

Kaa = K(a,a);     Kab = K(a,b);
Kba = K(b,a);     Kbb = K(b,b);

Pb  = P(b);

% Se resuelve el sistema de ecuaciones
Db   = Kbb\Pb;
D    = zeros(j*2,1);
D(b) = Db;

Pa   = Kab*Db;
P1   = zeros(6,1);    P1([1,2,4,5]) = Pa; 

%% Tensiones en cada elemento
sigma = zeros(1,ba);

for x = 1:ba
    ge       = D(gdl(x,:));
    sigmae   = E*[-eta_acu{x} -mu_acu{x} eta_acu{x} mu_acu{x}]*ge/le(x);
    sigma(x) = sigma(x)+sigmae;
end

%% Fuerzas al interior de cada elemento
[P, p]  = deal(cell(ba,1));

for e = 1:ba
   De   = D(gdl(e,:));
   p{e} = k_acu{e}*T_acu{e}*De;
   P{e} = K_acu{e}*De;
end

%% Graficar
 
XYdef = zeros(size(XY));
fac   = 1;
c     = 0;
for i = 1:j
    c          = c+1;
    XYdef(i,1) = XY(i,1)+fac*D(c);
    c          = c+1;
    XYdef(i,2) = XY(i,2)+fac*D(c);
end

% Se grafican la estructura y su deformada
figure
for e = 1:ba
    Q    = [XY(IJ(e,1),1)     XY(IJ(e,1),2);...
            XY(IJ(e,2),1)     XY(IJ(e,2),2)];
        
    Qdef = [XYdef(IJ(e,1),1)  XYdef(IJ(e,1),2);...
            XYdef(IJ(e,2),1)  XYdef(IJ(e,2),2)];
        
    if sigma(e)>0
        plot(Q(:,1),Q(:,2),'k--',Qdef(:,1),Qdef(:,2),'b')
    else
        plot(Q(:,1),Q(:,2),'k--',Qdef(:,1),Qdef(:,2),'r')
    end
    hold on
end
title('Estructura y Deformada: Uniones Articuladas')
xlabel('x [m]')
ylabel('y [m]  factor = 1')
axis equal

% Se grafican las fuerzas axiales
figure
for e = 1:ba
    subplot(3,7,e)
    plot([0 le(e)],[0 0],'lineWidth',5)
    title('Barra',e)
    hold on
    plot([0 0 le(e) le(e)],[0 p{e}(3) p{e}(3) 0])
    text(le(e)/5,   p{e}(3)/2,   num2str(p{e}(3)),'FontSize',18);
end
sgtitle('Fuerza Axial: Uniones Articuladas')

%% Extraigo tablas resumen
bar = 1:ba;
n   = 1:j;

[Axial, Cortante] = deal(zeros(ba,1));

for e = 1:ba
    Axial(e)    = p{e}(3);
    Cortante(e) = p{e}(2);
end

disp(' ')
disp('Resumen Uniones Articuladas-----------------------------------------')

T1_1  = table(bar', sigma', Axial, Cortante, zeros(21,1),...
       'VariableNames',{'Barra','Esfuerzo','F Axial','F Cortante','Momento'})

[Dx, Dy]  = deal(zeros(j,1));
c         = 0;

for e     = 1:2:2*j
    c     = c+1;
    Dx(c) = D(e);
    Dy(c) = D(e+1);
end

T1_2 = table(n',Dx, Dy,'VariableNames',{'Nodo','Desp. X','Desp. Y'})


%% ------------------------------------------------------------------------
% Uniones Rigidas
% ------------------------------------------------------------------------- 
%% Cargas y Grados de Libertad

% Se ubican las cargas en los nodos especificados en la dirección dada[en kN]
P  = zeros(j*3,1);
P([5 8 11 14 17],1)=-2;  P([23 35],1)=-5;  P([26 32],1)=-10;  P(29,1)=-15;

% Se define la matriz con los grados de libertad de cada barra
gdl=zeros(ba,6);

for e=1:ba
   gdl(e,:)=[3*IJ(e,1)-2 3*IJ(e,1)-1 3*IJ(e,1) 3*IJ(e,2)-2 3*IJ(e,2)-1 3*IJ(e,2)];
end

%% Calcular la matriz k local

% Se definen unos acumuladores para usar datos mas adelante
[eta_acu, mu_acu, K_acu, k_acu, T_acu] = deal(cell(1,ba));

% Se separa el espacio de la matriz K global
K = zeros(j*3);

% Se calculan los ke locales, se guardan en los acumuladores y se ensabla
% la matriz global K
for e  = 1:ba
    
    % Por legibilidad se definen
    Ea = E*A/le(e);    Ei=E*I/le(e)^3;
    ke = [  Ea           0             0  -Ea            0             0
             0       12*Ei    6*le(e)*Ei    0       -12*Ei    6*le(e)*Ei
             0  6*le(e)*Ei  4*le(e)^2*Ei    0  -6*le(e)*Ei  2*le(e)^2*Ei
           -Ea           0             0   Ea            0             0
             0      -12*Ei   -6*le(e)*Ei    0        12*Ei   -6*le(e)*Ei
             0  6*le(e)*Ei  2*le(e)^2*Ei    0  -6*le(e)*Ei  4*le(e)^2*Ei ];   
    k_acu{e} = ke;
    
    %Se extrae el angulo del elemento e del vector de angulos
    eta  = cosd(beta(e)); eta_acu{e} = eta;      
    mu   = sind(beta(e)); mu_acu{e}  = mu;
    
    %Se define la matriz de transformación
    Te   = [ eta  mu   0    0   0   0
             -mu eta   0    0   0   0
               0   0   1    0   0   0
               0   0   0  eta  mu   0
               0   0   0  -mu eta   0
               0   0   0    0   0   1]; T_acu{e} = Te;
       
    Ke   = Te'*ke*Te;K_acu{e} = Ke;
    ge   = gdl(e,:);
    
    K(ge,ge) = K(ge,ge)+Ke;
end

%% Grados de Libertad Restringidos y libres
a = [1 2 3 19 20 21];   b = setdiff(1:3*j,a);

Kaa = K(a,a);     Kab = K(a,b);
Kba = K(b,a);     Kbb = K(b,b);

Pb = P(b);

% Se resuelve el sistema de ecuaciones
Db   = Kbb\Pb;
D    = zeros(j*3,1);
D(b) = Db;

Pa   = Kab*Db;

%% Fuerzas al interior de cada elemento
[P, p]  = deal(cell(ba,1));     % Globales-Locales

for e   = 1:ba
   De   = D(gdl(e,:));
   p{e} =  k_acu{e}*T_acu{e}*De;
   P{e} =  K_acu{e}*De;
end

%% Esfuerzos al interior
sigma = zeros(ba,1);

for e = 1:ba
    sigma(e) = p{e}(4)/A;
end

%% Graficar
XYdef = zeros(size(XY));
fac   = 1;
c     = 0;
for i = 1:j
    c          = c+1;
    XYdef(i,1) = XY(i,1)+fac*D(c);
    c          = c+1;
    XYdef(i,2) = XY(i,2)+fac*D(c);
    c          = c+1;
end

% Se grafican la estructura y su deformada
figure
for e=1:ba
    Q    = [XY(IJ(e,1),1)     XY(IJ(e,1),2);...
            XY(IJ(e,2),1)     XY(IJ(e,2),2)];
        
    Qdef = [XYdef(IJ(e,1),1)  XYdef(IJ(e,1),2);...
            XYdef(IJ(e,2),1)  XYdef(IJ(e,2),2)];
        
    if sigma(e)>0
        plot(Q(:,1),Q(:,2),'k--',Qdef(:,1),Qdef(:,2),'b')
    else
        plot(Q(:,1),Q(:,2),'k--',Qdef(:,1),Qdef(:,2),'r')
    end
    hold on
end
title('Estructura y Deformada: Uniones Rigidas')
xlabel('x [m]')
ylabel('y [m]  factor = 1')
axis equal

% Se grafican las fuerzas axiales
figure
for e=1:ba
    subplot(3,7,e)
    plot([0 le(e)],[0 0],'lineWidth',5)
    title('Barra',e)
    hold on
    plot([0 0 le(e) le(e)],[0 p{e}(4) p{e}(4) 0])
    text(le(e)/5,   p{e}(4)/2,   num2str(p{e}(4)),'FontSize',18);
end
sgtitle('Fuerza Axial: Uniones Rigidas')

% Se grafican las fuerzas cortantes
figure
for e=1:ba
    subplot(3,7,e)
    plot([0 le(e)],[0 0],'lineWidth',5)
    title('Barra',e)
    hold on
    plot([0 0 le(e) le(e)],[0 p{e}(5) p{e}(5) 0])
    text(le(e)/5,   p{e}(5)/2,   num2str(p{e}(5)),'FontSize',18);
end
sgtitle('Fuerza Cortante: Uniones Rigidas')

% Se grafican los momentos
figure
for e=1:ba
    subplot(3,7,e)
    plot([0 le(e)],[0 0],'lineWidth',5)
    title('Barra',e)
    hold on
    plot([0 0 le(e) le(e)],[0 p{e}(3) -p{e}(6) 0])
    text(0,   p{e}(3)/2,   num2str(p{e}(3)),'FontSize',11);
    text(3*le(e)/5,   -p{e}(6)/2,   num2str(p{e}(6)),'FontSize',11);
end
sgtitle('Momento: Uniones Rigidas')

%% Extraigo tablas resumen
[Axial, Cortante, Momento] = deal(zeros(ba,1));

for e = 1:ba
    Axial(e)    = p{e}(4);
    Cortante(e) = p{e}(5);
    Momento(e)  = p{e}(6);
end

disp('Resumen Uniones Rigidas--------------------------------------------')

T2_1  = table(bar', sigma, Axial, Cortante, Momento,...
     'VariableNames',{'Barra','Esfuerzo','F Axial','F Cortante','Momento'})

[Dx, Dy]  = deal(zeros(j,1));
c         = 0;

for e     = 1:3:3*j
    c     = c+1;
    Dx(c) = D(e);
    Dy(c) = D(e+1);
end

T2_2 = table(n',Dx, Dy,'VariableNames',{'Nodo','Desp. X','Desp. Y'})

%% ------------------------------------------------------------------------
% Comparacion
% -------------------------------------------------------------------------
%% Reacciones

disp('Comparacion--------------------------------------------------------')

Reac = table(P1, Pa, 'VariableNames',{'Uniones Articuladas', 'Uniones Rigidas'})

%% Desplazamientos

Dx1 = table2array(T1_2(:,2));
Dx2 = table2array(T2_2(:,2));

Dxx = table(Dx1,Dx2,'VariableNames',{'Uniones Articuladas', 'Uniones Rigidas'})

Error = max(abs(Dx1-Dx2));
fprintf('\n La maxima diferencia en los desplazamientos \n en x fue de: %c m \n',Error)

Dy1 = table2array(T1_2(:,3));
Dy2 = table2array(T2_2(:,3));

Dyy = table(Dy1,Dy2,'VariableNames',{'Uniones Articuladas', 'Uniones Rigidas'})

Error = max(abs(Dy1-Dy2));
fprintf('\n La maxima diferencia en los desplazamientos \n en y fue de: %c m \n',Error)

%% Axial

Ax1 = table2array(T1_1(:,3));
Ax2 = table2array(T2_1(:,3));

Axx = table(Ax1,Ax2,'VariableNames',{'Uniones Articuladas', 'Uniones Rigidas'})

Error = max(abs(Ax1-Ax2));
fprintf('\n La maxima diferencia en la fuerza axial \n fue de: %c kN \n',Error)
