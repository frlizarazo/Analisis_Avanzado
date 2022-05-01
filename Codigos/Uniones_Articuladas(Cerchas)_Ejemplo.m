%% ANALISIS ESTRUCTURAL AVANZADO ANALISIS DE CERCHA
%
% Franklin Andres Lizarazo Muñoz
%
%-------------------------------------------------------------------------

%% Geometría
E  = 2e5;         % kPa
A  = 0.001536;       % m^2a  = [0 0 0 0 0 20.7 110.7 41.4 110.7 41.4 180-41.4 69.3 180-41.4 69.3 180-20.7 20.7 0 0 180-20.7 20.7 110.7 41.4 180-41.4 69.3 180-20.7 20.7 180-20.7];

ba = 27;          % Numero de barras
j  = 15;          % Numero de nodos

% Se define la longitud de cada barra manualmente [en m]
le = [2 2 6 2 2 1.87 0.71 2 1.41 2 2 1.41 2 0.71 1.87 1.87 2 2 1.87 1.87 0.71 2 2 0.71 1.87 1.87 1.87];

% Se define la inclinación de cada barra manualmente [en °]
a  = [0 0 0 0 0 20.7 110.7 41.4 110.7 41.4 180-41.4 69.3 180-41.4 69.3 180-20.7 20.7 0 0 180-20.7 20.7 110.7 41.4 180-41.4 69.3 180-20.7 20.7 180-20.7];

% Se ubican las cargas en los nodos especificados en la dirección dada[en kN]
%     x1 y1 ...
P  = [0 -60 0 0 0 0 0 0 0 0 0 -60 0 -120 -41.63 -112.55 0 -180 0 0 0 0 0 -180 84.83 -224.51 0 -240 0 -300]';

% Se define la topología de los elementos
IJ = [1 2; 2 3; 3 4; 4 5; 5 6; 1 7; 2 7; 2 9; 3 9; 3 10; 4 11; 4 12; 5 12; 5 8; 6 8; 7 9; 9 10; 11 12; 12 8; 9 13; 13 10; 10 15; 15 11; 11 14; 12 14; 13 15; 15 14];

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
for x=1:ba
    ke   = E*A/le(x)*[ 1 0 -1 0
                       0 0  0 0
                      -1 0  1 0
                       0 0  0 0]; k_acu{x} = ke;
                 
    %Se extrae el angulo del elemento e del vector de angulos
    beta = a(x);
    eta  = cosd(beta); eta_acu{x} = eta;      
    mu   = sind(beta); mu_acu{x}  = mu;
    
    %Se define la matriz de transformación
    Te   = [ eta  mu   0   0
             -mu eta   0   0
               0   0 eta  mu
               0   0 -mu eta]; T_acu{x} = Te;
           
    Ke   = Te'*ke*Te; K_acu{x} = Ke;
    ge   = gdl(x,:);
    
    K(ge,ge) = K(ge,ge)+Ke;   
end

%% Grados de Libertad Restringidos y libres
a = [1 2 12];   b = setdiff(1:2*j,a);

Kaa = K(a,a);     Kab = K(a,b);
Kba = K(b,a);     Kbb = K(b,b);

Pb  = P(b);

% Se resuelve el sistema de ecuaciones
Db   = Kbb\Pb;
D    = zeros(j*2,1);
D(b) = Db;

Pa   = Kab*Db;

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

% Se define la posición de cada uno de los nodos
XY = [ 0 0; 2 0; 4 0; 10 0; 12 0; 14 0; 1.75 0.66; 14-1.75 0.66; 3.5 1.32; 5.5 1.32; 14-5.5 1.32; 14-3.5 1.32; 3*1.75 1.98; 14-3*1.75 1.98; 7 2.65];
 
XYdef = zeros(size(XY));
fac   = 0.01;
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
xlabel('x [m]')
ylabel('y [m]  factor = 1')
axis equal

% Se grafican las fuerzas axiales
figure
for e = 1:ba
    subplot(3,9,e)
    plot([0 le(e)],[0 0],'lineWidth',5)
    title('Barra',e)
    hold on
    plot([0 0 le(e) le(e)],[0 p{e}(3) p{e}(3) 0])
    text(le(e)/5,   p{e}(3)/2,   num2str(p{e}(3)),'FontSize',18);
end

%% Extraigo tablas resumen
bar = 1:ba;
a   = 1:j;

[Axial, Cortante] = deal(zeros(ba,1));

for e = 1:ba
    Axial(e)    = p{e}(3);
    Cortante(e) = p{e}(2);
end

T1  = table(bar', sigma', Axial, Cortante, zeros(27,1),...
     'VariableNames',{'Barra','Esfuerzo','F Axial','F Cortante','Momento'});

[Dx, Dy]  = deal(zeros(j,1));
c         = 0;

for e     = 1:2:2*j
    c     = c+1;
    Dx(c) = D(e);
    Dy(c) = D(e+1);
end

T2 = table(a',Dx, Dy,'VariableNames',{'Nodo','Desp. X','Desp. Y'});
