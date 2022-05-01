%% TALLER 1 ANALISIS ESTRUCTURAL AVANZADO
%
% Franklin Andres Lizarazo Muñoz
%
%-------------------------------------------------------------------------

%% Geometría
E  = 2e5;         % kPa
A  = 0.005;       % m^2
ba = 21;          % Numero de barras
j  = 12;          % Numero de nodos

% Se define la longitud de cada barra manualmente [en m]
le = [5 5 5 5 5 5 6.4 5.5 5.04 5.04 5.5 6.4 4 6.4 6.3 8.04 7 8.04 6.3 6.4 4];

% Se define la inclinación de cada barra manualmente [en °]
a  = [0 0 0 0 0 0 39 25 8 -8 -25 -39 90 -39 90 -52 90 52 90 39 -90];

% Se ubican las cargas en los nodos especificados en la dirección dada[en kN]
%     x1 y1 ...
P  = [0  0  0 -2 0 -2 0 -2 0 -2 0 -2 0 0 0 -5 0 -10 0 -15 0 -10 0 -5]';

% Se define la topología de los elementos
IJ = [ 1  2;  2  3;  3  4; 4  5; 5  6; 6  7; 1 12; 12 11; 11 10; 10  9; 
       9  8; 8  7; 2 12; 12 3; 3 11; 11  4; 4 10; 4  9; 5  9; 5  8; 8  6];

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
a = [1 2 13 14];   b = setdiff(1:2*j,a);

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
XY = [ 0   0;  5   0; 10   0; 15   0; 20   0; 25   0;
      30   0; 25   4; 20 6.3; 15   7; 10 6.3;  5   4];
 
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

%% Extraigo tablas resumen
bar = 1:ba;
a   = 1:j;

[Axial, Cortante] = deal(zeros(ba,1));

for e = 1:ba
    Axial(e)    = p{e}(3);
    Cortante(e) = p{e}(2);
end

T1  = table(bar', sigma', Axial, Cortante, zeros(21,1),...
     'VariableNames',{'Barra','Esfuerzo','F Axial','F Cortante','Momento'});

[Dx, Dy]  = deal(zeros(j,1));
c         = 0;

for e     = 1:2:2*j
    c     = c+1;
    Dx(c) = D(e);
    Dy(c) = D(e+1);
end

T2 = table(a',Dx, Dy,'VariableNames',{'Nodo','Desp. X','Desp. Y'});
