% -------------------------------------------------------------------------
% 12/20/2011 Blade optimization code for strength and cavitation concerns.
%
% Written in Scilab by Stefano Brizzolara et al. in code BB-OpenProp.sce
%
% Translated to Matlab on 12/20/2011 by Brenden Epps
%
% Reference:
%   Brizzolara, Tincani, and Grassi, "Design of contra-rotating propellers
%   for high-speed stern thrusters," Ships and Offshore Structures, 2:2 
%   p.169-182, 2007
% -------------------------------------------------------------------------
%
% Returns:
%       CoD     chord/diameter  distribution at RC
%       t0oc    thickness/chord distribution at RC
%       C_res   maximum residual CoD between iterations
%
% -------------------------------------------------------------------------


function [CoD, t0oc, C_res] = Chord_Brizzolara(rho,n,D,Vs,H,Mp,Z,KT,KQ, RC,G,VSTAR,TANBIC,CoD,t0oc)
                     

                     
% -------------------------------------------------------------------------
% Translate OpenProp input variables into local variable names:

T = rho * KT * n^2 * D^4; % [N]  thrust (spinta)
Q = rho * KQ * n^2 * D^5; % [NM] torque

RO   = rho; % ro H2O di mare [kg/m3] == seawater density

Vsms = Vs; % velocita nave [m/s] == ship speed

giri = n*60; % giri dell'elica [rpm] == propeller rotation speed

Hs   = H;     % immersione della linea d'asse [m] == shaft centerline depth

R    = D/2;     % raggio max [m] == propeller radius

S    = Mp;     % numero di sezioni di pala considerate, intero == number of blade sections

iangle1  = 0*RC;  % [rad] angolo di rake positivo per aft rake == rake angle at RC
skangle1 = 0*RC;  % [rad] angolo di skew, positivo se opposto al senso di rotazione == skew angle

ROmat = 8.3; % [kg/dm3] densita' del materiale dell'elica == blade material density

X     = RC; % r/R radii

l     =         CoD * D;  % [m] chord length at RC
tx    = t0oc .* CoD * D;  % [m] thickness    at RC

p     = 1 ./ (RC .* TANBIC); % by definition, where p == 2*pi*R/(pitch/D)


CLl   = 2*(2*pi*R)*G'./VSTAR;  % == CL * CoD * D == lift coefficient * chord

TGBETAI = TANBIC;  % tan(betaI) at RC


CoD0  = CoD;  % last value of CoD
t0oc0 = t0oc; % last value of t0oc

% -------------------------------------------------------------------------
                     

% -------------------------------------------------------------------------
% COSTANTI FISICHE == physical constants
g  = 9.81;    % accelerazione di gravita' [m/s^2] == acceleration due to gravity
Ha = 10.1;    % pressione atmosferica [m], in colonna d'acqua == atmospheric pressure head
Hv = 0.17;    % pressione di vapore [m], in colonna d'acqua   == water vapor pressure head


% INPUT PER ELICA == design inputs
P = 0.9;   % frazione di CL del profilo generata per camber elica == fraction of CL to generate using camber

% TIPO DI PROFILO == blade profile type
foil = 1;  % 1--> NACA16, a=0.8,  2--> NACA16, a=1.0,   3--> NACA66, a=0.8,  4--> NACA66, a=1.0, 5--> NACA16, NACA65,  6--> NACA66, NACA65
              

% FATTORI CONTROLLO, ROBUSTEZZA E MARGINE ALLA CAVITAZIONE                  == factors that control the margin on strength and cavitation
sigmam     = 500*10^6;  %  [Pa] carico di snervamento                       == material yield stress
Ksicurezza = 5.0;       %  fattore di sicurezza sulla sigma_amm dell'elica  == margin for fatigue stress loading
Krob       = 1.0;       %  margine di robustezza                            == margin for stress loading
Kcav1      = 0.85;      %  margine alla cavitazione elica1                  == cavitation margin
epsi       = 0.0001;    %  tolleranza sulla convergenza di tsucmin          == tolerance for convergence of tsucmin
tsucmax1   = 0.24;      %  valore MAXIMUM alla radice (non dovrebbe servire)== maximum allowable thickness/chord ratio at the root
txmin      = 0;  % [m] typically txmin = 0.0015*D  spessore minimo al tip                           == minimum thickness at tip [m]


% VALORI DI INPUT DERIVATI == derived inputs
LAMBs   = Vsms/(pi*n*D);              % coefficiente d'avanzo == absolute advance coefficient
CTs     = T/(RO/2*(D^2/4)*pi*Vsms^2); % coefficiente di spinta == thrust coefficient
ROmatI  = ROmat*2.205/1000/0.394^3;   % [lb/in^3] blade material density
sigma_r = sigmam/Ksicurezza/Krob;     % [Pa] allowable stress
H       = Ha+Hs-Hv-((X)*(D/2));       % total pressure head for blade at top dead center, at X==RC

SIGMA   = 2*g*H ./ (VSTAR*Vs).^2;     % == rho*g*H ./ (0.5*rho*(VSTAR*Vs).^2), local cavitation number
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% .............................
% .............................
% CALL OF THE FUNCTION THAT CALC CHORDS AND THICKNESS
% .............................
% .............................
% -------------------------------------------------------------------------
% Valuta il Ct^2 con il metodo di Conolly (Trans. RINA 1961)
% e risolve l'equazione di 3?Ǭ? per trovare la corda (c)
%  con il metodo dei profili sottili per la valutazione del Cpmin
%   da Castagento Maioli (Naval Hydrodynamics (1968))
% -------------------------------------------------------------------------

% h is indexed versus "foil" type:  
h1=[0.278 0.250 0.278 0.250 0.310 0.310];
h2=[1.132 1.132 1.270 1.270 1.132 1.270];
h3=[0.131 0.131 0.130 0.130 0.131 0.130];

% (xV,aV,AA1) and (xV,aV,AA2) data:  (i.e. each is size [9,13])
AA1 =[7.848,7.888,7.704,7.36,6.896,6.456,6.024,5.608,5.28,4.928,4.632,4.36,4.16;
      8.64,8.576,8.248,7.752,7.096,6.456,5.928,5.44,5.016,4.64,4.312,4.016,3.832;
      9.712,9.448,8.952,8.232,7.448,6.712,6.096,5.56,5.12,4.72,4.36,4.064,3.832;
      11.032,10.536,9.88,8.96,8.04,7.264,6.56,5.936,5.448,5.04,4.648,4.32,4.048;
      12.72,12.056,11.12,10.04,8.96,8.048,7.28,6.568,6.008,5.528,5.144,4.8,4.488;
      14.944,14.04,12.88,11.536,10.248,9.224,8.312,7.592,6.88,6.304,5.84,5.432,5.104;
      17.68,16.40,14.976,13.36,11.888,10.64,9.64,8.72,7.992,7.312,6.736,6.28,5.912;
      19.0,17.568,16.0,14.256,12.712,11.36,10.256,9.352,8.52,7.8,7.184,6.656,6.24;
 19.34,18.01,16.42,14.61,13.00,11.60,10.50,9.55,8.69,7.95,7.30,6.78,6.32];%,
%      00.00,00.00,00.00,00.00,00.00,00.00,00.00,0.00,0.00,0.00,0.00,0.00,0.00];
% ... %

AA2 =[ 70.88,69.0,66.2,62.68,59.04,55.32,51.68,48.64,45.6,42.88,40.52,38.32,36.08;
       53.80,51.2,48.12,44.24,40.32,36.64,33.28,30.60,28.4,26.12,24.32,22.48,21.0;
       46.44,43.24,39.92,36.08,32.24,28.96,25.96,23.68,21.24,20.0,18.4,17.12,15.84;
       43.28,40.0,36.4,32.48,28.88,25.68,22.80,20.68,18.84,17.48,16.2,15.0,14.0;
       43.0,39.52,35.6,31.64,28.08,24.8,22.12,20.08,18.40,16.84,15.6,14.4,13.56;
       45.0,40.88,36.8,32.48,28.8,25.6,22.88,20.8,19.2,17.52,16.12,15.0,14.0;
       48.04,43.8,39.48,35.24,31.04,27.6,24.76,22.6,20.68,19.16,17.6,16.36,15.2;
       47.88,43.8,39.6,35.24,31.2,27.8,25.0,22.8,20.8,19.28,18.0,16.72,15.44;
       43.9,40.5,36.7,32.7,28.9,25.7,22.8,20.6,18.8,17.3,16.2,15.1,14.1];
     %  00.00,00.00,00.00,00.00,00.00,00.00,00.00,0.00,0.00,0.00,0.00,0.00,0.00];

% KshapeV, CC1, and CC2 are each given versus xV
KshapeV = [0.1070 0.1076 0.0900 0.0647 0.0395 0.0195 0.0070 0.0013 0.0];    % size [1,9]
CC1     = [16.6 12.0 9.5 7.7 6.2 4.8 3.4 1.9 0];                            % size [1,9]
CC2     = [57.4 33.5 21.1 13.4 8.3 4.6 2.1 0.6 0];                          % size [1,9]

aV = [1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0];  % size [1,13]  p == 1./(RC.*TANBIC)
xV = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];                  % size [1,9]   X == RC

[aVmat,xVmat] = meshgrid(aV,xV);  % each is size [9,13]

           
sigma_c  = zeros(S,1); % centrifugal stress at each blade section
suCorda  = zeros(S,1); % suCorda == 1/(chord)^2

% -------------------------------------------------------------------------
gapMax = 1;  % residual
step   = 0;  % iteration

% fig(1000); plot(X, l/D,'k')
% fig(2000); plot(X,tx/D ,'k')

while gapMax > 0.01 & step < 1000
    % .............................
    % [SIGMA, suCorda, sigma_c, Q, ctsq, ctsq_NEW, tx, l, pax, gapMax]=C1mono();
    % .............................

    % Valuta il Ct^2 con il metodo di Conolly (Trans. RINA 1961)
    % e risolve l'equazione di 3?Ǭ? per trovare la corda (c)
    %  con il metodo dei profili sottili per la valutazione del Cpmin
    %   da Castagento Maioli (Naval Hydrodynamics (1968))           
                    
    % for each blade section
      for s=1:S,

          vacca11(s) = interp2(aVmat,xVmat,AA1, p(s),X(s), 'spline');
          vacca21(s) = interp2(aVmat,xVmat,AA2, p(s),X(s), 'spline');
          
          Kshape1(s) = interp1(xV,KshapeV,      X(s)     ,'spline');     
          
          VAL(s) =  R*Kshape1(s)/Z * ( vacca11(s)*p(s)*T + vacca21(s)*Q/R );  % uso T e non T/2 perche' non sono CR
          
          if  VAL(s) < 1e-9 
              VAL(s) = 0.0;
          end
          
          ctsq(s) = VAL(s)/ (sigma_r-sigma_c(s));

          if ctsq(s) < 0               % to enforce minimum thickness, use:   if ctsq(s) < txmin^2*l(s)
             ctsq(s) = txmin^2*l(s);
          end

          %calcolo dei coefficienti dell'equazione implicita descritta in CETENA n�?45 == calculate the coefficients of equaiton (31)
          cet11 = ( h2(foil)*sqrt(ctsq(s)) );
          cet21 = ( h1(foil)*P + h3(foil)*(1-P) ) * CLl(s);
          cet31 = ( 1 - sqrt(1+Kcav1*SIGMA(s)) );
          
          if any(isnan([cet11,cet21,0,cet31]))
              CoD   = 0*RC;
              t0oc  = 0*RC;
              C_res = 0;
              return
          end
          
          % Solve cubic equation (31): cet11*x^3 + cet21*x^2 + 0*x + cet31 = 0
          xroot1 = roots([cet11,cet21,0,cet31]);
          
          
          for i=1:length(xroot1),  % thee roots of a cubic equation
              if  isreal(xroot1(i))
                  suCorda(s) = max( [suCorda(s), xroot1(i)] );
              end
          end
          
          l(s) = suCorda(s)^(-2);   % l == chord [m]
      end

          tx   = sqrt( ctsq ./ l ); % tx == thickness [m]


    % Calculation of centrifugal stress sigma_c [Pa]
    % calcolo di sigma_c sollecitazione dovuta alla forza centrifuga, espressi in Pa
    % coeff. C1 e C2 vanno calcolati per elica 1 e 2
    % 1 in = 2.54 cm = 0.0254m => 1 m = 39.37 in
    % 1 tons/in^2 = 15.44 MPa da PNA vol.I
      for s=1:S,

          vacc_1CC1(s)= interp1(xV,CC1, X(s), 'spline');
          vacc_1CC2(s)= interp1(xV,CC2, X(s), 'spline');

          sigma_c(s)  = (((giri^2*(D/2*39.37)^2)/10^10)*(vacc_1CC1(s)+(tan(iangle1(s))*(D/2*39.37)*vacc_1CC2(s)*cos(atan(TGBETAI(s)))/(tx(s)*39.37))))/0.3 * ROmatI*15.44*10^6;

          ctsq_NEW(s) = VAL(s)/ (sigma_r-sigma_c(s));

          if ctsq_NEW(s) < 0           % to enforce minimum thickness, use:   if ctsq(s) < txmin^2*l(s)
             ctsq_NEW(s) = txmin^2 * l(s);
          end

          gap1(s) = abs( (ctsq_NEW(s)-ctsq(s))/ctsq(s) );

      end
          gapMax = max(gap1);
      % .............................
      % .............................
      
%     fig(1000); plot(X, l/D,'b')
%     fig(2000); plot(X,tx/D,'r')
    
    step = step+1;
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
  



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% .............................
% [lmod, txmod, min1, kk]=C2mono(l, tx)
% .............................

% modifica le corde alla radice per r/R < 0.4 imponendo un valore dello spessore e mantendo costante ct^2.
% Modify the chord at radii r/R < 0.4 to impose maximum t/c at hub, but maintain constant ct^2
txh = tx(1);  % blade thickness at the hub
lh  =  l(1);  % blad  chord     at the hub

if  tx(1)/l(1) > tsucmax1
    
    txh = (tsucmax1*ctsq(1))^(1/3);  % enforce maximum t/c at hub
else
    txh = tx(1);
end


for i=1:S-1,
    
     mtxi(i) = (tx(i+1)-tx(i))/(X(i+1)-X(i));  %elica1
     mtxh(i) = (tx(i+1)-txh  )/(X(i+1)-X(1));
    
    demtx(i) = mtxi(i) - mtxh(i);
    
end

[min1, kk] = min(abs(demtx));


txmod = tx;  % modified thickness
 lmod = l;   % modified chord
 
for i=1:kk 
     txmod(i) =  tx(kk) - (X(kk)-X(i))*mtxh(kk);
      lmod(i) = ctsq(i)/(txmod(i)^2);
end

% NOTE:
%   If you want to enforce minimum thickness, (i.e. use:   if ctsq(s) < txmin^2*l(s) ... above)
%   then you now need to fair the tx distribution near the tip to blend this modification with the nominal blade...
%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%     fig(1000); plot(X, lmod/D,'g--')
%     fig(2000); plot(X,txmod/D,'g--')
    

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Translate local variables into outputs:

CoD   = lmod / D;       % chord/diameter

t0oc  = txmod ./ lmod;  % thickness/chord

C_res = max([abs(CoD - CoD0),abs(t0oc - t0oc0)]);  % residual chord/diameter difference
% -------------------------------------------------------------------------

% Prevent EppsOptimizer code from crashing
if any(isnan([CoD,t0oc,C_res])) | any(isinf([CoD,t0oc,C_res])) | any([CoD,t0oc,C_res] < 0)
   CoD   = 0*RC;
   t0oc  = 0*RC;
   C_res = 0;
end
     