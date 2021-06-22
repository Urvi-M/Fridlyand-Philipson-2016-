%-------------------------------------------------------------------------%
% Author          :   Athil Althaf A Z, Urvi Mishra and Seeraja Konnur    %
% Date            :   15-June-2018                                        %
% Version         :   2.0                                                 %
% Description     :   Pancreatic beta cell Second Messenger interactions  %
% Copyright       :   IISER Kolkata,NMIMS Mumbai,IBAB Bangalore           %
%-------------------------------------------------------------------------%
Ts = 6000;              % Total simulation time (s)
dt = .01;                 % Timestep (s)
t = 0:dt:Ts;            % Defining the array for time
fNc = 1;          % Conversion factor for converting microm^-2 to microM (microM micro m^2)   Doubt
fNc1 = 21.1;
a=zeros(1,length(t));
t1 = 3000;
t2 = 5000;
t3 = 5000;
t4 = 3000;
tpde=4000; % same as above
t6 = 5000;              % AM3 input time 
t7 = 1000;              % AR7 input time
t9 = 3000;              % kATP Blocker time
t10 =5000;              % VmPLP input time
%-------------------------------------------------------------------------%
%                   Glucose Metabolism                                    %
%-------------------------------------------------------------------------%
Cac1=0.35;
ATD_in = 1.5;           % Initial value of ATD eq.2
tATD = 100;             % Time constant (s) eq.2
ATDm = 32;              % Saturated ADT value eq.4
KGF = 5200;             % Coefficient of half saturated glucose (microM) eq.4
kGFr = 7.4;             % Coefficient of coversion FFA to glucose eq.4
hgla = 5;               % Hill coefficient eq.4
FFA = 0;                % Free Fatty acid concentration (microM) eq.4
AT = 4000;              % Total free ATP and free ADP concentration (microM) eq.5
Glui = 3000;            % Initial low glucose concentration (microM) eq.7
Glu1 = 8000;            % Increased Glucose concentration (microM) eq.7
Glu = Glui*ones(1,length(t)); % Initial Glucose level array
ATD = a;   % Defining the Array for ATP -ADP ratio eq.2
ATD(1) = ATD_in;            % Initializing the ATD array
ADPf = a;  % Initializing the free ADP array
ATP = a;  % Initializing the free ATP array
ADPf(1) = AT/(1+ATD(1));   % Initializing first element of ADPf
ATP(1) = AT - ADPf(1);     % Initializing first element of ATP
%-------------------------------------------------------------------------%
%                      Membrane Potential,channels and Pumps              %
%-------------------------------------------------------------------------%
Cm = 6158;              % Membrane capacitance (fF) eq.8
Vp_in = -62;            % Initial membrane potential (mV) eq.8
gmKATPi = 24000;        % Maximum conductance for IKATP (pS) eq.9
gmKATP1 = 9000;       % Action of KATP blocker at t9 (pS)  eq.14
gmKATP = gmKATPi*ones(1,length(t));       % Normal level of KATP conductance   eq.14
EK = -75;               % Reversal potential of pottassium (mV) eq.9
KKPI2 = 1125*fNc1;       % Activation constant (microM  converted by fNc) eq.10     Doubt
kdd = 17;               % Dissosiation coefficient (microM) eq.11
ktd = 26;               % Dissosiation coefficient (microM) eq.11
ktt = 50;               % Dissosiation coefficient (microM) eq.11
kKPKa = 1;              % Inhibition constant (microM^-1) eq.13                  doubt
gmKr = 45000;           % Maximum conductance for IKr (pS) eq.15
VdKr = -9;              % Half activaton potential (mV) eq.16
kdKr = 5;               % Slope of half activation potential (mV) eq.16
fKri = 1;               % Inactivating gating variable eq.17
gmCa = 900;             % Maximum conductance for ICa (pS) eq.18
ECa = 100;              % Reversal potential for calcium (mV) eq.18
VdCa = -19;             % Half activation potential (mV) eq.19
kdCa = 9.5;             % Slope of half activation potential (mV) eq.19
kdCap = .01;            % Coefficient of constitutive channel activity eq.19
VfCa = -9;              % Half inactivation potential (mV) eq.20
kfCa = 8;               % Slope of half activation potential (mV) eq.20
gmSOC = 10;             % Maximum whole cell conductance (pS) eq.21
KNS = 200;              % Calcium ER inhibition constant (microM) eq.22
PmCa = 6000;            % Maximum IPCa current (fA) eq.23
KpCa = .2;              % Cyosolic Calcium activation constant eq.23
ENa = 70;               % Reversal potential of sodium (mV) eq.24
gmNb = 10;              % Permanent conductance for INab (pS) eq.25
gmNM3 = 0;              % Maximum conductance of NALCN channels (pS) eq.25
Vp = a;      % Defining the Array for membrane potential eq.8
Vp(1) = Vp_in;          % Adding initial value of the Membrane potential


MgADPf = a; % Initializing array for MgADP eq.12
OKATP = a;  % Initializing array for OKATP eq.11
IKATP = a;  % Initializing array for IKATP eq.9
fKPI = a;   % Initialising array for fKPI eq.9
dKri = a;   % Initializing array for dKri eq.16
IKr = a;    % Initializing array for IKr eq.15
dCai = a;   % Initializing array for dCai eq.19
fCai = a;   % Initializing array for fCai eq.19
ICa = a;    % Initializing array for ICa eq.19
fSOC = a;   % Initializing array for fSOC eq.22
ISOC = a;   % Initializing array for ISOC eq.21
IPCa = a;   % Initializing array for IPCa eq.23
INab = a;   % Initializing array for INab eq.24
gNab =a;    % Initializing array for gNab eq.24

%-------------------------------------------------------------------------%
%                      Calcium Dynamics                                   %
%-------------------------------------------------------------------------%
Cac_in = .09;           % Initial cytosolic calcium concentration (microM) eq.86
ksg = .00001;           % Coefficient of calcium sequestration rate (s^-1) eq.86
CaER_in = 160;          % Initial ER calcium concentration (microM) eq.87
PCaER = .026;              % Maximum SERCA pump rate (microMs^-1) eq.88
Kser = .5;              % Calcium activation constant (microM) eq.88
kmIP = 7;               % Maximum permeability of IP3 activated channel (microMs^-1) eq.89
kleak = .002;           % Coefficient of passive leak from ER (s^-1) eq.89
KRPCa = .35;            % Cytosolic calcium activation constant (microM) eq.90
dinact = .4;            % Cytosolic inactivation constant (microM) eq.90 %% obtained from ref.203
KIP3 = 3.2;             % IP3 activation constant (microM) eq.90
KIP3R = .3;             % PKA activation coefficient 90
Cac = a;     % Defining array for cytosolic calcium eq.86
fcf = .01;              % Fraction of free calcium in cytoplasm eq.87
fER = .03;              % Fraction of free calcium at ER eq.87
F = 96.487;             % Faraday constant CmicroM^-1 eq.86 
Cac(1) = Cac_in;        % Adding initial value for cytosolic calcium
CaER  = a;   % Defining array for ER calcium eq.87
VER = 280;                    % Effective volume of ER microm^3 eq.87
Vc = 764;                     % Effective volume of cytoplasm eq.86
CaER(1) = CaER_in;      % Adding initial value for ER calcium
Jser = a;    % Initializing Jser array eq.88
Jrel = a;    % Initializing Jrel array eq.89
fIPKA = a;   % Initializing fIPKA array eq.90
PRIP3 = a;   % Initializing Jrel array eq.90



%-----------------------------------------------49--------------------------------------------------------%

AM3i = 0.0022; %Ln (microM)
RG6_in = 0.01*fNc; %RnGN(#per micro m^2)
ReM3d_in = 0.1*fNc; %Rnd(#per micro m^2)
LRG6_in = 0.01*fNc; %LnRnGn(#per micro m^2)
GaT6_in = 0.01*fNc; %Gan(GTP)(#per micro m^2)
GaD6_in = 0.01*fNc; %Gan(GDT)(#per micro m^2)
PLCM3_in = 0.1*fNc; %Gan(GTP)En(#per micro m^2)
ReM3t = 2*fNc; %Rnt(#per micro m^2)
G6t = 20*fNc; %Gnt(#per micro m^2)
KAM3 = 0.22; %K(M)(microM)
R61 = 10; %k1n(microM)
R62 = 0.68/fNc; %k2n(per# microm^2 per s)
R62r = 6.8; %k2nr(per s)
R63 = 0.65; %k3n(per s)
R64 = 1/fNc; %k4n(per# microm^2 per s)
R65 = 0.026; %k5n(per s)
R66 = 0.4; %k6n(per s)
R67 = 1/fNc; %k7n(per# microm^2 per s)
R68 = 0.007; %k8n(per s)
R69 = 0.00105; %k9n(per s)
AM31 = 2.2;

%----------------------------------------------------------50-----------------------------------------------%

AR7i = 0.0506; %Ln (microM)
ReR40d_in = 0.1*fNc; %Rnd(#per micro m^2)
LR7_in = 0.01*fNc; %LnRn(#per micro m^2)
LRG7_in = 0.01*fNc; %LnRnGn(#per micro m^2)
GaT7_in = 0.01*fNc; %Gan(GTP)(#per micro m^2)
GaD7_in = 0.01*fNc; %Gan(GDT)(#per micro m^2)
PLCR40_in = 0.1*fNc; %Gan(GTP)En(#per micro m^2)
ReR40t = 2*fNc; %Rnt(#per micro m^2)
G7t = 20*fNc; %Gnt(#per micro m^2)
PLCpt = 3*fNc; %Ent(#per micro m^2)
KAR7 = 4.6; %K(M)(microM)
R71 = 7; %k1n(microM)
R72 = 1/fNc; %k2n(per# microm^2 per s)
R72r = 0.68; %k2nr(per s)
R73 = 0.7; %k3n(per s)
R74 = 1/fNc; %k4n(per# microm^2 per s)
R75 = 0.026; %k5n(per s)
R76 = 0.4; %k6n(per s)
R77 = 1/fNc; %k7n(per# microm^2 per s)
R78 = 5*10^-5; %k8n(per s)
R79 = 2.83*10^-4; %k9n(per s)
PLCR40 = a;
ReR40 = a;
GaD7 = a;
GaT7 = a;
LRG7 = a;
G7 = a;
Gbg7 = a;
AR7 = AR7i*ones(1,length(t));
ReR40d = a;
LR7 = a;
PLCR40(1) = PLCR40_in;
GaD7(1) = GaD7_in ;
GaT7(1) = GaT7_in ;
LRG7(1) = LRG7_in ;
%AR7(1) = AR7i;
ReR40d(1) = ReR40d_in ;
LR7(1) = LR7_in ;
AR71 = 82;

%--------------------------------------------------73-86-----------------------------------------------------%

P4P_in = 4000*fNc; % #per micro m^2
kPI = 0.0015; %per s
kPIr = 0.006; %per s
kP4P = 0.02; %per s
kP4Pr = 0.014; %per s
PI = 140000*fNc; % #per micro m^2
KCaPI = 0.3;
KP4PK = 0.5;
kcpi = 0.2;
PIP2_in = 4200*fNc; %#per micro m^2
KPIP2 = 2370*fNc; % #per micro m^2 converted to micro M
kpPL = 15;% micro mol per s
PKCa_in = 0.1;
KCaPL = 0.4; % micro M
KCCaPL = 0.2; % micro M
IP3_in = 1; % micro M
kdIP3  = 0.04; % per s
DAG_in = 23*fNc; % #per micro m^2
kdDAG = 0.05; % per s
kPKC = 3*10^-6; % per s
kPKCr = 0.0034;% per s
PKCt = 1;

%-----------------------------------------------------M3R6-------------------------------------%

GaD6 = a;
PLCM3 = a;
GaT6 = a;
ReM3d = a;
LRG6 = a;
RG6 = a;
Gbg6 = a;
G6 = a;
AM3 = AM3i*ones(1,length(t));
GaD6(1) = GaD6_in;
PLCM3(1) = PLCM3_in;
GaT6(1) = GaT6_in;
ReM3d(1) = ReM3d_in;
LRG6(1) = LRG6_in;
RG6(1) = RG6_in;
%AM3(1) = AM3i;

P4P = a;
PIP2 = a;
IP3 = a;
DAG = a;
PKCa = a;
PLCpa = a;
fPIP2 = a;
PKCi = a;
fPI = a;
PLCpf = a;

VPLP = a;
VPLC = a;
VPL = a;
P4P(1)=P4P_in;
PIP2(1)=PIP2_in;
Cac_in = 0.09;% micro M
PKCa(1)= PKCa_in;
VmPLPi = 700; % micro mol per s
VmPLP1 = 100000;  % micro mol per s
VmPLP =VmPLPi*ones(1,length(t));
VmPLC = 50; % micro mol per s
IP3(1)=IP3_in;
DAG(1) = DAG_in; 
%----------------------------------incretin-------------------------------%
% for GLP1
GLP1i=(3.2)*(10^-7);
GLP11=(6.2)*(10^-4);
ReGLd_in=0.1*fNc;
LR1_in=0.01*fNc;
LRG1_in=0.01*fNc;
GaT1_in=0.01*fNc;
GaD1_in=0.01*fNc;
ACGLP_in=0.3*fNc;
ReGLt=2*fNc;
G1t=20*fNc;
ACpt=3*fNc;
KGLP1=(3.1)*(10^-5);
R11=7;
R12=1/fNc;
R12r=0.68;
R13=0.7;
R14=1/fNc;
R15=0.026;
R16=0.4;
R17=1/fNc;
R18=0.00005;
R19=0.000283;

LR1=a;
LR1(1)=LR1_in;

LRG1=a;
LRG1(1)=LRG1_in;

ReGLd=a;
ReGLd(1)=ReGLd_in;

GaT1=a;
GaT1(1)=GaT1_in;

ACGLP=a;

ACGLP(1)=ACGLP_in;

GaD1=a;
GaD1(1)=GaD1_in;

GLP1=GLP1i*ones(1,length(t));
ACpf=a;


ReGL=a;
G1=a;
Gbg1=a;

%.....for GIP...%

GIPi=(1.37)*(10^-6);
GIP1=(3.42)*(10^-3);
ReGld_in=0.1*fNc;
LR2_in=0.01*fNc;
LRG2_in=0.01*fNc;
GaT2_in=0.01*fNc;
GaD2_in=0.01*fNc;
ACGIP_in=0.3*fNc;
ReGlt=2.5*fNc;
G2t=25*fNc;
KGIP=(1.71)*(10^-4);
R21=7;
R22=1/fNc;
R22r=0.68;
R23=0.7;
R24=1/fNc;
R25=0.026;
R26=0.4;
R27=1/fNc;
R28=0.00005;
R29=0.000283;
%...for catecholamines....%
AR3i=0.002;
AR31=3.0;
ReARd_in=0.1*fNc;
LR3_in=0.01*fNc;
LRG3_in=0.01*fNc;
GaT3_in=0.01*fNc;
GaD3_in=0.01*fNc;
ACAR_in=0.3*fNc;
ReARt=3*fNc;
G3t=20*fNc;
KAR3=0.3;
R31=7;
R32=1/fNc;
R32r=0.68;
R33=0.7;
R34=1/fNc;
R35=0.026;
R36=0.4;
R37=1/fNc;
R38=0.00005;
R39=0.000283;

LR2=a;
LR2(1)=LR2_in;

LRG2=a;
LRG2(1)=LRG2_in;

ReGld=a;
ReGld(1)=ReGld_in;

GaT2=a;
GaT2(1)=GaT2_in;

ACGIP=a;

ACGIP(1)=ACGIP_in;

GaD2=a;
GaD2(1)=GaD2_in;

GIP=GIPi*ones(1,length(t));



ReGl=a;
G2=a;
Gbg2=a;

%...for catev...%

LR3=a;
LR3(1)=LR3_in;

LRG3=a;
LRG3(1)=LRG3_in;

ReARd=a;
ReARd(1)=ReARd_in;

GaT3=a;
GaT3(1)=GaT3_in;

ACAR=a;

ACAR(1)=ACAR_in;

GaD3=a;
GaD3(1)=GaD3_in;

AR3=AR3i*ones(1,length(t));



ReAR=a;
G3=a;
Gbg3=a;
%...calmodulin...%

CaCaM_in=0.42; %uM;
% 10^3 factor is not considered here!!!
k1f=(2.3); % uM-1 s-1
k1b=(2.4); % ms-1; 
k2f=(2.3); % (units:uM-1 s-1) (forward rate constant eq 52)
k2b=(2.4); %(units: ms-1)  (backward rate constant eq 52)
k3f=(160); %(units:uM-1 s-1)  (forward rate constant eq 53)
k3b=(405); %(units:ms-1)    (backward rate constant eq 53)
k4f=(160) ; %(units:uM-1 s-1)  (forward rate constant eq 54)
k4b=(405); %(units:ms-1)    (backward rate constant eq 54)
CaM0=11.25; % (units:uM)          (total amount of calmodulin eq 55)

%...for cAMP pka Epac.....%
kda=0.01; %(units:  (kda is the coefficient of PDE-independent cAMP degradation)
Vmcam=2; % (units: umol/s) (maximum ACp acitivity)
Kpcam=0.348; %  (Units:uM)(CaM(active) activation constant)
KNCa=75; % (Units:uM)([Ca2+]c inhibition constant.)
VmACc=0.2 ;%  (Units:umol/s)(maximum ACc activity)
KmAACS=1030; %  (Units:uM) ATP activation constant
KmCACS=0.5 ;%  (Units:uM) [Ca2+]c activation constant
kACS=0.01 ; %  (Units:umol/s)coefficient of the constitutive ACc activity.
kipdei=1; % (activation or inhibition coefficient)
kipde1=-0.4; % don't know value dummy input
Kdpe=0.348; %  (Units:uM)CaMa activation constant
Kpde=3; %  (Units:uM)cAMP activation constant
PKAa_in=0.24; % conc of PKA;
kak=1; % scaling factor,
tpka=900; %  (Units:s)time constant;
Kpcm=2.9; %  (Units:uM)cAMP activation coefficient
hpca=1.4; % Hill coefficient
EPa_in=0.01;
tep=900; %  (Units:s)time constant
Kmep=20.2; %  (Units:uM)cAMP activation constant for Epac,
hce=2; %  Hill coefficient


CaCaM=a;
CaCaM(1)=CaCaM_in;
Ca2CaM=a;
Ca3CaM=a;
Ca4CaM=a;
CaM=a;
CaMa=a;

%...cAMP..%
ACpa=a;
ACpar=a;
facca=a;
VACp=a;
VACc=a;
Vpde=a;
VAC=a;
cAMP=a;
PKAi=a;
PKAa=a;
PKAa(1)=PKAa_in;
fKPKa = a;
EPi=a;
EPa=a;
EPa(1)=EPa_in;
fpca=a;
fEcA=a;
kipde=kipdei*ones(1,length(t));
Vgpde=0.04*ones(1,length(t));
Vcpde=1.4*ones(1,length(t));
%-------------------------------------------------------------------------%
%                      Modeling Part                                      %
%-------------------------------------------------------------------------%
for i = 1:length(t)-1                     
    if t(i)>t1
        % Place for inputs and dummy variables
        % Dummy inputs mimicks the behavior of those variables
        Glu(i) =  Glu1;          % Increasing the Glucose concentration eq.7 
        Cac(i) = Cac1; 
        
    end
    if t(i)>t2
      GLP1(i)=GLP11;
    end
%     if t(i)>t3
%           GIP(i)=GIP1;
%     end
%     if t(i)>t4
%         AR3(i)=AR31;
%     end
%     if t(i)> t6
%          AM3(i) = AM31;
%     end
%     if t(i)> t7;
%          AR7(i) = AR71;
%    end
%    if t(i)>t9                       % Blocking KATP eq.14
%         gmKATP(i) = gmKATP1;
%     end
%    Membrane potential part
    if t(i)> t10;
        VmPLP(i+1) = VmPLP(i+1) + VmPLPi;
    end
    if t(i)>tpde
    kipde(i+1)=kipde(i+1)+kipde1;
    end
           
    % Place for solving the differential equations and the allied parameter

    % Glucose metabolism part
    ATD0 = ATDm*(Glu(i)+kGFr*FFA)^hgla/((Glu(i)+kGFr*FFA)^hgla+KGF^hgla);     % eq.4   
    %dATD(i) = dt*(ATD0 - ATD(i))/tATD;                                  % eq.2
    %ATD(i) = ATD(i) + dATD(i);      % Appending the ATD by euler method eq.2
    ATD(i+1) = ATD(i) + dt*(ATD0 - ATD(i))/tATD;
    ADPf(i+1) = AT/(1+ATD(i));        % eq.6
    ATP(i+1) = AT - ADPf(i);          % eq.5
 
    fKPI(i) = PIP2(i)/(PIP2(i)+KKPI2);       % PIP2 is dummy eq.10  %% Membrane Mechanism
    fKPKa(i) = 1/(1+kKPKa*PKAa(i));       % PKAa is dummy eq.13
    MgADPf(i) = 0.055*fKPKa(i)*ADPf(i);                 % eq.12 
    OKATP(i) = (.088*(1+2*MgADPf(i)/kdd)+.89*(MgADPf(i)/kdd)^2)/((1+MgADPf(i)/kdd)^2*(1+.45*MgADPf(i)/ktd +ATP(i)/ktt));    % eq.11
    IKATP(i) = gmKATP(i)*fKPI(i)*OKATP(i)*(Vp(i)-EK);      % eq.9
    dKri(i) = 1/(1+exp((VdKr - Vp(i))/kdKr));        % eq.16
    IKr(i)  = gmKr*dKri(i)^2*fKri*(Vp(i)-EK);        % eq.15
    dCai(i) = 1/(1+exp((VdCa-Vp(i))/kdCa))+kdCap;    % eq.19
    fCai(i) = 1/(1+exp(-(VfCa-Vp(i))/kfCa));         % eq.20
    ICa(i) = gmCa*dCai(i)*fCai(i)*(Vp(i)-ECa);       % eq.18
    fSOC(i) = KNS^4/(KNS^4+CaER(i)^4);               % eq.22
    ISOC(i) = gmSOC*fSOC(i)*(Vp(i)-ECa);             % eq.21
    IPCa(i) = PmCa*Cac(i)^2/(KpCa^2+Cac(i)^2);       % eq.23
    gNab(i) = gmNb + gmNM3*LRG6(i);                     % eq.24
    INab(i) = gNab(i)*(Vp(i)-ENa);                   % eq.24  
%     dVp(i) = -dt*(IKATP(i)+IKr(i)+ICa(i)+ISOC(i)+IPCa(i)+INab(i))/Cm;   % eq.8
    
    Vp(i+1) = Vp(i)-dt*(IKATP(i)+IKr(i)+ICa(i)+ISOC(i)+IPCa(i)+INab(i))/Cm;  % eq.8                       % Appending the Vp array
    % Part for Calcium dynamics
    Jser(i) = PCaER*Cac(i)^2/(Kser^2 + Cac(i)^2);    % eq.88
    fIPKA(i) = 1/(PKAa(i) + KIP3R);                      % eq.90
    PRIP3(i) = (Cac(i)/(KRPCa+Cac(i)))^3*(IP3(i)/(fIPKA(i)*KIP3+IP3(i)))^3*(dinact/(dinact+Cac(i)))^3; %  eq.90
    Jrel(i) = (kmIP*PRIP3(i) + kleak)*(CaER(i) - Cac(i));      % eq.89
%    Cac(i+1) = Cac(i) + dt*(fcf*((-ICa(i)-ISOC(i)-2*IPCa(i))/(2*F*Vc)-Jser(i)+Jrel(i)/Vc) - ksg*Cac(i)); % eq.86 
    CaER(i+1)= CaER(i) + dt*(fER*(Jser(i)*Vc-Jrel(i))/VER);      % eq.87
    
    ACpf(i) = ACpt- ACGLP(i)-ACGIP(i)-ACAR(i);
    ReGL(i)=ReGLt-ReGLd(i)-LR1(i)-LRG1(i);
    G1(i)=G1t-LRG1(i)-GaT1(i)-GaD1(i)-ACGLP(i);
    
    Gbg1(i)=GaT1(i)+GaD1(i)+ACGLP(i);
    
    LR1(i+1)=LR1(i) + dt*(R11*ReGL(i)*(GLP1(i)/(KGLP1+GLP1(i))) + R12r*LRG1(i)-R12*G1(i)*LR1(i)-R19*LR1(i));
    LRG1(i+1)=LRG1(i)+ dt*(R12*G1(i)*LR1(i)-R12r*LRG1(i)-R13*LRG1(i));
    ReGLd(i+1)=ReGLd(i) + dt*(R19*LR1(i)-R18*ReGLd(i));
    GaT1(i+1)=GaT1(i) + dt*(R13*LRG1(i)-R14*GaT1(i)*ACpf(i)-R15*GaT1(i));
    ACGLP(i+1)= ACGLP(i)+dt*(R14*GaT1(i)*ACpf(i)-R16*ACGLP(i));
    GaD1(i+1)=GaD1(i) + dt*(R15*GaT1(i)+R16*ACGLP(i)-R17*GaD1(i)*Gbg1(i));
    
    %......for GIP....%
    
      
    ReGl(i)=ReGlt-ReGld(i)-LR2(i)-LRG2(i);
    G2(i)=G2t-LRG2(i)-GaT2(i)-GaD2(i)-ACGIP(i);
    
    Gbg2(i)=GaT2(i)+GaD2(i)+ACGIP(i);
    
    LR2(i+1)=LR2(i) + dt*(R21*ReGl(i)*(GIP(i)/(KGIP+GIP(i))) + R22r*LRG2(i)-R22*G2(i)*LR2(i)-R29*LR2(i));
    LRG2(i+1)=LRG2(i)+ dt*(R22*G2(i)*LR2(i)-R22r*LRG2(i)-R23*LRG2(i));
    ReGld(i+1)=ReGld(i) + dt*(R29*LR2(i)-R28*ReGld(i));
    GaT2(i+1)=GaT2(i) + dt*(R23*LRG2(i)-R24*GaT2(i)*ACpf(i)-R25*GaT2(i));
    ACGIP(i+1)= ACGIP(i)+dt*(R24*GaT2(i)*ACpf(i)-R26*ACGIP(i));
    GaD2(i+1)=GaD2(i) + dt*(R25*GaT2(i)+R26*ACGIP(i)-R27*GaD2(i)*Gbg2(i));
    
    
    %.....for catecholamines...%
    
    ReAR(i)=ReARt-ReARd(i)-LR3(i)-LRG3(i);
    G3(i)=G3t-LRG3(i)-GaT3(i)-GaD3(i)-ACAR(i);
    
    Gbg3(i)=GaT3(i)+GaD3(i)+ACAR(i);
    
    LR3(i+1)=LR3(i) + dt*(R31*ReAR(i)*(AR3(i)/(KAR3+AR3(i))) + R32r*LRG3(i)-R32*G3(i)*LR3(i)-R39*LR3(i));
    LRG3(i+1)=LRG3(i)+ dt*(R32*G3(i)*LR3(i)-R32r*LRG3(i)-R33*LRG3(i));
    ReARd(i+1)=ReARd(i) + dt*(R39*LR3(i)-R38*ReARd(i));
    GaT3(i+1)=GaT3(i) + dt*(R33*LRG3(i)-R34*GaT3(i)*ACpf(i)-R35*GaT3(i));
    ACAR(i+1)= ACAR(i)+dt*(R34*GaT3(i)*ACpf(i)-R36*ACAR(i));
    GaD3(i+1)=GaD3(i) + dt*(R35*GaT3(i)+R36*ACAR(i)-R37*GaD3(i)*Gbg3(i));
    %....calmodulin...%
  
  Ca2CaM(i)=(k2f/k2b)*Cac(i)*CaCaM(i);
  Ca3CaM(i)=(k3f/k3b)*Cac(i)*Ca2CaM(i);
  Ca3CaM(i)=(k4f/k4b)*Cac(i)*Ca3CaM(i);
  
  CaM(i)=CaM0-CaCaM(i)-Ca2CaM(i)-Ca3CaM(i)-Ca4CaM(i);
  CaMa(i)=Ca3CaM(i)+Ca4CaM(i);
  
  CaCaM(i+1)=CaCaM(i) + dt*(k1f*Cac(i)*CaM(i)-k1b*CaCaM(i));

%...for cAMP...%

ACpf(i)=ACpt-ACGLP(i)-ACGIP(i)-ACAR(i);
ACpa(i)=ACGLP(i)+ACGIP(i);
ACpar(i)=ACpa(i)/ACpt;

facca(i)=((CaMa(i)/CaMa(i)+Kpcam))*(KNCa/(KNCa+Cac(i)));
VACp(i)=(Vmcam)*(facca(i))*ACpar(i);
VACc(i)=(((VmACc)*ATP(i))/(ATP(i)+KmAACS))*((Cac(i))/(Cac(i)+KmCACS))+kACS;
Vpde(i)=kipde(i)*(Vgpde(i)+Vcpde(i)*(CaMa(i)/(CaMa(i)+Kdpe)))*(cAMP(i)/(cAMP(i)+Kpde));

VAC(i)=VACp(i) + VACc(i);
cAMP(i+1)=cAMP(i) + dt*(VAC(i)-Vpde(i)-kda*cAMP(i));

    %...for PKA...%
    PKAi(i)=1-PKAa(i);
    fpca(i)=(cAMP(i)^1.4)/((Kpcm^1.4)+(cAMP(i)^1.4));
    PKAa(i+1)=PKAa(i) + (dt*(kak*(fpca(i)*PKAi(i)-PKAa(i))/tpka));

%...for Epac...%

EPi(i)=1-EPa(i);
fEcA(i)=(cAMP(i)^1.4)/((Kmep^1.4)+(cAMP(i)^1.4));
EPa(i+1)=EPa(i) + (dt*(kak*(fEcA(i)*EPi(i)-EPa(i))/tep));

    
    % IP3,PIP2 and PLC pathways
    PKCi(i) = PKCt - PKCa(i);                                                                          %% eq 85
    PKCa(i+1) =PKCa(i)+ dt*(kPKC*DAG(i)*PKCi(i) - kPKCr*PKCa(i));                                      %% eq 84
    fPI(i) = (Cac(i)^2/((Cac(i)^2)+(KCaPI^2)))*(PKCa(i)^2/((PKCa(i)^2)+(KP4PK^2)))+kcpi;               %% eq 74
    P4P(i+1) = P4P(i)+ dt*(kPI*fPI(i)*PI+kP4Pr*PIP2(i) - kPIr*P4P(i) - kP4P*P4P(i));                   %% eq 73
    fPIP2(i) = PIP2(i)/(KPIP2+PIP2(i));                                                                %% eq 76
    PLCpa(i) = PLCM3(i) + PLCR40(i);                                                                   %% eq 79
    PLCpf(i) = PLCpt - PLCpa(i);                                                                       %% eq 80
    VPLP(i) = (VmPLP(i))*(PLCpa(i)/PLCpt)*(Cac(i)^2/((KCaPL^2)+Cac(i)^2));                             %% eq 78
    VPLC(i) = VmPLC*(Cac(i)^2)/(KCCaPL^2 + Cac(i)^2);                                                  %% eq 81
    VPL(i) = VPLP(i) + VPLC(i) + kpPL;                                                                 %% eq 77
    PIP2(i+1) = PIP2(i)+ dt*(kP4P*P4P(i) - fPIP2(i)*VPL(i) - kP4Pr*PIP2(i));                           %% eq 75
    IP3(i+1) =IP3(i)+ dt*(VPL(i)*fNc*fPIP2(i) - kdIP3*IP3(i));                                         %% eq 82
    DAG(i+1) =DAG(i)+ dt*(VPL(i)*fPIP2(i) - kdDAG*DAG(i));   
%------------------------------------------------------------ M3R ----------------------------------------------------------%
    G6(i) = G6t - RG6(i) - LRG6(i) - GaT6(i) - GaD6(i) + PLCM3(i);                                  %% eq 49
    Gbg6(i) = GaT6(i) + GaD6(i) + PLCM3(i);                                                         %% '' 
    RG6(i+1) = RG6(i)+ dt*(((R62*ReM3t*G6(i)) -(R61*RG6(i)*AM3(i)))/(KAM3+AM3(i)) - R62r*RG6(i));   %% ''
    LRG6(i+1) = LRG6(i)+ dt*(R61*RG6(i)*AM3(i)/(KAM3+AM3(i)) - R63*LRG6(i) - R69*LRG6(i));          %% ''
    ReM3d(i+1) = ReM3d(i) + dt*(R69*LRG6(i) - R68*ReM3d(i));                                        %% ''
    GaT6(i+1) = GaT6(i)+ dt*(R63*LRG6(i) - R64*GaT6(i)*(PLCpf(i)) - R65*GaT6(i));                      %% ''
    PLCM3(i+1) = PLCM3(i)+ dt*(R64*GaT6(i)*PLCpf(i) - R66*PLCM3(i));                                   %% ''
    GaD6(i+1) = GaD6(i)+ dt*(R65*GaT6(i) +R66*PLCM3(i) - R67*GaD6(i)*Gbg6(i));                      %% ''
    ReM3 = ReM3t - ReM3d(i) - RG6(i) - LRG6(i); 

%----------------------------------------------------------FFAR1/GPR40-----------------------------------------------------%
                                                                                 
    ReR40(i) = ReR40t - ReR40d(i) - LRG7(i) - LR7(i);                                              %% eq 50
    G7(i) = G7t - LRG7(i) - GaT7(i) - GaD7(i) - PLCR40(i);                                         %% ''  
    Gbg7(i+1) = GaT7(i) + GaD7(i) + PLCR40(i);                                                     %% ''
    LR7(i+1) = LR7(i) + dt*(R71*ReR40(i)*AR7(i)/(KAR7 + AR7(i)) + R72r*LRG7(i) - R72*G7(i)*LR7(i) - R79*LR7(i));  %% ''
    LRG7(i+1) = LRG7(i)+ dt*(R72*G7(i)*LR7(i)*R72r*LRG7(i) - R73*LRG7(i));                         %% ''
    ReR40d(i+1) = ReR40d(i)+ dt*(R79*LR7(i) - R78*ReR40d(i));                                      %% ''
    GaT7(i+1) = GaT7(i)+ dt*(R73*LRG7(i) - R74*GaT7(i)*PLCpf(i) - R75*GaT7(i));                       %% ''
    PLCR40(i+1) = PLCR40(i)+ dt*(R74*GaT7(i)*PLCpf(i) - R76*PLCR40(i));                               %% ''
    GaD7(i+1) = GaD7(i)+ dt*(R75*GaT7(i) + R76*PLCR40(i) - R77*GaD7(i)*Gbg7(i));                   %% ''
end

%~PLOTS~%
% Section of plots
% figure(1);
% plot(t,ATD0)
% xlabel('time s')
% ylabel('Units')
% title('ATD_0')
% 
% figure(2);
% plot(t,ATD)
% xlabel('time s')
% ylabel('Units')
% title('ATD')
% 
% figure(3);
% plot(t,ADPf)
% xlabel('time s')
% ylabel('Concentration \muM')
% title('free ADP (ADPf)')
% 
% figure(4);
% plot(t,ATP)
% xlabel('time s')
% ylabel('Concentration \muM')
% title('ATP')
% 
% figure(5);
% plot(t,fKPI)
% xlabel('time s')
% ylabel('Units')
% title('f_{KPI}')
% 
% figure(6);
% plot(t,fKPKa)
% xlabel('time s')
% ylabel('Concentration^{-1} \muM^{-1}')        %Check this part 
% title('f_{KPKa}')
% 
% figure(7);
% plot(t,MgADPf)
% xlabel('time s')
% ylabel('Concentration \muM')
% title('Mg Bound ADP MgADP_f')
% 
% figure(8);
% plot(t,OKATP)
% xlabel('time s')
% ylabel('Probability')
% title('O_{KATP} fraction of open channels')
% 
% figure(9);
% plot(t,IKATP)
% xlabel('time s')
% ylabel('Current fA')
% title('ATP Dependant  pottassium channel current I_{KATP}')
% 
% figure(10);
% plot(t,dKri)
% xlabel('time s')
% ylabel('Probability')
% title('d_{Kri} channel activation probability')
% 
% figure(11);
% plot(t,IKr)
% xlabel('time s')
% ylabel('Current fA')
% title('Voltage Dependant  pottassium channel current  I_{Kr}')
% 
% figure(12);
% plot(t,dCai)
% xlabel('time s')
% ylabel('Probability')
% title('d_{Cai} channel activation probability')
% 
% figure(13);
% plot(t,fCai)
% xlabel('time s')
% ylabel('Probability')
% title('f_{Cai} channel inactivation probability')
% 
% figure(14);
% plot(t,ICa)
% xlabel('time s')
% ylabel('Current fA')
% title('Voltage Dependant  pottassium channel current I_{Ca}')
% 
% figure(15);
% plot(t,fSOC)
% xlabel('time s')
% ylabel('Probability')
% title('f_{SOC} channel inactivation probability')
% 
% figure(16);
% plot(t,ISOC)
% xlabel('time s')
% ylabel('Current fA')
% title('Store operated channel current I_{SOC}')
% 
% figure(17);
% plot(t,IPCa)
% xlabel('time s')
% ylabel('Current fA')
% title('Calcium pump current I_{PCa}')
% 
% figure(18);
% plot(t,gNab)
% xlabel('time s')
% ylabel(' Conductance pS')
% title('Sodium Conductance g_{Nab}')
% 
% figure(19);
% plot(t,INab)
% xlabel('time s')
% ylabel('Current fA')
% title('Sodium background current I_{Nab}')
% 
% figure(20);
% plot(t,Vp)
% xlabel('time s')
% ylabel('Voltage mV')
% title('Membrane potential V_p')
% 
% figure(21);
% plot(t,Jser)
% xlabel('time s')
% ylabel('flux \mumols^-1 ')
% title('Flux of Calcium from cytosol to ER J_{SER}')
% 
% figure(22);
% plot(t,fIPKA)
% xlabel('time s')
% ylabel('fraction')
% title(' f_{IPKA}')
% 
% figure(23);
% plot(t,PRIP3)
% xlabel('time s')
% ylabel('Probability')
% title('Channel open probability P_{RIP3}')
% 
% figure(24);
% plot(t,Jrel)
% xlabel('time s')
% ylabel('flux \mumols^-1 ')
% title('Flux of Calcium from ER to cytosol J_{rel}')
% 
% figure(25);
% plot(t,Cac)
% xlabel('time s')
% ylabel('Concentration \muM')
% title('Cytosolic calcium concentration')
% 
% 
% figure(26);
% plot(t,CaER)
% xlabel('time s')
% ylabel('Concentration \muM')
% title('ER calcium concentration')

%....Incretin catecholamines graph...%

% % figure(27);
% % plot(t,LR1)
% % xlabel('time s')
% % ylabel('LR1')
% % title('LR1')
% % 
% % figure(28);
% % plot(t,LR2)
% % xlabel('time s')
% % ylabel('LR2')
% % title('LR2')
% % 
% % figure(29);
% % plot(t,LR3)
% % xlabel('time s')
% % ylabel('LR3')
% % title('LR3')
% % 
% % figure(30);
% % plot(t,LRG1)
% % xlabel('time s')
% % ylabel('LRG1')
% % title('LRG1')
% % 
% % figure(31);
% % plot(t,LRG2)
% % xlabel('time s')
% % ylabel('LRG2')
% % title('LRG2')
% % 
% % figure(32);
% % plot(t,LRG3)
% % xlabel('time s')
% % ylabel('LRG3')
% % title('LRG3')
% % 
% % figure(33);
% % plot(t,G1)
% % xlabel('time s')
% % ylabel('G1')
% % title('G1')
% % 
% % figure(34);
% % plot(t,G2)
% % xlabel('time s')
% % ylabel('G2')
% % title('G2')
% % 
% % figure(35);
% % plot(t,G3)
% % xlabel('time s')
% % ylabel('G3')
% % title('G3')
% % 
% % figure(36);
% % plot(t,GaD1)
% % xlabel('time s')
% % ylabel('GaD1')
% % title('GaD1')
% % 
% % figure(37);
% % plot(t,GaD2)
% % xlabel('time s')
% % ylabel('GaD2')
% % title('GaD2')
% % 
% % figure(38);
% % plot(t,GaD3)
% % xlabel('time s')
% % ylabel('GaD3')
% % title('GaD3')
% % 
% % figure(39);
% % plot(t,ACGLP)
% % xlabel('time s')
% % ylabel('ACGLP')
% % title('ACGLP')
% % 
% % figure(40);
% % plot(t,ACGIP)
% % xlabel('time s')
% % ylabel('ACGIP')
% % title('ACGIP')
% % 
% % figure(41);
% % plot(t,ACAR)
% % xlabel('time s')
% % ylabel('ACAR')
% % title('ACAR')
% %
% %figure(42);
% %plot(t,PKAa)
% %xlabel('time s')
% %ylabel('PKAa')
% %title('PKAa')
% %
% %figure(43);
% %plot(t,CaMa)
% %xlabel('time s')
% %ylabel('CaMa')
% %title('CaMa')



















