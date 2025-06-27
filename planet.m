function outstruc=planet(data,ustar_flag,depo_flag,params)
%% PLAnET Microbial Model Description
% The PLAnET model describes both phyllospheric population dynamics and net
% fluxes of microorganisms on the basis of simple meteorological variables.
%
% SYNTAX:
%
% outstruc = planet(data,ustar_flag,depo_flag,params);
%
% PROGRAM STRUCTURE:
%
% i) INPUTS: 
%
% The model needs at least three input arguments. A data matrix
% (data), a flag indicating if wind speed or friction velocity are included
% in the data matrix (ustar_flag) and a flag indicating if depositional
% effects are to be computed or not (depo_flag). A fourth argument
% (params)can be provided by the user if there's the need to use different
% model parameters from the ones derived by its original formulation.
% In the latter case, params is a vector of 11 values specifying:
% 1) Minimum Growth Temperature (°C)
% 2) Maxmimum Growth Temperature (°C)
% 3) Optimal Growth Temperature (°C)
% 4) Lower population boundary (CFU m-2)
% 5) Upper population boundary (CFU m-2)
% 6-8) three coefficients for the Gompertz equation governing upward fluxes
% 9) A scaling coefficient multiplying the growth of microorganisms
% 10-11) Slope and offset of the function relating LAI and the background
% airborne microbial concentration
%
% The input data matrix has a specific format: it has a number of rows
% corresponding to the number of observations and four columns each one
% corresponding to a driving variable. The four columns are:
% 1) Temperature (°C)
% 2) Pressure (Pa)
% 3) Wind Speed / u* (m/s)
% 4) Leaf area index
% 5) Wind Speed (if only u* is provided in column 3, see below). 
%
% The ustar flag (ustar_flag) governs how the data matrix is treated by the
% model. If the flag == 1, it assumes that the third column in the matrix
% is already representing a friction velocity; if the flag == 0 it assumes
% that, instead, wind speed is given. The model converts wind speed in
% friction velocity on the basis of the height of the wind measurement in
% meter (z, in m) and the roughness length of the surface (z0, in m). These
% parameters are fixed at 10 and 0.15 m respectively but may be modified in
% the "MODEL CONSTANTS" section of the code by the user if needed. If ustar
% is directly given, the model needs wind speed as well to compute
% deposition. Wind speed should be added as a fifth column to the input
% matrix if deposition wants to be used. 
%
% The third input of the model is the flag governing the deposition. When
% depo_flag==1 the model calculates settling velocity and depositional
% fluxes using a diameter fixed at 3.3 microns. If the user wants to change
% it, it is possible to modify "dia" in the "MODEL CONSTANTS" section. When
% depo_flag==0 no depositional fluxes are computed and the net flux is
% never corrected. Besides gravitational settling velocity the model
% computes also the deposition due to interception/impaction (see function
% "slinnmod" at lines 325ff). General impaction/interception parameters can be
% tweaked under the function's "model parameters" (lines 336-340)
%
% ii) OUTPUTS: 
%
% To avoid having the user input a dozen of output arguments,
% the model exits with only a single argument (outstruc). The latter is a
% structure containing all the main model variables:
% 1) outstruc.population = the microbial population in the phyllosphere (CFU
% m-2)
% 2) outstruc.growth = the microbial growth rate
% 3) outstruc.kmax = the variation in population cap due to leaf senescence 
% 4) outstruc.gross_out = gross outward flux (CFU m-2 s-1)
% 5) outstruc.gross_in = gross inward flux (CFU m-2 s-1)
% 6) outstruc.conc = background airborne concentration (CFU m-3)
% 7) outstruc.net_flux = net microbial flux (CFU m-2 s-1)
% 8) outstruc.vg = computed settling velocity due to gravity (m s-1)
% 9) outstruc.vimp = computed settling velocity due to interception /impaction (m s-1)
% 10) oustruc.vdep = overall settling velocity (vg + vimp) (m s-1)
% 11) outstruc.dflag = deposition flag. If it's equal to 1 deposition is active at the given timestep.
% 12) outstruc.dieout = dieout flow (CFU m-2)
%
% iii) GENERAL GUIDELINES: 
%
% 1) The model assumes an half-hourly timestep in the input, and that's why
% it divides the output values in fluxes by d=1800 in order to output
% fluxes in CFU m-2 s-1. If the user has other needs, the value of d in
% "MODEL CONSTANTS".
% 2) The model assumes to be starting from a moment when the microbial
% population is at its minimum, so it's advisable to have a dataset starting
% in winter months when LAI is at its minimum. 
%
% PLAnET PATCH 1 (24/09/2018):
% Corrected some small overflow errors in population computation (see
% comments atl line 259)


%% MODEL CONSTANTS
% set model constants
z=10; % height of the anemometer for deriving u* from wind speed (m)
z0=0.15; % roughness length (m)
d=1800; % integration time (seconds)
dia=3.3; % diameter of the biological aerosol (micrometers)
zc=0.2; % canopy average height (m)

%% INITIALIZATION
% Check number of arguments
n=nargin;
if n==4; 
Tmin=params(1);
Tmax=params(2);
Topt=params(3);
kmin=params(4);
kmax=params(5);
slp=params(6);
slp2=params(7);
slp3=params(8);
c=params(9);
lai1=params(10);
lai2=params(11);

elseif n==3;
Tmin=12.959980999122827;
Tmax=30.159885612670195;
Topt=21.559933305896510;
kmin=5.003398245516982e+04;
kmax=4.816110269787611e+05;
slp=30;
slp2=256.26;
slp3=19;
c=0.132575241852830;
lai1=26.99;
lai2=115.9;
else
    error('Wrong numer of input arguments specified! The model needs a meterological data matrix, a friction velocity flag, a depositional flag and (optionally) a set of model parameters');
end

% read data from input
T=data(:,1);
P=data(:,2);
wspeed=data(:,3);
laivals=data(:,4);

% if wind speed is provided, convert to ustar
if ustar_flag==1;
    ust=wspeed;
    % if only ustar flag is provided a 5th column with wind speed is needed
    % if using deposition
    if depo_flag==1;
        wsp=data(:,5);
    else
    end
    
else
    ust=0.4.*wspeed./log(z/z0);
    if depo_flag==1;
        wsp=wspeed;
    else
    end
end

% Set beginning population
N0=kmin; % we set it to kmin assuming model starts when there is low LAI 
         % and T
%% SET MICROBIAL POOL
for k=1:length(data(:,1));
    if k==1;
        popul=N0; % to avoid index error when p(k-1)
        kmaxbuffer=kmax;
    else
        popul=pop(k-1); % to avoid index error when p(k-1)
    end
    % correct kmax for LAI
    kmax=kmin+((kmaxbuffer-kmin).*laivals(k));
    kmaxout(k)=kmax;

%% GROWTH
if T(k)<Tmin || T(k)>Tmax;
    growth(k)=0;
else
    growth(k)=(((Tmax-T(k))./(Tmax-Topt)).*((T(k)-Tmin)./(Topt-Tmin)).^((Topt-Tmin)./(Tmax-Topt))).*c;
end

%% GROSS UP
removal(k)=((slp.*exp(-slp2.*exp(-slp3*ust(k)))).*d).*(popul./kmax); % removal is multiplied by the time-step d here (CFU m-2 s-1 * s = CFU m-2)
if removal(k)<0; % we can't have NEGATIVE upward fluxes, 
                  %deposition is treated separatedly
    removal(k)=0;
else
end
if isnan(removal(k));
    error(['At k = ',num2str(k),' removal is nan']);
else
end

%% GROSS DOWN
if depo_flag==1; % active only if deposition module is activated
    vg(k)=settvel(dia,T(k)+273.15,P(k));
    % if ustar/wind are 0 gravity works, but turbulence not so much. check
    % it. 
    if ust(k)>0 & wsp(k)>0;
    vimp(k)=slinnmod(wsp(k),ust(k),dia,vg(k),zc,z0);
    else
        vimp(k)=0;
    end
    vdep(k)=vg(k)+vimp(k);
    conc(k)=(lai1*laivals(k)+lai2);
    Fd(k)=(vdep(k).*conc(k)).*d; % deposition flux is here multiplied by the timestep (CFU m-2 s-1 * s = CFU m-2)
    dflag(k)=1;
%% NET FLUX
% For arithmetic reasons, net flux has to be negative when microorganisms
% fly away (subtracts from population) and positive when deposited (add to
% population). This arithmetic will be converted to the standard
% "environmental" one at the end of the processing.
    Fn(k)=-removal(k)+Fd(k);
else
    Fn(k)=-removal(k);
    Fd(k)=0;
end

%% ACTUAL POPULATION CHANGE
pop(k)=popul+(popul*growth(k))+Fn(k);

if pop(k)>kmax;
       % provided that we have growth, we may artificially hinder growth to
       % reduce population under kmax: too much cells, too few food. Stop
       % growing. 
       if (kmax-popul-Fn(k))./popul>=0; % check if stunted growth > 0
           if growth(k)==0;
               error('stunting growth with growth = 0 makes no sense...');
           else
           end
        growth(k)=(kmax-popul-Fn(k))./popul;
        pop(k)=popul+(popul*growth(k))+Fn(k);
        dieout(k)=0;
       else
           % if we can't stunt growth, we set it to 0 and see if that's
           % enough for not going above kmax (no space and substrate, see
           % above).
           growth(k)=0;
           pop(k)=popul+(popul*growth(k))+Fn(k);
           dieout(k)=0;
           if pop(k)>kmax;
               % if not even setting growth to 0 solves the problem, then
               % it's a LAI driven dieout. We reduce population with a
               % lateral flow. This is a dieout of microorganisms that
               % happens because we have reduced their habitat (fallen
               % leaves)
               dieout(k)=pop(k)-kmax; % this is the lateral flow
               pop(k)=popul+(popul*growth(k))+Fn(k)-dieout(k);
           else
           end
       end
elseif pop(k)<kmin;
    % no outgoing flux can happen when population is below kmin
    rem(k)=0;Fn(k)=-rem(k)+Fd(k);
        if popul+(popul.*growth(k))+Fn(k)>=kmin;
            % if the deposition flux is enough to go above threshold, then...
            pop(k)=popul+(popul.*growth(k))+Fn(k);
            dieout(k)=0;
            % in the rare cases (<0.05% in the original dataset) that LAI
            % goes very low decreasing kmax a lot and there is no growth,
            % but a substantial deposition flux, this can overshoot kmax,
            % we adjust here
            if pop(k)> kmax;
                pop(k)=kmax;
            else
            end
        elseif (kmin-popul-Fn(k))./popul>0;
           % if we have growth, we will adjust the growth rate to reach
           % kmin: we can make the assumption that when more space and 
           % resources are available, microorganisms get a growth boost...
          growth(k)=(kmin-popul-Fn(k))./popul;
          pop(k)=popul+(popul*growth(k))+Fn(k);
          dieout(k)=0;
        else
            error('Removal zeroing and boost growth not ok. There are probably major issues in input parameters. The model is stopping');
        end
        
else
        dieout(k)=0;
end

% logic checks for errors in population
ovflw(k)=0;
undflw(k)=0;
if pop(k)>kmax;
    if abs(kmax-pop(k))<1e-4*(min(abs(kmax),abs(pop(k))));
        ovflw(k)=0;
    else
        warning(['Population overflow @ k = ',num2str(k)]);
        ovflw(k)=1;
    end
elseif pop(k)<kmin;
    if abs(kmin-pop(k))<1e-4*(min(abs(kmin),abs(pop(k))));
        undflw(k)=0;
    else
        warning(['Population underflow @ k = ',num2str(k)]);
        undflw(k)=1;
    end
else
end

end
%% OUTPUTS
% let's put back fluxes in CFU m-2 s-1 by dividing it by the timestep d
conc=conc';
removal=removal./d;
removal=removal';
% let's change sign to net flux to respect conventions, negative fluxes are
% downward, positive fluxes upward. 
Fn=-1*(Fn./d);Fn=Fn';
Fd=Fd./d;Fd=Fd';



outstruc.population=pop;
outstruc.growth=growth;
outstruc.kmax=kmaxout;
outstruc.gross_out=removal;
outstruc.gross_in=Fd;
outstruc.conc=conc;
outstruc.net_flux=Fn;
outstruc.vg=vg;
outstruc.vimp=vimp;
outstruc.vdep=vdep;
outstruc.dflag=dflag;
outstruc.dieout=dieout;
outstruc.overflow=ovflw;
outstruc.underflow=undflw;
end

function vg=settvel(d,T,P)
% Settling velocity evaluation for spherical particles
%
% inputs:
% d = diameter (um)
% T = temperature (K)
% P = pressure (Pa);
%
% output:
% vg = settling velocity (m/s)

% set constants
rhop=1100; % bioaerosol particle density (kg/m3) from Cox CS & Wathes CM (1995), Bioaerosols Handbook
eta=1.833e-5; % air viscosity (Pa*s)
a1=1.257; % Cunningham slip correction factor 1
a2=0.4; % Cunningham slip correction factor 2
a3=1.1; % Cunningham slip correction factor 3
lr= 0.0665; % air reference mean free path (um)
g=9.78; % gravity accel. (m/s^2)

% calculations
kP=P./1000;
l=lr.*(101./kP).*(T./293).*((1+(110/293))./(1+(110./T)));
Cc=1+((2.*l)./d).*(a1+a2.*exp((-a3.*d)./(2.*l)));
dm=d./1e6;
vg=(g.*rhop.*dm.^2.*Cc)./(18.*eta);
end

function vimp=slinnmod(wsp,ust,dp,vg,zc,z0)

% inputs
% wind speed (m/s)
% ustar (m/s)
% particle diameter (micrometers)
% gravitational settling velocity (m/s)
% canopy height (m)
% roughness length (z0)

% model parameters
asmall=10e-06; % diameter of biggest interceptor (m)
abig=1e-03; % diameter of smaller interceptor (m)
cvcd=1/3; % ratio between the portion of cd arising from viscous drag (cv) and the average drag coefficient for vegetation at height z (cd)
F=0.01; % fraction of the total collected momentum collected by smallest collectors (e.g.: vegetative hairs)
b=2; % coefficient for the rebound factor R

% constants
vonk=0.4;

% computations
l=zc;
gam=(l.*100)^(1/2);
if dp<1;
    char1=char(8805);
    char2=char(956);
    error(['Slinn Model has been implemented without diffusional effects and is reliable for particles ',char1,' 1 ',char2,'m']);
else
    dp=dp.*10^-6;
end

Cd=ust.^2./wsp.^2;
uhur=(ust./(vonk.*wsp)).*log(l./z0);
St=(vg.*ust)./(abig.*9.78);
Ein=cvcd.*(F.*(dp./(dp+asmall))+(1-F).*(dp./(dp+abig)));
Eim=St./(1+St.^2);
R=exp(-b.*sqrt(St));
epsilon=(Ein+Eim).*R;
vimp=Cd.*wsp.*(1+uhur.*((1-epsilon)./(epsilon+sqrt(epsilon).*tanh(gam).*sqrt(epsilon)))).^-1;
end
