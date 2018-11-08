%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
%%%% Original 88 line code, adapted for heat conductivity problem
%%%% Adapted to solve for the electrical problem
%%%% Added IV Curves and bus voltage
%%%% Added shading
%%%% INITIAL DESIGN CHANGE / Start with electrode / PUM holes
%%%% Solve for problem Lx = 2*Lx -> matrix bad conditioned (Je<0 / Je<-J0)
%%%% Including correct sensitivity data
%%%% Use gauss sample points/averaged current/averaged voltage
%%%% Adapted to use MMA / OC
%%%% New objective implemented (Power)
%%%% Sensitivity Check using FD / Central FD
%top88PVCell(300,100,0.1,3,1.5,1)

function top88PVCell(nelx,nely,volfrac,penal,rmin,ft)
%% MATERIAL PROPERTIES

filename = 'deepak_test3';               %mapnaam

sigmaTL = 1e5;
sigmaEL = 1e7;
hTL = 200e-9;
hEL = 10e-6;
E0 = hEL*sigmaEL;               %refers the electrode
Emin = hTL*sigmaTL;             %refers the transparent conducting layer       
Lx = 15e-3;                     %length in x direction added
%Lx = 2*Lx;

%Optimisation options
shading = 1;                    %control variable for sensitivity
jecorrect = 1;                 %correct current density -> negative
senscheck = 0;                  %control variable for sensitivity analysis
xpert = 0.0001;                 %pertubation for Finite difference check
penalS = 3;                     %penalty factor for shading
Qmethod  = 3;                   %1,2 or 3
elec = 2;                       %control variable for initial electrode
nrofe = 2;                      %number of electrodes
objective = 2;                  %average voltage 1) or power 2);
useMMA = 1;                     %controls use of MMA Optimizer
nroffixeddofs = 2;              %number of fixed points in the grid  

%MMA CONTROL
f0fac = 100;          % objective scaling --> 1
f0add = 0;          % constant addition to objective

if useMMA
    % Initialize MMA parameters
    m = 1;                  % n of constraints / volume constraint
    n = nelx*nely;          % number of variables
    xmin = 0;               % minimum value of density (now 0)

    xminvec  = xmin*ones(n,1);  %Column vector with the lower bounds for the variables x_j.
    xmaxvec  = ones(n,1);       %Column vector with the upper bounds for the variables x_j.
    low   = xminvec;            %Column vector with the lower asymptotes from the previous  iteration (provided that iter>1).
    upp   = xmaxvec;            %Column vector with the upper asymptotes from the previous
    
    c = 10000*ones(m,1);        %Column vector with the constants c_i in the terms c_i*y_i.
    d = zeros(m,1);             %Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
    a0 = 1;                     %The constants a_0 in the term a_0*z.
    a = zeros(m,1);             %Column vector with the constants a_i in the terms a_i*z
end

maxiter = 250;

%history of values as function of the loop
objhist = NaN*ones(1,maxiter);
volhist = NaN*ones(1,maxiter);


%IV-Curve properties
j0 = 310;
jV = -0.006;
beta = 16.4;

Vbus = 0.5;                    %V, fixed voltage for testing...
% Vbus = 0.5239;
% Vbus = 0.52;


%including domain size
elh = Lx/nelx;                  % element size
Ly = nely*elh;                  % Length in y-direction
elA = elh^2;                    % element area

%% Compute shape functions


Xa = [0.211324865405187 0.788675134594813 0.788675134594813 0.211324865405187];
Yb = [0.211324865405187 0.211324865405187 0.788675134594813 0.788675134594813];
a = 1 ; b = 1;

phi1 = (1-Xa./a).*(Yb./b);
phi2 = (1-Xa./a).*(1-Yb./b);
phi3 = (Xa./a).*(Yb./b);
phi4 = (Xa./a).*(1-Yb./b);

%put them in the order 2 4 3 1!!
phi = [phi2;phi4;phi3;phi1];               %U*phi gives voltages....

%% CREATE ELECTRODES
ewidth = 14;
elnr = nely/(nrofe+1);
evec = [((-1/2*ewidth)+1):(1/2*ewidth)]+elnr;
if nrofe == 2;
elnr = nely/(nrofe*2);
evec = [((-1/2*ewidth)+1):(1/2*ewidth)]+elnr;
elnr2 = 3*nely/(nrofe*2);
evec2 = [((-1/2*ewidth)+1):(1/2*ewidth)]+elnr2;
evec = [evec evec2];
end
%% PREPARE FINITE ELEMENT ANALYSIS
A = [ 2/3 -1/6
     -1/6  2/3];
B = [-1/3 -1/6 
     -1/6 -1/3]; 
KE = [A B;B A];

Qvec = elA/4*ones(4,1);
Qmat = elA/16*ones(4);
    
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,4)+repmat([0 nely+[1 0] -1],nelx*nely,1);
%explained in Andreassen et al. 2011

iK = reshape(kron(edofMat,ones(4,1))',16*nelx*nely,1);

jK = reshape(kron(edofMat,ones(1,4))',16*nelx*nely,1);


% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
% F = sparse((nely+1)*(nelx+1),1);
% F(:,1) = 0.01;
R = zeros((nely+1)*(nelx+1),1);
Q = zeros((nely+1)*(nelx+1),1);
U = zeros((nely+1)*(nelx+1),1);

if nroffixeddofs == 1
fixeddofs = [1:nely+1];                     % left column
end
% fixeddofs = [nely-3:nely+1];                  % single point corner
% fixeddofs = [round(0.95*nely):nely+1];    %single point corner

if nroffixeddofs ==2
fixeddofs = round((nely+1)/2) + (nely+1)*round((nelx+1)/2);
end

if nroffixeddofs == 4
    fixeddofs = round((nely+1)/4) + (nely+1)*round((nelx+1)/4);
    fixeddofs = [fixeddofs fixeddofs(1)+round((nely+1)/2)];
    fixeddofs = [fixeddofs fixeddofs + (nely+1)*round((nelx+1)/2)];
end

if nroffixeddofs == 5 %one fixed busbar point at the left centre
    fixeddofs = [round((nely+1)/2)-10:round((nely+1)/2)+10];
end

if nroffixeddofs == 6
    fixeddofs = [nely-10:nely+1];
end

alldofs = [1:(nely+1)*(nelx+1)];            % can be faster by using nodenrs(:)' 
freedofs = setdiff(alldofs,fixeddofs);

U(fixeddofs) = Vbus;                        % set voltage

%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0.001,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION

%FOR INITIAL DESIGN!
if elec == 1
x = zeros(nely,nelx);
x(evec,:) = 1;
else
x = repmat(volfrac,nely,nelx);
end

xPhys = x;
loop = 0;
change = 1;
if senscheck
c_old = 0;  % needed for FD-Analysis
x_old = zeros(nely,nelx);  % needed for FD-Analysis
x_old2 = zeros(nely,nelx);
c_old2 = 0;
cstar_old = 0;
cstar_old2 = 0;
dc_old = ones(nely,nelx);
method1 = 1;   %pertubation on 1 value of x
p1 = 15;
p2 = 30;
end

if useMMA
    xold1 = x;
    xold2 = x;
end

NRiterMax = 30;
NRtol = (nelx*nely)*1e-10;
%% START ITERATION

% compute PV currents:
Vavg = .25*sum(U(edofMat),2);
if shading == 1
  Xs = 1-reshape(xPhys,nelx*nely,1);          %densities for shading
else
  Xs = ones(nelx*nely,1);                     %no density correction
end

if Qmethod == 1
    if shading
    je = Xs.^penalS*j0+jV*(exp(beta*Vavg)-1);                 %current density for 1 node
    else
    je = j0+jV*(exp(beta*Vavg)-1);
    end
elseif Qmethod == 2
    if shading
    je = Xs.^penalS*j0+1/4*sum(jV*(exp(beta*U(edofMat))-1),2);                 %current density for 1 node
    else
    je = j0+1/4*sum(jV*(exp(beta*U(edofMat))-1),2);                 %current density for 1 node
    end
else
    if shading
    je = Xs.^penalS*j0+1/4*sum(jV*(exp(beta*U(edofMat)*phi)-1),2);                 %current density for 1 node
    else
    je = j0+1/4*sum(jV*(exp(beta*U(edofMat)*phi)-1),2);                 %current density for 1 node
    end
end

je_check = je < -j0;

if any(je_check)
    disp('Voltage above upper limit, solution not correct')
    if jecorrect == 1
        je(je_check) = 0;
        disp('Correction on current density applied')
    end
end


for i=1:size(edofMat,1)
    Q(edofMat(i,:))=Q(edofMat(i,:))+je(i)*Qvec;
end




while change > 0.01 && loop < maxiter
  loop = loop + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),16*nelx*nely,1);

  %start newton iteration

  NRiter = 0;
  Qnorm = norm(Q);
  
  Kv = sparse(iK,jK,sK); Kv = (Kv+Kv')/2;
  
  while 1
    NRiter = NRiter +1;
    Qnorm = norm(Q);
%     
%     if (NRiter > NRiterMax) && (loop > 10)
%         error('Maximum nr of iterations reached'); 
%     elseif NRiter > NRiterMax
%         break
%     end

    if NRiter > NRiterMax, break; end
   
    
    % compute residual
        R = Kv*U - Q;
    
    %norm(R(freedofs))
    if norm(R(freedofs))/Qnorm<(NRtol/100), 
%     if norm(R)/Qnorm<(NRtol/100),     
        disp(['Convergence reached in ',num2str(NRiter),' iterations']);
        break 
    else
        disp(['convergence: ' num2str(norm(R(freedofs))/Qnorm)])
    end
    % build dRdV, solve for voltage update:
    
    if Qmethod == 1
    dJdV = repmat(1/4*(jV*beta*exp(beta*Vavg(:)')),4,1);
    elseif Qmethod == 2;
    dJdV = 1/4*(jV*beta*exp(beta*U(edofMat)))';
    else
    dJdV = 1/4*phi*(jV*beta*exp(beta*U(edofMat)*phi))';
    end
    
    if jecorrect == 1
    dJdV(:,je_check) = 0;                                       %makes derivative zero if je < 0  trick!S       
    end
    
    
%     sKQ = reshape(Qmat(:)*dJdV,16*nelx*nely,1);  
    sKQ = reshape(Qvec(:)*dJdV(:)',16*nelx*nely,1);    
    dQdv = sparse(iK,jK,sKQ);
    K = Kv - dQdv; %K = (K+K')/2;           %This K equals dR/dU
    
    % solve increment and update U vector;
    deltaU = K(freedofs,freedofs)\R(freedofs);
%     deltaU = K\R;
    U(freedofs) = U(freedofs) - deltaU;
    clear deltaU
    
    % compute PV currents:
    Q = zeros((nely+1)*(nelx+1),1); 
    Vavg = .25*sum(U(edofMat),2);
    % je = j0+jV*(exp(beta*Vavg)-1);
    
if Qmethod == 1
    if shading
    je = Xs.^penalS*j0+jV*(exp(beta*Vavg)-1);                 %current density for 1 node
    else
    je = j0+jV*(exp(beta*Vavg)-1);
    end
elseif Qmethod == 2
    if shading
    je = Xs.^penalS*j0+1/4*sum(jV*(exp(beta*U(edofMat))-1),2);                 %current density for 1 node
    else
    je = j0+1/4*sum(jV*(exp(beta*U(edofMat))-1),2);                 %current density for 1 node
    end
else
    if shading
    je = Xs.^penalS*j0+1/4*sum(jV*(exp(beta*U(edofMat)*phi)-1),2);                 %current density for 1 node
    else
    je = j0+1/4*sum(jV*(exp(beta*U(edofMat)*phi)-1),2);                 %current density for 1 node
    end
end
    
    disp(['Highest voltage found during iteration: ' num2str(max(U)) ]) 
    disp(['Lowest current found during iteration: ' num2str(min(je)) ]) 
    
    je_check = je < -j0;
    
    if any(je_check)
    disp(['Voltage above upper limit, newton iteration not correct, iteration'...
           num2str(NRiter)])
        if jecorrect ==1
        je(je_check) = 0;
        disp('Correction on current density applied')
        end
    end

    
    
    for i=1:size(edofMat,1)
    Q(edofMat(i,:))=Q(edofMat(i,:))+je(i)*Qvec;
    end
   
    % KLEINE TEST!!!!
%     Q(fixeddofs) = 0;
    
  end    
  
  
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  lambda = zeros((nely+1)*(nelx+1),1);
  
  if objective == 1; 
  lambda(freedofs) = -K(freedofs,freedofs)\(dQdv(freedofs,freedofs)*U(freedofs)+Q(freedofs));   %(compute lambda) K = dRdU
  ce = reshape(sum((U(edofMat)*KE).*lambda(edofMat),2),nely,nelx);
  ce2 = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  cstar = (Emin+xPhys.^penal*(E0-Emin)).*ce2;
  c = Q.'*U;
%   c = U.'*Kv*U;
  dc = penal*(E0-Emin)*xPhys.^(penal-1).*ce;
  
  if shading
      dshading = - penalS*Xs.^(penalS-1);
      dc =  -(dc + reshape(j0.*dshading.*((U(edofMat) - lambda(edofMat))*Qvec),nely,nelx));
  end
  
  dv = ones(nely,nelx);
  
  else
  % c = Power = Vbus*elA*(sum of all currents);
  
  lambda(freedofs) = -K(freedofs,freedofs)\(Vbus*sum(dQdv(freedofs,freedofs),1).');   %(compute lambda) K = dRdU
  
  ce = reshape(sum((U(edofMat)*KE).*lambda(edofMat),2),nely,nelx);
  ce2 = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  cstar = (Emin+xPhys.^penal*(E0-Emin)).*ce2;
  
  c = Vbus*sum(Q);

  dc = penal*(E0-Emin)*xPhys.^(penal-1).*ce;
  
  if shading
      dshading = - penalS*Xs.^(penalS-1);
      dc =  -(dc + reshape(j0.*dshading.*((-lambda(edofMat))*Qvec),nely,nelx) + Vbus*reshape(j0.*dshading,nely,nelx)*elA);
  end
  
  dv = ones(nely,nelx);    
  end
  
  
    %% finite difference check on sensitivities
  if senscheck
  fd = cstar - cstar_old;    % finite difference
  step = xPhys-x_old;
  maxstepsize = max(abs(step(:)));
  minstepsize = min(abs(step(:)));
  dc_fd = fd./step; 
  dc_fdcentral = (cstar-cstar_old2)./(xPhys-x_old2);
  
  Exact = dc(p1,p2)
  Exactprevious = dc_old(p1,p2)
  
  if method1
  finiteDiff = (c-c_old)/(x(p1,p2) - x_old(p1,p2))
  finiteCentral = (c-c_old2)/(x(p1,p2) - x_old2(p1,p2))
  disp('methode 1 toegepast')
  else
  finiteDiff = dc_fd(p1,p2)
  finiteCentral = dc_fdcentral(p1,p2)
  end
  
  ratioDiff = finiteDiff/Exact
  ratiocentral = finiteCentral/Exactprevious
  
  x_old2 = x_old;
  c_old2 = c_old;
  cstar_old2 = cstar_old;
  
  dc_old = dc;
  x_old = xPhys;    % new value for x_old
  c_old = c;       % new value for c_old
  cstar_old = cstar;
  
  
  if method1
  xPhys(p1,p2) = x(p1,p2) + xpert;
  else
  xPhys = x + xpert;
  end     
  x = xPhys;
  pause
  else
  
  
  %TEMPORARY

  if max(dc(:))>0 && ~ useMMA
          disp(sprintf(' >>>>>  Artificial sensitivity correction for OC: %e\n %5i',max(dc(:))/(sum(abs(dc(:)))/length(dc)),length(find(dc>0))));
          dc = dc - 2*max(0,max(dc(:)));
  end
  
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  if ft == 1
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  elseif ft == 2
    dc(:) = H*(dc(:)./Hs);
    dv(:) = H*(dv(:)./Hs);
  end
  
  %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
  if ~useMMA
  l1 = 0; l2 = 100000; move = 0.2;
  while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    xnew = max(Emin,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));  %% CHANGED Emin io 0
%     xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));  %% CHANGED Emin io 0
    if ft == 1
      xPhys = xnew;
    elseif ft == 2
      xPhys(:) = (H*xnew(:))./Hs;
    end
    if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
  end
  change = max(abs(xnew(:)-x(:)));
  x = xnew;
  end
  
  %% MMA - Method of Moving Asymptotes
  if useMMA
      f0val = c*f0fac + f0add;
      df0dx = dc(:)*f0fac;
      df0dx2 = 0*df0dx;
      fval = sum(xPhys(:))/(n*volfrac) - 1 ;
      dfdx = dv(:)/(n*volfrac);
%       fval = 0;
%       dfdx = 0;
      dfdx2 = 0*dfdx';
      iter = loop;
      xval = x(:);
      
      [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
          mmasub(m,n,iter,xval,xminvec,xmaxvec,xold1,xold2, ...
          f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a,c,d);
      
      xold2 = xold1;
      xold1 = xval;      
      xnew = reshape(xmma,nely,nelx);
      
      if ft == 1
          xPhys = xnew;
      elseif ft == 2
          xPhys(:) = (H*xnew(:))./Hs;
      end
      change = max(abs(xnew(:)-x(:)));
      x = xnew;
  end
  


  %% Compute cell power CHANGE!!!
  P = sum(elA*je)*Vbus;               % Power of the cell = sum of all currents*voltage
  Pe = je.*Vavg*elA;    % Local Power
  Pe_eff = je*Vbus*elA; % Contribution to total power // Effective power
  Pe_loss = Pe-Pe_eff;                % Power loss
  
  objhist(loop) = P;
  volhist(loop)=mean(xPhys(:));
  
  %% loop condition
  if loop > 1
  change = 500*abs((P-objhist(loop-1))/(P));
  else
  change = 1;
  end
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f Pow.:%7.3f\n',loop,c, ...
    mean(xPhys(:)),change,P);
  %% PLOT DENSITIES
ss=get(0,'screensize');
figure(1); %set(1,'position',[6 2*ss(4)/4+20 ss(3)/3 ss(4)/3]);
colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off;
title(['Iteration: ',num2str(loop),';  Power: ', num2str(P)]); drawnow;


figure(2); set(2,'position',[ss(3)/3 2*ss(4)/4+20 ss(3)/3 ss(4)/3]);   
subplot(2,1,1);
[mx,my]=meshgrid(0:nelx,0:nely);
imagesc(reshape(Vavg,nely,nelx)); axis equal; axis tight;colorbar; %shading interp;%interpolatie --> shape function;
subplot(2,1,2);
imagesc(reshape(je,nely,nelx)); axis equal; axis tight;colorbar;drawnow

figure(3);set(3,'position',[2*ss(3)/3 2*ss(4)/4+20 ss(3)/3 ss(4)/3]); 
subplot(2,1,1);
[mx,my]=meshgrid(0:nelx,0:nely);
imagesc(reshape(Pe,nely,nelx)); axis equal; axis tight;colorbar; axis([0 nelx 0 nely]);drawnow
subplot(2,1,2);
imagesc(reshape(Pe_loss,nely,nelx)); axis equal; axis tight;colorbar; axis([0 nelx 0 nely]);drawnow;

figure(4); set(4,'position',[6 50 ss(3)/3 ss(4)/3]);
plot(objhist(1:loop))

figure(5); set(5,'position',[ss(3)/3 50 ss(3)/3 ss(4)/3]); 
plot(volhist(1:loop))

fig=figure(1);
set(gcf,'DoubleBuffer','on'); % double buffer
set(gcf,'PaperPositionMode','auto');
M(loop)=getframe(gcf);

  end
  
  if shading == 1
    Xs = 1-reshape(xPhys,nelx*nely,1);          %densities for shading
  else
    Xs = ones(nelx*nely,1);                     %no density correction
  end
  
end

 %% post
 jemax = j0+jV*(exp(beta*Vbus)-1);   % Current density if cell has zero resistance
 Pmax = elA*jemax*(nelx*nely)*Vbus   % Maximum power possible (withouth shading etc)
 Umax = max(U);                      % maximum voltage found
 je_Umax = j0+jV*(exp(beta*Umax)-1);
  
 disp(['Maximum possible power from cell: Pmax: ' num2str(Pmax)])
 disp(['Fraction obtained: P/Pmax :  ' num2str(P/Pmax)])
 disp(['Fraction obtained: Pshading/Pmax :  ' num2str(c/Pmax)])
 disp(['Maximum Voltage found: ' num2str(max(U))])
 disp(['Corresponding current: ' num2str(je_Umax)])
 
%  save figures
save([filename 'workspace'])
figure(1); set(1,'position',[6 2*ss(4)/4+20 ss(3)/3 ss(4)/3]);
% print -depsc [filename 'design']
print(1, '-depsc', [filename 'design'])
figure(2)
print(2, '-depsc', [filename 'voltage_current'])
% print -depsc [filename 'voltage_current']
figure(3)
print(3, '-depsc', [filename 'Power'])
% print -depsc [filename 'Power']
figure(4)
title('Objective as function of iteration step')
print(4, '-depsc', [filename 'Convergence'])
% print -depsc [filename 'Convergence']
figure(5)
title('Density as function of iteration step')
print(5, '-depsc', [filename 'density'])

% print -depsc [filename 'density']

movie2avi(M,'filename.avi', 'compression', 'none','fps',3);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Andreassen, A. Clausen, M. Schevenels,%
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,           %
%  Technical University of Denmark,                                        %
%  DK-2800 Lyngby, Denmark.                                                %
% Please sent your comments to: sigmund@fam.dtu.dk                         %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
% This version is based on earlier 99-line code                            %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

