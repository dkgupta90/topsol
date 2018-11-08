%##########################################################################
%TopSol: Topology Optimizimation for front electrode pattern of solar cells
%**************************************************************************
%MATLAB code developed at Structural Optimization Group (SOM), Dept. of
%Precision and Microsystems Engineering, TU Delft
%**************************************************************************
% This code is developed my modifying the code published in:
% E. Andreassen, A. Clausen, M. Schevenels,                                
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010
%##########################################################################
%##########################################################################

function topsol(nelx, nely, volfrac, penal, rmin, ft)

%input parameters
%nelx:  no. of rows
%nely:  no. of columns
%volfrac:   max. allowed volume fraction of space for filing material
%penal: penalty factor for total conductivity matrix
%rmin:  
%ft:    control variable for filter type

%declaration of global variables

%checking the no. of input arguments
clc; 
if nargin > 6
    disp('Wrong no. of inputs');
    exit();
elseif nargin == 0
    nelx = 299;
    nely = nelx;
    volfrac = 1;
    penal = 3.0;
    rmin = 1.5;
    ft = 2;
elseif nargin == 1
    nely = nelx;
    volfrac = 1;
    penal = 3.0;
    rmin = 1.5;
    ft = 2;
elseif nargin == 2
    volfrac = 1;
    penal = 3.0;
    rmin = 1.5;
    ft = 2;
elseif nargin == 3
    penal = 3.0;
    rmin = 1.5;
    ft = 2;
elseif nargin == 4
    rmin = 1.5;
    ft = 2;
elseif nargin == 5
    ft = 2;
end
penalS = 3;                     %penalty factor for shading

disp('START');
 %penal = str2num(getenv('iter1'))/10;
 %penalS = str2num(getenv('iter2'))/10;

mkdir('output/output_penal_penalS');
output_path = 'output/output_penal_penalS/';
folder_new = ['Output_', num2str(nelx), '_', num2str(nely), '_penal', num2str(penal), '_penalS', num2str(penalS)];
mkdir([output_path, folder_new]);
output_path = [output_path, folder_new, '/'];
%MATERIAL PROPERTIES
%TCO layer
sigmaTL = 1e5;
hTL = 200e-9;
Emin = sigmaTL * hTL;
%electrode layer
sigmaEL = 1e7;
hEL = 10e-6;
E0 = sigmaEL * hEL;

%OPTIMIZATION OPTIONS
shading = 1;                    %control variable for sensitivity
jecorrect = 1;                  %correct current density -> negative
senscheck = 0;                  %control variable for sensitivity analysis
xpert = 0.0001;                 %pertubation for Finite difference check

Qmethod  = 3;                   %1,2 or 3
elec = 2;                       %control variable for initial electrode
nrofe = 991;                      %number of electrodes
objective = 2;                  %1-> average voltage 2-> power
optimizer = 3;                  %1-> MMA 2-> Optimality criteria
nroffixeddofs = 1;              %number of fixed points in the grid
maxiter = 100;                  %maximum number of iterations
objhist = NaN*ones(1,maxiter);  %saves all values of objective function
volhist = NaN*ones(1,maxiter);  %saves all values of the volume fraction
senscheck1 = 0;

%Experimental IV-characteristics (obtained from TNO Eindhoven)
jL = 353;       %Photoillumination current at zero shading (in mA/cm2)
j0 = -0.0029;     %reverse bias current (in mA/cm2)
beta = 19.31;    %describes the p-n junction properties
Uref = 0.6;     %fixed voltage for testing
s_Vbus = 0.79;  %factor that defines Vbus as a function of Uref
Vbus = Uref * s_Vbus;   %voltage for fixed points in the space

%Reading the free form structure
xleaf = load('leaf/leaf1.dat', '-ascii');
%xleaf = load('leaf/leaf_maple2.dat', '-ascii');
xleaf = xleaf(1:860, 225:650);
xleaf(830:860, 195:220) = 100;
%xleaf(348:350, 302:304) = 100;
imagesc(xleaf);
nely = 860 - 1;
nelx = 426 - 1;
%nely = 461 - 1;
%nelx = 615 - 1;
fixeddofs = find(xleaf == 100);
nodofs = find(xleaf == 255);



%Defining design space dimension parameters
Lx = 15e-3;                     %Length in x-direction
elh = Lx/nelx;                  % element size
Ly = nely*elh;                  % Length in y-direction
elA = elh^2;                    % element area

%Create electrodes
ewidth = 14;
evec = create_electrode(ewidth, nrofe, nely);

% Compute shape functions
phi = compute_shape_fn();       %U*phi gives voltages....

%Prepare Finite Element Analysis
KE = finite_elem();

%Generation of nodal matrix (%explained in Andreassen et al. 2011)
Qvec = elA/4*ones(4,1);
Qmat = elA/16*ones(4);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,4)+repmat([0 nely+[1 0] -1],nelx*nely,1);

iK = reshape(kron(edofMat,ones(4,1))',16*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,4))',16*nelx*nely,1);

% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
R = zeros((nely+1)*(nelx+1),1);
Q = zeros((nely+1)*(nelx+1),1);
U = zeros((nely+1)*(nelx+1),1);


%Defining the fixed degrees of freedom
%fixeddofs = fixed_dofs(nelx, nely, nroffixeddofs);
%fixeddofs = [fixeddofs evec+301];
alldofs = [1:(nely+1)*(nelx+1)];
actualdofs = setdiff(alldofs, nodofs);
freedofs = setdiff(alldofs, fixeddofs);
freedofs = setdiff(freedofs, nodofs);

area_factor = length(actualdofs)/length(alldofs);

% %computing the indices of the free elements
% freelems = [];
% for i = 1:1:nelx*nely
%     if(sum(ismember(edofMat(i, :), actualdofs)) == 4)
%         freelems = [freelems i];
%     end
% end
% save freelems_maple.mat freelems;
load freelems.mat;
%load freelems_maple.mat;

noelems = setdiff(1:nely*nelx, freelems);
if senscheck1 == 1
    nsens = 3;
else
    nsens = 1;
end

perturbx = 25;
perturby = 25;
perturb = 1e-4;


%% PREPARE FILTER
[H, Hs] = prepare_filter(nelx, nely, rmin);

%% INITIALIZE ITERATION

%FOR INITIAL DESIGN!
x = initial_design(nelx, nely, elec, 0.5, evec);
Hbeta = 1;
if ft == 1 || ft == 2
    xPhys = x;
elseif ft == 3
    xTilde = x;
    xPhys = 1 - exp(-Hbeta * xTilde) + xTilde*exp(-Hbeta);
end
%Initializing the optimizer parameters
if optimizer == 1
%MMA control parameters
    [f0fac, f0add, m, n, xminvec, xmaxvec, low, upp, c1, d, a0, a] = ...
        intialize_MMA(nelx, nely);
    xold1 = x(:);
    xold2 = x(:);
    xold1(n+1) = s_Vbus;
    xold2(n+1) = s_Vbus;
end

%% START ITERATION

% compute PV currents:
Vavg = .25*sum(U(edofMat),2);
if shading == 1
  Xs = 1-reshape(xPhys,nelx*nely,1);          %densities for shading
else
  Xs = ones(nelx*nely,1);                     %no density correction
end

%Xs(noelems) = 0.0;

loop = 0;
Hbetaloop = 0;
[je, je_check] = compute_je(Qmethod, shading, Xs, penalS, j0, jL, beta, Vavg,...
    edofMat, phi, U, loop, jecorrect);

for i=1:size(edofMat,1)
    Q(edofMat(i,:))=Q(edofMat(i,:))+je(i)*Qvec;
end

change = 1;
changemin = 0.01;
while change > changemin && loop < maxiter
  loop = loop + 1;
  Hbetaloop = Hbetaloop + 1;
%   if penal  > 1.2 && loop > 20
%       penal = penal - 0.01;
%   end
  for i1 = 1:1:nsens
%Initializing the optimizer parameters
    if i1 == 1
        xPhys(perturbx, perturby) = xPhys(perturbx, perturby) - perturb;
        Q11 = Q;
        U1 = U;
    elseif i1 == 2
        xPhys(perturbx, perturby) = xPhys(perturbx, perturby) + 2*perturb;
        Q = Q11;
        U = U1;
    else
        Q = Q11;
        U = U1;
        xPhys(perturbx, perturby) = xPhys(perturbx, perturby) - perturb;
    end
    if shading == 1
      Xs = 1-reshape(xPhys,nelx*nely,1);          %densities for shading
    else
      Xs = ones(nelx*nely,1);                     %no density correction
    end
  %% FE-ANALYSIS
  Edist = Emin+xPhys(:)'.^penal*(E0-Emin);
  Edist(noelems) = 1e-6;
  sK = reshape(KE(:)*(Edist),16*nelx*nely,1);

  %start newton iteration
  NRiter = 0;
  Qnorm = norm(Q);
  NRiterMax = 10;
  NRtol = (nelx*nely)*1e-10;
  Kv = sparse(iK,jK,sK); Kv = (Kv+Kv')/2;
  
  while 1
    NRiter = NRiter +1;
    Qnorm = norm(Q);
    if NRiter > NRiterMax, break; end    %loop control statement
    U(fixeddofs) = Vbus; %setting boundary values
    % compute residual
        R = Kv*U - Q;
    
    %norm(R(freedofs))
    if norm(R(freedofs))/Qnorm<(NRtol/100), 
%     if norm(R)/Qnorm<(NRtol/100),     
        disp(['Convergence reached in ',num2str(NRiter),' iterations']);
        break 
    else
        %disp(['convergence: ' num2str(norm(R(freedofs))/Qnorm)]);
    end
    dJdV = compute_dJdV(Qmethod, j0, beta, Vavg, edofMat, phi,...
        jecorrect, je_check, U, []); 
    sKQ = reshape(Qvec(:)*dJdV(:)',16*nelx*nely,1);    
    dQdv = sparse(iK,jK,sKQ);
    K = Kv - dQdv;           %This K equals dR/dU
    
    % solve increment and update U vector;
    deltaU = K(freedofs,freedofs)\R(freedofs);
    U(freedofs) = U(freedofs) - deltaU;
    clear deltaU
    
    % compute PV currents:
    Q = zeros((nely+1)*(nelx+1),1); 
    Vavg = .25*sum(U(edofMat),2);
    % je = j0+jV*(exp(beta*Vavg)-1);
    
    [je, je_check] = compute_je(Qmethod, shading, Xs, penalS, j0, jL, beta, Vavg,...
    edofMat, phi, U, loop, jecorrect);
     
    for i=1:size(edofMat,1)
    Q(edofMat(i,:))=Q(edofMat(i,:))+je(i)*Qvec;
    end  
  end  
  if senscheck1 == 1
    fd(i1) = -(Uref*s_Vbus)*sum(Q);
  end
  end
  %%computing dQfdvf and dQfdvp

  if senscheck1 == 1
    fdsens = (fd(2) - fd(1))/(2*perturb)
  end  
  
%   if loop == 201
%       disp('HALT');
%   end
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
[c, dc, dv, lambda] = sens_analysis(objective, freedofs, U, K, Q, jL, edofMat, KE, nelx, nely,...
     Emin, xPhys, penal, E0, dQdv, shading, penalS, Xs, Qvec, Vbus, elA);
 
 dc(noelems) = 0;
  
dc_Vbus = sens_analysis_Vbus(lambda, Qmethod, j0, beta, Vavg, edofMat,...
    phi,jecorrect, je_check, U, Uref, Q, s_Vbus, fixeddofs, freedofs,...
    Kv, Qvec, nelx, nely, iK, jK);
  dc(25, 25)
    %% finite difference check on sensitivities
  if senscheck
  %do finite difference analysis
  else
  %TEMPORARY

  if max(dc(:))>0 && optimizer == 0
          disp(sprintf(' >>>>>  Artificial sensitivity correction for OC: %e\n %5i',max(dc(:))/(sum(abs(dc(:)))/length(dc)),length(find(dc>0))));
          dc = dc - 2*max(0,max(dc(:)));
  end
  
%% FILTERING/MODIFICATION OF SENSITIVITIES
if ft == 1 || ft == 2
    [dc, dv] = filtering(ft, H, x, dc, dv, Hs, 0.0);
elseif ft == 3
    [dc, dv] = filtering(ft, H, xTilde, dc, dv, Hs, Hbeta);
end
  %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
  if optimizer == 0
    [x, xPhys] = optim_criteria(nelx, nely, x, Emin, dc, dv,...
        xPhys, H, Hs, volfrac, ft);
  end
  %% MMA - Method of Moving Asymptotes
  if optimizer == 1
  f0val = c*f0fac + f0add;
  df0dx = dc(:)*f0fac;
  df0dx(n+1) = dc_Vbus*f0fac;
  df0dx2 = 0*df0dx;
  fval = sum(xPhys(:))/(n*volfrac) - 1 ;
  dfdx = dv(:)/(n*volfrac);
  dfdx(n+1) = 0;
  dfdx2 = 0*dfdx';
  iter = loop;
  xval = x(:);
  xval(n+1) = s_Vbus;
  xminvec(n+1) = 0.30;
  xmaxvec(n+1) = 0.8;
  [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
      mmasub(m,n+1,iter, xval, xminvec, xmaxvec, xold1, xold2, ...
      f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a,c1,d);
  xmma(noelems) = 0;
  
  if loop == maxiter
      xmma(noelems) = 0.5;
  end
  xold2 = xold1;
  xold1 = xval;      
  xnew = reshape(xmma(1:n),nely,nelx);
  

  if ft == 1
      xPhys = xnew;
  elseif ft == 2
      xPhys(:) = (H*xnew(:))./Hs;
  elseif ft == 3
      xTilde(:) = (H*xnew(:))./Hs;
      xPhys = 1 - exp(-Hbeta*xTilde) + xTilde*exp(-Hbeta);
  end
   x = xnew;
  end
 
  if optimizer == 2
      xTilde = (H*xPhys(:))./Hs;
      xTilde = reshape(xTilde, nely, nelx);
      xPhys = 1 - exp(-Hbeta*xTilde) + xTilde*exp(-Hbeta);
  end
  
% if loop == 200
%       xPhys(xPhys > 0.1) = 1;
%       xPhys(xPhys < 0.1) = 0;
%       x(x > 0.1) = 1;
%       x(x < 0.1) = 0;
%   end

        if optimizer == 1
          s_Vbus = xmma(n+1);
        end
  end

  %% Compute cell power CHANGE!!!
  P = sum(elA*je)*Vbus;               % Power of the cell = sum of all currents*voltage
  Pe = je.*Vavg*elA;    % Local Power
  Pe_eff = je*Vbus*elA; % Contribution to total power // Effective power
  Pe_loss = Pe-Pe_eff;                % Power loss
  Eff = ((P/(Lx*Ly))/1000)*100*area_factor;
   Vbushist(loop) = Vbus;
  objhist(loop) = Eff;
  volhist(loop)=mean(xPhys(:));
  Vavghist(loop) = sum(U)/length(U);
  
  Vbus = s_Vbus*Uref;
  
  
  %% loop condition
  if loop > 1
  change = 50*abs((P-objhist(loop-1))/(P));
  else
  change = 1;
  end
  x_low = 0.001;
  gray_factor = get_gray_factor(xPhys, nelx, nely, x_low);
  
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f GF.:%7.3f ch.:%7.3f Pow.:%7.3f\n',loop,c, ...
    mean(xPhys(:)),gray_factor, change,P);
%   %% PLOT DENSITIES
   Mov(loop) = plot_densities(nelx, nely, xPhys, loop, P, Pe,...
       Pe_loss, Vavg, je, objhist, volhist, Vbus, Eff, Vbushist);
    
    if ft == 3 && Hbeta < 2048 && (Hbetaloop >=50 || change <= changemin)
        Hbeta = 2*Hbeta;
        Hbetaloop = 0;
        change = 1;
        Hbeta
  end
  
  if shading == 1
    Xs = 1-reshape(xPhys,nelx*nely,1);          %densities for shading
  else
    Xs = ones(nelx*nely,1);                     %no density correction
  end
  
end

  save_output(nelx, nely, U, Q, xPhys, jL, j0, P, Eff,  beta, Vbus, elA, c, output_path); 
 %% final output display and saving output in files
%show_output(nelx, nely, U, jL, j0, P, beta, Vbus, elA, c, output_path, Mov);

