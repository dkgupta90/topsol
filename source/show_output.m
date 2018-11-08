%FUNCTION TO SHOW FINAL OUTPUT AND SAVE DATA IN FILES
function dummy = show_output(nelx, nely, U, j0, jV, P, beta, Vbus, elA, c, filename, Mov)
 ss=get(0,'screensize');
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
%title('Objective as function of iteration step')
xlabel('Iterations')
ylabel('Efficiency %')
print(4, '-depsc', [filename 'Convergence'])
% print -depsc [filename 'Convergence']
figure(5)
%title('Density as function of iteration step')
xlabel('Iterations')
ylabel('V_{bus}');
print(5, '-depsc', [filename 'density'])
movie2avi(Mov, 'filename.avi', 'compression', 'none','fps',3);