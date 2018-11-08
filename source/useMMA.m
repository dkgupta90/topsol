function [x, xPhys, xold1, xold2, iter] = useMMA(x, xPhys, xold1,...
    xold2, nelx, nely, volfrac, loop, dc, dv, f0fac, f0add, H, Hs, ...
    m, n, xminvec, xmaxvec, low, upp, a0, a, c, d, ft)
  
  f0val = c*f0fac + f0add;
  df0dx = dc(:)*f0fac;
  df0dx2 = 0*df0dx;
  fval = sum(xPhys(:))/(n*volfrac) - 1 ;
  dfdx = dv(:)/(n*volfrac);
  dfdx2 = 0*dfdx';
  iter = loop;
  xval = x(:);

  [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
      mmasub(m,n,iter, xval, xminvec, xmaxvec, xold1, xold2, ...
      f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a,c,d);

  xold2 = xold1;
  xold1 = xval;      
  xnew = reshape(xmma,nely,nelx);

  if ft == 1
      xPhys = xnew;
  elseif ft == 2
      xPhys(:) = (H*xnew(:))./Hs;
  end
  x = xnew;
  
end