function [x, xPhys] = optim_criteria(nelx, nely, x, Emin, dc, dv, xPhys, H, Hs,...
    volfrac, ft)

l1 = 0; l2 = 100000; move = 0.2;
while (l2-l1)/(l1+l2) > 1e-3
    lmid = 0.5*(l2+l1);
    xnew = max(Emin,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));  %% CHANGED Emin io 0
    
    if ft == 1
        xPhys = xnew;
    elseif ft == 2
        xPhys(:) = (H*xnew(:))./Hs;
    end
    
    if sum(xPhys(:)) > volfrac*nelx*nely
        l1 = lmid; 
    else
        l2 = lmid; 
    end
end
x = xnew;

end