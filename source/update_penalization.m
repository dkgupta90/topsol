function [penal, penalS] = update_penalization(penalmax, loop, maxiter)

dpenal = penalmax - 1.5;
dmaxiter = maxiter - 900;
pfactor = min(1, loop/dmaxiter);

penal = 1.5 + (dpenal*pfactor);
penalS = penal;

end