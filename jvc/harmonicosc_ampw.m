function A=harmonicosc_ampw(w,m,c,k,g,w0,f);

A=f./m./sqrt((w0.^2-w.^2).^2 + 4.*g.^2.*w.^2);

