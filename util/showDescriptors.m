close all
dloc = 200;
fmap = Fmaps{1}.ppMaps{1};
fpushed = transferFunction(f, basis1, basis2, fmap);
h1=showDescriptor(mesh1, dloc, f, 'f'); c1 = caxis;
h2=showDescriptor(mesh2, dloc, g, 'g'); c2 = caxis;
h3=showDescriptor(mesh2, dloc, fpushed, 'f pushed'); c3 = caxis;
c = min([c1;c2;c3]);
C = max([c1;c2;c3]);
cnew = [c(1) C(2)];
%figure(h1); caxis(cnew); figure(h2); caxis(cnew);figure(h3); caxis(cnew);

