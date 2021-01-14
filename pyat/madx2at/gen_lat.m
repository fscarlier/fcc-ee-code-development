clear
format long

global NUMDIFPARAMS
NUMDIFPARAMS.XYStep = 1e-9;
NUMDIFPARAMS.DPStep = 1e-9;

v = load('fcch.mat');
ring = v.ring;

cav = findcells(ring,'Class','RFCavity');
for i=length(cav):-1:1
    if ring{cav(i)}.Voltage==0
        ring{cav(i)}=atmarker('CAV_marker');
    end
end

cav = findcells(ring,'Class','RFCavity');
dip = findcells(ring,'Class','Bend');
quad = findcells(ring,'Class','Quadrupole');
sext = findcells(ring,'Class','Sextupole');
ring(dip) = atsetfieldvalues(ring(dip),'NumIntSteps',20);
ring(quad) = atsetfieldvalues(ring(quad),'NumIntSteps',20);
ring(sext) = atsetfieldvalues(ring(sext),'NumIntSteps',20);

ring(dip) = atsetfieldvalues(ring(dip),'EntranceAngle',0);
ring(dip) = atsetfieldvalues(ring(dip),'ExitAngle',0);

volt = sum(atgetfieldvalues(ring(cav),'Voltage'));
freq = mean(atgetfieldvalues(ring(cav),'Frequency'));
len = findspos(ring,length(ring)+1);
harm = floor(freq/2.99792e8*len);

idpass = findcells(ring,'PassMethod','IdentityPass');
ring(idpass) = atsetfieldvalues(ring(idpass),'Length',0);

%ring w.o. radiations
ring = atsetcavity(ring,volt,0,harm);
ring = atradoff(ring,'IdentityPass','auto','auto');
save fcch_norad.mat ring

%ring with radiations
ring = atsetcavity(ring,volt,1,harm);
ring = atradon(ring,'CavityPass','auto','auto');
save fcch_rad.mat ring

%ring with radiations and tapering
co = findorbit6(ring,1:length(ring)+1);
cop_i = co(5,:);
bends = findcells(ring,'Class','Bend');
k0 = atgetfieldvalues(ring(bends),'BendingAngle')./atgetfieldvalues(ring(bends),'Length');
ring(bends) = atsetfieldvalues(ring(bends),'PolynomB',{1},k0.*cop_i(bends)');
save fcch_rad_tapered.mat ring