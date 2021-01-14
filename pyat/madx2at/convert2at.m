clear
%% FCC lattice
mod = py.importlib.import_module('prepare_fcc_seq');
fname = 'Optics_FCCee_h_217_nosol_3.seq';
mod.reformat_seq(fname);
mod.save_seq(fname);
seqfilemadX='Optics_FCCee_h_217_nosol_3_for_at.seq';
seqname = 'L000013';
E0=120e9;
%t: ttbar, Ebeam = 175 GeV or 182.5 GeV h: ZH, Ebeam = 120 GeV w: W+W-, Ebeam = 80 GeV z: Z, Ebeam = 45.6 GeV
% execute conversion. file ..._AT_LATTICE.mat will be created.
atfrommadx(seqfilemadX,E0);
load([seqfilemadX(1:end-4) '_AT_LATTICE.mat']); 
ring = L000013;
save fcch.mat ring