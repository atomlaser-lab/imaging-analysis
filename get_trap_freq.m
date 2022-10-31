function [f,D,trap] = get_trap_freq(P1,P2)
%
% Use to calculate optical trap frequencies
%
trap = optical_trap('Rb87',[gaussian_beam(1064e-9,P1,100e-6,0),gaussian_beam(1090e-9,P2,150e-6,0)]);
trap.lasers(1).set_rotations('x',90);trap.lasers(2).set_rotations('x',90,'z',30);
trap.ext_force.Fx = @(~,~,~) 0;
trap.ext_force.Fy = @(~,~,~) 0;
trap.ext_force.Fz = @(~,~,~) const.mRb*const.g;

f = trap.freq(0,0,0,1e-6);
D = trap.depth(0,0,0,'z')/const.kb*1e6;
