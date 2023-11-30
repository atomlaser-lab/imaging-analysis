function [f,D,trap] = get_trap_freq(P1,P2,w1,w2)
%
% Use to calculate optical trap frequencies
%

if nargin < 3
    w1 = 106e-6;
    w2 = 153e-6;
elseif nargin < 4
    w1 = 106e-6;
end
% trap = optical_trap('Rb87',[gaussian_beam(1064e-9,P1,112e-6,0),gaussian_beam(1090e-9,P2,150e-6,0)]);
trap = optical_trap('Rb87',[gaussian_beam(1064e-9,P1,w1,0),gaussian_beam(1090e-9,P2,w2,0)]);
trap.lasers(1).set_rotations('x',90);trap.lasers(2).set_rotations('x',90,'z',30);
trap.ext_force.Fx = @(~,~,~) 0;
trap.ext_force.Fy = @(~,~,~) 0;
trap.ext_force.Fz = @(~,~,~) const.mRb*const.g;

f = trap.freq(0,0,0,1e-6);
D = trap.depth(0,0,0,'z')/const.kb*1e6;

f = f(:)';