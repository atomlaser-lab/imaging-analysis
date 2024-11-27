function r = double_resonance(x)

A1 = 0.8;
w1 = 8.9;
x1 = 4.2;
A2 = 3.8;
w2 = 6.9;
x2 = 0.9;

x = x + 1.0784; %This is the detuning at which the maximum response occurs.
r = A1./(1 + 4*(x - x1).^2./w1.^2) + A2./(1 + 4*(x - x2).^2./w2.^2);
r = r/4.326;    %This normalises the response to 1 at the detuning of maximum response.