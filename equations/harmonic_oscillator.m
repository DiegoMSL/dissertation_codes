function dy = harmonic_oscillator(t,y,param)

dy = [
param(4);
y(3);
- param(1)*y(3) - param(2)*y(2) + param(3)*cos(y(1)); 
];
