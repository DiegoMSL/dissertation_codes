function dy = van_der_pol(t,y,param)

dy = [
y(2);
param(1)*y(2) - param(1)*(y(1)^2)*y(2) - param(2)*y(1) + param(3)*cos(y(3));
param(4); 
];