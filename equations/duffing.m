function dy = duffing(t,y,param)

dy = [
 y(2);
 -param(1)*y(2) - param(2)*y(1) - param(3)*y(1)^3 + param(4)*cos(y(3));
 param(5);
];