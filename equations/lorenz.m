function dy = lorenz(t,y,param)

dy = [
param(1)*(y(2)-y(1));
y(1)*(param(3)-y(3))-y(2);
y(1)*y(2)-param(2)*y(3);
];