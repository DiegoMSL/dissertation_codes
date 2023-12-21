function dy = rossler(t,y,param)

dy = [
-y(2)-y(3);
y(1) + param(1)*y(2);
param(2) + y(3)*(y(1) - param(3));
];
