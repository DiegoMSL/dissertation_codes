function dy = simple_pendulum(t, theta, param)

dy = [
theta(2);
-(param(1)/param(2))*sin(theta(1));
];

end