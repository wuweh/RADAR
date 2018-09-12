

function z = MstEq(x)


z = [0.8*cos(x(1,:))-0.2*cos(x(1,:)+x(2,:)); 0.8*sin(x(1,:))-0.2*sin(x(1,:)+x(2,:))];
