clear all;
m = load('boundary_view.txt');
x = m(:,1);
y = m(:,2);
z = m(:,3);

plot3(x,y,z,'.')
grid on
xlabel('x')
ylabel('y')
zlabel('z')
