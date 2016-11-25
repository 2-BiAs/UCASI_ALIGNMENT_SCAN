//Generate uCasi Mirror Scan Boundary

alpha = 54.7*%pi/180;
theta = [linspace(%pi/2 - alpha, %pi/2 + alpha, 8)' ; linspace(3*%pi/2 - alpha, 3*%pi/2 + alpha, 8)'];
R = 26;

mPoints = [R * cos(theta), R * sin(theta)];
plot(mPoints(:,1), mPoints(:,2));
