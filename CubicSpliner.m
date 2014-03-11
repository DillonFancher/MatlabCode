function   [yvals, xvals]= CubicSpliner(x, y)
%Routine to find the Coefficients for the cubic spline
%--------------------------------------------------------------------------
xpoints = x;
ypoints = y;

n = length(xpoints);
h = zeros(size(xpoints)-1);

b = zeros(1,n-1);
c = zeros(1,n-1);
d = zeros(1,n-1);

for i = 1:n-1
    h(i) = xpoints(i+1)- xpoints(i);
end

for i = 2:n-1
    alpha(i) = (3/h(i))*(ypoints(i+1)-ypoints(i)) - (3/h(i-1))*(ypoints(i)-ypoints(i-1));
end

lnot = zeros(1,n-1);
munot = zeros(1,n-1);
znot = zeros(1,n-1);

lnot(1) = 1;
munot(1) = 0;
znot(1) = 0;

for i = 2:n-1
    lnot(i) = 2*(xpoints(i+1)-xpoints(i-1))-h(i-1)*munot(i-1);
    munot(i) = h(i)/lnot(i);
    znot(i) = (alpha(i)-h(i-1)*znot(i-1))/lnot(i);
end

lnot(n) = 1;
munot(n) = 0;
c(n) = 0;

for i = 1:n-1
    c(n-i) = znot(n-i)-munot(n-i)*c(n-i+1);
    b(n-i) = (ypoints(n-i+1)-ypoints(n-i))/h(n-i) - h(n-i)*(c(n-i+1)+2*c(n-i))/3;
    d(n-i) = (c(n-i+1)-c(n-i))/(3*h(n-i));
end

a = ypoints;

for i = 1:n-1
disp('For spline');
disp(i-1);
disp('a =');
disp(a(i));
disp('b =');
disp(b(i));
disp('c =');
disp(c(i));
disp('d = ');
disp(d(i));
end

%end of the routine for the coefficients
%---------------------------------------------------------------------------

%routine to evaluate the cubic spline
%--------------------------------------------------------------------------
E = 100;
xvals = zeros(n-1, E);
for i = 1:n-1
    xvals(i, :) = linspace(xpoints(i), xpoints(i+1), E);
end

Spline = zeros(n-1, 4);
for i = 1:n-1
    Spline(i, 1) = a(i);
    Spline(i, 2) = b(i);
    Spline(i, 3) = c(i);
    Spline(i, 4) = d(i);
end

SplineEvals = zeros(n-1, E);
for i = 1:n-1
    for j = 1:E
        SplineEvals(i, j) = Spline(i, 1) + Spline(i, 2)*(xvals(i,j)-xpoints(i))...
            + Spline(i, 3)*(xvals(i,j) - xpoints(i))^2 + Spline(i, 4)*(xvals(i,j) - xpoints(i))^3;
    end
end

xvals = xvals;
yvals = SplineEvals;

%Commented out so LastNameLetter can clal CubicSpliner without interferance
%for i = 1:n-1
%    plot(xvals(i,:), SplineEvals(i,:))
%    hold on 
%end

%cubic spline evaluated
%---------------------------------------------------------------------------