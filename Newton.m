function y = Newton(fun, dfun, x0, a, b)

%Created by Dillon Fancher, 1/27/13
%--------------------------------------------------------------------------
%Set format in order to see all of the digits to determine the accuracy
format long

%Read in the functions from the user input as objects
f = inline(fun);
df = inline(dfun);

%Calculate the first two iterations through brute force in order to keep
%good organization in the while loop
x(1) = x0 - (f(x0)/df(x0));
x(2) = x(1) - (f(x(1))/df(x(1)));

%Create a set of values representing our calculated x value corresponding
%to the root with a re-ordered incices for a pretty plot without missing
%values. (i.e. x(0) is not allowed)
xvals(1) = x0;
xvals(2) = x(1);
xvals(3) = x(2);

%Store values of i in a set in order to plot our calculated value of x with
%its corresponding iteration number
ivals(1) = 0;
ivals(2) = 1;
ivals(3) = 2;

%Declare a starting point for our while loop
i = 2;

%While loop to iterate over the bisection method until the difference
%between the current and previous x values is less than or equal to 10^-6
%in order to garuntee at least 6 digits of accuracy in the fractional part
%of our calculated value.
while(abs(x(i)-x(i-1)) >= 10^(-6))
    
    x(i+1) = x(i) - (f(x(i))/df(x(i)));
    i = i+1;
    
    xvals(i+1) = x(i);
    ivals(i+1) = i;
end

%Display our final answer and plot the absolute error at each iteration
display(x(i))
plot(ivals, (xvals-x(i)))
%--------------------------------------------------------------------------