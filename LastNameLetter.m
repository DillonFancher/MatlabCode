function LastNameLetter(x, y)

  % [yvals, xvals] = CubicSpliner(x,y);
   xPos = x;
   yPos = y;
   n = length(xPos);
   xt = [1:n];
   yt = [1:n];
   
        [xvals, xtime] = CubicSpliner(xt, xPos);
        [yvals, ytime] = CubicSpliner(yt, yPos);
   
   for i = 1:n-1
        plot(xvals(i,:), yvals(i,:))
        hold on
   end