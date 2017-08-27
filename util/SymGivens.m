function [c,s,d] = SymGivens(a,b)
	% Symmetric Givens rotation from M. A. Saunders and S.-C. Choi.
	% http://www.stanford.edu/group/SOL/dissertations/sou-cheng-choi-thesis.pdf
  	if b == 0
    	if a == 0
      		c = 1;
    	else
      		c = sign(a);  % NOTE: sign(0) = 0 in MATLAB
    	end
    	s = 0;
    	d = abs(a);

  	elseif a == 0
    	c = 0;
    	s = sign(b);
    	d = abs(b);

  	elseif abs(b) > abs(a)
   		t = a/b;
    	s = sign(b) / sqrt(1 + t^2);
    	c = s*t;
    	d = b/s; % computationally better than d = a / c since |c| <= |s|

  	else
    	t = b/a;
    	c = sign(a) / sqrt(1 + t^2);
    	s = c*t;
    	d = a/c; % computationally better than d = b / s since |s| <= |c|
  	end