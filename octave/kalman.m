  ####################################################################################
#                                                                                  #
# The MIT License (MIT)                                                            #
#                                                                                  #
# Copyright (c) 2013 abellagonzalo                                                 #
#                                                                                  #
# Permission is hereby granted, free of charge, to any person obtaining a copy of  #
# this software and associated documentation files (the "Software"), to deal in    #
# the Software without restriction, including without limitation the rights to     #
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of #
# the Software, and to permit persons to whom the Software is furnished to do so,  #
# subject to the following conditions:                                             #
#                                                                                  #
# The above copyright notice and this permission notice shall be included in all   #
# copies or substantial portions of the Software.                                  #
#                                                                                  #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR       #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS #
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR   #
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER   #
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN          #
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.       #  
#                                                                                  #
####################################################################################

source classic_kalman.m

# Initial state
x = [100; 0; 0; 0];

P = [100,   0,    0,    0;
       0, 100,    0,    0;
       0,   0, 1000,    0;
       0,   0,    0, 1000];

for i = 1:100

  #disp("Iter "), disp(i)

  # Odometry. Robot does not move.
  u = [0; 0; 0; 30; 0; 0.9];
  w = [0.01; 0.01; 0.01; 0; 10; 0];

  # Predict
  x = f(x, u, [0,0,0,0,0,0]);
  P = A(x, u, w) * P * A(x, u, w).' + W(x, u, w) * Q(w) * W(x, u, w).';

  # Observation
  z = [100 + 10*i; 0];
  v = [0.1; 0.1];

  # Update
  gain = P * H(x,v).' * inv( H(x,v) * P * H(x,v).' + V(x,v) * R(v) * V(x,v).' );
  x = x + gain * ( z - h(x, [0,0,0,0]) );
  P = ( eye(4) - gain * H(x,v) ) * P;

endfor

x
P
