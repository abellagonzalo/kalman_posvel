
function retval = g_1(ix, iy, itheta, it, x, y, vx, vy)

  retval = [  x*cos(itheta) + vx*cos(itheta)*it + y*sin(itheta) - ix;
             -x*sin(itheta) - vy*sin(itheta)*it + y*cos(itheta) - iy;
              vx;
              vy];

endfunction

function retval = g(u_, mu_)

  retval = g_1(  u_(1),  u_(2),  u_(3),  u_(4),
                mu_(1), mu_(2), mu_(3), mu_(4));

endfunction


function retval = G_2(cosO, sinO, t)

  retval = [cosO   sinO  cosO*t        0;
            cosO  -sinO       0  -sinO*t;
               0      0       1        0;
               0      0       0        1];

endfunction


function retval = G_1(theta, t)

  retval = G_2(cos(theta), sin(theta), t);

endfunction

function retval = G(u)

  retval = G_1( u(3), u(4) );

endfunction

function retval = h(mu)
  
  retval = [ mu(1); mu(2)];

endfunction

function retval = H()

  retval = [ 1, 0, 0, 0;
             0, 1, 0, 0];

endfunction

function retval = predict_state(g)

  retval = g;

endfunction

function retval = predict_sigma(G, sigma, R)

  retval = G * sigma * G.' + R; #'

endfunction

function retval = K(H_, sigma, Q)

  retval = sigma * H_.' * (H_ * sigma * H_.' + Q);

endfunction

function retval = update_state(mu, K, z, h_)

  retval = mu + K * (z - h_);

endfunction

function retval = update_sigma(K, H, sigma)

  retval = (eye(4) - K * H) * sigma;

endfunction

mu = [0; 0; 0; 0];
sigma = eye(4) * 1000;
u = [0; 0; 0; 1];


for i = 1:5

  # Predict
  mu = predict_state( g(u, mu) );
  sigma = predict_sigma( G(u), sigma, eye(4));

  z = [i+10; 0;];

  #Update
  K_ = sigma * H().' * inv( H() * sigma * H().' + 0 );
  mu = mu + K_ * (z - h(mu))
  sigma = update_sigma(K_, H(), sigma)

endfor

mu;
sigma;


source classic_kalman.m
