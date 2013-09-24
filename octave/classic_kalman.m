disp ("Classic kalman")
disp ("--------------")

######
# f1 #
######

function ret = f1(x, u, w)
  ret = x(1)*cos(u(3)+w(3)) + x(3)*(u(4)+w(4))*cos(u(3)+w(3)) + x(2)*sin(u(3)+w(3)) - (u(1)+w(1)) + (u(5)+w(5));
endfunction

function ret = df1_dxx(x, u, w)
  ret = cos(u(3)+w(3));
endfunction

function ret = df1_dxy(x, u, w)
  ret = sin(u(3)+w(3));
endfunction

function ret = df1_dxvx(x, u, w)
  ret = (u(4)+w(4))*cos(u(3)+w(3));
endfunction

function ret = df1_dxvy(x, u, w)
  ret = 0;
endfunction

######
# f2 #
######

function ret = f2(x, u, w)
  ret = x(2)*cos(u(3)+w(3)) + x(4)*(u(4)+w(4))*cos(u(3)+w(3)) - x(1)*sin(u(3)+w(3)) - (u(2)+w(2)) + (u(5)+w(5));
endfunction

function ret = df2_dxx(x, u, w)
  ret = -sin(u(3)+w(3));
endfunction

function ret = df2_dxy(x, u, w)
  ret = cos(u(3)+w(3));
endfunction

function ret = df2_dxvx(x, u, w)
  ret = 0;
endfunction

function ret = df2_dxvy(x, u, w)
  ret = (u(4)+w(4)) * cos(u(3)+w(3));
endfunction

######
# f3 #
######

function ret = f3(x, u, w)
  ret = x(3)*cos(u(3)+w(3))*(u(6)+w(6)) + x(4)*sin(u(3)+w(3))*(u(6)+w(6));
endfunction

function ret = df3_dxx(x, u, w)
  ret = 0;
endfunction

function ret = df3_dxy(x, u, w)
  ret = 0;
endfunction

function ret = df3_dxvx(x, u, w)
  ret = cos(u(3)+w(3)) * (u(6)+w(6));
endfunction

function ret = df3_dxvy(x, u, w)
  ret = sin(u(3)+w(3)) * (u(6)+w(6));
endfunction

######
# f4 #
######

function ret = f4(x, u, w)
  ret = x(4)*cos(u(3)+w(3))*(u(6)+w(6)) - x(3)*sin(u(3)+w(3))*(u(6)+w(6));
endfunction

function ret = df4_dxx(x, u, w)
  ret = 0;
endfunction

function ret = df4_dxy(x, u, w)
  ret = 0;
endfunction

function ret = df4_dxvx(x, u, w)
  ret = -sin(u(3)+w(3)) * (u(6)+w(6));
endfunction

function ret = df4_dxvy(x, u, w)
  ret = cos(u(3)+w(3)) * (u(6)+w(6));
endfunction

#################################################
# Newtons difference quotient. Used for testing #
#################################################

function ret = lim_xx(f, x, u, w, h)
  ret = (feval(f, x+[h;0;0;0], u, w) - feval(f, x, u, w)) / h;
endfunction

function ret = lim_xy(f, x, u, w, h)
  ret = (feval(f, x+[0;h;0;0], u, w) - feval(f, x, u, w)) / h;
endfunction

function ret = lim_xvx(f, x, u, w, h)
  ret = (feval(f, x+[0;0;h;0], u, w) - feval(f, x, u, w)) / h;
endfunction

function ret = lim_xvy(f, x, u, w, h)
  ret = (feval(f, x+[0;0;0;h], u, w) - feval(f, x, u, w)) / h;
endfunction


#############################################
# To run tests from octave terminal execute #
#   > test classic_kalman.m                 #
#############################################

%!shared x, u, w, h
%! x = [800; 1200; 567; 100];
%! u = [832; 1358; 0.12; 1/30; 0.1; 0.9];
%! w = [100; 234; 234; 94; 245; 456];
%! h = 0.0001;

%!test assert ( df1_dxx(x, u, w), lim_xx(@f1, x, u, w, h), h );
%!test assert ( df1_dxy(x, u, w), lim_xy(@f1, x, u, w, h), h );
%!test assert ( df1_dxvx(x, u, w), lim_xvx(@f1, x, u, w, h), h );
%!test assert ( df1_dxvy(x, u, w), lim_xvy(@f1, x, u, w, h), h );

%!test assert ( df2_dxx(x, u, w), lim_xx(@f2, x, u, w, h), h );
%!test assert ( df2_dxy(x, u, w), lim_xy(@f2, x, u, w, h), h );
%!test assert ( df2_dxvx(x, u, w), lim_xvx(@f2, x, u, w, h), h );
%!test assert ( df2_dxvy(x, u, w), lim_xvy(@f2, x, u, w, h), h );

%!test assert ( df3_dxx(x, u, w), lim_xx(@f3, x, u, w, h), h );
%!test assert ( df3_dxy(x, u, w), lim_xy(@f3, x, u, w, h), h );
%!test assert ( df3_dxvx(x, u, w), lim_xvx(@f3, x, u, w, h), h );
%!test assert ( df3_dxvy(x, u, w), lim_xvy(@f3, x, u, w, h), h );

%!test assert ( df4_dxx(x, u, w), lim_xx(@f4, x, u, w, h), h );
%!test assert ( df4_dxy(x, u, w), lim_xy(@f4, x, u, w, h), h );
%!test assert ( df4_dxvx(x, u, w), lim_xvx(@f4, x, u, w, h), h );
%!test assert ( df4_dxvy(x, u, w), lim_xvy(@f4, x, u, w, h), h );
