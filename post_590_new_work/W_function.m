yalmip('clear')
 

%%%%%%%start here
% syms x 
% % x = sdpvar(1,1);
% y = piecewise(x<1,(1/2)*x*x+1,(1<=x)<2.5,2*x-(1/2), (2.5<=x)<6, -x + 7, x>=6, x-5);
% fplot(y)



f = [-inf(), -5, 0, 5;
    0, 0, 0, 0;
    -1, 1, -1, 1;
    -5, 5, 5, -5];

%  syms a1 a2 a3 a4 b1 b2 b3 b4 c1 c2 c3 c4
a1 = sdpvar(1,1);
a2 = sdpvar(1,1);
a3 = sdpvar(1,1);
a4 = sdpvar(1,1);
b1 = sdpvar(1,1);
b2 = sdpvar(1,1);
b3 = sdpvar(1,1);
b4 = sdpvar(1,1);
c1 = sdpvar(1,1);
c2 = sdpvar(1,1);
c3 = sdpvar(1,1);
c4 = sdpvar(1,1);

rho = [-inf(), -5, 0, 5;
    a1, a2, a3, a4;
    b1, b2, b3, b4;
    c1, c2, c3, c4];

% f = [-inf(), 1, 2.5, 6;
%     (1/2), 0, 0, 0;
%     0, 2, -1, 1;
%     1, -(1/2), 7, -5];
% 
% % syms a1 a2 a3 a4 b1 b2 b3 b4 c1 c2 c3 c4
% a1 = sdpvar(1,1);
% a2 = sdpvar(1,1);
% a3 = sdpvar(1,1);
% a4 = sdpvar(1,1);
% b1 = sdpvar(1,1);
% b2 = sdpvar(1,1);
% b3 = sdpvar(1,1);
% b4 = sdpvar(1,1);
% c1 = sdpvar(1,1);
% c2 = sdpvar(1,1);
% c3 = sdpvar(1,1);
% c4 = sdpvar(1,1);
% 
% rho = [-inf(), 1, 2.5, 6;
%     a1, a2, a3, a4;
%     b1, b2, b3, b4;
%     c1, c2, c3, c4];
% 
% % integral1= (2*a1 - 1)^2/20 - 2*c1 - (2*a1)/3 + b1*(c1 - 1) + b1^2/3 + c1^2 + (b1*(30*a1 - 15))/60 + (c1*(2*a1 - 1))/3 + (Inf*(2*a1 - 1)^2 - Inf*b1*(c1 - 1) - Inf*b1*(30*a1 - 15) + Inf*(20*b1^2 - 40*a1 + 20*c1*(2*a1 - 1) + 20) - 120*c1 + 60*c1^2 + 60)*Inf + 4/3;
% integral1=49750499500250*(2*a1 - 1)^2 - 2000*c1 - (1994006000*a1)/3 - 998000*b1*(c1 - 1) + (997003000*b1^2)/3 + 1000*c1^2 - (49800299800*b1*(30*a1 - 15))/3 + (997003000*c1*(2*a1 - 1))/3 + 997006000/3; %%taking -if to -999
% % integral1=(9999500009999900064399361*(2*a1 - 1)^2)/20 - 200000*c1 - (1999940000600000*a1)/3 - 9999800000*b1*(c1 - 1) + (999970000300000*b1^2)/3 + 100000*c1^2 - (4999800002999980000*b1*(30*a1 - 15))/3 + (999970000300000*c1*(2*a1 - 1))/3 + 999970000600000/3; %%taking -if to -99999
% % integral1 = (99950009999000040001*(2*a1 - 1)^2)/20 - 20000*c1 - (1999400060000*a1)/3 - 99980000*b1*(c1 - 1) + (999700030000*b1^2)/3 + 10000*c1^2 - (499800029998000*b1*(30*a1 - 15))/3 + (999700030000*c1*(2*a1 - 1))/3 + 999700060000/3; %%taking -if to -9999
% integral2 =(39*(b2 - 2)^2)/8 + (3*(2*c2 + 1)^2)/8 + (609*a2*(b2 - 2))/32 + (21*(2*b2 - 4)*(2*c2 + 1))/16 + (3093*a2^2)/160 + (39*a2*(2*c2 + 1))/8;
% integral3 =(1603*(b3 + 1)^2)/24 + (7*(c3 - 7)^2)/2 + (20111*a3*(b3 + 1))/32 + (1603*a3*(c3 - 7))/12 + (245707*a3^2)/160 + (119*(b3 + 1)*(c3 - 7))/4;
% % integral4 =Inf*(c4 + 5)^2 + Inf*((b4 - 1)^2/3 + (2*a4*(c4 + 5))/3) + Inf*a4*(b4 - 1) + Inf*a4^2 + Inf*(b4 - 1)*(c4 + 5) - 72*(b4 - 1)^2 - 6*(c4 + 5)^2 - 648*a4*(b4 - 1) - 144*a4*(c4 + 5) - (7776*a4^2)/5 - 36*(b4 - 1)*(c4 + 5);
% integral4 = 332334261*(b4 - 1)^2 + 993*(c4 + 5)^2 + (996005994705*a4*(b4 - 1))/2 + 664668522*a4*(c4 + 5) + (995009989997223*a4^2)/5 + 997965*(b4 - 1)*(c4 + 5); %%taking inf to 999
% % integral4 = 333323333433261*(b4 - 1)^2 + 99993*(c4 + 5)^2 + 49998000029999799672*a4*(b4 - 1) + 666646666866522*a4*(c4 + 5) + (9999500009999899966824864*a4^2)/5 + 9999799965*(b4 - 1)*(c4 + 5);  %%taking inf to 99999-huge values wont work in objective
% % integral4 = 333233343261*(b4 - 1)^2 + 9993*(c4 + 5)^2 + 4998000299979352*a4*(b4 - 1) + 666466686522*a4*(c4 + 5) + (99950009999000035744*a4^2)/5 + 99979965*(b4 - 1)*(c4 + 5);  %%taking inf to 9999


%%new
% integral1=199001998001000*(a1 - 1/2)^2 + 1000*(c1 - 1)^2 - 498002998000*b1*(a1 - 1/2) - 998000*b1*(c1 - 1) + (997003000*b1^2)/3 + (1994006000*(a1 - 1/2)*(c1 - 1))/3;
% integral2=(39*(b2 - 2)^2)/8 + (3*(c2 + 1/2)^2)/2 + (609*a2*(b2 - 2))/32 + (39*a2*(c2 + 1/2))/4 + (3093*a2^2)/160 + (21*(b2 - 2)*(c2 + 1/2))/4;
% integral3=(1603*(b3 + 1)^2)/24 + (7*(c3 - 7)^2)/2 + (20111*a3*(b3 + 1))/32 + (1603*a3*(c3 - 7))/12 + (245707*a3^2)/160 + (119*(b3 + 1)*(c3 - 7))/4;
% integral4=332334261*(b4 - 1)^2 + 993*(c4 + 5)^2 + (996005994705*a4*(b4 - 1))/2 + 664668522*a4*(c4 + 5) + (995009989997223*a4^2)/5 + 997965*(b4 - 1)*(c4 + 5);


% integral1=(997002874*(b1 + 1)^2)/3 + 994*(c1 + 5)^2 - 498002997688*a1*(b1 + 1) + (1994005748*a1*(c1 + 5))/3 + (995009990001874*a1^2)/5 - 997976*(b1 + 1)*(c1 + 5);
% integral1=(970174*(b1 + 1)^2)/3 + 94*(c1 + 5)^2 - 48029488*a1*(b1 + 1) + (1940348*a1*(c1 + 5))/3 + (9509897374*a1^2)/5 - 9776*(b1 + 1)*(c1 + 5);
% integral1=Inf*(c1 + 5)^2 + Inf*((b1 + 1)^2/3 + (2*a1*(c1 + 5))/3) + Inf*a1^2 - Inf*a1*(b1 + 1) - Inf*(b1 + 1)*(c1 + 5) - (125*(b1 + 1)^2)/3 - 5*(c1 + 5)^2 + (625*a1*(b1 + 1))/2 - (250*a1*(c1 + 5))/3 - 625*a1^2 + 25*(b1 + 1)*(c1 + 5);
integral1=(26875*(b1 + 1)^2)/3 + 25*(c1 + 5)^2 - (809375*a1*(b1 + 1))/2 + (53750*a1*(c1 + 5))/3 + 4859375*a1^2 - 875*(b1 + 1)*(c1 + 5);
integral2=(125*(b2 - 1)^2)/3 + 5*(c2 - 5)^2 - (625*a2*(b2 - 1))/2 + (250*a2*(c2 - 5))/3 + 625*a2^2 - 25*(b2 - 1)*(c2 - 5);
integral3=(125*(b3 + 1)^2)/3 + 5*(c3 - 5)^2 + (625*a3*(b3 + 1))/2 + (250*a3*(c3 - 5))/3 + 625*a3^2 + 25*(b3 + 1)*(c3 - 5);
% integral4=(997002874*(b4 - 1)^2)/3 + 994*(c4 + 5)^2 + 498002997688*a4*(b4 - 1) + (1994005748*a4*(c4 + 5))/3 + (995009990001874*a4^2)/5 + 997976*(b4 - 1)*(c4 + 5);
% integral4=(970174*(b4 - 1)^2)/3 + 94*(c4 + 5)^2 + 48029488*a4*(b4 - 1) + (1940348*a4*(c4 + 5))/3 + (9509897374*a4^2)/5 + 9776*(b4 - 1)*(c4 + 5);
% integral4=piecewise(a4 == 0 & imag(b4)^2 == (real(b4) - 1)^2, - (125*(b4 - 1)^2)/3 - 5*(c4 + 5)^2 - 25*(b4 - 1)*(c4 + 5), a4 == 0 & imag(b4)^2 < (real(b4) - 1)^2, Inf - (125*(b4 - 1)^2)/3 - 5*(c4 + 5)^2 - 25*(b4 - 1)*(c4 + 5), a4 == 0 & (real(b4) - 1)^2 < imag(b4)^2, - Inf - (125*(b4 - 1)^2)/3 - 5*(c4 + 5)^2 - 25*(b4 - 1)*(c4 + 5), imag(a4)^2 ~= real(a4)^2 & (real(a4) == 0 | in(a4, 'real')), - Inf*sign(imag(a4)^2/5 - real(a4)^2/5) - (125*(b4 - 1)^2)/3 - 5*(c4 + 5)^2 - (625*a4*(b4 - 1))/2 - (250*a4*(c4 + 5))/3 - 625*a4^2 - 25*(b4 - 1)*(c4 + 5), 0 < imag(a4)*real(a4) & imag(a4)^2 ~= real(a4)^2, Inf*1i - Inf*sign(imag(a4)^2/5 - real(a4)^2/5), imag(a4)*real(a4) < 0 & imag(a4)^2 ~= real(a4)^2, - Inf*1i - Inf*sign(imag(a4)^2/5 - real(a4)^2/5), 0 < imag(a4)*real(a4) & real(a4)*(real(b4) - 1) ~= imag(a4)*imag(b4) & imag(a4)^2 == real(a4)^2, Inf*1i + Inf*sign(real(a4)*(real(b4)/2 - 1/2) - (imag(a4)*imag(b4))/2), imag(a4)*real(a4) < 0 & real(a4)*(real(b4) - 1) ~= imag(a4)*imag(b4) & imag(a4)^2 == real(a4)^2, Inf*sign(real(a4)*(real(b4)/2 - 1/2) - (imag(a4)*imag(b4))/2) - Inf*1i, 0 < imag(a4)*real(a4) & real(a4)*(real(b4) - 1) == imag(a4)*imag(b4) & imag(a4)^2 == real(a4)^2 & imag(b4)^2 + 2*imag(a4)*imag(c4) == (real(b4) - 1)^2 + real(a4)*(2*real(c4) + 10), Inf*1i - (125*(b4 - 1)^2)/3 - 5*(c4 + 5)^2 - (625*a4*(b4 - 1))/2 - (250*a4*(c4 + 5))/3 - 625*a4^2 - 25*(b4 - 1)*(c4 + 5), imag(a4)*real(a4) < 0 & real(a4)*(real(b4) - 1) == imag(a4)*imag(b4) & imag(a4)^2 == real(a4)^2 & imag(b4)^2 + 2*imag(a4)*imag(c4) == (real(b4) - 1)^2 + real(a4)*(2*real(c4) + 10), - Inf*1i - (125*(b4 - 1)^2)/3 - 5*(c4 + 5)^2 - (625*a4*(b4 - 1))/2 - (250*a4*(c4 + 5))/3 - 625*a4^2 - 25*(b4 - 1)*(c4 + 5), 0 < imag(a4)*real(a4) & real(a4)*(real(b4) - 1) == imag(a4)*imag(b4) & imag(a4)^2 == real(a4)^2 & 2*real(b4) + imag(b4)^2 + 2*imag(a4)*imag(c4) < 10*real(a4) + real(b4)^2 + 2*real(a4)*real(c4) + 1, Inf + Inf*1i, 0 < imag(a4)*real(a4) & real(a4)*(real(b4) - 1) == imag(a4)*imag(b4) & imag(a4)^2 == real(a4)^2 & 10*real(a4) + real(b4)^2 + 2*real(a4)*real(c4) + 1 < 2*real(b4) + imag(b4)^2 + 2*imag(a4)*imag(c4), - Inf + Inf*1i, imag(a4)*real(a4) < 0 & real(a4)*(real(b4) - 1) == imag(a4)*imag(b4) & imag(a4)^2 == real(a4)^2 & 2*real(b4) + imag(b4)^2 + 2*imag(a4)*imag(c4) < 10*real(a4) + real(b4)^2 + 2*real(a4)*real(c4) + 1, Inf - Inf*1i, imag(a4)*real(a4) < 0 & real(a4)*(real(b4) - 1) == imag(a4)*imag(b4) & imag(a4)^2 == real(a4)^2 & 10*real(a4) + real(b4)^2 + 2*real(a4)*real(c4) + 1 < 2*real(b4) + imag(b4)^2 + 2*imag(a4)*imag(c4), - Inf - Inf*1i);
integral4=(26875*(b4 - 1)^2)/3 + 25*(c4 + 5)^2 + (809375*a4*(b4 - 1))/2 + (53750*a4*(c4 + 5))/3 + 4859375*a4^2 + 875*(b4 - 1)*(c4 + 5);
% obj = integral_1 + integral_2 + integral_3 + integral_4;
obj = integral1 + integral2 + integral3 + integral4;
 

%yy = int((-x-5-(a1*x*x+b1*x+c1))^2, [-999 -5])
 integral1 = (997002874*(b1 + 1)^2)/3 + 994*(c1 + 5)^2 - 498002997688*a1*(b1 + 1) + (1994005748*a1*(c1 + 5))/3 + (995009990001874*a1^2)/5 - 997976*(b1 + 1)*(c1 + 5);
% integral1 = (999700029874*(b1 + 1)^2)/3 + 9994*(c1 + 5)^2 - 4998000299979688*a1*(b1 + 1) + (1999400059748*a1*(c1 + 5))/3 + (99950009999000046874*a1^2)/5 - 99979976*(b1 + 1)*(c1 + 5);
% integral1 = 498988*b1 - (997002874*a1)/3 - 994*c1 + 494018;
% integral1=(997002874*a1)/3 - 498988*b1 + 994*c1 - 494018;

%yy = int((x+5-(a2*x*x+b2*x+c2))^2, [-5 0])
integral2=(125*(b2 - 1)^2)/3 + 5*(c2 - 5)^2 - (625*a2*(b2 - 1))/2 + (250*a2*(c2 - 5))/3 + 625*a2^2 - 25*(b2 - 1)*(c2 - 5);
% integral2 = (25*b2)/2 - (125*a2)/3 - 5*c2 + 25/2;
% integral2=(125*a2)/3 - (25*b2)/2 + 5*c2 - 25/2;


%yy = int((-x+5-(a3*x*x+b3*x+c3))^2, [0 5])
integral3=(125*(b3 + 1)^2)/3 + 5*(c3 - 5)^2 + (625*a3*(b3 + 1))/2 + (250*a3*(c3 - 5))/3 + 625*a3^2 + 25*(b3 + 1)*(c3 - 5);
% integral3 = 25/2 - (25*b3)/2 - 5*c3 - (125*a3)/3;
% integral3=(125*a3)/3 + (25*b3)/2 + 5*c3 - 25/2;

%yy = int((x-5-(a4*x*x+b4*x+c4))^2, [5 999])
integral4 = (997002874*(b4 - 1)^2)/3 + 994*(c4 + 5)^2 + 498002997688*a4*(b4 - 1) + (1994005748*a4*(c4 + 5))/3 + (995009990001874*a4^2)/5 + 997976*(b4 - 1)*(c4 + 5);
% integral4=(999700029874*(b4 - 1)^2)/3 + 9994*(c4 + 5)^2 + 4998000299979688*a4*(b4 - 1) + (1999400059748*a4*(c4 + 5))/3 + (99950009999000046874*a4^2)/5 + 99979976*(b4 - 1)*(c4 + 5);
% integral4 = 494018 - 498988*b4 - 994*c4 - (997002874*a4)/3;
% integral4=(997002874*a4)/3 + 498988*b4 + 994*c4 - 494018;

obj = integral1 + integral2 + integral3 + integral4;


Constraints = [];
for i=1:4
    Constraints = [Constraints, rho(2,i)>=0];
end

% 
% disp(Constraints)
% 
% Constraints = [Constraints, function1(rho(1,2)) == function2(rho(1,2))];
% Constraints = [Constraints, function2(rho(1,3)) == function3(rho(1,3))];
% Constraints = [Constraints, function3(rho(1,4)) == function4(rho(1,4))];

disp(Constraints)

% x=sdpvar(1,1); 
%rho_i(x_i+1) == rho_i+1(x_i+1)
Constraints = [Constraints, rho(2,1)*(rho(1,2)*rho(1,2)) + rho(3,1)*rho(1,2) + rho(4,1) == rho(2,2)*(rho(1,2)*rho(1,2)) + rho(3,2)*rho(1,2) + rho(4,2)];
Constraints = [Constraints, rho(2,2)*(rho(1,3)*rho(1,3)) + rho(3,2)*rho(1,3) + rho(4,2) == rho(2,3)*(rho(1,3)*rho(1,3)) + rho(3,3)*rho(1,3) + rho(4,3)];
Constraints = [Constraints, rho(2,3)*(rho(1,4)*rho(1,4)) + rho(3,3)*rho(1,4) + rho(4,3) == rho(2,4)*(rho(1,4)*rho(1,4)) + rho(3,4)*rho(1,4) + rho(4,4)];

disp(Constraints)


Constraints = [Constraints,2*rho(2,1)*rho(1,2) + rho(3,1) <= 2*rho(2,2)*rho(1,2) + rho(3,2)];
Constraints = [Constraints,2*rho(2,2)*rho(1,3) + rho(3,2) <= 2*rho(2,3)*rho(1,3) + rho(3,3)];
Constraints = [Constraints,2*rho(2,3)*rho(1,4) + rho(3,3) <= 2*rho(2,4)*rho(1,4) + rho(3,4)];
% Constraints = [Constraints,(integral1 + integral2 + integral3 + integral4)>=0];

% Constraints = [Constraints,(5-5==(a1*(-5)*(-5)+b1*(-5)+c1))];

Constraints = [Constraints,(999-5==(a1*(-999)*(-999)+b1*(-999)+c1))];
% Constraints = [Constraints,(555-5==(a1*(-555)*(-555)+b1*(-555)+c1))];
% Constraints = [Constraints,(5-5==(a4*5*5+b4*5+c4))];

Constraints = [Constraints,(999-5==(a4*999*999+b4*999+c4))];
% Constraints = [Constraints,(555-5==(a4*555*555+b4*555+c4))];
% Constraints = [Constraints,(998-5==(a4*998*998+b4*998+c4))];

disp(Constraints)

 
%   options = sdpsettings('solver','', 'debug', 1, 'verbose', 1);
  options = sdpsettings('solver','bmibnb', 'verbose', 1);
objec = obj;
% Constraints = [Constraints, (5>=obj)>=0]
sol = optimize(Constraints,objec, options)
% sol = optimize(Constraints,obj)
disp("a1="+value(a1))
disp("a2="+value(a2))
disp("a3="+value(a3))
disp("a4="+value(a4))

disp("b1="+value(b1))
disp("b2="+value(b2))
disp("b3="+value(b3))
disp("b4="+value(b4))

disp("c1="+value(c1))
disp("c2="+value(c2))
disp("c3="+value(c3))
disp("c4="+value(c4))

disp(value(obj))
% disp(cvx_optval)

%%%displaying constraints with calculated values
for i=1:4
    if value(rho(2,i))>=0
        disp("True");
    else
        disp("False");
    end
%     Constraints = [Constraints, rho(2,i)>=0];
end

% Constraints = [Constraints, rho(2,1)*(rho(1,2)*rho(1,2)) + rho(3,1)*rho(1,2) + rho(4,1) == rho(2,2)*(rho(1,2)*rho(1,2)) + rho(3,2)*rho(1,2) + rho(4,2)];
% Constraints = [Constraints, rho(2,2)*(rho(1,3)*rho(1,3)) + rho(3,2)*rho(1,3) + rho(4,2) == rho(2,3)*(rho(1,3)*rho(1,3)) + rho(3,3)*rho(1,3) + rho(4,3)];
% Constraints = [Constraints, rho(2,3)*(rho(1,4)*rho(1,4)) + rho(3,3)*rho(1,4) + rho(4,3) == rho(2,4)*(rho(1,4)*rho(1,4)) + rho(3,4)*rho(1,4) + rho(4,4)];
if value(rho(2,1))*(rho(1,2)*rho(1,2)) + value(rho(3,1))*rho(1,2) + value(rho(4,1)) == value(rho(2,2))*(rho(1,2)*rho(1,2)) + value(rho(3,2))*rho(1,2) + value(rho(4,2))
 disp("True");
else
    disp(value(rho(2,1))*(rho(1,2)*rho(1,2)) + value(rho(3,1))*rho(1,2) + value(rho(4,1)))
    disp(value(rho(2,2))*(rho(1,2)*rho(1,2)) + value(rho(3,2))*rho(1,2) + value(rho(4,2)))
    disp("False");
end

if value(rho(2,2))*(rho(1,3)*rho(1,3)) + value(rho(3,2))*rho(1,3) + value(rho(4,2)) == value(rho(2,3))*(rho(1,3)*rho(1,3)) + value(rho(3,3))*rho(1,3) + value(rho(4,3))
 disp("True");
else
    disp(value(rho(2,2))*(rho(1,3)*rho(1,3)) + value(rho(3,2))*rho(1,3) + value(rho(4,2)))
    disp(value(rho(2,3))*(rho(1,3)*rho(1,3)) + value(rho(3,3))*rho(1,3) + value(rho(4,3)))
    disp("False");
end

if value(rho(2,3))*(rho(1,4)*rho(1,4)) + value(rho(3,3))*rho(1,4) + value(rho(4,3)) == value(rho(2,4))*(rho(1,4)*rho(1,4)) + value(rho(3,4))*rho(1,4) + value(rho(4,4))
 disp("True");
else
    disp(value(rho(2,3))*(rho(1,4)*rho(1,4)) + value(rho(3,3))*rho(1,4) + value(rho(4,3)))
    disp(value(rho(2,4))*(rho(1,4)*rho(1,4)) + value(rho(3,4))*rho(1,4) + value(rho(4,4)))
    disp("False");
end
%%%%%%%%%%%%

% if value(a1)*(rho(1,2)*rho(1,2)) + value(b1)*rho(1,2) + value(c1) == value(a2)*(rho(1,2)*rho(1,2)) + value(b2)*rho(1,2) + value(c2)
%  disp("True");
% else
%     disp(value(a1)*(rho(1,2)*rho(1,2)) + value(b1)*rho(1,2) + value(c1))
%     disp(value(a2)*(rho(1,2)*rho(1,2)) + value(b2)*rho(1,2) + value(c2))
%     disp("False");
% end
% 
% if value(rho(2,1))*(rho(1,2)*rho(1,2)) + value(rho(3,1))*rho(1,2) + value(rho(4,1)) == value(rho(2,2))*(rho(1,2)*rho(1,2)) + value(rho(3,2))*rho(1,2) + value(rho(4,2))
%  disp("True");
% else
%     disp(value(rho(2,1))*(rho(1,2)*rho(1,2)) + value(rho(3,1))*rho(1,2) + value(rho(4,1)))
%     disp(value(rho(2,2))*(rho(1,2)*rho(1,2)) + value(rho(3,2))*rho(1,2) + value(rho(4,2)))
%     disp("False");
% end
% 
% if value(rho(2,3))*(rho(1,4)*rho(1,4)) + value(rho(3,3))*rho(1,4) + value(rho(4,3)) == value(rho(2,4))*(rho(1,4)*rho(1,4)) + value(rho(3,4))*rho(1,4) + value(rho(4,4))
%  disp("True");
% else
%     disp(value(rho(2,3))*(rho(1,4)*rho(1,4)) + value(rho(3,3))*rho(1,4) + value(rho(4,3)))
%     disp(value(rho(2,4))*(rho(1,4)*rho(1,4)) + value(rho(3,4))*rho(1,4) + value(rho(4,4)))
%     disp("False");
% end

% Constraints = [Constraints,2*rho(2,1)*rho(1,2) + rho(3,1) <= 2*rho(2,2)*rho(1,2) + rho(3,2)];
% Constraints = [Constraints,2*rho(2,2)*rho(1,3) + rho(3,2) <= 2*rho(2,3)*rho(1,3) + rho(3,3)];
% Constraints = [Constraints,2*rho(2,3)*rho(1,4) + rho(3,3) <= 2*rho(2,4)*rho(1,4) + rho(3,4)];

if 2*value(rho(2,1))*rho(1,2) + value(rho(3,1)) <= 2*value(rho(2,2))*rho(1,2) + value(rho(3,2))
    disp("True")
else
    disp("False")
end

if 2*value(rho(2,2))*rho(1,3) + value(rho(3,2)) <= 2*value(rho(2,3))*rho(1,3) + value(rho(3,3))
    disp("True")
else
    disp("False")
end


if 2*value(rho(2,3))*rho(1,4) + value(rho(3,3)) <= 2*value(rho(2,4))*rho(1,4) + value(rho(3,4))
    disp("True")
else
    disp("False")
end

%%%%%%%%%%%%%%%
% disp(sol)

function result = function1(x)
    result = (1/2)*x*x+1;
end


function result = function2(x)
    result = 2*x-(1/2);
end
 
 
 
function result = function3(x)
    result = -x + 7;
end
 
function result = function4(x)
    result = x-5;
end