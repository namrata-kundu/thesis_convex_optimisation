% syms a1 a2 a3 a4 b1 b2 b3 b4 c1 c2 c3 c4 x
% yy = int((-x-5-(a1*x*x+b1*x+c1))^2, [-999 -5]);
% disp(yy)
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
% % a1=4;
% % b1=2;
% % c1=3;
% subs(yy)
% 
% % disp(yy)



yalmip('clear')
 

%%%%%%%start here
syms a1 a2 a3 a4 b1 b2 b3 b4 c1 c2 c3 c4 x
int1 = int((-x-5-(a1*x*x+b1*x+c1))^2, [-999 -5])
int2 = int((x+5-(a2*x*x+b2*x+c2))^2, [-5 0])
int3 = int((-x+5-(a3*x*x+b3*x+c3))^2, [0 5])
int4 = int((x-5-(a4*x*x+b4*x+c4))^2, [5 999])

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

integral1 = eval(char(int1))
integral2 = eval(char(int2))
integral3 = eval(char(int3))
integral4 = eval(char(int4))



obj = integral1 + integral2 + integral3 + integral4;


Constraints = [];
for i=1:4
    Constraints = [Constraints, rho(2,i)>=0];
end

 

disp(Constraints)

 
Constraints = [Constraints, rho(2,1)*(rho(1,2)*rho(1,2)) + rho(3,1)*rho(1,2) + rho(4,1) == rho(2,2)*(rho(1,2)*rho(1,2)) + rho(3,2)*rho(1,2) + rho(4,2)];
Constraints = [Constraints, rho(2,2)*(rho(1,3)*rho(1,3)) + rho(3,2)*rho(1,3) + rho(4,2) == rho(2,3)*(rho(1,3)*rho(1,3)) + rho(3,3)*rho(1,3) + rho(4,3)];
Constraints = [Constraints, rho(2,3)*(rho(1,4)*rho(1,4)) + rho(3,3)*rho(1,4) + rho(4,3) == rho(2,4)*(rho(1,4)*rho(1,4)) + rho(3,4)*rho(1,4) + rho(4,4)];

disp(Constraints)


Constraints = [Constraints,2*rho(2,1)*rho(1,2) + rho(3,1) <= 2*rho(2,2)*rho(1,2) + rho(3,2)];
Constraints = [Constraints,2*rho(2,2)*rho(1,3) + rho(3,2) <= 2*rho(2,3)*rho(1,3) + rho(3,3)];
Constraints = [Constraints,2*rho(2,3)*rho(1,4) + rho(3,3) <= 2*rho(2,4)*rho(1,4) + rho(3,4)];
 
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