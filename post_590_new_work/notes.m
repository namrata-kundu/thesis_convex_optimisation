% In order to create dynamic variables:
for i=1:3
    eval(['A' num2str(i) '= i'])
end

% substitute syms with sdpvar:
syms a1 b1 c1 x
int1 = int((-x-5-(a1*x*x+b1*x+c1))^2, [-999 -5])

a1 = sdpvar(1,1);
b1 = sdpvar(1,1);
c1 = sdpvar(1,1);

integral1 = eval(char(int1))
