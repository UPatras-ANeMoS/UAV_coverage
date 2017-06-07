% Copyright 2016 Sotiris Papatheodorou
% 
% Licensed under the Apache License, Version 2.0 (the \"License\");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%    http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an \"AS IS\" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

function out1 = FJni(a,t,xi,xj,yi,yj)
%FJNI
%    OUT1 = FJNI(A,T,XI,XJ,YI,YJ)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    19-May-2016 02:04:54

t2 = yi-yj;
t3 = xi-xj;
t4 = 1.0./t3.^2;
t5 = t2.^2;
t6 = t4.*t5;
t7 = t6+1.0;
t8 = 1.0./sqrt(t7);
t9 = sinh(t);
t10 = a.*t8.*t9;
t11 = cosh(t);
t12 = 1.0./t3;
t13 = t3.^2;
t14 = t13.*(1.0./4.0);
t15 = t5.*(1.0./4.0);
t16 = a.^2;
t17 = t14+t15-t16;
t18 = sqrt(t17);
t19 = t2.*t8.*t11.*t12.*t18;
t20 = t10+t19;
t21 = abs(t20);
t23 = t8.*t11.*t18;
t24 = a.*t2.*t8.*t9.*t12;
t25 = t23-t24;
t22 = abs(t25);
t26 = t21.^2;
t27 = t22.^2;
t28 = t26+t27;
t29 = 1.0./sqrt(t28);
t30 = conj(t18);
t31 = a.*t8.*t11;
t33 = t2.*t8.*t9.*t12.*t30;
t34 = t31+t33;
t35 = t20.*t29.*t34;
t36 = t8.*t9.*t30;
t37 = a.*t2.*t8.*t11.*t12;
t38 = t36-t37;
t39 = t25.*t29.*t38;
t40 = t35+t39;
t42 = t20.*t29.*t40;
t43 = t2.*t8.*t9.*t12.*t18;
t44 = t31-t42+t43;
t32 = abs(t44);
t47 = t8.*t9.*t18;
t48 = t25.*t29.*t40;
t49 = t37-t47+t48;
t41 = abs(t49);
t45 = 1.0./t7.^(3.0./2.0);
t46 = t32.^2;
t50 = t41.^2;
t51 = t46+t50;
t52 = 1.0./sqrt(t51);
t53 = 1.0./t30;
t54 = xi.*(1.0./2.0);
t55 = t54-xj.*(1.0./2.0);
t56 = 1.0./t3.^3;
t57 = 1.0./t3.^4;
t58 = yi.*2.0;
t62 = yj.*2.0;
t59 = t58-t62;
t60 = yi.*(1.0./2.0);
t61 = t60-yj.*(1.0./2.0);
out1 = [-t44.*t52.*(a.*t5.*t11.*t45.*t56-t2.*t4.*t8.*t9.*t30+t2.*t8.*t9.*t12.*t53.*t55.*(1.0./2.0)+t2.*t5.*t9.*t30.*t45.*t57-1.0./2.0)+t49.*t52.*(t8.*t9.*t53.*t55.*(1.0./2.0)+a.*t2.*t4.*t8.*t11+t5.*t9.*t30.*t45.*t56-a.*t2.*t5.*t11.*t45.*t57);-t44.*t52.*(t8.*t9.*t12.*t30-a.*t4.*t11.*t45.*t59.*(1.0./2.0)+t2.*t8.*t9.*t12.*t53.*t61.*(1.0./2.0)-t2.*t9.*t30.*t45.*t56.*t59.*(1.0./2.0))+t49.*t52.*(-a.*t8.*t11.*t12+t8.*t9.*t53.*t61.*(1.0./2.0)-t4.*t9.*t30.*t45.*t59.*(1.0./2.0)+a.*t2.*t11.*t45.*t56.*t59.*(1.0./2.0)+1.0./2.0)];
