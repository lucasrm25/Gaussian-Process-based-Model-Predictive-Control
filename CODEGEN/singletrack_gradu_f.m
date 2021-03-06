function gradu = singletrack_gradu_f(in1,in2,in3)
%SINGLETRACK_GRADU_F
%    GRADU = SINGLETRACK_GRADU_F(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    09-Jan-2020 11:52:34

I_z = in3(2,:);
M = in3(1,:);
T = in2(2,:);
V_vx = in1(4,:);
V_vy = in1(5,:);
c_f = in3(8,:);
delta = in2(1,:);
deltamax = in3(5,:);
l_f = in3(3,:);
maxbrakeWForce = in3(6,:);
maxmotorWForce = in3(7,:);
psi_dot = in1(6,:);
t2 = 1.0./I_z;
t3 = 1.0./M;
t4 = T.*5.0e+1;
t5 = delta.*5.0e+1;
t6 = deltamax.*5.0e+1;
t8 = V_vy.*1i;
t9 = V_vx.*1.0e+2;
t13 = l_f.*psi_dot.*1i;
t7 = -t4;
t10 = -t5;
t11 = -t6;
t12 = tanh(t9);
t24 = V_vx+t8+t13;
t14 = t7+5.0e+1;
t15 = t7-5.0e+1;
t18 = t6+t10;
t21 = t10+t11;
t26 = angle(t24);
t16 = exp(t14);
t17 = exp(t15);
t19 = exp(t18);
t22 = exp(t21);
t20 = t16+1.0;
t23 = t17+1.0;
t25 = t19+1.0;
t29 = t22+1.0;
t27 = 1.0./t20;
t30 = 1.0./t23;
t32 = 1.0./t25;
t34 = 1.0./t29;
t28 = t27.^2;
t31 = t30.^2;
t33 = t32.^2;
t35 = t34.^2;
t36 = t27.*5.0e+1;
t37 = t27-1.0;
t38 = t30.*5.0e+1;
t40 = deltamax.*t32;
t41 = t32-1.0;
t43 = t34-1.0;
t39 = -t36;
t42 = -t38;
t44 = -t40;
t45 = deltamax.*t43;
t47 = t16.*t28.*5.0e+1;
t49 = t17.*t31.*5.0e+1;
t50 = t16.*t28.*2.5e+3;
t53 = t17.*t31.*2.5e+3;
t55 = t6.*t19.*t33;
t56 = deltamax.*t19.*t33.*-5.0e+1;
t57 = t6.*t22.*t35;
t58 = t30.*t37;
t59 = deltamax.*t22.*t35.*-5.0e+1;
t61 = t37.*t38;
t64 = t34.*t41;
t67 = t4.*t16.*t28.*t30;
t69 = t4.*t17.*t31.*t37;
t70 = t5.*t19.*t33.*t34;
t72 = t5.*t22.*t35.*t41;
t46 = -t45;
t48 = -t47;
t51 = -t49;
t52 = -t50;
t54 = -t53;
t60 = T.*t58;
t63 = t4.*t58;
t65 = delta.*t64;
t68 = T.*t30.*t50;
t71 = T.*t37.*t53;
t90 = t56+t59+t64+t70+t72;
t62 = -t60;
t66 = -t65;
t74 = t39+t42+t63+5.0e+1;
t86 = t26+t44+t46+t65;
t87 = t48+t51+t58+t67+t69;
t89 = t52+t54+t61+t68+t71;
t73 = t30+t37+t62;
t75 = exp(t74);
t77 = t40+t45+t66;
t76 = t75+1.0;
t78 = sin(t77);
t81 = cos(t77);
t79 = 1.0./t76;
t80 = t79.^2;
t82 = maxmotorWForce.*t79;
t83 = t79-1.0;
t84 = maxbrakeWForce.*t12.*t83;
t91 = maxmotorWForce.*t75.*t80.*t89;
t92 = maxbrakeWForce.*t12.*t75.*t80.*t89;
t85 = -t84;
t93 = -t92;
t88 = t82+t85;
t94 = t91+t93;
gradu = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t3.*(c_f.*t78.*t90-c_f.*t81.*t86.*t90+(t73.*t78.*t88.*t90)./2.0),-t3.*((t73.*t94)./2.0+(t87.*t88)./2.0+(t73.*t81.*t94)./2.0+(t81.*t87.*t88)./2.0),0.0,-t3.*(c_f.*t81.*t90+c_f.*t78.*t86.*t90+(t73.*t81.*t88.*t90)./2.0),-t3.*((t73.*t78.*t94)./2.0+(t78.*t87.*t88)./2.0),0.0,-t2.*(c_f.*l_f.*t81.*t90+c_f.*l_f.*t78.*t86.*t90+(l_f.*t73.*t81.*t88.*t90)./2.0),-t2.*((l_f.*t73.*t78.*t94)./2.0+(l_f.*t78.*t87.*t88)./2.0),0.0,0.0,0.0,1.0],[3,7]);
