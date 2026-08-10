function Field=readCeresdata(addsp1,addsp2,name)

Field.Tc = readmatrix([addsp1,name],'Sheet','T_c');
Field.phi= readmatrix([addsp1,name],'Sheet','phi_c');
Field.v_ol= readmatrix([addsp1,name],'Sheet','v_ol_c');
Field.C= readmatrix([addsp1,name],'Sheet','C_c');
H_c= readmatrix([addsp1,name],'Sheet','H_c');
Field.H=H_c';

Field.t2= readmatrix([addsp2,name],'Sheet','t');
Field.V_phi= readmatrix([addsp2,name],'Sheet','phi');
Field.V_phiol= readmatrix([addsp2,name],'Sheet','olphi');

t2_1m=Field.t2-4.568e9;
[~, Field.b]=min(abs(t2_1m));

Field.MOI=readmatrix([addsp2,name],'Sheet','IC');
