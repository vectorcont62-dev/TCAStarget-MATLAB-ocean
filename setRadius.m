function r = setRadius(r_rf, t_si, t_tox1, t_tn1, t_tox2, t_ctn, t_box, t_top_channel, CTN_slice)
r = struct();
if (r_rf == 0) r_rf = 0.01e-7; end
if (t_si == 0) t_si = 0.01e-7; end
if (t_tox1 == 0) t_tox1 = 0.01e-7; end
if (t_tn1 == 0) t_tn1 = 0.01e-7; end
if (t_tox2 == 0) t_tox2 = 0.01e-7; end
if (t_ctn == 0) t_ctn = 0.01e-7; end
if (t_box == 0) t_box = 0.01e-7; end
if (t_top_channel == 0) t_top_channel = 0.01e-7; end

r_0 = r_rf + t_si;

r.r_rf = r_rf;
r.t_si = t_si;
r.r_si = r_rf + t_si;

r.r_0 = r_0;

r.t_tox1 = t_tox1;
r.t_o1 = t_tox1;
r.r_tox1 = r_0 + t_tox1;
r.r_o1 = r_0 + t_tox1;

r.t_tn1 = t_tn1;
r.t_n1 = t_tn1;
r.r_n1 = r.r_o1 + t_tn1;
r.r_tn1 = r.r_tox1 + t_tn1;

r.t_o2 = t_tox2;
r.r_tox2 = r.r_tn1 + t_tox2;
r.r_tox = r.r_tox2;          % top = box, blo
r.r_o2 = r.r_tox2;

r.t_ctn = t_ctn;
r.r_ctn = r.r_tox + t_ctn;
r.t_box = t_box;
r.r_box = r.r_ctn + t_box;
r.t_top_channel = t_top_channel;
r.r_top_channel = r.r_box + t_top_channel;
R = r.r_top_channel;
n_num = 1e3;
r.n_num = n_num;
r.dr = (R-r_0)/n_num; dr = r.dr;
r.r = ((r_0-1e-7):dr:R)';
r.l_WL = 22e-7;

r.r_CTN=zeros(CTN_slice,1);
for i=1:length(r.r_CTN)
    r.r_CTN(i,1) = r.r_tox2 + i/CTN_slice*r.t_ctn;
end

end