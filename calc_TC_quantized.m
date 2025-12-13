function [TC_ch, TC_WL] = calc_TC_quantized(r, p, E_c0, E_v0, Field, cell_info)

[idx, radius_quantized, error] = getRadiusIndex(cell_info(1), r);
idx_WL = r.n_num;

% idx-- : TC_ch calculation, for TC_ch = TC_ch*TC_value(idx--)

TC_ch = 1;
TC_WL = 1;

for i=idx:-1:2
    TC_ch = TC_ch*TC_value(r.r(i), r.r(i-1), E_c0(i)-(E_c0(idx)-cell_info(2)), E_c0(i-1)-(E_c0(idx)-cell_info(2)), Field(i), 0.5, p);
end

for i=idx:idx_WL
    TC_WL = TC_WL*TC_value(r.r(i), r.r(i+1), E_c0(i)-(E_c0(idx)-cell_info(2)), E_c0(i+1)-(E_c0(idx)-cell_info(2)), Field(i), 0.5, p);
end

end