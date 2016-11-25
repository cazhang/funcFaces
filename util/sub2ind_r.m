% sub2ind, row major
function loc = sub2ind_r(s,i,j)
loc = sub2ind([s(2),s(1)],j,i);