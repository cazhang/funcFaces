function S = seg_colormap(n)
% Colormap from: Brewer, C.A., G.W. Hatchard, and M.A. Harrower, (2003), 
% Colorbrewer in print: A catalog of color schemes for maps, 
% Cartography and Geographic Information Science, 30, 5-32
% (http://www.ColorBrewer.org).
S1 = [
    228 26 28
    55 126 184
    77 175 74
    152 78 163
    255 127 0
    255 255 51
    166 86 40
    247 129 191   
    153 153 153
    ];
S1 = S1/256;

S2 = [
    141 211 199
    255 255 179
    190 186 218
    251 128 114
    128 177 211
    253 180 98
    179 222 105
    252 205 229
    217 217 217 
    188 128 189
    204 235 197
    255 237 111
    ];
S2 = S2/256;

if n <= size(S1,1)
    S = S1(1:n,:);
elseif n <= size(S2,1);
    S = S2(1:n,:);
else
    nc = size(S2,1);
    jetCol = hsv(n - nc + 1);
    S = [S2; jetCol(1:end-1,:)];
end