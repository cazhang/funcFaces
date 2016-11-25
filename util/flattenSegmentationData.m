% Generate a bunch of txt files with segmentations from the .mat file
function flattenSegmentationData(segDataFile, outdir)

load(segDataFile);
%      1 -  20  Human
%     21 -  40  Cup
%     41 -  60  Glasses
%     61 -  80  Airplane
%     81 - 100  Ant
%    101 - 120  Chair
%    121 - 140  Octopus
%    141 - 160  Table
%    161 - 180  Teddy
%    181 - 200  Hand
%    201 - 220  Plier
%    221 - 240  Fish
%    241 - 260  Bird
%   (261 - 280) Spring (excluded from our study)
%    281 - 300  Armadillo
%    301 - 320  Bust
%    321 - 340  Mech
%    341 - 360  Bearing
%    361 - 380  Vase
%    381 - 400  Fourleg 
groups = {'Human', 'Cup', 'Glass', 'Airplane', 'Ant', 'Chair', 'Octopus', ...
    'Desk', 'Teddy', 'Hand', 'Plier', 'Fish', 'Bird', 'Spring', 'Armadillor', 'Bust', 'Mech', ...
    'Bearing', 'Utensil', 'FourLegs'};
for groupid = 1:length(groups)
    % no springs
    if groupid == 14
        continue;
    end
    for meshid = 1:20
        fid = (groupid-1)*20 + meshid;
        cellname = ['C_' sprintf('%02d',groupid) '_' groups{groupid} '_Smoothed'];
        seg = eval([cellname '{' num2str(meshid) '}']);
        fname = [outdir '/' num2str(fid) '_joint_seg.txt'];
        if size(seg,2) > 1, seg = seg'; end;
        dlmwrite(fname,seg,'newline','pc');
    end
end
