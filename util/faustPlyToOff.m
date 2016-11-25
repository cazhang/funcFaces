% Ply to Off

clc;
addpath('/local/sda7/cz679/myToolbox/toolbox_graph');

pathToPly = '/local/sda7/cz679/myData/MPI-FAUST/training/registrations/';
plyList = dir( [pathToPly, '*.ply']);

for i = 11:length(plyList)
    i
    plyName = plyList(i).name;
    offName = [plyName(1:end-4), '.off'];
    plyFullName = [pathToPly, plyName];
    offFullName = [pathToPly, offName];
    [pts, tri] = read_ply( plyFullName );
    
    write_off( offFullName, pts, tri );
end