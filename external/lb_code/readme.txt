To use this directory, you'll need to compile MeshLP code. Generally, just
"cd MeshLP"
"make"
should work. You may have to modify the Matlab paths in the "makefile" to point to the right files. Sometimes Matlab complains about missing libraries. In that case, in Ubuntu, the following works:

launch Maltab using the following command line:

LD_PRELOAD=/usr/lib/libstdc++.so.6 /opt/matlab2010a/bin/matlab
