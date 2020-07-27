


%generate initial tracks roughly in the SLF region using the WM fod template
unix(['tckgen wmfod_template.mif -seed sphere 0,18,57,2 -initdir 0,0,1 fornix.tck']);

unix(['tckgen wmfod_template.mif -seed sphere 0,18,57,2 -initdir 0,0,1 fornix.tck -include -28,-7,42,10']);