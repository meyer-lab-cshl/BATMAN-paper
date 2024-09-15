clear

epitope_data=readtable('TCR_epitope_database.xlsx');

mhcs=string(unique(epitope_data.mhc));

for imhc=1:size(mhcs,1)
    epitopes=unique(epitope_data(epitope_data.mhc==mhcs(imhc),17));
    filename=erase(['epitopes_',char(mhcs(imhc)),'.xlsx'],["*",":"]);
    writetable(epitopes,filename)
end