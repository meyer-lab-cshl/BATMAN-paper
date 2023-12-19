clear

TCR_list=["TCR3","18A2","NYE-S1","868Z11","TCR6","TCR7","A6","A23","TCR2-T"];
%% Read data
epitope_database=readtable('TCR_epitope_database.xlsx');

for tcr=1:size(TCR_list,2)
  
  Peptide=string(epitope_database.peptide(epitope_database.tcr_name==TCR_list(tcr)));
  CDR3b=repmat(string(unique(epitope_database.cdr3b(epitope_database.tcr_name==TCR_list(tcr)))),size(Peptide,1),1);
  
  writetable(table(CDR3b,Peptide),'input-sequences.csv','Delimiter',';','WriteMode','append','WriteRowNames',false)
end
