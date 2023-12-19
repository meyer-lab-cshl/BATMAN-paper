clear

TCR_list=["TCR3","18A2","NYE-S1","868Z11","TCR6","TCR7","A6","A23","TCR2-T"];
HLA_list=["A*02:01","B*81:01","A*02:01","A*02:01","A*02:01","B*27:05","A*02:01","A*02:01","A*02:01"];

%% Read data
epitope_database=readtable('TCR_epitope_database.xlsx');

for tcr=1:size(TCR_list,2)
  
  Antigen=string(epitope_database.peptide(epitope_database.tcr_name==TCR_list(tcr)));
  CDR3=repmat(string(unique(epitope_database.cdr3b(epitope_database.tcr_name==TCR_list(tcr)))),size(Antigen,1),1);
  HLA=repmat(HLA_list(tcr),size(Antigen,1),1);
  disp([TCR_list(tcr),size(Antigen,1)])
  writetable(table(CDR3,Antigen,HLA),'pmtnet_input.csv','WriteMode','append','WriteRowNames',false)
end