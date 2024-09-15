clear

TCR_list=["TCR3","18A2","NYE-S1","868Z11","TCR6","TCR7","A6","A23","TCR2-T"];
HLA_list=["A*02:01:01:01","B*81:01:01:01","A*02:01:01:01","A*02:01:01:01","A*02:01:01:01","B*27:05:02:01",...
    "A*02:01:01:01","A*02:01:01:01","A*02:01:01:01"];
MHC_list=repmat("YFAMYGEKVAHTHVDTLYVRYHYYTWAVLAYTWY",1,size(TCR_list,2));%populate with A0201 first
MHC_list(2)="YYSEYRNIYAQTDESNLYLSYNYYSLAVLAYEWY";
MHC_list(6)="YHTEYREICAKTDEDTLYLNYHDYTWAVLAYEWY";

%% Read data
epitope_database=readtable('TCR_epitope_database.xlsx');

for tcr=1:size(TCR_list,2)
  
  epitope=string(epitope_database.peptide(epitope_database.tcr_name==TCR_list(tcr)));
  CDR3b=repmat(string(unique(epitope_database.cdr3b(epitope_database.tcr_name==TCR_list(tcr)))),size(epitope,1),1);
  HLA=repmat(HLA_list(tcr),size(epitope,1),1);
  MHC=repmat(MHC_list(tcr),size(epitope,1),1);
  
  binder=randi([0,1],size(epitope));%mock binder data
  
  writetable(table(CDR3b,epitope,HLA,binder,MHC),'input.csv','WriteMode','append','WriteRowNames',false)
end
