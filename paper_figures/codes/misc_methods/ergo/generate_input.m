clear

TCR_list=["TCR3","18A2","NYE-S1","868Z11","TCR6","TCR7","A6","A23","TCR2-T"];

clear

TCR_list=["TCR3","18A2","NYE-S1","868Z11","TCR6","TCR7","A6","A23","TCR2-T"];
HLA_list=["A-02:01","B-81:01","A-02:01","A-02:01","A-02:01","B-27:05","A-02:01","A-02:01","A-02:01"];

%% Read data
epitope_database=readtable('TCR_epitope_database.xlsx');

for tcr=1:size(TCR_list,2)
  
  Peptide=string(epitope_database.peptide(epitope_database.tcr_name==TCR_list(tcr)));
  TRA=repmat(string(unique(epitope_database.cdr3a(epitope_database.tcr_name==TCR_list(tcr)))),size(Peptide,1),1);
  TRB=repmat(string(unique(epitope_database.cdr3b(epitope_database.tcr_name==TCR_list(tcr)))),size(Peptide,1),1);
  TRAV=repmat(append("TRAV",string(unique(epitope_database.trav(epitope_database.tcr_name==TCR_list(tcr))))),size(Peptide,1),1);
  TRAJ=repmat(append("TRAJ",string(unique(epitope_database.traj(epitope_database.tcr_name==TCR_list(tcr))))),size(Peptide,1),1);
  TRBV=repmat(append("TRBV",string(unique(epitope_database.trbv(epitope_database.tcr_name==TCR_list(tcr))))),size(Peptide,1),1);
  TRBJ=repmat(append("TRBJ",string(unique(epitope_database.trbj(epitope_database.tcr_name==TCR_list(tcr))))),size(Peptide,1),1);
  T_Cell_Type=repmat("CD8",size(Peptide,1),1);
  MHC=repmat(HLA_list(tcr),size(Peptide,1),1);
  writetable(table(TRA,TRB,TRAV,TRAJ,TRBV,TRBJ,T_Cell_Type,Peptide,MHC),'ergo_input.csv','WriteMode','append','WriteRowNames',false)
end

%% Remember to remove "See Paper" cells in the csv if present

