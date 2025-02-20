clear

%% Common across all epitopes
tcr_name="T3";
cdr3a="See paper";
cdr3b="See paper";
trav="See paper";
traj="See paper";
trbv="See paper";
trbd="See paper";
trbj="See paper";
assay="ELISA";
tcr_source_organism="Human";
index_peptide='ALYDKTKRIFL';
mhc="HLA-A*02:01";
pmid="34873326";
peptide_source_organism="Human";
peptide_type="cancer (TdT)";
scan_type="Complete 1 Hamming";

%% Read peptide activity data
aaorder='ACDEFGHIKLMNPQRSTVWY';

V=readmatrix('41587_2021_1089_MOESM4_ESM.xlsx','sheet','Fig. 2g','NumHeaderLines',0);
V=V(14:24,2:21);%read IFN values

[~,index_aa]=ismember(index_peptide,aaorder);
index_peptide_activity=mean(V(sub2ind([size(index_peptide,2),20],1:size(index_peptide,2),index_aa)),'all');%average over replicates

%% Start by writing the index peptide
a=table(cellstr(tcr_name),cellstr(cdr3a),cellstr(cdr3b),cellstr(trav),cellstr(traj),...
cellstr(trbv),cellstr(trbd),cellstr(trbj),cellstr(assay),...
cellstr(tcr_source_organism),cellstr(index_peptide),index_peptide_activity,...
cellstr(mhc),cellstr(pmid),cellstr(peptide_source_organism),cellstr(peptide_type),...
cellstr(index_peptide),index_peptide_activity,1,cellstr(scan_type));

writetable(a,'epitope_data_T1_T3.xlsx','WriteMode','Append','WriteVariableNames',false,'WriteRowNames',false);

for i=1:size(index_peptide,2)
for j=1:size(aaorder,2)
        disp([i,j])
epitope=index_peptide;
if epitope(i)~=aaorder(j) %skip index peptide
    
epitope(i)=aaorder(j);
value=V(i,j);
a=table(cellstr(tcr_name),cellstr(cdr3a),cellstr(cdr3b),cellstr(trav),cellstr(traj),...
cellstr(trbv),cellstr(trbd),cellstr(trbj),cellstr(assay),...
cellstr(tcr_source_organism),cellstr(index_peptide),index_peptide_activity,...
cellstr(mhc),cellstr(pmid),cellstr(peptide_source_organism),cellstr(peptide_type),...
cellstr(epitope),V(i,j),V(i,j)/index_peptide_activity,cellstr(scan_type));

writetable(a,'epitope_data_T1_T3.xlsx','WriteMode','Append','WriteVariableNames',false,'WriteRowNames',false);
end

end
end


