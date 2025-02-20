clear

%% Common across all epitopes
tcr_name="1G4";
cdr3a="See Notes";
cdr3b="See Notes";
trav="See Notes";
traj="See Notes";
trbv="See Notes";
trbd="See Notes";
trbj="See Notes";
assay="ELISA";
tcr_source_organism="Human";
index_peptide='SLLMWITQC';
mhc="HLA-A*02:01";
pmid="37607971";
peptide_source_organism="Human";
peptide_type="cancer";
scan_type="Complete 1 Hamming";

%% Read peptide activity data

f=readtable('processed_data.xlsx');
f=f(172:end,:);
V=table2array(f(:,2));%read IFN values

index_peptide_activity=100;
%% write epitope from list

for i=1:size(f,1)
        disp(i)    
epitope=char(table2array(f(i,1)));

a=table(cellstr(tcr_name),cellstr(cdr3a),cellstr(cdr3b),cellstr(trav),cellstr(traj),...
cellstr(trbv),cellstr(trbd),cellstr(trbj),cellstr(assay),...
cellstr(tcr_source_organism),cellstr(index_peptide),index_peptide_activity,...
cellstr(mhc),cellstr(pmid),cellstr(peptide_source_organism),cellstr(peptide_type),...
cellstr(epitope),V(i),V(i)/index_peptide_activity,cellstr(scan_type));

writetable(a,'epitope_data_A23_1G4.xlsx','WriteMode','Append','WriteVariableNames',false,'WriteRowNames',false);

end


