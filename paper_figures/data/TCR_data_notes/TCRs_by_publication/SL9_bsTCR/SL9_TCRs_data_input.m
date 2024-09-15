clear

%% Common across all epitopes
tcr_name="868Z11";
cdr3a="See Notes";
cdr3b="See Notes";
trav="See Notes";
traj="See Notes";
trbv="See Notes";
trbd="See Notes";
trbj="See Notes";
assay="NFAT-luc2 luminescence";
tcr_source_organism="Human";
index_peptide='SLYNTVATL';
mhc="DS-A*02:01";
pmid="31324691";
peptide_source_organism="HIV";
peptide_type="viral";
scan_type="Complete 1 Hamming";

%% Read peptide activity data

f=readtable('aav0860_table_s5.xlsx','sheet','Figure 6b');
V=table2array(f(:,5));%read TCR activation values at 100pM TCR concentration

index_peptide_activity=V(1);
%% write epitope from list

for i=1:size(f,1)
        disp(i)    
epitope=char(table2array(f(i,1)));

a=table(cellstr(tcr_name),cellstr(cdr3a),cellstr(cdr3b),cellstr(trav),cellstr(traj),...
cellstr(trbv),cellstr(trbd),cellstr(trbj),cellstr(assay),...
cellstr(tcr_source_organism),cellstr(index_peptide),index_peptide_activity,...
cellstr(mhc),cellstr(pmid),cellstr(peptide_source_organism),cellstr(peptide_type),...
cellstr(epitope),V(i),V(i)/index_peptide_activity,cellstr(scan_type));

writetable(a,'epitope_data_SL9.xlsx','WriteMode','Append','WriteVariableNames',false,'WriteRowNames',false);

end


