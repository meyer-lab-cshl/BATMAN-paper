clear

%% Common across all epitopes
tcr_name="B7";
cdr3a="See Notes";
cdr3b="See Notes";
trav="See Notes";
traj="See Notes";
trbv="See Notes";
trbd="See Notes";
trbj="See Notes";
assay="minigene depletion";
tcr_source_organism="Human";
index_peptide='LLFGYPVYV';
mhc="HLA-A*02:01";
pmid="32184297";
peptide_source_organism="HTLV(Tax)";
peptide_type="viral";
scan_type="Complete 1 Hamming";

%% Read peptide activity data

f=readtable('B7_position_depletion_matrix.csv');
V=table2array(f(1:171,3));%read depletion values
V(V>1.5)=1.5;%clip values (B7 only)

pos=table2array(f(:,1));%mutation position
aa=char(table2array(f(:,2)));%mutated aa

%normalize
V(:)=normalize(-log2(V(:)),"range");

%identify index peptide
index_peptide_activity=1.0001;

%write index peptide
a=table(cellstr(tcr_name),cellstr(cdr3a),cellstr(cdr3b),cellstr(trav),cellstr(traj),...
cellstr(trbv),cellstr(trbd),cellstr(trbj),cellstr(assay),...
cellstr(tcr_source_organism),cellstr(index_peptide),index_peptide_activity,...
cellstr(mhc),cellstr(pmid),cellstr(peptide_source_organism),cellstr(peptide_type),...
cellstr(index_peptide),1.0001,1,cellstr(scan_type));

writetable(a,'epitope_data_A6_B7.xlsx','WriteMode','Append','WriteVariableNames',false,'WriteRowNames',false);


%% write epitope from list (take index peptide only once)

for i=1:171
        disp(i) 

epitope=index_peptide; %start with index peptide
epitope(pos(i))=aa(i);%put the mutation


a=table(cellstr(tcr_name),cellstr(cdr3a),cellstr(cdr3b),cellstr(trav),cellstr(traj),...
cellstr(trbv),cellstr(trbd),cellstr(trbj),cellstr(assay),...
cellstr(tcr_source_organism),cellstr(index_peptide),index_peptide_activity,...
cellstr(mhc),cellstr(pmid),cellstr(peptide_source_organism),cellstr(peptide_type),...
cellstr(epitope),V(i),V(i),cellstr(scan_type));

writetable(a,'epitope_data_A6_B7.xlsx','WriteMode','Append','WriteVariableNames',false,'WriteRowNames',false);
end


