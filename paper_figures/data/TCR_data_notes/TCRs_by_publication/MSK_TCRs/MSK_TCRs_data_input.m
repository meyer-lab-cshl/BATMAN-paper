clear

%% Common across all epitopes
tcr_name="TCR7";
cdr3a="See Notes";
cdr3b="See Notes";
trav="See Notes";
traj="See Notes";
trbv="See Notes";
trbd="See Notes";
trbj="See Notes";
assay="CD137 expression";
tcr_source_organism="Human";
index_peptide='GRLKALCQR';
mhc="HLA-B*27:05";
pmid="35589842";
peptide_source_organism="Human";
peptide_type="cancer";
scan_type="Complete 1 Hamming";

%% Read peptide activity data

f=readtable('41586_2022_4735_MOESM13_ESM.xlsx');
mutation=string(table2array(f(:,31)));%extract peptide mutations
V=table2array(f(:,32));%read fitted EC50 values
%clip values
V(V==-Inf|V<1E-4)=1E-4;
V(V==Inf|V>1E4)=1E4;
%take log and normalize
V=normalize(-log10(V),"range");

%identify index peptide
index_peptide_index=6;
index_peptide_activity=V(index_peptide_index);

%write the index peptide
a=table(cellstr(tcr_name),cellstr(cdr3a),cellstr(cdr3b),cellstr(trav),cellstr(traj),...
cellstr(trbv),cellstr(trbd),cellstr(trbj),cellstr(assay),...
cellstr(tcr_source_organism),cellstr(index_peptide),index_peptide_activity,...
cellstr(mhc),cellstr(pmid),cellstr(peptide_source_organism),cellstr(peptide_type),...
cellstr(index_peptide),V(index_peptide_index),V(index_peptide_index)/index_peptide_activity,cellstr(scan_type));
writetable(a,'epitope_data_MSK.xlsx','WriteMode','Append','WriteVariableNames',false,'WriteRowNames',false);


%% write epitope from list

for i=1:size(f,1)
        disp(i) 

epitope=index_peptide; %start with index peptide
p=char(mutation(i));

if p(1)~=p(3) %skip the index peptide
epitope(str2num(p(2)))=p(3);%put the mutation


a=table(cellstr(tcr_name),cellstr(cdr3a),cellstr(cdr3b),cellstr(trav),cellstr(traj),...
cellstr(trbv),cellstr(trbd),cellstr(trbj),cellstr(assay),...
cellstr(tcr_source_organism),cellstr(index_peptide),index_peptide_activity,...
cellstr(mhc),cellstr(pmid),cellstr(peptide_source_organism),cellstr(peptide_type),...
cellstr(epitope),V(i),V(i)/index_peptide_activity,cellstr(scan_type));

writetable(a,'epitope_data_MSK.xlsx','WriteMode','Append','WriteVariableNames',false,'WriteRowNames',false);
end

end


