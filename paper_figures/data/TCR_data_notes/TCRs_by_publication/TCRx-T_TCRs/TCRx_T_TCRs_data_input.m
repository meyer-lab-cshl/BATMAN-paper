clear

%% Common across all epitopes
tcr_name="TCR3-T";
cdr3a="see paper";
cdr3b="see paper";
trav="see paper";
traj="see paper";
trbv="see paper";
trbd="see paper";
trbj="see paper";
assay="ELISA";
tcr_source_organism="Human";
index_peptide='FMNKFIYEI';
index_peptide_activity=1;
mhc="HLA-A*02:01";
pmid="32395117";
peptide_source_organism="Human";
peptide_type="cancer(AFP)";
scan_type="Complete 1 Hamming";

%% Read peptide activity data
aaorder='ACDEFGHIKLMNPQRSTVWY';

f=readmatrix('extracted_heights_of_bars.xlsx');%read full data

l_height=f(:,6:10:end);%extract lower bar height from mutational scan data for one TCR
u_height=f(:,11:10:end);%extract upper bar height from mutational scan data for one TCR

l_scale=f(1,2:10:end)./f(1,3:10:end);%extract scale for upper bar at 9 positions (unit: value/px height)
l_scale=repmat(l_scale,[20,1]);%same scale for all aa
u_scale=f(1,7:10:end)./f(1,8:10:end);%extract scale for upper bar at 9 positions (unit: value/px height)
u_scale=repmat(u_scale,[20,1]);%same scale for all aa

V=l_height.*l_scale+u_height.*u_scale;%add heights and scales of upper and lower bar to get peptide activity values to store
V=V';
%% Start by writing the index peptide
a=table(cellstr(tcr_name),cellstr(cdr3a),cellstr(cdr3b),cellstr(trav),cellstr(traj),...
cellstr(trbv),cellstr(trbd),cellstr(trbj),cellstr(assay),...
cellstr(tcr_source_organism),cellstr(index_peptide),index_peptide_activity,...
cellstr(mhc),cellstr(pmid),cellstr(peptide_source_organism),cellstr(peptide_type),...
cellstr(index_peptide),1,1,cellstr(scan_type));

writetable(a,'epitope_data_TCRx_T.xlsx','WriteMode','Append','WriteVariableNames',false,'WriteRowNames',false);

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
cellstr(epitope),V(i,j),V(i,j),cellstr(scan_type));

writetable(a,'epitope_data_TCRx_T.xlsx','WriteMode','Append','WriteVariableNames',false,'WriteRowNames',false);
end

end
end


