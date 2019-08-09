% write_all_tables

ci = [3.0862    0.0009];
fn = 'TADPOLE_Submission_Tohka_Ciszek_method2.xlsx';
subtab = writeTADPOLEtables(srids,squeeze(D2pred(:,1,:)),squeeze(D2pred(:,2,:)),permute(D2prob_scores,[1 3 2]),ci,fn,'ID1');
subtab2 = writeTADPOLEtables(D3_rids,squeeze(D3pred(:,1,:)),squeeze(D3pred(:,2,:)),permute(D3prob_scores,[1 3 2]),ci,fn,'ID5');