nthred = 24;
nassign = 5500;
nparamset = 130918;

parfor ithred = 1:nthred
    output = FunctionalConnectome_recalcForTau_forParallel(ithred,nassign);
end

output = zeros(nparamset,64);
pyfinal= zeros(nparamset,49);
for ithred = 1:nthred
    nstart = (ithred-1)*nassign+1;
    nend   = min(ithred*nassign, nparamset);
    filename=strcat("./Matrics_addMeasures",num2str(ithred),".csv");
    output(nstart:nend,:) = readmatrix(filename);
    filename=strcat("./Updated_yfinals",num2str(ithred),".csv");
    pyfinal(nstart:nend,:)=readmatrix(filename);
end
npassed = sum(output(:,3));
disp(['npassed = ', num2str(npassed), '          outof ',num2str(nparamset),' paramsets']);

writematrix(pyfinal(nstart:nend,:),"./Updated_yfinals_all.csv");

output_head = {'Num' 'ifconverge' 'ifpass' 'pCa' 'kfca' 'IbarNCX' 'ks' 'koCa' 'ec50SR' 'Vmax_SRCaP' ...
               'Kmf' 'Kmr'            'SRCamax'   'SRCaRel'   'SRCaVmx'   'SRCaTau'   'Camin'   'Camax'   'CaTau' 'amp' ...
               'SRCaRfit' 'CaTfitR' 'M-SRCamax' 'M-SRCaRel' 'M-SRCaVmx' 'M-SRCaTau' 'M-Camin' 'M-Camax' 'M-CaTau' '#passed' ...
               'APD90' 'APD80' 'APD50' 'APD25' 'APmax' 'APmin'  ...
               'sumICa' 'sumIncx' 'sumJsrrel' 'sumJserca' 'sumIki' 'sumItof' 'sumItos' 'sumIks' 'sumIkr' 'sumIna' 'sumIcabk' 'sumIpca' 'sumSRleak' ...
               'peakICa' 'peakIncx' 'peakJsrrel' 'peakJserca' 'peakIki' 'peakItof' 'peakItos' 'peakIks' 'peakIkr' 'peakIna' 'peakIcabk' 'peakIpca' 'peakSRleak' ...
               'Naconc' 'CICRamp'};
output_table = array2table(output,'VariableNames',output_head);
writetable(output_table,'Matrics_addMeasures_all.csv','Delimiter',',');
