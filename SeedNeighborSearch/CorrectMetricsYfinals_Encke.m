nparamset = 100000;
py_output = zeros(nparamset,49);
metrics_out = zeros(nparamset,33);
nexist = 0;
npassed= 0;
for iparamset = 1:nparamset
    filename1=strcat("./result/metrics",num2str(iparamset),".txt");
    filename2=strcat("./result/yfinal_",num2str(iparamset),".mat");
    try
        outputA = load(filename1);
        load(filename2);
        nexist  = nexist+1;
        params  = outputA(4:12);
        outline = [params, yfinal];
        py_output(nexist,:)= outline;
        metrics_out(nexist,:)=outputA;
        if(outputA(3)==1)
            npassed = npassed + 1;
        end
    catch
%        disp(["noyfinal" num2str(iparamset)]);
    end
end
disp(["passed param num = ",num2str(npassed)," outof ", num2str(nexist)]);
output_head = {'pCa' 'kfca' 'IbarNCX' 'ks' 'koCa' 'ec50SR' 'Vmax_SRCaP' 'Kmf' 'Kmr' ...
                'm'       'h'       'j'       'd'       'f'       'fcaBj'   'fcaBsl'   'xtos'    'ytos'    'xtof'    'ytof'    'xkr'     'xks'   ...
                'RyRr'    'RyRo'    'RyRi'    'NaBj'    'NaBsl'   'TnCL'    'TnCHc'   'TnCHm'   'CaM'     'Myoc'    'Myom'  ...
                'SRB'     'SLLj'   'SLLsl'    'SLHj'    'SLHsl'  'Csqnb'...
                'Ca_sr'   'Naj'     'Nasl'    'Nai'     'Ki'      'Caj'     'Casl'    'Cai'     'Vm'    'rtos'  ...
               };
output_table = array2table(py_output(1:nexist,:),'VariableNames',output_head);
writetable(output_table,'params_yfinal_SeedNeighbor_Encke.txt','Delimiter','\t');
    
output_head = {'Num' 'ifconverge' 'ifpass' 'pCa' 'kfca' 'IbarNCX' 'ks' 'koCa' 'ec50SR' 'Vmax_SRCaP' ...
               'Kmf' 'Kmr'            'SRCamax'   'SRCaRel'   'SRCaVmx'   'SRCaTau'   'Camin'   'Camax'   'CaTau' 'amp' ...
               'SRCaRfit' 'CaTfitR' 'M-SRCamax' 'M-SRCaRel' 'M-SRCaVmx' 'M-SRCaTau' 'M-Camin' 'M-Camax' 'M-CaTau' '#passed' ...
               'ICamax' 'ICaTau' 'APD'};
output_table = array2table(metrics_out(1:nexist,:),'VariableNames',output_head);
writetable(output_table,'TauModified_SeedNeighbor_Encke.csv','Delimiter',',');
