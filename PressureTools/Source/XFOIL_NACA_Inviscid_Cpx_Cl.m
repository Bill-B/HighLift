function [ x, Cp ] = XFOIL_NACA_Inviscid_Cpx_Cl( NACA, Cl )
%XFOIL_ Summary of this function goes here
%   Detailed explanation goes here
%%


NACAs = sprintf('%04i', NACA);


exePath = 'C:\Users\Bill\Desktop\XFOIL6.99\xfoil.exe';


   
    fileID = fopen('XFOIL_Commands.txt','w');
    fprintf(fileID,'plop\n');
    fprintf(fileID,'g\n');
    fprintf(fileID,'\n');
    fprintf(fileID,'naca\n');
    fprintf(fileID,'%s\n', NACAs);
    fprintf(fileID,'oper\n');

    fprintf(fileID,'init\n');
    fprintf(fileID,'Cl\n');
    fprintf(fileID,'%.3f\n', Cl);
    
    fprintf(fileID,'cpwr\n');
    fprintf(fileID,'%s\n', 'XFOIL_Cpx_Output.dat');
    fprintf(fileID,'\n');
    
    fclose(fileID);
    
    %% RUN XFOIL
    [~,~] = system(sprintf('%s <XFOIL_Commands.txt> output.txt',exePath));
    
    %% Read XFOIL_Cpx_Output.dat
    fileID = fopen('XFOIL_Cpx_Output.dat','r');
    formatSpec = '%f%f%f%[^\n\r]';
    dataArray = textscan(fileID, formatSpec, 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'HeaderLines' ,3, 'ReturnOnError', false);
    fclose(fileID);
    
    %% Try to delete output
    try
        delete('XFOIL_Cpx_Output.dat');
        delete('XFOIL_Commands.txt');
        delete('output.txt');
    end
    
    %% Reformat data
    x = [dataArray{:,1}];
    y = [dataArray{:,2}];
    Cp = [dataArray{:,3}];

    


end













