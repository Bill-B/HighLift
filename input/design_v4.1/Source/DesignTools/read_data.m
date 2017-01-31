function [ method, des ] = read_data(  )
%Read initial data from design.txt (CAD or FreeWake input method)
%Uses: input method (CAD or FreeWake)
%Returns: runcond struct, timecond struct, aerocond struct, check (y/n), number of wings wings, current wingcount
loc = 0;

%% READ DATA FROM design.txt %%

%open up design.txt

des = fopen('design.txt', 'r');


%% READ MODEL DESIGN DATA FROM design.txt %%

%% FIND MODEL DESIGN SECTION %%
while strcmp('Model-Design',loc) == 0    
    loc = fscanf(des, '%s',1);
end

%ok, found Model-Design Section, continue....
%which design method are we using?
while loc ~='='    
    loc = fscanf(des, '%s', 1);
end
method = fscanf(des,'%s',1); % CAD or FreeWake?
if strcmp('CAD',method) == 0 && strcmp('FreeWake',method) == 0
   error('Design Method must be either CAD or FreeWake');
end
loc = 0;

% %do we want to run checkmesh?
% while loc ~='='    
%     loc = fscanf(des, '%s', 1);
% end
% check = fscanf(des,'%s',1); % run checkmesh? y/n
% loc = 0;

%read number of wings
% while loc ~='='    
%     loc = fscanf(des, '%s', 1);
% end

%store data
% wings = fscanf(des,'%d',1); %number of wings
wingcount = 0;

end

