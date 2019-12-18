isempty
while (1)

    i=i+1;    
survey{i}.title = input('Title = ','s');    
survey{i}.index = input('Index *,* = ','s');    
survey{i}.journal = input('Journal, Volume, Number, Month, Year = ','s');    
survey{i}.methods_for_analysis = input('Methods for analysis (FEA,HA,Exp) = ','s');    
survey{i}.analysis_assumptions = input('Analysis assumptions (FEA,HA,Exp) = ','s');    
survey{i}.modeling_domains = input('Modeling domains (electronic, mechanical, electrostatic, magnetic) = ','s');    
survey{i}.sugar_now = input('Sugar now = ','s');    
survey{i}.minor_modifications = input('Minor modification = ','s');    
survey{i}.major_modifications = input('Major modification = ','s');    
survey{i}.device_appications = input('Device applications = ','s');    
survey{i}.institution = input('Institution = ','s');    
end
