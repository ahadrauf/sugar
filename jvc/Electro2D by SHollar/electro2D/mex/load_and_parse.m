function [variable] = load_and_parse(FNAME)

load(FNAME);

seg_flag=0;
for x = 1:size(FNAME,2);
     if (FNAME(x)=='.') & (~seg_flag)
	F_SHORT = FNAME(1:(x-1));
	seg_flag=1;
      end
end
if (~seg_flag)
     F_SHORT = FNAME;
end

variable = eval(F_SHORT);
clear eval(F_SHORT);

