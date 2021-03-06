%x1 y1  x2  y2  voltage
%left
-18 -5.0  -16 -5.0  -55
-16 -5.0  -16 5.0   -55
-16 5.0   -18 5.0   -55
-18 5.0   -18 -5.0  -55

-14 -5.0  -6  -5.0  55
-6  -5.0  -6  -1.0   55
-6  -1.0    -14 -1.0    55
-14   -1.0  -14 -5.0    55

-14 1.0 -6  1.0 55
-6  1.0 -6  5.0 55
-6  5.0   -14 5.0   55
-14 5.0   -14 1.0  55

%right
16  -5.0  18  -5.0  -55
18  -5.0  18  5.0   -55
18  5.0   16  5.0   -55
16  5.0   16  -5.0  -55

6   -5.0  14  -5.0  55
14  -5.0  14  -1.0   55
14  -1.0    6   -1.0    55
6   -1.0    6   -5.0    55

6   1.0 14  1.0 55
14  1.0 14  5.0 55
14  5.0 6   5.0 55
6   5.0 6   1.0 55

%bottom
-5.0  -14 -1.0   -14 55
-1.0    -14 -1.0    -6  55
-1.0    -6  -5.0    -6  55
-5.0    -6  -5.0    -14 55

1.0 -14 5.0 -14 55
5.0 -14 5.0 -6  55
5.0 -6  1.0 -6  55
1.0 -6  1.0 -14  55

-5.0  -18 5.0   -18 -55
5.0   -18 5.0   -16 -55
5.0   -16 -5.0  -16 -55
-5.0  -16 -5.0  -18 -55

%top
-5.0  6   -1.0   6   55
-1.0    6   -1.0    14  55
-1.0    14  -5.0    14  55
-5.0    14  -5.0    6   55

1.0 6   5.0 6   55
5.0 6   5.0 14  55
5.0 14  1.0 14  55
1.0 14  1.0 6   55

-5.0  16  5.0   16  -55
5.0   16  5.0   18  -55
5.0   18  -5.0  18  -55
-5.0  18  -5.0  16  -55

%edge bottom right
6   -14 14  -14 55
14  -14 14  -6  55
14  -6  6   -6  55
6   -6  6   -14 55

6   -18 14  -18 -155
14  -18 14  -16 -155
14  -16 6   -16 -155
6   -16 6   -18 -155

%edge bottom left
-14 -14 -6  -14 55
-6  -14 -6  -6  55
-6  -6  -14 -6 55
-14 -6  -14 -14 55

-14 -18 -6  -18 -155
-6  -18 -6  -16 -155
-6  -16 -14 -16 -155
-14 -16 -14 -18 -155 


%edge top right
6   6   14  6   55
14  6   14  14  55
14  14  6   14  55
6   14  6   6   55

6   16  14  16  -155
14  16  14  18  -155
14  18  6   18  -155
6   18  6   16  -155

%edge top left
-14  6  -6  6   55
-6  6   -6  14  55
-6  14  -14 14  55
-14 14  -14 6   55

-14 16  -6  16  -155
-6  16  -6  18  -155 
-6  18  -14 18  -155
-14 18  -14 16  -155

