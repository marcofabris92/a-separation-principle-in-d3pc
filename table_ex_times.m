% Run this script after executing MAIN and saving the workspace in order
% to get the numerical results collected in Table 1 of this paper.


clear all
clc

pathw6 = '..\model1_SNR6_wrkspc.mat'; % set your path here
load(pathw6)
close all
etw6 = ave_exe_times;

pathw20 = '..\model1_SNR20_wrkspc.mat'; % set your path here
load(pathw20)
close all
etw20 = ave_exe_times;

pathc6 = '..\model3_SNR6_wrkspc.mat'; % set your path here
load(pathc6)
close all
etc6 = ave_exe_times;

pathc20 = '..\model3_SNR20_wrkspc.mat'; % set your path here
load(pathc20)
close all
etc20 = ave_exe_times;

a = (etw6+etw20+etc6+etc20)/4;
s = sqrt(((a-etw6).^2+(a-etw20).^2+(a-etc6).^2+(a-etc20).^2)/4);

training = cell(1,9);
offline_search = cell(1,9);
optimization = cell(1,9);


% MPC
i = 1;
training{i} = [NaN];
offline_search{i} = [NaN];
optimization{i} = [a(1) s(1)];

% 2 offline
i = i + 1;
atr = a(13)+a(18);
str = sqrt(((atr-etw6(13)-etw6(18)).^2+...
    (atr-etw20(13)-etw20(18)).^2+...
    (atr-etc6(13)-etc6(18)).^2+...
    (atr-etc20(13)-etc20(18)).^2)/4);
training{i} = [atr str];
offline_search{i} = [a(2)/2 s(2)/2];
optimization{i} = [a(5) s(5)];

% 3 offline
i = i + 1;
training{i} = [atr str];
offline_search{i} = [a(2)/2 s(2)/2];
optimization{i} = [a(6) s(6)];

% 23 offline
i = i + 1;
training{i} = [atr str];
offline_search{i} = [a(3) s(3)];
optimization{i} = [a(7) s(7)];

% DeePC
i = i + 1;
atr = a(17)+a(18);
str = sqrt(((atr-etw6(17)-etw6(18)).^2+...
    (atr-etw20(17)-etw20(18)).^2+...
    (atr-etc6(17)-etc6(18)).^2+...
    (atr-etc20(17)-etc20(18)).^2)/4);
training{i} = [atr str];
offline_search{i} = [a(4) s(4)];
optimization{i} = [a(8) s(8)];

% opt
i = i + 1;
atr = a(14)+a(18);
str = sqrt(((atr-etw6(14)-etw6(18)).^2+...
    (atr-etw20(14)-etw20(18)).^2+...
    (atr-etc6(14)-etc6(18)).^2+...
    (atr-etc20(14)-etc20(18)).^2)/4);
training{i} = [atr str];
offline_search{i} = [NaN];
optimization{i} = [a(9) s(9)];

% thm2 (subopt)
i = i + 1;
atr = a(15)+a(18);
str = sqrt(((atr-etw6(15)-etw6(18)).^2+...
    (atr-etw20(15)-etw20(18)).^2+...
    (atr-etc6(15)-etc6(18)).^2+...
    (atr-etc20(15)-etc20(18)).^2)/4);
training{i} = [atr str];
offline_search{i} = [NaN];
optimization{i} = [a(10) s(10)];

% 2 online
i = i + 1;
atr = a(16)+a(18);
str = sqrt(((atr-etw6(16)-etw6(18)).^2+...
    (atr-etw20(16)-etw20(18)).^2+...
    (atr-etc6(16)-etc6(18)).^2+...
    (atr-etc20(16)-etc20(18)).^2)/4);
training{i} = [atr str];
offline_search{i} = [NaN];
optimization{i} = [a(11) s(11)];

% 3 online
i = i + 1;
training{i} = [atr str];
offline_search{i} = [NaN];
optimization{i} = [a(12) s(12)];







% final table in [s]
scale = 1e3;
for i = 1:9
    training{i} = scale*training{i};
    offline_search{i} = scale*offline_search{i};
    optimization{i} = scale*optimization{i};
end
training
offline_search
optimization