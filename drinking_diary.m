BrAC=zeros(length(MergedData.Breathlyzer_Report.RESULT),2);
BrAC(:,2)=transpose(MergedData.Breathlyzer_Report.RESULT);
TAC=zeros(length(MergedData.SCRAM_Report(1).TAC),2);
TAC(:,2)=transpose(MergedData.SCRAM_Report(1).TAC);
BrAC_t=datetime(MergedData.Breathlyzer_Report.TESTDATENUM,'ConvertFrom','datenum');
TAC_t=datetime(MergedData.SCRAM_Report(1).DateNum,'ConvertFrom','datenum');

BrAC(:,1)=[29 55 65 75 85 95 105 116 125 135 145 155 165 175 185 195 205 215 225 235 245 255 265 275 285 295 305 315]./60;
TAC(:,1)=[0 30 31 36 42 47 53 56 62 66 67 72 76 82 86 87 92 96 106 117 122 126 136 145 156 166 176 186 186 195 205 215 221 225 235 245 255 256 265 275 285 285 295 305 315 315 346 376 407 437 468 498 529 559 590 620 651 681 712 742 773 803 834 864 895 925 956 986 1017 1047 1077]./60;
