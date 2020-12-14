import pandas as pd
import sys
from scipy import stats

SAMPLE_TMP = sys.argv[1]
DD = sys.argv[2]
sgRNA_tmp = sys.argv[3]
TF = sys.argv[4]

FILE_tmp = DD + '/' + sgRNA_tmp + '/'+ SAMPLE_TMP + "_" + sgRNA_tmp + '_' + TF + '_raw_data.txt'
data_tmp = pd.read_table(FILE_tmp)
res_WMW  = stats.mannwhitneyu(data_tmp['deviation'][data_tmp['sgRNA'] == 'sgNTC'], data_tmp['deviation'][data_tmp['sgRNA'] == sgRNA_tmp])
Pvalue_tmp = res_WMW.pvalue
print(sgRNA_tmp + '__' + TF, res_WMW.pvalue)

