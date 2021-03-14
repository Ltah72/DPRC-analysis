#!/apps/psyc/anaconda3/bin

import pandas as pd
import sys

#print('Arguments:', len(sys.argv))
#print('List:', str(sys.argv))

if len(sys.argv) < 3:
    print('Too few arguments, please specify a filename for input and output')

filename = sys.argv[1]
outname = sys.argv[2]
df=pd.read_table(filename)
df=df[['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z']]
df.to_csv(outname,header=False,index=False)


#Example usage: fmriprep2hcp_mvmt_regressors <prefix>_desc-confounds_regressors.tsv <output>.par

#from FSL forum here: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;6e71614b.1911
