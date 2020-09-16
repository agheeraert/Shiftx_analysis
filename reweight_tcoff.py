import pandas as pd 
from Parser import ShiftXParser, expVt
import argparse
import mdtraj as md 
import numpy as np
from os import listdir
from os.path import join as jn
from tqdm import tqdm
import pickle as pkl

parser = argparse.ArgumentParser(description='Gets the table of aa temperature coefficients')
parser.add_argument('f',  type=str, nargs=2,
                    help='Two folders with SHIFTX2 results')
parser.add_argument('t',  type=str, nargs=3,
                    help='Two file containing trajectories and one file containing topology')
parser.add_argument('e',  type=str,
                    help='Experimental data')
parser.add_argument('o',  type=str,
                    help='output folder')                                       
parser.add_argument('-T',  type=float, default=20,
                    help='Difference of temperature during the computation')
parser.add_argument('-s', action='store_true', help='Computes also standard deviation')


args = parser.parse_args()
# N_CHUNKS = 200
# rmsd = []
# for chunk1, chunk2 in tqdm(zip(md.iterload(args.t[0], top=args.t[2], chunk=N_CHUNKS), md.iterload(args.t[1], top=args.t[2], chunk=N_CHUNKS))):
#     for i in range(N_CHUNKS):
#         rmsd.append(md.rmsd(chunk1, chunk2, frame=i)[i])
# rmsd = np.array(rmsd)

# np.save(jn(args.o, 'rmsd.npy'), rmsd)
# weights = np.exp(-rmsd)/np.sum(np.exp(-rmsd))
# np.save(jn(args.o, 'weights.npy'), weights)

weights = np.load(jn(args.o, 'weights.npy'))

# def get_frames(d):
#     n = len(listdir(d))
#     return [jn(d, 'frame_%s.pdb.cs'% i) for i in range(n-10000, n)]

# ts30 = ShiftXParser(get_frames(args.f[0]), ['H'], range(253), avg=False)
# pkl.dump(ts30, open(jn(args.o, 'ts30_reweight.spy'), 'wb'))
# ts50 = ShiftXParser(get_frames(args.f[1]), ['H'], range(253), avg=False)
# pkl.dump(ts50, open(jn(args.o, 'ts50_reweight.spy'), 'wb'))
ts30 = pkl.load(open(jn(args.o, 'ts30_reweight.spy'), 'rb'))
ts50 = pkl.load(open(jn(args.o, 'ts50_reweight.spy'), 'rb'))
tcoff = []
for i in range(ts50.shifts.shape[0]):
    tcoff.append(np.sum(np.multiply(weights, ts50.shifts[i] - ts30.shifts[i])))

if args.s:
    std = np.std((ts50.shifts-ts30.shifts)/args.T*10**3, axis=-1)
    
tcoff = np.array(tcoff)/args.T*10**3
np.save(jn(args.o, 'tcoff.npy'), tcoff)

# tcoff = np.load(jn(args.o, 'tcoff.npy'))
dH = ts30.df.copy()
dH['TCOFF'] = tcoff
if args.s:
    dH['STD'] = std

edh = np.genfromtxt(args.e, delimiter=',')
edh[edh== 0] = np.nan
indexes = np.where(~np.isnan(edh)==True)[0]
edh = edh[~np.isnan(edh)]
dh_exp = pd.DataFrame({'NUM': indexes, 'TCOFF': edh})
expVt(dh_exp, dH, jn(args.o, 'dH.png'), name='TCOFF', axisstring='Temp. coefficient', offset=0, inf=5, sup=-13)

L=[40, 41, 71, 74, 122, 123, 161, 195]
for elt in L:
    line = dH.loc[dH['NUM'] == elt-1]
    if not args.s:
        print(line['RES'].values[0]+str(line['NUM'].values[0]+1),line['TCOFF'].values[0])
    else:
        print(line['RES'].values[0]+str(line['NUM'].values[0]+1),line['TCOFF'].values[0], line['STD'].values[0])

