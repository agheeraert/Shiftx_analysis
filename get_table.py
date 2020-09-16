import pandas as pd 
from Parser import Shifts, temperature_coefficient
import argparse

parser = argparse.ArgumentParser(description='Gets the table of aa temperature coefficients')
parser.add_argument('f1',  type=str, 
                    help='First .pdy file with shifts')
parser.add_argument('f2',  type=str,
                    help='Second .pdy file with shifts')
parser.add_argument('-T',  type=float, default=20,
                    help='Difference of temperature during the computation')
parser.add_argument('-L',  type=int, default=[40, 41, 71, 74, 122, 123, 161, 195], nargs='+',
                    help='List of residues to output')
parser.add_argument('-s', action='store_true', help='Gets the standard deviation (if applicable)')

args = parser.parse_args()
ts30 = Shifts(args.f1)
ts50 = Shifts(args.f2)
dH = temperature_coefficient(ts30.h, ts50.h, dT=args.T)

for elt in args.L:
    line = dH.loc[dH['NUM'] == elt-1]
    print(line['RES'].values[0]+str(line['NUM'].values[0]+1),line['TCOFF'].values[0])

if args.s:
    print('Standard deviations')
    for elt in args.L:
        line1 = ts30.h.loc[dH['NUM'] == elt-1]
        line2 = ts50.h.loc[dH['NUM'] == elt-1]
        print(line1['RES'].values[0]+str(line1['NUM'].values[0]+1),line1['STD_SHIFT'].values[0])
        print(line2['RES'].values[0]+str(line2['NUM'].values[0]+1),line2['STD_SHIFT'].values[0])
        


