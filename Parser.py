import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import linregress
mpl.style.use('classic')
from os.path import join as jn
from collections import Counter
from operator import itemgetter
import pickle as pkl

class ShiftXParser():
    """Class parsing results from one or multiple ShiftX file"""
    def __init__(self, files, atoms, residues, empty=False, avg=True):
        """ Parses files (can be list or str) only for given atoms and residues
        """
        self.atoms = atoms
        self.residues = residues
        if type(files) == str:
            self.parse_struct(files)
            self.shifts = Shifts(self.parse_shifts(files))
        elif type(files) == list:
            L_shifts = []
            self.parse_struct(files[0])
            for f in tqdm(files):
                L_shifts.append(self.parse_shifts(f))
            self.shifts = np.stack(L_shifts, axis=-1)
        # self.avg_discard_wrong()
        if avg:
            self.avg_shifts()
    
    def parse_struct(self, f):
        df = pd.read_csv(f)
        self.indexes = np.where((df['ATOMNAME'].isin(self.atoms)) & (df['NUM'].isin(self.residues)))
        df = df.loc[self.indexes]
        # self.idx2res = pd.Series(df['RES'].values, df.index).to_dict()
        # self.idx2atom = pd.Series(df['ATOMNAME'].values, df.index).to_dict()
        # self.indexes = df.index
        # self.columns = df.columns
        # self.res = df['RES']
        # self.atoms = df['ATOMNAME']
        self.df = df

    def parse_shifts(self, f):
        df = pd.read_csv(f)
        return df.loc[self.indexes]['SHIFT']

    def plot_shifts_time(self, output='/home/aria/landslide/RESULTS/IGPS_TDEP/TIMEPLOT/30/'):
        for i in range(self.shifts.shape[0]):
            line = self.df.loc[self.df['NUM']==i//2]
            try:
                string = line['RES'].values[0]+str(line['NUM'].values[0]+1)
            except IndexError:
                print(line)
                print(i)
                string = str(i)
            if i%2==0:
                string+='_H'
            else:
                string+='_N'
            f = plt.figure()
            plt.xlim(0,len(self.shifts[i]))
            plt.xlabel('Frame')
            plt.ylabel('$^1H$ chemical shift')
            plt.scatter(range(self.shifts.shape[1]), self.shifts[i], marker='+')
            plt.savefig(jn(output, '%s.png' %string))
            plt.close()

    # def avg_discard_wrong(self):
    #     for i in range(self.shifts.shape[0]):
    #         c = dict(Counter(self.shifts[i]))
    #         # wrongval = np.argmax(np.array([v for k,v in sorted(c.items())]))
    #         wrongval = max(c.items(), key=itemgetter(1))[0]
    #         self.shifts[i][self.shifts[i]==wrongval] = np.nan
    #     self.df['SHIFT'] = np.nanmean(self.shifts, axis=-1)
    #     self.shifts = Shifts(self.df)

    def avg_shifts(self):
        self.df['SHIFT'] = np.mean(self.shifts, axis=-1)
        self.df['STD_SHIFT'] = np.std(self.shifts, axis=-1)
        self.shifts = Shifts(self.df)

class ExpParser():
    """Class parsing experimental results"""
    def __init__(self, f):
        table = np.genfromtxt(f, delimiter=',')
        self.residues = table[:,0].astype(int)
        shifts = np.zeros(table.shape[0]*2)
        shifts[::2] = table[:,1]
        shifts[1::2] = table[:,2]
        numtwice = np.zeros(self.residues.shape[0]*2, dtype=int)
        numtwice[::2] = self.residues
        numtwice[1::2] = self.residues
        h = ['H' for elt in range(self.residues.shape[0])]
        n = ['N' for elt in range(self.residues.shape[0])]
        atoms = ['a' for elt in range(self.residues.shape[0]*2)]
        atoms[::2] = h
        atoms[1::2] = n
        # self.res2idx = dict(zip(self.residues, range(self.residues.shape[0])))
        self.shifts = Shifts(pd.DataFrame({'NUM': numtwice, 'ATOMNAME': atoms, 'SHIFT':shifts}))

class Shifts():
    """class to handle shifts"""
    def __init__(self, df):
        if type(df) == str:
            self.df = pd.read_pickle(df)
        else:
            self.df = df
        self.h = self.df[self.df['ATOMNAME'] == 'H']
        self.n = self.df[self.df['ATOMNAME'] == 'N']
    def save(self, output):
        self.df.to_pickle(output)

def temperature_coefficient(shifts1, shifts2, dT=20):
    df = shifts1
    S1 = np.zeros(len(shifts1))
    S2 = np.zeros_like(S1)
    nums = np.zeros_like(S1)
    for i, elt in enumerate(shifts1.itertuples()):
        num = elt[1]
        nums[i] = num
        S1[i] = elt[-1]
        S2[i] = shifts2.loc[shifts2['NUM'] == num]['SHIFT'].values[0]
    df = pd.DataFrame({'NUM': nums.astype(int), 'RES': shifts1['RES'], 'TCOFF':(S2-S1)/dT*10**3})
    return df

def getXY(exp_values, t_values, offset=1, name='SHIFT'):
    f = plt.figure()
    X, Y = [], []
    for elt in exp_values.itertuples():
        num = elt[1] 
        eshift = elt[-1]
        try:
            sshift = t_values.loc[t_values['NUM'] == num-offset][name].values[0]
            X.append(eshift)
            Y.append(sshift)
        except IndexError:
            pass
    return X, Y

def expVt(exp_values, t_values, output, name='SHIFT', axisstring='$^1$H Chemical shift', text=None, offset=1, inf=6, sup=10):
    X, Y = getXY(exp_values, t_values, offset=offset, name=name)
    plt.style
    plt.scatter(X, Y, marker='s', color='k')
    plt.plot([inf, sup], [inf, sup], color='r')
    plt.xlabel('%s from Exp. (ppm)' %axisstring)
    plt.ylabel('%s from SHIFTX2 (ppm)' %axisstring)
    plt.xlim(sup+1, inf-0.5)
    plt.ylim(sup+1, inf-0.5)
    X = np.array(X).reshape(1, -1)
    Y = np.array(Y).reshape(1, -1)
    if text:
        plt.text(sup-0.05, inf-0.2, text, fontsize=14)
    plt.savefig(output)
    plt.close()
    return X, Y

# def expVt_bar(exp_values, t_values, output, labels, offset=0):
#     """ only for temperature coefficient (not good)"""
#     exp, th = getXY(exp_values, t_values, offset=offset, name='TCOFF')
#     x = np.arange(len(labels))  # the label locations
#     width = 0.5  # the width of the bars
#     fig, ax = plt.subplots()
#     exp = np.array(exp)
#     th = np.array(th)
#     exp = exp[~np.isnan(exp)]
#     th = th[~np.isnan(th)]
#     rects1 = ax.bar(x - width/2, exp, width, label='Experimental')
#     rects2 = ax.bar(x + width/2, th, width, label='Theory')
#     ax.set_ylabel('Temperature coefficient')
#     ax.set_title('Comparison between experimental and theoretical temperature coefficient')
#     ax.set_xticks(x)
#     ax.set_xticklabels(labels)
#     ax.legend()
#     autolabel(rects1, ax)
#     autolabel(rects2, ax)
#     fig.tight_layout()
#     plt.savefig(output) 
#   

def expVt_lines(exp_values, t_values, output, labels, dT=20.5, offset=0, moredata=None):
    exp, th = getXY(exp_values, t_values, offset=offset, name='TCOFF')
    if moredata:
        exp2, th2 = getXY(moredata[0], moredata[1], offset=offset)
    if not moredata:
        Y1 = np.array(exp)*dT/10**3
        Y2 = np.array(th)*dT/10**3
    if moredata:
        Y1 = np.array(exp)[:-4]/np.array(exp2)*dT/10**3
        Y2 = np.array(th)[:-4]/np.array(th2)*dT/10**3
    print(np.mean(Y1), np.mean(Y2))
    x = np.arange(len(labels))
    fig, ax = plt.subplots()
    # ax.plot(x[:-4], Y1*100, color='g')
    # ax.plot(x[:-4], Y2*100, color='r')
    ax.plot(x[:-4], (Y2-Y1)*100, color='purple')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    fig.tight_layout()
    plt.savefig(output)


def expVt_bar(exp_values, t_values, output, labels, offset=0):
    """ only for temperature coefficient (not good)"""
    exp, th = getXY(exp_values, t_values, offset=offset, name='TCOFF')
    perc = (np.array(th)-np.array(exp))/np.array(exp)*100
    verbose_expVt(perc)
    x = np.arange(len(labels))  # the label locations
    fig, ax = plt.subplots()
    rects = ax.bar(x, perc, width=1)
    ax.set_ylabel('Difference between experience and theory')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    autolabel(rects, ax)
    fig.tight_layout()
    plt.savefig(output)

def verbose_expVt(perc):
    d = pkl.load(open('/home/aria/landslide/RESULTS/IGPS_TDEP/DICT/np2res.p', 'rb'))
    for elt in perc.argsort()[-20:][::-1]:
        print(elt*2, d[elt])

def rate(exp_values, t_values):
    """ only for temperature coefficient (not good)"""
    exp, th = getXY(exp_values, t_values, offset=1, name='TCOFF')
    def inrange(num):
        return num <= -1 and num >= -4
    truepos, falsepos, falseneg, trueneg = 0, 0, 0, 0
    for i in range(len(exp)):
        if inrange(th[i]):
            if inrange(exp[i]):
                truepos += 1
            else:
                falsepos += 1
        else:
            if inrange(exp[i]):
                falseneg += 1
            else:
                trueneg += 1
    # print(truepos, falsepos, falseneg, trueneg)
    print(truepos/(truepos+falsepos))


def autolabel(rects, ax):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = round(rect.get_height(), 1) 
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')





        

if __name__ == '__main__':
    DD_PATH = '/home/aria/landslide/RESULTS/IGPS_TDEP/'
    es25 = ExpParser(jn(DD_PATH, 'exp_nmr/25.csv'))
    es45 = ExpParser(jn(DD_PATH, 'exp_nmr/45.csv'))    
    ts30 = ShiftXParser([jn(DD_PATH, 'shiftx/30/frame_%s.pdb.cs') %i for i in range(1, 10025)], ['H', 'N'], range(253))
    # ts30.plot_shifts_time()
    ts50 = ShiftXParser([jn(DD_PATH, 'shiftx/50/frame_%s.pdb.cs') %i for i in range(1, 10052)], ['H', 'N'], range(253))
    # ts50.plot_shifts_time(output='/home/aria/landslide/RESULTS/IGPS_TDEP/TIMEPLOT/50/')
    ts30.shifts.save(jn(DD_PATH, 'ts30.pdy'))
    ts50.shifts.save(jn(DD_PATH, 'ts50.pdy'))
    # expVt(es25.shifts.h, ts30.shifts.h, jn(DD_PATH, 'VSTHEORY/25v30_2.png'), text='297.78K')
    # expVt(es45.shifts.h, ts50.shifts.h, jn(DD_PATH, 'VSTHEORY/45v50_2.png'), text='317.48K')
    # dH = temperature_coefficient(ts30.shifts.h, ts50.shifts.h)
    # dN = temperature_coefficient(ts30.shifts.n, ts50.shifts.n)

    # ts30 = Shifts(jn(DD_PATH, 'ts30.pdy'))
    # ts50 = Shifts(jn(DD_PATH, 'ts50.pdy'))
    # ts30.h.to_csv(jn(DD_PATH, 'VSTHEORY/30H_disc.csv'))
    # ts50.h.to_csv(jn(DD_PATH, 'VSTHEORY/50H_disc.csv'))
    # expVt(es25.shifts.h, ts30.h, jn(DD_PATH, 'VSTHEORY/25v30.png'), text='297.78K')
    # expVt(es45.shifts.h, ts50.h, jn(DD_PATH, 'VSTHEORY/45v50.png'), text='317.48K')
    # dH = temperature_coefficient(ts30.h, ts50.h, dT=20.5)
    # dN = temperature_coefficient(ts30.n, ts50.n, dT=20.5)

    # enh1 = np.genfromtxt(jn(DD_PATH, 'exp_nmr/nh1.csv'), delimiter=',')
    # enh1[enh1== 0] = np.nan
    # edh = np.genfromtxt(jn(DD_PATH, 'exp_nmr/dh.csv'), delimiter=',')
    # edh[edh== 0] = np.nan
    # indexes = np.where(~np.isnan(edh)==True)[0]
    # edh = edh[~np.isnan(edh)]
    # dh_exp = pd.DataFrame({'NUM': indexes, 'TCOFF': edh})
    # expVt(dh_exp, dH, jn(DD_PATH, 'VSTHEORY/dH_3.png'), name='TCOFF', axisstring='Temp. coefficient', offset=0, inf=5, sup=-13)
    # expVt_lines(dh_exp, dH, jn(DD_PATH, 'VSTHEORY/dH_lines_3.png'), indexes, offset=0, moredata=[es25.shifts.h, ts30.h])
    # expVt_bar(dh_exp, dH, jn(DD_PATH, 'VSTHEORY/dH_bar2_normalized.png'), indexes, offset=0)
    # rate(dh_exp, dH)
    # nh1_exp = pd.DataFrame({'NUM': range(253), 'TCOFF': enh1})
    # nh1_shifty = pd.DataFrame({'NUM': dH['NUM'], 'TCOFF':(dH['TCOFF']+0.2*dN['TCOFF'])})
    # dH.to_csv(jn(DD_PATH, 'VSTHEORY/dhshifty.csv'))
    # nh1_exp.to_csv(jn(DD_PATH, 'VSTHEORY/nh1exp.csv'))
    # nh1_shifty.to_csv(jn(DD_PATH, 'VSTHEORY/nh1shifty.csv'))
    # expVt(nh1_exp, nh1_shifty, jn(DD_PATH, 'VSTHEORY/nh1.png'), name='TCOFF', axisstring='Temp. coefficient', offset=0, inf=4, sup=-12)

            



        

        




        
