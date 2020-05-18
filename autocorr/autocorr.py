import numpy as np

def autocorr(data, npoints):
    y_vals = data[:,1]
    y_avg = y_vals.mean()
    dy = y_vals - y_avg
    AC = np.zeros(npoints)
    for i in range(AC.shape[0]):
        perc = (i+1)/AC.shape[0] * 100.0
        if perc % 10.0 == 0:
            print(f'{perc:5.2f}% complete')
        yy0 = dy*np.roll(dy, i)
        AC[i] = yy0[i:(i+1+int(data.shape[0]/2))].mean()/dy.std()**2
    return AC


if __name__ == '__main__':
    pass
