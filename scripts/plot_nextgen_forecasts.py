import numpy as np

files = ['Downloads/GFDL_A06_1982_2018_CCA_hindcast.txt', 'Downloads/CMC2_1982_2018_CCA_hindcast.txt', 'Downloads/CMC1_1982_2018_CCA_hindcast.txt', 'Downloads/GFDL_B01_1982_2018_CCA_hindcast.txt', 'Downloads/CFSv2_1_1982_2018_CCA_hindcast.txt', 'Downloads/COLA_1982_2018_CCA_hindcast.txt', 'Downloads/GFDL_aer04_1982_2018_CCA_hindcast.txt']

lats, flag = [], 0
ensemble_data = []
for filename in files:
    f = open(filename, 'r')
    ndx = 0
    vals, model_data = [], []
    prev_ndx = 1
    for line in f:
        if ndx == 0:
            xlmns = line
        elif ndx == 1:
            field, T, nrow, ncol, row, col, units, missing = line.strip().split(',')
        elif ndx == 2:
            long = line.split()
        else:
            if line[0:4] == 'cpt:' and ndx > 1:
                model_data.append(vals)
                vals = []
                prev_ndx = ndx
                flag = 1
            else:
                if ndx == prev_ndx + 1:
                    pass
                else:
                    line = line.strip().split()
                    vals.append([ float(i) for i in line[1:]])
                    if flag == 0:
                        lats.append(line[0])
        ndx += 1
    f.close()
    model_data.append(vals)
    model_data = np.asarray(model_data)
    ensemble_data.append(model_data)

ensemble_data = np.asarray(ensemble_data)
ensemble_data[np.where(ensemble_data == -999.0)] = np.nan

ensemble_mean = np.nanmean(ensemble_data, axis=0)
print(ensemble_mean.shape)

f = open('kjch_ensemble_mean.txt', 'w')
f.write('xmlns:cpt=http://iri.columbia.edu/CPT/v10/\n')
f.write('cpt:field=prcp, cpt:T=1982-06/09, cpt:nrow=32, cpt:ncol=28, cpt:row=Y, cpt:col=X, cpt:units=mm/day, cpt:missing=-999.000000000\n')
for i in range(ensemble_mean.shape[2]):
    f.write('\t{}'.format(long[i]))
f.write('\n')

for i in range(ensemble_mean.shape[1]):
    fmt_str = '\t{}\t' + ''.join(['{}\t' for j in range(ensemble_mean.shape[2])])
    t = fmt_str.format(*ensemble_mean[0, i]) + '\n'
    f.write(t)

for year in range(1, ensemble_mean.shape[0]):
    f.write('cpt:T={}-06/09\n'.format(1982+year))
    for i in range(ensemble_mean.shape[2] + 1):
        f.write('\t{}'.format(long[i]))
    f.write('\n')

    for i in range(ensemble_mean.shape[1]):
        fmt_str = '\tlat_placeholder\t' + ''.join(['{}\t' for j in range(ensemble_mean.shape[2])])
        f.write(fmt_str.format(*ensemble_mean[year, i, :]) + '\n')

f.close()
