import pystan
import sys
import numpy as np

def main():
    test_data, test_labels = load_data('digitstest.txt')
    train_data, train_labels = load_data('digitstrain.txt')
    nb_dat = {'N': len(train_labels),
              'N_test': len(test_labels),
              'D': train_data.shape[1],
              'label': train_labels,
              'train': train_data,
              'test': test_data}
    pystan.misc.stan_rdump(nb_dat, 'test_nb.dat')
    """fit=pystan.stan(file='./test_nb.stan', data=nb_dat, iter=100)
    print(fit.extract()['mu'])
    means = np.asarray(fit.extract()['mu'])
    print means.shape
    means = np.mean(means, axis=0)
    means1 = np.mean(np.asarray([x for (x, y) in zip(train_data, train_labels) if y == 1]))
    means0 = np.mean(np.asarray([x for (x, y) in zip(train_data, train_labels) if y == 0]))
    delta0 = means[0] - means0
    delta1 = means[1] - means1
    print delta0 + delta1
    print np.mean(delta0 + delta1)"""



def load_data(filepath):
    data = []
    labels = []
    with open(filepath) as file:
        for line in file:
            tokens = line.strip().split(',')
            line = np.asarray(map(float, tokens[:-1]))
            data.append(line)
            labels.append(int(tokens[-1]))


    labels = np.array(labels)
    data = np.array(data)
    return data, labels

if __name__ == '__main__':
    sys.exit(main())
