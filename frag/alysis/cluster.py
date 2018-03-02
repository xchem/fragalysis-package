#!/usr/bin/env python
import sys
import random
import math


class KMeans(object):

    def __init__(self, _X, _k, _xVal = 0, _stop=False):
        # X is sample size lists of dim length
        #
        # _xVal is the number of records to hold out cross-validation.
        # To use this you must randomize input data!
        #
        # Setting _stop=True causes iteration to stop when out of cross-validate
        # error starts to rise.
        #
        self.nFeatures = len(_X[0])
        self.xValSize = _xVal
        self.allSize = len(_X)
        self.size = self.allSize - self.xValSize
        self.X = _X
        self.k = _k
        self.stop = _stop
        # Initialize group memebership
        self.dataClusterId = [-1 for i in range(0, self.allSize)] # index of group for each data pair
        self.clusters = {}
        idx = 0
        # initialize to k random data points
        # don't assign x-val as a strat center
        for i in random.sample(list(range(0, self.size)), self.k):
            self.clusters[idx] = self.X[i]
            idx += 1
        # output records
        self.record = []
        self.errorRecord = []

    def dSquared(self, x, y):
        dist2 = 0.0
        for j,k in zip(x,y):
            dist2 += (j - k)**2
        return dist2

    def error(self):
        res = 0.0
        for i in range(0, self.size):
            res += self.dSquared(self.X[i], self.clusters[self.dataClusterId[i]])
        # error on non training data
        res1 = 0.0
        err1 = 0.0
        for i in range(self.size, self.allSize):
            res1 += self.dSquared(self.X[i], self.clusters[self.dataClusterId[i]])
        if res1 > 0.0:
            err1 = res1/self.xValSize
        return res/self.size, err1

    def nearestCluster(self, x):
        cmin = sys.maxsize
        cidx = -sys.maxsize
        for j in self.clusters:
            dist = math.sqrt(self.dSquared(x, self.clusters[j]))
            if dist < cmin:  # record closest centroid
                cmin = dist
                cidx = j
        return cidx, cmin

    def assign(self):
        for i in range(0, self.allSize):
            self.dataClusterId[i], dmin = self.nearestCluster(self.X[i])

    def updateClusters(self):
        ctemp = {} # dim sums by cluster
        for j in range(0, self.k):
            ctemp[j] = []
            for k in range(0, self.nFeatures):
                ctemp[j].append(0.0) # init sums
            ctemp[j].append(0) # init counter
        # only calculate clusters on training, not cross-validation set
        for i in range(0,self.size):
            for j in range(0, self.nFeatures):
                ctemp[self.dataClusterId[i]][j] += self.X[i][j]
            ctemp[self.dataClusterId[i]][self.nFeatures] += 1 # count
        for c in self.clusters:
            if ctemp[c][self.nFeatures] != 0:
                self.clusters[c] = [ ctemp[c][k]/ctemp[c][self.nFeatures] for k in range(0,self.nFeatures)]
            else:
                # no members in this cluster
                pass
        return

    def run(self, nmax = 100, eps = 1e-7):
        prev = 0.0
        prevXVal = float(sys.maxsize)
        for iter in range(0,nmax):
            # update assignments
            self.assign()
            # calculate error
            err, errXVal = self.error()
            #
            if self.stop and errXVal - prevXVal >= 0.0:
                sys.stderr.write("Cross-validation error increasing at step %d\n"%iter)
                break
            prevXVal = errXVal
            #
            if abs(err-prev) < eps:
                sys.stderr.write("Tolerance reached at step %d\n"%iter)
                break
            prev = err
            # going on...
            self.errorRecord.append((iter, err, errXVal))
            self.output(str(iter))
            self.updateClusters()
        sys.stderr.write("Iterations completed: %d\n"%iter)
        sys.stderr.write("Final error: %f\n"%prev)
        sys.stderr.write("Final cross-validation error: %f\n"%prevXVal)
        # This is a step past stop if using cross-validation...
        self.output("Final")
        return err, errXVal

    def output(self, iter):
        for i in range(0,self.size):
            self.record.append([str(y) for y in self.X[i]] + [str(self.dataClusterId[i])] + ["Iter-%s"%iter])
        for i in range(self.size, self.allSize):
            self.record.append([str(y) for y in self.X[i]] + [str(self.dataClusterId[i])] + ["Xval-Iter-%s"%iter])
        for k in self.clusters:
            self.record.append([str(y) for y in self.clusters[k]] + [str(k)] + ["Cent-Iter-%s"%iter])

    def getOutput(self):
        for x in self.record:
            yield x

    def getErrors(self):
        for x in self.errorRecord:
            yield x


class DPMeans(KMeans):
    def __init__(self, _X, _lam = 1, _xVal = 0, _stop=False):
        # init k-means with 1 cluster
        KMeans.__init__(self, _X, 1, _xVal, _stop)
        self.lam = _lam

    def assign(self):
        for i in range(0, self.size):
            cidx, dmin = self.nearestCluster(self.X[i])
            if dmin > self.lam:
                self.k += 1
                self.clusters[self.k-1] = self.X[i]
                self.dataClusterId[i] = self.k - 1
            else:
                self.dataClusterId[i] = cidx
        # don't create new clusters on cross-validation data
        for i in range(self.size, self.allSize):
            self.dataClusterId[i], dmin = self.nearestCluster(self.X[i])

    def error(self):
        err, xValErr = KMeans.error(self)
        return err + self.lam * self.k, xValErr + self.lam * self.k


def dp_means(data, lamb):
    """
    Carry out DP Means clustering on a vector of vectors
    :param data: vector of vectors (all same length)
    :param lamb: the parameter - the distance between clusters
    :return:
    """
    dp = DPMeans(data, lamb, 0, False)
    dp.run()
    return dp
