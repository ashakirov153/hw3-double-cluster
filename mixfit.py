#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from scipy import stats


tau = 0.5
mu1 = 0.5
sigma1 = 0.2
mu2 = 1.5
sigma2 = 0.7
n = 100000

x_1 = np.random.normal(mu1, sigma1, int(tau*n))

x_2 = np.random.normal(mu2, sigma2, int((1-tau)*n))

x = np.concatenate((x_1, x_2))

def t(x, tau, mu1, mu2, sigma12, sigma22):
    tau0 = tau
    tau1 = 1 - tau
    t1 = tau0 / np.sqrt(2*np.pi * sigma12) * np.exp(-0.5*((x - mu1)**2)/sigma12)
    t2 = tau1 / np.sqrt(2*np.pi * sigma22) * np.exp(-0.5*((x - mu2)**2)/sigma22)
    t_norm = t1 + t2
    t1 = np.divide(t1, t_norm, out=np.full_like(t_norm, 0.5), where=t_norm!=0)
    t2 = np.divide(t2, t_norm, out=np.full_like(t_norm, 0.5), where=t_norm!=0)
    return np.vstack((t1,t2))

def theta(x, *old):
    t1, t2 = t(x, *old)
    tau = np.sum(t1) / np.sum(t1+t2)
    mu1 = np.sum(x * t1) / np.sum(t1)
    mu2 = np.sum(x * t2) / np.sum(t2)
    sigma12 = np.sum((x - mu1)**2 * t1) / np.sum(t1)
    sigma22 = np.sum((x - mu2)**2 * t2) / np.sum(t2)
    return tau, mu1, mu2, sigma12, sigma22

def max_likelihood(x, tau, mu1, sigma1, mu2, sigma2, rtol=1e-3):
    pass


def em_double_gauss(x, tau, mu1, sigma1, mu2, sigma2, rtol=1e-3):
    th = (tau, mu1, mu2, sigma1**2, sigma2**2)
    for i in range(100):
        if (np.abs(np.array(th) - theta(x, *th)) < 10**(-5)).all:
            break
        th = theta(x, *th)
    return (th[0], th[1], th[3]**0.5, th[2], th[4]**0.5)


def em_double_cluster(x, tau1, tau2, muv, mu1, mu2, sigma02, sigmax2, sigmav2, rtol=1e-5):
    pass


if __name__ == "__main__":
    pass
