import numpy as np 
import matplotlib.pyplot as plt 
from scipy.fft import fft, fftshift, fftfreq
from scipy.optimize import curve_fit
sigma = 1000.0
kappa = 1000.0
prefactor = 1
prefactor = 2/(np.pi)
# Compute the size of each sample using a single data file

datadir = "data"
# data = np.loadtxt(f"{datadir}/curvradii_1.txt", delimiter=" ")
# data = np.loadtxt(f"{datadir}/elastic_1.txt", delimiter=" ")
# data = np.loadtxt(f"{datadir}/radii2_1.txt", delimiter=" ")
data = np.loadtxt(f"{datadir}/radii_1.txt", delimiter=" ")
# data = np.loadtxt(f"{datadir}/curv_1.txt", delimiter=" ")
hmean = data[:,1].mean()
dh = data[:,1] - hmean 
# plt.plot(dh)
# plt.show()
N = dh.size
xf = fftshift(fftfreq(N))
fft_dh = fftshift(fft(dh))
N = fft_dh.size
samples = 1000
normsqs = np.zeros([samples, N])
for i in range(samples):
# i=0
    data = np.loadtxt(f"{datadir}/radii_{i}.txt", delimiter=" ")
    # data = np.loadtxt(f"{datadir}/radii2_{i}.txt", delimiter=" ")
    # data = np.loadtxt(f"{datadir}/curv_{i}.txt", delimiter=" ")
    # data = np.loadtxt(f"data/curvradii_{i}.txt", delimiter=" ")
    # data = np.loadtxt(f"data/elastic_{i}.txt", delimiter=" ")
    hmean = data[:,1].mean()
    dh = data[:,1] - hmean 
    # plt.plot(dh)
    # plt.show()
    N = dh.size
    xf = 2*np.pi*fftshift(fftfreq(N))
    fft_dh = 1/np.sqrt((2*np.pi))*fftshift(fft(dh))
    normsq = np.abs(fft_dh)**2
    normsqs[i] = normsq

# normsqs = np.log(normsqs)
avg = np.mean(normsqs, axis=0)
std = np.std(normsqs, axis=0)
# plt.plot(xf, normsq)
N = xf.size
xf = xf[N//2+1:]
avg = avg[N//2+1:]
std = std[N//2+1:]
# N = xf.size
# xf = xf[N//2+1:]
# avg = avg[N//2+1:]
# plt.loglog(xf, avg,'o-', label="data")


[fw, fh] = plt.rcParams["figure.figsize"]
plt.figure(figsize=(fw/2, fh/2))

plt.plot(xf, avg,'o-', color='g', label="data")
plt.fill_between(xf, avg-std/2, avg+std/2, color='g', alpha=0.3)

plt.xscale("log")
plt.yscale("log")
def func1(x, sigma):
    return 1 / (sigma * np.abs(x)**2)

def func2(x, kappa):
    return 1 / (kappa * np.abs(x)**4)

def func3(x, sigma, kappa):
    return 1 / (sigma * np.abs(x)**2 + kappa * np.abs(x)**4)

popt1, pcov1 = curve_fit(func1, xf, avg, p0=[1e4])
popt2, pcov2 = curve_fit(func2, xf, avg, p0=[1e6])
# popt3, pcov3 = curve_fit(func3, xf, avg, p0=[1e4, 1e6])
# plt.plot(xf, func1(xf, *popt1), '--',
#          label="fit: 1 / (%.1e "% tuple(popt1) + r"$q^{-2}$" + ")")
# plt.plot(xf, func2(xf, *popt2), '--',
#          label="fit: 1 / (%.1e "% tuple(popt2) + r"$q^{-4}$" + ")")
# plt.plot(xf, func3(xf, *popt3), '--',
#          label="fit: 1 / (%.1e "% tuple(popt3)[0] + r"$q^{-2}$" + " + %.1e "% tuple(popt3)[1] + r"$q^{-4}$" + ")")
# plt.loglog(xf, 1e-6/np.abs(xf)**4, '--', label=r"$q^{-4}$")
# plt.loglog(xf, 1e-3/np.abs(xf)**2, '--', label=r"$q^{-2}$")

plt.plot(xf, prefactor/(sigma*np.abs(xf)**2), '--', label=r"$k_B T/ (\sigma q^2)$")
plt.plot(xf, prefactor/(kappa*np.abs(xf)**4), '--', label=r"$k_B T/ (\kappa q^4)$")
plt.xlabel(r"$q$")
plt.ylabel(r"$\langle|\delta h(q)|^2\rangle$")
plt.plot(xf, prefactor/(sigma*np.abs(xf)**2 + kappa*np.abs(xf)**4), 'k--', label=r"$(\sigma q^2 + \kappa q^4)^{-1}$")
# plt.plot(xf[:N//2], avg[0:N//2])
plt.legend()
plt.tight_layout()
# plt.savefig("CurvSpectrum.png")
# plt.savefig("CurvSpectrum.pdf")
# plt.savefig("LineTensionSpectrum.png")
# plt.savefig("LineTensionSpectrum.pdf")
plt.show()
