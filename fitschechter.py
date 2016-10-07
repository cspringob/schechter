#!/Users/christopherspringob/anaconda/bin/python

import matplotlib as mpl
import matplotlib.pyplot as pp
import matplotlib.ticker as ticker
import numpy as np

inputfile = "lumfunc.txt"
f1 = open("fitschechter.txt","w")

absmagd=[]
logphid=[]
error=[]
absmagarr=[]
phiarr=[]
logphiarr=[]

#The program reads in a data file that tells us the frequency of absolute magnitudes (the luminosity function), and then fits a Schechter function to that luminosity function.


for line in open(inputfile):
    if line[0] == "#":
        continue
    entries = list(line.strip().split())
    absmagd.append(float(entries[0]))
    logphid.append(float(entries[1]))
    error.append(float(entries[3]))

minchi2 = 10000.0
for j in range(0,10):
    alpha=-1.2+(0.02*j)
    for k in range(0,10):
        mstar=-23.5+(0.1*k)
        for l in range(0,10):
            phistar=10.0**(-3.5+(0.1*l))
            chi2sum=0
            for i in range(0,len(absmagd)):
                absmag=0.25*((1.0*i)+0.5)-25.75
                phi=0.4*np.log(10.0)*phistar*(10.0**(-0.4*(absmag-mstar)*(alpha+1.0)))*np.exp(-1.0*(10.0**(-0.4*(absmag-mstar))))
                logphi=np.log10(phi)
                chi2sum=chi2sum+((logphi-logphid[i])*(logphi-logphid[i])/(error[i]*error[i]))
            print(j,k,l,chi2sum,file=f1)
            if(chi2sum<minchi2):
                minchi2=chi2sum
                print(alpha,mstar,phistar,chi2sum)
                bestalpha=alpha
                bestmstar=mstar
                bestphistar=phistar


f1.close()

for i in range(0,30):
    absmagarr.append(0.25*((1.0*i)+0.5)-25.75)
    phiarr.append(0.4*np.log(10.0)*bestphistar*(10.0**(-0.4*(absmagarr[i]-bestmstar)*(bestalpha+1.0))*np.exp(-1.0*(10.0**(-0.4*(absmagarr[i]-bestmstar))))))
    logphiarr.append(np.log10(phiarr[i]))


#Write the entire array for the best fit Schechter function to the file bestfit.txt:
np.savetxt('bestfit.txt', np.transpose([absmagarr,phiarr,logphiarr]),fmt='%7.3f %8.6f %8.5f')

#Plot the best fit Schechter function, and save it to the file fitschechter.eps

fig = pp.figure()
ax = fig.add_subplot(111)
plot = ax.scatter(absmagarr,logphiarr)

pp.plot(absmagarr,logphiarr)

pp.savefig('fitschechter.eps',bbox_inches='tight',pad_inches=0.03)

