############################## IMPORT
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
import math
from math import log
import cmath
import matplotlib.ticker as tck
import numpy
import scipy
from scipy import interpolate
from scipy.optimize import curve_fit
from pymongo import MongoClient
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker
from pylab import figure, axes, pie, title, show
import sys, os
import matplotlib.ticker as mtick
#### Path
path='/home/thomas/CLL2Species/Database_online/pictures_summary/'
################################ class CALCULATION
class calculation:
####################### GLOBAL VARIABLES
    db=None
    def __init__(self):
        client=MongoClient('localhost',27017)
      	uri = "mongodb://antonello:diocancaro@ds062339.mlab.com:62339/two_species_db"
      	client=MongoClient(uri)
        db=client.two_species_db
        dati=db.Result2Species
        self.db=dati
########### SUB-FUNCTION
##### Max List
    def pi_val(self, x,pos): # this function give the number for the x-axis
        return str(int(x/cmath.pi)) + "$\pi$"

    def get_max_list(self, infornate, parametro_gruppo, parametro_max):
        """
        this function return the list of results for a maximum value of a parameter given a parameter to group by
        """
        inf = dict()
        inf_list=list()
        for i in infornate:
            val = inf.get(str(i[parametro_gruppo]))
            if (val is not None and \
                    val[parametro_max]<i[parametro_max]) or\
                    val is None:
                inf[str(i[parametro_gruppo])] = i
        for i in inf.values():
            inf_list.append(i)
        return sorted(inf_list, key=lambda k:k["L"])

    def get_max_list_float(self, infornate, parametro_gruppo, parametro_max):
        """
        this function return the list of results for a maximum value of a parameter given a parameter to group by
        here there is a convertion the parameter in a float to have a good sort of the values
        """
        inf = dict()
        inf_list=list()
        for i in infornate:
            val = inf.get(str(i[parametro_gruppo]))
            if (val is not None and \
                    float(val[parametro_max])<float(i[parametro_max])) or\
                    val is None:
                inf[str(i[parametro_gruppo])] = i
        for i in inf.values():
            inf_list.append(i)
        return sorted(inf_list, key=lambda k:k["L"])

##### First Brillouin Zone
    def FBZ(self,L):
        '''
        Here we find the values of k from the FBZ and imposing the  boundary conditions.

        FBZ:
                    [-pi,pi]
        ----|----|--(--|--)--|----|----|----|----|-> k
                       0     1

        |k| in/ [-pi,pi] we can shift |k| /in [0,2pi].

        boundary conditions:

        ____|____|____|_|_|_|_|____|____|____|____|
            |         |       |         |    |
                      0      2pi

        |k| /in [0,2pi]

                      k = 2pim/L                m in/ Z

        '''
        k=list()
        lung=range(L)
        for l in lung:
            k.append(2*math.pi*l/L)
        return k
#### Number of Cluster
    def cluster(self,L,density,rc):
        M = (1.0-density)*float(L)/float(rc)
        return M.is_integer()


####################### FUNCTIONS #######################
    ##### ENTROPY #####
    def entropy_convi_fit(self,L,sweep,v,u,rc): ## here I compute for only one value
        '''
        This function plot the entropy in function of dmrg states.
        The argument "L" is a LIST
        '''
        infornate=self.db.find({"sweep":sweep,"L":L,"v":v,"u":u,"rc":rc})
        list_entropy=list()
        list_sts=list()
        for l in infornate:
            list_entropy.append(l["entropy"])
            list_sts.append(l["n_sts"])
        print list_entropy
        print list_sts
        fig=plt.figure()
        plt.plot(list_sts,list_entropy) # improve the plotting
        fig.suptitle('L_'+str(L)+'_v_'+str(v)+'_rc_'+str(rc), fontsize=20)
        plt.xlabel('DMRG states', fontsize=18)
        plt.ylabel('entropy', fontsize=16)
        plt.show()

    def entropy_convi_(self,L,sweep,v,u,rc,rho1,rho2): ## here I compute for only one value
        '''
        This function plot the entropy in function of dmrg states.
        The argument "L" is a LIST
        '''
        infornate=self.db.find({"sweep":sweep,"L":L,"v":v,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
        list_entropy=list()
        list_sts=list()
        ee=list()
        for l in infornate:
            list_entropy.append(l["entropy"][0])
            list_sts.append(l["n_sts"])
        #list_sts=sorted(list_sts, key=lambda n_sts:n_sts)
        print list_entropy
        print list_sts
        #for p in range(0,len(list_entropy)):
        #    ee.append(list_entropy[p][0])
        #print list_sts
        fig=plt.figure()
        plt.plot(list_sts,list_entropy,'r8',linewidth=3) # improve the plotting
        fig.suptitle('L_'+str(L)+'_v_'+str(v)+'_rc_'+str(rc), fontsize=20)
        plt.xlabel('DMRG states', fontsize=18)
        plt.ylabel('entropy', fontsize=16)
        plt.show()

    def entropy_convi_list(self,L,sweep,v,u,rc): ## here I compute for only one value
        '''
        This function plot the entropy in function of dmrg states.
        The argument "L" is a LIST
        '''
        fig=plt.figure()
        c=['r','k','b','c','g','m','y']
        marker=['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X','.']
        coppia=list()
        for colore in c:
            for m in marker:
                coppia.append([colore,m])
        k=0
        for  LL in L:
            infornate=self.db.find({"sweep":sweep,"L":LL,"v":v,"u":u,"rc":rc})
            list_entropy=list()
            list_sts=list()
            axLL=fig.add_subplot(1, 1, 1)
            for l in infornate:
                list_entropy.append(l["entropy"])
                list_sts.append(l["n_sts"])
            name="Lenght_"+str(LL)
            axLL.scatter(list_sts,list_entropy,c=coppia[k][0],marker=coppia[k][1],label=name)
            k+=1
        print list_sts
        fig.suptitle('L_'+str(L)+'_v_'+str(v)+'_rc_'+str(rc), fontsize=20)
        plt.xlabel('DMRG states', fontsize=18)
        plt.ylabel('entropy', fontsize=16)
        plt.show()
        plt.legend(loc='upper left')
        plt.show()

    def entropy_convi_list_size(self,sweep,v,u,rc,L,rho1,rho2,fit=False):
        '''
        Computation of the ground State Entropy versus L (system size)
        '''
        fig=plt.figure()
        c=['r','k','b','c','g','m','y']
        marker=['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X','.']
        coppia=list()
        for colore in c:
            for m in marker:
                coppia.append([colore,m])

        list_sts=list()
        for  LL in range(1,L/8+1):
            list_sts2=list()
            ls=8*LL
            infornate=self.db.find({"sweep":sweep,"L":ls,"v":v,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
            #infornate=self.get_max_list(inf,"L","n_sts")
            for l in infornate:
                list_sts2.append(l["n_sts"])
            list_sts.append(list_sts2)
        #print list_sts
        d={}
        for x in range(0,len(list_sts[0])):
            d["list_entropy_{0}".format(x)]=list()
            d["list_Lsize_{0}".format(x)]=list()
            d["plot_{0}".format(x)]=plt
        for  Ll in range(1,L/8+1):
            LL=8*Ll
            for i in range(0,len(list_sts[0])):
                    inf=[]
                    inf=self.db.find({"sweep":sweep,"v":v,"u":u,"rc":rc,"L":LL,"rho1":rho1,"rho2":rho2,"n_sts":list_sts[Ll-1][i]})
                    for l in inf:
                        d["list_entropy_{0}".format(i)].append(l["entropy"][0])
                        d["list_Lsize_{0}".format(i)].append(l["L"])

        for i in range(0,len(list_sts[0])):
                d["list_Lsize_{0}".format(i)]=[y*3.0/(math.pi) for y in d["list_Lsize_{0}".format(i)]]
                d["list_Lsize_{0}".format(i)]=[log(z) for z in d["list_Lsize_{0}".format(i)]]
        fig=plt.figure()
        for i in range(0,len(list_sts[0])):
            d["plot_{0}".format(i)]=plt.plot(d["list_Lsize_{0}".format(i)],d["list_entropy_{0}".format(i)],'s--',label='#states'+str(list_sts[3][i]))

        if fit:
            d["list_entropy_0"][3]=float(d["list_entropy_0"][3])
            popt, pcov = curve_fit(fit, d["list_Lsize_0"], d["list_entropy_0"])
            xfine = numpy.linspace( d["list_Lsize_0"][0], d["list_Lsize_0"][-1],100 )  # define values to plot the function for
            plt.plot(xfine, fit(xfine, popt[0], popt[1]), 'r-', linewidth=2)
            # fit plot
            ax=fig.add_subplot(1,1,1)
            box1 = TextArea("a = "+str(popt[0])+" +/- "+str(math.sqrt(pcov[0,0]))+"\n"+"b = "+str(popt[1])+" +/- "+str(math.sqrt(pcov[1,1])),textprops=dict(color="k"))

            box = HPacker(children=[box1],align="center", pad=5, sep=5)

            anchored_box = AnchoredOffsetbox(loc=3,child=box, pad=0.,frameon=True,bbox_to_anchor=(0., 1.02),bbox_transform=ax.transAxes,borderpad=0.)

            ax.add_artist(anchored_box)

            fig.subplots_adjust(top=0.8)

        fig.suptitle('L_'+str(L)+'_v_'+str(v)+'_rc_'+str(rc), fontsize=20)
        plt.xlabel('k(L)', fontsize=18)
        plt.ylabel('3*S(L/2)', fontsize=16)
        plt.legend(loc=2)
        plt.show()

    ##### GROUND STATE ENERGY #####
    def GroundStateSizeDstate(self,sweep,v,u,rc,L,Lmax,rho1,rho2):
        '''
        Ground state energy per site versus L (system size)
        '''
        list_sts=list()
        for  LL in range(1,L/8+1):
            list_sts2=list()
            ls=8*LL
            infornate=self.db.find({"sweep":sweep,"L":ls,"v":v,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
            #infornate=self.get_max_list(inf,"L","n_sts")
            for l in infornate:
                list_sts2.append(l["n_sts"])
            list_sts.append(list_sts2)
        d={}
        for x in range(0,len(list_sts[0])):
            d["list_energy_{0}".format(x)]=list()
            d["list_Lsize_{0}".format(x)]=list()
            d["plot_{0}".format(x)]=plt
        for  Ll in range(1,L/8+1):
            LL=8*Ll
            for i in range(0,len(list_sts[0])):
                    inf=self.db.find({"sweep":sweep,"v":v,"u":u,"rc":rc,"L":LL,"rho1":rho1,"rho2":rho2,"n_sts":list_sts[Ll-1][i]})
                    for l in inf:
                        if isinstance(l["E"][0], float):
                            d["list_energy_{0}".format(i)].append(l["E"][0]/l["L"])
                            d["list_Lsize_{0}".format(i)].append(l["L"])
                        else:
                            d["list_energy_{0}".format(i)].append(float(l["E"][0][0])/float(l["L"]))
                            d["list_Lsize_{0}".format(i)].append(l["L"])

        for i in range(0,len(list_sts[0])):
                d["list_Lsize_{0}".format(i)]=[1./y for y in d["list_Lsize_{0}".format(i)]]
                #d["list_energy_{0}".format(i)]=[z/L for z in d["list_energy_{0}".format(i)]]
        fig=plt.figure()
        for i in range(0,len(list_sts[0])):
            d["plot_{0}".format(i)]=plt.plot(d["list_Lsize_{0}".format(i)],d["list_energy_{0}".format(i)],'s--',label='#states'+str(list_sts[3][i]))

        fig.suptitle('L_'+str(L)+'_v_'+str(v)+'_rc_'+str(rc), fontsize=20)
        plt.xlabel('1/L', fontsize=18)
        plt.ylabel('E0/L', fontsize=16)
        plt.legend(loc=5)
        plt.show()


    def SpectraGap(self,sweep,v,u,rc,L,rho1,rho2): ## here I compute for only one value
        '''
        This function plot the entropy in function of L system size for 1 sweep only
        '''
        fig=plt.figure()
        c=['r','k','b','c','g','m','y']
        marker=['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X','.']
        coppia=list()
        for colore in c:
            for m in marker:
                coppia.append([colore,m])

        list_sts=list()
        for  LL in range(1,L/8+1):
            list_sts2=list()
            ls=8*LL
            infornate=self.db.find({"sweep":sweep,"L":ls,"v":v,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
            #infornate=self.get_max_list(inf,"L","n_sts")
            for l in infornate:
                list_sts2.append(l["n_sts"])
            list_sts.append(list_sts2)
        print list_sts
        d={}
        for x in range(0,len(list_sts[0])):
            d["list_energy_{0}".format(x)]=list()
            d["list_Lsize_{0}".format(x)]=list()
            d["plot_{0}".format(x)]=plt
        for  Ll in range(1,L/8+1):
            LL=8*Ll
            for i in range(0,len(list_sts[0])):
                    inf=self.db.find({"sweep":sweep,"v":v,"u":u,"rc":rc,"L":LL,"rho1":rho1,"rho2":rho2,"n_sts":list_sts[Ll-1][i]})
                    for l in inf:
                        if isinstance(l["E"][0], float):
                            d["list_energy_{0}".format(i)].append(math.fabs(l["E"][0]-l["E"][1]))
                            d["list_Lsize_{0}".format(i)].append(l["L"])
                        else:
                            d["list_energy_{0}".format(i)].append(math.fabs(float(l["E"][0][0])-float(l["E"][1][0])))
                            d["list_Lsize_{0}".format(i)].append(l["L"])

        #print list_entropy3
        #print d["list_entropy_2"]
        #print d["list_Lsize_3"]

        for i in range(0,len(list_sts[0])):
                d["list_Lsize_{0}".format(i)]=[1./y for y in d["list_Lsize_{0}".format(i)]]
                #d["list_energy_{0}".format(i)]=[z*1./L for z in d["list_energy_{0}".format(i)]]
        fig=plt.figure()
        for i in range(0,len(list_sts[0])):
            d["plot_{0}".format(i)]=plt.plot(d["list_Lsize_{0}".format(i)],d["list_energy_{0}".format(i)],'s--',label='#states'+str(list_sts[3][i]))
        #p1=plt.plot(d["list_Lsize_0"],d["list_entropy_0"],'rh--',label='#states'+str(list_sts[3][0]))
        #p2=plt.plot(d["list_Lsize_1"],d["list_entropy_1"],'g^--',label='#states'+str(list_sts[3][1]))
        #p3=plt.plot(d["list_Lsize_2"],d["list_entropy_2"],'ks--',label='#states'+str(list_sts[3][2]))
        #p4=plt.plot(d["list_Lsize_3"],d["list_entropy_3"],'b8--',label='#states'+str(list_sts[3][3]))

        fig.suptitle('L_'+str(L)+'_v_'+str(v)+'_rc_'+str(rc), fontsize=20)
        plt.xlabel('1/L', fontsize=18)
        plt.ylabel('$\Delta_{E0}$', fontsize=16)
        plt.legend(loc=5)
        plt.show()

    def energy_convi_list_size(self,sweep,v,u,rc,L): ## here I compute for only one value
        '''
        This function plot the entropy in function of L system size for 1 sweep only
        '''
        fig=plt.figure()
        c=['r','k','b','c','g','m','y']
        marker=['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X','.']
        coppia=list()
        for colore in c:
            for m in marker:
                coppia.append([colore,m])
        list_entropy=list()
        list_Lsize=list()
        entropy=list()
        for  Ll in range(1,L/8+1):
            LL=8*Ll
            inf=self.db.find({"sweep":sweep,"v":v,"u":u,"rc":rc,"L":LL})
            infornate=self.get_max_list(inf, "L","n_sts")
            for l in infornate:
                list_entropy.append(l["E"])
                list_Lsize.append(l["L"])
        for p in range(0,len(list_entropy)):
            entropy.append(list_entropy[p][0])
        fig=plt.figure()
        plt.plot(list_Lsize,entropy,'rh--') # improve the plotting
        fig.suptitle('L_'+str(L)+'_v_'+str(v)+'_rc_'+str(rc), fontsize=20)
        plt.xlabel('L', fontsize=18)
        plt.ylabel('energy', fontsize=16)
        plt.show()

    def entropy_convi(self,L,sweep,v,u,rc): ## ENTROPY VS DMRG STATES
        '''
        for entropy function I will use the parameters:
           sweep *
           L(single value)
           v
         * this is linked to the data precision(we can choose te better value to use studying the oters parameters)
        we want build two methods that:
           1)Plot in terms of number of states
        '''
        if isinstance(L,list) and isinstance(L[0],int):
            entropy_convi(self,L,sweep,v,u,rc)
        elif isinstance(L,int):
            entropy_convi(self,L,sweep,v,u,rc)
        else:
            print "Error: the value of (L) inserted is not correct are accepted only INT or LIST of int "
    ####
    def entropy(self,sweep,rc,v,u): ## ENTROPY VS L
        inf=self.db.find({"sweep":sweep,"v":v,"u":u,"rc":rc}) ## TODO CHECK
        infornata=self.get_max_list(inf, "L","n_sts")
        list_entropy=list()
        list_L=list()
        for l in infornata:
            list_entropy.append(l["entropy"])
            list_L.append(l["L"])
        fig=plt.figure()
        plt.plot(list_L,list_entropy, "ro")
        fig.suptitle('v_'+str(v)+'rc'+str(rc),fontsize=20)
        plt.xlabel('sites', fontsize=18)
        plt.ylabel('entropy', fontsize=16)
        plt.show()
############ ENERGY
    def energy_convi_single(self,L,sweep,v,u,rc,rho1,rho2): ## here I compute for only one value
        '''
        This function plot the ground state energy in function of dmrg states.
        The argument "L" is a LIST

        To plot a fit you have to put a function in "f" defined like a function that you think could be the best
        '''
        find=self.db.find({"sweep":sweep,"L":L,"v":v,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
        infornate = list()
        for i in find:
            infornate.append({"n_sts":int(i["n_sts"]),"E": i["E"]})
        infornate=sorted(infornate, key=lambda x:x["n_sts"])
        list_energytmp=list()
        list_energy=list()
        list_sts=list()
        for l in infornate:
            if isinstance(l["E"][0], float):
                list_energy.append(l["E"][0])
                list_sts.append(l["n_sts"])
            else:
                list_energy.append(l["E"][0][0])
                list_sts.append(l["n_sts"])

        # here I define an order inside the dmrg states
        #for p in range(0,len(list_energytmp)):
        #        list_energy.append(list_energytmp[p][0])
        print list_energy
        print list_sts
        fig=plt.figure()
        ax=fig.add_subplot(1,1,1)
        #ax.set_yticks(list_energy)
        ax.yaxis.set_major_formatter(tck.ScalarFormatter())
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        ax.get_yaxis().get_major_formatter().set_scientific(False)
        plt.plot(list_sts,list_energy,marker='o', color='r',linestyle='--') # improve the plotting
#compute fit
        """  if f:
            print "fit"
            popt, pcov = curve_fit(f, list_sts, list_entropy)
            ax.annotate("a="+str(popt[0])+"+/-"+str(math.sqrt(pcov[0,0])),xy=(list_sts[1],list_entropy[1]))
            ax.annotate("b="+str(popt[1])+"+/-"+str(math.sqrt(pcov[1,1])),xy=(list_sts[2],list_entropy[1]))
# fit plot
            xfine = np.linspace(float(list_sts[0]), float(list_sts[-1]), 10)  # define values to plot the function for
            plt.plot(xfine, f(xfine, popt[0], popt[1],popt[2]), 'b--')
        """
        fig.suptitle('L_'+str(L)+'_v_'+str(v)+'_rc_'+str(rc), fontsize=20)
        plt.xlabel('DMRG states', fontsize=18)
        plt.ylabel('E0', fontsize=16)
        plt.show()
##
    def energy_convi_list(self,L,sweep,v,u,rc,rho1,rho2): ## here I compute for only one value
            '''
            This function plot the ground state energy in function of dmrg states.
            The argument "L" is a LIST
            '''
            fig=plt.figure()
            c=['r','k','b','c','g','m','y']
            marker=['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X','.']
            coppia=list()
            for colore in c:
                for m in marker:
                    coppia.append([colore,m])
            k=0
            for  LL in L:
                infornate=self.db.find({"sweep":sweep,"L":LL,"v":v,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
                list_entropy=list()
                list_sts=list()
                axLL=fig.add_subplot(1, 1, 1)
                for l in infornate:
                    list_entropy.append(l["E0"])
                    list_sts.append(l["n_sts"])
                name="Lenght_"+str(LL)
                axLL.scatter(list_sts,list_entropy,c=coppia[k][0],marker=coppia[k][1],label=name)
                k+=1
            fig.suptitle('L_'+str(L)+'_v_'+str(v)+'_rc_'+str(rc), fontsize=20)
            plt.xlabel('states', fontsize=18)
            plt.ylabel('E0', fontsize=16)
            plt.show()
            plt.legend(loc='upper left')
            plt.show()

    def calabrese_cardy_scal(self,L,sweep,rc,v,u,rho1,rho2,see,cut,fit=False):## ENTROPY VS K(l)
        '''
        This function plot the results of entropy in terms K(l) a variable definite from the calabrese-cardy formula of
        entanglement entropy
        cut is to select the last block size and avoid strong size effect
        '''
        inf=self.db.find({"sweep":sweep,"L":L,"v":v,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
        infornata=self.get_max_list_float(inf, "L","n_sts")[0]
        print infornata["n_sts"]
        k=[]
        entrop=[]
        for l in infornata["entangl_entr"]:
            k.append(l[1])
            entrop.append(l[2])
        k1=k[-cut:]
        #I print k1 to see the size of the block I have selected
        print k1
        k11=[]
        kk=[]
        for q in k1:
            k11.append(numpy.log2(numpy.sin(q*math.pi/L)*L/math.pi))
        entrop1=entrop[-cut:]

        fig=plt.figure()
        if fit:
            popt, pcov = curve_fit(fit, k11, entrop1)
            xfine = numpy.linspace( min(k11), max(k11),100 )  # define values to plot the function for
            plt.plot(xfine, fit(xfine, popt[0], popt[1]), 'k-', linewidth=3)
            ax=fig.add_subplot(1,1,1)

            #multipication by 3 to obtain centrale charge (following Calabrese Formula)
            box = TextArea("c = "+str(popt[0]*3)+" +/- "+str(math.sqrt(pcov[0,0]))+"\n"+"a0 = "+str(popt[1])+" +/- "+str(math.sqrt(pcov[1,1])),textprops=dict(color="k",size=20))
            box = HPacker(children=[box],align="center", pad=5, sep=5)
            anchored_box = AnchoredOffsetbox(loc=3,child=box, pad=0.,frameon=True,bbox_to_anchor=(0., 1.02),bbox_transform=ax.transAxes,borderpad=0.)
            ax.add_artist(anchored_box)

            fig.subplots_adjust(top=0.8)


        plt.plot(k11,entrop1,"ro",markersize=10)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.xlabel('$\kappa (\ell)$', fontsize=30)
        plt.ylabel('$S_L(\ell)$', fontsize=30)
        if see==0:
            plt.show()
        else:
            plt.savefig('/home/thomas/CLL2Species/scanU20/pictures/entropy_rc'+str(rc)+'_V'+str(v)+'_U'+str(u)+'_na'+str(rho1)+'_L'+str(L)+'.pdf',bbox_inches='tight')

##
    def calabrese_cardy_scal_compa(self,L,sweep,rc,u,rho1,rho2,see,inset,cut,fit=False):## ENTROPY VS K(l)
        '''
        This function will extract by a linear fitting the central charge for each value of interaction V
        cut is to select number of block size you will take (see previous function calabrese_cardy_scal)
        inset is to plot in certain range of V the central charge, here between 5.4 and 6.5
        '''
        list_v=list()
        infornate=self.db.find({"sweep":sweep,"L":L,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
        for l in infornate:
            list_v.append(l["v"])
        list_v=sorted(list_v) #sorte for the value of interaction
        d={}
        c=list()
        v=list()
        for x in range(0,len(list_v)):
            d["list_entropy_{0}".format(x)]=list()
            d["list_Lsize_{0}".format(x)]=list()
            d["plot_{0}".format(x)]=plt
        for i in range(0,len(list_v)):
                inf=self.db.find({"sweep":sweep,"u":u,"rc":rc,"L":L,"rho1":rho1,"rho2":rho2,"v":list_v[i]})
                infornata=self.get_max_list_float(inf, "L","n_sts")[0]
                for l in infornata["entangl_entr"]:
                        d["list_entropy_{0}".format(i)].append(l[2])
                        d["list_Lsize_{0}".format(i)].append(l[1])

                d["list_Lsize_{0}".format(i)]=d["list_Lsize_{0}".format(i)][-cut:]
                d["list_entropy_{0}".format(i)]=d["list_entropy_{0}".format(i)][-cut:]

        for i in range(0,len(list_v)):
                d["list_Lsize_{0}".format(i)]=[numpy.log2(numpy.sin(y*math.pi/L)*L/math.pi) for y in d["list_Lsize_{0}".format(i)]]

        if fit:
            for i in range(0,len(list_v)):
                popt, pcov = curve_fit(fit,d["list_Lsize_{0}".format(i)] , d["list_entropy_{0}".format(i)])
                if inset==1:
                        if list_v[i]<=6.5 and list_v[i]>=5.4:
                            c.append(popt[0]*3)
                            v.append(list_v[i])
                else:
                        c.append(popt[0]*3)
                        v.append(list_v[i])
        fig=plt.figure()
        plt.plot(v,c,"k^--",markersize=10,linewidth=3)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.xlabel('V/t', fontsize=30)
        plt.ylabel('$c_1$', fontsize=30)
        if see==0:
            plt.show()
        else:
            plt.savefig('/home/thomas/CLL2Species/scanU20/pictures/C1_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'_L'+str(L)+'.pdf',bbox_inches='tight')
    def centralcharge_vs_V_L(self,choice,L,rc,u,v,rho1,rho2,inset,cut,see,fit=False):## ENTROPY VS K(l)
        '''
        This function will extract by a linear fitting the central charge for each value of interaction V
        cut is to select number of block size you will take (see previous function calabrese_cardy_scal)
        inset is to plot in certain range of V the central charge, here between 5.4 and 6.5
        choice parameter allow to compute central charge for L or U
        '''
        d={}

        if choice==1:
            cutt=[7 , 12]
            infornata=self.db.find({"v":v,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
            list_L=list()
            for l in infornata:
                list_L.append(l["L"])
            list_L=sorted(list_L)
            lab1="L="+str(list_L[0])
            lab2="L="+str(list_L[1])
        else:
            cutt=[7 , 7]
            infornata=self.db.find({"v":v,"L":L,"rc":rc,"rho1":rho1,"rho2":rho2})
            list_L=list()
            for l in infornata:
                list_L.append(l["u"])
                list_L=sorted(list_L)
            lab1="U="+str(list_L[0])
            lab2="U="+str(list_L[1])
        #print list_L


        for y in range(0,len(list_L)):
            d["list_c_{0}".format(y)]=list()
            d["list_inter_{0}".format(y)]=list()
            d["max_{0}".format(y)]=[0,0]
            d["plot_{0}".format(y)]=plt

        for j in range(0,len(list_L)):
            list_v=list()
            if choice==1:
                infornate=self.db.find({"L":list_L[j],"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
            else:
                infornate=self.db.find({"L":L,"u":list_L[j],"rc":rc,"rho1":rho1,"rho2":rho2})
            for l in infornate:
                list_v.append(l["v"])
            list_v=sorted(list_v) #sorte for the value of interaction
            c=list()
            v=list()
            for x in range(0,len(list_v)):
                d["list_entropy_{0}".format(x)]=list()
                d["list_Lsize_{0}".format(x)]=list()
            if choice==1:
                for i in range(0,len(list_v)):
                        inf=self.db.find({"u":u,"rc":rc,"L":list_L[j],"rho1":rho1,"rho2":rho2,"v":list_v[i]})
                        infornata=self.get_max_list_float(inf, "L","n_sts")[0]
                        for l in infornata["entangl_entr"]:
                                d["list_entropy_{0}".format(i)].append(l[2])
                                d["list_Lsize_{0}".format(i)].append(l[1])
                        d["list_Lsize_{0}".format(i)]=d["list_Lsize_{0}".format(i)][-cutt[j]:]
                        d["list_entropy_{0}".format(i)]=d["list_entropy_{0}".format(i)][-cutt[j]:]
                for i in range(0,len(list_v)):
                        d["list_Lsize_{0}".format(i)]=[numpy.log2(numpy.sin(y*math.pi/list_L[j])*list_L[j]/math.pi) for y in d["list_Lsize_{0}".format(i)]]
            else:
                for i in range(0,len(list_v)):
                        inf=self.db.find({"u":list_L[j],"rc":rc,"L":L,"rho1":rho1,"rho2":rho2,"v":list_v[i]})
                        infornata=self.get_max_list_float(inf, "L","n_sts")[0]

                        for l in infornata["entangl_entr"]:
                                d["list_entropy_{0}".format(i)].append(l[2])
                                d["list_Lsize_{0}".format(i)].append(l[1])
                        d["list_Lsize_{0}".format(i)]=d["list_Lsize_{0}".format(i)][-cutt[j]:]
                        d["list_entropy_{0}".format(i)]=d["list_entropy_{0}".format(i)][-cutt[j]:]

                for i in range(0,len(list_v)):
                        d["list_Lsize_{0}".format(i)]=[numpy.log2(numpy.sin(y*math.pi/L)*L/math.pi) for y in d["list_Lsize_{0}".format(i)]]



            if fit:
                for i in range(0,len(list_v)):
                    popt, pcov = curve_fit(fit,d["list_Lsize_{0}".format(i)] , d["list_entropy_{0}".format(i)])
                    if inset==1:
                            if list_v[i]<=6.5 and list_v[i]>=5.4:
                                c.append(popt[0]*3)
                                v.append(list_v[i])
                    else:
                            d["list_c_{0}".format(j)].append(popt[0]*3)
                            d["list_inter_{0}".format(j)].append(list_v[i])
                            if d["max_{0}".format(j)][1] < popt[0]*3 :
                                        #max_xy = [k, sommatoria_a.real]
                                        d["max_{0}".format(j)]=[list_v[i], popt[0]*3]
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.annotate("{}".format(d["max_{0}".format(0)][0]), xy=d["max_{0}".format(0)],fontsize=25)
        ax.annotate("{}".format(d["max_{0}".format(1)][0]), xy=d["max_{0}".format(1)],fontsize=25)
        plt.plot(d["list_inter_{0}".format(1)],d["list_c_{0}".format(1)],"g^--",markersize=10,linewidth=3, label=lab2 )
        plt.plot(d["list_inter_{0}".format(0)],d["list_c_{0}".format(0)],"rD--",markersize=10,linewidth=3 ,label=lab1)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.xlabel('V/t', fontsize=30)
        plt.ylabel('$c_1$', fontsize=30)
        plt.legend(loc='upper left',fontsize=25)
        if see==0:
            plt.show()
        else:
            plt.savefig(os.path.join(path,'C1_rc'+str(rc)+'_choice_'+str(choice)+'_U'+str(u)+'_na'+str(rho1)+'_L'+str(L)+'.png'),bbox_inches='tight')

    def entropy_scal_compa(self,L,sweep,rc,u,rho1,rho2,see,cut,fit=False):## ENTROPY VS K(l)

        '''
        This function will extract by a linear fitting the central charge for each value of interaction V
        cut is to select number of block size you will take (see previous function calabrese_cardy_scal)
        inset is to plot in certain range of V the central charge, here between 5.4 and 6.5
        HERE different formula to extract see note project
        '''

        entrotest=list()
        list_v=list()
        infornate=self.db.find({"sweep":sweep,"L":L,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
        #infornate=self.get_max_list_float(inf,"L","n_sts")
        for l in infornate:
            list_v.append(l["v"])
        list_v=sorted(list_v)
        d={}
        c=list()
        SL=list()
        SL1=list()
        v=list()
        for x in range(0,len(list_v)):
            d["list_entropy_{0}".format(x)]=list()
            d["list_Lsize_{0}".format(x)]=list()
            d["plot_{0}".format(x)]=plt
        for i in range(0,len(list_v)):
                inf=self.db.find({"sweep":sweep,"u":u,"rc":rc,"L":L,"rho1":rho1,"rho2":rho2,"v":list_v[i]})
                infornata=self.get_max_list_float(inf, "L","n_sts")[0]
                for l in infornata["entangl_entr"]:
                        d["list_entropy_{0}".format(i)].append(l[2])
                        d["list_Lsize_{0}".format(i)].append(l[1])
                for p in range(2,len(infornata["entangl_entr"])):
                    if d["list_Lsize_{0}".format(i)][p]==L/2:
                        SL.append(d["list_entropy_{0}".format(i)][p])
                    if d["list_Lsize_{0}".format(i)][p]==L/2-1:
                        SL1.append(d["list_entropy_{0}".format(i)][p])
                v.append(list_v[i])
        for i in range(0,len(list_v)):
                c.append((3*(SL1[i]-SL[i]))/numpy.log2(numpy.cos(numpy.pi/L)))
        fig=plt.figure()
        plt.plot(v,c,"gD--",markersize=10,linewidth=3)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.xlabel('V/t', fontsize=30)
        plt.ylabel('$c_2$', fontsize=30)
        if see==0:
            plt.show()
        else:
            plt.savefig('/home/thomas/CLL2Species/scanU20/pictures/C2_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'_L'+str(L)+'.pdf',bbox_inches='tight')

    def entropy2_scal_compa(self,L,sweep,rc,u,rho1,rho2,see,cut,fit=False):## ENTROPY VS K(l)
        '''
        This function will extract by a linear fitting the central charge for each value of interaction V
        cut is to select number of block size you will take (see previous function calabrese_cardy_scal)
        inset is to plot in certain range of V the central charge, here between 5.4 and 6.5
        HERE different formula to extract see note project
        '''
        entrotest=list()
        list_v=list()
        infornate=self.db.find({"sweep":sweep,"L":L,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
        for l in infornate:
            list_v.append(l["v"])
        list_v=sorted(list_v)
        d={}
        c=list()
        SL=list()
        SL1=list()
        v=list()
        for x in range(0,len(list_v)):
            d["list_entropy_{0}".format(x)]=list()
            d["list_Lsize_{0}".format(x)]=list()
            d["plot_{0}".format(x)]=plt
        for i in range(0,len(list_v)):
                inf=self.db.find({"sweep":sweep,"u":u,"rc":rc,"L":L,"rho1":rho1,"rho2":rho2,"v":list_v[i]})
                infornata=self.get_max_list_float(inf, "L","n_sts")[0]
                for l in infornata["entangl_entr"]:
                        d["list_entropy_{0}".format(i)].append(l[2])
                        d["list_Lsize_{0}".format(i)].append(l[1])
                for p in range(2,len(infornata["entangl_entr"])):
                    if d["list_Lsize_{0}".format(i)][p]==L/2:
                        SL.append(d["list_entropy_{0}".format(i)][p])
                    if d["list_Lsize_{0}".format(i)][p]==L/2-1:
                        SL1.append(d["list_entropy_{0}".format(i)][p])
                v.append(list_v[i])
        for i in range(0,len(list_v)):
                c.append((6*L*L*(SL[i]-SL1[i]))/(numpy.pi*numpy.pi))
        fig=plt.figure()
        plt.plot(v,c,"bs--",markersize=10,linewidth=3)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.xlabel('V/t', fontsize=30)
        plt.ylabel('$c_3$', fontsize=30)
        if see==0:
            plt.show()
        else:
            plt.savefig('/home/thomas/CLL2Species/scanU20/pictures/C3_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'_L'+str(L)+'.pdf',bbox_inches='tight')


    def energy_convi(self,L,sweep,v,u,rc,rho1,rho2): ## ENERGY VS DMRG STATES
        '''
        parameters:
           sweep *
           L(single value)
           v
         * this is linked to the data precision(we can choose te better value to use studying the oters parameters)
        we want observe the GS energy behavior in function of DMRG states:
           1)Plot in terms of number of states
        '''
        if isinstance(L,list) and isinstance(L[0],int): # LIST
            self.energy_convi_list(L,sweep,v,u,rc,rho1,rho2)
        elif isinstance(L,int):                         # SINGLE
            self.energy_convi_single(L,sweep,v,u,rc,rho1,rho2)
        else:
            print "Error: the value of (L) inserted is not correct are accepted only INT or LIST of int "
####
    def gsLL(self,sweep,rc,v,u,rho1,rho2,fit=False): ## GSEnergy VS 1/L
        '''
        This function return a plot of the energy in terms of lenght of chain
        '''
        inf=self.db.find({"sweep":sweep,"v":v,"u":u,"rc":rc,'rho1':rho1,'rho2':rho2})## TODO CHECK
        infornata=self.get_max_list(inf,"L","n_sts")
        list_entropy=list()
        list_L=list()
        fig=plt.figure()
        ax=fig.add_subplot(1,1,1)
        cou=0
        for l in infornata:
            if cou > 0: ## with this counter i consider only the points over this value of cou.
                list_entropy.append(l["E"][0]/l['L'])
                list_L.append(1.0/(l["L"]))
            cou+=1
# fit function            I put a linear fit only because from the theory I know that have to be linear

        plt.plot(list_L,list_entropy,"ro--")
        if fit:
            plt.plot(xfine, fit(xfine, popt[0], popt[1]), 'b--')
            popt, pcov = curve_fit(fit, list_L, list_entropy)
            y=(list_entropy[1]-list_entropy[2])/6 + list_entropy[1]
# fit plot
            xfine = numpy.linspace(list_L[0], list_L[-1],100 )  # define values to plot the function for
            box1 = TextArea("a = "+str(popt[0])+" +/- "+str(math.sqrt(pcov[0,0]))+"\n"+"b = "+str(popt[1])+" +/- "+str(math.sqrt(pcov[1,1])),textprops=dict(color="k"))

            box = HPacker(children=[box1],align="center", pad=5, sep=5)

            anchored_box = AnchoredOffsetbox(loc=3,child=box, pad=0.,frameon=True,bbox_to_anchor=(0., 1.02),bbox_transform=ax.transAxes,borderpad=0.)

            ax.add_artist(anchored_box)

            fig.subplots_adjust(top=0.8)



        fig.suptitle('v='+str(v)+',u='+str(u)+r'$,r_{c}=$'+str(rc)+r"$,\rho_{1}=$"+str(rho1)+r"$,\rho_{2}=$"+str(rho2),fontsize=20)
        plt.xlabel('1/L', fontsize=18)
        plt.ylabel('E0/L', fontsize=16)
        plt.show()
    def gsL(self,sweep,rc,v,u,rho1,rho2,fit=False): ## GSEnergy VS 1/L^2
            '''
            This function return a plot of the energy in terms of lenght of chain
            '''
            inf=self.db.find({"sweep":sweep,"v":v,"u":u,"rc":rc,'rho1':rho1,'rho2':rho2})## TODO CHECK
            infornata=self.get_max_list(inf,"L","n_sts")
            list_entropy=list()
            list_L=list()
            fig=plt.figure()
            ax=fig.add_subplot(1,1,1)
            cou=0
            for l in infornata:
                if cou > 0:
                    list_entropy.append(l["E"][0]/l['L'])
                    list_L.append(1.0/(l["L"]*l["L"]))
                cou+=1
# fit function            I put a linear fit only because from the theory I know that have to be linear
            plt.plot(list_L,list_entropy,"ro--")
            if fit:
                popt, pcov = curve_fit(fit, list_L, list_entropy)
                vi=popt[0]*6/math.pi #here I compute the value of the fermi velocity
                y=(list_entropy[1]-list_entropy[2])/6 + list_entropy[1]
# fit plot
                xfine = numpy.linspace(list_L[0], list_L[-1],100 )  # define values to plot the function for
                plt.plot(xfine, fit(xfine, popt[0], popt[1]), 'b--')
                box1 = TextArea("a = "+str(popt[0])+" +/- "+str(math.sqrt(pcov[0,0]))+"\n"+"b = "+str(popt[1])+" +/- "+str(math.sqrt(pcov[1,1]))+"\n""v = "+str(vi),textprops=dict(color="k"))

                box = HPacker(children=[box1],align="center", pad=5, sep=5)

                anchored_box = AnchoredOffsetbox(loc=3,child=box, pad=0.,frameon=True,bbox_to_anchor=(0., 1.02),bbox_transform=ax.transAxes,borderpad=0.)

                ax.add_artist(anchored_box)

                fig.subplots_adjust(top=0.8)


            fig.suptitle('v='+str(v)+',u='+str(u)+r'$,r_{c}=$'+str(rc)+r"$,\rho_{1}=$"+str(rho1)+r"$,\rho_{2}=$"+str(rho2),fontsize=20)
            ax.tick_params(axis='both', which='major', labelsize=15)
            plt.xlabel(r'$1/L^2$', fontsize=25)
            plt.ylabel('E0/L', fontsize=25)
            plt.show()

    def gsLc(self,sweep,rc,v,u,rho1,rho2,fit=False): ## GSEnergy cluster VS 1/L^2
        '''
        This function return a plot of the energy in terms of lenght of chain
        '''
        inf=self.db.find({"sweep":sweep,"v":v,"u":u,"rc":rc,'rho1':rho1,'rho2':rho2})## TODO CHECK
        infornata=self.get_max_list(inf,"L","n_sts")
        list_entropy=list()
        list_L=list()
        fig=plt.figure()
        ax=fig.add_subplot(1,1,1)
        for l in infornata:
            if self.cluster(l["L"],l["rho1"],l["rc"]): ## I putted
                list_entropy.append(l["E"][0]/l['L'])
                list_L.append(1.0/(l["L"]*l["L"]))
# fit function            I put a linear fit only because from the theory I know that have to be linear
        print list_entropy
        print list_L
        plt.plot(list_L,list_entropy,"ro--")
        try:
            e2=list_entropy[2]
            l2=list_L[2]
        except:
            e2=None
            l2=None
        if fit and e2 and l2:
            popt, pcov = curve_fit(fit, list_L, list_entropy)
            y=(list_entropy[1]-list_entropy[2])/9 + list_entropy[1]
            #ax.annotate("a="+str(popt[0])+"+/-"+str(math.sqrt(pcov[0,0])),xy=(list_L[-1],list_entropy[1]))
            #ax.annotate("b="+str(popt[1])+"+/-"+str(math.sqrt(pcov[1,1])),xy=(list_L[-1],y))
# fit plot
            xfine = numpy.linspace(list_L[0], list_L[-1],100 )  # define values to plot the function for
            plt.plot(xfine, fit(xfine, popt[0], popt[1]), 'b--')
            fig=plt.figure(1, figsize=(3,3))
            ax = plt.subplot(111)

            box1 = TextArea("a = "+str(popt[0])+" +/- "+str(math.sqrt(pcov[0,0]))+"\n"+"b = "+str(popt[1])+" +/- "+str(math.sqrt(pcov[1,1])),textprops=dict(color="k"))

            box = HPacker(children=[box1],align="center", pad=5, sep=5)

            anchored_box = AnchoredOffsetbox(loc=3,child=box, pad=0.,frameon=True,bbox_to_anchor=(0., 1.02),bbox_transform=ax.transAxes,borderpad=0.)

            ax.add_artist(anchored_box)

            fig.subplots_adjust(top=0.8)
            fig.suptitle('v='+str(v)+'u='+str(u)+'$,r_{c}=$'+str(rc)+"$,\rho_{1}=$"+str(rho1)+"$,\rho_{2}=$"+str(rho2),fontsize=20)
            plt.xlabel('1/L^2', fontsize=18)
            plt.ylabel('E0/L', fontsize=16)
            plt.show()

    def energyL_list(self,sweep,rc,v,u,n_en,rho1,rho2): ## ENERGY VS L
        '''
        This function return a plot of the energy in terms of lenght of chain
        '''
        list_e=range(n_en)
        c=['r','k','b','c','g','m','y']
        marker=['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X','.']
        coppia=list()
        for colore in c:
            for m in marker:
                coppia.append([colore,m])

#inf=self.db.find({"sweep":sweep,"v":v,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})## TODO CHECK
        inf=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v,"u":u,"rc":rc})
        infornata=self.get_max_list(inf,"L","n_sts")
        list_entropy=list()
        list_L=list()
        for i in range(5):
            list_entropy.append(list())
            list_L.append(list())
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        for l in infornata:
            i=0
            for e in l["E"]:
                list_entropy[i].append(e)
                list_L[i].append(1.0/l["L"])
                i+=1
        for i in range(len(list_entropy)):
            list_entropy[i]=list_entropy[i][3:]
            list_L[i]=list_L[i][3:]
        for i in range(n_en):
            ax.scatter(list_L[i],list_entropy[i],c=coppia[i*19][0],marker=coppia[i][1],label="E" + str(i),s=70)
        fig.suptitle('v_'+str(v)+'rc'+str(rc),fontsize=20)
        plt.xlabel('1/L', fontsize=18)
        plt.ylabel('E', fontsize=16)
        plt.legend(loc='center right')
        plt.show()
    def energyscal(self,sweep,rc,v,u,n_en,rho1,rho2): ## ENERGY VS L TODO CHECK scaling
        '''
        This function return a plot of the energy in terms of lenght of chain
        '''
        list_en=range(n_en)
        c=['r','k','b','c','g','m','y']
        marker=['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X','.']
        coppia=list()
        for colore in c:
            for m in marker:
                coppia.append([colore,m])

        inf=self.db.find({"sweep":sweep,"v":v,"u":u,"rc":rc,"rho1":rho1,'rho2':rho2})## TODO CHECK
        infornata=self.get_max_list(inf,"L","n_sts")
        list_entropy=list()
        list_1ovL=list()
        for i in range(5):
            list_entropy.append(list())
            list_1ovL.append(list())
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        for l in infornata:
            i=0
            for e in reversed(l["E"]):
                list_entropy[i].append(e)
                list_1ovL[i].append(1./l["L"])
                i+=1
        for i in range(n_en):
            ax.scatter(list_1ovL[i],list_entropy[i],c=coppia[i*19][0],marker=coppia[i][1],label="E" + str(i))
        fig.suptitle('v_'+str(v)+'rc'+str(rc),fontsize=20)
        plt.xlabel('sites', fontsize=18)
        plt.ylabel('E', fontsize=16)
        plt.legend(loc='upper right')
        plt.show()

    def charge_gap(self,sweep,rc,v,u,rho1,rho2,rho11,rho22): # ENERGY DIFFERENCE IN DENSITY VS SITES
        ''' GS gap
        this function return the plot of values of energy gap in terms of L
        we choose rho1,rho2 = initial density
                  rho11,rho22 = final density
        '''
        inf1=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v,"u":u,"rc":rc}) ## TODO CHECK
        inf2=self.db.find({"rho1":rho11,'rho2':rho22,"sweep":sweep,"v":v,"u":u,"rc":rc}) ## TODO CHECK
        gap=list()
        lenght=list()
        infornata1=sorted(self.get_max_list(inf1,"L","n_sts"),key=lambda k :k["L"])
        infornata2=sorted(self.get_max_list(inf2,"L","n_sts"),key=lambda k :k["L"])

        if len(infornata1)<len(infornata2):
            k=0
            for l in infornata1:
                gap.append(math.fabs(l["E"][0]-infornata2[k]["E"][0])) # TODO CHECK
                lenght.append(l["L"])
                k+=1
        else:
            k=0
            for l in infornata2:
                gap.append(math.fabs(l["E0"]-infornata1[k]["E0"]))
                lenght.append(l["L"])
                k+=1
        fig=plt.figure()
        plt.plot(lenght,gap)
        fig.suptitle('v_'+str(v)+'rc'+str(rc)+'_rho1_'+str(rho2)+'rho2'+str(rho2),fontsize=20)
        plt.xlabel('sites', fontsize=18)
        plt.ylabel('gap', fontsize=16)
        plt.show()


    def single_particle_gap(self,sweep,Lmax,rc,v,u,rho1,rho2,rho11,rho22,rho111,rho222,see,fit=False): # ENERGY DIFFERENCE IN DENSITY VS SITES
        ''' GS gap
        this function return the plot of values of energy gap in terms of L
        we choose rho1,rho2 = initial density
                  rho11,rho22 = final density
        '''
        inf1=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v,"u":u,"rc":rc}) ## TODO CHECK
        inf2=self.db.find({"rho1":rho11,'rho2':rho22,"sweep":sweep,"v":v,"u":u,"rc":rc}) ## TODO CHECK
        inf3=self.db.find({"rho1":rho111,'rho2':rho222,"sweep":sweep,"v":v,"u":u,"rc":rc}) ## TODO CHECK
        gap=list()
        lenght=list()
        infornata1=sorted(self.get_max_list_float(inf1,"L","n_sts"),key=lambda k :k["L"])
        infornata2=sorted(self.get_max_list_float(inf2,"L","n_sts"),key=lambda k :k["L"])
        infornata3=sorted(self.get_max_list_float(inf3,"L","n_sts"),key=lambda k :k["L"])
        #if len(infornata1)<len(infornata2):

        for l in infornata1:
            for l2 in infornata2:
                for l3 in infornata3:
                    if l["L"]==l2["L"] and l2["L"]==l3["L"] and l["L"]==l3["L"]:
                        gap.append(math.fabs(1.0/2*(2*float(l["E"][0][0])-float(l2["E"][0][0])-float(l3["E"][0][0])))) # TODO CHECK
                        lenght.append(l["L"])
        '''
        k=0
        for l in infornata1:
            #if l["L"]==Lmax:
            print l["n_sts"]
            print infornata2[k]["n_sts"]
            print infornata3[k]["n_sts"]
            gap.append(math.fabs(1.0/2*(2*float(l["E"][0][0])-float(infornata2[k]["E"][0])-float(infornata3[k]["E"][0])))) # TODO CHECK
            lenght.append(l["L"])
            #k+=1
            #else:
                #print l["L"]

                print l["n_sts"]
                print infornata2[k]["n_sts"]
                print infornata3[k]["n_sts"]
                gap.append(math.fabs(1.0/2*(2*l["E"][0]-infornata2[k]["E"][0]-infornata3[k]["E"][0]))) # TODO CHECK
                lenght.append(l["L"])
                k+=1

        '''
	lenght=[1.0/y for y in lenght]
        fig=plt.figure()
        if fit:
            '''
            print  d["list_Lsize_0"]
            print  d["list_entropy_5"]
            print  d["list_entropy_0"][3]
            print  type(d["list_entropy_0"][3])
	    '''
            #d["list_entropy_0"][3]=float(d["list_entropy_0"][3])

            '''
	    print  d["list_entropy_8"]
            print  d["list_entropy_8"][0:3]
            print  type(d["list_entropy_0"][3])
            '''
            popt, pcov = curve_fit(fit, lenght,gap)
            xfine = numpy.linspace( 0, 0.06,100 )  # define values to plot the function for
            #y=(list_entropy[1]-list_entropy[2])/6 + list_entropy[1]
            plt.plot(xfine, fit(xfine, popt[0], popt[1],popt[2]), 'g-', linewidth=3)
# fit plot

            ax=fig.add_subplot(1,1,1)
            box1 = TextArea("a1 = "+str(popt[0])+" +/- "+str(math.sqrt(pcov[0,0]))+"\n"+"a2 = "+str(popt[1])+" +/- "+str(math.sqrt(pcov[1,1]))+"\n"+"a0 = "+str(popt[2])+" +/- "+str(math.sqrt(pcov[2,2])),textprops=dict(color="k",size=20))

            box = HPacker(children=[box1],align="center", pad=5, sep=5)

            anchored_box = AnchoredOffsetbox(loc=3,child=box, pad=0.,frameon=True,bbox_to_anchor=(0., 1.02),bbox_transform=ax.transAxes,borderpad=0.)

            ax.add_artist(anchored_box)

            fig.subplots_adjust(top=0.8)

        #print k
        #else:
         #   k=0
          #  for l in infornata2:
          #      gap.append(math.fabs(l["E0"]-infornata1[k]["E0"]))
            #     lenght.append(l["L"])
        #        k+=1
	#lenght=[log(y,10) for y in lenght]
	#gap=[log(y,10) for y in gap]
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.plot(lenght,gap,'rh',markersize=12)
        #fig.suptitle('v_'+str(v)+'rc'+str(rc)+'_rho1_'+str(rho2)+'rho2'+str(rho2),fontsize=25)
        plt.xlabel('1/L', fontsize=30)
        plt.ylabel('$\Delta_{sp}$', fontsize=30)
        if see==0:
            plt.show()
        else:
            plt.savefig('/home/thomas/CLL2Species/scanU20/pictures/SingleGap_V'+str(v)+'_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'.pdf',bbox_inches='tight')

    def single_particle_gap_compa(self,sweep,Lmax,rc,v1,v2,v3,v4,u,rho1,rho2,rho11,rho22,rho111,rho222,see,fit=False): # ENERGY DIFFERENCE IN DENSITY VS SITES

        inf1=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v1,"u":u,"rc":rc}) ## TODO CHECK
        inf2=self.db.find({"rho1":rho11,'rho2':rho22,"sweep":sweep,"v":v1,"u":u,"rc":rc}) ## TODO CHECK
        inf3=self.db.find({"rho1":rho111,'rho2':rho222,"sweep":sweep,"v":v1,"u":u,"rc":rc}) ## TODO CHECK
        gap=list()
        lenght=list()
        infornata1=sorted(self.get_max_list_float(inf1,"L","n_sts"),key=lambda k :k["L"])
        infornata2=sorted(self.get_max_list_float(inf2,"L","n_sts"),key=lambda k :k["L"])
        infornata3=sorted(self.get_max_list_float(inf3,"L","n_sts"),key=lambda k :k["L"])
        #if len(infornata1)<len(infornata2):

        for l in infornata1:
            for l2 in infornata2:
                for l3 in infornata3:
                    if l["L"]==l2["L"] and l2["L"]==l3["L"] and l["L"]==l3["L"]:
                        gap.append(math.fabs(1.0/2*(2*float(l["E"][0][0])-float(l2["E"][0][0])-float(l3["E"][0][0])))) # TODO CHECK
                        lenght.append(l["L"])

	lenght=[1.0/y for y in lenght]
        #fig=plt.figure()
        if fit:
            popt, pcov = curve_fit(fit, lenght,gap)
            xfine = numpy.linspace( 0, 0.06,100 )  # define values to plot the function for
            #y=(list_entropy[1]-list_entropy[2])/6 + list_entropy[1]
            #plt.plot(xfine, fit(xfine, popt[0], popt[1]), 'g-', linewidth=3)

########
        inf12=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v2,"u":u,"rc":rc}) ## TODO CHECK
        inf22=self.db.find({"rho1":rho11,'rho2':rho22,"sweep":sweep,"v":v2,"u":u,"rc":rc}) ## TODO CHECK
        inf32=self.db.find({"rho1":rho111,'rho2':rho222,"sweep":sweep,"v":v2,"u":u,"rc":rc}) ## TODO CHECK
        gap2=list()
        lenght2=list()
        infornata12=sorted(self.get_max_list_float(inf12,"L","n_sts"),key=lambda k :k["L"])
        infornata22=sorted(self.get_max_list_float(inf22,"L","n_sts"),key=lambda k :k["L"])
        infornata32=sorted(self.get_max_list_float(inf32,"L","n_sts"),key=lambda k :k["L"])
        #if len(infornata1)<len(infornata2):

        for ll in infornata12:
            for ll2 in infornata22:
                for ll3 in infornata32:
                    if ll["L"]==ll2["L"] and ll2["L"]==ll3["L"] and ll["L"]==ll3["L"]:
                        gap2.append(math.fabs(1.0/2*(2*float(ll["E"][0][0])-float(ll2["E"][0][0])-float(ll3["E"][0][0])))) # TODO CHECK
                        lenght2.append(ll["L"])
	lenght2=[1.0/y for y in lenght2]
        #fig=plt.figure()
        if fit:
            popt2, pcov2 = curve_fit(fit, lenght2,gap2)
            xfine = numpy.linspace( 0, 0.06,100 )  # define values to plot the function for
            #y=(list_entropy[1]-list_entropy[2])/6 + list_entropy[1]

########
        inf13=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v3,"u":u,"rc":rc}) ## TODO CHECK
        inf23=self.db.find({"rho1":rho11,'rho2':rho22,"sweep":sweep,"v":v3,"u":u,"rc":rc}) ## TODO CHECK
        inf33=self.db.find({"rho1":rho111,'rho2':rho222,"sweep":sweep,"v":v3,"u":u,"rc":rc}) ## TODO CHECK
        gap3=list()
        lenght3=list()
        infornata13=sorted(self.get_max_list_float(inf13,"L","n_sts"),key=lambda k :k["L"])
        infornata23=sorted(self.get_max_list_float(inf23,"L","n_sts"),key=lambda k :k["L"])
        infornata33=sorted(self.get_max_list_float(inf33,"L","n_sts"),key=lambda k :k["L"])
        #if len(infornata1)<len(infornata2):

        for p in infornata13:
            for p2 in infornata23:
                for p3 in infornata33:
                    if p["L"]==p2["L"] and p2["L"]==p3["L"] and p["L"]==p3["L"]:
                        gap3.append(math.fabs(1.0/2*(2*float(p["E"][0][0])-float(p2["E"][0][0])-float(p3["E"][0][0])))) # TODO CHECK
                        lenght3.append(p["L"])

	lenght3=[1.0/y for y in lenght3]
        fig=plt.figure()
        if fit:
            popt3, pcov3 = curve_fit(fit, lenght3,gap3)
            xfine = numpy.linspace( 0, 0.06,100 )  # define values to plot the function for
            #y=(list_entropy[1]-list_entropy[2])/6 + list_entropy[1]

########
        inf14=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v4,"u":u,"rc":rc}) ## TODO CHECK
        inf24=self.db.find({"rho1":rho11,'rho2':rho22,"sweep":sweep,"v":v4,"u":u,"rc":rc}) ## TODO CHECK
        inf34=self.db.find({"rho1":rho111,'rho2':rho222,"sweep":sweep,"v":v4,"u":u,"rc":rc}) ## TODO CHECK
        gap4=list()
        lenght4=list()
        infornata14=sorted(self.get_max_list_float(inf14,"L","n_sts"),key=lambda k :k["L"])
        infornata24=sorted(self.get_max_list_float(inf24,"L","n_sts"),key=lambda k :k["L"])
        infornata34=sorted(self.get_max_list_float(inf34,"L","n_sts"),key=lambda k :k["L"])
        #if len(infornata1)<len(infornata2):

        for p in infornata14:
            for p2 in infornata24:
                for p3 in infornata34:
                    if p["L"]==p2["L"] and p2["L"]==p3["L"] and p["L"]==p3["L"]:
                        gap4.append(math.fabs(1.0/2*(2*float(p["E"][0][0])-float(p2["E"][0][0])-float(p3["E"][0][0])))) # TODO CHECK
                        lenght4.append(p["L"])
        lenght4=[1.0/y for y in lenght4]
        fig=plt.figure()
        if fit:
            popt4, pcov4 = curve_fit(fit, lenght4,gap4)
            xfine = numpy.linspace( 0, 0.06,100 )  # define values to plot the function for
            #y=(list_entropy[1]-list_entropy[2])/6 + list_entropy[1]

        plt.plot(lenght,gap,'gh',markersize=12,label='V='+str(v1))
        plt.plot(xfine, fit(xfine, popt[0], popt[1],popt[2]), 'g-', linewidth=3)
        plt.plot(lenght2,gap2,'k^',markersize=12,label='V='+str(v2))
        plt.plot(xfine, fit(xfine, popt2[0], popt2[1],popt2[2]), 'k-', linewidth=3)
        plt.plot(lenght3,gap3,'ro',markersize=12,label='V='+str(v3))
        plt.plot(xfine, fit(xfine, popt3[0], popt3[1],popt3[2]), 'r-', linewidth=3)
        plt.plot(lenght4,gap4,'bp',markersize=12,label='V='+str(v4))
        plt.plot(xfine, fit(xfine, popt4[0], popt4[1],popt4[2]), 'b-', linewidth=3)
        plt.legend(fontsize=15)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        #fig.suptitle('v_'+str(v)+'rc'+str(rc)+'_rho1_'+str(rho2)+'rho2'+str(rho2),fontsize=25)
        plt.xlabel('1/L', fontsize=30)
        plt.ylabel('$\Delta_{sp}$', fontsize=30)
        if see==0:
            plt.show()
        else:
            plt.savefig('/home/thomas/CLL2Species/scanU20/pictures/SingleGap_compa_V1'+str(v1)+'_V4'+str(v4)+'_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'.pdf',bbox_inches='tight')

    def SingleGap_vs_V(self,L,sweep,rc,u,rho1,rho2,rho3,rho4,rho5,rho6,see):## ENTROPY VS K(l)

        list_v=list()
        infornate=self.db.find({"sweep":sweep,"L":L,"u":u,"rc":rc,"rho1":rho3,"rho2":rho4})
        #infornate=self.get_max_list_float(inf,"L","n_sts")
        for l in infornate:
            list_v.append(l["v"])
        list_v=sorted(list_v)
        #print list_v
        d={}
        gap=list()
        v=list()
        for i in range(0,len(list_v)):
                inf1=self.db.find({"sweep":sweep,"L":L,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2,"v":list_v[i]})
                inf2=self.db.find({"sweep":sweep,"L":L,"u":u,"rc":rc,"rho1":rho3,"rho2":rho4,"v":list_v[i]})
                inf3=self.db.find({"sweep":sweep,"L":L,"u":u,"rc":rc,"rho1":rho5,"rho2":rho6,"v":list_v[i]})
                #infornata1=self.get_max_list_float(inf1, "L","n_sts")[0]
                #infornata2=self.get_max_list_float(inf2, "L","n_sts")[0]
                #infornata3=self.get_max_list_float(inf3, "L","n_sts")[0]
                for l1 in inf1:
                        for l2 in inf2:
                                for l3 in inf3:
                                        gap.append(math.fabs(1.0/2*(2*float(l1["E"][0][0])-float(l2["E"][0][0])-float(l3["E"][0][0])))) # TODO CHECK

                v.append(list_v[i])
        fig=plt.figure()
        plt.plot(v,gap,"rs--",markersize=10,markeredgecolor='k',linewidth=3)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.xlabel('V/t', fontsize=30)
        plt.ylabel('$\Delta_{singlegap}$', fontsize=30)
        if see==0:
            plt.show()
        else:
            plt.savefig('/home/thomas/CLL2Species/scanU20/pictures/SingleGap_vs_V_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'_L'+str(L)+'.pdf',bbox_inches='tight')

    def Singlegap_vs_V_compaU(self,L,sweep,rc,u,v,rho1,rho2,rho3,rho4,rho5,rho6,see):## ENTROPY VS K(l)

        list_u=list()
        infornate=self.db.find({"sweep":sweep,"L":L,"v":v,"rc":rc,"rho1":rho3,"rho2":rho4})
        #infornate=self.get_max_list_float(inf,"L","n_sts")
        for l in infornate:
                list_u.append(l["u"])
                list_u=sorted(list_u)
        #print list_u
        label=list()
        label.append(list_u[0])
        label.append(list_u[1])

        d={}

        for x in range(0,len(list_u)):
            d["list_sp_{0}".format(x)]=list()
            d["list_L_{0}".format(x)]=list()
            d["plot_{0}".format(x)]=plt

        for j in range(0,len(list_u)):
            list_v=list()
            infornate=self.db.find({"sweep":sweep,"L":L,"u":list_u[j],"rc":rc,"rho1":rho3,"rho2":rho4})
            #infornate=self.get_max_list_float(inf,"L","n_sts")
            for l in infornate:
                list_v.append(l["v"])
            list_v=sorted(list_v)
            #print list_v
            gap=list()
            v=list()
            for i in range(0,len(list_v)):
                    inf1=self.db.find({"sweep":sweep,"L":L,"u":list_u[j],"rc":rc,"rho1":rho1,"rho2":rho2,"v":list_v[i]})
                    inf2=self.db.find({"sweep":sweep,"L":L,"u":list_u[j],"rc":rc,"rho1":rho3,"rho2":rho4,"v":list_v[i]})
                    inf3=self.db.find({"sweep":sweep,"L":L,"u":list_u[j],"rc":rc,"rho1":rho5,"rho2":rho6,"v":list_v[i]})
                    #infornata1=self.get_max_list_float(inf1, "L","n_sts")[0]
                    #infornata2=self.get_max_list_float(inf2, "L","n_sts")[0]
                    #infornata3=self.get_max_list_float(inf3, "L","n_sts")[0]
                    for l1 in inf1:
                            for l2 in inf2:
                                for l3 in inf3:
                                                gap.append(math.fabs(1.0/2*(2*float(l1["E"][0][0])-float(l2["E"][0][0])-float(l3["E"][0][0]))))
                                                #gap.append(math.fabs(float(l1["E"][0][0])-float(l2["E"][0][0]))) # TODO CHECK
                    v.append(list_v[i])
                    d["list_L_{0}".format(j)]=v
                    d["list_sp_{0}".format(j)]=gap
        mark=['s','8','D','o','x','*','P','^']
        linetype=['-','--',':','-.','-']
        fig=plt.figure()
        for i in range(0,len(list_u)):
                plt.plot(d["list_L_{0}".format(i)],d["list_sp_{0}".format(i)],linestyle=linetype[i],marker=mark[i],markersize=10,markeredgecolor='k',linewidth=3,label='U='+str(label[i]))
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.xlabel('V/t', fontsize=30)
        plt.ylabel('$\Delta_{chargegap}$', fontsize=30)
        plt.legend(fontsize=25)
        if see==0:
            plt.show()
        else:
            plt.savefig(os.path.join(path,'Singlegap_vs_V_compaU_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'_L'+str(L)+'.png'),bbox_inches='tight')

    def SpinGap_vs_V(self,L,sweep,rc,u,rho1,rho2,rho3,rho4,see):## ENTROPY VS K(l)

        list_v=list()
        infornate=self.db.find({"sweep":sweep,"L":L,"u":u,"rc":rc,"rho1":rho3,"rho2":rho4})
        #infornate=self.get_max_list_float(inf,"L","n_sts")
        for l in infornate:
            list_v.append(l["v"])
        list_v=sorted(list_v)
        del list_v[1:3] #problem for L=30 data missing
        #print list_v
        d={}
        gap=list()
        v=list()
        for i in range(0,len(list_v)):
                inf1=self.db.find({"sweep":sweep,"L":L,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2,"v":list_v[i]})
                inf2=self.db.find({"sweep":sweep,"L":L,"u":u,"rc":rc,"rho1":rho3,"rho2":rho4,"v":list_v[i]})
                #infornata1=self.get_max_list_float(inf1, "L","n_sts")[0]
                #infornata2=self.get_max_list_float(inf2, "L","n_sts")[0]
                #infornata3=self.get_max_list_float(inf3, "L","n_sts")[0]
                for l1 in inf1:
                        for l2 in inf2:
                                    gap.append(math.fabs(float(l1["E"][0][0])-float(l2["E"][0][0]))) # TODO CHECK

                v.append(list_v[i])
        fig=plt.figure()
        plt.plot(v,gap,"b8--",markersize=10,linewidth=3)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.xlabel('V/t', fontsize=30)
        plt.ylabel('$\Delta_{spingap}$', fontsize=30)
        if see==0:
            plt.show()
        else:
            plt.savefig(os.path.join(path,'SpinGap_vs_V_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'_L'+str(L)+'.pdf'),bbox_inches='tight')

    def SpinGap_vs_V_compa(self,L,sweep,rc,u,v,rho1,rho2,rho3,rho4,see):## ENTROPY VS K(l)

        list_u=list()
        infornate=self.db.find({"sweep":sweep,"L":L,"v":v,"rc":rc,"rho1":rho3,"rho2":rho4})
        #infornate=self.get_max_list_float(inf,"L","n_sts")
        for l in infornate:
                list_u.append(l["u"])
                list_u=sorted(list_u)
        #print list_u
        label=list()
        label.append(list_u[0])
        label.append(list_u[1])

        d={}

        for x in range(0,len(list_u)):
            d["list_sp_{0}".format(x)]=list()
            d["list_L_{0}".format(x)]=list()
            d["plot_{0}".format(x)]=plt

        for j in range(0,len(list_u)):
            list_v=list()
            infornate=self.db.find({"sweep":sweep,"L":L,"u":list_u[j],"rc":rc,"rho1":rho3,"rho2":rho4})
            #infornate=self.get_max_list_float(inf,"L","n_sts")
            for l in infornate:
                list_v.append(l["v"])
            list_v=sorted(list_v)
            #print list_v
            gap=list()
            v=list()
            for i in range(0,len(list_v)):
                    inf1=self.db.find({"sweep":sweep,"L":L,"u":list_u[j],"rc":rc,"rho1":rho1,"rho2":rho2,"v":list_v[i]})
                    inf2=self.db.find({"sweep":sweep,"L":L,"u":list_u[j],"rc":rc,"rho1":rho3,"rho2":rho4,"v":list_v[i]})
                    #infornata1=self.get_max_list_float(inf1, "L","n_sts")[0]
                    #infornata2=self.get_max_list_float(inf2, "L","n_sts")[0]
                    #infornata3=self.get_max_list_float(inf3, "L","n_sts")[0]
                    for l1 in inf1:
                            for l2 in inf2:
                                        gap.append(math.fabs(float(l1["E"][0][0])-float(l2["E"][0][0]))) # TODO CHECK
                    v.append(list_v[i])
                    d["list_L_{0}".format(j)]=v
                    d["list_sp_{0}".format(j)]=gap
        mark=['s','8','D','o','x','*','P','^']
        linetype=['-','--',':','-.','-']
        fig=plt.figure()
        for i in range(0,len(list_u)):
                plt.plot(d["list_L_{0}".format(i)],d["list_sp_{0}".format(i)],linestyle=linetype[i],marker=mark[i],markersize=10,linewidth=3,label='U='+str(label[i]))
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.xlabel('V/t', fontsize=30)
        plt.ylabel('$\Delta_{spingap}$', fontsize=30)
        plt.legend(fontsize=25)
        if see==0:
            plt.show()
        else:
            plt.savefig(os.path.join(path,'SpinGap_vs_V_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'_L'+str(L)+'.png'),bbox_inches='tight')


    def Oss(self,sweep,rc,v,u,rho1,rho2,L,see,observable1='ss_ab'):
        inf1=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v,"u":u,"rc":rc,"L":L})
        infornata1=self.get_max_list_float(inf1,"L","n_sts")
        ssab=list()
        size=list()
        for i in range(0,L):
            for l in infornata1:
                   ssab.append(l[observable1][i])
                   size.append(i+1)
        fig=plt.figure()
        ax = fig.add_subplot(111)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.xlabel('j', fontsize=30)
        plt.ylabel('$O_{SS}$', fontsize=30)
        plt.plot(size,ssab,'bo--')
        if see==0:
                plt.show()
        else:
                plt.savefig('/home/thomas/CLL2Species/scanU20/pictures/Oss_V'+str(v)+'_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'_L'+str(L)+'.pdf',bbox_inches='tight')

    def Tss0(self,sweep,rc,v,u,rho1,rho2,L,see,observable1='Tss_0a',observable2='Tss_0b',observable3='Tss_0aa',observable4='Tss_0bb'):
        inf1=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v,"u":u,"rc":rc,"L":L})
        infornata1=self.get_max_list_float(inf1,"L","n_sts")
        Tss0a=list()
        Tss0b=list()
        Tss0=list()
        size=list()
        for i in range(0,L-1):
            for l in infornata1:
                   Tss0.append( l[observable1][i] + l[observable2][i] + l[observable3][i] + l[observable4][i])
                   size.append(i+1)

        fig=plt.figure()
        ax = fig.add_subplot(111)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
        plt.xlabel('j', fontsize=30)
        plt.ylabel('$O_{TS_0}$', fontsize=30)
        plt.plot(size,Tss0,'bo--')
        if see==0:
               plt.show()
        else:
                plt.savefig('/home/thomas/CLL2Species/scanU20/pictures/Tss0_V'+str(v)+'_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'_L'+str(L)+'.pdf',bbox_inches='tight')

    def BCDW(self,sweep,rc,v,u,rho1,rho2,L,see,observable1='bcdw_a',observable2='bcdw_b',observable3='Ibcdw_b',observable4='Ibcdw_b'):
        inf1=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v,"u":u,"rc":rc,"L":L})
        infornata1=self.get_max_list_float(inf1,"L","n_sts")
        BCDWa=list()
        BCDWb=list()
        IBCDWa=list()
        IBCDWb=list()
        BCDWz=list()
        size=list()
        for i in range(0,L-1):
            for l in infornata1:
                   BCDWa.append(l[observable1][i])
                   BCDWb.append(l[observable2][i])
                   IBCDWa.append(l[observable3][i])
                   IBCDWb.append(l[observable4][i])
                   size.append(i+1)
        for i in range(0,L-1):
            BCDWz.append(BCDWa[i]+BCDWb[i]+IBCDWa[i]+IBCDWb[i])

        fig=plt.figure()
        ax = fig.add_subplot(111)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.xlabel('j', fontsize=30)
        plt.ylabel('$O_{BCDWz}$', fontsize=30)
        plt.plot(size,BCDWz,'bo--')
        if see==0:
               plt.show()
        else:
                plt.savefig('/home/thomas/CLL2Species/scanU20/pictures/BCDW_V'+str(v)+'_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'_L'+str(L)+'.pdf',bbox_inches='tight')

    def SCDW(self,sweep,rc,v,u,rho1,rho2,L,see,observable1='bcdw_a',observable2='bcdw_b',observable3='Ibcdw_b',observable4='Ibcdw_b'):
        inf1=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v,"u":u,"rc":rc,"L":L})
        infornata1=self.get_max_list_float(inf1,"L","n_sts")
        BCDWa=list()
        BCDWb=list()
        IBCDWa=list()
        IBCDWb=list()
        SCDWz=list()
        size=list()
        for i in range(0,L-1):
            for l in infornata1:
                   BCDWa.append(l[observable1][i])
                   BCDWb.append(l[observable2][i])
                   IBCDWa.append(l[observable3][i])
                   IBCDWb.append(l[observable4][i])
                   size.append(i+1)
        for i in range(0,L-1):
            SCDWz.append((BCDWa[i]-BCDWb[i])+(IBCDWa[i]-IBCDWb[i]))

        fig=plt.figure()
        ax = fig.add_subplot(111)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.xlabel('j', fontsize=30)
        plt.ylabel('$O_{SCDWz}$', fontsize=30)
        plt.plot(size,SCDWz,'bo--')
        if see==0:
                plt.show()
        else:
                plt.savefig('/home/thomas/CLL2Species/scanU20/pictures/SCDW_V'+str(v)+'_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'_L'+str(L)+'.pdf',bbox_inches='tight')

    def spin_gap(self,sweep,Lmax,rc,v,u,rho1,rho2,rho3,rho4,see,fit=False): # ENERGY DIFFERENCE IN DENSITY VS SITES


        list_L=list()
        infornate=self.db.find({"v":v,"u":u,"rc":rc,"rho1":rho3,"rho2":rho4})
        for l in infornate:
                list_L.append(l["L"])

        list_L=sorted(list_L) #sorte for the value of interaction
        print list_L

        d={}
        gap=list()
        size=list()
        for x in range(0,len(list_L)):
            d["list_gap_{0}".format(x)]=list()
            d["list_Lsize_{0}".format(x)]=list()
            d["plot_{0}".format(x)]=plt
        for i in range(0,len(list_L)):
                    inf=self.db.find({"u":u,"rc":rc,"L":list_L[i],"rho1":rho1,"rho2":rho2,"v":v})
                    inf2=self.db.find({"u":u,"rc":rc,"L":list_L[i],"rho1":rho3,"rho2":rho4,"v":v})
                    infornata=self.get_max_list_float(inf, "L","n_sts")[0]
                    infornata2=self.get_max_list_float(inf2, "L","n_sts")[0]

                    gap.append(math.fabs(float(infornata["E"][0][0])-float(infornata2["E"][0][0])))
                    size.append(list_L[i])
        size=[1.0/y for y in size]
        fig=plt.figure()
        if fit:
            popt, pcov = curve_fit(fit,size, gap, maxfev=1000)
            xfine = numpy.linspace( 0, 0.06,100 )  # define values to plot the function for
            #y=(list_entropy[1]-list_entropy[2])/6 + list_entropy[1]
            plt.plot(xfine, fit(xfine, popt[0], popt[1]), 'k-', linewidth=3)

            ax=fig.add_subplot(1,1,1)
            box2 = TextArea("a1 = "+str(popt[0])+" +/- "+str(math.sqrt(pcov[0,0]))+"\n"+"a2 = "+str(popt[1])+" +/- "+str(math.sqrt(pcov[1,1])),textprops=dict(color="k",size=20))
            box2 = HPacker(children=[box2],align="center", pad=5, sep=5)
            #+"\n"+"a0 = "+str(popt[2])+" +/- "+str(math.sqrt(pcov[2,2])
            anchored_box2 = AnchoredOffsetbox(loc=3,child=box2, pad=0.,frameon=True,bbox_to_anchor=(0., 1.02),bbox_transform=ax.transAxes,borderpad=0.)

            ax.add_artist(anchored_box2)

            fig.subplots_adjust(top=0.8)


        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.plot(size,gap,'bs--',markersize=12)
        #fig.suptitle('v_'+str(v)+'rc'+str(rc)+'_rho1_'+str(rho2)+'rho2'+str(rho2),fontsize=20)
        plt.xlabel('1/L', fontsize=30)
        plt.ylabel('$\Delta_{SpinGap}$', fontsize=30)
        if see==0:
            plt.show()
        else:
            plt.savefig('/home/thomas/CLL2Species/scanU20/pictures/SpinGap_V'+str(v)+'_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'.pdf',bbox_inches='tight')

    def spin_gap_compa2(self,rc,list_v,u,rho1,rho2,rho3,rho4,see,fit=False): # ENERGY DIFFERENCE IN DENSITY VS SITES

        d={}
        for x in range(0,len(list_v)):
                        d["list_gap_{0}".format(x)]=list()
                        d["list_Lsize_{0}".format(x)]=list()
                        d["list_popt_{0}".format(x)]=list()
                        d["list_pcov_{0}".format(x)]=list()
                        d["plot_{0}".format(x)]=plt
        for j in range(0,len(list_v)):
            list_L=list()
            infornate=self.db.find({"v":list_v[j],"u":u,"rc":rc,"rho1":rho3,"rho2":rho4})
            for l in infornate:
                    list_L.append(l["L"])

            list_L=sorted(list_L) #sorte for the value of interaction
            #print list_L
            gap=list()
            size=list()
            for i in range(0,len(list_L)):
                        inf=self.db.find({"u":u,"rc":rc,"L":list_L[i],"rho1":rho1,"rho2":rho2,"v":list_v[j]})
                        inf2=self.db.find({"u":u,"rc":rc,"L":list_L[i],"rho1":rho3,"rho2":rho4,"v":list_v[j]})
                        infornata=self.get_max_list_float(inf, "L","n_sts")[0]
                        infornata2=self.get_max_list_float(inf2, "L","n_sts")[0]

                        gap.append(math.fabs(float(infornata["E"][0][0])-float(infornata2["E"][0][0])))
                        size.append(list_L[i])
            size=[1.0/y for y in size]
            print size
            d["list_gap_{0}".format(j)]=gap
            d["list_Lsize_{0}".format(j)]=size
            print d["list_Lsize_{0}".format(j)]
            fig=plt.figure()
            if fit:
                d["list_popt_{0}".format(j)],d["list_pcov_{0}".format(j)] = curve_fit(fit,size, gap, maxfev=1000)
                xfine = numpy.linspace( 0, 0.06,100 )  # define values to plot the function for
                #y=(list_entropy[1]-list_entropy[2])/6 + list_entropy[1]
                #plt.plot(xfine, fit(xfine, popt[0], popt[1]), 'k-', linewidth=3)
                '''
                ax=fig.add_subplot(1,1,1)
                box2 = TextArea("a1 = "+str(popt[0])+" +/- "+str(math.sqrt(pcov[0,0]))+"\n"+"a2 = "+str(popt[1])+" +/- "+str(math.sqrt(pcov[1,1])),textprops=dict(color="k",size=20))
                box2 = HPacker(children=[box2],align="center", pad=5, sep=5)
                #+"\n"+"a0 = "+str(popt[2])+" +/- "+str(math.sqrt(pcov[2,2])
                anchored_box2 = AnchoredOffsetbox(loc=3,child=box2, pad=0.,frameon=True,bbox_to_anchor=(0., 1.02),bbox_transform=ax.transAxes,borderpad=0.)

                ax.add_artist(anchored_box2)

                fig.subplots_adjust(top=0.8)
                '''
        xfine = numpy.linspace( 0, 0.06,100 )
        mark=['s','8','D','o','x','*','P','^']
        linetype=['-','--',':','-.','-']
        col=['r','k','b','g','c','y']
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        for j in range(0,len(list_v)):
                plt.plot(d["list_Lsize_{0}".format(j)],d["list_gap_{0}".format(j)],marker=mark[j],color=col[j],markersize=12,label='V'+str(list_v[j]))
                plt.plot(xfine,fit(xfine,d["list_popt_{0}".format(j)][0],d["list_popt_{0}".format(j)][1]), '-',color=col[j],linewidth=3)
                #plt.plot(xfine,fit(xfine,d["list_popt_{0}".format(j)][0]), '-',color=col[j],linewidth=3)
        #fig.suptitle('v_'+str(v)+'rc'+str(rc)+'_rho1_'+str(rho2)+'rho2'+str(rho2),fontsize=20)
        plt.xlabel('1/L', fontsize=30)
        plt.ylabel('$\Delta_{SpinGap}$', fontsize=30)
        plt.legend(fontsize=25)
        if see==0:
            plt.show()
        else:
            plt.savefig('/home/thomas/CLL2Species/scanU20/pictures/SpinGap_V'+str(v)+'_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'.pdf',bbox_inches='tight')


    def singlegap_compa2(self,rc,list_v,u,rho1,rho2,rho3,rho4,rho5,rho6,see,fit=False): # ENERGY DIFFERENCE IN DENSITY VS SITES

        d={}
        for x in range(0,len(list_v)):
                        d["list_gap_{0}".format(x)]=list()
                        d["list_Lsize_{0}".format(x)]=list()
                        d["list_popt_{0}".format(x)]=list()
                        d["list_pcov_{0}".format(x)]=list()
                        d["plot_{0}".format(x)]=plt
        for j in range(0,len(list_v)):
            list_L=list()
            infornate=self.db.find({"v":list_v[j],"u":u,"rc":rc,"rho1":rho3,"rho2":rho4})
            for l in infornate:
                    list_L.append(l["L"])

            list_L=sorted(list_L) #sorte for the value of interaction
            #print list_L
            gap=list()
            size=list()
            for i in range(0,len(list_L)):
                        inf=self.db.find({"u":u,"rc":rc,"L":list_L[i],"rho1":rho1,"rho2":rho2,"v":list_v[j]})
                        inf2=self.db.find({"u":u,"rc":rc,"L":list_L[i],"rho1":rho3,"rho2":rho4,"v":list_v[j]})
                        inf3=self.db.find({"u":u,"rc":rc,"L":list_L[i],"rho1":rho5,"rho2":rho6,"v":list_v[j]})
                        infornata=self.get_max_list_float(inf, "L","n_sts")[0]
                        infornata2=self.get_max_list_float(inf2, "L","n_sts")[0]
                        infornata3=self.get_max_list_float(inf3, "L","n_sts")[0]

                        gap.append(math.fabs(0.5*(2*float(infornata["E"][0][0])-float(infornata2["E"][0][0])-float(infornata3["E"][0][0]))))
                        size.append(list_L[i])
            size=[1.0/y for y in size]
            print size
            d["list_gap_{0}".format(j)]=gap
            d["list_Lsize_{0}".format(j)]=size
            print d["list_Lsize_{0}".format(j)]
            fig=plt.figure()
            if fit:
                d["list_popt_{0}".format(j)],d["list_pcov_{0}".format(j)] = curve_fit(fit,size, gap, maxfev=1000)
                xfine = numpy.linspace( 0, 0.06,100 )  # define values to plot the function for
                #y=(list_entropy[1]-list_entropy[2])/6 + list_entropy[1]
                #plt.plot(xfine, fit(xfine, popt[0], popt[1]), 'k-', linewidth=3)
                '''
                ax=fig.add_subplot(1,1,1)
                box2 = TextArea("a1 = "+str(popt[0])+" +/- "+str(math.sqrt(pcov[0,0]))+"\n"+"a2 = "+str(popt[1])+" +/- "+str(math.sqrt(pcov[1,1])),textprops=dict(color="k",size=20))
                box2 = HPacker(children=[box2],align="center", pad=5, sep=5)
                #+"\n"+"a0 = "+str(popt[2])+" +/- "+str(math.sqrt(pcov[2,2])
                anchored_box2 = AnchoredOffsetbox(loc=3,child=box2, pad=0.,frameon=True,bbox_to_anchor=(0., 1.02),bbox_transform=ax.transAxes,borderpad=0.)

                ax.add_artist(anchored_box2)

                fig.subplots_adjust(top=0.8)
                '''
        xfine = numpy.linspace( 0, 0.06,100 )
        mark=['s','8','D','o','x','*','P','^']
        linetype=['-','--',':','-.','-']
        col=['r','k','b','g','c','y']
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        for j in range(0,len(list_v)):
                plt.plot(d["list_Lsize_{0}".format(j)],d["list_gap_{0}".format(j)],marker=mark[j],color=col[j],markersize=12,label='V'+str(list_v[j]))
                plt.plot(xfine,fit(xfine,d["list_popt_{0}".format(j)][0],d["list_popt_{0}".format(j)][1]), '-',color=col[j],linewidth=3)
        #fig.suptitle('v_'+str(v)+'rc'+str(rc)+'_rho1_'+str(rho2)+'rho2'+str(rho2),fontsize=20)
        plt.xlabel('1/L', fontsize=30)
        plt.ylabel('$\Delta_{SpinGap}$', fontsize=30)
        plt.legend(fontsize=25)
        if see==0:
            plt.show()
        else:
            plt.savefig('/home/thomas/CLL2Species/scanU20/pictures/SingleGap_V'+str(v)+'_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'.pdf',bbox_inches='tight')

    def spin_gap_compa(self,sweep,Lmax,rc,v1,v2,v3,u,rho1,rho2,rho11,rho22,see,fit=False): # ENERGY DIFFERENCE IN DENSITY VS SITES

        inf1=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v1,"u":u,"rc":rc}) ## TODO CHECK
        inf2=self.db.find({"rho1":rho11,'rho2':rho22,"sweep":sweep,"v":v1,"u":u,"rc":rc}) ## TODO CHECK
        gap=list()
        lenght=list()
        infornata1=sorted(self.get_max_list_float(inf1,"L","n_sts"),key=lambda k :k["L"])
        infornata2=sorted(self.get_max_list_float(inf2,"L","n_sts"),key=lambda k :k["L"])

        for l in infornata1:
            for l2 in infornata2:
                    if l["L"]==l2["L"]:
                        gap.append(math.fabs(float(l["E"][0][0])-float(l2["E"][0][0]))) # TODO CHECK
                        lenght.append(l["L"])

        lenght=[1.0/y for y in lenght]
        #fig=plt.figure()
        if fit:
            popt, pcov = curve_fit(fit, lenght,gap)
            xfine = numpy.linspace( 0, 0.06,100 )  # define values to plot the function for
            #y=(list_entropy[1]-list_entropy[2])/6 + list_entropy[1]
            #plt.plot(xfine, fit(xfine, popt[0], popt[1]), 'g-', linewidth=3)
# fit plot
            '''
            ax=fig.add_subplot(1,1,1)
            box1 = TextArea("a = "+str(popt[0])+" +/- "+str(math.sqrt(pcov[0,0]))+"\n"+"b = "+str(popt[1])+" +/- "+str(math.sqrt(pcov[1,1])),textprops=dict(color="k",size=12))

            box = HPacker(children=[box1],align="center", pad=5, sep=5)

            anchored_box = AnchoredOffsetbox(loc=3,child=box, pad=0.,frameon=True,bbox_to_anchor=(0., 1.02),bbox_transform=ax.transAxes,borderpad=0.)

            ax.add_artist(anchored_box)

            fig.subplots_adjust(top=0.8)
            '''
##################
        inf12=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v2,"u":u,"rc":rc}) ## TODO CHECK
        inf22=self.db.find({"rho1":rho11,'rho2':rho22,"sweep":sweep,"v":v2,"u":u,"rc":rc}) ## TODO CHECK
        gap2=list()
        lenght2=list()
        infornata12=sorted(self.get_max_list_float(inf12,"L","n_sts"),key=lambda k :k["L"])
        infornata22=sorted(self.get_max_list_float(inf22,"L","n_sts"),key=lambda k :k["L"])

        for ll in infornata12:
            for ll2 in infornata22:
                    if ll["L"]==ll2["L"]:
                        gap2.append(math.fabs(float(ll["E"][0][0])-float(ll2["E"][0][0]))) # TODO CHECK
                        lenght2.append(ll["L"])
        lenght2=[1.0/y for y in lenght2]
        if fit:
            popt2, pcov2 = curve_fit(fit, lenght2,gap2)
            xfine = numpy.linspace( 0, 0.06,100 )  # define values to plot the function for
############

        inf13=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v3,"u":u,"rc":rc}) ## TODO CHECK
        inf23=self.db.find({"rho1":rho11,'rho2':rho22,"sweep":sweep,"v":v3,"u":u,"rc":rc}) ## TODO CHECK
        gap3=list()
        lenght3=list()
        infornata13=sorted(self.get_max_list_float(inf13,"L","n_sts"),key=lambda k :k["L"])
        infornata23=sorted(self.get_max_list_float(inf23,"L","n_sts"),key=lambda k :k["L"])
        #if len(infornata1)<len(infornata2):

        for p in infornata13:
            for p2 in infornata23:
                    if p["L"]==p2["L"]:
                        gap3.append(math.fabs(float(p["E"][0][0])-float(p2["E"][0][0]))) # TODO CHECK
                        lenght3.append(p["L"])

        lenght3=[1.0/y for y in lenght3]
        fig=plt.figure()
        if fit:
            popt3, pcov3 = curve_fit(fit, lenght3,gap3)
            xfine = numpy.linspace( 0, 0.06,100 )  # define values to plot the function for
            #y=(list_entropy[1]-list_entropy[2])/6 + list_entropy[1]

        plt.plot(lenght,gap,'rh--',markersize=12,label='V='+str(v1))
        plt.plot(xfine, fit(xfine, popt[0], popt[1]), 'r-', linewidth=3)
        plt.plot(lenght2,gap2,'b^--',markersize=12,label='V='+str(v2))
        plt.plot(xfine, fit(xfine, popt2[0], popt2[1]), 'b-', linewidth=3)
        plt.plot(lenght3,gap3,'kp--',markersize=12,label='V='+str(v3))
        plt.plot(xfine, fit(xfine, popt3[0], popt3[1]), 'k-', linewidth=3)
        plt.legend(fontsize=15)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.xlabel('1/L', fontsize=30)
        plt.ylabel('$\Delta_{SpinGap}$', fontsize=30)
        if see==0:
            plt.show()
        else:
            plt.savefig('/home/thomas/CLL2Species/scanU20/pictures/SpinGap_compaV1'+str(v1)+'_V3'+str(v3)+'_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'.pdf',bbox_inches='tight')



    def spectral_gap(self,sweep,Lmax,rc,vvv,uuu,rho1,rho2,ns,see,fit=False): # ENERGY DIF SPECTRUM VS SITES
        ''' GS gap
        this function return the plot of values of energy gap in terms of L
        '''
        inf1=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":vvv,"u":uuu,"rc":rc}) ## TODO CHECK
        gap=list()
        lista=list()
        lenght=list()
        for i in range(ns):
            lista.append(list())
            lenght.append(list())
        s0=list()
        infornata1=sorted(self.get_max_list_float(inf1,"L","n_sts"),key=lambda k :k["L"])
        u=0
        for l in infornata1:
            i=0
            for k in range(1,ns):
                #print l["n_sts"]
                #if l["L"]==Lmax:
                lista[i].append(float(l["E"][k][0])-float(l["E"][0][0])) # TODO CHECK
                lenght[i].append(1.0/l["L"])
                i+=1
                #else:
                #lista[i].append(float(l["E"][k])-float(l["E"][0])) # TODO CHECK
                #lenght[i].append(1.0/l["L"])
                #i+=1
        '''
        for i in range(len(lista)):
            lista[i]=lista[i][3:]
            lenght[i]=lenght[i][3:]
            #if p: # if I want to see this I have to put at the last function variable the value (True)
                #print str(l["E1"])+"-"+str(l["E0"])+"="+str(l["E1"] - l["E0"])
            i+=1
        '''
        fig=plt.figure()
        ax=fig.add_subplot(111)
### labels data
        list_e=range(ns)
        c=['r','k','b','c','g','m','y']
        marker=['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X','.']
        listes=['--','-',':','--']
        coppia=list()
        for l in listes:
            for colore in c:
                for m in marker:
                    coppia.append([colore,m,l])
###
        v=1
        for i in range(0,ns-1):
            plt.plot(lenght[i],lista[i],c=coppia[i*16][0],marker=coppia[i+1][1],markersize=14,linestyle=coppia[i*33][2],label="E0-E" + str(v))
            v+=1
        #i=0
        #for xy in zip(lenght,gap):
         #   xy=(xy[0]+0.05,xy[1]+0.05)
          #  ax.annotate('(%s)' % s0[i] , xy=xy, textcoords='data')
           # i+=1

        if fit:
           popt, pcov = curve_fit(fit, lenght[0], lista[0])
           xfine = numpy.linspace(0, lenght[0][0],100 )  # define values to plot the function for
           plt.plot(xfine, fit(xfine, popt[0], popt[1]), 'g-',linewidth=3)

           ax=fig.add_subplot(1,1,1)
           box2 = TextArea("a = "+str(popt[0])+" +/- "+str(math.sqrt(pcov[0,0]))+"\n"+"b = "+str(popt[1])+" +/- "+str(math.sqrt(pcov[1,1])),textprops=dict(color="k",size=20))

           box2 = HPacker(children=[box2],align="center", pad=5, sep=5)

           anchored_box2 = AnchoredOffsetbox(loc=3,child=box2, pad=0.,frameon=True,bbox_to_anchor=(0., 1.02),bbox_transform=ax.transAxes,borderpad=0.)

           ax.add_artist(anchored_box2)

           fig.subplots_adjust(top=0.8)
                             #fig.suptitle('v='+str(v)+',u='+str(uuu)+r'$,r_{c}=$'+str(rc)+r"$,\rho_{1}=$"+str(rho1)+r"$,\rho_{2}=$"+str(rho2),fontsize=20)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.xlabel(r"$1/L$", fontsize=25)
        plt.ylabel('gap', fontsize=25)
        #plt.legend().draggable()
        ax.legend().draggable()
        if see==0:
            plt.show()
        else:
            plt.savefig('/home/thomas/CLL2Species/scanU20/pictures/SpectralGap_V'+str(vvv)+'_rc'+str(rc)+'_U'+str(uuu)+'_na'+str(rho1)+'.pdf',bbox_inches='tight')

    def spec_gap_cluster(self,sweep,rc,v,u,rho1,rho2,ns): # ENERGY DIF SPECTRUM VS SITES
        ''' GS gap
        this function return the plot of values of energy gap in terms of L
        '''
        inf1=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v,"u":u,"rc":rc}) ## TODO CHECK
        gap=list()
        lista=list()
        lenght=list()
        for i in range(ns):
            lista.append(list())
            lenght.append(list())
        s0=list()
        infornata1=sorted(self.get_max_list(inf1,"L","n_sts"),key=lambda k :k["L"])
        i=0
        for l in infornata1:
            i=0
            for k in range(ns):
                if k>0:
                    if self.cluster(l["L"],l["rho"],l["rc"]):
                        lista[i].append(l["E"][k]-l["E"][0])
                        lenght[i].append(1.0/l["L"])
                        i+=1
        #s0.append(l["s0"])
        """
        for i in range(len(lista)):# here I cut the value under 3 place
            lista[i]=lista[i][3:]
            lenght[i]=lenght[i][3:]
        """
        print lista
        print lenght
        i+=1
        fig=plt.figure()
        ax=fig.add_subplot(111)
### labels data
        list_e=range(ns)
        c=['r','k','b','c','g','m','y']
        marker=['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X','.']
        listes=['--','-',':','--']
        coppia=list()
        for l in listes:
            for colore in c:
                for m in marker:
                    coppia.append([colore,m,l])
###
        v=1
        for i in range(ns):
            ax.plot(lenght[i],lista[i],c=coppia[i*16][0],marker=coppia[i+1][1],linestyle=coppia[i*33][2],label="E0-E" + str(v))
            v+=1
        #i=0
        #for xy in zip(lenght,gap):
         #   xy=(xy[0]+0.05,xy[1]+0.05)
          #  ax.annotate('(%s)' % s0[i] , xy=xy, textcoords='data')
           # i+=1
        fig.suptitle('v_'+str(v)+'rc'+str(rc)+'_rho_'+str(rho),fontsize=20)
        plt.xlabel('1/L', fontsize=18)
        plt.ylabel('gap', fontsize=16)
        plt.legend(loc='upper right')
        plt.show()

    ##### CORRELATION ####
    def LocObser(self,rc,v,u,rho1,rho2,L,sweep,reduce,see,observable="aa"):
        '''
        Local observable
        '''
        inf=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v,"u":u,"rc":rc,"L":L})
        infornata=self.get_max_list_float(inf,"L","n_sts")[0]
        print infornata["n_sts"]
        list_sts=list()
        infornate=self.db.find({"sweep":sweep,"L":L,"v":v,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
        for l in infornate:
                list_sts.append(float(l["n_sts"]))
        #print list_sts
        #print max(list_sts)
        Mn_sts=max(list_sts)
        test=list()
        inf2=self.db.find({"sweep":sweep,"L":L,"v":v,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2,"n_sts":list_sts[0]})
        for i in inf2:
            test.append(i["E"][0])
        #print test
        aa=infornata[observable]
        for j in range(0,len(aa)):
            aa[j]=aa[j]*1.0
        #print aa
        #print len(aa)
        lenght=range(L-reduce)
        #print lenght
        fig=plt.figure()
        plt.xticks(fontsize = 25)
        ax=fig.add_subplot(1,1,1)
        plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.5f'))
        plt.yticks(fontsize = 25)
        plt.xlabel("$j$", fontsize=25)
        plt.ylabel(str(observable), fontsize=25)
        plt.plot(lenght,aa,'r8--',linewidth=3, markersize=10)
        if see==0:
            plt.show()
        else:
            plt.savefig(os.path.join(path,'Obs_'+str(observable)+'_V'+str(v)+'_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'.pdf'),bbox_inches='tight')
    def LocObser_compa(self,rc,list_v,u,rho1,rho2,L,sweep,logoption,see,observable="aa"):
        '''
        Local observable
        '''
        '''
        list_v=list()
        infornate=self.db.find({"sweep":sweep,"L":L,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
        for l in infornate:
            list_v.append(l["v"])
        list_v=sorted(list_v) #sorte for the value of interaction

        c=list()
        v=list()
        '''
        d={}
        for x in range(0,len(list_v)):
            d["list_entropy_{0}".format(x)]=list()
            d["list_Lsize_{0}".format(x)]=list()
            d["plot_{0}".format(x)]=plt
        for i in range(0,len(list_v)):
                inf=self.db.find({"sweep":sweep,"u":u,"rc":rc,"L":L,"rho1":rho1,"rho2":rho2,"v":list_v[i]})
                infornata=self.get_max_list_float(inf, "L","n_sts")[0]
                for p in range(0,len(infornata[observable])/2+1):
                    d["list_entropy_{0}".format(i)].append(abs(infornata[observable][p]))#[0:len(infornata[observable])/2+1]
                d["list_Lsize_{0}".format(i)]=range(0,len(infornata[observable])/2+1)

        fig=plt.figure()
        plt.xticks(fontsize = 25)
        plt.xlim(10**0,10)
        plt.ylim(10**-2,10**0)
        ax=fig.add_subplot(1,1,1)
        plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
        plt.yticks(fontsize = 25)
        plt.xlabel("$j$", fontsize=25)
        #plt.ylabel("$<c^{\dagger}_{j,a} c_{k,a}> $", fontsize=25)
        plt.ylabel( '$< O_{cdw}(j) O_{cdw} (k)>$', fontsize=25)
        mark=['s','8','D','o','x','*','P','^']
        linetype=['-','--',':','-.','-']
        if logoption:
            for i in range(0,len(list_v)):
                    plt.loglog(d["list_Lsize_{0}".format(i)],d["list_entropy_{0}".format(i)],linestyle=linetype[i],marker=mark[i],markersize=10,linewidth=3,label='V'+str(list_v[i]))
        else:
            for i in range(0,len(list_v)):
                    plt.plot(d["list_Lsize_{0}".format(i)],d["list_entropy_{0}".format(i)],linestyle=linetype[i],marker=mark[i],markersize=10,linewidth=3,label='V'+str(list_v[i]))
        #plt.plot(lenght,aa,'r8--',linewidth=3, markersize=10)

        plt.legend(fontsize=25)
        if see==0:
            plt.show()
        else:
            plt.savefig(os.path.join(path,'Obscompa_'+str(observable)+'_L'+str(L)+'_V'+str(list_v[0])+'_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'.png'),bbox_inches='tight')

    def MeanObser(self,rc,v,u,rho1,rho2,L,sweep,reduce,see,observable="ss_ab"):

        '''
        Average of Observable per site versus V/t
        '''

        AveOb=list()
        list_v=list()
        infornate=self.db.find({"sweep":sweep,"L":L,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
        #infornate=self.get_max_list_float(inf,"L","n_sts")
        for l in infornate:
            list_v.append(l["v"])
        list_v=sorted(list_v)
        #print list_v
        d={}
        c=list()
        v=list()
        for i in range(0,len(list_v)):
                l=[]
                inf=self.db.find({"sweep":sweep,"u":u,"rc":rc,"L":L,"rho1":rho1,"rho2":rho2,"v":list_v[i]})
                infornata=self.get_max_list_float(inf, "L","n_sts")[0]
                l=infornata[observable]
                #print sum(l)/(L-reduce)
                c.append(sum(l)/(L-reduce))
        fig=plt.figure()
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.5f'))
        plt.xlabel("$V$", fontsize=25)
        plt.ylabel('Mean_'+str(observable), fontsize=25)
        plt.plot(list_v,c,'r8--',linewidth=3, markersize=10)
        if see==0:
            plt.show()
        else:
            plt.savefig('/home/thomas/CLL2Species/scanU20/pictures/MeanObs_'+str(observable)+'_V'+str(v)+'_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'.pdf',bbox_inches='tight')

    def nn_concorr(self,rc,v,u,rho1,rho2,L,sweep):
        '''
        NOT TEST YET
        pair pair correlation
        '''
        inf=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v,"u":u,"rc":rc,"L":L})
        infornata=self.get_max_list_float(inf,"L","n_sts")[0]
        print infornata["n_sts"]
        length=range(L)
        g2=[]
        for i in len(length):
                g2.append(infornata['naa'][i]+infornata['nab'][i]+infornata['nbb'][i]+infornata['nba'][i]-(infornata['naa'][0]+infornata['nbb'][0])**2)

        fig=plt.figure()
        plt.plot(length, g2, 'ro--', linewidth=3, markersize=12)
        plt.legend()
        fig.suptitle('v='+str(v)+r'$,r_{c}=$'+str(rc)+r'$,\rho_{1}=$'+str(rho1)+r'$,\rho_{2}$='+str(rho2)+r'$,L$='+str(L),fontsize=20)
        plt.xlabel('j', fontsize=25)
        plt.ylabel('$g_2 (j)$', fontsize=25)
        plt.show()

    def nn_concorr_log(self,rc,v,u,rho1,rho2,L,sweep):
        '''
        NOT TEST YET
        pair pair correlation comparison logscale
        '''
        inf=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v,"u":u,"rc":rc,"L":L})
        infornata=self.get_max_list(inf,"L","n_sts")[0]
        aa=infornata[observable]
        bb=infornata[observable2]
        lenght=range(L)
        g2=[]
        for i in len(length):
                g2.append(infornata['naa'][i]+infornata['nab'][i]+infornata['nbb'][i]+infornata['nba'][i]-(infornata['naa'][0]+infornata['nbb'][0])**2)

        fig=plt.figure()
        plt.yscale('log')
        plt.xscale('log')
        plt.plot(length, g2, 'ro--', linewidth=3, markersize=12)
        plt.legend()
        plt.legend()
        fig.suptitle('v='+str(v)+r'$,r_{c}=$'+str(rc)+r'$,\rho_{1}=$'+str(rho1)+r'$,\rho_{2}$='+str(rho2)+r'$,L$='+str(L),fontsize=20)
        plt.xlabel('$\log[j]$', fontsize=25)
        plt.ylabel('$\log[g_2 (j)]$', fontsize=25)
        plt.show()

    def nn_precision(self,rc,v1,v2,v3,v4,u,rho1,rho2,L,sweep,n,see,observable='cdw_a'):

        '''
        Function to compute local deviation of an observable
        '''

        list_v=list()
        nststemp=0
        infornate=self.db.find({"sweep":sweep,"L":L,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
        #infornate=self.get_max_list_float(inf,"L","n_sts")
        for l in infornate:
            list_v.append(l["v"])
            nststemp +=float(l["n_sts"])
        nststemp=nststemp/len(list_v)
        list_v=sorted(list_v)
        #print list_v
        d={}
        nave=list()
        v=list()
        for x in range(0,len(list_v)):
            d["list_preci_{0}".format(x)]=list()
            d["list_Lsize_{0}".format(x)]=list()
            d["plot_{0}".format(x)]=plt
        for i in range(0,len(list_v)):
                ntemp=0
                inf=self.db.find({"sweep":sweep,"u":u,"rc":rc,"L":L,"rho1":rho1,"rho2":rho2,"v":list_v[i]})
                infornata=self.get_max_list_float(inf, "L","n_sts")[0]
                for l in infornata[observable]:
                        d["list_preci_{0}".format(i)].append(l-n)
                        ntemp +=abs(l-n)
                        #print len(d["list_preci_{0}".format(i)])
                        d["list_Lsize_{0}".format(i)].append(len(d["list_preci_{0}".format(i)]))
                nave.append(ntemp/L)

        fig=plt.figure()

        for i in range(0,len(list_v)):
           d["plot_{0}".format(i)]=plt.plot(d["list_Lsize_{0}".format(i)],d["list_preci_{0}".format(i)],'s--',label='V'+str(list_v[i]))

        #plt.plot(list_v,nave,'ko--',linewidth=3,markersize=12,label='#states'+str(nststemp))
        plt.legend(loc=3,fontsize=15)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.xlabel('L', fontsize=30)
        plt.ylabel('$n_{j}-n$', fontsize=30)
        if see==0:
            plt.show()
        else:
            plt.savefig('/home/thomas/CLL2Species/scanU20/pictures/nn_precision.pdf',bbox_inches='tight')

    def nn_precision_average(self,rc,u,rho1,rho2,L,sweep,n,see,observable='cdw_a'):

        '''
        Function to compute truncation error or local deviation (of the density) as function of the strength interaction
        Function will be separated in two part when I have time
        '''

        list_v=list()
        nststemp=0
        infornate=self.db.find({"sweep":sweep,"L":L,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
        #infornate=self.get_max_list_float(inf,"L","n_sts")
        for l in infornate:
            list_v.append(l["v"])
            nststemp +=float(l["n_sts"])
        nststemp=nststemp/len(list_v)
        list_v=sorted(list_v)
        d={}
        nave=list()
        trunc=list()
        v=list()
        for x in range(0,len(list_v)):
            d["list_preci_{0}".format(x)]=list()
            d["list_Lsize_{0}".format(x)]=list()
            d["plot_{0}".format(x)]=plt
        for i in range(0,len(list_v)):
                ntemp=0
                inf=self.db.find({"sweep":sweep,"u":u,"rc":rc,"L":L,"rho1":rho1,"rho2":rho2,"v":list_v[i]})
                infornata=self.get_max_list_float(inf, "L","n_sts")[0]
                try:
                    for l in infornata[observable]:
                        d["list_preci_{0}".format(i)].append(l-n)
                        ntemp +=abs(l-n)
                        d["list_Lsize_{0}".format(i)].append(len(d["list_preci_{0}".format(i)]))
                except:
                        trunc.append(infornata[observable])
                        d["list_Lsize_{0}".format(i)].append(len(d["list_preci_{0}".format(i)]))
                nave.append(ntemp/L)

        fig=plt.figure()

        #for i in range(0,len(list_v)):
        #       d["plot_{0}".format(i)]=plt.plot(list_v[i],d["list_preci_{0}".format(i)],'s--',label='V'+str(list_v[i]))
        plt.semilogy(list_v,trunc,'ko--',linewidth=3,markersize=12,label='#states'+str(nststemp))
        ax=fig.add_subplot(1,1,1)
        ax=fig.add_subplot(111)
        ax.set_yticks([10**-5, 10**-6,10**-7])
        plt.legend(loc=2,fontsize=20)
        plt.xlabel('V', fontsize=30)
        #ax.set_yscale('log')
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 35)
        plt.ylabel('truncation error', fontsize=30)
        #plt.ylabel('$1/L*\sum_{j}(n_{j}-n)$', fontsize=30)
        if see==0:
            plt.show()
        else:
            plt.savefig(os.path.join(path,'precision_average_'+str(observable)+'_L_'+str(L)+'_sweep_'+str(sweep)+'_u_'+str(u)+'_rc_'+str(u)+'.png'),bbox_inches='tight')

    def MaxSq(self,rc,v,u,rho1,rho2,L,sweep,reduce,see,observable='naa',fit=None):
        ''' pair pair correlation
        I take all the "infornata" with 2 sweep because my code compute for this value
        '''

        list_v=list()
        nststemp=0
        infornate=self.db.find({"sweep":sweep,"L":L,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
        #infornate=self.get_max_list_float(inf,"L","n_sts")
        for l in infornate:
            list_v.append(l["v"])
            nststemp +=float(l["n_sts"])
        nststemp=nststemp/len(list_v)
        list_v=sorted(list_v)
        d={}
        nave=list()
        trunc=list()
        v=list()
        lenght=range(L-reduce)
        kk=list()
        ka2=list()
        ka=self.FBZ(L)
        numpy.fft.rfft(ka)
        mannaia=0
        max_xy=[0,0]
        Sq=list()
        for l in range(0,len(list_v)):
                inf=self.db.find({"sweep":sweep,"u":u,"rc":rc,"L":L,"rho1":rho1,"rho2":rho2,"v":list_v[l]})
                infornata=self.get_max_list_float(inf, "L","n_sts")[0]
                aa=infornata[observable]
                for k in ka: # fourier transform
                    sommatoria_a = 0 # I have to put the division by L only at the end to improve the calculation
                    for p in lenght:
                        i=1j
                        sommatoria_a += ((aa[p])*cmath.exp(i*p*k)/(L-reduce)) #TODO this is true only if the translational invariance is respected.
                        kk.append(sommatoria_a.real)
                        mannaia+=1
                        if max_xy[1] < sommatoria_a.real:
                            max_xy = [k, sommatoria_a.real]
                        ka2.append(k)
                Sq.append(max_xy[1])
                v.append(list_v[l])
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(v,Sq, 'rh--',linewidth=3 ,markersize=12)
        '''
        ax.annotate("{}".format(max_xy[0]), xy=max_xy)
        ax.set_xticks([0.,.25*cmath.pi,.5*cmath.pi, .75*cmath.pi , cmath.pi])
        ax.set_xticklabels(["$0$", r"$\frac{1}{4}\pi$",r"$\frac{1}{2}\pi$", r"$\frac{3}{4}\pi$", r"$\pi$"])
        ax.tick_params(axis='both', which='major', labelsize=15)
        fig.suptitle('v='+str(v)+',u='+str(u)+r'$,r_{c}=$'+str(rc)+r"$,\rho_{1}=$"+str(rho1)+r"$,\rho_{2}=$"+str(rho2)+','+observable,fontsize=20)
        '''
        plt.xlabel('V/t', fontsize=25)
        plt.ylabel('Max(S(k))_'+str(observable), fontsize=25)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        if see==0:
            plt.show()
        else:
            plt.savefig(os.path.join(path,'MaxSq_'+str(observable)+'.png'),bbox_inches='tight')

    def MaxSq_cdw_sdw(self,rc,v,u,rho1,rho2,L,sweep,reduce,see,observable='naa',fit=None):
        ''' pair pair correlation
        I take all the "infornata" with 2 sweep because my code compute for this value
        '''

        list_v=list()
        nststemp=0
        infornate=self.db.find({"sweep":sweep,"L":L,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
        #infornate=self.get_max_list_float(inf,"L","n_sts")
        for l in infornate:
            list_v.append(l["v"])
            nststemp +=float(l["n_sts"])
        nststemp=nststemp/len(list_v)
        list_v=sorted(list_v)
        d={}
        nave=list()
        trunc=list()
        v=list()
        lenght=range(L-reduce)

        ka=self.FBZ(L)
        numpy.fft.rfft(ka)
        mannaia=0
        max_xy=[0,0]
        Sq=list()
        for l in range(0,len(list_v)):
                kk=list()
                ka2=list()
                inf=self.db.find({"sweep":sweep,"u":u,"rc":rc,"L":L,"rho1":rho1,"rho2":rho2,"v":list_v[l]})
                infornata=self.get_max_list_float(inf, "L","n_sts")[0]
                aa=infornata[observable]
                naa=infornata['naa'][0]
                nbb=infornata['nbb'][0]
                nn=(naa+nbb)**2
                for k in ka: # fourier transform
                    sommatoria_a = 0 # I have to put the division by L only at the end to improve the calculation
                    for p in lenght:
                        i=1j
                        sommatoria_a += ((aa[p]-nn)*cmath.exp(i*p*k)/L) #TODO this is true only if is respected the traslational invariance.
                        mannaia+=1
                        if max_xy[1] < sommatoria_a.real:
                            max_xy = [k, sommatoria_a.real]
                    kk.append(sommatoria_a.real)
                    ka2.append(k)

                #print max(kk)
                Sq.append(max_xy[1])
                v.append(list_v[l])
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(v,Sq, 'rh--',linewidth=3, markersize=10)
        #ax.plot(ka2,kk,'bo--',linewidth=3, markersize=10)
        '''
        ax.annotate("{}".format(max_xy[0]), xy=max_xy)
        ax.set_xticks([0.,.25*cmath.pi,.5*cmath.pi, .75*cmath.pi , cmath.pi])
        ax.set_xticklabels(["$0$", r"$\frac{1}{4}\pi$",r"$\frac{1}{2}\pi$", r"$\frac{3}{4}\pi$", r"$\pi$"])
        ax.tick_params(axis='both', which='major', labelsize=15)
        fig.suptitle('v='+str(v)+',u='+str(u)+r'$,r_{c}=$'+str(rc)+r"$,\rho_{1}=$"+str(rho1)+r"$,\rho_{2}=$"+str(rho2)+','+observable,fontsize=20)
        '''
        plt.xlabel('V/t', fontsize=25)
        plt.ylabel('$Max(S_{CDW}(k))$', fontsize=25)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        if see==0:
            plt.show()
        else:
            plt.savefig(os.path.join(path,'MaxSq_cdw_'+str(observable)+'.png'),bbox_inches='tight')

    def MaxSq_TSs(self,rc,v,u,rho1,rho2,L,sweep,see,fit=None):
        ''' pair pair correlation
        I take all the "infornata" with 2 sweep because my code compute for this value
        '''

        list_v=list()
        nststemp=0
        infornate=self.db.find({"sweep":sweep,"L":L,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
        #infornate=self.get_max_list_float(inf,"L","n_sts")
        for l in infornate:
            list_v.append(l["v"])
            nststemp +=float(l["n_sts"])
        nststemp=nststemp/len(list_v)
        list_v=sorted(list_v)
        d={}
        nave=list()
        trunc=list()
        v=list()
        lenght=range(L-1)
        kk=list()
        ka2=list()
        ka=self.FBZ(L)
        numpy.fft.rfft(ka)
        mannaia=0
        max_xy=[0,0]
        Sq=list()
        for l in range(0,len(list_v)):
                inf=self.db.find({"sweep":sweep,"u":u,"rc":rc,"L":L,"rho1":rho1,"rho2":rho2,"v":list_v[l]})
                infornata=self.get_max_list_float(inf, "L","n_sts")[0]
                Ta=infornata['Tss_0a']
                Taa=infornata['Tss_0aa']
                Tb=infornata['Tss_0b']
                Tbb=infornata['Tss_0bb']
                for k in ka: # fourier transform
                    sommatoria_a = 0 # I have to put the division by L only at the end to improve the calculation
                    for p in lenght:
                        i=1j
                        sommatoria_a += ((Ta[p] +Tb[p] +Taa[p] +Tbb[p] )*cmath.exp(i*p*k)/(L)) #TODO this is true only if is respected the traslational invariance.
                        if max_xy[1] < sommatoria_a.real:
                            max_xy = [k, sommatoria_a.real]
                    ka2.append(k)
                    kk.append(sommatoria_a.real)
                #print max(kk)
                Sq.append(max_xy[1])
                v.append(list_v[l])
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(v,Sq, 'rh--',linewidth=3, markersize=12)
        '''
        ax.annotate("{}".format(max_xy[0]), xy=max_xy)
        ax.set_xticks([0.,.25*cmath.pi,.5*cmath.pi, .75*cmath.pi , cmath.pi])
        ax.set_xticklabels(["$0$", r"$\frac{1}{4}\pi$",r"$\frac{1}{2}\pi$", r"$\frac{3}{4}\pi$", r"$\pi$"])
        ax.tick_params(axis='both', which='major', labelsize=15)
        fig.suptitle('v='+str(v)+',u='+str(u)+r'$,r_{c}=$'+str(rc)+r"$,\rho_{1}=$"+str(rho1)+r"$,\rho_{2}=$"+str(rho2)+','+observable,fontsize=20)
        '''
        plt.xlabel('V/t', fontsize=25)
        plt.ylabel('$S_{TS0}(k=0)$', fontsize=25)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        if see==0:
            plt.show()
        else:
            plt.savefig(os.path.join(path,'MaxSq_TSs.pdf'),bbox_inches='tight')

    def MaxSq_compa(self,rc,v,u,rho1,rho2,L,sweep,reduce,see,fit=None):
        ''' pair pair correlation
        I take all the "infornata" with 2 sweep because my code compute for this value
        '''

        list_v=list()
        nststemp=0
        infornate=self.db.find({"sweep":sweep,"L":L,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
        #infornate=self.get_max_list_float(inf,"L","n_sts")
        for l in infornate:
            list_v.append(l["v"])
            nststemp +=float(l["n_sts"])
        nststemp=nststemp/len(list_v)
        list_v=sorted(list_v)
        d={}
        nave=list()
        trunc=list()
        v=list()
        lenght=range(L)
        lenght2=range(L-1)
        kk1=list()
        ka1=list()
        kk2=list()
        ka2=list()
        kk3=list()
        ka3=list()
        ka=self.FBZ(L)
        numpy.fft.rfft(ka)
        mannaia=0
        cc=list()
        max_xy1=[0,0]
        max_xy2=[0,0]
        max_xy3=[0,0]
        Sq1=list()
        Sq2=list()
        Sq3=list()
        for l in range(0,len(list_v)):
                inf=self.db.find({"sweep":sweep,"u":u,"rc":rc,"L":L,"rho1":rho1,"rho2":rho2,"v":list_v[l]})
                infornata=self.get_max_list_float(inf, "L","n_sts")[0]
                aa=infornata["ss_ab"]
                bb=infornata["Tss_1a"]
                cc1=infornata["Tss_0a"]
                cc2=infornata["Tss_0aa"]
                cc3=infornata["Tss_0b"]
                cc4=infornata["Tss_0bb"]
                for y in range(0,len(cc1)):
                    cc.append(cc1[y]+cc2[y]+cc3[y]+cc4[y])
                for k in ka: # fourier transform
                    sommatoria_a = 0 # I have to put the division by L only at the end to improve the calculation
                    sommatoria_b = 0
                    sommatoria_c = 0
                    for p in lenght:
                        i=1j
                        sommatoria_a += ((aa[p])*cmath.exp(i*p*k)/(L)) #TODO this is true only if the translational invariance is respected.
                        mannaia+=1
                        if max_xy1[1] < sommatoria_a.real:
                            max_xy1 = [k, sommatoria_a.real]
                    kk1.append(sommatoria_a.real)
                    ka1.append(k)
                    for p in lenght2:
                        i=1j
                        sommatoria_b += ((bb[p])*cmath.exp(i*p*k)/(L))
                        sommatoria_c += ((cc[p])*cmath.exp(i*p*k)/(L)) #TODO this is true only if the translational invariance is respected.
                        if max_xy2[1] < sommatoria_b.real:
                            max_xy2 = [k, sommatoria_b.real]
                        if max_xy3[1] < sommatoria_c.real:
                            max_xy3= [k, sommatoria_c.real]
                    kk3.append(sommatoria_c.real)
                    ka3.append(k)
                    kk2.append(sommatoria_b.real)
                    ka2.append(k)
                Sq1.append(max_xy1[1])
                Sq2.append(max_xy2[1])
                Sq3.append(max_xy3[1])
                v.append(list_v[l])
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(v,Sq1, 'rh--',linewidth=3 ,markersize=12,label='$SS$')
        ax.plot(v,Sq2, 'gD--',linewidth=3 ,markersize=12,label='$TS_1$')
        ax.plot(v,Sq3, 'bo--',linewidth=3 ,markersize=12,label='$TS_0$')
        '''
        ax.annotate("{}".format(max_xy[0]), xy=max_xy)
        ax.set_xticks([0.,.25*cmath.pi,.5*cmath.pi, .75*cmath.pi , cmath.pi])
        ax.set_xticklabels(["$0$", r"$\frac{1}{4}\pi$",r"$\frac{1}{2}\pi$", r"$\frac{3}{4}\pi$", r"$\pi$"])
        ax.tick_params(axis='both', which='major', labelsize=15)
        fig.suptitle('v='+str(v)+',u='+str(u)+r'$,r_{c}=$'+str(rc)+r"$,\rho_{1}=$"+str(rho1)+r"$,\rho_{2}=$"+str(rho2)+','+observable,fontsize=20)
        '''
        plt.xlabel('V/t', fontsize=25)
        plt.ylabel('Max(S(k))', fontsize=25)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.legend(fontsize=25)
        if see==0:
            plt.show()
        else:
            plt.savefig(os.path.join(path,'MaxSq_compa_U'+str(u)+'_L_'+str(L)+'.png'),bbox_inches='tight')


    def MaxSq_compa2(self,list_L,rc,v,list_u,rho1,rho2,n,see,observable='cdw',fit=None):
        ''' pair pair correlation
        I take all the "infornata" with 2 sweep because my code compute for this value
        '''
        d={}
        count=0
        for x in range(0,len(list_L)*len(list_u)):
                        d["list_Sq_{0}".format(x)]=list()
                        d["list_v_{0}".format(x)]=list()
                        d["list_max_xy_{0}".format(x)]=[0,0]
                        d["plot_{0}".format(x)]=plt
        for y in range(0,len(list_u)):
            for j in range(0,len(list_L)):
                list_v=list()
                infornate=0
                try:
                    infornate=self.db.find({"u":list_u[y],"L":list_L[j],"rc":rc,"rho1":rho1,"rho2":rho2})

                    for l in infornate:
                        list_v.append(l["v"])
                    list_v=sorted(list_v)
                    lenght=range(list_L[j])
                    kk1=list()
                    ka1=list()
                    ka=self.FBZ(list_L[j])
                    numpy.fft.rfft(ka)
                    mannaia=0
                    cc=list()
                    Sq1=list()
                    for l in range(0,len(list_v)):
                            inf=self.db.find({"u":list_u[y],"rc":rc,"L":list_L[j],"rho1":rho1,"rho2":rho2,"v":list_v[l]})
                            infornata=self.get_max_list_float(inf, "L","n_sts")[0]
                            aa=infornata[observable]

                            cc1=infornata["aa"][0]
                            cc2=infornata["bb"][0]
                            nn=(cc1+cc2)**2
                            for k in ka: # fourier transform
                                sommatoria_a = 0 # I have to put the division by L only at the end to improve the calculation
                                for p in lenght:
                                    i=1j
                                    if n:
                                        sommatoria_a += ((aa[p]-nn)*cmath.exp(i*p*k)/(list_L[j])) #TODO this is true only if the translational invariance is respected.
                                        mannaia+=1
                                        if d["list_max_xy_{0}".format(count)][1] < sommatoria_a.real:
                                            d["list_max_xy_{0}".format(count)] = [k, sommatoria_a.real]
                                    else:
                                        sommatoria_a += ((aa[p])*cmath.exp(i*p*k)/(list_L[j])) #TODO this is true only if the translational invariance is respected.
                                        mannaia+=1
                                        if d["list_max_xy_{0}".format(count)][1] < sommatoria_a.real:
                                            d["list_max_xy_{0}".format(count)] = [k, sommatoria_a.real]
                            kk1.append(sommatoria_a.real)
                            ka1.append(k)
                            d["list_Sq_{0}".format(count)].append(d["list_max_xy_{0}".format(count)][1])
                            d["list_v_{0}".format(count)].append(list_v[l])

                except:
                        print "no data"
                count=count+1
        fig=plt.figure()
        ax=fig.add_subplot(111)
        mark=['s','8','D','o','x','*','P','^']
        linetype=['-','--',':','-.','-']
        col=['r','b','g','c','k','y']
        count=0
        #print d["list_Sq_{0}".format(1)]
        for y in range(0,len(list_u)):
            for j in range(0,len(list_L)):
                if not d["list_Sq_{0}".format(count)]:
                    print "no data for L="+str(list_L[j])+' and U='+str(list_u[y])
                else:
                        ax.plot(d["list_v_{0}".format(count)],d["list_Sq_{0}".format(count)], '--',marker=mark[j+y],color=col[j+y],linewidth=3 ,markersize=10,markeredgecolor='k',label='L='+str(list_L[j])+', U='+str(list_u[y]))
                count=count+1
        plt.xlabel('V/t', fontsize=25)
        plt.ylabel('Max(S(k))_'+str(observable), fontsize=25)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.legend(fontsize=25)
        if see==0:
            plt.show()
        else:
            plt.savefig(os.path.join(path,'MaxSq_'+str(observable)+'.png'),bbox_inches='tight')


    def MaxSq_compa_dw(self,rc,v,u,rho1,rho2,L,sweep,see,fit=None):
        ''' pair pair correlation
        I take all the "infornata" with 2 sweep because my code compute for this value
        '''
        list_v=list()
        nststemp=0
        infornate=self.db.find({"sweep":sweep,"L":L,"u":u,"rc":rc,"rho1":rho1,"rho2":rho2})
        #infornate=self.get_max_list_float(inf,"L","n_sts")
        for l in infornate:
            list_v.append(l["v"])
            nststemp +=float(l["n_sts"])
        nststemp=nststemp/len(list_v)
        list_v=sorted(list_v)

        nave=list()
        trunc=list()
        v=list()
        lenght=range(L)
        lenght2=range(L)
        kk1=list()
        ka1=list()
        kk2=list()
        ka2=list()
        kk3=list()
        ka3=list()
        ka=self.FBZ(L)
        numpy.fft.rfft(ka)
        mannaia=0
        cc=list()
        max_xy1=[0,0]
        max_xy2=[0,0]
        max_xy3=[0,0]
        Sq1=list()
        Sq2=list()
        Sq3=list()
        for l in range(0,len(list_v)):
                inf=self.db.find({"sweep":sweep,"u":u,"rc":rc,"L":L,"rho1":rho1,"rho2":rho2,"v":list_v[l]})
                infornata=self.get_max_list_float(inf, "L","n_sts")[0]
                aa=infornata["cdw"]
                bb=infornata["sdw"]
                cc1=infornata["aa"][0]
                cc2=infornata["bb"][0]
                nn=(cc1+cc2)**2
                for k in ka: # fourier transform
                    sommatoria_a = 0 # I have to put the division by L only at the end to improve the calculation
                    sommatoria_b = 0
                    sommatoria_c = 0
                    for p in lenght:
                        i=1j
                        sommatoria_a += ((aa[p]-nn)*cmath.exp(i*p*k)/(L)) #TODO this is true only if the translational invariance is respected.
                        mannaia+=1
                        if max_xy1[1] < sommatoria_a.real:
                            max_xy1 = [k, sommatoria_a.real]
                    kk1.append(sommatoria_a.real)
                    ka1.append(k)
                    for p in lenght2:
                        i=1j
                        sommatoria_b += ((bb[p])*cmath.exp(i*p*k)/(L))
                        if max_xy2[1] < sommatoria_b.real:
                            max_xy2 = [k, sommatoria_b.real]
                    kk2.append(sommatoria_b.real)
                    ka2.append(k)
                Sq1.append(max_xy1[1])
                Sq2.append(max_xy2[1])
                v.append(list_v[l])
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(v,Sq1, 'rh--',linewidth=3 ,markersize=12,label='CDW $< O_{cdw}(j) O_{cdw}(k)> $')
        ax.plot(v,Sq2, 'gD--',linewidth=3 ,markersize=12,label='SDW $< O_{sdw}(j) O_{sdw}(k)> $')
        '''
        ax.annotate("{}".format(max_xy[0]), xy=max_xy)
        ax.set_xticks([0.,.25*cmath.pi,.5*cmath.pi, .75*cmath.pi , cmath.pi])
        ax.set_xticklabels(["$0$", r"$\frac{1}{4}\pi$",r"$\frac{1}{2}\pi$", r"$\frac{3}{4}\pi$", r"$\pi$"])
        ax.tick_params(axis='both', which='major', labelsize=15)
        fig.suptitle('v='+str(v)+',u='+str(u)+r'$,r_{c}=$'+str(rc)+r"$,\rho_{1}=$"+str(rho1)+r"$,\rho_{2}=$"+str(rho2)+','+observable,fontsize=20)
        '''
        plt.xlabel('V/t', fontsize=25)
        plt.ylabel('Max(S(k))', fontsize=25)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.legend(fontsize=25)
        if see==0:
            plt.show()
        else:
            plt.savefig(os.path.join(path,'MaxSq_compa_dw_U'+str(u)+'_L_'+str(L)+'.png'),bbox_inches='tight')


    def str_fact(self,rc,v,u,rho1,rho2,L,sweep,reduce,see,observable='naa',fit=None):
        ''' pair pair correlation
        I take all the "infornata" with 2 sweep because my code compute for this value
        '''
        inf=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v,"u":u,"rc":rc,"L":L})
        infornata=self.get_max_list(inf,"L","n_sts")[0]
        aa=infornata[observable]
        nanb=[]
        naa=infornata['naa'][0]
        nbb=infornata['nbb'][0]
        nn=(naa-nbb)**2
        if observable[1] != observable[2]:
            print "eccolo!"
            nanb.append(infornata['naa'][0])
            nanb.append(infornata['nbb'][0])
        else:
            nanb.append(aa[0])
            nanb.append(aa[0])
        lenght=range(L-reduce)
        kk=list()
        ka=self.FBZ(L)
        numpy.fft.rfft(ka)
        mannaia=0
        ka2=list()
        max_xy=[0,0]
        for k in ka: # fourier transform
            sommatoria_a = 0 # I have to put the division by L only at the end to improve the calculation
            for p in lenght:
                i=1j
                sommatoria_a += ((aa[p])*cmath.exp(i*p*k)/L) #TODO this is true only if is respected the traslational invariance.

                mannaia+=1
                if max_xy[1] < sommatoria_a.real:
                    max_xy = [k, sommatoria_a.real]
            ka2.append(k)
            kk.append(sommatoria_a.real)

        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(ka2,kk, 'rh--',linewidth=3, markersize=12)

        if fit: # I try to obtain a good value of LUTTINGER PARAMETER fitting the lowest point of structure factor
           popt, pcov = curve_fit(fit, ka2[:5], kk[:5])
           xfine = numpy.linspace(ka2[0], ka2[9],100 )  # define values to plot the function for
           plt.plot(xfine, fit(xfine, popt[0], popt[1]), 'g-')
           box1 = TextArea("a = "+str(popt[0])+" +/- "+str(math.sqrt(pcov[0,0]))+"\n"+"b = "+str(popt[1])+" +/- "+str(math.sqrt(pcov[1,1]))+"\n"+"K = "+str(popt[0]*math.pi/2),textprops=dict(color="k"))

           box = HPacker(children=[box1],align="center", pad=5, sep=5)

           anchored_box = AnchoredOffsetbox(loc=3,child=box, pad=0.,frameon=True,bbox_to_anchor=(0., 1.02),bbox_transform=ax.transAxes,borderpad=0.)

           ax.add_artist(anchored_box)


        ax.annotate("{}".format(max_xy[0]), xy=max_xy)
        ax.set_xticks([0.,.25*cmath.pi,.5*cmath.pi, .75*cmath.pi , cmath.pi, 1.25*cmath.pi, 1.5*cmath.pi, 1.75*cmath.pi ,2.0*cmath.pi ])
        ax.set_xticklabels(["$0$", r"$\frac{1}{4}\pi$",r"$\frac{1}{2}\pi$", r"$\frac{3}{4}\pi$", r"$\pi$",  r"$\frac{5}{4}\pi$",r"$\frac{3}{2}\pi$", r"$\frac{7}{4}\pi$",r"$2\pi$"])
        ax.tick_params(axis='both', which='major', labelsize=15)
        fig.suptitle('v='+str(v)+',u='+str(u)+r'$,r_{c}=$'+str(rc)+r"$,\rho_{1}=$"+str(rho1)+r"$,\rho_{2}=$"+str(rho2)+','+observable,fontsize=20)
        plt.xlabel('k', fontsize=25)
        plt.ylabel('S(k)_'+str(observable), fontsize=25)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        if see==0:
            plt.show()
        else:
            plt.savefig(os.path.join(path,'Sq_'+str(observable)+'_V'+str(v)+'_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'_L'+str(L)+'.pdf'),bbox_inches='tight')

    def str_facT(self,rc,v,u,rho1,rho2,L,sweep,see,fit=None):
        ''' pair pair correlation
        I take all the "infornata" with 2 sweep because my code compute for this value
        '''
        inf=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v,"u":u,"rc":rc,"L":L})
        infornata=self.get_max_list_float(inf,"L","n_sts")[0]
        print infornata["n_sts"]
        aa=[]
        for i in range(L):
            aa.append(infornata['naa'][i] + infornata['nab'][i]+infornata['nbb'][i] + infornata['nba'][i])
            cost=infornata['naa'][0]**2 + 2*infornata['naa'][0]*infornata['nbb'][0] + infornata['nbb'][0]**2
        lenght=range(L)
        kk=list()
        ka=self.FBZ(L)
        numpy.fft.rfft(ka)
        mannaia=0
        ka2=list()
        max_xy=[0,0]
        for k in ka: # fourier transform
            sommatoria_a = 0 # I have to put the division by L only at the end to improve the calculation
            for p in lenght:
                i=1j
                sommatoria_a += (1.0/L)*((aa[p]-cost)*cmath.exp(i*p*k)) #TODO this is true only if is respected the traslational invariance.
            if mannaia < L/2:
                kk.append(sommatoria_a.real)
                mannaia+=1
                if max_xy[1] < sommatoria_a.real:
                    max_xy = [k, sommatoria_a.real]
                ka2.append(k)
            else:
                break
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(ka2,kk, 'ro--', linewidth = 3, markersize =12)

        if fit: # I try to obtain a good value of LUTTINGER PARAMETER fitting the lowest point of structure factor
           popt, pcov = curve_fit(fit, ka2[:5], kk[:5])
           xfine = numpy.linspace(ka2[0], ka2[9],100 )  # define values to plot the function for
           plt.plot(xfine, fit(xfine, popt[0], popt[1]), 'g-')
           box1 = TextArea("a = "+str(popt[0])+" +/- "+str(math.sqrt(pcov[0,0]))+"\n"+"b = "+str(popt[1])+" +/- "+str(math.sqrt(pcov[1,1]))+"\n"+"K = "+str(popt[0]*math.pi/2),textprops=dict(color="k"))

           box = HPacker(children=[box1],align="center", pad=5, sep=5)

           anchored_box = AnchoredOffsetbox(loc=3,child=box, pad=0.,frameon=True,bbox_to_anchor=(0., 1.02),bbox_transform=ax.transAxes,borderpad=0.)

           ax.add_artist(anchored_box)


        ax.annotate("{}".format(max_xy[0]), weight='bold', xy=max_xy) #xytext=(2.5,0.009)
        ax.set_xticks([0.,.25*cmath.pi,.5*cmath.pi, .75*cmath.pi , cmath.pi])
        ax.set_xticklabels(["$0$", r"$\frac{1}{4}\pi$",
                     r"$\frac{1}{2}\pi$", r"$\frac{3}{4}\pi$", r"$\pi$"])
        ax.tick_params(axis='both', which='major', labelsize=30)
        ax = fig.add_subplot(111)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)
        plt.xlabel('k', fontsize=30)
        plt.ylabel('S(k)', fontsize=30)
        if see==0:
            plt.show()
        else:
            plt.savefig('/home/thomas/CLL2Species/scanU20/pictures/S(k)_V'+str(v)+'_rc'+str(rc)+'_U'+str(u)+'_na'+str(rho1)+'_L'+str(L)+'.pdf',bbox_inches='tight')

    def str_facT_compa(self,choice,rc1,rc2,rc3,rc4,v1,v2,v3,v4,u,rho1,rho2,L,sweep,see,fit=None):

        # Function to plot the comparison between different parameters of the Bose Hubbard model extended
        inf=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v1,"u":u,"rc":rc1,"L":L})
        infornata=self.get_max_list_float(inf,"L","n_sts")[0]
        print infornata["n_sts"]
        aa=[]
        for i in range(L):
            aa.append(infornata['naa'][i] + infornata['nab'][i]+infornata['nbb'][i] + infornata['nba'][i])
            cost=infornata['naa'][0]**2 + 2*infornata['naa'][0]*infornata['nbb'][0] + infornata['nbb'][0]**2
        lenght=range(L)
        kk=list()
        ka=self.FBZ(L)
        numpy.fft.rfft(ka)
        mannaia=0
        ka2=list()
        max_xy=[0,0]
        for k in ka: # fourier transform
            sommatoria_a = 0 # I have to put the division by L only at the end to improve the calculation
            for p in lenght:
                i=1j
                sommatoria_a += ((aa[p]-cost)*cmath.exp(i*p*k)/L) #TODO this is true only if is respected the traslational invariance.
            if mannaia < L/2:
                kk.append(sommatoria_a.real)
                mannaia+=1
                if max_xy[1] < sommatoria_a.real:
                    max_xy = [k, sommatoria_a.real]
                ka2.append(k)
            else:
                break

        inf2=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v2,"u":u,"rc":rc2,"L":L})
        infornata2=self.get_max_list_float(inf2,"L","n_sts")[0]
        print infornata2["n_sts"]
        aa2=[]
        for i in range(L):
            aa2.append(infornata2['naa'][i] + infornata2['nab'][i]+infornata2['nbb'][i] + infornata2['nba'][i])
            cost2=infornata2['naa'][0]**2 + 2*infornata2['naa'][0]*infornata2['nbb'][0] + infornata2['nbb'][0]**2
        lenght=range(L)
        kk2=list()
        ka=self.FBZ(L)
        numpy.fft.rfft(ka)
        mannaia=0
        ka22=list()
        max_xy=[0,0]
        for k in ka2: # fourier transform
            sommatoria_a = 0 # I have to put the division by L only at the end to improve the calculation
            for p in lenght:
                i=1j
                sommatoria_a += ((aa2[p]-cost2)*cmath.exp(i*p*k)/L) #TODO this is true only if is respected the traslational invariance.
            if mannaia < L/2:
                kk2.append(sommatoria_a.real)
                mannaia+=1
                if max_xy[1] < sommatoria_a.real:
                    max_xy = [k, sommatoria_a.real]
                ka22.append(k)
            else:
                break

        if fit: # I try to obtain a good value of LUTTINGER PARAMETER fitting the lowest point of structure factor
           popt, pcov = curve_fit(fit, ka2[:5], kk[:5])
           xfine = numpy.linspace(ka2[0], ka2[9],100 )  # define values to plot the function for
           plt.plot(xfine, fit(xfine, popt[0], popt[1]), 'g-')
           box1 = TextArea("a = "+str(popt[0])+" +/- "+str(math.sqrt(pcov[0,0]))+"\n"+"b = "+str(popt[1])+" +/- "+str(math.sqrt(pcov[1,1]))+"\n"+"K = "+str(popt[0]*math.pi/2),textprops=dict(color="k"))

           box = HPacker(children=[box1],align="center", pad=5, sep=5)

           anchored_box = AnchoredOffsetbox(loc=3,child=box, pad=0.,frameon=True,bbox_to_anchor=(0., 1.02),bbox_transform=ax.transAxes,borderpad=0.)

           ax.add_artist(anchored_box)


        inf3=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v3,"u":u,"rc":rc3,"L":L})
        infornata3=self.get_max_list_float(inf3,"L","n_sts")[0]
        print infornata3["n_sts"]
        aa3=[]
        for i in range(L):
            aa3.append(infornata3['naa'][i] + infornata3['nab'][i]+infornata3['nbb'][i] + infornata3['nba'][i])
            cost3=infornata3['naa'][0]**2 + 2*infornata3['naa'][0]*infornata3['nbb'][0] + infornata3['nbb'][0]**2
        lenght=range(L)
        kk3=list()
        ka=self.FBZ(L)
        numpy.fft.rfft(ka)
        mannaia=0
        ka23=list()
        max_xy=[0,0]
        for k in ka: # fourier transform
            sommatoria_a = 0 # I have to put the division by L only at the end to improve the calculation
            for p in lenght:
                i=1j
                sommatoria_a += ((aa3[p]-cost3)*cmath.exp(i*p*k)/L) #TODO this is true only if is respected the traslational invariance.
            if mannaia < L/2:
                kk3.append(sommatoria_a.real)
                mannaia+=1
                if max_xy[1] < sommatoria_a.real:
                    max_xy = [k, sommatoria_a.real]
                ka23.append(k)
            else:
                break

        inf4=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v4,"u":u,"rc":rc4,"L":L})
        infornata4=self.get_max_list_float(inf4,"L","n_sts")[0]
        print infornata3["n_sts"]
        aa4=[]
        for i in range(L):
            aa4.append(infornata4['naa'][i] + infornata4['nab'][i]+infornata4['nbb'][i] + infornata4['nba'][i])
            cost4=infornata4['naa'][0]**2 + 2*infornata4['naa'][0]*infornata4['nbb'][0] + infornata4['nbb'][0]**2
        lenght=range(L)
        kk4=list()
        ka=self.FBZ(L)
        numpy.fft.rfft(ka)
        mannaia=0
        ka24=list()
        max_xy=[0,0]
        for k in ka: # fourier transform
            sommatoria_a = 0 # I have to put the division by L only at the end to improve the calculation
            for p in lenght:
                i=1j
                sommatoria_a += ((aa4[p]-cost4)*cmath.exp(i*p*k)/L) #TODO this is true only if is respected the traslational invariance.
            if mannaia < L/2:
                kk4.append(sommatoria_a.real)
                mannaia+=1
                if max_xy[1] < sommatoria_a.real:
                    max_xy = [k, sommatoria_a.real]
                ka24.append(k)
            else:
                break
        fig=plt.figure()
        ax2=fig.add_subplot(111)
        ax2.plot(ka2,kk, 'rh--')

        if fit: # I try to obtain a good value of LUTTINGER PARAMETER fitting the lowest point of structure factor
           popt, pcov = curve_fit(fit, ka2[:5], kk[:5])
           xfine = numpy.linspace(ka2[0], ka2[9],100 )  # define values to plot the function for
           plt.plot(xfine, fit(xfine, popt[0], popt[1]), 'g-')
           box1 = TextArea("a = "+str(popt[0])+" +/- "+str(math.sqrt(pcov[0,0]))+"\n"+"b = "+str(popt[1])+" +/- "+str(math.sqrt(pcov[1,1]))+"\n"+"K = "+str(popt[0]*math.pi/2),textprops=dict(color="k"))

           box = HPacker(children=[box1],align="center", pad=5, sep=5)

           anchored_box = AnchoredOffsetbox(loc=3,child=box, pad=0.,frameon=True,bbox_to_anchor=(0., 1.02),bbox_transform=ax.transAxes,borderpad=0.)

           ax2.add_artist(anchored_box)

        ax2.annotate("{}".format(max_xy[0]), xy=max_xy)
        ax2.set_xticks([0.,.25*cmath.pi,.5*cmath.pi, .75*cmath.pi , cmath.pi])
        ax2.set_xticklabels(["$0$", r"$\frac{1}{4}\pi$",r"$\frac{1}{2}\pi$", r"$\frac{3}{4}\pi$", r"$\pi$"])
        ax2.tick_params(axis='both', which='major', labelsize=15)
        plt.xlabel('k', fontsize=30)
        plt.ylabel('S(k)', fontsize=30)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)

	if choice==1: # as function of V
        	plt.plot(ka2,kk,'rh--',linewidth=3,markersize=12,label='V='+str(v1))
        	plt.plot(ka22,kk2,'bs--',linewidth=3,markersize=12,label='V='+str(v2))
        	plt.plot(ka23,kk3,'g8--',linewidth=3,markersize=12, label='V='+str(v3))
        	plt.plot(ka24,kk4,'kv--',linewidth=3,markersize=12, label='V='+str(v4))
        	plt.legend(fontsize=15)
	elif choice==0: # as function of rc
        	plt.plot(ka2,kk,'rh--',linewidth=3,markersize=12,label='rc='+str(rc1))
        	plt.plot(ka22,kk2,'bs--',linewidth=3,markersize=12,label='rc='+str(rc2))
        	plt.plot(ka23,kk3,'g8--',linewidth=3,markersize=12, label='rc='+str(rc3))
        	plt.plot(ka24,kk4,'kv--',linewidth=3,markersize=12, label='rc='+str(rc4))
        	plt.legend(fontsize=15)

        if see==0:
            plt.show()
        else:
            plt.savefig(os.path.join(path,'S(K)_compa_U'+str(u)+'_V'+str(v1)+'_na'+str(rho1)+'_L'+str(L)+'.png'),bbox_inches='tight')
    def str_facT_sdw_compa(self,choice,rc1,rc2,rc3,rc4,v1,v2,v3,v4,u,rho1,rho2,L,sweep,see,fit=None):

        # Function to plot the comparison between different parameters of the Bose Hubbard model extended
        inf=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v1,"u":u,"rc":rc1,"L":L})
        infornata=self.get_max_list_float(inf,"L","n_sts")[0]
        print infornata["n_sts"]
        aa=[]
        aa=infornata["sdw"]
        lenght=range(L)
        kk=list()
        ka=self.FBZ(L)
        numpy.fft.rfft(ka)
        mannaia=0
        ka2=list()
        max_xy=[0,0]
        for k in ka: # fourier transform
            sommatoria_a = 0 # I have to put the division by L only at the end to improve the calculation
            for p in lenght:
                i=1j
                sommatoria_a += ((aa[p])*cmath.exp(i*p*k)/L) #TODO this is true only if is respected the traslational invariance.

                mannaia+=1
                if max_xy[1] < sommatoria_a.real:
                    max_xy = [k, sommatoria_a.real]
            ka2.append(k)
            kk.append(sommatoria_a.real)

        inf2=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v2,"u":u,"rc":rc2,"L":L})
        infornata2=self.get_max_list_float(inf2,"L","n_sts")[0]
        print infornata2["n_sts"]
        aa2=[]
        aa2=infornata2["sdw"]
        lenght=range(L)
        kk2=list()
        ka=self.FBZ(L)
        numpy.fft.rfft(ka)
        mannaia=0
        ka22=list()
        max_xy=[0,0]
        for k in ka2: # fourier transform
            sommatoria_a = 0 # I have to put the division by L only at the end to improve the calculation
            for p in lenght:
                i=1j
                sommatoria_a += ((aa2[p])*cmath.exp(i*p*k)/L) #TODO this is true only if is respected the traslational invariance.

                mannaia+=1
                if max_xy[1] < sommatoria_a.real:
                    max_xy = [k, sommatoria_a.real]
            ka22.append(k)
            kk2.append(sommatoria_a.real)

        if fit: # I try to obtain a good value of LUTTINGER PARAMETER fitting the lowest point of structure factor
           popt, pcov = curve_fit(fit, ka2[:5], kk[:5])
           xfine = numpy.linspace(ka2[0], ka2[9],100 )  # define values to plot the function for
           plt.plot(xfine, fit(xfine, popt[0], popt[1]), 'g-')
           box1 = TextArea("a = "+str(popt[0])+" +/- "+str(math.sqrt(pcov[0,0]))+"\n"+"b = "+str(popt[1])+" +/- "+str(math.sqrt(pcov[1,1]))+"\n"+"K = "+str(popt[0]*math.pi/2),textprops=dict(color="k"))

           box = HPacker(children=[box1],align="center", pad=5, sep=5)

           anchored_box = AnchoredOffsetbox(loc=3,child=box, pad=0.,frameon=True,bbox_to_anchor=(0., 1.02),bbox_transform=ax.transAxes,borderpad=0.)

           ax.add_artist(anchored_box)


        inf3=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v3,"u":u,"rc":rc3,"L":L})
        infornata3=self.get_max_list_float(inf3,"L","n_sts")[0]
        print infornata3["n_sts"]
        aa3=[]
        aa3=infornata3["sdw"]
        lenght=range(L)
        kk3=list()
        ka=self.FBZ(L)
        numpy.fft.rfft(ka)
        mannaia=0
        ka23=list()
        max_xy=[0,0]
        for k in ka: # fourier transform
            sommatoria_a = 0 # I have to put the division by L only at the end to improve the calculation
            for p in lenght:
                i=1j
                sommatoria_a += ((aa3[p])*cmath.exp(i*p*k)/L) #TODO this is true only if is respected the traslational invariance.
                mannaia+=1
                if max_xy[1] < sommatoria_a.real:
                    max_xy = [k, sommatoria_a.real]
            ka23.append(k)
            kk3.append(sommatoria_a.real)
        inf4=self.db.find({"rho1":rho1,'rho2':rho2,"sweep":sweep,"v":v4,"u":u,"rc":rc4,"L":L})
        infornata4=self.get_max_list_float(inf4,"L","n_sts")[0]
        print infornata3["n_sts"]
        aa4=[]
        aa4=infornata4["sdw"]
        lenght=range(L)
        kk4=list()
        ka=self.FBZ(L)
        numpy.fft.rfft(ka)
        mannaia=0
        ka24=list()
        max_xy=[0,0]
        for k in ka: # fourier transform
            sommatoria_a = 0 # I have to put the division by L only at the end to improve the calculation
            for p in lenght:
                i=1j
                sommatoria_a += ((aa4[p])*cmath.exp(i*p*k)/L) #TODO this is true only if is respected the traslational invariance.

                mannaia+=1
                if max_xy[1] < sommatoria_a.real:
                    max_xy = [k, sommatoria_a.real]
            ka24.append(k)
            kk4.append(sommatoria_a.real)
        fig=plt.figure()
        ax2=fig.add_subplot(111)
        ax2.plot(ka2,kk, 'rh--')

        if fit: # I try to obtain a good value of LUTTINGER PARAMETER fitting the lowest point of structure factor
           popt, pcov = curve_fit(fit, ka2[:5], kk[:5])
           xfine = numpy.linspace(ka2[0], ka2[9],100 )  # define values to plot the function for
           plt.plot(xfine, fit(xfine, popt[0], popt[1]), 'g-')
           box1 = TextArea("a = "+str(popt[0])+" +/- "+str(math.sqrt(pcov[0,0]))+"\n"+"b = "+str(popt[1])+" +/- "+str(math.sqrt(pcov[1,1]))+"\n"+"K = "+str(popt[0]*math.pi/2),textprops=dict(color="k"))

           box = HPacker(children=[box1],align="center", pad=5, sep=5)

           anchored_box = AnchoredOffsetbox(loc=3,child=box, pad=0.,frameon=True,bbox_to_anchor=(0., 1.02),bbox_transform=ax.transAxes,borderpad=0.)

           ax2.add_artist(anchored_box)

        ax2.annotate("{}".format(max_xy[0]), xy=max_xy)
        ax2.set_xticks([0.,.25*cmath.pi,.5*cmath.pi, .75*cmath.pi , cmath.pi])
        ax2.set_xticklabels(["$0$", r"$\frac{1}{4}\pi$",r"$\frac{1}{2}\pi$", r"$\frac{3}{4}\pi$", r"$\pi$",  r"$\frac{5}{4}\pi$",r"$\frac{3}{2}\pi$", r"$\frac{7}{4}\pi$",r"$2\pi$"])
        ax2.tick_params(axis='both', which='major', labelsize=15)
        plt.xlabel('k', fontsize=30)
        plt.ylabel('S(k)', fontsize=30)
        plt.xticks(fontsize = 25)
        plt.yticks(fontsize = 25)

	if choice==1: # as function of V
        	plt.plot(ka2,kk,'rh--',linewidth=3,markersize=12,label='V='+str(v1))
        	plt.plot(ka22,kk2,'bs--',linewidth=3,markersize=12,label='V='+str(v2))
        	plt.plot(ka23,kk3,'g8--',linewidth=3,markersize=12, label='V='+str(v3))
        	plt.plot(ka24,kk4,'kv--',linewidth=3,markersize=12, label='V='+str(v4))
        	plt.legend(fontsize=15)
	elif choice==0: # as function of rc
        	plt.plot(ka2,kk,'rh--',linewidth=3,markersize=12,label='rc='+str(rc1))
        	plt.plot(ka22,kk2,'bs--',linewidth=3,markersize=12,label='rc='+str(rc2))
        	plt.plot(ka23,kk3,'g8--',linewidth=3,markersize=12, label='rc='+str(rc3))
        	plt.plot(ka24,kk4,'kv--',linewidth=3,markersize=12, label='rc='+str(rc4))
        	plt.legend(fontsize=15)

        if see==0:
            plt.show()
        else:
            plt.savefig(os.path.join(path,'S(K)_sdw_compa_U'+str(u)+'_V'+str(v1)+'_na'+str(rho1)+'_L'+str(L)+'.png'),bbox_inches='tight')
