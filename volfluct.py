"""
Authors: A.Rustamov

Analytic formulas for volume fluctuations and for their corrections

R. Holzmann, V. Koch, A. Rustamov, J. Stroth

arXiv:

If you use the code to produce results, we ask you to be fair and cite the above paper!

"""

import sympy as sp
from sympy import simplify
from sympy.printing.cxx import cxxcode

class Volfluct:
    def __init__(self):
        self.snamevolume = []
        self.snamecum_n = []
        self.snamecum_nA = []
        self.snamecum_nB = []
        self.resultscum_N = []
        self.snamecum_N = []
        self.resultscum_W = []
        self.x = sp.Symbol('x')
        self.y = sp.Symbol('y')
        self.f = sp.Function('f')(self.x)
        self.g = sp.Function('g')(self.f)
        self.f2 = sp.Function('f2')(self.x, self.y)
        self.g2 = sp.Function('g2')(self.f2)

    def getCumulant(self, orderA, fact = "fact", name = ["n", "N"]):
        cumname = 'k_'
        if fact == "fact":
            cumname = 'C_'
        for i in range(orderA):
            self.snamevolume.extend(['k_'+str(i+1)+"[W]"])
            self.snamecum_n.extend([cumname+str(i+1)+'['+name[0]+']'])
            self.snamecum_N.extend([cumname+str(i+1)+'['+name[1]+']'])
            self.snamevolume[0] = '<W>'
            self.snamecum_n[0] = '<'+name[0]+'>'
        pcum_volume = sp.symbols(self.snamevolume)
        pcum_n = sp.symbols(self.snamecum_n)
        for j in range(orderA):
            self.resultscum_N.extend([self.g.diff(self.x,j+1)])
            for i in range(orderA, 0, -1):
                self.resultscum_N[j] = self.resultscum_N[j].subs(sp.Derivative(self.g,(self.f,i)), pcum_volume[i-1])
                self.resultscum_N[j] = self.resultscum_N[j].subs(sp.Derivative(self.f,(self.x,i)), pcum_n[i-1])
        return [self.resultscum_N, self.snamecum_N]
    
    def getCrossCumulant(self, orderA, fact = "fact", name1 = ["a", "A"], name2 = ["b","B"]):
        cumname = 'k_'
        covname = 'cov'
        if fact == "fact":
            cumname = 'C_'
        for i in range(orderA*2):
            self.snamevolume.extend(['k_'+str(i+1)+"[W]"])

        for i in range(orderA):
            self.snamecum_nA.extend([cumname+str(i+1)+'['+name1[0]+']'])
            self.snamecum_nB.extend([cumname+str(i+1)+'['+name2[0]+']'])
            for j in range(orderA):
                self.snamecum_n.extend([covname+str(i+1)+str(j+1)+'['+name1[0]+name2[0]+']'])
                self.snamecum_N.extend([covname+str(i+1)+str(j+1)+'['+name1[1]+name2[1]+']'])
        self.snamecum_nA[0] = '<'+name1[0]+'>'
        self.snamecum_nB[0] = '<'+name2[0]+'>' 
        self.snamevolume[0] = '<W>'
        pcum_volume = sp.symbols(self.snamevolume)
        pcum_n = sp.symbols(self.snamecum_n)
        pcum_nA = sp.symbols(self.snamecum_nA)
        pcum_nB = sp.symbols(self.snamecum_nB)
        
        for j in range(orderA):
            for k in range(orderA):
                self.resultscum_N.extend([sp.expand(self.g2.diff(self.x,j+1, self.y,k+1))])

        for j in range(len(self.resultscum_N)):
            for i in range(orderA*2, 0, -1):
                self.resultscum_N[j] = self.resultscum_N[j].subs(sp.Derivative(self.g2,(self.f2,i)), pcum_volume[i-1])
        
        for j in range(len(self.resultscum_N)):
            count = orderA*orderA
            for i in range(orderA, 0, -1):
                for k in range(orderA, 0, -1):
                    count -= 1
                    self.resultscum_N[j] = self.resultscum_N[j].subs(sp.Derivative(self.f2,self.x,i,self.y,k), pcum_n[count])

        for j in range(len(self.resultscum_N)):
            for i in range(orderA, 0, -1):
                for k in range(orderA, 0, -1):
                    self.resultscum_N[j] = self.resultscum_N[j].subs(sp.Derivative(self.f2,self.x,i), pcum_nA[i-1])
                
        for j in range(len(self.resultscum_N)):
            for i in range(orderA, 0, -1):
                for k in range(orderA, 0, -1):
                    self.resultscum_N[j] = self.resultscum_N[j].subs(sp.Derivative(self.f2,self.y,k), pcum_nB[k-1])
    
        return [self.resultscum_N, self.snamecum_N]
    
    def test(self):
        [r1, r2] = self.getCumulant(10, '')
        for i in range(len(r1)):
            print(r2[i], " == ", r1[i])
        
     
    def getCumulant1(self, orderA, fact = "fact", name = ["n", "N"]):
        cumname = 'k_'
        if fact == "fact":
            cumname = 'C_'
        snamevolume = []
        snamecum_n = []
        for i in range(orderA):
            snamevolume.extend(['k_'+str(i+1)+"[W]"])
            snamecum_n.extend([cumname+str(i+1)+'['+name[0]+']'])
            snamecum_N = cumname+str(i+1)+'['+name[1]+']'
        snamevolume[0] = '<W>'
        snamecum_n[0] = '<'+name[0]+'>'
        pcum_volume = sp.symbols(snamevolume)
        pcum_n = sp.symbols(snamecum_n)
        resultscum_N = self.g.diff(self.x, orderA)
        for i in range(orderA, 0, -1):
            resultscum_N = resultscum_N.subs(sp.Derivative(self.g,(self.f,i)), pcum_volume[i-1])
            resultscum_N = resultscum_N.subs(sp.Derivative(self.f,(self.x,i)), pcum_n[i-1])
        return(resultscum_N, snamecum_N)
    
    def getCrossCumulant1(self, orderA, orderB, fact = "fact", name1 = ["a", "A"], name2 = ["b","B"]):
        cumname = 'k_'
        covname = 'k_'
        if fact == "fact":
            cumname = 'C_'
        snamevolume = []
        snamecum_nA = []
        snamecum_nB = []
        snamecum_n = []
        for i in range(orderA+orderB):
            snamevolume.extend(['k_'+str(i+1)+"[W]"])

        for i in range(orderA):
            snamecum_nA.extend([cumname+str(i+1)+'['+name1[0]+']'])
        for j in range(orderB):
            snamecum_nB.extend([cumname+str(j+1)+'['+name2[0]+']'])
        
        for i in range(orderA):
            for j in range(orderB):
                snamecum_n.extend([covname+str(i+1)+str(j+1)+'['+name1[0]+name2[0]+']'])
                
        snamecum_N  =  covname+str(orderA)+str(orderB)+'['+name1[1]+name2[1]+']'
        snamecum_nA[0] = '<'+name1[0]+'>'
        snamecum_nB[0] = '<'+name2[0]+'>' 
        snamevolume[0] = '<W>'
        pcum_volume = sp.symbols(snamevolume)
        pcum_n = sp.symbols(snamecum_n)
        pcum_nA = sp.symbols(snamecum_nA)
        pcum_nB = sp.symbols(snamecum_nB)
        
        
        resultscum_N = sp.expand(self.g2.diff(self.x,orderA, self.y,orderB))
        
        for i in range(orderA+orderB, 0, -1):
            resultscum_N = resultscum_N.subs(sp.Derivative(self.g2,(self.f2,i)), pcum_volume[i-1])
        
        count = orderA*orderB
        for i in range(orderA, 0, -1):
            for k in range(orderB, 0, -1):
                count -= 1
                resultscum_N = resultscum_N.subs(sp.Derivative(self.f2,self.x,i,self.y,k), pcum_n[count])

        for i in range(orderA, 0, -1):
            resultscum_N = resultscum_N.subs(sp.Derivative(self.f2,self.x,i), pcum_nA[i-1])
                
        for k in range(orderB, 0, -1):
            resultscum_N = resultscum_N.subs(sp.Derivative(self.f2,self.y,k), pcum_nB[k-1])
    
        return [resultscum_N, snamecum_N]
    
    def getCumW2(self, orderA, ismix = 0):
        snamevolume = []
        snamecum_N  = [] # kappa[N]
        snamecum_n  = [] # kappa[n]

        snamecum_A  = [] # kappa[N]
        snamecum_a  = [] # kappa[n]

        snamecum_B  = [] # kappa[N]
        snamecum_b  = [] # kappa[n]

        snamecum_M  = [] # C[M]
        snamecum_m  = [] # C[m]
        snamecum_N_bar = [] #kappa_bar
        snamecum_A_bar = [] #kappa_bar
        snamecum_B_bar = [] #kappa_bar
        snamecum_M_bar = [] #C_bar
        snamecummix_n = []
        snamecummix_N = []
        snamecummix_N_bar = []

        for i in range(orderA):
            snamecum_N.extend(['k_'+str(i+1)+'[N]'])
            snamecum_n.extend(['k_'+str(i+1)+'[n]'])
  
            snamecum_A.extend(['k_'+str(i+1)+'[N1]'])
            snamecum_a.extend(['k_'+str(i+1)+'[n1]'])

            snamecum_B.extend(['k_'+str(i+1)+'[N2]'])
            snamecum_b.extend(['k_'+str(i+1)+'[n2]'])

            snamecum_m.extend(["C_"+str(i+1)+'[m]'])
            snamecum_N_bar.extend(['k_'+'bar_'+str(i+1)+'[N]'])
            snamecum_A_bar.extend(['k_'+'bar_'+str(i+1)+'[N1]'])
            snamecum_B_bar.extend(['k_'+'bar_'+str(i+1)+'[N2]'])
            snamecum_M_bar.extend(['C_'+'bar_'+str(i+1)+'[M]'])   
            snamecum_M.extend(['C_'+str(i+1)+'[M]'])
            snamevolume.extend(['k_'+str(i+1)+'[W]'])
            
        snamevolume[0] = '<W>'
        pcum_N = sp.symbols(snamecum_N)
        pcum_n = sp.symbols(snamecum_n)
        pcum_A = sp.symbols(snamecum_A)
        pcum_a = sp.symbols(snamecum_a)
        pcum_B = sp.symbols(snamecum_B)
        pcum_b = sp.symbols(snamecum_b)

        pcum_m = sp.symbols(snamecum_m)
        pcum_N_bar = sp.symbols(snamecum_N_bar)
        pcum_A_bar = sp.symbols(snamecum_A_bar)
        pcum_B_bar = sp.symbols(snamecum_B_bar)
        pcum_M_bar = sp.symbols(snamecum_M_bar)
        pcum_M = sp.symbols(snamecum_M)
        pcum_volume = sp.symbols(snamevolume)
        mean_m, mean_M = sp.symbols('<m>,<M>')
        mean_n, mean_N = sp.symbols('<n>,<N>')
        mean_a,mean_b, mean_A, mean_B = sp.symbols('<n1>, <n2>, <N1>, <N2>')

        
        factcumresults = []
        cumresults = []
        mixcums = []
        volk2results = []
        volk2results.extend([0])

        for i in range(orderA):
            [r,n]=self.getCumulant1(i+1, 'fact', ['m', 'M'])
            factcumresults.extend([r])
            [r,n]=self.getCumulant1(i+1, '', ['n', 'N'])
            cumresults.extend([r])
        if ismix:
            for i in range(orderA):
                for j in range(orderA):
                    if (i+j+2) > orderA:
                        continue
                    snamecummix_n.extend(['k_'+str(i+1)+str(j+1)+'[n1n2]'])
                    snamecummix_N.extend(['k_'+str(i+1)+str(j+1)+'[N1N2]'])
                    snamecummix_N_bar.extend(['mxbar'+str(i+1)+str(j+1)+'[AB]'])
                    [r, n] = self.getCrossCumulant1(i+1,j+1, '', ['n1','N1'], ['n2','N2'])
                    mixcums.extend([r])
        
        pcummix_N_bar = sp.symbols(snamecummix_N_bar)
        pcummix_N = sp.symbols(snamecummix_N)
        pcummix_n = sp.symbols(snamecummix_n)
    
        for i in range(1, orderA):
            a = sp.solve(factcumresults[i] - pcum_M[i], pcum_volume[i])
            volk2results.extend([a[0]])

        for i in range(1, orderA):
            for j in range(orderA-1, 0, -1):
                cumresults[i] = cumresults[i].subs(pcum_volume[j],volk2results[j], simultaneous=True)
            cumresults[i] = cumresults[i].subs(mean_m,mean_M/pcum_volume[0], simultaneous=True)
            cumresults[i] = cumresults[i].subs(mean_n,mean_N/pcum_volume[0], simultaneous=True)        
            for j in range(orderA-1, 0, -1):  
                cumresults[i] = cumresults[i].subs(pcum_volume[0]*pcum_n[j],pcum_N_bar[j], simultaneous=True)
                cumresults[i] = cumresults[i].subs(pcum_volume[0]*pcum_m[j],pcum_M_bar[j], simultaneous=True)

        for i in range(0, len(mixcums)):
            for j in range(len(volk2results)-1, 0, -1):
                mixcums[i] = mixcums[i].subs(pcum_volume[j],volk2results[j], simultaneous=True)
            mixcums[i] = mixcums[i].subs(mean_m,mean_M/pcum_volume[0], simultaneous=True)
            mixcums[i] = mixcums[i].subs(mean_n,mean_N/pcum_volume[0], simultaneous=True)
            mixcums[i] = mixcums[i].subs(mean_a,mean_A/pcum_volume[0], simultaneous=True)
            mixcums[i] = mixcums[i].subs(mean_b,mean_B/pcum_volume[0], simultaneous=True)
            for j in range(len(pcum_N_bar)-1, 0, -1):
                mixcums[i] = mixcums[i].subs(pcum_volume[0]*pcum_n[j],pcum_N_bar[j], simultaneous=True)
                mixcums[i] = mixcums[i].subs(pcum_volume[0]*pcum_m[j],pcum_M_bar[j], simultaneous=True)
                mixcums[i] = mixcums[i].subs(pcum_volume[0]*pcum_a[j],pcum_A_bar[j], simultaneous=True)
                mixcums[i] = mixcums[i].subs(pcum_volume[0]*pcum_b[j],pcum_B_bar[j], simultaneous=True)
            for j in range(len(pcummix_N_bar), 0, -1):
                mixcums[i] = mixcums[i].subs(pcum_volume[0]*pcummix_n[j-1],pcummix_N_bar[j-1], simultaneous=True)      

        kn_bar_mix = []
        kn_corr_mix = []
        deltan_mix = []

        kn_corr_mix_name = []
        deltan_mix_name = []

        for i in range(orderA):
            for j in range(orderA):
                if (i+j+2) > orderA:
                    continue
                kn_corr_mix_name.extend(['k_'+str(i+1)+str(j+1)+'_corr[N1N2]'])
                deltan_mix_name.extend(['Delta'+str(i+1)+str(j+1)])
        
        for i in range(0, len(mixcums)):
            b = sp.solve(mixcums[i]-pcummix_N[i],pcummix_N_bar[i])
            kn_bar_mix.extend([b[0]])
            kn_bar_mix[i] = sp.expand(kn_bar_mix[i])
            for j in range(0, i):
                kn_bar_mix[i] = kn_bar_mix[i].subs(pcummix_N_bar[j],kn_bar_mix[j])
            tesrCorrMix =  kn_bar_mix[i]
            for j in range(0, len(pcum_M_bar)):
                tesrCorrMix = tesrCorrMix.subs(pcum_M_bar[j],0., simultaneous=True)
            kn_corr_mix.extend([tesrCorrMix])
            deltan_mix.extend([kn_bar_mix[i] - kn_corr_mix[i]])
            kn_corr_mix[i] = sp.expand(sp.simplify(kn_corr_mix[i]))
            deltan_mix[i] = sp.expand(sp.simplify(deltan_mix[i]))
        kn_bar = []
        kn_corr = []
        deltan = []
        kn_bar.extend([0])
        kn_corr.extend([0])
        deltan.extend([0])

        kn_corr_name = []
        deltan_name = []

        kn_corr_name.extend([0])
        deltan_name.extend([0])

        for i in range(1, orderA):
            kn_corr_name.extend(["k_corr["+str(i+1)+"]"])
            deltan_name.extend(["delta["+str(i+1)+"]"])
            b = sp.solve(cumresults[i]-pcum_N[i],pcum_N_bar[i])
            kn_bar.extend([b[0]])
            kn_bar[i] = sp.expand(kn_bar[i])
            for j in range(1, i):
                kn_bar[i] = kn_bar[i].subs(pcum_N_bar[j],kn_bar[j])
            tesrCorr =  kn_bar[i]
            for j in range(1, orderA):
                tesrCorr = tesrCorr.subs(pcum_M_bar[j],0, simultaneous=True)
            kn_corr.extend([tesrCorr])
            deltan.extend([kn_bar[i] - kn_corr[i]])
            kn_corr[i] = sp.expand(sp.simplify(kn_corr[i]))
            deltan[i] = sp.expand(sp.simplify(deltan[i]))

        results = [kn_corr, deltan, kn_corr_mix, deltan_mix]
        resultsname = [kn_corr_name, deltan_name, kn_corr_mix_name, deltan_mix_name]
        return(results, resultsname)
  



        

