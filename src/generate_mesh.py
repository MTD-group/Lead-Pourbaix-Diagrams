import numpy as np
import math

#constants
R=.008314 #kJ/mol-K
T=298.15 #K
F=96.485 #kJ/V-g

#Gibbs Free Energy from PBEsol, kJ/mol
GPb = 0
GPbO = -182.179911180077
GPbO2 = -225.05892068112
GPb2O3 = -409.660823950954
GPb3O4 = -606.321927722997
PbCO3 = -636.8927823795264
hc = -1680.6559541150234
pbphos = -2265.1179660491134

#experimental Gibbs Free Energy values for ions, kJ/mol
GPb4 = 302.50074
GPb2 = -23.9743
GHPbO2 = -338.904
GPbO3 = -277.56656
GPbO4 = -282.08946
GPbH2 = 290.788
OH = -157.33514
CO3 = -527.97896
H2CO3 = -623.2068
HCO3 = -586.93152

#set up pH range (consider doing this in input)
pHstart = -2
pHend = 16
dpH = 0.045
pHvec = np.arange(pHstart, pHend, dpH)

#set up potential range(consider doing this in input)
Vstart = 4
Vend = -2
dV = -0.045
Vvec = np.arange(Vstart, Vend, dV)

#set up grid
pH_, V_ = np.meshgrid(pHvec, Vvec)
Z = np.empty(np.shape(pH_), dtype = object)

def generate(conc = 1.5e-8, carbon = False, phosphate = False):
    print("calculating stable species")
    print("Pb concentration: ", conc)
    i = 0
    j = 0
    for pH in pHvec:
        i = 0
        for V in Vvec:
            lowpot = 10000000
            stable = "check"
            ue = -F * V
            uH = -R*T*math.log(10)*pH
            uPb = GPb
            uH2O = -237.18
            
            #Pb
            lowpot = uPb
            stable = "Pb"
            
            #Pb2+
            pot = GPb2 + R*T*math.log(conc)
            urxn = pot + 2*ue - uPb
            #print("Pb2+: ", urxn)
            #pb2plus[j] == urxn
            if (urxn <= lowpot) :
                lowpot = urxn
                stable = "Pb++"
                
            #PbO
            urxn = GPbO + 2*ue + 2*uH - uPb - uH2O
            #print("PbO: ", urxn)
            if (urxn <= lowpot) :
                lowpot = urxn
                stable = "PbO"
                
            #PbO2
            urxn = GPbO2 + 4*ue + 4*uH - uPb - 2*uH2O
            if (urxn <= lowpot) :
                lowpot = urxn
                stable = "PbO2"
                
            #Pb3O4
            urxn = (GPb3O4 + 8*ue + 8*uH - 3*uPb -4*uH2O)/3
            if (urxn <= lowpot) :
                lowpot = urxn
                stable = "Pb3O4"
                
            #Pb4+
            pot = GPb4 + R*T*math.log(conc)
            urxn = pot + 4*ue - uPb
            if (urxn <= lowpot) :
                lowpot = urxn
                stable = "Pb++++"
                
            #HPbO2
            pot = GHPbO2 + R*T*math.log(conc)
            urxn = pot + 3*uH + 2*ue - 2*uH2O - uPb
            #print("HPbO2: ", urxn)
            if (urxn <= lowpot) : 
                lowpot = urxn
                stable = "HPbO2-"
                
            #PbO3--
            pot = GPbO3 + R*T*math.log(conc)
            urxn = pot + 6*uH + 4*ue - 3*uH2O - uPb
            #print("PbO3: ", urxn)
            if (urxn <= lowpot) :
                lowpot = urxn
                stable = "PbO3--"
                
            if (carbon):
                #Hc   
                urxn = (hc + 2*uH + 6*ue - 2*uH2O - 2*(CO3 + R*T*math.log(carbon)))/3 
                #print(urxn)
                if (urxn <= lowpot):
                    lowpot = urxn
                    stable = "Hc"
                
                #PbCO3

                border1 = 6.35
                border2 = 10.33
                if pH < border1:
                    urxn = PbCO3 + 2*ue + 2*uH - (H2CO3 + R*T*math.log(carbon))
                elif pH < border2:
                    urxn = PbCO3 + 2*ue + uH - (HCO3 + R*T*math.log(carbon))
                else:
                    urxn = PbCO3 + 2*ue - (CO3 + R*T*math.log(carbon))
                #print("PbCO3: ",urxn )
                if (urxn <= lowpot) :
                    lowpot = urxn
                    stable = "PbCO3"
                    
            if (phosphate): 
                #ortho
                
                h3po4 = -1147.2 + R*T*math.log(phosphate)
                h2po4 = -1135.119 + R*T*math.log(phosphate)
                hpo4 = -1094.116 + R*T*math.log(phosphate)
                po4 = -1025.498 + R*T*math.log(phosphate)
                
                border1 = 2.03
                border2 = 7.19
                border3 = 12.03
                
                if (pH<border1):
                    urxn = (pbphos + 6*uH + 6*ue - 2*h3po4)/3
                elif (pH<border2):
                    urxn = (pbphos + 4*uH + 6*ue - 2*h2po4)/3
                elif (pH<border3):
                    urxn = (pbphos + 2*uH + 6*ue - 2*hpo4)/3
                else:
                    urxn = (pbphos + 6*ue - 2*po4)/3
                    
                if (urxn <= lowpot): 
                    lowpot = urxn
                    stable = "Pb3(PO4)2"
                
            Z[i,j] = stable
            i+=1
        j+=1
    return pH_, V_, Z