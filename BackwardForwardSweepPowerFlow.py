# -*- coding: utf-8 -*-

import numpy as np
import math
from numpy.lib.function_base import append
import bibliotecaLerDados_SweepMethod
import SweepMethodBib 
import cmath

# Variaveis globais
Sbase = 1000 # KVA
Vbase = 11.00 # kV
#Zbase = ((Vbase*1000)**2)/(Sbase*1000) #Ohms


arquivo = open('G:\\O meu disco\\Biblioteca de Programas em Python\\FluxoMetodoBackwardForwardSweep\\Teste6barras_Sweep.txt')
DadosBarra,DadosLinhas,barra,Sbase,nb,nl,branch, ntd, ntn = bibliotecaLerDados_SweepMethod.LeitorDados.ler_dados(arquivo)

# Impões o valor da tensão base

Vbase = 11.00 # 11 kV

###################################################################################################################
# Importando os parâmetros do sistema em pu
###################################################################################################################

#Zpu, Zbase = SweepMethodBib.putranformer.Z_pu(DadosBarra,DadosLinhas,Sbase,Vbase,nl)

###################################################################################################################
# Inicia Snode, Inode, Ibranch e V0
###################################################################################################################

Inode = np.zeros(nb-1,dtype=complex)
Snode = np.zeros(nb-1,dtype=complex)
Ibranch = np.zeros(nb-1,dtype=complex)
V0 = np.ones(nb-1,dtype=complex)

for i in range(nb-1):
    V0[i] = cmath.rect((DadosBarra['tensao'][i]),DadosBarra['angulo'][i])  #transformando as tensões na forma polar para retangular
for i in range(nb-1):
    Snode[i] =  ((DadosBarra['Pger'][i+1] - DadosBarra['Pcarga'][i+1])) + 1j*((DadosBarra['Qger'][i+1]-DadosBarra['Qcarga'][i+1])) #Preciso verificar se não preciso dividir ou retirars as potências pelo valor do Sbase
for i in range(nb-1):
    Inode[i] = np.conj((Snode[i]/V0[i]))

Vnode = np.copy(V0)

#print(V0)
#print('^=+'*30)
#print(Snode)
#print('^=+'*30)
#print(len(Inode))

###################################################################################################################
# Criação de ponteiros para a etapa backward
###################################################################################################################

PositionBackward = []
accumulator = [] 
for n in range(nl):
    accumulator = [] 
    for m in range(nl):
        if DadosLinhas['Para'][n] == DadosLinhas['De'][m]:
            accumulator.append(m)
    PositionBackward.append(accumulator) 

#print(PositionBackward)
###################################################################################################################
# Inicializa as variaveis de controle do processo iterativo
###################################################################################################################

k = 0
contK = []
mismatche = []
contK.append(k)
E = 0.0001
dic_deltVnode = []
cop = np.copy(V0)
dic_deltVnode.append(cop)
cop = np.copy(V0)
MAX_DELTVBUS = 1
convergence_report = {'iteration':None,'mismatche':None}
mismatche.append(MAX_DELTVBUS)
convergence_report['iteration'] = contK
convergence_report['mismatche'] = mismatche

while MAX_DELTVBUS > E and k < 100:

    #Calculo das correntes correntes nos nós

    for i in range(nb-1):
        Inode[i] = np.conj((Snode[i]/Vnode[i]))
        
    # Etapa backward para calculo das correntes nos ramos

    summatory = 0
    for n in reversed(range(nl)):
        for m in range(len(PositionBackward[n])):
            summatory = summatory + Ibranch[PositionBackward[n][m]]
        Ibranch[n] = Inode[n] + summatory
        summatory = 0 

    # Etapa Forward para atualização das tensões nos nós

    for i in range(nl):
        if i == 0:
            Vnode[i] = (V0[i]) + ((((DadosLinhas['R'][i]) + (1j*DadosLinhas['Xl'][i]))*Ibranch[i])/1000)
        else:
            Vnode[DadosLinhas['Para'][i]-1] = (Vnode[DadosLinhas['De'][i]-1]) + ((((DadosLinhas['R'][i]) + (1j*DadosLinhas['Xl'][i]))*Ibranch[i])/1000)
        #print(Vnode[i])

    cop = np.copy(Vnode)
    dic_deltVnode.append(cop)

    # Incrementa a iteração

    k = k + 1
    contK.append(k)

    # Calcula os erros

    MAX_DELTVBUS = np.abs(np.max(((np.abs(dic_deltVnode[k]) - (np.abs(dic_deltVnode[k-1]))))))
    mismatche.append(MAX_DELTVBUS)
    convergence_report['iteration'] = contK
    convergence_report['mismatche'] = mismatche

#print(dic_deltVnode)
#print(k+1)

###################################################################################################################
# Tranforma as tensões nos nós em fasores
###################################################################################################################

V_mod = []
V_mod_pu = []
theta_graus = []
theta_rad = []

for i in range(nb):
    if i == 0: #Para adicionar o valor da primeira barra do sistema com o valor de 1 pu com o angulo de 0 
        V_mod.append(Vbase)
        V_mod_pu.append(Vbase/Vbase)
        theta_graus.append(0)
        theta_rad.append(0)
    else:
        V_mod.append(np.sqrt((np.real(Vnode[i-1])**2)+np.imag(Vnode[i-1])**2))
        V_mod_pu.append(V_mod[i]/Vbase)
        theta_graus.append(((np.arctan((np.imag(Vnode[i-1])/np.real(Vnode[i-1]))))/np.pi)*180)
        theta_rad.append((np.arctan((np.imag(Vnode[i-1])/np.real(Vnode[i-1])))))


#print(V_mod)
#print('=+'*30)
#print(V_mod_pu)

###################################################################################################################
# Calculo dos Fluxos de Potência e injeções de potência nos nós
###################################################################################################################

Power_Flow = SweepMethodBib.powewrflow_calculation.powerflow(DadosBarra, DadosLinhas, nl, nb, V_mod, theta_rad, branch)

Pk = []
Qk = []
Snode = []

for i in range(nb):
    if i == 0: #Para adicionar o valor da primeira barra do sistema com o valor de 1 pu com o angulo de 0 
        Snode.append(Power_Flow[0]['Pkm'][i] + 1j* Power_Flow[0]['Qkm'][i])
        Pk.append(Power_Flow[0]['Pkm'][i])
        Qk.append(Power_Flow[0]['Qkm'][i])
    else:
        Snode.append(Vnode[i-1]*(np.conj(Inode[i-1])))
        Pk.append(np.real(Snode[i])/1000)
        Qk.append(np.imag(Snode[i])/1000)

#print(Pk)
#print('=+'*30)
#print(Qk)

###################################################################################################################
# Gera o relatório de convergência
###################################################################################################################

path = 'G:\\O meu disco\\Biblioteca de Programas em Python\\FluxoMetodoBackwardForwardSweep\\Relatório_SweepMethod.txt'

SweepMethodBib.report_generation.report(path, arquivo, ntd,ntn, nb, nl, convergence_report, Power_Flow, DadosBarra,DadosLinhas, E, V_mod_pu, V_mod, theta_rad, theta_graus, Pk, Qk)

# Observações:
# Na entrada de dados, o valor das tensões iniciais devem estar no valor da tensão base (20 kV, 11 kV...)
# Os valores de Potência gerada e demandada devem ser colocadas em kW e kVAr 
# Aparentemente, os métodos apresentam inconsisteências quando as potências são muito pequenas (não consigo entender o motivo)



       

    


        