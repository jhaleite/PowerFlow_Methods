from os import path
import numpy as np
from numpy.lib.function_base import append
from bibliotecaLerDados_DirectMethod import LeitorDados
import cmath
import DirectMethodBib                  # Importa a classe com as funções que montam as matrizes do metodo direto e o calculo do fluxo de potência e as perdas

arquivo = open('G:\\O meu disco\\Biblioteca de Programas em Python\\FluxoMétodoDireto\\Teste6barras_DirectMethod.txt')
DadosBarra,DadosLinhas,barra,Sbase,nb,nl,branch, ntd, ntn = LeitorDados.ler_dados(arquivo)
Vbase = 11.00 # 11 kV

###################################################################################################################
# Importando as matrizes BIBC e BCBV
###################################################################################################################

BIBC_c = DirectMethodBib.matrixbiulder.BIBC(DadosBarra,DadosLinhas,nl,nb,Sbase) # Matriz Bus injection to branch current
BCBV_c = DirectMethodBib.matrixbiulder.BCBV(DadosBarra,DadosLinhas,nl,nb,Sbase) # Matriz branch current to bus voltage
DLF_temp = BCBV_c@BIBC_c[0]
DLF_c = DLF_temp[0]                                                                  # Matriz Distribution Load Flow


###################################################################################################################
# Inicia Snode, Inode e V0
###################################################################################################################
Inode = np.zeros(nb-1,dtype=complex)
Snode = np.zeros(nb-1,dtype=complex)
V0 = np.ones(nb-1,dtype=complex)

for i in range(nb-1):
    V0[i] = cmath.rect((DadosBarra['tensao'][i]),DadosBarra['angulo'][i])  #transformando as tensões na forma polar para retangular
for i in range(nb-1):
    Snode[i] =  ((DadosBarra['Pger'][i+1] - DadosBarra['Pcarga'][i+1])) + 1j*((DadosBarra['Qger'][i+1]-DadosBarra['Qcarga'][i+1])) #Preciso verificar se não preciso dividir ou retirars as potências pelo valor do Sbase
for i in range(nb-1):
    Inode[i] = np.conj((Snode[i]/V0[i]))

Ibranch = BIBC_c@Inode
###################################################################################################################
# Calculando as correntes nos ramos Ibranch e a diferença de tensão nas barras deltVbus mais as tensões nas barras
###################################################################################################################

deltVnode = DLF_c@Inode.T
deltVnode = deltVnode/1000 #Para transformar o resultado em kV
Vnode = V0 + deltVnode

#print('+='*30)
#print(Vbus)

###############################################################
# Inicializa as variaveis de controle do processo iterativo
###############################################################

k = 0
it = 100
contK = []
mismatche = []
contK.append(k)
E = 0.0001
dic_deltVnode = []
dic_deltVnode.append(Vnode)
MAX_DELTVBUS = 1
convergence_report = {'iteration':None,'mismatche':None}
mismatche.append(MAX_DELTVBUS)
convergence_report['iteration'] = contK
convergence_report['mismatche'] = mismatche


while MAX_DELTVBUS > E and k < 500:

    # 1 - Calcula o valor das injeções de correntes nos barramentos

    for i in range(nb-1):
        Inode[i] = np.conj((Snode[i]/Vnode[i]))
    #print(Inode)
    Ibranch = BIBC_c@Inode
    # 2 - Calcula o valor dos mismatches de tensão nos barramentos

    deltVnode = DLF_c@Inode.T
    deltVnode = deltVnode/1000 #tranforma para kV

    # 3 - Atualiza o vetor das tensões nos barramentos

    Vnode = V0 + deltVnode
    dic_deltVnode.append(Vnode)

    # 4 - Incrementa a iteração

    k += 1
    contK.append(k)

    # 5 - Calcula o valor do erro

    MAX_DELTVBUS = np.abs(np.max(((np.abs(dic_deltVnode[k]) - (np.abs(dic_deltVnode[k-1]))))))
    mismatche.append(MAX_DELTVBUS)
    convergence_report['iteration'] = contK
    convergence_report['mismatche'] = mismatche


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

Power_Flow = DirectMethodBib.powewrflow_calculation.powerflow(DadosBarra, DadosLinhas, nl, nb, V_mod, theta_rad, branch)

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

path = 'G:\\O meu disco\Biblioteca de Programas em Python\\FluxoMétodoDireto\\Relatório_DirectMethod.txt'

DirectMethodBib.report_generation.report(path, arquivo, ntd,ntn, nb, nl, convergence_report, Power_Flow, DadosBarra,DadosLinhas, E, V_mod_pu, V_mod, theta_rad, theta_graus, Pk, Qk)

# Observações:
# Na entrada de dados, o valor das tensões iniciais devem estar no valor da tensão base (20 kV, 11 kV...)
# Os valores de Potência gerada e demandada devem ser colocadas em kW e kVAr 
# Aparentemente, os métodos apresentam inconsisteências quando as potências são muito pequenas (não consigo entender o motivo)
