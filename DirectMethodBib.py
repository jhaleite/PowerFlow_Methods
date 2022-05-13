from os import path
from types import CoroutineType
import numpy as np

#arquivo = open('C:\\Users\\Joao Henrique Leite\\Biblioteca de Programas em Python\\FluxoMétodoDireto\\CIGRE.txt')
#DadosBarra,DadosLinhas,barra,Sbase,nb,nl, branch = LeitorDados.ler_dados(arquivo)

class matrixbiulder: 
    def BIBC(bus_data, branch_data, nl,nb, Sbase):

        BIBC = np.zeros((nl,(nb)),dtype=int) #Essa matriz possui dimensão m x (n-1) onde m é o numero de ramos e n é o numero de barras do sistema
        BIBC[0,0] = 1 

        for k in range(nl): #ramo k 
            BIBC[:,branch_data['Para'][k]] = np.copy(BIBC[:,branch_data['De'][k]])
            BIBC[k,branch_data['Para'][k]] = 1 
        
        BIBC = BIBC[:,1:nb] #Essa matriz possui dimensão m x (n-1) onde m é o numero de ramos e n é o numero de barras do sistema

        return[BIBC]

    def BCBV(bus_data, branch_data, nl, nb, Sbase):

        BCBV = np.zeros((nb,nl),dtype=complex)

        BCBV[0,0] = (branch_data['R'][0]) + (1j*branch_data['Xl'][0])

        for k in range(nl): #ramo k 
            BCBV[branch_data['Para'][k],:] = np.copy(BCBV[branch_data['De'][k],:])
            BCBV[branch_data['Para'][k],k] = (branch_data['R'][k]) + (1j*branch_data['Xl'][k])
        BCBV = BCBV[1:nb,:] #Essa matriz possui dimensão (m-1) x n onde m é o numero de ramos e n é o numero de barras do sistema

        return[BCBV]

class powewrflow_calculation:
    def powerflow(bus_data, branch_data, nl, nb, V, theta, branch):
        Fluxo = {'Pkm':None,'Pmk':None,'Qkm':None,'Qmk':None,'PerdaAtiva_km':None,'PerdaReativa_km':None}
        temp = [] 
        temp1 = []
        temp2 = []
        temp3 = []
        temp4 = []
        temp5 = []
        for J in range(nl):
            k = int(branch_data['De'][J])
            m = int(branch_data['Para'][J])
            gkm = branch['g'][J]
            bkm = branch['b'][J]
            bkmsh = branch['bsh'][J]
            akm = branch_data['tap'][J]
            fikm = branch_data['fi'][J]
            temp.append((((1/akm)*V[k])**2)*gkm - ((1/akm)*V[k]*V[m]*(gkm*np.cos(theta[k]-theta[m]+fikm) + bkm*np.sin(theta[k]-theta[m]+fikm))))
            temp1.append(((V[m]**2))*gkm - (1/akm)*V[k]*V[m]*(gkm*np.cos(theta[k]-theta[m]+fikm) - bkm*np.sin(theta[k]-theta[m]+fikm)))
            temp2.append(-(((1/akm)*V[k])**2)*(bkm+bkmsh) - (1/akm)*V[k]*V[m]*(gkm*np.sin(theta[k]-theta[m]+fikm) - bkm*np.cos(theta[k]-theta[m]+fikm)))
            temp3.append(-(V[m]**2)*(bkm+bkmsh) + (1/akm)*V[k]*V[m]*(gkm*np.sin(theta[k]-theta[m]+fikm) + bkm*np.cos(theta[k]-theta[m]+fikm)))
            temp4.append(np.abs((np.abs(temp[J])) - np.abs(temp1[J])))
            temp5.append(np.abs((np.abs(temp2[J])) - np.abs(temp3[J])))
            Fluxo['Pkm'] = temp[:]
            Fluxo['Pmk'] = temp1[:]
            Fluxo['Qkm'] = temp2[:]
            Fluxo['Qmk'] = temp3[:]
            Fluxo['PerdaAtiva_km'] = temp4[:]
            Fluxo['PerdaReativa_km'] = temp5[:]
        del(temp,temp1,temp2,temp3,temp4,temp5)

        return[Fluxo]

class report_generation:
    def report(path,arquivo,ntd,ntn, nb, nl, convergence_report, Power_Flow, bus_data, branch_data, E, V_mod_pu, V_mod, theta_rad, theta_graus, Pk, Qk):

        arq = open(path,'w')
        arq.write('=='*60)
        arq.write('\n')
        arq.write('                                               RELATÓRIO DE SAÍDA                    \n')
        arq.write('=='*60)
        arq.write('\n')
        arq.write('                                           Informações Gerais do Sistema                 \n')
        arq.write('--'*60)
        arq.write('\n')
        arq.write('  Arquivo de Entrada: {}\n'.format(arquivo))    
        arq.write('  Numero de Barras: {}\n'.format(nb))
        arq.write('  Numero de Linhas: {}\n'.format(nl))
        arq.write('  Numero de trafos convencionais: {}\n'.format(ntn))
        arq.write('  Numero de trafos defasadores: {}\n'.format(ntd))
        arq.write('=='*60)
        arq.write('\n')
        arq.write('                                                 DADOS DAS BARRAS                               \n')
        arq.write('=='*60)
        arq.write('\n')
        arq.write('  Barra   ')
        arq.write('  Tipo    ')
        arq.write('    Tensão  ')
        arq.write('    Ângulo (º)  ')
        arq.write('    Ângulo (rad) ')
        arq.write('    Pk (MW)  ')
        arq.write('    Qk (MVAr ')
        #arq.write('    dP (kW)  ')
        #arq.write('    dQ (kVAr) ')
        arq.write('\n')
        arq.write('--'*60) 
        arq.write('\n')
        for i in range(nb):
            arq.write('{:^10.4f}{:^10.4f}  {:^10.4f}   {:^10.4f}        {:^10.4f}      {:^10.4f}   {:^10.4f}\n'.format(float(bus_data['num'][i]+1),float(bus_data['tipo'][i]),float(V_mod[i]),float(theta_graus[i]),float(theta_rad[i]),float(Pk[i]),float(Qk[i])))
        arq.write('\n')
        arq.write('=='*60)
        arq.write('\n')
        arq.write('                                                 DADOS DOS RAMOS                               \n')
        arq.write('=='*60)
        arq.write('\n')
        arq.write(' '*69)
        arq.write('                         PERDAS\n')    
        arq.write('   k    ')
        arq.write('      m    ')
        arq.write('      Pkm (MW)   ')
        arq.write('   Pmk (MW)   ')
        arq.write('   Qkm (MVAr)   ')
        arq.write('   Qmk (MVAr)   ')
        arq.write('   Ativa (MW)   ')
        arq.write('   Reativa (MVAr)   ')
        arq.write('\n')
        arq.write('--'*60)
        arq.write('\n') 
        for i in range(nl):
            arq.write('{:^10.4f}{:^10.4f}    {:^10.4f}   {:^10.4f}      {:^10.4f}      {:^10.4f}     {:^10.4f}        {:^10.4f}\n'.format(float((branch_data['De'][i])+1),float((branch_data['Para'][i])+1),float(Power_Flow[0]['Pkm'][i]),float(Power_Flow[0]['Pmk'][i]),float(Power_Flow[0]['Qkm'][i]),float(Power_Flow[0]['Qmk'][i]),float(Power_Flow[0]['PerdaAtiva_km'][i]),float(Power_Flow[0]['PerdaReativa_km'][i])))
        arq.write('                                                                                      ------            ------')
        arq.write('\n')
        arq.write('                                                                             Total')
        arq.write('  {:^10.4f}        {:^10.4f}'.format(float(sum(Power_Flow[0]['PerdaAtiva_km'])),float(sum(Power_Flow[0]['PerdaReativa_km']))))
        arq.write('\n')
        arq.write('--'*60) 
        arq.write('\n')
        arq.write('=='*60)
        arq.write('\n')
        arq.write('                                               RELATÓRIO DE CONVERGÊNCIA                    \n')
        arq.write('=='*60)
        arq.write('\n')
        arq.write('  Direct Method\n')
        arq.write('  Convergence Tolerance: {:0.4f}\n'.format(float(E)))
        arq.write('--'*60) 
        arq.write('\n')
        arq.write('  Iteration  ')
        arq.write('  Mismatche  ')
        arq.write('\n')
        arq.write('--'*60) 
        arq.write('\n')
        for i in range(len(convergence_report['iteration'])):
            arq.write('    {:^5.4f}      {:^5.4f} \n'.format(float(convergence_report['iteration'][i]),float(convergence_report['mismatche'][i])))
        arq.write('\n')
        arq.write('=='*60)
        arq.write('\n')
        arq.write('  FIM')
        arq.close()

        return[]
