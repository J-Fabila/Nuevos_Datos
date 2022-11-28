import numpy as np
import pandas as pd

Nat=36
Nconf=4457
fuerzas=pd.read_csv("fuerzas.csv")
posiciones=pd.read_csv("posiciones.csv")
energias=pd.read_csv("energias_torch.csv")
simbolos=pd.read_csv("atomos_simbolos.csv")

f_np=fuerzas.to_numpy()
e_np=energias.to_numpy()
p_np=posiciones.to_numpy()
s_np=simbolos.to_numpy()
#f_np[:,0]

f_listo=f_np.reshape((Nconf, Nat, 3),order='C')
p_listo=p_np.reshape((Nconf, Nat, 3),order='C')
e_listo=e_np.reshape((Nconf,1),order='C')
simbolos_listo=s_np.reshape((Nat),order='C')
simbolos_listo=simbolos_listo.astype('uint8')

np.savez("Ti18C18", E=e_listo,F=f_listo,R=p_listo,z=simbolos_listo)

