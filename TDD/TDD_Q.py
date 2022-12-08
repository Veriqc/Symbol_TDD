import numpy as np
from TDD.TN import Index,Tensor,TensorNetwork
from qiskit.quantum_info.operators import Operator
from sympy import *
from qiskit.circuit.parameterexpression import ParameterExpression

def is_diagonal(U):
    i, j = np.nonzero(U)
    return np.all(i == j)

def add_hyper_index(var_list,hyper_index):
    for var in var_list:
        if not var in hyper_index:
            hyper_index[var]=0
            
def reshape(M):
    v=M.shape[0]
    h=M.shape[1]
    M=np.array(np.vsplit(M,v/2))
    M=np.array(np.dsplit(M,h/2))
    if M.shape[0]!=2 and M.shape[1]!=2:
        M=reshape(M)
    # else:
    #     M=M.T
    return M

    # U=M
    # if U.shape==(1,1):
    #     return U
    
    # if U.shape[0]==U.shape[1]:
    #     split_U=np.split(U,2,1)
    # else:
    #     split_U=np.split(U,2,0)
    # split_U[0]=reshape(split_U[0])
    # split_U[1]=reshape(split_U[1]) 
    # return np.array([split_U])[0]     
            
def get_real_qubit_num(cir):
    """Calculate the real number of qubits of a circuit"""
    gates=cir.data
    q=0
    for k in range(len(gates)):
        q=max(q,max([qbit.index for qbit in gates[k][1]]))
    return q+1

def cir_2_tn(cir,input_s=[],output_s=[],cong=False):
    """return the dict that link every quantum gate to the corresponding index"""
   
    hyper_index=dict()
    qubits_index = dict()
    start_tensors= dict()
    end_tensors = dict()
    
    qubits_num=get_real_qubit_num(cir)

    for k in range(qubits_num):
        qubits_index[k]=0
        
    tn=TensorNetwork([],tn_type='cir',qubits_num=qubits_num)

        
    if input_s:
        U0=np.array([1,0])
        U1=np.array([0,1])
        for k in range(qubits_num):
            if input_s[k]==0:
                ts=Tensor(U0,[Index('x'+str(k))],'in',[k])
            elif input_s[k]==1:
                ts=Tensor(U1,[Index('x'+str(k))],'in',[k])
            else:
                print('Only support computational basis input')
            tn.tensors.append(ts)
                
    gates=cir.data
    for k in range(len(gates)):
        from qiskit import version 
        if tuple(int(a) for a in version.get_version_info().split('.')) >= (0, 21, 0):
            g=(gates[k].operation,gates[k].qubits)
        else:
            g=gates[k]
        nam=g[0].name
        # print(g,'\n',nam)
        q = [q.index for q in g[1]]
        
        var=[]
        # print(g,nam,q)

        if nam=='cx':
#             print(Operator(g[0]).data)
            var_con='x'+ str(q[0])+'_'+str(qubits_index[q[0]])
            var_tar_in='x'+ str(q[1])+'_'+str(qubits_index[q[1]])
            var_tar_out='x'+ str(q[1])+'_'+str(qubits_index[q[1]]+1)
            add_hyper_index([var_con,var_tar_in,var_tar_out],hyper_index)
            var+=[Index(var_con,hyper_index[var_con]),Index(var_con,hyper_index[var_con]+1),Index(var_tar_in,hyper_index[var_tar_in]),Index(var_tar_out,hyper_index[var_tar_out])]
            U=np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
            U=reshape(U)
            ts=Tensor(U,var,nam,q)
            tn.tensors.append(ts)
            if qubits_index[q[0]]==0 and hyper_index[var_con]==0:
                start_tensors[q[0]]=ts
            if qubits_index[q[1]]==0 and hyper_index[var_tar_in]==0:
                start_tensors[q[1]]=ts
            end_tensors[q[0]]=ts
            end_tensors[q[1]]=ts
            hyper_index[var_con]+=1
            qubits_index[q[1]]+=1            
            continue
        ts=Tensor([],[],nam,q)
        
        if len(g[0].params)>0 :
            assert len(g[0].params)==1,'Not Support yet'
            param_expr = g[0].params[0] 
            # print(param_expr,'\n',type(param_expr))
            
            if isinstance(param_expr,ParameterExpression) and len(param_expr.parameters) != 0:
                s=symbols(list(param_expr.parameters)[0].name.replace("[","").replace("]",""))
                sp_expr=sympify(str(param_expr).replace("[","").replace("]",""))
                const=sp_expr.replace(s,0)
                coeff=int(format(float(sp_expr.replace(s,1)-const),".0f"))
                new_expr=s+const/coeff/2

                from sympy.simplify.fu import TR9, TR11
                # sin_expr=nsimplify(TR11(TR9(fourier_series(sin(sp_expr)).truncate())),tolerance=1e-3)
                # cos_expr=nsimplify(TR11(TR9(fourier_series(cos(sp_expr)).truncate())),tolerance=1e-3)

                sin_expr=nsimplify(TR11(TR9(fourier_series(sin(new_expr)).truncate())),tolerance=1e-3)
                cos_expr=nsimplify(TR11(TR9(fourier_series(cos(new_expr)).truncate())),tolerance=1e-3)

                # if nam=='ry':
                #     U=np.array([[cos_expr,-sin_expr],[sin_expr,cos_expr]])
                #     # print(U)
                # if nam=='rz':
                #     U=np.array([[cos_expr-1j*sin_expr,0],[0,cos_expr+1j*sin_expr]])
                #     # print(U)
                #     if cong==True:
                #         U=np.array([[cos_expr+1j*sin_expr,0],[0,cos_expr-1j*sin_expr]])
                # if nam=='rx':
                #     U=np.array([[cos_expr,-1j*sin_expr],[-1j*sin_expr,cos_expr]])
                #     # print(U)
                #     if cong==True:
                #         U=np.array([[cos_expr,1j*sin_expr],[1j*sin_expr,cos_expr]])

                if coeff > 0:
                    sins=sin_expr
                    coss=cos_expr
                if coeff < 0:
                    sins=-sin_expr
                    coss=cos_expr


                if nam=='ry':
                    U=np.linalg.matrix_power(np.array([[coss,-sins],[sins,coss]]),abs(coeff))
                if nam=='rz':
                    U=np.linalg.matrix_power(np.array([[coss-1j*sins,0],[0,coss+1j*sins]]),abs(coeff))
                    if cong==True:
                        U=np.linalg.matrix_power(np.array([[coss+1j*sins,0],[0,coss-1j*sins]]),abs(coeff))
                if nam=='rx':
                    U=np.linalg.matrix_power(np.array([[coss,-1j*sins],[-1j*sins,coss]]),abs(coeff))
                    if cong==True:
                        U=np.linalg.matrix_power(np.array([[coss,1j*sins],[1j*sins,coss]]),abs(coeff))
            
            else:
                U=Operator(g[0]).data 
                # U=g[0].to_matrix()       

        else:
            U=Operator(g[0]).data 
            # U=g[0].to_matrix()

        if is_diagonal(U):
            for k in q:
                var_in='x'+ str(k)+'_'+str(qubits_index[k])
                add_hyper_index([var_in],hyper_index)
                var+=[Index(var_in,hyper_index[var_in]),Index(var_in,hyper_index[var_in]+1)]
                if qubits_index[k]==0 and hyper_index[var_in]==0:
                    start_tensors[k]=ts
                end_tensors[k]=ts             
                hyper_index[var_in]+=1
        else:
            for k in q:
                var_in='x'+ str(k)+'_'+str(qubits_index[k])
                var_out='x'+ str(k)+'_'+str(qubits_index[k]+1)
                add_hyper_index([var_in,var_out],hyper_index)
                var+=[Index(var_in,hyper_index[var_in]),Index(var_out,hyper_index[var_out])]
                if qubits_index[k]==0 and hyper_index[var_in]==0:
                    start_tensors[k]=ts
                end_tensors[k]=ts                
                qubits_index[k]+=1
        if len(q)>1:
            U=reshape(U)
        if qubits_num>1:
            if len(q)==1:
                U=U.T
        ts.data=U
        ts.index_set=var
        tn.tensors.append(ts)

    for k in range(qubits_num):
        if k in start_tensors:
            last1=Index('x'+str(k)+'_'+str(0),0)
            new1=Index('x'+str(k),0)            
            start_tensors[k].index_set[start_tensors[k].index_set.index(last1)]=new1
        if k in end_tensors:
            last2=Index('x'+str(k)+'_'+str(qubits_index[k]),hyper_index['x'+str(k)+'_'+str(qubits_index[k])])
            new2=Index('y'+str(k),0)            
            end_tensors[k].index_set[end_tensors[k].index_set.index(last2)]=new2
               
    for k in range(qubits_num):
        U=np.eye(2)
        if qubits_index[k]==0 and not 'x'+str(k)+'_'+str(0) in hyper_index:
            var_in='x'+str(k)
            var=[Index('x'+str(k),0),Index('y'+str(k),0)]
            ts=Tensor(U,var,'nu_q',[k])
            tn.tensors.append(ts)            

    if output_s:
        U0=np.array([1,0])
        U1=np.array([0,1])
        for k in range(qubits_num):
            if input_s[k]==0:
                ts=Tensor(U0,[Index('y'+str(k))],'out',[k])
            elif input_s[k]==1:
                ts=Tensor(U1,[Index('y'+str(k))],'out',[k])
            else:
                print('Only support computational basis output')
            tn.tensors.append(ts)
    
    all_indexs=[]
    for k in range(qubits_num):
        all_indexs.append('x'+str(k))
        for k1 in range(qubits_index[k]+1):
            all_indexs.append('x'+str(k)+'_'+str(k1))
        all_indexs.append('y'+str(k))

    return tn,all_indexs

def add_inputs(tn,input_s,qubits_num):
    U0=np.array([1,0])
    U1=np.array([0,1])
    if len(input_s)!= qubits_num:
        print("inputs is not match qubits number")
        return 
    for k in range(qubits_num-1,-1,-1):
        if input_s[k]==0:
            ts=Tensor(U0,[Index('x'+str(k))],'in',[k])
        elif input_s[k]==1:
            ts=Tensor(U1,[Index('x'+str(k))],'in',[k])
        else:
            print('Only support computational basis input')
        tn.tensors.insert(0,ts)
            
def add_outputs(tn,output_s,qubits_num):
    U0=np.array([1,0])
    U1=np.array([0,1])
    if len(output_s)!= qubits_num:
        print("outputs is not match qubits number")
        return 
    for k in range(qubits_num):
        if output_s[k]==0:
            ts=Tensor(U0,[Index('y'+str(k))],'out',[k])
        elif output_s[k]==1:
            ts=Tensor(U1,[Index('y'+str(k))],'out',[k])
        else:
            print('Only support computational basis output')
        tn.tensors.append(ts)       

def add_trace_line(tn,qubits_num):
    U=np.eye(2)
    for k in range(qubits_num-1,-1,-1):
        var_in='x'+str(k)
        var=[Index('x'+str(k),0),Index('y'+str(k),0)]
        ts=Tensor(U,var,'tr',[k])
        tn.tensors.insert(0,ts)
        