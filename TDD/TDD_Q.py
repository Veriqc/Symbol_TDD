import numpy as np
from TDD.TN import Index,Tensor,TensorNetwork
from qiskit.quantum_info.operators import Operator
from sympy import *
from qiskit.circuit.parameterexpression import ParameterExpression
import math


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

# def cir_2_tn(cir,input_s=[],output_s=[],cong=False):
#     """return the dict that link every quantum gate to the corresponding index"""
   
#     hyper_index=dict()
#     qubits_index = dict()
#     start_tensors= dict()
#     end_tensors = dict()
    
#     qubits_num=get_real_qubit_num(cir)

#     for k in range(qubits_num):
#         qubits_index[k]=0
        
#     tn=TensorNetwork([],tn_type='cir',qubits_num=qubits_num)

        
#     if input_s:
#         U0=np.array([1,0])
#         U1=np.array([0,1])
#         for k in range(qubits_num):
#             if input_s[k]==0:
#                 ts=Tensor(U0,[Index('x'+str(k))],'in',[k])
#             elif input_s[k]==1:
#                 ts=Tensor(U1,[Index('x'+str(k))],'in',[k])
#             else:
#                 print('Only support computational basis input')
#             tn.tensors.append(ts)
                
#     gates=cir.data
#     for k in range(len(gates)):
#         from qiskit import version 
#         if tuple(int(a) for a in version.get_version_info().split('.')) >= (0, 21, 0):
#             g=(gates[k].operation,gates[k].qubits)
#         else:
#             g=gates[k]
#         nam=g[0].name
#         # print(g,'\n',nam)
#         q = [q.index for q in g[1]]
        
#         var=[]
#         # print(g,nam,q)

#         if nam=='cx':
# #             print(Operator(g[0]).data)
#             var_con='x'+ str(q[0])+'_'+str(qubits_index[q[0]])
#             var_tar_in='x'+ str(q[1])+'_'+str(qubits_index[q[1]])
#             var_tar_out='x'+ str(q[1])+'_'+str(qubits_index[q[1]]+1)
#             add_hyper_index([var_con,var_tar_in,var_tar_out],hyper_index)
#             var+=[Index(var_con,hyper_index[var_con]),Index(var_con,hyper_index[var_con]+1),Index(var_tar_in,hyper_index[var_tar_in]),Index(var_tar_out,hyper_index[var_tar_out])]
#             U=np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]])
#             U=reshape(U)
#             ts=Tensor(U,var,nam,q)
#             tn.tensors.append(ts)
#             if qubits_index[q[0]]==0 and hyper_index[var_con]==0:
#                 start_tensors[q[0]]=ts
#             if qubits_index[q[1]]==0 and hyper_index[var_tar_in]==0:
#                 start_tensors[q[1]]=ts
#             end_tensors[q[0]]=ts
#             end_tensors[q[1]]=ts
#             hyper_index[var_con]+=1
#             qubits_index[q[1]]+=1            
#             continue
#         ts=Tensor([],[],nam,q)



#         '''
#         TrDD part start
#         '''
#         # if len(g[0].params)>0 :
#         if g[0].is_parameterized():
#             assert len(g[0].params)==1,'Not Support yet'
#             param_expr = g[0].params[0] 
#             # print(param_expr,'\n',type(param_expr))
            
#             if isinstance(param_expr,ParameterExpression) and len(param_expr.parameters) != 0:
#                 '''
#                 由於旋轉矩陣的三角函數為sin((coeff*theta+const)/2)，分數的角度不好處理，因此實際上是把theta/2當作一個變數。
#                 下方取出const與coeff，再求出新的表達式new_expr。
#                 '''
#                 s=symbols(list(param_expr.parameters)[0].name.replace("[","").replace("]",""))
#                 sp_expr=sympify(str(param_expr).replace("[","").replace("]",""))
                
#                 #取得旋轉角度所加的常數
#                 const=sp_expr.replace(s,0)
#                 #取得旋轉角度變數所加的係數
#                 # coeff=int(format(float(sp_expr.replace(s,1)-const),".0f"))
#                 coeff=sp_expr.replace(s,1)-const
#                 coeff_self=int(format(float(sp_expr.replace(s,1)-const),".0f"))
#                 coeff_plus=int(format(float(sp_expr.replace(s,1)-const),".0f"))+1
#                 if math.isclose(coeff , coeff_self, rel_tol = 1e-3):
#                     coeff = coeff_self
#                 elif math.isclose(coeff , coeff_plus, rel_tol = 1e-3):
#                     coeff = coeff_plus
#                 elif math.isclose(coeff+1 , coeff_plus, rel_tol = 1e-3):
#                     coeff = coeff_plus-1
#                 else:
#                     print(sp_expr)
#                     print(const,coeff)
#                     print(coeff_self,coeff_plus)
#                     assert isinstance(coeff,int), 'coeff is %s and the value is %i not int'%(type(coeff),coeff)
#                 #求出新的變數表達式
#                 # new_expr=s+const/coeff/2
#                 new_expr=coeff*s+const/2

#                 # from sympy.simplify.fu import TR9, TR11
#                 # # sin_expr=nsimplify(TR11(TR9(fourier_series(sin(sp_expr)).truncate())),tolerance=1e-3)
#                 # # cos_expr=nsimplify(TR11(TR9(fourier_series(cos(sp_expr)).truncate())),tolerance=1e-3)

#                 # sin_expr=nsimplify(TR11(TR9(fourier_series(sin(new_expr)).truncate())),tolerance=1e-3)
#                 # cos_expr=nsimplify(TR11(TR9(fourier_series(cos(new_expr)).truncate())),tolerance=1e-3)

#                 sin_expr=sin(new_expr).expand(trig=True)
#                 cos_expr=cos(new_expr).expand(trig=True)

#                 # if coeff > 0:
#                 #     sins=sin_expr
#                 #     coss=cos_expr
#                 # if coeff < 0:
#                 #     sins=-sin_expr
#                 #     coss=cos_expr

#                 if nam=='ry':
#                     U=np.array([[cos_expr,-sin_expr],[sin_expr,cos_expr]])
#                 if nam=='rz':
#                     U=np.array([[cos_expr-1j*sin_expr,0],[0,cos_expr+1j*sin_expr]])
#                     if cong==True:
#                         U=np.array([[cos_expr+1j*sin_expr,0],[0,cos_expr-1j*sin_expr]])
#                 if nam=='rx':
#                     U=np.array([[cos_expr,-1j*sin_expr],[-1j*sin_expr,cos_expr]])
#                     if cong==True:
#                         U=np.array([[cos_expr,1j*sin_expr],[1j*sin_expr,cos_expr]])

#                 # if nam=='ry':
#                 #     U=np.linalg.matrix_power(np.array([[coss,-sins],[sins,coss]]),abs(coeff))
#                 # if nam=='rz':
#                 #     U=np.linalg.matrix_power(np.array([[coss-1j*sins,0],[0,coss+1j*sins]]),abs(coeff))
#                 #     if cong==True:
#                 #         U=np.linalg.matrix_power(np.array([[coss+1j*sins,0],[0,coss-1j*sins]]),abs(coeff))
#                 # if nam=='rx':
#                 #     U=np.linalg.matrix_power(np.array([[coss,-1j*sins],[-1j*sins,coss]]),abs(coeff))
#                 #     if cong==True:
#                 #         U=np.linalg.matrix_power(np.array([[coss,1j*sins],[1j*sins,coss]]),abs(coeff))
            
#             else:
#                 # U=Operator(g[0]).data 
#                 U=g[0].to_matrix()       
#             '''
#             TrDD part end
#             '''


#         else:
#             # U=Operator(g[0]).data 
#             U=g[0].to_matrix()

#         if is_diagonal(U):
#             for k in q:
#                 var_in='x'+ str(k)+'_'+str(qubits_index[k])
#                 add_hyper_index([var_in],hyper_index)
#                 var+=[Index(var_in,hyper_index[var_in]),Index(var_in,hyper_index[var_in]+1)]
#                 if qubits_index[k]==0 and hyper_index[var_in]==0:
#                     start_tensors[k]=ts
#                 end_tensors[k]=ts             
#                 hyper_index[var_in]+=1
#         else:
#             for k in q:
#                 var_in='x'+ str(k)+'_'+str(qubits_index[k])
#                 var_out='x'+ str(k)+'_'+str(qubits_index[k]+1)
#                 add_hyper_index([var_in,var_out],hyper_index)
#                 var+=[Index(var_in,hyper_index[var_in]),Index(var_out,hyper_index[var_out])]
#                 if qubits_index[k]==0 and hyper_index[var_in]==0:
#                     start_tensors[k]=ts
#                 end_tensors[k]=ts                
#                 qubits_index[k]+=1
#         if len(q)>1:
#             U=reshape(U)
#         if qubits_num>1:
#             if len(q)==1:
#                 U=U.T
#         ts.data=U
#         ts.index_set=var
#         tn.tensors.append(ts)

#     for k in range(qubits_num):
#         if k in start_tensors:
#             last1=Index('x'+str(k)+'_'+str(0),0)
#             new1=Index('x'+str(k),0)            
#             start_tensors[k].index_set[start_tensors[k].index_set.index(last1)]=new1
#         if k in end_tensors:
#             last2=Index('x'+str(k)+'_'+str(qubits_index[k]),hyper_index['x'+str(k)+'_'+str(qubits_index[k])])
#             new2=Index('y'+str(k),0)            
#             end_tensors[k].index_set[end_tensors[k].index_set.index(last2)]=new2
               
#     for k in range(qubits_num):
#         U=np.eye(2)
#         if qubits_index[k]==0 and not 'x'+str(k)+'_'+str(0) in hyper_index:
#             var_in='x'+str(k)
#             var=[Index('x'+str(k),0),Index('y'+str(k),0)]
#             ts=Tensor(U,var,'nu_q',[k])
#             tn.tensors.append(ts)            

#     if output_s:
#         U0=np.array([1,0])
#         U1=np.array([0,1])
#         for k in range(qubits_num):
#             if input_s[k]==0:
#                 ts=Tensor(U0,[Index('y'+str(k))],'out',[k])
#             elif input_s[k]==1:
#                 ts=Tensor(U1,[Index('y'+str(k))],'out',[k])
#             else:
#                 print('Only support computational basis output')
#             tn.tensors.append(ts)
    
#     all_indexs=[]
#     for k in range(qubits_num):
#         all_indexs.append('x'+str(k))
#         for k1 in range(qubits_index[k]+1):
#             all_indexs.append('x'+str(k)+'_'+str(k1))
#         all_indexs.append('y'+str(k))

#     return tn,all_indexs


def cir_2_tn(cir,input_s=[],output_s=[],cong=False):
    """return the dict that link every quantum gate to the corresponding index"""
   
    hyper_index=dict()
    qubits_index = dict()
    start_tensors= dict()
    end_tensors = dict()
    
    # from qiskit import transpile
    # cir=transpile(cir,basis_gates=['cx','u'],optimization_level=0)
    qubits_num=get_real_qubit_num(cir)
    parameters_dict=dict((v,k) for k, v in dict(enumerate(cir.parameters)).items())
    # parameter_num=len(parameters_dict)
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


        if g[0].is_parameterized():
            '''
            待修改：control-U 帶變數的目前這方法還需再修正。
            '''
            from qiskit import QuantumCircuit,transpile
            from qiskit import QuantumRegister
            from qiskit.circuit.quantumregister import Qubit
            
            g2=gates[k].copy()
            g2.qubits=(Qubit(QuantumRegister(1, 'q'), 0),)
            
            temp=QuantumCircuit(1)
            temp.append(g2)
            temp=transpile(temp,basis_gates='u')
            params=temp.data[0].operation.params
            '''
            global phase的問題還是還沒解決
            '''
            
            def u3(theta,phi,lam):
                half_theta= theta/2
                cos_half_theta=exp(str(1j*half_theta))/2+exp(str(-1j*half_theta))/2
                sin_half_theta=exp(str(1j*half_theta))/2j-exp(str(-1j*half_theta))/2j
                # return np.array([[cos_half_theta,               -exp(str(1j*lam))*sin_half_theta],
                #                 [exp(str(1j*phi))*sin_half_theta,      exp(str(1j*(phi+lam)))*cos_half_theta]])
                return np.array([[cos_half_theta,               exp(str(-1j*half_theta+1j*lam))/2j-exp(str(1j*half_theta+1j*lam))/2j],
                                [-exp(str(-1j*half_theta+1j*phi))/2j+exp(str(1j*half_theta+1j*phi))/2j,      exp(str(1j*half_theta+1j*phi+1j*lam))/2+exp(str(-1j*half_theta+1j*phi+1j*lam))/2]])
            U=u3(*params)
            # def extract_expr(param):
            #     if len(param.parameters)==0:
            #         return None, None, float(param)
            #     p=tuple(param.parameters)[0]
            #     const=param.bind({p:0})
            #     const=const.sympify().simplify()
            #     coeff=param.bind({p:1}).sympify()-const
            #     coeff=coeff.simplify()
            #     return float(coeff), p, float(const)
            # expr= [extract_expr(item) for item in params]

            # def construct_exp_expr(coeff,p,const):
            #     r=np.abs(const)
            #     theta = np.angle(const)
            #     temp=[0]*parameter_num+[theta]
            #     if p:
            #         pos=parameters_dict[p]
            #         temp[pos]=coeff
            #     from TDD.Exp.EXP import BDD
            #     return BDD({tuple(temp):r})
            
            # def angle_mul(theta, scale):
            #     return [theta[0]*scale,theta[1],theta[2]*scale]

            # def u3_exp(theta,phi,lam):
            #     half_theta=angle_mul(theta,1j/2)
            #     exp_i_half_theta=construct_exp_expr(*half_theta)
            #     minus_half_theta= angle_mul(theta,-1j/2)
            #     exp_minus_i_half_theta=construct_exp_expr(*minus_half_theta)

            #     from TDD.Exp.EXP import add, mul

            #     cos_half_theta=mul(0.5,add(exp_i_half_theta,exp_minus_i_half_theta))
            #     sin_half_theta=mul(-0.5j,add(exp_i_half_theta,mul(-1,exp_minus_i_half_theta)))

            #     exp_i_lam=construct_exp_expr(*lam)
            #     exp_i_phi=construct_exp_expr(*phi)
                

            #     return np.array([[cos_half_theta,              mul(-1.0,exp_i_lam)*sin_half_theta],
            #                     [exp_i_phi*sin_half_theta,      exp_i_lam*exp_i_phi*cos_half_theta]])
            # U=u3_exp(*expr)
            
        else:
            U=g[0].to_matrix()

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
        