from __future__ import annotations
import sympy as sp
import numpy as np

def generate_qubit_symbols(x_index):
    x = ['']*2
    for i in range(2):
        n_str = 'n' if i == 0 else ''
        x[i] = sp.symbols('x%s%i' % (n_str, x_index))
        
    return x

class PBF():
    def __init__(self, expr,tolerance=6):
        tolerance1=10**-tolerance
        self.expr = sp.nsimplify(expr,tolerance=tolerance1,rational=False)
        
    def __repr__(self):
        # return str(self.expr)
        return r'%s'%sp.latex(self.expr)
    
    def __str__(self):
        return str(self.expr)

    def __eq__(self, g):   # defind == 
        # return self.expr==g.expr
        return self.expr.equals(g.expr)
        
    def get_coeffs_poly(self):
        return [self.expr] if self.is_constant() else sp.Poly(self.expr).coeffs()

    def get_qubit_idx_set(self):
        sym_set = self.expr.free_symbols
        num_set = set()
        for sym in sym_set:
            num = int(sym.name.replace('x','').replace('n',''))
            num_set.add(num)
        return num_set
    
    def is_constant(self):
        return self.expr.is_constant()
    
    def xreplace(self, rule):
        #print(rule)
        return PBF(self.expr.xreplace(rule))
        # return PBF(self.expr.evalf(subs=rule,chop=True))

class NormalForm():
    def __init__(self, pbf: PBF, weight=1, qubit_idx_set =set(),tolerance=6):

        tolerance1=10**-tolerance
        
        self.qubit_idx_set = qubit_idx_set.copy()

        if weight == 0 or pbf.expr.is_zero:
            self.pbf = PBF(0)
            self.weight = 0
            return
        if pbf.is_constant():
            pbf_weight = pbf.expr
        elif len(pbf.expr.args) > 0:
            pbf_weight = pbf.expr.args[0]
            while not pbf_weight.is_constant():
                pbf_weight=NormalForm(PBF(pbf_weight)).weight
        else:
            pbf_weight = sp.Integer(1)
        
        if pbf_weight.is_constant() and not pbf_weight.is_zero:
            self.pbf = PBF(pbf.expr / pbf_weight)
            x=sp.Symbol('x')
            self.weight= sp.nsimplify((weight * pbf_weight)*x,tolerance=tolerance1,rational=False)/x
        else:
            self.pbf = pbf
            self.weight = weight

        
    def __repr__(self): 
        return r'%s' % sp.latex(self.weight*self.pbf.expr)
    
    def __str__(self):
        if self.weight==complex(0):
            return '0' 
        elif self.pbf.expr==1:
             return '%s' % str(self.weight)
        elif self.weight==complex(1):
            return ' %s ' % str(self.pbf)
        else:
            return '%s * [%s] ' % (str(self.weight), str(self.pbf))
    @property
    def _repr_latex_(self):
        from IPython.display import Math
        return Math(r'%s'%self.__repr__())

    def copy(self):
        return NormalForm(self.pbf, self.weight, self.qubit_idx_set.copy())
    
    def is_constant(self):
        return self.pbf.is_constant()

    def is_zero(self):
        
        return self.weight==0 or self.pbf.expr.is_zero
    
    def xreplace(self, rule, x_index):
        new_qubit_idx_set = self.qubit_idx_set.copy()
        if x_index in new_qubit_idx_set:
            new_qubit_idx_set.remove(x_index)
            
        new_pbf = PBF(self.pbf.xreplace(rule))
        if new_pbf.is_constant():
            new_qubit_idx_set = set()
        return NormalForm(new_pbf, weight=self.weight, qubit_idx_set=new_qubit_idx_set)

    def __add__(self, g):
        if self.is_constant() and g.is_constant():
            return NormalForm(PBF(self.weight*self.pbf.expr + g.weight*g.pbf.expr))
        if self.pbf.expr.is_zero:
            return g.copy()
        if g.pbf.expr.is_zero:
            return self.copy()
        if self.pbf==g.pbf:
            return NormalForm(weight=self.weight+g.weight,pbf=self.pbf)

        return NormalForm.op_subroutine(self, g, NormalForm.__add__)
    
    def __mul__(self, g:NormalForm):
        if self.is_constant() or g.is_constant():
            qubit_idx_set = self.qubit_idx_set.union(g.qubit_idx_set)
            return NormalForm(PBF(self.pbf.expr * g.pbf.expr), weight=self.weight*g.weight, qubit_idx_set=qubit_idx_set)
        elif self.pbf.expr.is_zero or g.pbf.expr.is_zero:
            return NormalForm(PBF(0))
        '''
        new start
        '''
        n=min(self.find_fork_idx,g.find_fork_idx)
        if n>0:
            coeff_dict1=self.mul_coeff_dict(n)
            coeff_dict2=g.mul_coeff_dict(n)
            def coeff_dict_mul(dict1,dict2):
                output=1
                dict1=dict1.copy()
                dict2=dict2.copy()
                qubit_idx_set = set(dict1.keys()).union(set(dict2.keys()))
                for x_index in dict1.keys():
                    from TDD.normal_form_v2 import generate_qubit_symbols
                    x = generate_qubit_symbols(x_index)
                    if x_index in dict2.keys():
                        output*=dict1[x_index][0]*dict2[x_index][0]*x[1]+sp.simplify(dict1[x_index][1]*dict2[x_index][1])*x[0]
                        del dict2[x_index]
                    else:
                        output*=dict1[x_index][0]*x[1]+dict1[x_index][1]*x[0]
                for x_index in dict2.keys():
                    x = generate_qubit_symbols(x_index)
                    output*=dict2[x_index][0]*x[1]+dict2[x_index][1]*x[0]
                return PBF(output), qubit_idx_set
            mul=coeff_dict_mul(coeff_dict1,coeff_dict2)
            h=NormalForm(mul[0],weight=self.weight*g.weight,qubit_idx_set=mul[1])
        else:
            coeff_dict1=set()
            coeff_dict2=set()
            h=NormalForm(PBF(self.weight*g.weight))
        # print(h)

        pf=self.remaining_term(coeff_dict1,n)
        pg=g.remaining_term(coeff_dict2,n)
        
        '''mul part start'''
        total_qubit_idx_set = pf.qubit_idx_set.union(pg.qubit_idx_set)

        if len(total_qubit_idx_set)!=0:
            x_index = total_qubit_idx_set.pop()
        else:
            return h

        x = generate_qubit_symbols(x_index)
        
        sub_nf = [0]*2
        ws = [0]*2
        
        for i in range(2):
            sub_nf_f = pf.xreplace({x[i]: 1, x[1-i]: 0}, x_index)
            sub_nf_g = pg.xreplace({x[i]: 1, x[1-i]: 0}, x_index)
            sub_nf[i] = sub_nf_f*sub_nf_g
            ws[i] = sub_nf[i].weight
            sub_nf[i].weight = sp.Integer(1) if not sub_nf[i].pbf.expr.is_zero else sp.Integer(0)            
        
        new_qubit_idx_set = sub_nf[0].qubit_idx_set.union(sub_nf[1].qubit_idx_set)
        new_qubit_idx_set.add(x_index)


        #NormalForm.shannon_expansion(sub_nf, ws, x, new_qubit_idx_set, x_index)
        fs=sub_nf
        xs=x
        qubit_idx_set=new_qubit_idx_set.union(h.qubit_idx_set)
        '''shannon_expansion part'''

        if complex(ws[1]) == complex(0):
            result =  (ws[0]*fs[0].weight, xs[0]*fs[0].pbf.expr)
        elif fs[0].pbf == fs[1].pbf:
            if complex(ws[0]) == complex(ws[1]):
                result =  (ws[0]*fs[0].weight, fs[0].pbf.expr)
                #qubit_idx_set.pop()
                qubit_idx_set.remove(x_index)
            else:
                result =  (ws[1]*fs[0].weight, (xs[1]+ws[0]/ws[1]*xs[0])*fs[0].pbf.expr)
        else:
            result = (ws[1]*fs[1].weight, (xs[1]*fs[1].pbf.expr+ws[0]/ws[1]*fs[0].weight/fs[1].weight*xs[0]*fs[0].pbf.expr))


        result1=h.pbf.expr*result[1]
        result0=h.weight*result[0]
        return NormalForm(PBF(result1), weight=result0, qubit_idx_set=qubit_idx_set)
        # return NormalForm(PBF(result[1]), weight=result[0], qubit_idx_set=qubit_idx_set)
        
        '''
        new end
        '''
        # return NormalForm.op_subroutine(self, g, NormalForm.__mul__)
    
    def __truediv__(self, g):
        if g.pbf.expr.is_zero:
            return self.copy()
        elif g.is_constant():
            new_f = self.copy()
            new_f.weight = self.weight / g.weight
            return new_f

        return NormalForm.op_subroutine(self, g, NormalForm.__truediv__)

    def __eq__(self, g):   # defind '==' 
        return self.weight==g.weight and self.pbf==g.pbf
    
    def normalise(self, g: NormalForm):
        if self.is_zero():
            return (g.copy(), NormalForm(PBF(0)), NormalForm(PBF(1)))
        if self.is_constant():
            return (self.copy(), NormalForm(PBF(1)), g/self)

        total_qubit_idx_set = self.qubit_idx_set.union(g.qubit_idx_set)
        x_index = total_qubit_idx_set.pop()
        x = generate_qubit_symbols(x_index)

        w=[0]*2

        def biggest_coeff(f):
            return  max(f.pbf.get_coeffs_poly(), key=abs)*f.weight
        
        def grab_sub_and_coeff(f, g , rule):
            sub_f=f.xreplace(rule,x_index)
            sub_g_comp=False
            if sub_f.is_zero():
                sub=g.xreplace(rule,x_index)
                sub_g_comp=True
            else:
                sub=sub_f
                
            return sub_f, sub ,biggest_coeff(sub),sub_g_comp

        rule=[{x[0]:1,x[1]:0},{x[0]:0,x[1]:1}]
        
        sub_f0, sub0, w[0], sub_g0_comp = grab_sub_and_coeff(self,g ,rule[0])
        sub_f1, sub1, w[1], sub_g1_comp = grab_sub_and_coeff(self,g ,rule[1])
        
        w[0]=NormalForm(PBF(w[0]))
        w[1]=NormalForm(PBF(w[1]))

        N0=NormalForm(PBF(0))
        Nx=NormalForm.normal_form_init(PBF(x[1]))
        Nxn=NormalForm.normal_form_init(PBF(x[0]))
        
        sub_g0=sub0 if sub_g0_comp else g.xreplace(rule[0],x_index)
        sub_g1=sub1 if sub_g1_comp else g.xreplace(rule[1],x_index)

        (h0,f0,g0)=NormalForm.normalise(sub_f0/w[0],sub_g0/w[0]) if not w[0].is_zero() else (N0,N0,N0) 
        (h1,f1,g1)=NormalForm.normalise(sub_f1/w[1],sub_g1/w[1]) #w1=0 will??
        
        h = w[0]*Nxn*h0+w[1]*Nx*h1
        f = Nxn*f0+Nx*f1
        g = Nxn*g0+Nx*g1
        return (h ,f ,g )

    @staticmethod
    def op_subroutine(f, g, opfunc):
        total_qubit_idx_set = f.qubit_idx_set.union(g.qubit_idx_set)
        x_index = total_qubit_idx_set.pop()
        x = generate_qubit_symbols(x_index)
        
        sub_nf = [0]*2
        ws = [0]*2
        
        for i in range(2):
            sub_nf_f = f.xreplace({x[i]: 1, x[1-i]: 0}, x_index)
            sub_nf_g = g.xreplace({x[i]: 1, x[1-i]: 0}, x_index)
            sub_nf[i] = opfunc(sub_nf_f, sub_nf_g)
            ws[i] = sub_nf[i].weight
            sub_nf[i].weight = sp.Integer(1) if not sub_nf[i].pbf.expr.is_zero else sp.Integer(0)            
        
        new_qubit_idx_set = sub_nf[0].qubit_idx_set.union(sub_nf[1].qubit_idx_set)
        new_qubit_idx_set.add(x_index)
        
        return NormalForm.shannon_expansion(sub_nf, ws, x, new_qubit_idx_set, x_index)

    @staticmethod
    def create_from_pbf(pbf):
        NormalForm.normal_form_init(pbf)
        
    @staticmethod
    def normal_form_init(pbf, qubit_idx_set=set(),tolerance=6):

        #evalf
        tolerance1=10**-tolerance
        pbf=PBF(sp.nsimplify(pbf.expr,tolerance=tolerance1,rational=False))

        if pbf.is_constant():
            return NormalForm(pbf)
        
        if len(qubit_idx_set) == 0:
            qubit_idx_set = pbf.get_qubit_idx_set()
            
        sub_qubit_idx_set = qubit_idx_set.copy()
        x_index = sub_qubit_idx_set.pop()
        x = generate_qubit_symbols(x_index)
        
        w = [0]*2
        f = [0]*2
        sub_pbf = [0]*2
        
        for i in range(2):
            sub_pbf[i] = pbf.xreplace({x[i]: 1, x[1-i]: 0})
            if sub_pbf[i].is_constant() or len(qubit_idx_set) == 0:
                w[i] = sub_pbf[i].expr
            else:
                w[i] = max(sub_pbf[i].get_coeffs_poly(), key=abs)
            
        def sub_init(weight,f, idx_set):
            new_pbf = PBF( f.expr/weight if weight!=complex(0) else complex(0) )
            return NormalForm.normal_form_init(new_pbf, idx_set)


        for i in range(2):
            f[i]=sub_init(w[i],sub_pbf[i], sub_qubit_idx_set)


        return NormalForm.shannon_expansion(f, w, x, qubit_idx_set, x_index)
    @staticmethod
    def shannon_expansion(fs, ws, xs, qubit_idx_set, x_index):

        if complex(ws[1]) == complex(0):
            result =  (ws[0]*fs[0].weight, xs[0]*fs[0].pbf.expr)
        elif fs[0].pbf == fs[1].pbf:
            if complex(ws[0]) == complex(ws[1]):
                result =  (ws[0]*fs[0].weight, fs[0].pbf.expr)
                #qubit_idx_set.pop()
                qubit_idx_set.remove(x_index)
            else:
                result =  (ws[1]*fs[0].weight, (xs[1]+ws[0]/ws[1]*xs[0])*fs[0].pbf.expr)
        else:
            result = (ws[1]*fs[1].weight, (xs[1]*fs[1].pbf.expr+ws[0]/ws[1]*fs[0].weight/fs[1].weight*xs[0]*fs[0].pbf.expr))

        return NormalForm(PBF(result[1]), weight=result[0], qubit_idx_set=qubit_idx_set)

    @property
    def find_fork_idx(self):
        qubit_idx_set=self.qubit_idx_set.copy()
        if len(qubit_idx_set)==0 or len(qubit_idx_set)==1:
            return sp.oo
        remaining_qubits_num=len(qubit_idx_set)
        expr=self.pbf.expr 
        first_alphabet=sp.srepr(expr)[0] #get first_alphabet 
        # print(sp.srepr(expr))
        if first_alphabet=='A': #fork
            fork_idx=qubit_idx_set.pop()
        elif first_alphabet=='M':
            f=len(expr.as_coeff_mul()[1])
            # print(f)
            if f==remaining_qubits_num: # all no fork
                fork_idx=sp.oo
            else: # fork at len
                fork_idx=list(qubit_idx_set)[f-1]
        return fork_idx

    def mul_coeff_dict(self,fork_idx=sp.oo):
        f=self.pbf.expr.as_coeff_mul()[1]
        idx_set=set(filter(lambda i:i<fork_idx,self.qubit_idx_set))
        coeff=dict()
        count=np.arange(len(idx_set))
        count2=len(idx_set)
        def in_idx_set(g,qubit_idx):
            if len(g)==1:
                output=(0,1) if g[0].name[1]=='n' else (1,0)
            else:
                # print('g=',g[1])
                output=(1,g[1].as_poly().coeffs()[0])
            coeff[qubit_idx]=output
        def out_idx_set(count2):
            g=f[count2].as_coeff_add()[1]
            count2+=1
            qubit_idx=int(g[0].name.replace('x','').replace('n',''))
            if qubit_idx in idx_set:
                in_idx_set(g,qubit_idx)
            else:
                out_idx_set(count2)
        for i in count:
            g=f[i].as_coeff_add()[1]
            qubit_idx=int(g[0].name.replace('x','').replace('n',''))
            if qubit_idx in idx_set:
                in_idx_set(g, qubit_idx)
            else:
                out_idx_set(count2)
        return coeff

    def remaining_term(self,coeff_dict,fork_idx):
        replace_set=set(filter(lambda i:i<fork_idx,self.qubit_idx_set))
        remaining_nf=self.copy()
        remaining_nf.weight=1
        for key in replace_set:
            x =generate_qubit_symbols(key)
            rule=[{x[0]:1,x[1]:0},{x[0]:0,x[1]:1}]
            if coeff_dict[key][0]==1:
                remaining_nf=remaining_nf.xreplace(rule[1],key)
            else:
                remaining_nf=remaining_nf.xreplace(rule[0],key)
        return remaining_nf

#############################
#
# Some preserved functions
#

def preorder(expr):
    cls = expr.class_key()[2]
    if cls in ['Add', 'Mul']:
        print(cls, end=' ')
    else:
        print(expr, end=' ')
    for arg in expr.args:
        preorder(arg)
        
def inorder(expr):
    op_dict = {'Add': '+', 'Mul': '*'}
    cls = expr.class_key()[2]
    if cls in ['Add', 'Mul']:
        print('(', end='')
        inorder(expr.args[0])
        for arg in expr.args[1:]:
            print(op_dict[cls], end=' ')
            inorder(arg)
        print(')', end='')
    else:
        print(expr, end=' ')

def max_value_from_dict(d):
    d_items = d.items()
    return max(d_items, key=lambda k: abs(operator.itemgetter(1)(k)) )

def max_value_from_expr(expr):
    return max(sp.Poly(expr).coeffs(), key=abs)


#
# End of some preserved functions
#
#############################