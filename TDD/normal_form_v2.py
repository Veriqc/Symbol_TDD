from unittest import FunctionTestCase
import sympy as sp

def generate_qubit_symbols(x_index):
    x = ['']*2
    for i in range(2):
        n_str = 'n' if i == 0 else ''
        x[i] = sp.symbols('x%s%i' % (n_str, x_index))
        
    return x

class PBF():
    def __init__(self, expr):
        self.expr = expr if isinstance(expr, sp.Basic) else sp.simplify(expr)
        
    def __repr__(self):
        return str(self.expr)
        
    def get_coeffs_poly(self):
        return sp.Poly(self.expr).coeffs()

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

class NormalForm():
    def __init__(self, pbf, weight=1, qubit_idx_set=set()):
        #print('==init==', pbf, weight)
        self.qubit_idx_set = qubit_idx_set.copy()
        
        if weight == 0 or pbf.expr.is_zero:
            self.pbf = PBF(0)
            self.weight = 0
            return
        
        if pbf.is_constant():
            pbf_weight = pbf.expr
        elif len(pbf.expr.args) > 0:
            pbf_weight = pbf.expr.args[0]
        else:
            pbf_weight = sp.Integer(1)
        
        if pbf_weight.is_constant() and not pbf_weight.is_zero:
            self.pbf = PBF(pbf.expr / pbf_weight)
            self.weight = weight * pbf_weight
        else:
            self.pbf = pbf
            self.weight = weight
        
    def __repr__(self):
        return '(%s) * [ %s ], set=%s' % (str(self.weight), str(self.pbf), str(self.qubit_idx_set))

    def copy(self):
        return NormalForm(self.pbf, self.weight, self.qubit_idx_set.copy())
    
    def is_constant(self):
        return self.pbf.is_constant()
    
    def xreplace(self, rule, x_index):
        new_qubit_idx_set = self.qubit_idx_set.copy()
        if x_index in new_qubit_idx_set:
            new_qubit_idx_set.remove(x_index)
            
        #print('xreplace_remove:', self.qubit_idx_set, new_qubit_idx_set)
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
        
        return NormalForm.op_subroutine(self, g, NormalForm.__add__)
    
    def __mul__(self, g):
        if self.is_constant() or g.is_constant():
            qubit_idx_set = self.qubit_idx_set.union(g.qubit_idx_set)
            return NormalForm(PBF(self.pbf.expr * g.pbf.expr), weight=self.weight*g.weight, qubit_idx_set=qubit_idx_set)
        elif self.pbf.expr.is_zero or g.pbf.expr.is_zero:
            return NormalForm(PBF(0))

        return NormalForm.op_subroutine(self, g, NormalForm.__mul__)
    
    def __truediv__(self, g):
        if g.pbf.expr.is_zero:
            return self.copy()
        elif g.is_constant():
            new_f = self.copy()
            new_f.weight = self.weight / g.weight
            return new_f

        return NormalForm.op_subroutine(self, g, NormalForm.__truediv__)
    
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
            # seperate sub_nf[i] to (weight, primary(weight=1) normal form)
            ws[i] = sub_nf[i].weight
            sub_nf[i].weight = sp.Integer(1) if not sub_nf[i].pbf.expr.is_zero else sp.Integer(0)            
        
        new_qubit_idx_set = sub_nf[0].qubit_idx_set.union(sub_nf[1].qubit_idx_set)
        new_qubit_idx_set.add(x_index)
        
        return NormalForm.shannon_expansion(sub_nf, ws, x, new_qubit_idx_set, x_index)

    @staticmethod
    def create_from_pbf(pbf):
        NormalForm.normal_form_init(pbf)
        
    @staticmethod
    def normal_form_init(pbf, qubit_idx_set=set()):
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
            #new_pbf = PBF( sp.simplify(f.expr/weight) if weight!=complex(0) else complex(0) )
            new_pbf = PBF( f.expr/weight if weight!=complex(0) else complex(0) )
            return NormalForm.normal_form_init(new_pbf, idx_set)
        
        for i in range(2):
            f[i]=sub_init(w[i],sub_pbf[i], sub_qubit_idx_set)
        
        return NormalForm.shannon_expansion(f, w, x, qubit_idx_set, x_index)
    
    @staticmethod
    def shannon_expansion(fs, ws, xs, qubit_idx_set, x_index):
        if complex(ws[1]) == complex(0):
            result =  (ws[0]*fs[0].weight, xs[0]*fs[0].pbf.expr)
        elif str(fs[0].pbf) == str(fs[1].pbf):
            if complex(ws[0]) == complex(ws[1]):
                result =  (ws[0]*fs[0].weight, fs[0].pbf.expr)
                #qubit_idx_set.pop()
                qubit_idx_set.remove(x_index)
            else:
                result =  (ws[1]*fs[0].weight, (xs[1]+ws[0]/ws[1]*xs[0])*fs[0].pbf.expr)
        else:
            result = (ws[1]*fs[1].weight, (xs[1]*fs[1].pbf.expr+ws[0]/ws[1]*xs[0]*fs[0].weight/fs[1].weight*fs[0].pbf.expr))
            
        #print(result)
            
        return NormalForm(PBF(result[1]), weight=result[0], qubit_idx_set=qubit_idx_set)

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