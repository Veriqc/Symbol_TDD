from matplotlib.font_manager import _Weight
import sympy

class primary_normal_form:
    def __init__(self,
                    pnf_list:list):
        self.__pnf_list: list=pnf_list
    @property
    def get_pnf_list(self) -> list:
        return self.__pnf_list

class normal_form:
    def __init__(self,
                    weight:int,
                    primary_normal_form:list):
        self.__weight=self.weight
        self.__primary_normal_form=self.primary_normal_form

    @property
    def weights(self) -> int:
        return self.__weight
    @property
    def primary_normal_form(self) -> primary_normal_form:
        return self.primary_normal_form
    
    @staticmethod
    def normal_form_init(Bool_Poly):
        '''
        The goal of this code is want to turn 'Psudo Boolean fuanction' into 'Normal_form'.
        
        input: 
            Bool_Poly: Psudo Boolean fuanction.
        output:
            res: Normal_form 
                [weight,[x_1+c_1*x1_n,...,x_n+c_n*xn_n]]
        '''
        symbol_list=list(Bool_Poly.free_symbols)
        order_dic={}
        symbol_dic={}
        using_qubit_list=[]
        res=[1,[]]
        for i , element in enumerate(symbol_list):
            order_dic[element.name]=i 
            symbol_dic[element.name]=element
            qubit_label=int(element.name.replace('x','').replace('n',''))

            if qubit_label not in using_qubit_list:
                using_qubit_list.append(qubit_label)
                # print(qubit_label)
            using_qubit_list=sorted(using_qubit_list)
        # print(order_dic)
        # print(symbol_dic)
        # print(using_qubit_list)    
        for i in using_qubit_list:
            if 'x%i'%i in order_dic.keys() and 'xn%i'%i in order_dic.keys():
                # print('x11:'+str(i))
                fx1=sympy.Poly(Bool_Poly,symbol_dic['x%i'%i]).coeffs()[0]
                fx0=sympy.Poly(Bool_Poly,symbol_dic['xn%i'%i]).coeffs()[0]
                c=sympy.simplify(fx0/fx1)
                Factor=symbol_dic['x%i'%i]+c*symbol_dic['xn%i'%i]
                Bool_Poly=sympy.simplify(Bool_Poly/Factor)
                if c!=1:
                    res[1].append(Factor)
            elif 'x%i'%i in order_dic.keys():
                # print('x1:'+str(i))
                Factor=symbol_dic['x%i'%i]
                res[1].append(Factor)
                Bool_Poly=sympy.simplify(Bool_Poly/Factor)
            elif 'xn%i'%i in order_dic.keys():
                # print('x0:'+str(i))
                Factor=symbol_dic['xn%i'%i]
                res[1].append(Factor)
                Bool_Poly=sympy.simplify(Bool_Poly/Factor)
        res[0]=Bool_Poly
        weight=res[0]
        primary_normal_form=res[1]
        # print(res)
        return normal_form(weight,primary_normal_form)

        