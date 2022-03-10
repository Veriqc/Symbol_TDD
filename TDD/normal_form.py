from __future__ import annotations
import sympy

class normal_form:
    '''
        Introduce:
            Normal_form is ....
    '''
    def __init__(self,
                weight: int,
                primary_normal_form: list,
                symbol_information:list
                ):
  
        self.__weight: int =weight
        self.__primary_normal_form: list =primary_normal_form
        self.__data=[self.weight,self.primary_normal_form]
        self.__symbol_information=symbol_information
        '''
            weight:
            primary_normal_form:
            data:
            symbol_information:

        '''

    @property
    def weight(self) -> int:
        return self.__weight

    @property
    def primary_normal_form(self) -> list:
        return self.__primary_normal_form

    @property
    def data(self) -> list:
        return self.__data
    @property
    def symbol_information(self) -> list:
        return self.__symbol_information 


    def __repr__(self):
        return "__repr__"

    def __str__(self):
        return "__str__"
        
    @staticmethod
    def normal_form_init(Bool_Poly:sympy.core.mul.Mul) ->normal_form:
        '''
        The goal of this code is want to turn 'Psudo Boolean fuanction' into 'Normal_form'.
        
        input: 
            Bool_Poly: Psudo Boolean fuanction.
        output:
            res: Normal_form 
                [weight,[x_0+c_0*x0_n,...,x_n+c_n*xn_n]]
        '''

        [symbol_list,order_dic,symbol_dic,using_qubit_list]=normal_form.find_symbol_information(Bool_Poly)

        res=[1,[]]
        del_list=[]
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
                else:
                    del_list.append(i)
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

        normal_form.symbol_information_deleter(del_list,symbol_list,order_dic,symbol_dic,using_qubit_list)

        return normal_form(weight=res[0],primary_normal_form=res[1],symbol_information=[symbol_list,order_dic,symbol_dic,using_qubit_list])

    @staticmethod
    def add (nf1, nf2):
        b=nf1+nf2
        res=[1,[]]
        res[0]=1
        res[1]=[]
        return normal_form(weight=res[0],primary_normal_form=res[1])

    @staticmethod
    def find_symbol_information(Bool_Poly:normal_form|sympy.core.mul.Mul) -> list:
        '''
            Input Bool_Poly or normal_form to get the symbol_information. Inclued symbol_list,order_dic,symbol_dic,using_qubit_list.

            input : Bool_Poly|normal_form

            output : symbol_information that be ordered [symbol_list,order_dic,symbol_dic,using_qubit_list]
            
            symbol_list : Type = list 
                        List of the symbol item.  
            order_dic : Type = dictionary 
                        The order of each each symbol item in symbol_list.
                        e.x. order_dic['xi']= xi order in symbol_list (Type = int)
            symbol_dic : Type = dictionary 
                        The symbol item.
                        e.x. order_dic['xi']= symbol item
            using_qubit_list : Type = list 
                        The sorted using qubit list.  
                        e.x. xn5, x1, x2, xn2 are in used. using_qubit_list=[1,2,5]

        '''             
        if type(Bool_Poly)==normal_form:
            return Bool_Poly.symbol_information

        symbol_list=list(Bool_Poly.free_symbols)

        order_dic={}
        symbol_dic={}
        using_qubit_list=[]

        for i , element in enumerate(symbol_list):
            order_dic[element.name]=i 
            symbol_dic[element.name]=element
            qubit_label=int(element.name.replace('x','').replace('n',''))

            if qubit_label not in using_qubit_list:
                using_qubit_list.append(qubit_label)
                # print(qubit_label)
            using_qubit_list=sorted(using_qubit_list)

        return [symbol_list,order_dic,symbol_dic,using_qubit_list]

    @staticmethod
    def symbol_information_deleter(del_list,symbol_list,order_dic,symbol_dic,using_qubit_list) -> list:
        '''
            Delete useless symbol_information after finded normal form.
        '''
        [symbol_list,order_dic,symbol_dic,using_qubit_list]=[symbol_list,order_dic,symbol_dic,using_qubit_list]

        for i in del_list:
            using_qubit_list.remove(i)
            symbol_list.remove(symbol_dic['x%i'%i])
            symbol_list.remove(symbol_dic['xn%i'%i])
            del order_dic['x%i'%i]
            del order_dic['xn%i'%i]
            del symbol_dic['x%i'%i]
            del symbol_dic['xn%i'%i]
        for i , element in enumerate(symbol_list):
            order_dic[element.name]=i 
            symbol_dic[element.name]=element
            qubit_label=int(element.name.replace('x','').replace('n',''))

            if qubit_label not in using_qubit_list:
                using_qubit_list.append(qubit_label)
                # print(qubit_label)
            using_qubit_list=sorted(using_qubit_list)

        return [symbol_list,order_dic,symbol_dic,using_qubit_list]
    