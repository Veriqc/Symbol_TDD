from __future__ import annotations
from threading import local
from numpy import poly
from scipy.fft import fft2
import sympy

class normal_form:
    '''
        Introduce:
            Normal_form is ....
    '''
    def __init__(self,
                weight: int|complex|float,
                primary_normal_form: list,
                symbol_information:list
                ):

        self.__weight: int =complex(weight)
        self.__primary_normal_form: list =primary_normal_form
        self.__data=[weight,primary_normal_form]
        self.__symbol_information=symbol_information
        '''
            weight:weight
            primary_normal_form:[x_0+c_0*x0_n,...,x_n+c_n*xn_n]
            data:[weight,[primary_normal_form]]
            symbol_information:[symbol_dic,using_qubit_list]
        '''

    @property
    def weight(self) -> int|complex|float:
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
        return str(self.data)

    # def __str__(self):
    #     return "__str__"
        
    @staticmethod
    def normal_form_init(Bool_Poly:sympy.core.mul.Mul) ->normal_form:
        '''
        The goal of this code is want to turn 'Psudo Boolean fuanction' into 'Normal_form'.
        
        input: 
            Bool_Poly: Psudo Boolean fuanction.
        output:
            res: Normal_form 
                [weight,primary_normal_form,symbol_information]

                primary_normal_form:[x_0+c_0*x0_n,...,x_n+c_n*xn_n]
                symbol_information:[symbol_dic,using_qubit_list]
        '''

        [symbol_dic,using_qubit_list]=normal_form.find_symbol_information(Bool_Poly)
        # print(symbol_dic.keys())

        if len(using_qubit_list)==0:
            return normal_form(weight=complex(Bool_Poly),primary_normal_form=complex(1),symbol_information=[symbol_dic,using_qubit_list])

        first_using_qubit=using_qubit_list[0]

        have_x=False
        have_xn=False
         #Initialized f0 weight
        if 'xn%i'%first_using_qubit in symbol_dic:
            have_xn=True
            xn=symbol_dic['xn%i'%first_using_qubit]

            f0=sympy.Poly(Bool_Poly,xn).coeffs()[0]
            if len(using_qubit_list)!=1:
                weight0=max(sympy.Poly(f0).coeffs(), key=abs) #the coefficient of the biggest term in f|_{x=0}
                # print('weight0a=',weight0)
            else:
                weight0=complex(f0)
                # weight0=complex(0)
                f0=complex(1)
                # print('weight0b=',weight0)
        else:
            weight0=complex(0)
            # print('weight1c=',weight0)

        #Initialized f1 weight
        if 'x%i'%first_using_qubit in symbol_dic:
            have_x=True
            x=symbol_dic['x%i'%first_using_qubit]        
            f1=sympy.Poly(Bool_Poly,x).coeffs()[0]
            if len(using_qubit_list)!=1:
                weight1=max(sympy.Poly(f1).coeffs(), key=abs) #the coefficient of the biggest term in f|_{x=1}
                # print('weight1a=',weight1)
            else:
                weight1=complex(f1)
                # weight1=complex(0)
                f1=complex(1)
                # print('weight1b=',weight1)
        else:
            weight1=complex(0)


        #Check weight0=0 or not    
        if weight0==complex(0):
            f0=complex(0)
        else:
            f0=normal_form.normal_form_init(sympy.sympify(f0/weight0))
    
        if weight1==complex(0):
            f1=complex(0)
        else:
            f1=normal_form.normal_form_init(sympy.sympify(f1/weight1))


        # if type(f0)==normal_form:
        #     print('f0_data',f0.data)
        # if type(f1)==normal_form:
        #     print('f1_data',f1.data)
        # print('f0=',f0)
        # print('f1=',f1)


        
        if normal_form.is_equal(f0,f1): #check if f0, f1 equall.
            if weight1==weight0:
                # print('weight1=',weight1)
                if type(f0)==normal_form:
                    if f0.primary_normal_form==complex(0): #set weight0=0 when primary_normal_form=0  
                        weight0=complex(0)
                    return normal_form(weight=weight0,primary_normal_form=f0.primary_normal_form,symbol_information=[symbol_dic,using_qubit_list])
                else:
                    if f0==complex(0):
                        weight0=complex(0)
                    return normal_form(weight=weight0,primary_normal_form=f0,symbol_information=[symbol_dic,using_qubit_list])
            else:
                if type(f0)==normal_form:
                    f0_primary_normal_form=f0.primary_normal_form
                else:
                    f0_primary_normal_form=f0
                if have_x:
                    if have_xn:
                        if f0_primary_normal_form!=complex(0):
                            return normal_form(weight=weight1,primary_normal_form=[[x+weight0/weight1*xn],f0_primary_normal_form],symbol_information=[symbol_dic,using_qubit_list])
                        else:
                            return normal_form(weight=weight1,primary_normal_form=[x+weight0/weight1*xn],symbol_information=[symbol_dic,using_qubit_list])
                    else:
                        if f0_primary_normal_form!=complex(0):
                            return normal_form(weight=weight1,primary_normal_form=[x,f0_primary_normal_form],symbol_information=[symbol_dic,using_qubit_list]) 
                        else:
                            return normal_form(weight=weight1,primary_normal_form=[x],symbol_information=[symbol_dic,using_qubit_list]) 
                else:
                    if f0_primary_normal_form!=complex(0):
                        return normal_form(weight=weight1,primary_normal_form=[weight0/weight1*xn,f0_primary_normal_form],symbol_information=[symbol_dic,using_qubit_list])     
                    else:
                        return normal_form(weight=weight1,primary_normal_form=[weight0/weight1*xn],symbol_information=[symbol_dic,using_qubit_list])     
        else:
            if weight1==complex(0):
                return normal_form(weight=weight0,primary_normal_form=[xn,f0.primary_normal_form],symbol_information=[symbol_dic,using_qubit_list])
            else:
                if type(f1)==normal_form and type(f0)==normal_form :
                    print(f1.primary_normal_form,f0.primary_normal_form)
                    if f1.primary_normal_form==[]:
                        if f0.primary_normal_form==[]:
                            output=[x+weight0/weight1*xn]
                        else:
                            output=[x,'+',[weight0/weight1*xn,f0.primary_normal_form]]
                    elif f0.primary_normal_form==[]:
                        output=[[x,f1.primary_normal_form],'+',xn]
                    else:
                        output=[[x,f1.primary_normal_form],'+',[weight0/weight1*xn,f0.primary_normal_form]]

                    primary_normal_form=output
                    print('e',f1.primary_normal_form)

                        
                else:
                    if weight0==complex(0):
                        if f1.primary_normal_form==complex(1):
                            primary_normal_form=[x]
                            print('c',f1.primary_normal_form)
                        else:
                            primary_normal_form=[x,f1.primary_normal_form]
                            print('a',f1.primary_normal_form)
                    else:
                        primary_normal_form=[x*f1+weight0/weight1*xn*f0]
                        print('b',f1.primary_normal_form)


                return normal_form(weight=weight1,primary_normal_form=primary_normal_form,symbol_information=[symbol_dic,using_qubit_list])
        

        # for i in using_qubit_list:
        #     if 'x%i'%i in symbol_dic.keys() and 'xn%i'%i in symbol_dic.keys():
        #         # print('x11:'+str(i))
        #         fx1=sympy.Poly(Bool_Poly,symbol_dic['x%i'%i]).coeffs()[0]
        #         fx0=sympy.Poly(Bool_Poly,symbol_dic['xn%i'%i]).coeffs()[0]
        #         c=sympy.simplify(fx0/fx1)
        #         Factor=symbol_dic['x%i'%i]+c*symbol_dic['xn%i'%i]
        #         Bool_Poly=sympy.simplify(Bool_Poly/Factor)
        #         if c!=1:
        #             res[1].append(Factor)
        #         else:
        #             del_list.append(i)
        #     elif 'x%i'%i in symbol_dic.keys():
        #         # print('x1:'+str(i))
        #         Factor=symbol_dic['x%i'%i]
        #         res[1].append(Factor)
        #         Bool_Poly=sympy.simplify(Bool_Poly/Factor)
        #     elif 'xn%i'%i in symbol_dic.keys():
        #         # print('x0:'+str(i))
        #         Factor=symbol_dic['xn%i'%i]
        #         res[1].append(Factor)
        #         Bool_Poly=sympy.simplify(Bool_Poly/Factor)
          
        # res[0]=Bool_Poly

        # normal_form.symbol_information_deleter(del_list,symbol_dic,using_qubit_list)

        # return normal_form(weight=res[0],primary_normal_form=res[1],symbol_information=[symbol_dic,using_qubit_list])

    @staticmethod
    def add(nf1,nf2)->normal_form:
        if len(nf1.data)==1 and len(nf2.data)==1:
            # print(nf1.weight+nf2.weight)
            return nf.normal_form(weight=nf1.weight+nf2.weight,primary_normal_form=[],symbol_information=[{},[]])

        if nf1.weight==0:
            # print(nf2)
            return nf2

        if nf2.weight==0:
            # print(nf1)
            return nf1
            
        [symbol_dic1,using_qubit_list1]=nf1.symbol_information
        [symbol_dic2,using_qubit_list2]=nf2.symbol_information

        first_using_qubit=sorted(list(set(using_qubit_list1+using_qubit_list2)))[0]
        # print(first_using_qubit)
        symbol_dic={**symbol_dic1,**symbol_dic2}
        # print(symbol_dic)
        x=symbol_dic['x%i'%first_using_qubit]
        xn=symbol_dic['xn%i'%first_using_qubit]

        # print('x%i'%first_using_qubit,'xn%i'%first_using_qubit)

    #Get f property
        # f1 area
        nf11=nf1.sub_normal_form(first_using_qubit,x)
        nf21=nf2.sub_normal_form(first_using_qubit,x)
        f1=normal_form.add(nf11,nf21)
        # f0 area
        nf10=nf1.sub_normal_form(first_using_qubit,xn)
        nf20=nf2.sub_normal_form(first_using_qubit,xn)
        f0=normal_form.add(nf10,nf20)

        if f1.weight==0:
            #initialize
            primary_normal_form=list(f0.primary_normal_form)
            symbol_dic=dict(f0.symbol_information[0])
            using_qubit_list=list(f0.symbol_information[1])
            #adjust
            weight=f0.weight
            primary_normal_form.insert(xn,0)
            symbol_dic['xn%i'%first_using_qubit]=xn
            using_qubit_list.insert(first_using_qubit,0)

            nf= nf.normal_form(weight=weight,primary_normal_form=primary_normal_form,symbol_information=[symbol_dic,using_qubit_list])
            print(nf.data,nf.symbol_information)  
            return nf
        if f0.primary_normal_form==f1.primary_normal_form:
            if f0.weight==f1.weight:
                return f0
            else:
                #initialize
                primary_normal_form=list(f0.primary_normal_form)
                symbol_dic=dict(f0.symbol_information[0])
                using_qubit_list=list(f0.symbol_information[1])
                #adjust
                weight=f1.weight
                primary_normal_form.insert(x+f0.weight/f1.weight*xn,0)

                symbol_dic['x%i'%first_using_qubit]=x
                symbol_dic['xn%i'%first_using_qubit]=xn
                using_qubit_list.insert(first_using_qubit,0)

                nf= nf.normal_form(weight=weight,primary_normal_form=primary_normal_form,symbol_information=[symbol_dic,using_qubit_list])
                print(nf.data,nf.symbol_information)  
                return nf
        else:
            
            #initialize

            primary_normal_form_1=list(f1.primary_normal_form)
            primary_normal_form_0=list(f0.primary_normal_form)

            #adjust
            weight=f1.weight
            # insert=
            primary_normal_form.insert(x+f0.weight/f1.weight*xn,0)

            symbol_dic['x%i'%first_using_qubit]=x
            symbol_dic['xn%i'%first_using_qubit]=xn
            using_qubit_list.insert(first_using_qubit,0)

            nf= nf.normal_form(weight=weight,primary_normal_form=primary_normal_form,symbol_information=[symbol_dic,using_qubit_list])
            print(nf.data,nf.symbol_information)  
            return nf

    @staticmethod
    def find_symbol_information(Bool_Poly:normal_form|sympy.core.mul.Mul|int|complex|float) -> list:
        '''
            Input Bool_Poly or normal_form to get the symbol_information. Inclued symbol_dic,using_qubit_list.

            input : Bool_Poly|normal_form

            output : symbol_information that be ordered [symbol_dic,using_qubit_list]
            
            symbol_dic : Type = dictionary 
                        Collect the symbol item with symbol lebel.
                        e.x. symbol_dic['xi']= symbol item
            using_qubit_list : Type = list 
                        The sorted using qubit list.  
                        e.x. xn5, x1, x2, xn2 are in used. using_qubit_list=[1,2,5]

        '''             

        if type(Bool_Poly)==normal_form:
            return Bool_Poly.symbol_information
        

        Bool_Poly=sympy.simplify(Bool_Poly)

        symbol_dic={}
        using_qubit_list=[]
        symbol_list=list(Bool_Poly.free_symbols)

        if len(symbol_list)==0:
            return [symbol_dic,using_qubit_list]

        for element in (symbol_list):
            symbol_dic[element.name]=element
            qubit_label=int(element.name.replace('x','').replace('n',''))

            if qubit_label not in using_qubit_list:
                using_qubit_list.append(qubit_label)
                # print(qubit_label)
            using_qubit_list=sorted(using_qubit_list)

        return [symbol_dic,using_qubit_list]

    @staticmethod
    def symbol_information_deleter(del_list,symbol_dic,using_qubit_list) -> list:
        '''
            Delete useless symbol_information after finded normal form.
        '''
        [symbol_dic,using_qubit_list]=[symbol_dic,using_qubit_list]

        for i in del_list:
            using_qubit_list.remove(i)
            del symbol_dic['x%i'%i]
            del symbol_dic['xn%i'%i]
        using_qubit_list=sorted(using_qubit_list)

        return [symbol_dic,using_qubit_list]


#
    #@staticmethod
    def is_equal (nf1,nf2)->bool:

        if type(nf1)==normal_form:
            nf1=nf1.primary_normal_form
        if type(nf2)==normal_form:
            nf2=nf2.primary_normal_form
        
        # print('nf0=',nf1)
        # print('nf1=',nf2)
        # print(nf1==nf2)
        if nf1==nf2:
            return True
        else:
            return False

    @staticmethod
    def sub_normal_form(Bool_Poly:normal_form|sympy.core.mul.Mul,symbol)-> normal_form: 

        if type(Bool_Poly)==normal_form:
            primary_normal_form=list(Bool_Poly.primary_normal_form)
            symbol_dic=dict(Bool_Poly.symbol_information[0])
            using_qubit_list=list(Bool_Poly.symbol_information[1])

            if symbol.name in symbol_dic.keys():
                weight=Bool_Poly.weight*sympy.Poly(primary_normal_form[0],symbol_dic[symbol.name]).coeffs()[0]
                # print(primary_normal_form[0])
                # print(sympy.Poly(primary_normal_form[0],symbol_dic[symbol]))
                # print(sympy.Poly(primary_normal_form[0],symbol_dic[symbol]).coeffs())
                # print(weight)
                del primary_normal_form[0]
                del symbol_dic[symbol.name]
                del using_qubit_list[0]

                nf= normal_form(weight=weight,primary_normal_form=primary_normal_form,symbol_information=[symbol_dic,using_qubit_list])
                print(nf.data,nf.symbol_information)  
                return nf

            else:
                nf= normal_form(weight=0,primary_normal_form=[],symbol_information=[{},[]])
                print(nf.data,nf.symbol_information)
                return nf


        # else:
        #     fx=sympy.Poly(Bool_Poly,symbol).coeffs()[0]
        #     [symbol_dic,using_qubit_list]=normal_form.find_symbol_information(fx)
        #     weight=max(sympy.Poly(fx).coeffs())
            
        #     nf=normal_form.normal_form_init(fx)

        #     return nf


    