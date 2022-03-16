from __future__ import annotations
from threading import local
from numpy import poly, ubyte
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
        return str(self.__data)

    # def __str__(self):
    #     return "__str__"
        
    @staticmethod
    def normal_form_init(Bool_Poly:sympy.core.mul.Mul|sympy.core.add.Add|int|float|complex) ->normal_form:
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
            if complex(Bool_Poly)==0:
                return normal_form(weight=complex(Bool_Poly),primary_normal_form=[complex(0)],symbol_information=[symbol_dic,using_qubit_list])
            else:
                return normal_form(weight=complex(Bool_Poly),primary_normal_form=[complex(1)],symbol_information=[symbol_dic,using_qubit_list])

        first_using_qubit=using_qubit_list[0]

        '''
        If first_using_qubit have x or xn initialized x or xn.

        '''

        have_x=False
        have_xn=False

        if 'xn%i'%first_using_qubit in symbol_dic:
            xn=symbol_dic['xn%i'%first_using_qubit]
            have_xn=True
        if 'x%i'%first_using_qubit in symbol_dic:
            have_x=True
            x=symbol_dic['x%i'%first_using_qubit]   
        
        # print(have_x,have_xn)

        '''
        Initialized f0 weight and f0(weight0!=0) 
        if x=1 xn=0, vice versa.
        If first_using_qubit have xn, do upper fuction.
        '''

        if have_xn:
            f0=Bool_Poly.xreplace({xn:1})
            # print('have_xn',sympy.Poly(Bool_Poly,xn).coeffs()[0],type(sympy.Poly(Bool_Poly,xn).coeffs()[0]))
            # f0=[sympy.Poly(Bool_Poly,xn).coeffs()[0]]
            if have_x:
                f0=f0.xreplace({x:0})
            # print(f0)
            f0=[sympy.simplify(f0)]
            [symbol_dic0,using_qubit_list0]=normal_form.find_symbol_information(f0[0])
            # print(using_qubit_list0)
            if len(using_qubit_list0)!=0:
                weight0=max(sympy.Poly(f0[0]).coeffs(), key=abs) #the coefficient of the biggest term in f|_{x=0}
                # print('weight0a=',weight0)
            else:
                # print(f0)
                # print(using_qubit_list)
                if len(using_qubit_list)==1:
                    weight0=f0[0]
                    f0=[complex(1)]
                else:
                    weight0=f0[0]
        else:
            weight0=complex(0)

        '''
        Initialized f1 weight and f1(weight1!=0)
        If first_using_qubit have x, do upper fuction.
        '''
        if have_x:     
            f1=Bool_Poly.xreplace({x:1})
            if have_xn:
                f1=f1.xreplace({xn:0})
            # print(f1,type(f1))
            f1=[sympy.simplify(f1)]
            [symbol_dic1,using_qubit_list1]=normal_form.find_symbol_information(f1[0])
            # print(using_qubit_list1)
            if len(using_qubit_list1)!=0:
                weight1=max(sympy.Poly(f1[0]).coeffs(), key=abs) #the coefficient of the biggest term in f|_{x=1}
            else:
                # print(f1)
                # print(using_qubit_list)
                if len(using_qubit_list)==1:
                    weight1=f1[0]
                    f1=[complex(1)]
                else:
                    weight1=f1[0]
        else:
            weight1=complex(0)

        print('weight0=',weight0)
        print('weight1=',weight1)

        '''
        Check weight=0 or not

        '''   
        if weight0==complex(0):
            f0=normal_form.normal_form_init(0)
        else:
            f0=normal_form.normal_form_init(sympy.sympify(f0[0]/weight0))
    
        if weight1==complex(0):
            f1=normal_form.normal_form_init(0)
        else:
            f1=normal_form.normal_form_init(sympy.sympify(f1[0]/weight1))

        # print('f0=',f0,type(f0))
        # print('f1=',f1,type(f1))

        '''
        In this stage, f0,f1 may be normal_form.

        type(weight)=complex,
        primary_normal_form=[primary_normal_form] 
            p.s. In this stage, type(primary_normal_form) also may be list.

        symbol_information=[symbol_dic,using_qubit_list] that remove first using qubit.

        Then belowing is return stage. 
        
        1. If weight1==complex(0):
        2. If f1==f0:
            If weight1==weight0:
           else: 
        3. else:

        '''
        f0_out=f0.primary_normal_form
        f1_out=f1.primary_normal_form

        if weight1==complex(0):
            # print('weight1==complex(0)',f0)
            # print(type(f0))
            if f0_out[0]==complex(0):
                return normal_form(weight=complex(0),primary_normal_form=[complex(0)],symbol_information=[symbol_dic,using_qubit_list])
            elif f0_out[0]==complex(1):
                return normal_form(weight=weight0,primary_normal_form=[xn],symbol_information=[symbol_dic,using_qubit_list])
            else:
                return normal_form(weight=weight0,primary_normal_form=[xn,*f0_out],symbol_information=[symbol_dic,using_qubit_list])

        print('f0=',f0,type(f0))
        print('f1=',f1,type(f1))

        if str(f0)==str(f1): 
            if weight1==weight0:
                return normal_form(weight=weight0,primary_normal_form=[*f0_out],symbol_information=[symbol_dic,using_qubit_list])
            else:
                if have_xn:
                    if f0_out[0]==complex(1): 
                        return normal_form(weight=weight1,primary_normal_form=[x+weight0/weight1*xn],symbol_information=[symbol_dic,using_qubit_list])
                    else:
                        return normal_form(weight=weight1,primary_normal_form=[x+weight0/weight1*xn,*f0_out],symbol_information=[symbol_dic,using_qubit_list])
                else:
                    if f0_out[0]==complex(1):
                        return normal_form(weight=weight1,primary_normal_form=[x],symbol_information=[symbol_dic,using_qubit_list])
                    else:
                        return normal_form(weight=weight1,primary_normal_form=[x,*f0_out],symbol_information=[symbol_dic,using_qubit_list])
        else:
            print('is_not_equal',f1)
            if have_xn:
                if f1_out[0]==complex(1):
                    if f0_out[0]==complex(1):
                        return normal_form(weight=weight1,primary_normal_form=[x+weight0/weight1*xn],symbol_information=[symbol_dic,using_qubit_list])
                    else:
                        return normal_form(weight=weight1,primary_normal_form=[x,'+',weight0/weight1*xn,*f0_out],symbol_information=[symbol_dic,using_qubit_list])
                else:
                    if f0_out[0]==complex(1):
                        return normal_form(weight=weight1,primary_normal_form=[x,*f1_out,'+',weight0/weight1*xn],symbol_information=[symbol_dic,using_qubit_list])  
                    else:
                        return normal_form(weight=weight1,primary_normal_form=[x,*f1_out,'+',weight0/weight1*xn,*f0_out],symbol_information=[symbol_dic,using_qubit_list])  
            elif f1_out[0]==complex(1):
                return normal_form(weight=weight1,primary_normal_form=[x],symbol_information=[symbol_dic,using_qubit_list])
            else:
                # print('a',f1)
                return normal_form(weight=weight1,primary_normal_form=[x,*f1_out],symbol_information=[symbol_dic,using_qubit_list])

    @staticmethod
    def normal_form_init_2(Bool_Poly:sympy.core.mul.Mul|sympy.core.add.Add|int|float|complex) ->normal_form:
        '''
        The goal of this code is want to turn 'Psudo Boolean fuanction' into 'Normal_form'.
        
        input: 
            Bool_Poly: Psudo Boolean fuanction.
        output:
            res: Normal_form 
                [weight,primary_normal_form,symbol_information]

                primary_normal_form:[x_0+c_0,...,x_n+c_n]
                symbol_information:[symbol_dic,using_qubit_list]
        '''

        [symbol_dic,using_qubit_list]=normal_form.find_symbol_information(Bool_Poly)
        # print(symbol_dic.keys())

        if len(using_qubit_list)==0:
            return normal_form(weight=complex(Bool_Poly),primary_normal_form=[complex(1)],symbol_information=[symbol_dic,using_qubit_list])

        first_using_qubit=using_qubit_list[0]

        '''
        Setup first_using_qubit (i) symble xi
        and find the biggest term weight_i in f.
        '''
        x=symbol_dic['x%i'%first_using_qubit]   
        # weight_i=max(sympy.Poly(Bool_Poly).coeffs(), key=abs)

        weight_i=sympy.Poly(Bool_Poly).coeffs()[0]

        '''
        Initialized f1 weight and f1(weight1!=0)
        If first_using_qubit have x, do upper fuction.
        '''
     
        f1=Bool_Poly.xreplace({x:1})
        f1=[sympy.simplify(f1)]
        [symbol_dic1,using_qubit_list1]=normal_form.find_symbol_information(f1[0])
        # print(using_qubit_list1)
        if len(using_qubit_list1)!=0:
            # weight1=max(sympy.Poly(f1[0]).coeffs(), key=abs) #the coefficient of the biggest term in f|_{x=1}
            weight1=sympy.Poly(f1[0]).coeffs()[0]
        else:
            # print(f1)
            # print(using_qubit_list)
            if len(using_qubit_list)==1:
                weight1=f1[0]
                f1=[complex(1)]
            else:
                weight1=f1[0]
        '''
        Initialized f0 weight and f0(weight1!=0)
        '''
        f0=Bool_Poly.xreplace({x:0})
        f0=[sympy.simplify(f0)]
        [symbol_dic1,using_qubit_list0]=normal_form.find_symbol_information(f0[0])
        # print(using_qubit_list1)
        if len(using_qubit_list0)!=0:
            # weight0=max(sympy.Poly(f0[0]).coeffs(), key=abs) #the coefficient of the biggest term in f|_{x=0}
            weight0=sympy.Poly(f0[0]).coeffs()[0]

        else:
            # print(f1)
            # print(using_qubit_list)
            if len(using_qubit_list)==1:
                weight0=f0[0]
                f0=[complex(1)]
            else:
                weight0=f0[0]

        weight1/=weight_i
        weight0/=weight_i
        print('weight_in=',weight_i)
        print('weight1=',weight1)
        print('weight0=',weight0)


        '''
        Check weight=0 or not
        '''   
        if weight0==complex(0):
            f0=normal_form.normal_form_init(0)
        else:
            f0=normal_form.normal_form_init(sympy.sympify(f0[0]/weight0))
    
        if weight1==complex(0):
            f1=normal_form.normal_form_init(0)
        else:
            f1=normal_form.normal_form_init(sympy.sympify(f1[0]/weight1))

        # print('f0=',f0,type(f0))
        # print('f1=',f1,type(f1))

        '''
        In this stage, f0,f1 may be normal_form.

        type(weight)=complex,
        primary_normal_form=[primary_normal_form] 
            p.s. In this stage, type(primary_normal_form) also may be list.

        symbol_information=[symbol_dic,using_qubit_list] that remove first using qubit.

        Then belowing is return stage. 
        
        1. If weight1==complex(0):
        2. If f1==f0:
            If weight1==weight0:
           else: 
        3. else:

        '''

        

        f1_out=f1.primary_normal_form        
        f0_out=f0.primary_normal_form


        if weight1==complex(0):
            # print('weight1==complex(0)',f0)
            # print(type(f0))
            if f0_out[0]==complex(0):
                return normal_form(weight=complex(0),primary_normal_form=[complex(0)],symbol_information=[symbol_dic,using_qubit_list])
            elif f0_out[0]==complex(1):
                return normal_form(weight=-weight_i*weight0,primary_normal_form=[x-1],symbol_information=[symbol_dic,using_qubit_list])
            else:
                return normal_form(weight=-weight_i*weight0,primary_normal_form=[x-1,*f0_out],symbol_information=[symbol_dic,using_qubit_list])

        print('f0=',f0,type(f0))
        print('f1=',f1,type(f1))

        if str(f0)==str(f1): 
            if weight1==weight0:
                return normal_form(weight=weight_i*weight1,primary_normal_form=[*f1_out],symbol_information=[symbol_dic,using_qubit_list])
            else:
                weight_out=weight_i*weight1**2/(weight1-weight0)
                const=weight0/(weight1-weight0)
                if f1_out[0]==complex(1):
                    return normal_form(weight=weight_out,primary_normal_form=[x+const],symbol_information=[symbol_dic,using_qubit_list])
                else:
                    return normal_form(weight=weight_out,primary_normal_form=[x+const,*f1_out],symbol_information=[symbol_dic,using_qubit_list])
        else:
            print('is_not_equal f1=',f1)
            print('is_not_equal f0=',f0)
            if weight0==complex(0):
                if f1_out[0]==complex(1):
                    return normal_form(weight=weight_i*weight1,primary_normal_form=[x],symbol_information=[symbol_dic,using_qubit_list])
                else:
                    return normal_form(weight=weight_i*weight1,primary_normal_form=[x,*f1_out],symbol_information=[symbol_dic,using_qubit_list])
            else:
                if f1_out[0]==complex(1):
                    if f0_out[0]==complex(1):
                        return normal_form(weight=weight_i*(weight1-weight0),primary_normal_form=[x+weight0/(weight1-weight0)],symbol_information=[symbol_dic,using_qubit_list])
                    else:
                        return normal_form(weight=weight_i*weight1,primary_normal_form=[x,'-',weight0/weight1,x-1],symbol_information=[symbol_dic,using_qubit_list])
                else:
                    return normal_form(weight=weight_i*weight1,primary_normal_form=[x,*f1_out,'-',weight0/weight1,x-1,*f0_out],symbol_information=[symbol_dic,using_qubit_list])

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


    