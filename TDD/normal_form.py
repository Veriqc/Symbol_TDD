from __future__ import annotations
from typing_extensions import Self
import weakref
import sympy



class normal_form:
    '''
        Introduce:
            Normal_form is ....
    '''
    def __init__(self,
                weight: int|complex|float,
                primary_normal_form: primary_normal_form,
                using_qubit_list:list
                ):

        self.__weight: complex =complex(weight)
        self.__primary_normal_form: primary_normal_form =primary_normal_form
        self.__data=[weight,primary_normal_form]
        self.__using_qubit_list=using_qubit_list
        '''
            weight:weight
            primary_normal_form:[xf1 + w0/w1 ·x′f0]
            data:[weight,[primary_normal_form]]
            using_qubit_list:The order of using_qubit_list.
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
    def using_qubit_list(self) -> list:
        return self.__using_qubit_list


    def __repr__(self):
        return str(self.__data)

    def __str__(self):
        return str(self.__data)



    @staticmethod
    def normal_form_init(Bool_Poly:sympy.core.mul.Mul|sympy.core.add.Add|int|float|complex,tolerance=6) ->normal_form:
        '''
        The goal of this code is want to turn 'Psudo Boolean fuanction' into 'Normal_form'.
        
        input: 
            Bool_Poly: Psudo Boolean fuanction.
        output:
            res: Normal_form 
                [weight,primary_normal_form,using_qubit_list]

                primary_normal_form:[x_0+c_0*x0_n,...,x_n+c_n*xn_n]
                using_qubit_list:The order of using_qubit_list.
        '''
        
        
        #evalf
        tolerance1=10**-tolerance
        # Bool_Poly=sympy.simplify(Bool_Poly)
        Bool_Poly=sympy.nsimplify(Bool_Poly,tolerance=tolerance1,rational=False)
        # .evalf(n=tolerance)

        using_qubit_list=normal_form.find_using_qubit_list(Bool_Poly)
        # print(symbol_dic.keys())
        

        if len(using_qubit_list)==0:
            if complex(Bool_Poly)==0:
                primary_normal_form_out=primary_normal_form(weight=complex(1),primary_normal_form=[complex(0)],using_qubit_list=using_qubit_list)
                nf=normal_form(weight=complex(Bool_Poly),primary_normal_form=primary_normal_form_out,using_qubit_list=using_qubit_list)
                print('nf=',nf)
                return nf
            else:
                primary_normal_form_out=primary_normal_form(weight=complex(1),primary_normal_form=[complex(1)],using_qubit_list=using_qubit_list)
                nf=normal_form(weight=complex(Bool_Poly),primary_normal_form=primary_normal_form_out,using_qubit_list=using_qubit_list)
                print('nf=',nf)
                return nf
        

        '''
        Initialized x and xn at first_using_qubit  .

        '''
        first_using_qubit=using_qubit_list[0]
        
        xn=sympy.symbols('xn%i'%first_using_qubit)
        x=sympy.symbols('x%i'%first_using_qubit) 
        f0=sympy.simplify(Bool_Poly)
        f1=sympy.simplify(Bool_Poly)

        '''
        Initialized f0 weight and f0(weight0!=0) 
        if x=1 xn=0, vice versa.
        If first_using_qubit have xn, do upper fuction.
        '''
        f0=f0.xreplace({x:0,xn:1})
        f0=[sympy.simplify(f0)]
        using_qubit_list0=normal_form.find_using_qubit_list(f0[0])
        # print(using_qubit_list0)
        if len(using_qubit_list0)!=0:
            weight0=max(sympy.Poly(f0[0]).coeffs(), key=abs) #the coefficient of the biggest term in f|_{x=0}
        else:
            weight0=f0[0]

        '''
        Initialized f1 weight and f1(weight1!=0)
        If first_using_qubit have x, do upper fuction.
        '''
        f1=f1.xreplace({x:1,xn:0})
        f1=[sympy.simplify(f1)]
        using_qubit_list1=normal_form.find_using_qubit_list(f1[0])
        # print(using_qubit_list1)
        if len(using_qubit_list1)!=0:
            weight1=max(sympy.Poly(f1[0]).coeffs(), key=abs) #the coefficient of the biggest term in f|_{x=1}
        else:
            weight1=f1[0]

        print('weight0=',weight0)
        print('weight1=',weight1)

        '''
        Check weight=0 or not

        '''   
        if weight0==complex(0):
            f0=normal_form.normal_form_init(complex(0))
        else:
            f0=normal_form.normal_form_init(sympy.sympify(f0[0]/weight0))
    
        if weight1==complex(0):
            f1=normal_form.normal_form_init(complex(0))
        else:
            f1=normal_form.normal_form_init(sympy.sympify(f1[0]/weight1))

        # print('f0=',f0,type(f0))
        # print('f1=',f1,type(f1))

        '''
        In this stage, f0,f1 may be normal_form.

        type(weight)=complex,
        primary_normal_form=[primary_normal_form] 
            p.s. In this stage, type(primary_normal_form) also may be list.

        using_qubit_list=using_qubit_list that remove first using qubit.

        Then belowing is return stage. 
        
        1. If weight1==complex(0):
        2. If f1==f0:
            If weight1==weight0:
           else: 
        3. else:
        '''
        
        return normal_form.Shannon_expansion(f1,f0,weight1,weight0,x,xn,using_qubit_list)
        
    @staticmethod
    def add(nf1:normal_form|primary_normal_form,nf2:normal_form|primary_normal_form)->normal_form:
        
        # '''
        # Check if nf1 or nf2 is normal_form.
        # If not, do normal_form_init.
        # '''
        # print(type(nf1),type(nf2))
        # if type(nf1)!=primary_normal_form or type(nf1)!=normal_form:
        #     nf1=normal_form.normal_form_init(nf1)
        # if type(nf2)!=primary_normal_form or type(nf2)!=normal_form:
        #     nf2=normal_form.normal_form_init(nf2)
        '''
        Check if nf1 or nf2 is constant.
        '''
        if len(nf1.using_qubit_list)==0 and len(nf2.using_qubit_list)==0:
            primary_normal_form_out=primary_normal_form(weight=complex(1),primary_normal_form=[complex(1)],using_qubit_list=[])
            nf=normal_form(weight=nf1.weight+nf2.weight,primary_normal_form=primary_normal_form_out,using_qubit_list=[])
            print('nf=',nf)
            return nf
        '''
        Check nf1 or nf2 is 0.
        '''
        if nf1.weight==complex(0):
            # print(nf2)
            return nf2

        if nf2.weight==complex(0):
            # print(nf1)
            return nf1
        
        ''''
        Find first_using_qubit.
        '''
        using_qubit_list1=nf1.using_qubit_list
        using_qubit_list2=nf2.using_qubit_list

        using_qubit_list=sorted(list(set(using_qubit_list1+using_qubit_list2)))
        first_using_qubit=using_qubit_list[0]

        xn=sympy.symbols('xn%i'%first_using_qubit)
        x=sympy.symbols('x%i'%first_using_qubit) 

        '''
        Do nf_xreplace, and add
        '''
        nf1_1=nf1.nf_xreplace({x:1,xn:0})
        nf2_1=nf2.nf_xreplace({x:1,xn:0})
        nf1_0=nf1.nf_xreplace({x:0,xn:1})
        nf2_0=nf2.nf_xreplace({x:0,xn:1})
        
        f1=normal_form.add(nf1_1,nf2_1)
        f0=normal_form.add(nf1_0,nf2_0)
        
        '''
        Set Shannon_expansion input.
        '''

        weight1=f1.weight
        weight0=f0.weight
        f1=f1.primary_normal_form
        f0=f0.primary_normal_form
        using_qubit_list=using_qubit_list.remove(first_using_qubit)

        nf=normal_form.Shannon_expansion(f1,f0,weight1,weight0,x,xn,using_qubit_list)
        return nf

    @staticmethod
    def find_using_qubit_list(Bool_Poly:normal_form|sympy.core.mul.Mul|int|complex|float) -> list:
        '''
            Input Bool_Poly or normal_form to get the using_qubit_list.

            input : Bool_Poly|normal_form

            output : using_qubit_list that be sorted. 
            
            using_qubit_list : Type = list 
                        The sorted using qubit list.  
                        e.x. xn5, x1, x2, xn2 are in used. using_qubit_list=[1,2,5]

        '''             
        if type(Bool_Poly)==normal_form:
            return Bool_Poly.using_qubit_list
        

        if type(Bool_Poly)==list: #primary normal form list
            symbol_list=[]
            for item in Bool_Poly:
                symbol_list+=[*list(item.free_symbols)]
            # print(symbol_list)
        else:
            Bool_Poly=sympy.simplify(Bool_Poly)
            symbol_list=list(Bool_Poly.free_symbols)

        using_qubit_list=[]
        if len(symbol_list)==0:
            return using_qubit_list
        
        for element in (symbol_list):
            qubit_label=int(element.name.replace('x','').replace('n',''))
            if qubit_label not in using_qubit_list:
                using_qubit_list.append(qubit_label)
                # print(qubit_label)
            using_qubit_list=sorted(using_qubit_list)

        return using_qubit_list
      
    @staticmethod
    def Shannon_expansion(f1,f0,weight1,weight0,x,xn,using_qubit_list)->normal_form:
        f1_out=f1.get_primary_normal_form_list()
        f0_out=f0.get_primary_normal_form_list()
        weight1_out=f1.weight
        weight0_out=f0.weight

        if weight1==complex(0): #w0*(xn*f0)
            if f0_out[0]==complex(0): #0
                weight=complex(0)
                primary_normal_form_out=primary_normal_form(weight=complex(0),primary_normal_form=[complex(0)],using_qubit_list=[])
            elif f0_out[0]==complex(1): #w0*w0_out(xn)
                weight=weight0*weight0_out
                primary_normal_form_out=primary_normal_form(primary_normal_form=[xn],using_qubit_list=using_qubit_list)
            else: #w0*w0_out(xn*f0_out)
                weight=weight0*weight0_out
                primary_normal_form_out=primary_normal_form(primary_normal_form=[xn,*f0_out],using_qubit_list=using_qubit_list)
            nf=normal_form(weight=weight,primary_normal_form=primary_normal_form_out,using_qubit_list=using_qubit_list)
            print('nf=',nf)
            return nf
        if str(f0)==str(f1): #w1(x+w0/w1*xn)f0
            if weight1==weight0: #w0*w0_out(f0_out)
                weight=weight0*weight0_out
                del using_qubit_list[0]
                primary_normal_form_out=primary_normal_form(primary_normal_form=[*f0_out],using_qubit_list=using_qubit_list)
            else: # w1*w0_out (x+w0/w1*xn)*f0_out
                weight=weight1*weight0_out
                if f0_out[0]==complex(1): # w1*w0_out (x+w0/w1*xn)
                    primary_normal_form_out=primary_normal_form(primary_normal_form=[x+weight0/weight1*xn],using_qubit_list=using_qubit_list) 
                else: # w1*w0_out (x+w0/w1*xn)*f0_out
                    primary_normal_form_out=primary_normal_form(primary_normal_form=[x+weight0/weight1*xn,*f0_out],using_qubit_list=using_qubit_list)
            nf=normal_form(weight=weight,primary_normal_form=primary_normal_form_out,using_qubit_list=using_qubit_list)
            print('nf=',nf)
            return nf

        else: #w1(x*f1+w0/w1*xn*f0)=w1*w1_out(x*f1_out+w0*weight0_out/w1/*w1_out*xn*f0_out)  
            weight=weight1*weight1_out #set weight
            print('is_not_equal','f1=',f1,'f0=',f0)
            if f1_out[0]==complex(1): #w1*w1_out(x+w0*weight0_out/w1/w1_out*(xn*f0_out))
                if f0_out[0]==complex(1): #w1*w1_out(x+w0*weight0_out/w1/w1_out*(xn))
                    primary_normal_form_out=primary_normal_form(primary_normal_form=[x+weight0*weight0_out/weight1/weight1_out*xn],using_qubit_list=using_qubit_list)
                elif f0_out[0]==complex(0):#w1*w1_out(x)
                    primary_normal_form_out=primary_normal_form(primary_normal_form=[x],using_qubit_list=using_qubit_list)
                else: #w1*w1_out(x+w0*weight0_out/w1/w1_out*(xn*f0_out))
                    if weight0*weight0_out/weight1/weight1_out==complex(1): #w1*w1_out(x+(xn*f0_out))
                        primary_normal_form_out2=primary_normal_form(primary_normal_form=[xn,*f0_out],using_qubit_list=using_qubit_list)
                        primary_normal_form_out=primary_normal_form(primary_normal_form=[[x],'+',primary_normal_form_out2],using_qubit_list=using_qubit_list)
                    else: #w1*w1_out(x+w0*weight0_out/w1/w1_out*(xn*f0_out))
                        primary_normal_form_out2=primary_normal_form(primary_normal_form=[xn,*f0_out],using_qubit_list=using_qubit_list)
                        nf_out=normal_form(weight=weight0*weight0_out/weight1/weight1_out,primary_normal_form=primary_normal_form_out2,using_qubit_list=using_qubit_list)
                        primary_normal_form_out=primary_normal_form(primary_normal_form=[[x],'+',nf_out],using_qubit_list=using_qubit_list)
            else: #w1(x*f1+w0/w1*xn*f0)=w1*w1_out(x*f1_out+w0*weight0_out/w1/*w1_out*xn*f0_out)  
                if f0_out[0]==complex(1): #w1*w1_out(x*f1_out+w0*weight0_out/w1/*w1_out*xn)      
                    if weight0*weight0_out/weight1/weight1_out ==complex(1): #w1*w1_out(x*f1_out+xn)
                        primary_normal_form_out=primary_normal_form(primary_normal_form=[[x,*f1_out],'+',[xn]],using_qubit_list=using_qubit_list)
                    else: #w1*w1_out(x*f1_out+w0*weight0_out/w1/*w1_out*xn) 
                        primary_normal_form_out2=primary_normal_form(primary_normal_form=[xn],using_qubit_list=using_qubit_list[0])
                        nf_out=normal_form(weight=weight0*weight0_out/weight1/weight1_out,primary_normal_form=primary_normal_form_out2,using_qubit_list=using_qubit_list[0])
                        primary_normal_form_out=primary_normal_form(primary_normal_form=[[x,*f1_out],'+',nf_out],using_qubit_list=using_qubit_list)
                elif f0_out[0]==complex(0): # w1*w1_out(x*f1_out)
                    primary_normal_form_out=primary_normal_form(primary_normal_form=[x,*f1_out],using_qubit_list=using_qubit_list) 
                else: #w1*w1_out(x*f1_out+w0*weight0_out/w1/*w1_out*xn*f0_out)  
                    if weight0*weight0_out/weight1/weight1_out ==complex(1): #w1*w1_out(x*f1_out+xn*f0_out)  
                        primary_normal_form_out=primary_normal_form(primary_normal_form=[[x,*f1_out],'+',[xn,*f0_out]],using_qubit_list=using_qubit_list)
                    else: #w1*w1_out(x*f1_out+w0*weight0_out/w1/*w1_out*xn*f0_out)  
                        primary_normal_form_out2=primary_normal_form(primary_normal_form=[xn,*f0_out],using_qubit_list=using_qubit_list)
                        nf_out=normal_form(weight=weight0*weight0_out/weight1/weight1_out,primary_normal_form=primary_normal_form_out2,using_qubit_list=using_qubit_list)
                        primary_normal_form_out=primary_normal_form(primary_normal_form=[[x,*f1_out],'+',nf_out],using_qubit_list=using_qubit_list)
            nf=normal_form(weight=weight,primary_normal_form=primary_normal_form_out,using_qubit_list=using_qubit_list)
            print('nf=',nf)
            return nf  

    def get_primary_normal_form_list(cls)->list:
        pnf=cls.primary_normal_form
        if type(pnf)==list:
            # print('list')
            return pnf
        else:
            pnf=normal_form.get_primary_normal_form_list(pnf)
            return pnf
    
    @staticmethod
    def nf_xreplace(nf_in:normal_form|primary_normal_form,qubit_num:int,bool:bool)-> normal_form: 
        '''
        It will use sympy.xreplace. 
        The format is ...
        
        '''

        using_qubit_list=normal_form.find_using_qubit_list(nf_in)
        if qubit_num not in using_qubit_list:
            return nf_in 
        
        if type(nf_in)==normal_form:
            weight=nf_in.weight
            primary_normal_form_list=nf_in.get_primary_normal_form_list
        else:
            primary_normal_form_list=nf_in

        if '+' in primary_normal_form_list:
            nf1=normal_form.nf_xreplace(primary_normal_form_list[0],qubit_num=qubit_num,bool=bool)
            nf2=normal_form.nf_xreplace(primary_normal_form_list[3],qubit_num=qubit_num,bool=bool)
            weight=weight*nf1.weight
            weight2=nf2.weight/nf1.weight
            pnf1=nf1.primary_normal_form
            if weight2==complex(1):
                nf2=nf2.primary_normal_form
            elif weight2==complex(0):
                nf2=primary_normal_form(weight=complex(0),primary_normal_form=[complex(0)],using_qubit_list=[])
            else:
                nf2=normal_form(weight=weight2,primary_normal_form=nf2.primary_normal_form,using_qubit_list=nf2.using_qubit_list)
            
            if nf1.weight==0:
                primary_normal_form_out=nf2
            if nf2.weight==0:
                primary_normal_form_out=nf1
            else:    
                primary_normal_form_out=[pnf1,'+',nf2]
            using_qubit_list.remove(qubit_num)
            nf=normal_form(weight=weight,primary_normal_form=primary_normal_form_out,using_qubit_list=using_qubit_list)
        else:
            nf=normal_form.list_xreplace(pnf_in=primary_normal_form_list,qubit_num=qubit_num,bool=bool)
        
        return nf
        

    @staticmethod
    def list_xreplace(pnf_in:list,qubit_num:int,bool:bool)-> normal_form: 
        '''
        This input should be primary_normal_form.
        '''
        using_qubit_list=normal_form.find_using_qubit_list(pnf_in)
        if qubit_num not in using_qubit_list:
            nf=primary_normal_form(primary_normal_form=pnf_in,using_qubit_list=using_qubit_list)
            return nf
        else:
            index=using_qubit_list.index(qubit_num)
        
        xn=sympy.symbols('xn%i'%qubit_num)
        x=sympy.symbols('x%i'%qubit_num)
        
        if bool==True:
            rule={x:1,xn:0}
        else:
            rule={x:0,xn:1}

        weight=complex(pnf_in[index].xreplace(rule))
        del pnf_in[index]
        del using_qubit_list[index]

        if weight==complex(1):
            nf=primary_normal_form(primary_normal_form=pnf_in,using_qubit_list=using_qubit_list)
        elif weight==complex(0):
            nf=primary_normal_form(weight=weight,primary_normal_form=complex(0),using_qubit_list=[])
        else:
            nf=normal_form(weight=weight,primary_normal_form=pnf_in,using_qubit_list=using_qubit_list)
        return nf

    # @classmethod

class primary_normal_form(normal_form):
    def __init__(self,
                primary_normal_form: list,
                using_qubit_list:list,
                weight=complex(1)
                ):
        self.__weight: complex =complex(weight)
        self.__primary_normal_form: list =primary_normal_form
        self.__data=primary_normal_form
        self.__using_qubit_list=using_qubit_list
        '''
            primary_normal_form:[x_0+c_0*x0_n,...,x_n+c_n*xn_n]
            data:primary_normal_form
            using_qubit_list:The order of using_qubit_list
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
    def using_qubit_list(self) -> list:
        return self.__using_qubit_list


    def __repr__(self):
        return str(self.__data)
    
    def __str__(self):
        return str(self.__data)
