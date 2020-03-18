from sympy import symbols,Matrix,cos,sin,tan,ln,exp,log
from sympy.printing.cxxcode import cxxcode


def removeDelete(grid):
    ret = []
    loop = len(grid)
    temp = []
    for i in range (loop):
        for j in range (loop):
            if grid[i][j] != 'd':
                temp.append(grid[i][j])
        if temp != []:
            ret.append(temp)
            temp = []
    return ret

def deletevar(varlist,deleteindex):
    retlist = []
    for i in range (len(deleteindex)):
        varlist[deleteindex[i]] = 'd'
    for i in range (len(varlist)):
        if varlist[i]!= "d":
            retlist.append(varlist[i])
    return retlist

def mark(matrix,num,size):
    for i in range (size):
        for j in range (size):
            if i==num or j == num:
                matrix[i][j] = 'd'

def jacobian(difvar,diffun,loop):
    dvl = [x.strip() for x in difvar.split(',')]
    
    #print (diffun)
    new_diffun = []
    for fun in diffun:
        fun = fun.replace('pow','Pow')
        new_diffun.append(fun)
    diffun = new_diffun

    dfl = diffun

    varlen = len(dvl)
    funlen = len(dfl)

    mat_y = [[y] for y in dvl]
    mat_x = [[x] for x in dfl]

    Y = Matrix(mat_y)
    X = Matrix(mat_x)

    jac = X.jacobian(Y)
    Ljac = jac.tolist()
    print (Ljac)
   
    deleterow = []
    deletecol = []
    for i in range (varlen):
        deleteflag = 1
        for j in range (varlen):
            if Ljac[i][j] != 0:
                deleteflag = 0
        if deleteflag ==1:
            deleterow.append(i)

    for j in range (varlen):
        deleteflag = 1
        for i in range (varlen):
            if Ljac[i][j]!=0:
                deleteflag = 0
        if deleteflag ==1:
            deletecol.append(j)

    delete_element=[]
    if len(deletecol)!=0 and len(deleterow)!=0:
        delete_element = list(set(deletecol).intersection(deleterow))
        dvl = deletevar(dvl,delete_element)
        for i in range (len(delete_element)):
            mark(Ljac,delete_element[i],varlen)
        Ljac = removeDelete(Ljac)
    
    #print(Ljac)
    funlen -= len(delete_element)
    varlen -= len(delete_element)

    for i in range (funlen):
        for j in range (varlen):
            Ljac[i][j] = "Entry"+"="+cxxcode(Ljac[i][j],standard='C++11')

    naturestirng = ""
    naturestirng+=str(len(dvl))
    naturestirng+='\n'
    for i in range (len(dvl)):
        naturestirng+=dvl[i]
        naturestirng+='\n'
    #nodelist =[]
    naturestirng+=str(funlen*varlen)
    naturestirng+='\n'

    for i in range (funlen):
        for j in range (varlen):
            naturestirng+=Ljac[i][j]
            naturestirng+='\n'
    filename = ""
    filename +="../work-dir/jacobiannature"
    filename += str(loop+1)
    filename +=".txt"
    savefile=open(filename,'w')
    savefile.write(naturestirng)
    savefile.close()
    
    ComputeLDFstring = ""
    if loop == 0:
        # ComputeLDFstring += "from math import *\n"
        # ComputeLDFstring += "import numpy as np\n"
        # ComputeLDFstring += "import numpy.linalg as la\n"
        # #ComputeLDFstring += "import matplotlib.pyplot as plt\n"
        # ComputeLDFstring += "import sys\n"
        # ComputeLDFstring += "import time\n"
        # ComputeLDFstring += "from sympy import Derivative\n"
        # ComputeLDFstring += "def jcalc(listvalue,curstate):\n"
        # ComputeLDFstring += "    ret = []\n"
        ComputeLDFstring += "#include <vector>\n"
        ComputeLDFstring += "#include <cmath>\n"
        # ComputeLDFstring += "using namespace std;\n"
        ComputeLDFstring += 'extern "C" std::vector<double> jcalc(std::vector<double> listvalue, int curstate)\n'
        ComputeLDFstring += "{\n"
        ComputeLDFstring += "    std::vector<double> ret;\n"
    ComputeLDFstring += "    if (curstate == " + str(loop+1) + ")\n"
    ComputeLDFstring += "    {\n"
    for i in range (len(dvl)):
            tempstring = "        double "+dvl[i] + "= listvalue[" + str(i) + "];\n"
            ComputeLDFstring += tempstring
    ComputeLDFstring += "        double Entry = 0;\n"
    for i in range (funlen):
        for j in range (varlen):
            tempstring="        "+Ljac[i][j]+";\n"
            ComputeLDFstring+=tempstring
            ComputeLDFstring+="        ret.push_back(Entry);\n"
    ComputeLDFstring +="        return ret;\n"
    ComputeLDFstring += "    }\n"
    filename = ""
    filename +="../work-dir/jaThin.cpp"
    if loop ==0:
        savefile=open(filename,'w')
        savefile.write(ComputeLDFstring)
    else:
        savefile= open(filename, "a")
        savefile.write(ComputeLDFstring)
    savefile.close()

    return delete_element

def createCDFfunction(delete_element):
    ComputeLDFstring = ""
    ComputeLDFstring += "\n"
    with open("../work-dir/ThinVarProp",'r') as thinvarfile:
        thinprop = thinvarfile.readlines()
    thinvarfile.close()
    print("thinprop: "+str(thinprop))
    thinvarlist = []
    for i in range(len(thinprop)):
        if thinprop[i][0] == '1':
            thinvarlist.append(i)
    print(thinvarlist)
    ComputeLDFstring += "}\n\n"
    ComputeLDFstring += 'extern "C" void thinHandle(int state, std::vector<double> &notbloating, std::vector<double> &thin_indicate, std::vector<double> &bloating, std::vector<double> &not_thin)\n'
    ComputeLDFstring += "{\n"
    for i in range (len(delete_element)):
        ComputeLDFstring +="    if (int(state) == "
        ComputeLDFstring += str(i+1)
        ComputeLDFstring += ")\n"
        ComputeLDFstring += "    {\n"
        ComputeLDFstring += "        notbloating=std::vector<double>{"
        for j in range (len(delete_element[i])):
            ComputeLDFstring+=str(delete_element[i][j])
            if j!= len(delete_element[i])-1:
                ComputeLDFstring +=","
        ComputeLDFstring+="};\n" 

        ComputeLDFstring += "        thin_indicate=std::vector<double>{"
        emptyflag = 1
        for j in range (len(thinvarlist)):
            if thinvarlist[j] not in delete_element[i]:
                if emptyflag == 0:
                    ComputeLDFstring +=","
                ComputeLDFstring+=str(thinvarlist[j])
                emptyflag = 0
        ComputeLDFstring+="};\n"

        ComputeLDFstring += "        bloating=std::vector<double>{"
        emptyflag = 1
        for j in range (len(thinprop)):
            if j not in delete_element[i]:
                if emptyflag == 0:
                    ComputeLDFstring +=","
                ComputeLDFstring+=str(j)
                emptyflag = 0
        ComputeLDFstring+="};\n"

        ComputeLDFstring += "        not_thin=std::vector<double>{"
        emptyflag = 1
        for j in range (len(thinprop)):
            if (j not in thinvarlist) and (j not in delete_element[i]):
                if emptyflag == 0:
                    ComputeLDFstring +=","
                ComputeLDFstring+=str(j)
                emptyflag = 0
        ComputeLDFstring+="};\n"
        ComputeLDFstring+="        return;\n"
        ComputeLDFstring += "    }\n"
        
    ComputeLDFstring += "}\n"
    filename = ""
    filename +="../work-dir/jaThin.cpp"
    savefile=open(filename,'a')
    savefile.write(ComputeLDFstring)
    savefile.close()