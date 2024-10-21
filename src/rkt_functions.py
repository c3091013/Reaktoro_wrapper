import colorsys
import os
import subprocess
import math as m
import myfunctions.functions as mf
import myfunctions.molecule as ml
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import math as m
from reaktoro import *
from scipy.stats import linregress
from itertools import chain, combinations
import multiprocessing as mp

#wrapper to allow parralelisation of functions used in this script. input is a list with the name of the function followed by a list of the inputs 
def wrapper(lst):
    if lst[0] == "scan_1D_pH": return scan_1D_pH(*lst[1])
    if lst[0] == "scan_1D_Eh": return scan_1D_Eh(*lst[1])
    if lst[0] == "scan_1D_pC": return scan_1D_pC(*lst[1])
    if lst[0] == "scan_2D_xpH_yEh": return scan_2D_xpH_yEh(*lst[1])
    if lst[0] == "scan_2D_xpH_ypC": return scan_2D_xpH_ypC(*lst[1])
    if lst[0] == "scan_2D_xpC_yEh": return scan_2D_xpC_yEh(*lst[1])
    if lst[0] == "scan_2D_xpC_ypC": return scan_2D_xpC_ypC(*lst[1])
    return "missed"

#generates a 1D plot of data and casts it to an axis object from matplotlib
def get_1D_plot(ax,data,xlim,imask,title,xtit,ytit,colors):
    mask=list(imask & data.keys())
    for ii in mask:
        if data[ii].max() > 1E-2:
            ax.plot(data["x_vals"],data[ii],label=ii,color=colors[ii])
    ax.set_title(title,loc='left',fontsize=20)
    ax.set_xlabel(xtit,fontsize=20)
    ax.set_ylabel(ytit,fontsize=20)
    ax.set_xlim(*xlim)
    ax.set_ylim([0.0,1.01])
    ax.tick_params(axis='both', labelsize=20)
    return ax.get_legend_handles_labels()

#generates a 2D patches plotfor a given set of predominance data and casts it to an axis object.
def get_2D_plot(ax,data,xlim,ylim,imask,title,xtit,ytit,colors):
    mask=list(imask & data.keys())
    a=[]
    cl=[]
    for ii in range(0,len(data["y_vals"])):
        b=[]
        c=[]
        for jj in range(0,len(data["x_vals"])):
            tc=[]
            for kk in range(0,len(mask)):
                tc.append(data[mask[kk]][jj][ii])
            tc=np.array(tc)
            b.append(colors[mask[np.argmax(tc)]])
            c.append(np.argmax(tc))
        a.append(b)
        cl.append(c)
    cl=np.array(cl)
    values = np.unique(cl.ravel())
    patches = [ mpatches.Patch(color=colors[mask[values[i]]], label="{l}".format(l=mask[values[i]]) ) for i in range(len(values)) ]
    im = ax.imshow(a,extent=[data["x_vals"].min(), data["x_vals"].max(), data["y_vals"].min(), data["y_vals"].max()],interpolation='antialiased', origin ='lower',aspect=((data["x_vals"].max()-data["x_vals"].min())/(data["y_vals"].max()-data["y_vals"].min())))
    #ax.legend(handles=patches, borderaxespad=0.)#, bbox_to_anchor=(1.05, -0.05), loc='upper left')
    ax.set_title(title,loc='left',fontsize=20)
    ax.set_xlabel(xtit,fontsize=20)
    ax.set_ylabel(ytit,fontsize=20)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.tick_params(axis='both', labelsize=20)
    #ax.tick_params(axis='both',which='both',bottom=False,top=False,left=False,right=False,labelbottom=False,labelleft=False)
    return patches

def get_phases(s_dic):
    aq_spe=[]
    g_spe=[]
    s_spe=[]
    for ii in s_dic.keys():
        match s_dic[ii][4]:
            case "aq":
                aq_spe.append(ii)
            case "s":
                s_spe.append(ii)
            case "g":
                g_spe.append(ii)
            case _:
                g_spe.append(ii)
    #now i want to build each phases, first aq, then g, then solids as their own set of phases
    phases=[]
    if aq_spe:
        phases.append(AqueousPhase(aq_spe))
    if g_spe:
        phases.append(GaseousPhase(g_spe))
    #now the solids as seperate phases
    for ii in s_spe:
        phases.append(SolidPhase(ii))
    return phases


def equilibrate(s_dic,c_lst,TPEH):
    #build out database
    db=get_DB(s_dic,TPEH[0])
    phases = get_phases(s_dic)
    system = ChemicalSystem(db, *phases)
    state = ChemicalState(system)
    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.Eh()
    specs.pH()
    solver = EquilibriumSolver(specs)
    for ii in c_lst:
        state.set(*ii)
    conditions = EquilibriumConditions(specs)
    conditions.temperature(TPEH[0], "celsius")
    conditions.pressure(TPEH[1], "bar")
    conditions.Eh(TPEH[2])
    conditions.pH(TPEH[3])
    options = EquilibriumOptions()
    options.optima.maxiters=10000
    solver.setOptions(options)
    solver.solve(state, conditions)
    return state

def scan_2D_xpH_yEh(x_vals,y_vals,s_dic,c_lst,TP):
    concs={}
    concs["x_vals"]=np.array(x_vals)
    concs["y_vals"]=np.array(y_vals)
    for ii in s_dic.keys():
        concs[ii]=[]
    for ii in x_vals:
        tmp={}
        for jj in s_dic.keys():
            tmp[jj]=[]
        for jj in y_vals:
            state = equilibrate(s_dic,c_lst,[TP[0],TP[1],jj,ii]) 
            for kk in s_dic.keys():
                tmp[kk].append(float(state.speciesAmount(kk)))
        for jj in s_dic.keys():
            concs[jj].append(tmp[jj])
    for ii in s_dic.keys():
        concs[ii] = np.array(concs[ii])
    return concs

def scan_2D_xpH_ypC(x_vals,y_vals,nm,s_dic,c_lst,TPE):
    concs={}
    concs["x_vals"]=np.array(x_vals)
    concs["y_vals"]=np.array(y_vals)
    for ii in s_dic.keys():
        concs[ii]=[]
    for ii in x_vals:
        tmp={}
        for jj in s_dic.keys():
            tmp[jj]=[]
        for jj in y_vals:
            state = equilibrate(s_dic,c_lst+[[nm,m.pow(10,jj),'mol']],[TPE[0],TPE[1],TPE[2],ii])
            for kk in s_dic.keys():
                tmp[kk].append(float(state.speciesAmount(kk)))
        for jj in s_dic.keys():
            concs[jj].append(tmp[jj])
    for ii in s_dic.keys():
        concs[ii] = np.array(concs[ii])
    return concs

def scan_2D_xpC_yEh(x_vals,y_vals,nm,s_dic,c_lst,TPH):
    concs={}
    concs["x_vals"]=np.array(x_vals)
    concs["y_vals"]=np.array(y_vals)
    for ii in s_dic.keys():
        concs[ii]=[]
    for ii in x_vals:
        tmp={}
        for jj in s_dic.keys():
            tmp[jj]=[]
        for jj in y_vals:
            state = equilibrate(s_dic,c_lst+[[nm,m.pow(10,ii),'mol']],[TPH[0],TPH[1],jj,TPH[2]])
            for kk in s_dic.keys():
                tmp[kk].append(float(state.speciesAmounts(kk)))
        for jj in s_dic.keys():
            concs[jj].append(tmp[jj])
    for ii in s_dic.keys():
        concs[ii] = np.array(concs[ii])
    return concs

def scan_2D_xpC_ypC(x_vals,y_vals,nmx,nmy,s_dic,c_lst,TPEH):
    concs={}
    concs["x_vals"]=np.array(x_vals)
    concs["y_vals"]=np.array(y_vals)
    for ii in s_dic.keys():
        concs[ii]=[]
    for ii in x_vals:
        tmp={}
        for jj in s_dic.keys():
            tmp[jj]=[]
        for jj in y_vals:
            state = equilibrate(s_dic,c_lst+[[nmx,m.pow(10,ii),'mol'],[nmy,m.pow(10,jj),'mol']],[TPEH[0],TPEH[1],TPEH[2],TPEH[3]])
            for kk in s_dic.keys():
                tmp[kk].append(float(state.speciesAmount(kk)))
        for jj in s_dic.keys():
            concs[jj].append(tmp[jj])
    for ii in s_dic.keys():
        concs[ii] = np.array(concs[ii])
    return concs

#scan concentration change as a function of Eh
def scan_1D_Eh(x_vals,s_dic,c_lst,TPH):
    concs={}
    concs["x_vals"]=x_vals
    for ii in s_dic.keys():
        concs[ii]=[]
    for ii in x_vals:
        state = equilibrate(s_dic,c_lst,[TPH[0],TPH[1],ii,TPH[2]])
        for jj in s_dic.keys():
            concs[jj].append(float(state.speciesAmount(jj)))
    for ii in s_dic.keys():
        concs[ii]=np.array(concs[ii])
    return concs

#scan concentration change as a function of pH
def scan_1D_pH(x_vals,s_dic,c_lst,TPE):
    concs={}
    concs["x_vals"]=x_vals
    for ii in s_dic.keys():
        concs[ii]=[]
    for ii in x_vals:
        state = equilibrate(s_dic,c_lst,[TPE[0],TPE[1],TPE[2],ii])
        for jj in s_dic.keys():
            concs[jj].append(float(state.speciesAmount(jj)))
    for ii in s_dic.keys():
        concs[ii]=np.array(concs[ii])
    return concs

#scan concentration change as a function of some initial concentration
def scan_1D_pC(x_vals,nm,s_dic,c_lst,TPEH):
    concs={}
    concs["x_vals"]=x_vals
    for ii in s_dic.keys():
        concs[ii]=[]
    for ii in x_vals:
        state = equilibrate(s_dic,c_lst+[[nm,m.pow(10,ii),'mol']],[TPEH[0],TPEH[1],TPEH[2],TPEH[3]])
        for jj in s_dic.keys():
            concs[jj].append(float(state.speciesAmount(jj)))
    for ii in s_dic.keys():
        concs[ii]=np.array(concs[ii])
    return concs

def get_DB(specs,temp):
    species = SpeciesList()
    #we want to build a data base out of species information, including their names, and data. take input as a dictionary. first i want all the elements i will be dealing with.
    el_lst=[]
    for ii in specs.keys():
        if specs[ii][0]:
            for jj in specs[ii][1].keys():
                if jj not in el_lst:
                    el_lst.append(jj)

    elements=get_EL(el_lst) 
    #we have out element list, lets deal with primary species.
    for ii in specs.keys():
        if specs[ii][0]:
            els=[]
            for jj in specs[ii][1].keys():
                els.append((elements.get(jj),  specs[ii][1][jj]))
            species.append(Species()
                    .withName(ii)
                    .withFormula(ii)
                    .withCharge(specs[ii][2])
                    .withElements(ElementalComposition(els))
                    .withAggregateState(get_agstate(specs[ii][4]))
                    .withStandardGibbsEnergy(specs[ii][3]))

    #secondaries are actually going to be harder
    for ii in specs.keys():
        if not specs[ii][0]:
            eln={}
            els=[]
            rsp=[]
            for jj in specs[ii][1].keys():
                rsp.append((species.get(jj) , specs[ii][1][jj]))
                fac=specs[ii][1][jj]
                s_slt=species.get(jj).elements().symbols()
                c_clt=species.get(jj).elements().coefficients()
                for kk in range(0,len(s_slt)):
                    if s_slt[kk] in eln.keys():
                        eln[s_slt[kk]]=eln[s_slt[kk]]+fac*c_clt[kk]
                    else:
                        eln[s_slt[kk]]=fac*c_clt[kk]
            for jj in eln.keys():
                els.append((elements.get(jj),  eln[jj]))
            species.append(Species()
                    .withName(ii)
                    .withFormula(ii)
                    .withCharge(specs[ii][2])
                    .withElements(ElementalComposition(els))
                    .withAggregateState(get_agstate(specs[ii][4]))
                    .withFormationReaction(FormationReaction()
                        .withReactants(rsp)
                        .withEquilibriumConstant(specs[ii][3])))

    return Database(species.data())

def get_EL(lst):
    elements = ElementList()
    mol_mas = get_MM()
    for en in lst:
        elements.append(Element().withName(en).withSymbol(en).withMolarMass(mol_mas[en]))

    return elements

def get_MM():
    mass={}
    mass["H"]=1.0E-3
    mass["O"]=16.0E-3
    mass["C"]=12.0E-3
    mass["N"]=14.0E-3
    mass["S"]=32.0E-3
    mass["Co"]=58.9E-3
    mass["Ni"]=58.7E-3
    mass["Cu"]=63.5E-3
    mass["Zn"]=65.4E-3
    mass["Na"]=23.0E-3
    mass["Ca"]=40.1E-3
    mass["Mn"]=54.9E-3
    mass["Mg"]=24.3E-3
    mass["Au"]=197.0E-3
    mass["Fe"]=55.8E-3
    mass["As"]=74.9E-3
    mass["Hg"]=200.6E-3
    return mass

def get_agstate(nm):
    match nm:
        case "aq":
            st = AggregateState.Aqueous
        case "s":
            st = AggregateState.Solid
        case "g":
            st = AggregateState.Gas
        case _:
            st = AggregateState.Aqueous
    return st

def get_N_HexRGBCol(mask):
    N=len(mask)
    X=int(np.ceil(N**(1/3)))
    rgb_out=[[x*1.0/X,y*1.0/X,z*1.0/X] for x in range(X) for y in range(X) for z in range(X)]
    rm.shuffle(rgb_out)
    rgb={}
    for ii in range(0,len(mask)):
        rgb[mask[ii]]=rgb_out[ii]
    return rgb

def get_N_HexRGBCol_old2(N):
    X=int(np.ceil(N**(1/3)))
    rgb_out = [[x*1.0/X,y*1.0/X,z*1.0/X] for x in range(X) for y in range(X) for z in range(X)]
    hex_out=[]
    rm.shuffle(rgb_out)
    return hex_out, rgb_out

def get_N_HexRGBCol_old(N):
    HSV_tuples = [(x*1.0/N,0.4,0.8) for x in range(N)]
    hex_out = []
    rgb_out = []
    for rgb in HSV_tuples:
        rgbhex = map(lambda x: int(x*255),colorsys.hsv_to_rgb(*rgb))
        #hex_out.append('#%02x%02x%02x' % tuple(rgbhex))
        hex_out.append(list(rgbhex))
        rgbnum = map(lambda x: int(x*255),colorsys.hsv_to_rgb(*rgb))
        #rgb_out.append('%d %d %d' % tuple(rgbnum))
        rgb_out.append(list(rgbnum))
    #print(rgb_out)
    return hex_out, rgb_out


