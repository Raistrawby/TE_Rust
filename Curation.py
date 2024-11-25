#!/bin/python
import pandas as pd
import numpy as np
from Bio import SeqIO
import argparse
"""*********************************************************************
Copyright: 2024 INRAE https://www.inrae.fr
License:
  CeCILL: http://www.cecill.info/index.fr.html
  See the LICENCE file in the project's top-level directory for details.
Authors:
  * Emma CORRE, Ecogenomics of interactions team, IAM
************************************************************************"""
"""
Curation after MCHelper and CDHIT of REPET annotation
in : REPET FILE, MCHelper FILE, CDHIT FILE (in the same directory)
out: DataFrame of all specie corrected
"""
def argParser():
	"""
	Arguments managment
	:return: an ArgumentParser object
	:rtype: object
	.. sectionauthor:: Emma Corre
	"""
	parser = argparse.ArgumentParser()
	parser.add_argument("-R", "--summary_repet",
	    help = "path of the summary.csv corrected file built by REPET",
		required = True)
	parser.add_argument("-M", "--MCHelper",
	    help = "path of the MCHelper file",
		required = True)
	parser.add_argument("-C", "--CDHIT",
	    help = "path of the cdhit file",
		required = True)
	return parser.parse_args()



def readpd(summarfile):
    """
    Reads the summary file and returns a DataFrame.
    Parameters:
        summarfile (str): Path to the summary file.
    Returns:
        pd.DataFrame: DataFrame with the summary data.
    """
    summarfile = pd.read_csv(summarfile, header=0, sep=",")
    summarfile[['MC_class', 'MC_order', 'MC_SF', 'MC_info']] = "non"
    # print(summarfile.columns)
    return summarfile


def readMCHelper(MCHelperF):
    """
    Reads the MCHelper file and returns a list.
    Parameters:
        summarfile (str): Path to the summary file.
    Returns:
        list: list with recordid.
    """
    list_recordid=[]
    for record in SeqIO.parse(MCHelperF, "fasta"):
     list_recordid.append(record.id)
    return list_recordid


def parseinfoMC(infolenClass):
    """
    Reads the MCHelper file and returns a DataFrame.
    Parameters:
        summarfile (str): Path to the summary file.
    Returns:
        pd.DataFrame: DataFrame with the summary data.
    """
    match len(infolenClass):
        case 1:
           return 'RAS'
        case 2:
          if infolenClass[1].split("#")[0].isdigit():
             return "split"
          else:
             return infolenClass[1].split("#")[0]
        case _:
            return "toto"
    return None


def parserelou(relouID):
   match relouID:
        case _ if  'denovo' in relouID:
           return 'RAS'
        case _ if relouID.isdigit():
          return 'split'
        case _:
            return relouID
   pass


def parseMCHelperId(listID):
    """
    Reads the MCHelper file and returns a DataFrame.
    Parameters:
        summarfile (str): Path to the summary file.
    Returns:
        pd.DataFrame: DataFrame with the summary data.
    """
    dicoIDinfo = {}
    for i in listID:
       classification = i.split('#')[1]
       if len(i.split("Map")) == 1:
            testi=parserelou(i.split("#")[0].split("_")[-1])
            infoMC = testi
       else:
            infoWhy = i.split('Map')[1].split('_')
            infoWhyText=parseinfoMC(infoWhy)
       infolenClass = classification.split('/')
       ID= i.split('#')[0]
       match len(infolenClass):
            case 1:
                classMC = infolenClass[0]
                orderMC = 'Nan'
                SFMC = 'Nan'
                infoMC = infoWhyText
                dicoIDinfo[ID] = {'classMC': classMC, 'orderMC': orderMC, 'SFMC' : SFMC, 'infoMC': infoMC}
            case 2:
                classMC = infolenClass[0]
                orderMC = infolenClass[1]
                SFMC = 'Nan'
                infoMC = infoWhyText
                dicoIDinfo[ID] = {'classMC': classMC, 'orderMC': orderMC, 'SFMC' : SFMC, 'infoMC': infoMC}
            case 3:
                classMC = infolenClass[0]
                orderMC = infolenClass[1]
                SFMC = infolenClass[2]
                infoMC = infoWhyText
                dicoIDinfo[ID] = {'classMC': classMC, 'orderMC': orderMC, 'SFMC' : SFMC, 'infoMC': infoMC}
            case _:
                return "Aucun modèle ne correspond"
    return dicoIDinfo


def joininfo(dicoinfo, dfRepet):
    for i,j in dicoinfo.items():
        if str(i) in dfRepet['TE'].values:
            dfRepet.loc[dfRepet['TE']== i, 'MC_class']=j['classMC']
            dfRepet.loc[dfRepet['TE']== i, 'MC_order']=j['orderMC']
            dfRepet.loc[dfRepet['TE']== i, 'MC_SF']=j['SFMC']
            dfRepet.loc[dfRepet['TE']== i, 'MC_info']=j['infoMC']
        else: 
            mask = dfRepet['TE'].str.contains(i)
            if mask.any(): # P without the end
                dfRepet.loc[dfRepet['TE']== mask, 'MC_class']=j['classMC']
                dfRepet.loc[dfRepet['TE']== mask, 'MC_order']=j['orderMC']
                dfRepet.loc[dfRepet['TE']== mask, 'MC_SF']=j['SFMC']
                dfRepet.loc[dfRepet['TE']== mask, 'MC_info']=j['infoMC']
            else:
                if len(i.split('_s_')) >1 :
                    mask = dfRepet['TE'].str.contains(i.split('_s_')[0]) 
                    if mask.any(): # P without the end
                        dfRepet.loc[dfRepet['TE']== mask, 'MC_class']=j['classMC']
                        dfRepet.loc[dfRepet['TE']== mask, 'MC_order']=j['orderMC']
                        dfRepet.loc[dfRepet['TE']== mask, 'MC_SF']=j['SFMC']
                        dfRepet.loc[dfRepet['TE']== mask, 'MC_info']=j['infoMC']
                else:
                    presplit = i.split('_s_')[0].split('_')
                    result = '_'.join(presplit[:-1])
                    mask = dfRepet['TE'].str.contains(result) 
                    if mask.any(): # P without the end
                        dfRepet.loc[dfRepet['TE']== result, 'MC_class']=j['classMC']
                        dfRepet.loc[dfRepet['TE']== result, 'MC_order']=j['orderMC']
                        dfRepet.loc[dfRepet['TE']== result, 'MC_SF']=j['SFMC']
                        dfRepet.loc[dfRepet['TE']== result, 'MC_info']=j['infoMC']
                    elif str(result) in dfRepet['TE'].values:
                        dfRepet.loc[dfRepet['TE']== result, 'MC_class']=j['classMC']
                        dfRepet.loc[dfRepet['TE']== result, 'MC_order']=j['orderMC']
                        dfRepet.loc[dfRepet['TE']== result, 'MC_SF']=j['SFMC']
                        dfRepet.loc[dfRepet['TE']== result, 'MC_info']=j['infoMC']
                    else:
                        print(result)
                        raise Exception("oulalala relou bis")
    return dfRepet


def readpd2(summarfile):
    """
    Reads the summary file and returns a DataFrame.
    Parameters:
        summarfile (str): Path to the summary file.
    Returns:
        pd.DataFrame: DataFrame with the summary data.
    """
    summarfile = pd.read_csv(summarfile, header=0, sep=",")
    summarfile[['cdhit_class', 'cdhit_order', 'cdhit_SF']] = "non"
    return summarfile


def create_dicoCluster(FileIDSubset):
    dico_CDHit={}
    with open(FileIDSubset) as FileI:
        for line in FileI:
            if line[0]==">":
                identifiant=line.split("\n")[0].split(">Cluster ")[1]
                dico_CDHit[identifiant]=[]
            else:
                seqCl=line.split("\n")[0].split("\t")[1].split("...")[0].split("nt,")[1].split(">")[1]
                dico_CDHit[identifiant].append(seqCl)
    return dico_CDHit

def clustercompClass(info_cluster,TEID):
    CDHITClassif = 'non'
    query_string = 'TE in @info_cluster'
    result_df = dfRepet.query(query_string)
    result_dfRepet = result_df[~result_df['Classif'].isin(['Unclassified', 'Conflict','?', ''])]
    classifUniq = result_dfRepet.Classif.unique()
    if len(classifUniq) == 1:
        CDHITClassif = classifUniq[0]
    else:
        result_dfMC = result_df[~result_df['MC_class'].isin(['Unclassified', 'Conflict','non'])]
        classifUniqMC = result_dfMC.MC_class.unique()
        if len(classifUniqMC) == 1:
            CDHITClassif = classifUniqMC[0]
        else:
            CDHITClassif = 'nonDiffClass'
    return CDHITClassif


def clustercompOrder(info_cluster,TEID):
    CDHITorder = 'non'
    query_string = 'TE in @info_cluster'
    result_df = dfRepet.query(query_string)
    result_dfRepet = result_df[~result_df['Orderi'].isin(['Unclassified', 'Conflict','?',''])]
    orderUniq = result_dfRepet.Orderi.unique()
    if len(orderUniq) == 1:
        CDHITorder = orderUniq[0]
    else:
        result_dfMC = result_df[~result_df['MC_order'].isin(['Unclassified', 'Conflict','non'])]
        orderUniqMC = result_dfMC.MC_order.unique()
        if len(orderUniqMC) == 1:
            CDHITorder = orderUniqMC[0]
        else:
            CDHITorder = 'nonDiffOrder'
    return CDHITorder


def clustercompSF(info_cluster,TEID):
    CDHITSF = 'non'
    query_string = 'TE in @info_cluster'
    result_df = dfRepet.query(query_string)
    result_dfRepet = result_df[~result_df['SF'].isin(['Unclassified', 'Conflict','?',''])]
    SFUniq = result_dfRepet.SF.unique()
    if len(SFUniq) == 1 and SFUniq[0] != '':
        CDHITSF = SFUniq[0]
    else:
        result_dfMC = result_df[~result_df['MC_SF'].isin(['Unclassified', 'Conflict','non'])]
        SFUniqMC = result_dfMC.MC_SF.unique()
        if len(SFUniqMC) == 1 and SFUniqMC[0] != '':
            CDHITSF = SFUniqMC[0]
        else:
            CDHITSF = 'nonDiffSF'
    return CDHITSF


def clustercomp(info_cluster,TEID):
    CDHITClassif = clustercompClass(info_cluster,TEID)
    if CDHITClassif != 'non':
        CDHITorder = clustercompOrder(info_cluster,TEID)
        if CDHITorder != 'non':
            CDHITSF = clustercompSF(info_cluster,TEID)
            if CDHITSF != 'non':
                CDHITSF = clustercompSF(info_cluster,TEID)
            else:
                CDHITSF = 'non'
        else:
            CDHITorder = 'non'
            CDHITSF = 'non'
    else:
        CDHITClassif = 'non'
        CDHITorder = 'non'
        CDHITSF = 'non'
    return(CDHITClassif,CDHITorder,CDHITSF)


def annotCluster(dico_CDHit,dfRepet):
    for i,j in dico_CDHit.items():
        if len(dico_CDHit[i])==1:
            if not dfRepet.loc[(dfRepet['TE'] == j[0]) & (~dfRepet['Classif'].isin(['Conflict', 'Unclassified','?','']))].empty:
                dfRepet.loc[dfRepet['TE']== j[0], 'cdhit_class']=dfRepet.loc[dfRepet['TE'] == j[0], 'Classif']
            else:
                if not dfRepet.loc[(dfRepet['TE'] == j[0]) & (~dfRepet['MC_class'].isin(['Conflict', 'Unclassified','non']))].empty:
                    dfRepet.loc[dfRepet['TE']== j[0], 'cdhit_class']=dfRepet.loc[dfRepet['TE'] == j[0], 'MC_class']
            if not dfRepet.loc[(dfRepet['TE'] == j[0]) & (~dfRepet['Orderi'].isin(['Conflict', 'Unclassified','?','']))].empty:
                dfRepet.loc[dfRepet['TE']== j[0], 'cdhit_order']=dfRepet.loc[dfRepet['TE'] == j[0], 'Orderi']
            else:
                if not dfRepet.loc[(dfRepet['TE'] == j[0]) & (~dfRepet['MC_order'].isin(['Conflict', 'Unclassified','non','?','']))].empty:
                    dfRepet.loc[dfRepet['TE']== j[0], 'cdhit_class']=dfRepet.loc[dfRepet['TE'] == j[0], 'MC_order']
            if not dfRepet.loc[(dfRepet['TE'] == j[0]) & (~dfRepet['SF'].isin(['Conflict', 'Unclassified','','?']))].empty:
                dfRepet.loc[dfRepet['TE']== j[0], 'cdhit_SF']=dfRepet.loc[dfRepet['TE'] == j[0], 'SF']
            else:
                if not dfRepet.loc[(dfRepet['TE'] == j[0]) & (~dfRepet['MC_SF'].isin(['Conflict', 'Unclassified','non','?','']))].empty:
                    dfRepet.loc[dfRepet['TE']== j[0], 'cdhit_SF']=dfRepet.loc[dfRepet['TE'] == j[0], 'MC_SF']
        else:
            CDHITClassif, CDHITorder, CDHITSF = clustercomp(dico_CDHit[i],i)
            info_cluster = dico_CDHit[i]
            query_string = 'TE in @info_cluster'
            result_df = dfRepet.query(query_string).copy() # sert a eviter les erreurs et pas modif le df de base
            if not result_df.empty:
                dfRepet.loc[dfRepet['TE'].isin(info_cluster), 'cdhit_class'] = CDHITClassif
                dfRepet.loc[dfRepet['TE'].isin(info_cluster), 'cdhit_order'] = CDHITorder
                dfRepet.loc[dfRepet['TE'].isin(info_cluster), 'cdhit_SF'] = CDHITSF
    return dfRepet


def readpd3(summarfile):
    """
    Reads the summary file and returns a DataFrame.
    Parameters:
        summarfile (str): Path to the summary file.
    Returns:
        pd.DataFrame: DataFrame with the summary data.
    """
    summarfile = pd.read_csv(summarfile, header=0, sep=",")
    summarfile[['FinalAnnot']] = "non"
    return summarfile


def parse_info(dfRepetN):
    final_all_values = []
    final_Annot_values = []
    for index, row in dfRepetN.iterrows():
        if row['Classif'] not in ['', '?', 'Conflict', 'Unclassified','non']:
            finalClassif = row['Classif']
        elif row['MC_class'] not in ['', '?', 'Conflict', 'Unclassified','non']:
            finalClassif = row['MC_class']
        elif row['cdhit_class'] not in ['', '?', 'Conflict', 'Unclassified','non']:
            finalClassif = row['cdhit_class']
        else:
            finalClassif = 'Unclassified'

        if row['Orderi'] not in ['', '?', 'Conflict', 'Unclassified','non']:
            finalOrder = row['Orderi']
        elif row['MC_order'] not in ['', '?', 'Conflict', 'Unclassified','non']:
            finalOrder = row['MC_order']
        elif row['cdhit_order'] not in ['', '?', 'Conflict', 'Unclassified','non']:
            finalOrder = row['cdhit_order']
        else:
            finalOrder = 'Unclassified'

        if row['SF'] not in ['', '?', 'Conflict', 'Unclassified','non']:
            finalSF = row['SF']
        elif row['MC_SF'] not in ['', '?', 'Conflict', 'Unclassified','non']:
            finalSF = row['MC_SF']
        elif row['cdhit_SF'] not in ['', '?', 'Conflict', 'Unclassified','non']:
            finalSF = row['cdhit_SF']
        else:
            finalSF = 'Unclassified'

        final_all = f"{row['TE']}#{finalClassif}/{finalOrder}/{finalSF}"
        final_Annot=f"{finalClassif}/{finalOrder}/{finalSF}"
        final_all_values.append(final_all)
        final_Annot_values.append(final_Annot)


    print(f"Length of final_all_values: {len(final_all_values)}")  # Debug: Length check

    dfRepetN['Final_All'] = final_all_values
    dfRepetN['FinalAnnot'] = final_Annot_values
    return dfRepetN
    

   
def main(args):
    print("JOB START")
    info_REPET=args.summary_repet
    dfRepet = readpd(info_REPET)
    MCHelperFile=args.MCHelper
    list_recordid=readMCHelper(MCHelperFile)
    dicoinfo = parseMCHelperId(list_recordid)
    newdfMC=joininfo(dicoinfo, dfRepet)
    newdfMC.to_csv("/Users/ecorre/Documents/tmpMCHELPER.csv",index=False)
    FileIDSubset=args.CDHIT
    info_REPET=("/Users/ecorre/Documents/tmpMCHELPER.csv",index=False)
    dfRepet = readpd2(info_REPET)
    dico_CDHit = create_dicoCluster(FileIDSubset)
    dfRepetN = annotCluster(dico_CDHit,dfRepet)
    dfRepetN.to_csv("/Users/ecorre/Documents/tmpMCHELPERCDHIT.csv",index=False)
    print('ok')
    info_REPET=("/Users/ecorre/Documents/tmpMCHELPERCDHIT.csv")
    dfRepetN = readpd3(info_REPET)
    dfEnd = parse_info(dfRepetN)
    dfEnd.to_csv("/Users/ecorre/Documents/info_completeConflictCorrectedMCHELPERCDHITFINAL_09_10_2024.csv", index=False)
          
if __name__ == '__main__':
	# execute only if run as a script
    args = argParser()
    main(args)
