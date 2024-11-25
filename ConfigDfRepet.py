#!/usr/bin/python
"""*********************************************************************
Rusty
Copyright: 2024 INRAE https://www.inrae.fr
License:
  CeCILL: http://www.cecill.info/index.fr.html
  See the LICENCE file in the project's top-level directory for details.
Authors:
  * Emma CORRE, Ecogenomics of interactions team, IAM

************************************************************************"""
"""
Change the REPET_Summary in a proper output with multiple specie
in : specie list (in the same directory)
out: DataFrame of all specie 
"""
import pandas as pd

def readpd(summarfile, name, lengenome):
  columns = ["TE", "length","covg", "frags", "fullLgthFrags","copies", "fullLgthCopies", "meanId", "sdId", "clusMCL", "Seqname", "Wcode", "strand", "confused", "classi", "order", "sFamily", "CI", "untruc", "coding", "struct", "other"]
  df = pd.DataFrame()
  df=pd.read_csv(summarfile, names = columns, header=0,  sep="\t")
  print(name)
  rowshape = int(df.shape[0])
  print(df.shape)
  df.insert(22, "specie", [name for i in range(0,rowshape)], True)
  df.insert(23, "lengenome", [lengenome for i in range(0, rowshape)], True)
  info=""
  list_info=[]
  for i,j in zip(df.coding, df.confused):
    print(i,j)
    identifiant = i.split(":")
    if identifiant[0] == "coding=(NA)" or identifiant[0] == " coding=(NA)":
      info = "NA"
    elif identifiant[0] == "coding=(profiles":
      info="NA"
    else: 
      if j == "True" or j =="ok":
        print("jesuis l identifiant")
        print(identifiant)
        identifiant = identifiant[1]
        info = identifiant.split("-")[0]
      else:
        info="NA"
    list_info.append(info)
  list_order = []
  for i in df.order:
    if i in ['Unclassified', 'TIR', 'LTR', 'MITE', 'TRIM', 'Helitron', 'DIRS', 'LINE', 'SINE','Maverick','Ple','LARD','noCat']:
      list_order.append(i)
    else:
      list_order.append("Conflict")
      if len(str(i).split("|")) == 1:
        print(i)
  list_classi = []
  for i in df.classi:
    if i in ["I", "II", "Unclassified"]:
      list_classi.append(i)
    else:
      list_classi.append("Conflict")
  df.insert(24,"SF",[info for info in list_info],True)
  df.insert(24,"Classif",[info for info in list_classi],True)
  df.insert(24,"Orderi",[info for info in list_order],True)
  return df

def cheminacces():
  Pg9C = "/Users/ecorre/Documents/data4Zenodoo/REPET_output/TEsummary/Pg9C_annotTE_summary.txt"
  lengPg9C = 86818685
  Mlp23 = "/Users/ecorre/Documents/data4Zenodoo/REPET_output/TEsummary/Mellp2_annotTE_summary.txt"
  lengMlp23 = 109877997
  Pg9A = "/Users/ecorre/Documents/data4Zenodoo/REPET_output/TEsummary/Pg9A_annotTE_summary.txt"
  lenPg9A = 91841532
  MeliCB = "/Users/ecorre/Documents/data4Zenodoo/REPET_output/TEsummary/MeliCconcatecleaned2.txt"
  lenMeliCB = 239950640
  Pg201A = "/Users/ecorre/Documents/data4Zenodoo/REPET_output/TEsummary/Pg201A_annotTE_summary.txt"
  lenPg201A = 88783524
  Pg201B = "/Users/ecorre/Documents/data4Zenodoo/REPET_output/TEsummary/Pg201B_annotTE_summary.txt"
  lenPg201B = 89387476
  Micin = "/Users/ecorre/Documents/data4Zenodoo/REPET_output/TEsummary/Micin_annotTE_summary.txt"
  lenMicin = 23461035
  Lcre = "/Users/ecorre/Documents/data4Zenodoo/REPET_output/TEsummary/Lcre_annotTE_summary.txt"
  lenLcre = 26329764
  Mylid = "/Users/ecorre/Documents/data4Zenodoo/REPET_output/TEsummary/Mylid_annotTE_summary.txt"
  lenMylid = 33172729
  Roto = "/Users/ecorre/Documents/data4Zenodoo/REPET_output/TEsummary/Roto_annotTE_summary.txt"
  lenRoto =  20683845
  Mame = "/Users/ecorre/Documents/data4Zenodoo/REPET_output/TEsummary/Mame_annotTE_summary.txt"
  lenMame = 112350849
  MeliH = "/Users/ecorre/Documents/data4Zenodoo/REPET_output/TEsummary/MeliHconcatecleaned2.txt"
  lenMeliH = 242316735
  Melap = "/Users/ecorre/Documents/data4Zenodoo/REPET_output/TEsummary/Melap25_annotTE_summary.txt"
  lenMelap = 335730080
  AZ2A = "/Users/ecorre/Documents/data4Zenodoo/REPET_output/TEsummary/AZ2A_annotTE_summary.txt"
  lenAZ2A = 75333134
  AZ2B = "/Users/ecorre/Documents/data4Zenodoo/REPET_output/TEsummary/AZ2B_annotTE_summary.txt"
  lenAZ2B = 75649056
  Phord = "/Users/ecorre/Documents/data4Zenodoo/REPET_output/TEsummary/Phord_annotTE_summary.txt"
  lenPhord = 206919034
  Pt76A = "/Users/ecorre/Documents/data4Zenodoo/REPET_output/TEsummary/Pt76A_annotTE_summary.txt"
  lenPt76A=123569637
  Pt76B = "/Users/ecorre/Documents/data4Zenodoo/REPET_output/TEsummary/Pt76B_annotTE_summary.txt"
  lenPt76B=128930402
  OldMT2006 = "/Users/ecorre/Documents/data4Zenodoo/REPET_output/TEsummary/MT2006New.csv"
  lenOldMT2006=1057430280

  dfRoto = readpd(Roto, "Rhodosporidium toruloides IFO0880", lenRoto)
  dfLcre = readpd(Lcre, "Leucosporidiella creatinivora UCDFST 62-1032 v1.0", lenLcre)
  dfMylid = readpd(Mylid, "Microbotryum lychnidis-dioicae p1A1 Lamole", lenMylid)
  dfMicin = readpd(Micin, "Microbotryum intermedium 1389 BM 12 12", lenMicin)
  dfMlp23 = readpd(Mlp23, "Melampsora larici-populina v2.0",lengMlp23)
  dfMeliCB = readpd(MeliCB, "Melampsora lini haplotype C", lenMeliCB)
  dfPg9C = readpd(Pg9C, "Puccinia graminis f. sp. tritici Ug99 haplotype C",lengPg9C)
  dfPg9A = readpd(Pg9A, "Puccinia graminis f. sp. tritici Ug99 haplotype A", lenPg9A)
  dfPg201B = readpd(Pg201B, "Puccinia graminis f. sp. tritici 21-0 haplotype B", lenPg201B)
  dfPg201A = readpd(Pg201A, "Puccinia graminis f. sp. tritici 21-0 haplotype A", lenPg201A)
  dfMame = readpd(Mame, "Melampsora americana R15-033-03 v1.0", lenMame)
  dfMeliH = readpd(MeliH, "Melampsora lini haplotype H", lenMeliH)
  dfMelap = readpd(Melap,"Melampsora allii- populina 12AY07 v1.0", lenMelap)
  dfAZ2A = readpd(AZ2A," Puccinia striiformis f. sp. tritici haplotype A", lenAZ2A)
  dfAZ2B = readpd(AZ2B," Puccinia striiformis f. sp. tritici haplotype B", lenAZ2B)
  # dfPt76 = readpd(Pt76,"Puccinia triticina", lenPt76)
  dfPt76A = readpd(Pt76A,"Puccinia triticina A", lenPt76A)
  dfPt76B = readpd(Pt76B,"Puccinia triticina Btig", lenPt76B)
  dfPhord = readpd(Phord,"Puccinia hordei isolate 560", lenPhord)
  print("ok Phord")
  dfPOldMT2006 = readpd(OldMT2006,"Phakopsora pachyrhizi MT2006", lenOldMT2006)

  
  print("ok")
  dfconcat=pd.concat([dfRoto, dfMicin, dfMylid, dfLcre,dfAZ2A, dfAZ2B, dfPt76A,dfPt76B, dfPg201A, dfPg9A,dfPg201B, dfPg9C,dfPhord, dfMame, dfMlp23, dfMeliCB,dfMeliH, dfMelap, dfPOldMT2006])
  listdf=""
  return dfconcat

test = cheminacces()
list_order = (test.Orderi.unique()).tolist()
list_classi = (test.Classif.unique()).tolist()
list_sf = (test.SF.unique()).tolist()
print(list_sf)
test.to_csv("/Users/ecorre/Documents/info_complete_09_10_2024.csv")