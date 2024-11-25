#!/usr/bin/python
"""
*********************************************************************
Copyright: 2024 INRAE https://www.inrae.fr
License:
  CeCILL: http://www.cecill.info/index.fr.html
  See the LICENCE file in the project's top-level directory for details.
Authors:
  * Emma CORRE, Ecogenomics of interactions team, IAM

************************************************************************
Main goal is to classify conflicts
"""

import pandas as pd

def readpd(summarfile):
    """
    Reads the summary file and returns a DataFrame.
    Parameters:
        summarfile (str): Path to the summary file.
    Returns:
        pd.DataFrame: DataFrame with the summary data.
    """
    return pd.read_csv(summarfile, header=0, sep=",")

def extract_classification(annotation, threshold=50.0):
    """
    Extracts classification information from the annotation.
    Parameters:
        annotation (str): The annotation provided by REPET.
        threshold (float): The percentage threshold for considering a classification valid.
    Returns:
        tuple: (class, order, superfamily) if valid, otherwise (None, None, None).
    """
    try:
        step2 = annotation.split("TE_BLRtx:")[1]
        step3 = step2.split(";")[0]
        step4 = step3.split(",")
        for teannot in step4:
            fields = teannot.split(":")
            percentage = float(fields[5].split("%")[0])
            if percentage >= threshold:
                classif = fields[1].split("Class")[1]
                order = fields[2]
                superfamily = fields[3]
                return classif, order, superfamily
    except IndexError:
        pass
    return None, None, None

def update_classification(df, column_name):
    """
    Updates classification in the DataFrame based on a specified column.
    Parameters:
        df (pd.DataFrame): DataFrame containing the data to update.
        column_name (str): The column name to process (e.g., 'coding', 'other').
    Returns:
        pd.DataFrame: Updated DataFrame.
    """
    for index, annotation in df[column_name].dropna().items():
        classif, order, superfamily = extract_classification(annotation)
        if classif and order and superfamily:
            df.at[index, 'Classif'] = classif
            df.at[index, 'Orderi'] = order
            df.at[index, 'SF'] = superfamily
    return df

def select_conflict(dfRepet):
    """
    Processes and resolves conflicts in classification.
    Parameters:
        dfRepet (pd.DataFrame): DataFrame with the classification data.
    Returns:
        pd.DataFrame: Updated DataFrame with conflicts resolved.
    """
    # Process coding annotations
    conflictuals = dfRepet
    dfRepet = update_classification(conflictuals, 'coding')
    # Process other annotations
    conflictuals2 = dfRepet
    dfRepet = update_classification(conflictuals2, 'other')
    print("Conflicts remaining after processing:", dfRepet.query('Classif == "Conflict"').shape)
    return dfRepet

def confused_clean(dfRepet):
    """
    Cleans and updates the DataFrame for confused entries.
    Parameters:
        dfRepet (pd.DataFrame): DataFrame with the classification data.
    Returns:
        pd.DataFrame: Updated DataFrame with confused entries processed.
    """
    confusedTE = dfRepet
    print("Number of confused entries:", confusedTE.shape)
    # Process coding annotations
    dfRepet = update_classification(confusedTE, 'coding')
    # Process other annotations
    confusedTE2 = dfRepet
    dfRepet = update_classification(confusedTE2, 'other')
    print("Number of confused entries remaining:", dfRepet.query('confused == True').shape)
    return dfRepet

if __name__ == "__main__":
    Repet_filePath = "/Users/ecorre/Documents/info_complete_09_10_2024.csv"
    dfRepet = readpd(Repet_filePath)
    correctdfRepet = select_conflict(dfRepet)
    print('ok')
    correctdfRepet.drop(columns=correctdfRepet.columns[0], axis=1, inplace=True)
    doublecorrection = confused_clean(correctdfRepet)
    doublecorrection.to_csv("/Users/ecorre/Documents/info_completeConflictCorrected_09_10_2024.csv", index=False)

