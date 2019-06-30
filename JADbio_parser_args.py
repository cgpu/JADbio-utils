#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: cristina

Script for parsing all JAD output .txt files of an input_directory
to collect ref and equivalent signatures based on user defined threshold for ROC AUC.


User defined arguments:
    `input_dir`
    `output_dir`
    `filtering by (i) Average AUC, (ii) minimum AUC`
    `filtering threshold`

"""
# TASK: iterate over all files in a given directory
# @stack: https://stackoverflow.com/questions/10377998/how-can-i-iterate-over-files-in-a-given-directory
# =============================================================================
# # Set output dir
# output_dir      = input("Please enter the full path to the output folder that your summary file will be stored:\n(A path to a directory must end with a '/')\n")#'C:/Users/bruno/Dropbox/dotpy/'
# input_dir       = input("Please enter the full path of the input folder that your JAD output files are stored:\n(A path to a directory must end with a '/')\n")# "C:/Users/bruno/Desktop/50repeats/"
# 
# AUC_filtering = int(input(
# 
# '''
# Select AUC filtering. Type:\n1 for filtering by AVERAGE AUC\n2 for filtering by LOWER AUC
# Type 1 or 2:  
# '''
# ))
# 
# 
# AUC_threshold = float(input(
# '''
# Type prefered AUC lower allowed threshold and use dot ('.') as a decimal separator (eg.: 0.5 )
# In the above example, contrasts with AUC < 0.5 will not be taken into account
# Now, type your threshold of choice:
# '''
#    ))
# =============================================================================
def JADbio_parser(input_dir, output_dir, AUC_filtering, AUC_threshold):
    print("_____________S T A R T _____________")

    # Ze imports
    from pandas import DataFrame    
    from pathlib import Path
    import time
    import sys
    from pathlib import Path

    # the Path call is a generator and fetches the filepaths on the fly
    # To store the filepaths as strings, call the generator object `pathlist`    
    # Save prints in report file
    
    
    pathlist  =     Path(input_dir).glob('**/*.txt')  
    
    # STORE JAD FILENPATHS IN A LIST
    JAD_files_list = [str(path) for path in pathlist]
    JAD_files_list = [ path.replace("\\", "/")  for path in JAD_files_list]
    
    # WHY -4? to remove '.txt', storing the name of the file which should point to which ctrl/case contrast it has
    ctrl_case_contrast = [x[len(input_dir):-4] for x in JAD_files_list]  
          
    # LOOP OVER ALL FILES IN DIRECTORY; 1 file is 1 contrast
    for i in range(0,len(JAD_files_list)):
        print("\n\n>looping over files in the directory..\n")
        print("\n")
    
    # =============================================================================
    # STEP 1:
    # From now on, iterate over all files in the input dir:
    # Remeber to return for looping over all files
    # =============================================================================
        filename        = ctrl_case_contrast[i]
        #JAD_output_file     = "C:/Users/bruno/Dropbox/HVNT/50repeats/amide_pos_corrected_AD_controlAD_controlSQ_controlSCLC.txt"
        JAD_output_file = JAD_files_list[i]
        
        # READING THE .txt FILE LINE BY LINE IN A LIST
        # leave out 1st line and ---- [1:]
        # leave out empty lines  ---- if x!="\n"
        
        # LOOP OVER ALL FILES IN THE INPUT DIR
    
        # Read 1 file of directory
        ### for JAD_output_file in JAD_files_list
        
        print("Reading file with contrast:   ",  JAD_files_list[i], "\n")
        print("File " + str(i+1) + " of " + str(len(JAD_files_list)))
        print("\n")

        # Read file    
        with open(JAD_output_file) as f:
            lis = [x for x in f if x!="\n"][1:]  
                                                # [1:] first row is the header "selected" so leaving it out
                                                 # for passing our data in a numpy data structure.                                            
    
       
        index_of_Equivalent_Signatures = lis.index('Equivalent Signatures \n') # index = 3
        index_of_Metric_Results        = lis.index('Metric results\n')
        
        # Pass for now, not using this in printed dataframe but could
        selected_model_configuration = lis[1][:-1]
        COLUMN_1 = selected_model_configuration
                
        # Area_Under_ROC_Curve 
        AUC_of_ROC_row           =    lis[index_of_Metric_Results+1] 
        index_of_OPEN_parenth    =    AUC_of_ROC_row.find('(')
        index_of_CLOSE_parenth   =    AUC_of_ROC_row.find(')')
        index_of_colon           =    AUC_of_ROC_row.find('e:')
        index_of_CI_colon        =    AUC_of_ROC_row.find('I:')
        index_of_comma           =    AUC_of_ROC_row.find(',')
        
        # Aftre grepping anchor chars, assign AUC values based on string slicicng
        AUC_of_ROC               =    round(float(AUC_of_ROC_row[index_of_colon    + 3 : index_of_OPEN_parenth]) , 3)
        Lower_AUC_of_ROC_limit   =    round(float(AUC_of_ROC_row[index_of_CI_colon + 3 : index_of_comma]     )   , 3)
        Upper_AUC_of_ROC_limit   =    round(float(AUC_of_ROC_row[index_of_comma    + 2 : index_of_CLOSE_parenth]), 3)
        
    
        # Accuracy
        accuracy_row             =    lis[index_of_Metric_Results+2]  
        index_of_OPEN_parenth    =    accuracy_row.find('(')
        index_of_CLOSE_parenth   =    accuracy_row.find(')')
        index_of_colon           =    accuracy_row.find('y:')
        index_of_CI_colon        =    accuracy_row.find('I:')
        index_of_comma           =    accuracy_row.find(',')
        accuracy                 =    round(float(accuracy_row[index_of_colon    + 3 : index_of_OPEN_parenth]) , 3)
        Lower_accuracy_limit     =    round(float(accuracy_row[index_of_CI_colon + 3 : index_of_comma]     )   , 3)
        Upper_accuracy_limit     =    round(float(accuracy_row[index_of_comma    + 2 : index_of_CLOSE_parenth]), 3)

        # =============================================================================
        # # Creating a list with the ref selected featureset plus their selection frequency
        # =============================================================================
        index_of_ref_sign   = lis.index("Reference Signature \n")    
        ref_signature       = lis[(index_of_ref_sign + 1)].split()        # selecting the line that corresponds to the reference signature
       
        index_of_Equivalent_Signatures = lis.index('Equivalent Signatures \n') # index = 3
        index_of_Metric_Results        = lis.index('Metric results\n')
    
        # LEN(equiv) + 1 to ommit the 'Equivalent Signatures' line
        # 1 line is one string, the list elements are strings
        equiv_signatures          = lis[index_of_Equivalent_Signatures+1:index_of_Metric_Results]#lis[len(JAD_6_from_start)+1:-len(JAD_13_from_end)]     
        
        # Split the strings to isolate each feature, omit the "\n" from the string
        equiv_signatures_LIST     = [x.split() for x in equiv_signatures if x!="\n"]   
        
        
        
        #VINCENZO's dataframe! 
        
        # ref in first position of list, then equivalent follow
        # Create empty LISTS of sublists, one sublist breakdwon: sublists[0] == ref feat, all the rest elements, equiv of the former ref feat
        
        refequiv_every_ref_its_equiv        = [ [] for j in range(0,len(ref_signature)) ]    
        
        # LIST of sublists; needed because the JAD output includes the ref feat in the equivalents, so dupes created
        EQUIV_unique                        = [ [] for j in range(0,len(ref_signature)) ]    
    
        # Use rownames, one row = one ref feat; will transpose later
        rownames                            = ["ref_feat_"+ str(ref_signature.index(x)+1) for x in ref_signature]
    
        # len (list of ref feats) = number of rows
        
        # LOOP OVER FEATURES OF REFERENCE SIGNATURE
        
        print("Info from file '" +str(ctrl_case_contrast[i] + "' retrieved." ))
        print("Looping over reference signature features ..")
        print("\n")
        
        for j in range(0,len(ref_signature)):
            
            print("reference feature   " + str(j + 1) + "   of " + str(len(ref_signature)) )
            
            EQUIV_unique[j] = [x for x in equiv_signatures_LIST[j] if x not in ref_signature[j]]
            
            # Add each feature from the ref signature in a sublist
            refequiv_every_ref_its_equiv[j].append(ref_signature[j])
            refequiv_every_ref_its_equiv[j] = refequiv_every_ref_its_equiv[j] + EQUIV_unique[j]
            # Unique but retain order
            
            # Add header name (eg ref_feat_1)in first index of sublist
            refequiv_every_ref_its_equiv[j].insert(0, rownames[j])   
        # Fill in with 0's to make dataframes of equal rows:
        
        # Find length of longest sublist:
            length_of_maxlength_sublist = len(max(refequiv_every_ref_its_equiv,key=len))
       
        # Exit forLoop; check if sublists are of same len(); if not, fill in with "NA"s
        
        # Append NA's until length equal in all sublists
        for j in range(0,len(ref_signature)):
            
            # DANGER ZONE! while! but it's oke, I checked
            while len(refequiv_every_ref_its_equiv[j]) != length_of_maxlength_sublist:
                refequiv_every_ref_its_equiv[j].append("NA")
       
        # Create "ref" or "equiv" flag column to features
        ref_equiv_flag        = ["flag"]+["ref"] +  ['equiv'] * (len(refequiv_every_ref_its_equiv[0])-2)
    
        df_ref_equiv_assigned = DataFrame(refequiv_every_ref_its_equiv)
    
        df_ref_equiv_assigned = df_ref_equiv_assigned.T
    
        # Add "ref" or "equiv" flag to features
        df_ref_equiv_assigned.insert(0, "flag", ref_equiv_flag)
       
        # Add column to the end of dataframe: average AUC
        avgAUC = ["avgAUC"] + [AUC_of_ROC] * ((len(ref_equiv_flag)) - 1)
        df_ref_equiv_assigned.insert(len(df_ref_equiv_assigned.columns), "avgAUC", avgAUC)
        
        # Add column to the end of dataframe: min AUC
        minAUC = ["minAUC"] + [Lower_AUC_of_ROC_limit] * ((len(ref_equiv_flag)) - 1)
        df_ref_equiv_assigned.insert(len(df_ref_equiv_assigned.columns), "minAUC", minAUC )
            
        # Add column to the end of dataframe: MAX AUC
        maxAUC = ["maxAUC"] + [Upper_AUC_of_ROC_limit] * ((len(ref_equiv_flag)) - 1)
        df_ref_equiv_assigned.insert(len(df_ref_equiv_assigned.columns), "maxAUC", maxAUC )
        
        # Replace header with first row
        new_header = df_ref_equiv_assigned.iloc[0] #grab the first row for the header
        df_ref_equiv_assigned = df_ref_equiv_assigned[1:] #take the data less the header row
        df_ref_equiv_assigned.columns = new_header #set the header row as the df head
        
        
        if AUC_filtering == 1:
            filtering_threshold   = df_ref_equiv_assigned['avgAUC'][1]
                 
        elif AUC_filtering == 2:
            filtering_threshold = df_ref_equiv_assigned['minAUC'][1]
    
        if filtering_threshold > AUC_threshold:
            # Write dataframes in csv
            df_ref_equiv_assigned.to_csv(output_dir + filename + ".csv",    # exact filepath to the filename
                  
                  sep         = ','           ,    #.csv file
                  encoding    = 'utf-8'       , 
                  header      = True         ,
                  index       = False         
                  )        
            
            print(df_ref_equiv_assigned)
        else:
            print("Contrast   " +filename+ "   didn't make the AUC threshold")
            print("too low AUC at " + str(filtering_threshold))
            print("pass and continue with remaining contrasts")
            pass
    print("\n\nFinished!\n\nYour output files are now stored in the folder:\n" + str(output_dir))
    print("_____________ T h e     E N D _____________")
    
    
    
    #%% Example Usage
    
AUC_filtering = 2 # filter by miin AUC
AUC_threshold = 0.47

JADbio_parser("C:/Users/bruno/Dropbox/2_dataNov2016/JADBio_reports_96/Report_c18neg_minAUC0dot47_but_avgAUC_0dot647/", 
              "C:/Users/bruno/Dropbox/2_dataNov2016/JADBio_Signatures_metabo/extraContrast_NSCLC_minAUC0dot47_avgAUC0dot647/", 
              AUC_filtering, 
              AUC_threshold)