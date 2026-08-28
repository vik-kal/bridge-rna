from bridge_rna import cohorts as C
from bridge_rna.retrieval import run_cohort_retrieval
from bridge_rna.geo import _enrich_hits_from_ncbi_eutils
from bridge_rna.layout import samples_df
from functools import reduce
import pandas as pd
import time
import os
import chime

# Hardcode the study you want
  
topk = 10    
email_value = "kalahasthivikasni@gmail.com"      


def run_cohort_dataframing(study_id): #returns a dataframe with the topk hits for space and ground, combined into one df
    # Build all cohorts for this one study
    facets = ['study','spaceflight']
    cohort_list = C.build_cohorts(facets=facets, study=study_id)

    print(f"Found {len(cohort_list)} cohorts for study {study_id} with facets {facets}")


    df_list = []
    merged_df = pd.DataFrame()  # Initialize an empty DataFrame to hold merged results
    # Loop through each cohort and retrieve hits
    for cohort in cohort_list:
        if cohort.size < C.MIN_COHORT_SIZE:
            print(f"Skipping cohort {cohort.label}: only {cohort.size} sample(s)")
            continue

        members = list(cohort.members)
        #print(cohort.label)


        hits_df, rows, stability = run_cohort_retrieval(members, topk=topk)
        hits_df = _enrich_hits_from_ncbi_eutils(hits_df, email_value)
    

        altered_df = hits_df[['gsm','gse','geo_summary']].copy()
        altered_df['spaceflight'] = cohort.label



        df_list.append(altered_df)

        #print(altered_df.head())  # preview altered dataframe
    #merge the dataframes for all cohort_options while making the cohort facet grouping as a new column

    if df_list:
        merged_df = pd.concat(df_list,ignore_index=True)
        #print(merged_df)
    
    return merged_df
    


def loop_all_cohorts():

    all_studies = samples_df["study_id"].unique()

    for study in all_studies:
        merged_data = run_cohort_dataframing(study)
        if not merged_data.empty:
            merged_data.to_csv( f"archs4metadata_cohort/{study}_hits.csv")

def test_one_cohort(study):
       
    merged_data = run_cohort_dataframing(study)
    if not merged_data.empty:
        merged_data.to_csv( f"archs4metadata_cohort/{study}_hits.csv")

if __name__ == "__main__":
    start = time.time()
    loop_all_cohorts()
    #print('hi')
    #test_one_cohort("OSD-141")
    
    end = time.time()
    print("Execution time:", end - start, "seconds")
    