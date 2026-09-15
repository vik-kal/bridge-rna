from bridge_rna import cohorts as C
from bridge_rna.retrieval import run_cohort_retrieval
from bridge_rna.geo import _enrich_hits_from_ncbi_eutils
from bridge_rna.retrieval import run_cached_query_retrieval
from bridge_rna.layout import samples_df
from demo_osdr_top5 import fetch_archs4_metadata
import pandas as pd
from pathlib import Path
import time
import os
from collections import defaultdict



# Hardcode the study you want
  
topk = 10    

email_value = os.getenv("email_value")


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
        print(cohort.members)
        #print(cohort.label)


        hits_df, rows, stability = run_cohort_retrieval(members, topk=topk)
        hits_df = _enrich_hits_from_ncbi_eutils(hits_df, email_value)
    

        #altered_df = hits_df[['gsm','gse','geo_summary']].copy()
        #altered_df = hits_df.copy()
        hits_df['spaceflight'] = cohort.label



        df_list.append(hits_df)

        #print(altered_df.head())  # preview altered dataframe
    #merge the dataframes for all cohort_options while making the cohort facet grouping as a new column

    if df_list:
        merged_df = pd.concat(df_list,ignore_index=True)
        #print(merged_df)
    
    return merged_df

def run_alt_cohort_dataframing(study_id):
    merged_df = pd.DataFrame()
    df_list = []
    cohort_dict = get_cohort_members(study_id)

    #test_dict = {'space':['OSD-244|Mmus_C57-6T_TMS_BSL_ISS-T_Rep2_B2','OSD-244|Mmus_C57-6T_TMS_GC_LAR_Rep8_G9']}
    

    for k,v in cohort_dict.items():

        for sample_name in v:
            hits_df = run_cached_query_retrieval(sample_name,10)
            hits_df = _enrich_hits_from_ncbi_eutils(hits_df,email_value)
            hits_df['sample_name'] = sample_name
            hits_df['spaceflight'] = k

            #add to df_list to merge
            df_list.append(hits_df)
        
    if df_list:
            merged_df = pd.concat(df_list,ignore_index=True)    

    return merged_df




def add_metadata_to_hits():
    human_archs4_path = Path('/media/volume/H5-Files/archs4/human_gene_v2.latest.h5')
    mouse_archs4_path = Path('/media/volume/H5-Files/archs4/mouse_gene_v2.latest.h5')

    hit_path = Path('archs4metadata_cohort_noncbi')
    
    for f in hit_path.iterdir():
        if str(f).endswith('.csv'):
            daf = pd.read_csv(f) #test with one first
            geo_list = daf['gsm'].tolist()
            metadata_df = fetch_archs4_metadata(geo_list,human_archs4_path,mouse_archs4_path)
            daf['characteristics'] = metadata_df['characteristics_ch1']
           
            
            daf.to_csv(f, index=False)
            print(f"Saved:{f}")


def get_cohort_members(study_id): # returns a dict of cohort members in each group

        member_dict = defaultdict(list)
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
            member_dict[cohort.label] = members
        return member_dict



def loop_all_cohorts():

    all_studies = samples_df["study_id"].unique()

    for study in all_studies:
        merged_data = run_alt_cohort_dataframing(study)
        if not merged_data.empty:
            merged_data.to_csv( f"archs4data/archs4hitscohort_alt_pre/{study}_hits.csv")

def test_one_cohort(study):
       
    merged_data = run_alt_cohort_dataframing(study)

    if not merged_data.empty:
        merged_data.to_csv( f"archs4data/archs4hitscohort_alt_pre/{study}_hits.csv")




if __name__ == "__main__":
    start = time.time()
    loop_all_cohorts()
    #add_metadata_to_hits()
    #print('hi')
    #test_one_cohort("OSD-464")
    #print(get_cohort_members("OSD-244"))
    #print(run_alt_cohort_dataframing("OSD-244"))
    
    end = time.time()
    print("Execution time:", end - start, "seconds")
    