import os
import pandas as pd
import argparse

pathotype_path = "/g/data/xf3/ht5438/data/Pathotype_book/20210412_pathotype_book.xlsx"


parser = argparse.ArgumentParser(description="Filter out Yr genes with percentage of phenotype data lower than a specified threshold, e.g. 0.1.")
parser.add_argument("-p", "--pathotype_path", metavar="", required=True, help="Path to pathotype book excel sheet.")
parser.add_argument("-t", "--threshold", metavar="", required=True, type=float, help="Percentage threshold to filter out Yr genes with proportion of available phenotype data lower than this value.")
args = parser.parse_args()

def YrFilter(pathotype_path, filter_per):
    assert filter_per <= 1.0
    pb_df = pd.read_excel(pathotype_path, sheet_name="Pathotype_book", engine="openpyxl").iloc[:,3:]
    sum_list = list()
    pb_df.columns = pb_df.columns.astype(str)
    for col in pb_df.columns:
        sum_list.append(pb_df.loc[:,col].value_counts())
    threshold = len(pb_df.index)*float(filter_per)
    
    Filtered_Yr_List = [i.name for i in sum_list if sum(i) >= threshold]
    Removed_Yr_List = [i.name for i in sum_list if sum(i) < threshold]
    
    print(f"\n{'='*30}\n{len(Filtered_Yr_List)} Yr gene(s) passed the {filter_per*100}% threshold: {', '.join(Filtered_Yr_List)}")
    print("\nYAML format:")
    Yr_name = {}
    for yr in Filtered_Yr_List:
        Yr_name[f"Yr{yr}"] = f"Yr{yr}"
    for key,val in Yr_name.items():
        print(f'{key}:"{val}"')
    print("="*30)
    print(f"\n{len(Removed_Yr_List)} Yr gene removed: {', '.join(Removed_Yr_List)}\n")
        

if __name__ == "__main__":
    YrFilter(args.pathotype_path, args.threshold)
  
    
    

    
